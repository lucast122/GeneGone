import streamlit as st
import os
import tempfile
import pandas as pd
from utils import (
    ReadType,
    run_alignment,
    process_alignment,
    calculate_coverage,
    plot_gene_coverage,
    create_deletion_summary,
    analyze_deletion,
    CPU_CORES
)

st.set_page_config(page_title="GeneGone", layout="wide")

st.title("GeneGone: Gene Deletion Analysis Tool")

# Show system info in sidebar
with st.sidebar:
    st.header("System Information")
    st.info(f"Using {CPU_CORES} CPU cores for analysis")
    
    st.header("Analysis Configuration")
    read_type = st.selectbox(
        "Select Sequencing Type",
        [rt.value for rt in ReadType],
        format_func=lambda x: x.capitalize()
    )
    
    # Coverage threshold settings
    st.subheader("Coverage Settings")
    use_auto_threshold = st.checkbox(
        "Auto-calculate threshold",
        value=True,
        help="Automatically calculate threshold based on flanking coverage"
    )
    
    if not use_auto_threshold:
        coverage_threshold = st.number_input(
            "Coverage Threshold",
            min_value=1,
            value=20 if read_type == "illumina" else 10,
            help="Coverage below this value indicates deletion"
        )
    else:
        coverage_threshold = None
    
    min_mapping_quality = st.slider(
        "Minimum Mapping Quality",
        min_value=0,
        max_value=60,
        value=20 if read_type == "illumina" else 0,
        help="Higher values are more stringent. Recommended: 20 for Illumina, 0 for Nanopore/PacBio"
    )

# File upload section
st.header("Upload Files")
col1, col2 = st.columns(2)

with col1:
    reference_file = st.file_uploader("Upload Reference Genome (FASTA)", type=['fasta', 'fa'])
    if reference_file:
        st.success("✅ Reference genome uploaded")

with col2:
    reads_file = st.file_uploader("Upload Sequencing Reads (FASTQ)", type=['fastq', 'fq'])
    if reads_file:
        st.success("✅ Reads file uploaded")

# Region selection
st.header("Target Region")
col1, col2, col3 = st.columns(3)

with col1:
    chromosome = st.text_input("Chromosome Name", value="test_chr1")
with col2:
    start_pos = st.number_input("Start Position", value=4001, min_value=1)
with col3:
    end_pos = st.number_input("End Position", value=5000, min_value=1)

gene_name = st.text_input("Gene Name (optional)", value="Target Gene")

# Analysis section
if st.button("Run Analysis"):
    if not (reference_file and reads_file):
        st.error("Please upload both reference genome and reads files")
    else:
        with st.spinner("Processing files..."):
            # Create temporary directory
            with tempfile.TemporaryDirectory() as tmpdir:
                # Save uploaded files
                ref_path = os.path.join(tmpdir, "reference.fasta")
                reads_path = os.path.join(tmpdir, "reads.fastq")
                
                with open(ref_path, "wb") as f:
                    f.write(reference_file.getvalue())
                with open(reads_path, "wb") as f:
                    f.write(reads_file.getvalue())
                
                # Run alignment
                sam_path = os.path.join(tmpdir, "aligned.sam")
                bam_path = os.path.join(tmpdir, "aligned.bam")
                
                with st.status("Running analysis...", expanded=True) as status:
                    status.write("Aligning reads...")
                    if not run_alignment(reads_path, ref_path, sam_path, ReadType(read_type)):
                        st.error("Alignment failed")
                        st.stop()
                    
                    status.write("Processing alignment...")
                    if not process_alignment(sam_path, bam_path, ReadType(read_type)):
                        st.error("BAM processing failed")
                        st.stop()
                    
                    status.write("Calculating coverage...")
                    region = f"{chromosome}:{start_pos}-{end_pos}"
                    coverage_df = calculate_coverage(
                        bam_path, 
                        region, 
                        ReadType(read_type),
                        min_mapping_quality
                    )
                    
                    if coverage_df.empty:
                        st.error("Coverage calculation failed")
                        st.stop()
                    
                    # Analyze deletion
                    analysis_results = analyze_deletion(
                        coverage_df,
                        start_pos,
                        end_pos,
                        ReadType(read_type),
                        coverage_threshold
                    )
                    
                    status.write("Creating visualizations...")
                    status.update(label="Analysis complete!", state="complete")

        # Results section
        st.header("Analysis Results")
        
        # Show key metrics in a clean row
        col1, col2, col3 = st.columns(3)
        with col1:
            deletion_status = "DELETION DETECTED" if analysis_results['deletion_detected'] else "NO DELETION"
            status_color = "red" if analysis_results['deletion_detected'] else "green"
            st.markdown(f"<h3 style='text-align: center; color: {status_color};'>{deletion_status}</h3>", unsafe_allow_html=True)
        with col2:
            st.metric(
                "Coverage Reduction",
                f"{analysis_results['deletion_percentage']:.1f}%",
                delta=None
            )
        with col3:
            st.metric(
                "Mean Coverage",
                f"{analysis_results['gene_mean_coverage']:.1f}x",
                delta=f"{analysis_results['flank_mean_coverage']:.1f}x flanking"
            )
        
        # Show coverage plot
        fig = plot_gene_coverage(
            coverage_df,
            start_pos,
            end_pos,
            ReadType(read_type),
            gene_name,
            analysis_results
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # Optional: Add download buttons in a single row
        col1, col2 = st.columns(2)
        with col1:
            st.download_button(
                "Download Plot",
                fig.to_image(format="png"),
                "coverage_plot.png",
                "image/png"
            )
        with col2:
            # Export coverage data
            csv = coverage_df.to_csv(index=False)
            st.download_button(
                "Download Coverage Data",
                csv,
                "coverage_data.csv",
                "text/csv"
            )

# Help section
with st.expander("Help & Documentation"):
    st.markdown("""
    ### How to Use GeneGone
    1. **Upload Files**:
       - Reference genome in FASTA format
       - Sequencing reads in FASTQ format
    
    2. **Configure Analysis**:
       - Select your sequencing platform
       - Choose between automatic or manual coverage threshold
       - Adjust mapping quality if needed
       - Enter the target gene coordinates
    
    3. **Interpret Results**:
       - Coverage Plot: Shows read depth across the region with context
       - Red highlighted areas indicate potential deletions
       - Green dashed line shows mean flanking coverage
       - Red dashed line shows coverage threshold
    
    ### Coverage Settings
    - **Auto-calculate threshold**: Uses flanking regions to determine threshold
    - **Manual threshold**: Set a specific coverage cutoff
    - **Mapping Quality**:
      - Illumina: Use higher values (20+)
      - Nanopore/PacBio: Lower values (0-10) are acceptable
    """)
