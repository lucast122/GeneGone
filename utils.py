import subprocess
import os
import shutil
import pandas as pd
import numpy as np
from Bio import SeqIO
import plotly.graph_objects as go
from enum import Enum
from typing import Dict, Optional, Tuple
import multiprocessing

# Get number of available CPU cores
CPU_CORES = multiprocessing.cpu_count()

class ReadType(Enum):
    NANOPORE = "nanopore"
    PACBIO = "pacbio"
    ILLUMINA = "illumina"

def get_aligner_params(read_type: ReadType) -> Tuple[str, list]:
    """
    Get aligner command and parameters based on read type.
    
    Args:
        read_type: Type of sequencing reads
    
    Returns:
        Tuple of (command, parameters)
    """
    if read_type == ReadType.NANOPORE:
        return "minimap2", ["-ax", "map-ont", "-t", str(CPU_CORES)]
    elif read_type == ReadType.PACBIO:
        return "minimap2", ["-ax", "map-pb", "-t", str(CPU_CORES)]
    else:  # Illumina
        return "bwa", ["mem", "-t", str(CPU_CORES)]

def run_alignment(fastq_path: str, reference_path: str, output_sam: str, read_type: ReadType) -> bool:
    """
    Align reads to reference using appropriate aligner for read type.
    
    Args:
        fastq_path: Path to input FASTQ file
        reference_path: Path to reference FASTA file
        output_sam: Path for output SAM file
        read_type: Type of sequencing reads
    
    Returns:
        bool: True if alignment successful, False otherwise
    """
    try:
        aligner, params = get_aligner_params(read_type)
        
        if read_type == ReadType.ILLUMINA:
            # For BWA, we need to ensure the index exists
            if not os.path.exists(f"{reference_path}.bwt"):
                subprocess.run(["bwa", "index", reference_path], check=True)
        
        cmd = [aligner] + params + [reference_path, fastq_path]
        with open(output_sam, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        return True
    except subprocess.CalledProcessError:
        return False

def process_alignment(sam_path: str, output_bam: str, read_type: ReadType) -> bool:
    """
    Convert SAM to sorted BAM file with read-type specific processing.
    
    Args:
        sam_path: Path to input SAM file
        output_bam: Path for output BAM file
        read_type: Type of sequencing reads
    
    Returns:
        bool: True if processing successful, False otherwise
    """
    try:
        # Convert SAM to BAM using multiple threads
        subprocess.run([
            "samtools", "view", 
            "-@", str(CPU_CORES-1),  # Leave one core for system
            "-bS", sam_path, 
            "-o", "temp.bam"
        ], check=True)
        
        # Sort BAM using multiple threads
        subprocess.run([
            "samtools", "sort",
            "-@", str(CPU_CORES-1),
            "temp.bam",
            "-o", output_bam
        ], check=True)
        
        # Index BAM
        subprocess.run([
            "samtools", "index",
            "-@", str(CPU_CORES-1),
            output_bam
        ], check=True)
        
        # Additional processing for specific read types
        if read_type in [ReadType.NANOPORE, ReadType.PACBIO]:
            # For long reads, mark supplementary alignments
            subprocess.run([
                "samtools", "view",
                "-@", str(CPU_CORES-1),
                "-h", "-F", "0x800",
                output_bam,
                "-o", "temp_filtered.bam"
            ], check=True)
            shutil.move("temp_filtered.bam", output_bam)
            subprocess.run([
                "samtools", "index",
                "-@", str(CPU_CORES-1),
                output_bam
            ], check=True)
        
        # Clean up
        if os.path.exists("temp.bam"):
            os.remove("temp.bam")
        return True
    except subprocess.CalledProcessError:
        return False

def calculate_coverage(bam_path: str, region: str, read_type: ReadType, min_mapping_quality: int = None, 
                      include_flanking: bool = True, flank_size: int = 100) -> pd.DataFrame:
    """
    Calculate coverage for a specific genomic region with read-type specific parameters.
    
    Args:
        bam_path: Path to BAM file
        region: Region string (e.g., "chr1:1000-2000")
        read_type: Type of sequencing reads
        min_mapping_quality: Minimum mapping quality score (0-60). If None, uses platform-specific defaults:
            - Illumina: 20 (higher confidence due to shorter reads)
            - Nanopore/PacBio: 0 (length provides confidence despite errors)
        include_flanking: Whether to automatically include flanking regions
        flank_size: Size of flanking regions (in bp) to include on each side
    
    Returns:
        DataFrame with position and coverage
    """
    try:
        # Set default mapping quality based on read type if not provided
        if min_mapping_quality is None:
            min_mapping_quality = 20 if read_type == ReadType.ILLUMINA else 0
        
        # Convert mapping quality to string for samtools
        min_qual = str(min_mapping_quality)
        
        # Expand region to include flanking sequences if requested
        if include_flanking and ':' in region and '-' in region:
            # Parse region string (format: "chrom:start-end")
            chrom, pos_range = region.split(':')
            start, end = map(int, pos_range.split('-'))
            
            # Expand region to include flanking areas
            expanded_start = max(1, start - flank_size)  # Ensure start is at least 1
            expanded_end = end + flank_size
            
            # Create expanded region string
            region = f"{chrom}:{expanded_start}-{expanded_end}"
            
            print(f"Expanded coverage region to include flanking: {region}")
        
        # Use multiple threads for coverage calculation
        cmd = [
            "samtools", "depth",
            "-@", str(CPU_CORES-1),
            "-r", region,
            "-q", min_qual,
            bam_path
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        coverage_data = []
        for line in result.stdout.split('\n'):
            if line:
                chrom, pos, depth = line.split()
                coverage_data.append({
                    'position': int(pos),
                    'coverage': int(depth)
                })
        
        df = pd.DataFrame(coverage_data)
        
        # Add metadata about the analysis
        df.attrs['min_mapping_quality'] = min_mapping_quality
        df.attrs['read_type'] = read_type.value
        
        return df
    except subprocess.CalledProcessError:
        return pd.DataFrame()

def analyze_deletion(coverage_df: pd.DataFrame, gene_start: int, gene_end: int, 
                    read_type: ReadType, coverage_threshold: Optional[float] = None) -> dict:
    """
    Analyze deletion by comparing coverage in gene region vs flanking regions.
    
    Args:
        coverage_df: DataFrame with position and coverage data
        gene_start: Start position of the gene
        gene_end: End position of the gene
        read_type: Type of sequencing reads
        coverage_threshold: Optional custom coverage threshold
        
    Returns:
        Dictionary with deletion analysis results
    """
    if coverage_df.empty:
        return {
            'deletion_detected': False,
            'deletion_percentage': 0,
            'gene_mean_coverage': 0,
            'flank_mean_coverage': 0,
            'threshold_used': 0,
            'deletion_ratio': 1
        }

    # Use fixed 100bp flanking regions
    flank_size = 100
    
    # Get the actual ranges of positions in the data
    min_pos = min(coverage_df['position'])
    max_pos = max(coverage_df['position'])
    
    # Calculate gene region statistics
    gene_region = coverage_df[
        (coverage_df['position'] >= gene_start) & 
        (coverage_df['position'] <= gene_end)
    ]
    gene_mean = gene_region['coverage'].mean() if not gene_region.empty else 0
    
    # Calculate the left and right flanking region boundaries
    left_flank_start = max(min_pos, gene_start - flank_size)
    left_flank_end = max(min_pos, gene_start - 1)  # Stop before gene starts
    
    right_flank_start = min(max_pos, gene_end + 1)  # Start after gene ends
    right_flank_end = min(max_pos, gene_end + flank_size)
    
    # Extract the flanking regions from the coverage data
    left_flank = coverage_df[
        (coverage_df['position'] >= left_flank_start) & 
        (coverage_df['position'] <= left_flank_end)
    ]
    
    right_flank = coverage_df[
        (coverage_df['position'] >= right_flank_start) & 
        (coverage_df['position'] <= right_flank_end)
    ]
    
    # Calculate flanking statistics
    left_mean = left_flank['coverage'].mean() if not left_flank.empty else 0
    right_mean = right_flank['coverage'].mean() if not right_flank.empty else 0
    
    # Print debug information
    print(f"Debug - Gene region: {gene_start}-{gene_end}, Mean coverage: {gene_mean}")
    print(f"Debug - Left flank: {left_flank_start}-{left_flank_end}, Mean coverage: {left_mean}")
    print(f"Debug - Right flank: {right_flank_start}-{right_flank_end}, Mean coverage: {right_mean}")
    
    # Use available flanking regions for mean calculation
    if left_flank.empty and right_flank.empty:
        print("WARNING: No flanking regions found in data. Cannot calculate flanking coverage.")
        flank_mean = 0
    elif left_flank.empty:
        flank_mean = right_mean
    elif right_flank.empty:
        flank_mean = left_mean
    else:
        flank_mean = (left_mean + right_mean) / 2
    
    # Ensure we have a valid flanking mean
    if flank_mean == 0:
        print("WARNING: Zero flanking coverage. Using 1 as default to avoid division by zero.")
        flank_mean = 1  # Avoid division by zero
    
    # Set threshold based on read type if not provided
    if coverage_threshold is None:
        if read_type == ReadType.ILLUMINA:
            coverage_threshold = max(20, flank_mean * 0.2)  # At least 20x or 20% of flanking coverage
        else:
            coverage_threshold = max(10, flank_mean * 0.2)  # At least 10x or 20% of flanking coverage
    
    # Calculate deletion metrics
    deletion_ratio = gene_mean / flank_mean
    
    # Calculate reduction percentage
    # If gene coverage is lower than flanking, calculate reduction
    if gene_mean < flank_mean:
        deletion_percentage = ((flank_mean - gene_mean) / flank_mean) * 100
    else:
        deletion_percentage = 0
    
    print(f"Debug - Flanking mean: {flank_mean}, Gene mean: {gene_mean}")
    print(f"Debug - Deletion ratio: {deletion_ratio}, Deletion percentage: {deletion_percentage}")
    
    return {
        'deletion_detected': gene_mean < coverage_threshold,
        'deletion_percentage': deletion_percentage,
        'gene_mean_coverage': gene_mean,
        'flank_mean_coverage': flank_mean,
        'threshold_used': coverage_threshold,
        'deletion_ratio': deletion_ratio
    }

def plot_gene_coverage(coverage_df: pd.DataFrame, gene_start: int, gene_end: int, 
                      read_type: ReadType, gene_name: str = "Target Gene",
                      analysis_results: Optional[dict] = None) -> go.Figure:
    """
    Create a publication-quality coverage plot with gene annotation and deletion visualization.
    
    Args:
        coverage_df: DataFrame with position and coverage data
        gene_start: Start position of the gene
        gene_end: End position of the gene
        read_type: Type of sequencing reads
        gene_name: Name of the target gene
        analysis_results: Optional dictionary with deletion analysis results
        
    Returns:
        Plotly figure object configured for publication-quality output
    """
    # Calculate plot boundaries to show context
    gene_size = gene_end - gene_start
    context_size = max(gene_size * 2, 2000)  # Show at least 2kb or 2x gene size
    plot_start = max(min(coverage_df['position']) - context_size, 0)
    plot_end = max(coverage_df['position']) + context_size
    
    # Create figure with white background
    fig = go.Figure()
    
    # Define publication-quality color scheme
    colors = {
        'nanopore': '#1f77b4',  # Professional blue
        'illumina': '#2ca02c',  # Professional green
        'gene': '#d62728',      # Dark red
        'threshold': '#7f7f7f',  # Gray
        'flanking': '#17becf',  # Teal
        'deletion': 'rgba(214, 39, 40, 0.15)'  # Transparent red
    }
    
    # Select color based on read type
    coverage_color = colors['nanopore'] if read_type == ReadType.NANOPORE else colors['illumina']
    
    # Add coverage trace with refined styling
    fig.add_trace(go.Scatter(
        x=coverage_df['position'],
        y=coverage_df['coverage'],
        name='Coverage',
        line=dict(
            color=coverage_color,
            width=1.5
        ),
        fill='tozeroy',
        fillcolor=f'rgba{tuple(int(coverage_color.lstrip("#")[i:i+2], 16) for i in (0, 2, 4)) + (0.2,)}'
    ))
    
    # Calculate y position for gene annotation
    max_coverage = max(coverage_df['coverage'])
    gene_y = max_coverage * 1.1
    
    # Add gene box annotation with refined styling
    fig.add_shape(
        type="rect",
        x0=gene_start,
        x1=gene_end,
        y0=gene_y - max_coverage * 0.015,
        y1=gene_y + max_coverage * 0.015,
        line=dict(color=colors['gene'], width=1.5),
        fillcolor=colors['gene'],
    )
    
    # Add gene label with refined styling
    fig.add_annotation(
        x=(gene_start + gene_end) / 2,
        y=gene_y + max_coverage * 0.05,
        text=f"<b>{gene_name}</b>",
        showarrow=False,
        font=dict(size=12, color="black", family="Arial"),
        bgcolor="white",
        bordercolor=colors['gene'],
        borderwidth=1.5,
        borderpad=3
    )
    
    if analysis_results:
        # Add threshold line with refined styling
        fig.add_hline(
            y=analysis_results['threshold_used'],
            line_dash="dash",
            line_color=colors['threshold'],
            line_width=1,
            annotation=dict(
                text=f"Coverage Threshold ({analysis_results['threshold_used']:.1f}x)",
                font=dict(color=colors['threshold'], size=10, family="Arial"),
                bgcolor="white",
                bordercolor=colors['threshold'],
                borderwidth=1,
                borderpad=3
            ),
            annotation_position="right"
        )
        
        # Add mean coverage lines with refined styling
        fig.add_hline(
            y=analysis_results['flank_mean_coverage'],
            line_dash="dash",
            line_color=colors['flanking'],
            line_width=1,
            annotation=dict(
                text=f"Mean Flanking Coverage ({analysis_results['flank_mean_coverage']:.1f}x)",
                font=dict(color=colors['flanking'], size=10, family="Arial"),
                bgcolor="white",
                bordercolor=colors['flanking'],
                borderwidth=1,
                borderpad=3
            ),
            annotation_position="right"
        )
        
        # Add deletion indicator if detected
        if analysis_results['deletion_detected']:
            # Add semi-transparent red overlay with refined styling
            fig.add_shape(
                type="rect",
                x0=gene_start,
                x1=gene_end,
                y0=0,
                y1=max_coverage,
                fillcolor=colors['deletion'],
                line=dict(color=colors['gene'], dash="dash", width=1),
                layer="below"
            )
            
            # Add deletion text box with refined styling
            fig.add_annotation(
                x=(gene_start + gene_end) / 2,
                y=max_coverage * 0.5,
                text=f"<b>DELETION DETECTED</b><br>{analysis_results['deletion_percentage']:.1f}% reduction",
                showarrow=False,
                font=dict(size=14, color=colors['gene'], family="Arial"),
                bgcolor="rgba(255, 255, 255, 0.9)",
                bordercolor=colors['gene'],
                borderwidth=1.5,
                borderpad=4,
                align="center"
            )
    
    # Update layout with publication-quality styling
    fig.update_layout(
        title=dict(
            text=f"Gene Deletion Analysis ({read_type.value.capitalize()} Data)",
            x=0.5,
            font=dict(size=16, family="Arial", color="black")
        ),
        xaxis=dict(
            title=dict(
                text="Genomic Position",
                font=dict(size=12, family="Arial", color="black")
            ),
            range=[plot_start, plot_end],
            showgrid=True,
            gridcolor='rgba(0,0,0,0.1)',
            gridwidth=0.5,
            zeroline=False,
            tickfont=dict(size=10, family="Arial", color="black"),
            mirror=True,
            showline=True,
            linecolor='black',
            linewidth=1
        ),
        yaxis=dict(
            title=dict(
                text="Coverage Depth",
                font=dict(size=12, family="Arial", color="black")
            ),
            rangemode="tozero",
            showgrid=True,
            gridcolor='rgba(0,0,0,0.1)',
            gridwidth=0.5,
            zeroline=False,
            tickfont=dict(size=10, family="Arial", color="black"),
            mirror=True,
            showline=True,
            linecolor='black',
            linewidth=1
        ),
        showlegend=False,  # Hide legend for cleaner look
        hovermode='x unified',
        width=800,  # Fixed width for publication
        height=500,  # Fixed height for publication
        template="none",  # No template for clean look
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(t=60, r=120, b=60, l=60)  # Balanced margins
    )
    
    # Add watermark for publication source
    fig.add_annotation(
        x=1.0,
        y=-0.15,
        xref="paper",
        yref="paper",
        text="GeneGone Analysis",
        showarrow=False,
        font=dict(size=8, color="gray", family="Arial"),
        align="right",
    )
    
    return fig

def create_deletion_summary(coverage_df: pd.DataFrame, gene_start: int, gene_end: int,
                          threshold: int, read_type: ReadType) -> go.Figure:
    """
    Create a publication-quality summary visualization of the deletion analysis.
    
    Args:
        coverage_df: DataFrame with coverage data
        gene_start: Start position of the gene
        gene_end: End position of the gene
        threshold: Coverage threshold for deletion detection
        read_type: Type of sequencing reads
    
    Returns:
        Plotly figure object formatted for publication
    """
    # Calculate metrics
    gene_region = coverage_df[
        (coverage_df['position'] >= gene_start) & 
        (coverage_df['position'] <= gene_end)
    ]
    
    total_positions = len(gene_region)
    positions_below_threshold = len(gene_region[gene_region['coverage'] < threshold])
    deletion_percentage = (positions_below_threshold / total_positions) * 100
    
    mean_coverage = gene_region['coverage'].mean()
    max_coverage = gene_region['coverage'].max()
    
    # Define publication color scheme
    colors = {
        'gauge_low': '#f0f0f0',     # Light gray
        'gauge_medium': '#ffcc99',  # Light orange
        'gauge_high': '#c6e2ff',    # Light blue
        'gauge_bar': '#d62728',     # Dark red
        'text': '#333333',          # Dark gray
        'border': '#666666'         # Medium gray
    }
    
    # Create gauge chart with refined styling
    fig = go.Figure(go.Indicator(
        mode = "gauge+number",
        value = deletion_percentage,
        title = {
            'text': "Deletion Confidence", 
            'font': {'size': 16, 'family': 'Arial', 'color': colors['text']}
        },
        number = {
            'suffix': "%",
            'font': {'size': 20, 'family': 'Arial', 'color': colors['text']}
        },
        gauge = {
            'axis': {
                'range': [None, 100], 
                'tickwidth': 1, 
                'tickcolor': colors['text'],
                'tickfont': {'size': 10, 'family': 'Arial'}
            },
            'bar': {'color': colors['gauge_bar']},
            'bgcolor': "white",
            'borderwidth': 1,
            'bordercolor': colors['border'],
            'steps': [
                {'range': [0, 40], 'color': colors['gauge_low']},
                {'range': [40, 80], 'color': colors['gauge_medium']},
                {'range': [80, 100], 'color': colors['gauge_high']}
            ],
            'threshold': {
                'line': {'color': colors['text'], 'width': 2},
                'thickness': 0.75,
                'value': 80
            }
        }
    ))
    
    # Add metrics as annotations with refined styling
    annotations = [
        {
            'x': 0.15, 'y': -0.15,
            'text': f"Mean Coverage: {mean_coverage:.1f}x",
            'showarrow': False,
            'font': {'size': 10, 'family': 'Arial', 'color': colors['text']},
            'xref': 'paper', 'yref': 'paper'
        },
        {
            'x': 0.5, 'y': -0.15,
            'text': f"Max Coverage: {max_coverage:.1f}x",
            'showarrow': False,
            'font': {'size': 10, 'family': 'Arial', 'color': colors['text']},
            'xref': 'paper', 'yref': 'paper'
        },
        {
            'x': 0.85, 'y': -0.15,
            'text': f"Threshold: {threshold}x",
            'showarrow': False,
            'font': {'size': 10, 'family': 'Arial', 'color': colors['text']},
            'xref': 'paper', 'yref': 'paper'
        },
        {
            'x': 1.0, 'y': -0.25,
            'text': "GeneGone Analysis",
            'showarrow': False,
            'font': {'size': 8, 'family': 'Arial', 'color': 'gray'},
            'xref': 'paper', 'yref': 'paper',
            'align': 'right'
        }
    ]
    
    # Update layout with publication-quality styling
    fig.update_layout(
        width=600,
        height=350,
        margin=dict(t=80, b=80, l=40, r=40),
        paper_bgcolor='white',
        annotations=annotations
    )
    
    return fig

def export_publication_plot(fig: go.Figure, filename: str, format: str = "pdf", width: int = 800, height: int = 600, scale: int = 4) -> str:
    """
    Export a Plotly figure in high-resolution format suitable for publication.
    
    Args:
        fig: Plotly figure object
        filename: Output filename (without extension)
        format: Output format ('pdf', 'svg', 'png', or 'jpg')
        width: Width in pixels
        height: Height in pixels
        scale: Resolution scale factor (higher = better resolution)
    
    Returns:
        Path to the exported file
    """
    # Validate format
    valid_formats = ['pdf', 'svg', 'png', 'jpg']
    if format not in valid_formats:
        print(f"Warning: Invalid format '{format}'. Using 'pdf' instead.")
        format = 'pdf'
    
    # Construct full filename with extension
    full_filename = f"{filename}.{format}"
    
    # Set image export options
    config = {
        'toImageButtonOptions': {
            'format': format,
            'width': width,
            'height': height,
            'scale': scale
        }
    }
    
    # Export with appropriate settings for each format
    if format == 'pdf' or format == 'svg':
        # Vector formats
        fig.write_image(
            full_filename,
            width=width,
            height=height,
            scale=1,  # Scale doesn't apply to vector formats
        )
    else:
        # Raster formats (PNG, JPG)
        fig.write_image(
            full_filename,
            width=width,
            height=height,
            scale=scale,  # Higher scale for better resolution
        )
    
    print(f"Exported high-resolution plot to: {full_filename}")
    return full_filename
