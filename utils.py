import subprocess
import os
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
            os.rename("temp_filtered.bam", output_bam)
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

def calculate_coverage(bam_path: str, region: str, read_type: ReadType, min_mapping_quality: int = None) -> pd.DataFrame:
    """
    Calculate coverage for a specific genomic region with read-type specific parameters.
    
    Args:
        bam_path: Path to BAM file
        region: Region string (e.g., "chr1:1000-2000")
        read_type: Type of sequencing reads
        min_mapping_quality: Minimum mapping quality score (0-60). If None, uses platform-specific defaults:
            - Illumina: 20 (higher confidence due to shorter reads)
            - Nanopore/PacBio: 0 (length provides confidence despite errors)
    
    Returns:
        DataFrame with position and coverage
    """
    try:
        # Set default mapping quality based on read type if not provided
        if min_mapping_quality is None:
            min_mapping_quality = 20 if read_type == ReadType.ILLUMINA else 0
        
        # Convert mapping quality to string for samtools
        min_qual = str(min_mapping_quality)
        
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

    # Define flanking region size (same as gene size, but max 1kb)
    gene_size = gene_end - gene_start
    flank_size = min(gene_size, 1000)  # Use up to 1kb flanking regions
    
    # Get flanking regions, ensuring we don't go out of bounds
    min_pos = min(coverage_df['position'])
    max_pos = max(coverage_df['position'])
    
    left_start = max(min_pos, gene_start - flank_size)
    right_end = min(max_pos, gene_end + flank_size)
    
    # Get flanking regions
    left_flank = coverage_df[
        (coverage_df['position'] >= left_start) & 
        (coverage_df['position'] < gene_start)
    ]
    right_flank = coverage_df[
        (coverage_df['position'] > gene_end) & 
        (coverage_df['position'] <= right_end)
    ]
    gene_region = coverage_df[
        (coverage_df['position'] >= gene_start) & 
        (coverage_df['position'] <= gene_end)
    ]
    
    # Calculate statistics with safety checks
    gene_mean = gene_region['coverage'].mean() if not gene_region.empty else 0
    left_mean = left_flank['coverage'].mean() if not left_flank.empty else 0
    right_mean = right_flank['coverage'].mean() if not right_flank.empty else 0
    
    # Print debug information
    print(f"Debug - Gene region: {gene_start}-{gene_end}, Mean coverage: {gene_mean}")
    print(f"Debug - Left flank: {left_start}-{gene_start}, Mean coverage: {left_mean}")
    print(f"Debug - Right flank: {gene_end}-{right_end}, Mean coverage: {right_mean}")
    
    # Use available flanking regions for mean calculation
    if left_flank.empty and right_flank.empty:
        flank_mean = gene_mean  # If no flanking regions, use gene coverage
    elif left_flank.empty:
        flank_mean = right_mean
    elif right_flank.empty:
        flank_mean = left_mean
    else:
        flank_mean = (left_mean + right_mean) / 2
    
    # Ensure we have a valid flanking mean
    if flank_mean == 0:
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
    Create an enhanced coverage plot with gene annotation and deletion visualization.
    """
    # Calculate plot boundaries to show context
    gene_size = gene_end - gene_start
    context_size = max(gene_size * 2, 2000)  # Show at least 2kb or 2x gene size
    plot_start = max(min(coverage_df['position']) - context_size, 0)
    plot_end = max(coverage_df['position']) + context_size
    
    # Create figure
    fig = go.Figure()
    
    # Add coverage trace
    fig.add_trace(go.Scatter(
        x=coverage_df['position'],
        y=coverage_df['coverage'],
        name='Coverage',
        line=dict(
            color='blue' if read_type == ReadType.ILLUMINA else 'green',
            width=2
        ),
        fill='tozeroy',
        fillcolor='rgba(0, 0, 255, 0.1)'
    ))
    
    # Calculate y position for gene annotation
    max_coverage = max(coverage_df['coverage'])
    gene_y = max_coverage * 1.1
    
    # Add gene box annotation
    fig.add_shape(
        type="rect",
        x0=gene_start,
        x1=gene_end,
        y0=gene_y - max_coverage * 0.02,  # Small height for the gene box
        y1=gene_y + max_coverage * 0.02,
        line=dict(color="red", width=2),
        fillcolor="red",
    )
    
    # Add directional triangles to show gene orientation
    triangle_size = min((gene_end - gene_start) * 0.1, 50)  # Scale triangle with gene size
    triangle_spacing = (gene_end - gene_start) / 8
    
    for x_pos in [gene_start + triangle_spacing * i for i in range(1, 8, 2)]:
        fig.add_shape(
            type="path",
            path=f"M {x_pos} {gene_y} L {x_pos + triangle_size} {gene_y} L {x_pos + triangle_size/2} {gene_y + max_coverage * 0.04} Z",
            fillcolor="red",
            line=dict(color="red", width=1),
        )
    
    # Add gene label with background
    fig.add_annotation(
        x=(gene_start + gene_end) / 2,
        y=gene_y + max_coverage * 0.06,
        text=f"<b>{gene_name}</b>",
        showarrow=False,
        font=dict(size=14, color="black"),
        bgcolor="white",
        bordercolor="red",
        borderwidth=2,
        borderpad=4
    )
    
    if analysis_results:
        # Add threshold line
        fig.add_hline(
            y=analysis_results['threshold_used'],
            line_dash="dash",
            line_color="red",
            annotation=dict(
                text=f"Coverage Threshold ({analysis_results['threshold_used']:.1f}x)",
                font=dict(color="red"),
                bgcolor="white",
                bordercolor="red",
                borderwidth=1,
                borderpad=4
            ),
            annotation_position="right"
        )
        
        # Add mean coverage lines for context
        fig.add_hline(
            y=analysis_results['flank_mean_coverage'],
            line_dash="dash",
            line_color="green",
            annotation=dict(
                text=f"Mean Flanking Coverage ({analysis_results['flank_mean_coverage']:.1f}x)",
                font=dict(color="green"),
                bgcolor="white",
                bordercolor="green",
                borderwidth=1,
                borderpad=4
            ),
            annotation_position="right"
        )
        
        # Add deletion indicator if detected
        if analysis_results['deletion_detected']:
            # Add semi-transparent red overlay
            fig.add_shape(
                type="rect",
                x0=gene_start,
                x1=gene_end,
                y0=0,
                y1=max_coverage,
                fillcolor="rgba(255, 0, 0, 0.1)",
                line=dict(color="red", dash="dash", width=2),
                layer="below"
            )
            
            # Add deletion text with background
            fig.add_annotation(
                x=(gene_start + gene_end) / 2,
                y=max_coverage * 0.5,
                text=f"<b>DELETION DETECTED</b><br>{analysis_results['deletion_percentage']:.1f}% reduction",
                showarrow=False,
                font=dict(size=20, color="red"),
                bgcolor="rgba(255, 255, 255, 0.9)",
                bordercolor="red",
                borderwidth=2,
                borderpad=4,
                align="center"
            )
    
    # Update layout with improved styling
    fig.update_layout(
        title=dict(
            text=f"Gene Deletion Analysis ({read_type.value.capitalize()} Data)",
            x=0.5,
            font=dict(size=20)
        ),
        xaxis=dict(
            title=dict(
                text="Genomic Position",
                font=dict(size=14)
            ),
            range=[plot_start, plot_end],
            showgrid=True,
            gridcolor='rgba(0,0,0,0.1)',
            zeroline=False
        ),
        yaxis=dict(
            title=dict(
                text="Coverage Depth",
                font=dict(size=14)
            ),
            rangemode="tozero",
            showgrid=True,
            gridcolor='rgba(0,0,0,0.1)',
            zeroline=False
        ),
        showlegend=True,
        hovermode='x unified',
        height=600,
        template="plotly_white",
        plot_bgcolor='white',
        margin=dict(t=100, r=200)  # Extra right margin for annotations
    )
    
    return fig

def create_deletion_summary(coverage_df: pd.DataFrame, gene_start: int, gene_end: int,
                          threshold: int, read_type: ReadType) -> go.Figure:
    """
    Create a summary visualization of the deletion analysis.
    
    Args:
        coverage_df: DataFrame with coverage data
        gene_start: Start position of the gene
        gene_end: End position of the gene
        threshold: Coverage threshold for deletion detection
        read_type: Type of sequencing reads
    
    Returns:
        Plotly figure object
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
    
    # Create gauge chart
    fig = go.Figure(go.Indicator(
        mode = "gauge+number+delta",
        value = deletion_percentage,
        delta = {'reference': 80},
        title = {'text': "Deletion Confidence", 'font': {'size': 24}},
        gauge = {
            'axis': {'range': [None, 100], 'tickwidth': 1},
            'bar': {'color': "red"},
            'bgcolor': "white",
            'borderwidth': 2,
            'bordercolor': "gray",
            'steps': [
                {'range': [0, 40], 'color': 'lightgray'},
                {'range': [40, 80], 'color': 'orange'},
                {'range': [80, 100], 'color': 'lightgreen'}
            ],
            'threshold': {
                'line': {'color': "black", 'width': 4},
                'thickness': 0.75,
                'value': 80
            }
        }
    ))
    
    # Add metrics as annotations
    fig.add_annotation(
        x=0.1, y=-0.2,
        text=f"Mean Coverage: {mean_coverage:.1f}",
        showarrow=False,
        font=dict(size=14)
    )
    fig.add_annotation(
        x=0.5, y=-0.2,
        text=f"Max Coverage: {max_coverage:.1f}",
        showarrow=False,
        font=dict(size=14)
    )
    fig.add_annotation(
        x=0.9, y=-0.2,
        text=f"Threshold: {threshold}",
        showarrow=False,
        font=dict(size=14)
    )
    
    # Update layout
    fig.update_layout(
        height=400,
        margin=dict(t=100, b=100),
        template="plotly_white"
    )
    
    return fig
