#!/usr/bin/env python3
import os
import subprocess
import gzip
import shutil
from pathlib import Path

def download_file(url, output_path):
    """Download a file using wget."""
    subprocess.run(['wget', '-O', output_path, url], check=True)

def gunzip_file(gz_path, output_path):
    """Decompress a gzipped file."""
    with gzip.open(gz_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def main():
    # Create test data directories
    for platform in ['nanopore', 'pacbio', 'illumina']:
        os.makedirs(f'test_data/{platform}', exist_ok=True)

    # Download E. coli reference genome
    print("Downloading E. coli reference genome...")
    ref_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    ref_gz = "test_data/ecoli_ref.fna.gz"
    ref_fna = "test_data/ecoli_ref.fna"
    
    download_file(ref_url, ref_gz)
    gunzip_file(ref_gz, ref_fna)
    os.remove(ref_gz)

    # Download test data for each platform
    # Nanopore
    print("Downloading Nanopore test data...")
    nano_url = "https://sra-pub-src-1.s3.amazonaws.com/ERR3152364/16s_minion.fastq.gz"
    nano_gz = "test_data/nanopore/reads.fastq.gz"
    nano_fq = "test_data/nanopore/reads.fastq"
    
    download_file(nano_url, nano_gz)
    gunzip_file(nano_gz, nano_fq)
    os.remove(nano_gz)

    # PacBio
    print("Downloading PacBio test data...")
    pacbio_url = "https://sra-pub-src-1.s3.amazonaws.com/ERR3237140/pacbio.fastq.gz"
    pacbio_gz = "test_data/pacbio/reads.fastq.gz"
    pacbio_fq = "test_data/pacbio/reads.fastq"
    
    download_file(pacbio_url, pacbio_gz)
    gunzip_file(pacbio_gz, pacbio_fq)
    os.remove(pacbio_gz)

    # Illumina
    print("Downloading Illumina test data...")
    illumina_url = "https://sra-pub-src-1.s3.amazonaws.com/ERR022075/ERR022075_1.fastq.gz"
    illumina_gz = "test_data/illumina/reads.fastq.gz"
    illumina_fq = "test_data/illumina/reads.fastq"
    
    download_file(illumina_url, illumina_gz)
    gunzip_file(illumina_gz, illumina_fq)
    os.remove(illumina_gz)

    print("\nTest data preparation complete!")
    print("\nTest Gene Coordinates:")
    print("Chromosome: NC_000913.3")
    print("lacZ gene location:")
    print("Start: 365,438")
    print("End: 368,512")

if __name__ == "__main__":
    main()
