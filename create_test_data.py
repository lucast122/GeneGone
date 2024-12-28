#!/usr/bin/env python3
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

def create_reference_genome():
    """Create a simple reference genome with a target gene."""
    # Create a 10kb sequence with a 1kb gene in the middle
    backbone = ''.join(random.choice('ATCG') for _ in range(4000))
    gene = ''.join(random.choice('ATCG') for _ in range(1000))
    after_gene = ''.join(random.choice('ATCG') for _ in range(4000))
    
    full_sequence = backbone + gene + after_gene
    
    record = SeqRecord(
        Seq(full_sequence),
        id="test_chr1",
        description="Test reference genome"
    )
    
    os.makedirs('test_data', exist_ok=True)
    with open('test_data/reference.fasta', 'w') as f:
        SeqIO.write(record, f, 'fasta')
    
    return len(backbone), len(gene)

def create_nanopore_reads(gene_start, gene_length, deleted=False):
    """Create simulated Nanopore reads."""
    os.makedirs('test_data/nanopore', exist_ok=True)
    
    with open('test_data/reference.fasta') as f:
        reference = str(next(SeqIO.parse(f, 'fasta')).seq)
    
    reads = []
    # Create 100 reads
    for i in range(100):
        if deleted and random.random() < 0.8:  # 80% of reads show deletion
            # Create reads that skip the gene
            start = random.randint(0, gene_start - 1000)
            pre_deletion = reference[start:gene_start]
            post_start = gene_start + gene_length
            post_deletion = reference[post_start:post_start + 1000]
            seq = pre_deletion + post_deletion
        else:
            # Normal reads
            start = random.randint(0, len(reference) - 2000)
            seq = reference[start:start + 2000]
        
        # Add random errors (15% error rate)
        seq = list(seq)
        for j in range(len(seq)):
            if random.random() < 0.15:
                if random.random() < 0.6:  # substitution
                    seq[j] = random.choice('ATCG')
                elif random.random() < 0.8:  # deletion
                    seq[j] = ''
                else:  # insertion
                    seq[j] = seq[j] + random.choice('ATCG')
        seq = ''.join(filter(None, seq))
        
        # Create quality scores (Nanopore-like)
        quals = [max(4, min(35, int(random.gauss(15, 5)))) for _ in range(len(seq))]
        
        read = SeqRecord(
            Seq(seq),
            id=f"read_{i}",
            description=f"Simulated Nanopore read {i}"
        )
        read.letter_annotations["phred_quality"] = quals
        reads.append(read)
    
    with open('test_data/nanopore/reads_with_deletion.fastq', 'w') as f:
        SeqIO.write(reads, f, 'fastq')

def create_illumina_reads(gene_start, gene_length, deleted=False):
    """Create simulated Illumina reads."""
    os.makedirs('test_data/illumina', exist_ok=True)
    
    with open('test_data/reference.fasta') as f:
        reference = str(next(SeqIO.parse(f, 'fasta')).seq)
    
    reads = []
    # Create 1000 paired-end reads
    for i in range(1000):
        if deleted and random.random() < 0.8:  # 80% of reads show deletion
            # Create reads that span the deletion junction
            pre_start = random.randint(max(0, gene_start - 200), gene_start - 50)
            post_start = gene_start + gene_length
            seq = reference[pre_start:gene_start] + reference[post_start:post_start + 50]
            
            # Create 150bp paired-end reads
            read1_seq = seq[:150]
            read2_seq = seq[-150:]
        else:
            # Normal reads
            start = random.randint(0, len(reference) - 300)
            seq = reference[start:start + 300]
            read1_seq = seq[:150]
            read2_seq = seq[-150:]
        
        # Add random errors (1% error rate)
        for read_seq in [read1_seq, read2_seq]:
            read_seq = list(read_seq)
            for j in range(len(read_seq)):
                if random.random() < 0.01:
                    read_seq[j] = random.choice('ATCG')
            read_seq = ''.join(read_seq)
            
            # Create quality scores (Illumina-like)
            quals = [max(25, min(40, int(random.gauss(35, 3)))) for _ in range(len(read_seq))]
            
            read = SeqRecord(
                Seq(read_seq),
                id=f"read_{i}",
                description=f"Simulated Illumina read {i}"
            )
            read.letter_annotations["phred_quality"] = quals
            reads.append(read)
    
    with open('test_data/illumina/reads_with_deletion.fastq', 'w') as f:
        SeqIO.write(reads, f, 'fastq')

def main():
    print("Creating test dataset...")
    
    # Create reference genome
    gene_start, gene_length = create_reference_genome()
    print(f"\nReference genome created:")
    print(f"Target gene location: chr1:{gene_start+1}-{gene_start+gene_length}")
    
    # Create reads with deletion
    print("\nCreating Nanopore reads with deletion...")
    create_nanopore_reads(gene_start, gene_length, deleted=True)
    
    print("Creating Illumina reads with deletion...")
    create_illumina_reads(gene_start, gene_length, deleted=True)
    
    print("\nTest data creation complete!")
    print("\nTest files created:")
    print("- test_data/reference.fasta")
    print("- test_data/nanopore/reads_with_deletion.fastq")
    print("- test_data/illumina/reads_with_deletion.fastq")
    print(f"\nUse these coordinates in the app:")
    print(f"Chromosome: test_chr1")
    print(f"Start: {gene_start+1}")
    print(f"End: {gene_start+gene_length}")

if __name__ == "__main__":
    main()
