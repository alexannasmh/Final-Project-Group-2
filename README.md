# Final-Project-Group-2

Our project focuses on how gene expression differs in individuals with Alzheimer’s disease (AD) compared to healthy controls. To investigate this, we will use RNA-seq data from blood samples and apply tools such as HISAT2 for read alignment and DESeq2 for analyzing differential gene expression. Our approach centers on identifying genes that are significantly upregulated (overexpressed) or downregulated (underexpressed) in AD samples. We will then perform functional annotation and pathway analysis using EnrichR to interpret the biological significance of these changes. Through transcriptomic analysis, we aim to uncover key gene expression patterns associated with AD and gain insight into the molecular mechanisms driving disease progression.

## Downloading sequences
module load sra-toolkit


## Download reference genome
# Download the fasta file
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
- DONE!!! WOOHOO!

# Unzip the file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Rename it to reference.fasta
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa reference.fasta

# STAR MAPPING
- Create a script called `star.bh`
```
vi star.bh
```
- Type I to edit
```
#!/bin/bash

# Set paths
GENOME_DIR="star_index"
FASTQ_DIR="fastq_files"
OUTPUT_DIR="star_output"
THREADS=8

# Make output directory if it doesn't exist
module load STAR
mkdir -p "$OUTPUT_DIR"

# Loop over all *_1.fastq files to get sample names
for R1 in "$FASTQ_DIR"/*_1.fastq; do
    SAMPLE=$(basename "$R1" _1.fastq)
    R2="$FASTQ_DIR/${SAMPLE}_2.fastq"

    # Run STAR
    STAR --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1" "$R2" \
         --outFileNamePrefix "${OUTPUT_DIR}/${SAMPLE}_" \
         --runThreadN "$THREADS" \
         --outSAMtype BAM SortedByCoordinate
done
```

