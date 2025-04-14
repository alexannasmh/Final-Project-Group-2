# Identifying Gene Expressions in Alzheimer's Disease
Problem Statement
# Biological Question
What genes are differently expressed in Alzheimer's disease (AD) compared to healthy individuals, and how can bioinformatics tools help identify different gene expressions? 
# Hypothesis
There is no significant difference in gene expression between AD and healthy controls.
# Significance
Alzheimer's is a leading neurodegenerative disease with no cure. Identifying different gene expressions can inprove early diagnosis and therapeutic targets, as well as understanding them could reveal mechanisms of disease progression.
# Goal
Determine if there are similar differences in gene expressions between AD and control samples. 

SRR numbers - 8 sample size for project
(4 AD samples & 4 healthy samples)

# Tools and Software 
    #FastQC: Performs quality control checks on raw sequencing data, generating reports about read quality, GC content, adapter content, etc.

# NCBI Database: Public database where our SRR accessions were found. https://www.ncbi.nlm.nih.gov/sra 

## Downloading sequences
module load sra-toolkit

#!/bin/bash
module load anaconda3
conda activate sra-env

while read -r SRR; do
  echo "Downloading $SRR..."
  prefetch --max-size 100G $SRR
  fasterq-dump $SRR --split-files -O fastq_files/
done < srr_accessions.txt

## Download reference genome


# Download the fasta file
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Rename it to reference.fasta
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa reference.fasta

# STAR MAPPING
STAR --runMode genomeGenerate \
  --genomeDir star_index \
  --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf \
  --sjdbOverhang 100

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
# Check if you have subread environment
```
module load anaconda3
conda env list
```
# Files will be in

    star_output/SRRACCESSIONHEALTY_Aligned.sortedByCoord.out.bam
    star_output/SRRACCESSIONHEALTY_Aligned.sortedByCoord.out.bam
    star_output/SRRACCESSIONHEALTY_Aligned.sortedByCoord.out.bam
    star_output/SRRACCESSIONHEALTY_Aligned.sortedByCoord.out.bam
    star_output/SRRACCESSION_AD_Aligned.sortedByCoord.out.bam
    star_output/SRRACCESSION_AD_Aligned.sortedByCoord.out.bam
    star_output/SRRACCESSION_AD_Aligned.sortedByCoord.out.bam
    star_output/SRRACCESSION_AD_Aligned.sortedByCoord.out.bam

   
