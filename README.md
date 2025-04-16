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
```
module load sra-toolkit
```
```
#!/bin/bash
module load anaconda3
conda activate sra-env

while read -r SRR; do
  echo "Downloading $SRR..."
  prefetch --max-size 100G $SRR
  fasterq-dump $SRR --split-files -O fastq_files/
done < srr_accessions.txt
```

## Download reference genome


# Download the fasta file
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


# Unzip the file
'''
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
'''

# Rename it to reference.fasta

'''
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa reference.fasta
'''

# STAR MAPPING
'''
STAR --runMode genomeGenerate \
  --genomeDir star_index \
  --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf \
  --sjdbOverhang 100
'''

- Create a script called `star.bh`
```
vi star.slurm
```
- Type I to edit
```
#!/bin/bash
#SBATCH --job-name=star_alzheim
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=star_%A_%a.out
#SBATCH --array=0-7
#SBATCH --mail-user=aetaylo2@svsu.edu
#SBATCH --mail-type=ALL

# Load STAR module
module load STAR

# Set paths
GENOME_DIR="star_index"
FASTQ_DIR="fastq_files"
OUTPUT_DIR="star_output"
THREADS=8

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Get all sample names
SAMPLES=($(ls ${FASTQ_DIR}/*_1.fastq | sed 's|.*/||' | sed 's/_1.fastq//' | sort))

# Get current sample for this task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R1="${FASTQ_DIR}/${SAMPLE}_1.fastq"
R2="${FASTQ_DIR}/${SAMPLE}_2.fastq"

echo "Processing sample: $SAMPLE"

# Run STAR
STAR --genomeDir "$GENOME_DIR" \
     --readFilesIn "$R1" "$R2" \
     --outFileNamePrefix "${OUTPUT_DIR}/${SAMPLE}_" \
     --runThreadN "$THREADS" \
     --outSAMtype BAM SortedByCoordinate
```
# Download annotation file
```
# Download the GTF file
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

# Unzip the GTF file
gunzip Homo_sapiens.GRCh38.113.gtf.gz

# Rename it to annotation.gtf
mv Homo_sapiens.GRCh38.113.gtf annotation.gtf
```
# Check if you have subread environment
```
module load anaconda3
conda env list
```
# Files will be in

```
featureCounts -T 8 -a annotation.gtf -o counts.txt \
star_output/SRR16101430_Aligned.sortedByCoord.out.bam
star_output/SRR16101431_Aligned.sortedByCoord.out.bam
star_output/SRR16101432_Aligned.sortedByCoord.out.bam
star_output/SRR16101433_Aligned.sortedByCoord.out.bam
star_output/SRR16101435_Aligned.sortedByCoord.out.bam
star_output/SRR16101436_Aligned.sortedByCoord.out.bam
star_output/SRR16101437_Aligned.sortedByCoord.out.bam
star_output/SRR16101438_Aligned.sortedByCoord.out.bam


   
