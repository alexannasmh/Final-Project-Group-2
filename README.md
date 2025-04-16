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
- Run script
```
sbatch star.slurm
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
```
# Push the "counts.txt" to your repository
```
git add counts.txt
git commit -m "adding counts file"
git push origin main
```

# Use "counts.txt" file to make volcano plot
- This is the RStudio script to run to make your volcano plot
```
# Load necessary libraries
if (!requireNamespace("DESeq2", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(DESeq2)
library(ggplot2)

# Load counts
counts <- read.delim("counts.txt", comment.char = "#", stringsAsFactors = FALSE)
rownames(counts) <- counts$Geneid
counts <- counts[, grep("^SRR", colnames(counts))]

# Sample metadata
condition <- factor(c(rep("Healthy", 4), rep("Alzheimer", 4)))
coldata <- data.frame(row.names = colnames(counts), condition)

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)

# Add column for significance
res$padj[is.na(res$padj)] <- 1
res$sig <- ifelse(res$padj < 0.05 & res$log2FoldChange > 1, "Up",
           ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "Down", "NotSig"))


# Volcano plot
ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "gray")) +
  labs(title = "Volcano Plot: Alzheimer vs Healthy",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

```

/opt/packages/STAR/2.7.10b/STAR: line 26: 93377 Killed                  singularity exec $OPTIONS $DIRECTORY/singularity-$PACKAGE-$VERSION.sif /opt/STAR_2.7.10b/Linux_x86_64_static/$TOOL "$@"
slurmstepd: error: Detected 1 oom-kill event(s) in StepId=30472752.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.

Processing sample: SRR16101430
	/opt/STAR_2.7.10b/Linux_x86_64_static/STAR --genomeDir star_index --readFilesIn fastq_files/SRR16101430_1.fastq fastq_files/SRR16101430_2.fastq --outFileNamePrefix star_output/SRR16101430_ --runThreadN 8 --outSAMtype BAM SortedByCoordinate
	STAR version: 2.7.10b   compiled: 2022-11-01T09:53:26-04:00 :/home/dobin/data/STAR/STARcode/STAR.master/source
Apr 16 11:42:14 ..... started STAR run
Apr 16 11:42:14 ..... loading genome
/opt/packages/STAR/2.7.10b/STAR: line 26: 63661 Killed                  singularity exec $OPTIONS $DIRECTORY/singularity-$PACKAGE-$VERSION.sif /opt/STAR_2.7.10b/Linux_x86_64_static/$TOOL "$@"
slurmstepd: error: Detected 1 oom-kill event(s) in StepId=30472675.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.

   
