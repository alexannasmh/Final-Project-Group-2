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
```
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

# Rename it to reference.fasta

```
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa human.fasta
```
# Download annotation file
```
# Download the GTF file
```
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
```

# Unzip the GTF file
```
gunzip Homo_sapiens.GRCh38.113.gtf.gz
```

# Rename it to annotation.gtf
```
mv Homo_sapiens.GRCh38.113.gtf human.gtf
```
# STAR MAPPING
- Create logs directory
```
mkdir logs
```
- Create index script
```
vi index.slurm
```
- Type I to edit and paste:
```
#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=logs/star_index_%j.out
#SBATCH --mail-user=jparedes@svsu.edu
#SBATCH --mail-type=END,FAIL

# Load STAR
module load STAR

# Paths
GENOME_DIR="star_index"
FASTA="human.fa"
GTF="human.gtf"

mkdir -p "$GENOME_DIR"

# Run STAR genomeGenerate
STAR --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$FASTA" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang 100 \
     --runThreadN 32
  
```

- Create a script called `star.slurm`
```
vi star.slurm
```
- Type I to edit
```
#!/bin/bash
#SBATCH --job-name=star_map
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/star_map_%A_%a.out
#SBATCH --array=0-7
#SBATCH --mail-user=jparedes@svsu.edu
#SBATCH --mail-type=END,FAIL

# Load STAR module
module load STAR

# Paths
GENOME_DIR="star_index"
FASTQ_DIR="fastq_files"
OUTPUT_DIR="star_output"
THREADS=8

mkdir -p "$OUTPUT_DIR" logs

# Get sample names from _1.fastq files
SAMPLES=($(ls ${FASTQ_DIR}/*_1.fastq | sed 's|.*/||' | sed 's/_1.fastq//' | sort))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Define input FASTQ files
R1="${FASTQ_DIR}/${SAMPLE}_1.fastq"
R2="${FASTQ_DIR}/${SAMPLE}_2.fastq"

echo "[$(date)] Starting STAR alignment for sample: $SAMPLE"
```
- Run script
```
sbatch star.slurm
```

# Check if you have subread environment
```
module load anaconda3
conda env list
```

- If env not installed:
```
conda create -n subreads_env -c bioconda subread
```
- activate environment
```
conda activate subreads_env
```

# Files will be in

```
featureCounts -T 8 -a human.gtf -o counts.txt \
star_output/SRR16101430_Aligned.sortedByCoord.out.bam \
star_output/SRR16101431_Aligned.sortedByCoord.out.bam \
star_output/SRR16101432_Aligned.sortedByCoord.out.bam \
star_output/SRR16101433_Aligned.sortedByCoord.out.bam \
star_output/SRR16101435_Aligned.sortedByCoord.out.bam \
star_output/SRR16101436_Aligned.sortedByCoord.out.bam \
star_output/SRR16101437_Aligned.sortedByCoord.out.bam \
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
# Optional: set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Then reinstall DESeq2
BiocManager::install("DESeq2", force = TRUE)

install.packages("ggrepel")

library(DESeq2)
library(ggplot2)
library(ggrepel)

# Load count matrix
counts <- read.table("counts.txt", header = TRUE, row.names = 1, comment.char = "#", check.names = FALSE)
counts <- counts[, 6:ncol(counts)]  # Keep only sample columns

# Sample metadata
condition <- factor(c(rep("Healthy", 4), rep("Alzheimer", 4)))
coldata <- data.frame(row.names = colnames(counts), condition = condition)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Add significance labels
res$padj[is.na(res$padj)] <- 1  # Replace NAs
res_df <- as.data.frame(res)
res_df$sig <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Upregulated",
                     ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Downregulated", "Not Significant"))

# Label all significant genes
res_df$label <- ifelse(res_df$sig != "Not Significant", rownames(res_df), NA)

# Final volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6) +
  geom_text_repel(aes(label = label, color = sig), size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  xlim(-10, 10) +
  labs(title = "Differential Gene Expression in Alzheimer’s Disease",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()


ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = label, color = sig), size = 4, max.overlaps = Inf) +  # Bigger gene labels
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  xlim(-10, 10) +
  labs(title = "Differential Gene Expression in Alzheimer’s Disease",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", color = "Regulation") +
  theme_minimal(base_size = 16) +  # Bumps up all text sizes
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Filter to top N genes (or set your own logic)
label_df <- subset(res_df, !is.na(label) & label != "")

# Separate up/downregulated
label_up <- label_df[label_df$sig == "Upregulated", ]
label_down <- label_df[label_df$sig == "Downregulated", ]

# Sort by padj or any metric you want
label_up <- label_up[order(label_up$padj), ]
label_down <- label_down[order(label_down$padj), ]

# Evenly spaced y-axis label positions (custom)
label_up$label_y <- seq(from = max(-log10(res_df$padj), na.rm = TRUE),
                        to = max(-log10(res_df$padj), na.rm = TRUE) - 4,
                        length.out = nrow(label_up))

label_down$label_y <- seq(from = max(-log10(res_df$padj), na.rm = TRUE) - 4,
                          to = max(-log10(res_df$padj), na.rm = TRUE) - 6,
                          length.out = nrow(label_down))

# Set label x positions and hjust
label_up$label_x <- 7.5
label_up$hjust <- 0
label_up$color <- "red"

label_down$label_x <- -7.5
label_down$hjust <- 1
label_down$color <- "blue"

# Combine
label_data <- rbind(label_up, label_down)

# Base volcano
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6, size = 2) +
  
  # Line from point to label
  geom_segment(data = label_data,
               aes(x = log2FoldChange, xend = label_x,
                   y = -log10(padj), yend = label_y,
                   color = sig),
               size = 0.3, show.legend = FALSE) +
  
  # Label itself
  geom_text(data = label_data,
            aes(x = label_x, y = label_y, label = label, hjust = hjust),
            color = label_data$color,
            size = 4,
            show.legend = FALSE) +
  
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "gray")) +
  xlim(-12, 12) +
  labs(
    title = "Differential Gene Expression in Alzheimer’s Disease",
    x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", color = "Regulation"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

```

