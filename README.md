# Final-Project-Group-2

Our project focuses on how gene expression differs in individuals with Alzheimer’s disease (AD) compared to healthy controls. To investigate this, we will use RNA-seq data from blood samples and apply tools such as HISAT2 for read alignment and DESeq2 for analyzing differential gene expression. Our approach centers on identifying genes that are significantly upregulated (overexpressed) or downregulated (underexpressed) in AD samples. We will then perform functional annotation and pathway analysis using EnrichR to interpret the biological significance of these changes. Through transcriptomic analysis, we aim to uncover key gene expression patterns associated with AD and gain insight into the molecular mechanisms driving disease progression.

## Downloading sequences
module load sra-toolkit


## Download reference genome
# Download the fasta file
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Rename it to reference.fasta
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa reference.fasta
