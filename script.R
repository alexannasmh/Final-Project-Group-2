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

