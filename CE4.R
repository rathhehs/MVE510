# Q1
counts_file <- "16s_counts.txt"
annotation_file <- "16s_annotation.txt"

counts <- read.table(counts_file, header = TRUE, sep = "\t", quote = "", comment.char = "")
annotations <- read.table(annotation_file, header = TRUE, sep = "\t", quote = "", comment.char = "")

print("Counts Data:")
print(head(counts))
print("Annotations Data:")
print(head(annotations))

total_counts <- colSums(counts)
print("Total counts per sample:")
print(total_counts)

print("Annotations structure:")
print(str(annotations))

incomplete_annotations <- annotations[apply(annotations, 1, function(row) any(is.na(row))), ]
print("Incomplete annotations:")
print(incomplete_annotations)
print(paste("Number of incomplete annotations:", nrow(incomplete_annotations)))

counts <- counts[rowSums(counts) >= 5, ]
annotations <- annotations[rownames(counts), ]
print(paste("Remaining OTUs after filtering:", nrow(counts)))

# Q2
library(ggplot2)

counts_file <- "16s_counts.txt"
counts <- read.table(counts_file, header = TRUE, sep = "\t", quote = "", comment.char = "")

counts_t <- t(counts)

pca_raw <- prcomp(counts_t, scale. = TRUE)

pca_data_raw <- data.frame(pca_raw$x, Group = c("HC", "HC", "HC", "LC", "LC", "LC"))

ggplot(pca_data_raw, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "PCA", x = "Component1", y = "Component2") +
  theme_minimal()

counts_vst <- log(counts + 1)

counts_vst_t <- t(counts_vst)

pca_vst <- prcomp(counts_vst_t, scale. = TRUE)

pca_data_vst <- data.frame(pca_vst$x, Group = c("HC", "HC", "HC", "LC", "LC", "LC"))

ggplot(pca_data_vst, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "PCA", x = "Component1", y = "Component2") +
  theme_minimal()

# Q3
rarefy_sample <- function(OTUs, counts, depth) {
  reads <- rep(OTUs, times = counts)
  reads_sample <- sample(reads, size = depth, replace = FALSE)
  counts_sample <- as.data.frame(table(reads_sample))
  colnames(counts_sample) <- c("OTU", "Count")
  return(counts_sample)
}

depth <- min(colSums(counts))

rarefied_data_list <- list()
for (i in 1:ncol(counts)) {
  OTUs <- rownames(counts)
  counts_for_OTUs <- counts[, i]
  rarefied_data_list[[i]] <- rarefy_sample(OTUs, counts_for_OTUs, depth)
}

print("Rarefied Data for each sample:")
for (i in 1:length(rarefied_data_list)) {
  print(paste("Sample", i, ":"))
  print(rarefied_data_list[[i]])
}

# Q4
richness <- function(counts_sample) {
  richness_value <- sum(counts_sample$Count > 0)
  return(richness_value)
}

shannon_index <- function(counts_sample) {
  total_count <- sum(counts_sample$Count)
  p_i <- counts_sample$Count / total_count
  H_prime <- -sum(p_i * log(p_i))
  return(H_prime)
}

richness_values <- numeric(length(rarefied_data_list))
shannon_values <- numeric(length(rarefied_data_list))

for (i in 1:length(rarefied_data_list)) {
  counts_sample <- rarefied_data_list[[i]]
  richness_values[i] <- richness(counts_sample)
  shannon_values[i] <- shannon_index(counts_sample)
}

print("Richness Values for each sample:")
print(richness_values)
print("Shannon Index for each sample:")
print(shannon_values)

# Q5
library(DESeq2)

design.matrix <- data.frame(exposure = c(1, 1, 1, 0, 0, 0))
counts.ds <- DESeqDataSetFromMatrix(countData = counts, colData = design.matrix, design = ~exposure)
res.ds <- DESeq(counts.ds)
results_ds <- results(res.ds, independentFiltering = FALSE, cooksCutoff = FALSE)
result_df <- as.data.frame(results_ds)
result_df <- result_df[order(result_df$padj), ]

print("Differentially Abundant OTUs:")
print(result_df)

significant_OTUs <- result_df[result_df$padj < 0.05, ]
print("Significant OTUs:")
print(significant_OTUs)

num_significant_OTUs <- nrow(significant_OTUs)
print(paste("Number of Significant OTUs:", num_significant_OTUs))

# Q6
top10_OTUs <- head(result_df, 10)
print("Top 10 most significant OTUs:")
print(top10_OTUs)

# Q7
counts_file <- "gene_counts.txt"
annotation_file <- "gene_annotation.txt"

gene_counts <- read.table(counts_file, header = TRUE, sep = "\t", quote = "", comment.char = "")
gene_annotations <- read.table(annotation_file, header = TRUE, sep = "\t", quote = "", comment.char = "")

print("Gene Counts (head):")
print(head(gene_counts))

print("Gene Annotations (head):")
print(head(gene_annotations))

total_counts <- colSums(gene_counts)
print("Total reads per sample:")
print(total_counts)

print("Annotations structure:")
print(str(gene_annotations))

filtered_counts <- gene_counts[rowSums(gene_counts) >= 5, ]
filtered_annotations <- gene_annotations[rownames(filtered_counts), ]
print(paste("Remaining genes after filtering:", nrow(filtered_counts)))

library(ggplot2)

counts_t <- t(filtered_counts)
pca_raw <- prcomp(counts_t, scale. = TRUE)
pca_data_raw <- data.frame(pca_raw$x, Group = c("HC", "HC", "HC", "LC", "LC", "LC"))

ggplot(pca_data_raw, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "PCA (Raw Counts)", x = "PC1", y = "PC2") +
  theme_minimal()

counts_vst <- log(filtered_counts + 1)
counts_vst_t <- t(counts_vst)
pca_vst <- prcomp(counts_vst_t, scale. = TRUE)

pca_data_vst <- data.frame(pca_vst$x, Group = c("HC", "HC", "HC", "LC", "LC", "LC"))

ggplot(pca_data_vst, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "PCA (Log-Transformed Counts)", x = "PC1", y = "PC2") +
  theme_minimal()

# Q8
richness <- function(counts_sample) {
  richness_value <- sum(counts_sample$Count > 0)
  return(richness_value)
}

shannon_index <- function(counts_sample) {
  total_count <- sum(counts_sample$Count)
  p_i <- counts_sample$Count / total_count
  p_i <- p_i[p_i > 0]
  H_prime <- -sum(p_i * log(p_i))
  return(H_prime)
}

rarefied_gene_data_list <- list()
for (i in 1:ncol(filtered_counts)) {
  genes <- rownames(filtered_counts)
  counts_for_genes <- filtered_counts[, i]
  depth <- min(colSums(filtered_counts))
  rarefied_gene_data_list[[i]] <- rarefy_sample(genes, counts_for_genes, depth)
}

richness_values <- numeric(length(rarefied_gene_data_list))
shannon_values <- numeric(length(rarefied_gene_data_list))

for (i in 1:length(rarefied_gene_data_list)) {
  counts_sample <- rarefied_gene_data_list[[i]]
  richness_values[i] <- richness(counts_sample)
  shannon_values[i] <- shannon_index(counts_sample)
}

print("Richness Values for each sample:")
print(richness_values)

print("Shannon Index for each sample:")
print(shannon_values)

# Q9
library(DESeq2)

design.matrix <- data.frame(exposure = c(1, 1, 1, 0, 0, 0))
counts.ds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = design.matrix, design = ~exposure)
res.ds <- DESeq(counts.ds)
results_ds <- results(res.ds, independentFiltering = FALSE, cooksCutoff = FALSE)

result_df <- as.data.frame(results_ds)
result_df <- result_df[order(result_df$padj), ]

print("Differentially Abundant Genes:")
print(result_df)

significant_genes <- result_df[result_df$padj < 0.05, ]
print("Significant Genes:")
print(significant_genes)

num_significant_genes <- nrow(significant_genes)
print(paste("Number of Significant Genes:", num_significant_genes))

genes_less_than_zero <- sum(significant_genes$log2FoldChange < 0)
genes_greater_than_zero <- sum(significant_genes$log2FoldChange > 0)
print(paste("Number of genes with log2FoldChange < 0:", genes_less_than_zero))
print(paste("Number of genes with log2FoldChange > 0:", genes_greater_than_zero))

gene_pdxA <- grep("TIGR00557", rownames(significant_genes), ignore.case = TRUE, value = TRUE)
pdxA_result <- result_df[rownames(result_df) %in% gene_pdxA, ]
print("pdxA Gene Result:")
print(pdxA_result)
