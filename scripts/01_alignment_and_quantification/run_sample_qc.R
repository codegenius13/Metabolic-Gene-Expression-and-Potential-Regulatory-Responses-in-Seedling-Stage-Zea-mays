library(DESeq2)

count_file <- "/mnt/c/Users/Boluwatife Ikuerowo/Desktop/work/programming/library/content/src/assets/exerciseOne/04_alignment/counts/maize_gene_counts_fragments.txt"
qc_dir <- "/mnt/c/Users/Boluwatife Ikuerowo/Desktop/work/programming/library/content/src/assets/exerciseOne/04_alignment/qc"

x <- read.delim(count_file, comment.char = "#", check.names = FALSE)
counts <- x[, 7:ncol(x), drop = FALSE]
rownames(counts) <- x$Geneid
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

# Convert BAM-path column names to accession IDs
colnames(counts) <- sub("^.*/", "", colnames(counts))
colnames(counts) <- sub("\\.sorted\\.bam$", "", colnames(counts))

condition_map <- c(
  SRR3984708 = "Control",
  SRR3984749 = "Control",
  SRR3984762 = "Salinity",
  SRR3984771 = "Salinity",
  SRR3984782 = "Drought",
  SRR3984791 = "Drought",
  SRR3984794 = "Heat",
  SRR3984795 = "Heat",
  SRR3984802 = "Cold",
  SRR3984815 = "Cold"
)

condition <- factor(
  unname(condition_map[colnames(counts)]),
  levels = c("Control", "Cold", "Drought", "Heat", "Salinity")
)

if (anyNA(condition)) stop("Missing condition assignment for a sample ID.")

metadata <- data.frame(condition = condition, row.names = colnames(counts))

write.table(
  data.frame(sample = rownames(metadata), condition = metadata$condition),
  file.path(qc_dir, "sample_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

counts <- counts[rowSums(counts) > 0, , drop = FALSE]

write.table(
  data.frame(Sample = colnames(counts),
             Condition = metadata[colnames(counts), "condition"],
             Total_Counts = colSums(counts)),
  file.path(qc_dir, "library_sizes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)
dds <- estimateSizeFactors(dds)

write.table(
  data.frame(Sample = colnames(dds),
             Condition = metadata[colnames(dds), "condition"],
             SizeFactor = sizeFactors(dds)),
  file.path(qc_dir, "size_factors.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

vsd <- vst(dds, blind = TRUE)
vst_matrix <- assay(vsd)

cor_matrix <- cor(vst_matrix, method = "pearson")
write.table(cor_matrix, file.path(qc_dir, "sample_correlation_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

png(file.path(qc_dir, "sample_correlation_heatmap.png"),
    width = 1800, height = 1600, res = 200)
heatmap(cor_matrix, symm = TRUE, margins = c(10, 10),
        main = "Sample Pearson Correlation")
dev.off()

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
write.table(pca_data, file.path(qc_dir, "PCA_coordinates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(qc_dir, "PCA_samples.png"),
    width = 1800, height = 1600, res = 200)
plot(pca_data$PC1, pca_data$PC2, pch = 19,
     xlab = paste0("PC1 (", percent_var[1], "%)"),
     ylab = paste0("PC2 (", percent_var[2], "%)"),
     main = "PCA of Maize RNA-seq Samples")
text(pca_data$PC1, pca_data$PC2,
     labels = rownames(pca_data), pos = 3, cex = 0.7)
dev.off()

sample_dist_matrix <- as.matrix(dist(t(vst_matrix)))
write.table(sample_dist_matrix,
            file.path(qc_dir, "sample_distance_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

png(file.path(qc_dir, "sample_distance_heatmap.png"),
    width = 1800, height = 1600, res = 200)
heatmap(sample_dist_matrix, symm = TRUE, margins = c(10, 10),
        main = "Sample-to-Sample Distance")
dev.off()

cat("QC completed using corrected fragment counts.\n")
