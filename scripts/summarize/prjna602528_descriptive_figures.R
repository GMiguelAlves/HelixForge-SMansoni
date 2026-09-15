args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: prjna602528_descriptive_figures.R TPM METADATA QC OUTDIR")

tpm_path <- args[[1]]
metadata_path <- args[[2]]
qc_path <- args[[3]]
outdir <- args[[4]]

tpm <- read.delim(tpm_path, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)
qc <- read.delim(qc_path, check.names = FALSE)

colnames(tpm) <- sub("^PRJNA602528__", "", colnames(tpm))
metadata <- metadata[order(as.numeric(metadata$time_hours)), ]
sample_order <- metadata$sample_id
tpm <- tpm[, sample_order, drop = FALSE]
log_tpm <- log2(tpm + 1)

variances <- apply(log_tpm, 1, var)
keep <- names(sort(variances, decreasing = TRUE))[
  seq_len(min(2000, sum(is.finite(variances) & variances > 0)))
]
pca <- prcomp(t(log_tpm[keep, , drop = FALSE]), center = TRUE, scale. = TRUE)
pca_table <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
write.table(
  pca_table,
  file.path(outdir, "descriptive_pca_coordinates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

png(file.path(outdir, "descriptive_pca.png"), width = 1400, height = 1000, res = 160)
plot(
  pca$x[, 1], pca$x[, 2], pch = 19, col = "#2457A7",
  xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
  ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2]),
  main = "PRJNA602528 — descriptive PCA"
)
text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.75)
mtext("Descriptive only; no differential-expression inference", side = 1, line = 4, cex = 0.8)
dev.off()

correlation <- cor(log_tpm, method = "pearson")
write.table(
  correlation,
  file.path(outdir, "sample_correlation.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)
png(file.path(outdir, "sample_correlation.png"), width = 1400, height = 1200, res = 160)
heatmap(
  correlation, symm = TRUE,
  col = colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))(100),
  margins = c(10, 10), main = "PRJNA602528 — log2(TPM+1) correlation"
)
mtext("Descriptive only", side = 1, line = 8, cex = 0.8)
dev.off()

top_genes <- names(sort(variances, decreasing = TRUE))[seq_len(min(12, length(variances)))]
trajectory <- t(scale(t(log_tpm[top_genes, , drop = FALSE])))
png(file.path(outdir, "descriptive_trajectories.png"), width = 1500, height = 1000, res = 160)
matplot(
  as.numeric(metadata$time_hours), t(trajectory), type = "l", lty = 1, lwd = 2,
  col = rainbow(length(top_genes)), xlab = "Hours after praziquantel",
  ylab = "Gene-wise z-score of log2(TPM+1)",
  main = "Highly variable descriptive trajectories"
)
legend(
  "topright", legend = top_genes, col = rainbow(length(top_genes)),
  lty = 1, cex = 0.65, ncol = 2
)
mtext("Descriptive only; one resolved library per time point", side = 1, line = 4, cex = 0.8)
dev.off()

qc <- qc[match(sample_order, qc$sample_id), ]
png(file.path(outdir, "mapping_and_retention.png"), width = 1500, height = 1000, res = 160)
values <- rbind(qc$trim_retention_percent, qc$salmon_mapping_percent)
barplot(
  values, beside = TRUE, names.arg = sample_order, las = 2, ylim = c(0, 105),
  col = c("#4DAF4A", "#377EB8"), ylab = "Percent",
  main = "Read retention and Salmon mapping"
)
legend(
  "bottomleft", legend = c("Trim retention", "Salmon mapping"),
  fill = c("#4DAF4A", "#377EB8"), bty = "n"
)
dev.off()
