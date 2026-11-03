#!/usr/bin/env Rscript
# GSE245451 donor-aware differential expression starter analysis
# Run from the project root: Rscript scripts/01_deseq2_discovery.R
suppressPackageStartupMessages(library(DESeq2))

project <- normalizePath(getwd())
counts_path <- file.path(project, "data/processed/GSE245451/GSE245451_counts_matrix.tsv")
meta_path <- file.path(project, "metadata/GSE245451_sample_metadata.tsv")
outdir <- file.path(project, "results/differential_expression")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

count_df <- read.delim(counts_path, check.names = FALSE, row.names = 1)
meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$sample_accession
meta <- meta[colnames(count_df), , drop = FALSE]
stopifnot(identical(rownames(meta), colnames(count_df)))

# Primary paired contrast: lesion vs perilesional skin, blocked by donor.
paired <- meta[meta$group %in% c("lesion", "perilesion"), , drop = FALSE]
paired$group <- factor(paired$group, levels = c("perilesion", "lesion"))
paired$donor_id <- factor(paired$donor_id)
paired_counts <- round(as.matrix(count_df[, rownames(paired), drop = FALSE]))
dds_pair <- DESeqDataSetFromMatrix(paired_counts, paired, design = ~ donor_id + group)
dds_pair <- dds_pair[rowSums(counts(dds_pair) >= 10) >= 3, ]
dds_pair <- DESeq(dds_pair)
res_pair <- results(dds_pair, contrast = c("group", "lesion", "perilesion"), alpha = 0.05)
write.csv(as.data.frame(res_pair), file.path(outdir, "GSE245451_lesion_vs_perilesion.csv"))

# Secondary contrasts against healthy controls. These are unpaired comparisons.
for (contrast_group in c("lesion", "perilesion")) {
  submeta <- meta[meta$group %in% c(contrast_group, "healthy"), , drop = FALSE]
  submeta$group <- factor(submeta$group, levels = c("healthy", contrast_group))
  subcounts <- round(as.matrix(count_df[, rownames(submeta), drop = FALSE]))
  dds <- DESeqDataSetFromMatrix(subcounts, submeta, design = ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", contrast_group, "healthy"), alpha = 0.05)
  write.csv(as.data.frame(res), file.path(outdir, paste0("GSE245451_", contrast_group, "_vs_healthy.csv")))
}

# Save normalized counts for plots and downstream validation; do not pool datasets.
norm <- counts(dds_pair, normalized = TRUE)
write.csv(as.data.frame(norm), file.path(outdir, "GSE245451_paired_normalized_counts.csv"))
message("Primary DESeq2 contrasts complete. Review QC and design assumptions before manuscript reporting.")
