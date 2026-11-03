suppressPackageStartupMessages(library(DESeq2))

count_file <- "data/processed/GSE151243/GSE151243_Raw_gene_counts_matrix.txt.gz"
meta_file <- "metadata/GSE151243_sample_metadata.tsv"
old_file <- "results/validation/GSE151243_lesion_vs_perilesion.csv"
new_file <- "results/validation/GSE151243_lesion_vs_perilesion_reproduced.csv"

raw <- read.delim(gzfile(count_file), check.names = FALSE)
gene_ids <- raw[[1]]
count_matrix <- as.matrix(raw[, -1, drop = FALSE])
rownames(count_matrix) <- gene_ids
storage.mode(count_matrix) <- "numeric"

meta <- read.delim(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[meta$group %in% c("lesion", "perilesion"), , drop = FALSE]

sample_columns <- vapply(meta$sample_title, function(title) {
  matches <- which(startsWith(colnames(count_matrix), title))
  if (length(matches) != 1) stop("Expected one count column for: ", title)
  colnames(count_matrix)[matches]
}, character(1))

count_matrix <- count_matrix[, sample_columns, drop = FALSE]
colnames(count_matrix) <- meta$sample_accession
rownames(meta) <- meta$sample_accession
meta <- meta[colnames(count_matrix), , drop = FALSE]

meta$group <- factor(meta$group, levels = c("perilesion", "lesion"))
meta$donor_id <- factor(meta$donor_id)

pair_table <- table(meta$donor_id, meta$group)
stopifnot(all(pair_table[, "perilesion"] == 1))
stopifnot(all(pair_table[, "lesion"] == 1))

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = meta,
  design = ~ donor_id + group
)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
cat("Genes retained after count filter:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "lesion", "perilesion"), alpha = 0.05)
new <- as.data.frame(res)
write.csv(new, new_file)
cat("Saved reproduced result:", new_file, "\n")

if (file.exists(old_file)) {
  old <- read.csv(old_file, row.names = 1, check.names = FALSE)
  common <- intersect(rownames(old), rownames(new))

  cat("Genes compared with prior result:", length(common), "\n")
  cat("Spearman correlation of log2FC:",
      cor(old[common, "log2FoldChange"], new[common, "log2FoldChange"],
          method = "spearman", use = "complete.obs"), "\n")
  cat("Median absolute log2FC difference:",
      median(abs(old[common, "log2FoldChange"] -
                 new[common, "log2FoldChange"]), na.rm = TRUE), "\n")
  cat("Prior FDR < 0.05:", sum(old$padj < 0.05, na.rm = TRUE), "\n")
  cat("Reproduced FDR < 0.05:", sum(new$padj < 0.05, na.rm = TRUE), "\n")
}
