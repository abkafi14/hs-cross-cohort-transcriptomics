suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

discovery_file <- "results/differential_expression/GSE245451_lesion_vs_perilesion.csv"
validation_file <- "results/validation/GSE151243_lesion_vs_perilesion_reproduced.csv"
old_file <- "results/validation/concordant_shared_DEGs.csv"

discovery <- read.csv(discovery_file, row.names = 1, check.names = FALSE)
discovery$symbol <- rownames(discovery)

symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(read.csv(validation_file, row.names = 1, check.names = FALSE)),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "asNA"
)

validation <- read.csv(validation_file, row.names = 1, check.names = FALSE)
validation$symbol <- unname(symbols[rownames(validation)])

validation <- validation[!is.na(validation$symbol) & !is.na(validation$padj), ]
validation <- validation[
  order(validation$padj, -abs(validation$log2FoldChange)),
]
validation <- validation[!duplicated(validation$symbol), ]

discovery <- discovery[!is.na(discovery$padj), ]
discovery <- discovery[
  order(discovery$padj, -abs(discovery$log2FoldChange)),
]
discovery <- discovery[!duplicated(discovery$symbol), ]

discovery_sig <- discovery[
  discovery$padj <= 0.05 & abs(discovery$log2FoldChange) >= 1,
  c("symbol", "log2FoldChange", "padj")
]
validation_sig <- validation[
  validation$padj <= 0.05 & abs(validation$log2FoldChange) >= 1,
  c("symbol", "log2FoldChange", "padj")
]

names(discovery_sig)[2:3] <- c("log2FoldChange_discovery", "padj_discovery")
names(validation_sig)[2:3] <- c("log2FoldChange_validation", "padj_validation")

shared <- merge(discovery_sig, validation_sig, by = "symbol")
shared$same_direction <-
  sign(shared$log2FoldChange_discovery) == sign(shared$log2FoldChange_validation)
shared$direction <- ifelse(shared$log2FoldChange_discovery > 0, "Up", "Down")
shared$worst_padj <- pmax(shared$padj_discovery, shared$padj_validation)

concordant <- shared[shared$same_direction, ]
concordant <- concordant[order(concordant$worst_padj), ]

new_file <- "results/validation/concordant_shared_DEGs_reproduced.csv"
write.csv(concordant, new_file, row.names = FALSE)

cat("Reproduced concordant genes:", nrow(concordant), "\n")
cat("Up:", sum(concordant$direction == "Up"),
    "| Down:", sum(concordant$direction == "Down"), "\n")
cat("Saved:", new_file, "\n")

if (file.exists(old_file)) {
  old <- read.csv(old_file, check.names = FALSE)
  cat("Prior consensus genes:", nrow(old), "\n")
  cat("Gene symbols in both lists:",
      length(intersect(old$symbol, concordant$symbol)), "\n")
}
