suppressPackageStartupMessages({
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

read_result <- function(file, cohort, id_type) {
  x <- read.csv(file, row.names = 1, check.names = FALSE)
  x$feature_id <- rownames(x)

  if (id_type == "ENSEMBL") {
    symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = x$feature_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    x$gene_label <- unname(symbols)
    missing <- is.na(x$gene_label) | x$gene_label == ""
    x$gene_label[missing] <- x$feature_id[missing]
  } else {
    x$gene_label <- x$feature_id
  }

  x$cohort <- cohort
  x$neglog10_FDR <- -log10(pmax(x$padj, 1e-50, na.rm = FALSE))
  x$neglog10_FDR[is.na(x$padj)] <- 0

  x$status <- "Not significant"
  sig <- !is.na(x$padj) & x$padj < 0.05
  x$status[sig & x$log2FoldChange >= 1] <- "Up"
  x$status[sig & x$log2FoldChange <= -1] <- "Down"
  x$status[sig & abs(x$log2FoldChange) < 1] <- "FDR only"

  x
}

discovery <- read_result(
  "results/differential_expression/GSE245451_lesion_vs_perilesion.csv",
  "GSE245451", "SYMBOL"
)

validation <- read_result(
  "results/validation/GSE151243_lesion_vs_perilesion.csv",
  "GSE151243", "ENSEMBL"
)

dat <- rbind(discovery, validation)

labels <- do.call(rbind, lapply(split(dat, dat$cohort), function(x) {
  x <- x[!is.na(x$padj), ]
  up <- head(x[x$status == "Up", ][order(x[x$status == "Up", "padj"]), ], 6)
  down <- head(x[x$status == "Down", ][order(x[x$status == "Down", "padj"]), ], 6)
  rbind(up, down)
}))

p <- ggplot(dat, aes(x = log2FoldChange, y = neglog10_FDR, color = status)) +
  geom_point(alpha = 0.65, size = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_text(
    data = labels,
    aes(label = gene_label),
    size = 2.5,
    check_overlap = TRUE,
    vjust = -0.5,
    show.legend = FALSE
  ) +
  facet_wrap(~cohort, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c(
    "Not significant" = "grey75",
    "FDR only" = "#7B8DA6",
    "Up" = "#C83E4D",
    "Down" = "#3B6FB6"
  )) +
  labs(
    x = "log2 fold change: lesion vs perilesional",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/Figure_primary_volcano.pdf", p, width = 12, height = 6)
ggsave("results/figures/Figure_primary_volcano.png", p, width = 12, height = 6, dpi = 300)

cat("Saved: results/figures/Figure_primary_volcano.pdf and .png\n")
cat("\nDiscovery categories:\n")
print(table(discovery$status))
cat("\nValidation categories:\n")
print(table(validation$status))
