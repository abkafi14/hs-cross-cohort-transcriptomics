suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
})

discovery_file <- "results/differential_expression/GSE245451_lesion_vs_perilesion.csv"
validation_file <- "results/validation/GSE151243_lesion_vs_perilesion.csv"
outdir <- "results/gsea"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_de <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  x <- read.csv(file, row.names = 1, check.names = FALSE)
  if (!"stat" %in% names(x)) stop("No stat column in: ", file)
  x
}

make_rank <- function(x, keytype) {
  stats <- x$stat
  names(stats) <- rownames(x)
  stats <- stats[is.finite(stats) & !is.na(names(stats)) & names(stats) != ""]

  ids <- unique(names(stats))
  ann <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ids,
    keytype = keytype,
    columns = "ENTREZID"
  )
  ann <- ann[!is.na(ann$ENTREZID) & !is.na(ann[[keytype]]), ]

  mapped <- data.frame(
    ENTREZID = as.character(ann$ENTREZID),
    stat = unname(stats[match(ann[[keytype]], names(stats))])
  )
  mapped <- mapped[is.finite(mapped$stat), ]
  mapped <- mapped[order(abs(mapped$stat), decreasing = TRUE), ]
  mapped <- mapped[!duplicated(mapped$ENTREZID), ]

  gene_list <- mapped$stat
  names(gene_list) <- mapped$ENTREZID
  set.seed(20261103)
  gene_list <- gene_list + stats::runif(length(gene_list), min = -1e-10, max = 1e-10)
  sort(gene_list, decreasing = TRUE)
}

run_gsea <- function(gene_list, label) {
  reactome <- ReactomePA::gsePathway(
    geneList = gene_list,
    organism = "human",
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    seed = TRUE,
    verbose = FALSE
  )

  go_bp <- clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    seed = TRUE,
    verbose = FALSE
  )

  reactome_df <- as.data.frame(reactome)
  go_df <- as.data.frame(go_bp)

  write.csv(reactome_df, file.path(outdir, paste0(label, "_Reactome.csv")), row.names = FALSE)
  write.csv(go_df, file.path(outdir, paste0(label, "_GO_BP.csv")), row.names = FALSE)
  write.csv(
    data.frame(ENTREZID = names(gene_list), stat = unname(gene_list)),
    file.path(outdir, paste0(label, "_ranked_genes.csv")),
    row.names = FALSE
  )

  cat(
    label, "| mapped ranked genes:", length(gene_list),
    "| Reactome FDR<0.05:", sum(reactome_df$p.adjust < 0.05, na.rm = TRUE),
    "| GO BP FDR<0.05:", sum(go_df$p.adjust < 0.05, na.rm = TRUE), "\n"
  )

  list(reactome = reactome_df, go = go_df)
}

compare_gsea <- function(a, b, name, label_a, label_b) {
  if (!all(c("ID", "NES", "p.adjust") %in% names(a)) ||
      !all(c("ID", "NES", "p.adjust") %in% names(b))) {
    cat("Could not compare", name, "- missing result columns.\n")
    return(invisible(NULL))
  }

  names(a)[names(a) %in% c("NES", "p.adjust")] <- c(
    paste0("NES_", label_a), paste0("FDR_", label_a)
  )
  names(b)[names(b) %in% c("NES", "p.adjust")] <- c(
    paste0("NES_", label_b), paste0("FDR_", label_b)
  )

  out <- merge(a, b, by = "ID", suffixes = c("_discovery", "_validation"))
  out$direction_concordant <- sign(out[[paste0("NES_", label_a)]]) ==
    sign(out[[paste0("NES_", label_b)]])

  write.csv(out, file.path(outdir, paste0(name, "_cohort_comparison.csv")), row.names = FALSE)
  cat(
    name, "| pathways found in both:", nrow(out),
    "| same NES direction:", sum(out$direction_concordant, na.rm = TRUE), "\n"
  )
}

discovery <- read_de(discovery_file)
validation <- read_de(validation_file)

rank_discovery <- make_rank(discovery, "SYMBOL")
rank_validation <- make_rank(validation, "ENSEMBL")

gsea_discovery <- run_gsea(rank_discovery, "GSE245451")
gsea_validation <- run_gsea(rank_validation, "GSE151243")

compare_gsea(
  gsea_discovery$reactome, gsea_validation$reactome,
  "Reactome", "GSE245451", "GSE151243"
)
compare_gsea(
  gsea_discovery$go, gsea_validation$go,
  "GO_BP", "GSE245451", "GSE151243"
)

cat("GSEA outputs saved in:", outdir, "\n")
