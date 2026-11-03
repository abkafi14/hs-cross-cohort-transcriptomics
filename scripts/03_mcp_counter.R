suppressPackageStartupMessages({
  library(DESeq2)
  library(MCPcounter)
})

outdir <- "results/cell_composition"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prepare_vst <- function(counts, meta, accession_col, id_col_type) {
  meta <- meta[meta$group %in% c("lesion", "perilesion"), , drop = FALSE]
  meta$group <- factor(meta$group, levels = c("perilesion", "lesion"))
  meta$donor_id <- factor(meta$donor_id)

  counts <- counts[, meta[[accession_col]], drop = FALSE]
  rownames(meta) <- meta[[accession_col]]
  stopifnot(identical(colnames(counts), rownames(meta)))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = meta,
    design = ~ donor_id + group
  )
  dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 3, ]
  vst_obj <- DESeq2::vst(dds, blind = FALSE)
  list(expression = SummarizedExperiment::assay(vst_obj), meta = meta)
}

run_mcp <- function(expression, meta, id_type, label) {
  scores <- MCPcounter::MCPcounter.estimate(
    expression = expression,
    featuresType = id_type
  )

  if (anyNA(scores)) {
    warning(label, ": MCP-counter returned NA scores for one or more populations.")
  }

  score_file <- file.path(outdir, paste0(label, "_MCPcounter_scores.csv"))
  write.csv(t(scores), score_file)

  tests <- lapply(rownames(scores), function(population) {
    dat <- data.frame(
      sample_accession = colnames(scores),
      score = as.numeric(scores[population, ]),
      stringsAsFactors = FALSE
    )
    dat <- merge(
      dat,
      meta[, c("sample_accession", "group", "donor_id")],
      by = "sample_accession"
    )

    lesion <- dat[dat$group == "lesion", c("donor_id", "score")]
    perilesion <- dat[dat$group == "perilesion", c("donor_id", "score")]
    names(lesion)[2] <- "lesion_score"
    names(perilesion)[2] <- "perilesion_score"
    paired <- merge(lesion, perilesion, by = "donor_id")
    difference <- paired$lesion_score - paired$perilesion_score

    p <- if (length(difference) >= 2) {
      stats::wilcox.test(
        paired$lesion_score,
        paired$perilesion_score,
        paired = TRUE,
        exact = FALSE
      )$p.value
    } else {
      NA_real_
    }

    data.frame(
      cell_population = population,
      n_donor_pairs = length(difference),
      median_lesion_minus_perilesion = median(difference, na.rm = TRUE),
      mean_lesion_minus_perilesion = mean(difference, na.rm = TRUE),
      p_value = p
    )
  })

  tests <- do.call(rbind, tests)
  tests$FDR <- p.adjust(tests$p_value, method = "BH")
  test_file <- file.path(outdir, paste0(label, "_MCPcounter_paired_tests.csv"))
  write.csv(tests, test_file, row.names = FALSE)

  cat("\n", label, "| score matrix:", nrow(scores), "populations x",
      ncol(scores), "samples\n", sep = " ")
  print(tests[order(tests$FDR), ], row.names = FALSE)
  cat("Saved:", score_file, "and", test_file, "\n")

  invisible(list(scores = scores, tests = tests))
}

discovery_meta <- read.delim(
  "metadata/GSE245451_sample_metadata.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
discovery_counts <- read.delim(
  "data/processed/GSE245451/GSE245451_counts_matrix.tsv",
  row.names = 1,
  check.names = FALSE
)
discovery <- prepare_vst(
  discovery_counts,
  discovery_meta,
  "sample_accession",
  "HUGO_symbols"
)
run_mcp(
  discovery$expression,
  discovery$meta,
  "HUGO_symbols",
  "GSE245451"
)

validation_meta <- read.delim(
  "metadata/GSE151243_sample_metadata.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
validation_raw <- read.delim(
  gzfile("data/processed/GSE151243/GSE151243_Raw_gene_counts_matrix.txt.gz"),
  check.names = FALSE
)

validation_gene_ids <- validation_raw[[1]]
validation_counts <- as.matrix(validation_raw[, -1, drop = FALSE])
rownames(validation_counts) <- validation_gene_ids
storage.mode(validation_counts) <- "numeric"

sample_columns <- vapply(validation_meta$sample_title, function(title) {
  matches <- which(startsWith(colnames(validation_counts), title))
  if (length(matches) != 1) {
    stop("Expected one validation count column for sample title: ", title)
  }
  colnames(validation_counts)[matches]
}, character(1))

validation_counts <- validation_counts[, sample_columns, drop = FALSE]
colnames(validation_counts) <- validation_meta$sample_accession

validation <- prepare_vst(
  validation_counts,
  validation_meta,
  "sample_accession",
  "ENSEMBL_ID"
)
run_mcp(
  validation$expression,
  validation$meta,
  "ENSEMBL_ID",
  "GSE151243"
)

cat("\nMCP-counter analysis complete. Results saved in:", outdir, "\n")
