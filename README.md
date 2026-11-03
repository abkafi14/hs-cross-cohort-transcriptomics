# Cross-cohort transcriptomic analysis of hidradenitis suppurativa

R analysis of two public hidradenitis suppurativa (HS) RNA-seq cohorts: paired differential expression, cross-cohort validation, pathway enrichment, and MCP-counter cell-abundance estimates.

## Data sources

- **Discovery:** [GSE245451](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245451) — five matched lesion/perilesional pairs and six healthy samples.
- **Validation:** [GSE151243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151243) — 20 matched lesion/perilesional pairs.

## Repository contents

- `data/processed/GSE245451/GSE245451_counts_matrix.tsv` — discovery count matrix.
- `data/processed/GSE151243/GSE151243_Raw_gene_counts_matrix.txt.gz` — validation count matrix.
- `scripts/01_deseq2_discovery.R` — discovery differential expression.
- `scripts/02_ranked_gsea.R` — ranked Reactome and GO Biological Process enrichment.
- `scripts/03_mcp_counter.R` — MCP-counter cohort analysis.
- `scripts/04_volcano_primary.R` — volcano plot generation.
- `scripts/05_deseq2_validation_audit.R` and `scripts/06_reproduce_shared_DEGs.R` — validation reproduction and shared-gene analysis.
- `metadata/` — cohort sample sheets, sequencing run manifests, and software versions.
- `results/` — selected gene and pathway tables, MCP-counter outputs, QC summary, and figure PDFs.

## Reproducing the analysis

Run commands from the repository root. The count matrices are included at the paths expected by the scripts. The workflows use R/Bioconductor packages; versions recorded for the analysis are listed in `metadata/software_versions.txt`. Check each script for its required inputs and run order.

## Summary of current results

The analysis identified 595 shared, direction-concordant genes (392 higher and 203 lower in lesions), 50 replicated Reactome pathways, and 531 replicated GO Biological Process terms. MCP-counter comparisons showed no FDR-significant population differences in the five-pair discovery cohort; eight of ten populations had significantly higher lesion scores in the 20-pair validation cohort.

These bulk-expression and inferred-abundance results do not establish cell-intrinsic mechanisms or therapeutic effects.

## Repository status

The repository contains the analysis scripts, processed count matrices, selected derived result tables, and figure files. The versioned v1.0.0 release is archived on Zenodo: https://doi.org/10.5281/zenodo.23913669.

## Citation

Please cite the versioned archive:

Kafi AH. Cross-cohort transcriptomic analysis of hidradenitis suppurativa. Version 1.0.0. Zenodo. 2026. https://doi.org/10.5281/zenodo.23913669.
