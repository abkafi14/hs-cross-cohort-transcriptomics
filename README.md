# Cross-cohort transcriptomic analysis of hidradenitis suppurativa

Reproducible R analysis of two public hidradenitis suppurativa (HS) RNA-seq cohorts, covering paired differential expression, cross-cohort validation, ranked pathway enrichment, and MCP-counter cell-abundance estimates.

## Data sources

- **Discovery:** [GSE245451](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245451) — five matched lesion/perilesional pairs and six healthy samples.
- **Validation:** [GSE151243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151243) — 20 matched lesion/perilesional pairs.

Expression count matrices remain hosted by GEO and are not included in this repository.

## Repository contents

- `scripts/01_deseq2_discovery.R`: donor-blocked lesion-versus-perilesional analysis and secondary comparisons against healthy samples.
- `metadata/`: sample metadata for both cohorts.
- `results/GSE245451_count_qc.tsv`: sample-level count quality-control summary.
- `results/concordant_shared_DEGs.csv`: 595 genes with significant, direction-concordant changes in both cohorts under the project's stated thresholds.

The discovery script expects the processed count matrix at `data/processed/GSE245451/GSE245451_counts_matrix.tsv`. Obtain the public processed counts from GEO before running it. The script requires R and DESeq2.

## Analysis overview

- Donor-blocked differential expression for matched lesion and perilesional samples.
- Independent analysis of each cohort, followed by cross-cohort comparison.
- Preranked Reactome and Gene Ontology Biological Process enrichment.
- MCP-counter estimates of immune and stromal cell abundance, analyzed separately by cohort.

## Summary of current results

The current analysis identified 595 shared, direction-concordant genes (392 higher and 203 lower in lesions), 50 replicated Reactome pathways, and 531 replicated GO Biological Process terms. MCP-counter comparisons showed no FDR-significant population differences in the five-pair discovery cohort; eight of ten populations had significantly higher lesion scores in the 20-pair validation cohort.

These bulk-expression and inferred-abundance results do not establish cell-intrinsic mechanisms or therapeutic effects.

## Repository status

Additional analysis scripts, software/session information, other derived result tables, and figure files still need to be added and checked before this repository can be cited as a complete reproducibility archive. No versioned release or DOI has been issued yet.

## Citation

A citation and versioned archive DOI will be added after the repository contents and manuscript are finalized.
