# 3_S-PrediXcan: Transcriptome-Wide Association and Cross-Tissue Integration

This module performs transcriptome-informed association testing using summary statistics. It links GWAS SNP effects to genetically predicted gene expression, then aggregates evidence across tissues and merges with other post-GWAS signals for gene prioritization.

## Biological significance

GWAS loci often implicate non-coding regulation rather than coding sequence. S-PrediXcan helps answer:

- Which genes show association through predicted expression changes?
- In which tissues is expression-mediated signal strongest?
- Are associations consistent across tissues, or tissue-specific?

For social isolation traits, this can prioritize genes with plausible neurobiological or systemic regulatory roles even when lead SNPs are intergenic.

## Statistical framework

### 1) S-PrediXcan (single tissue)

S-PrediXcan uses pre-trained expression prediction models (typically elastic-net weights from eQTL reference data) and GWAS summary statistics to test each gene:

\[
Z_g \approx \frac{\sum_j w_{gj} Z_j}{\sqrt{w_g^\top \Sigma w_g}}
\]

- \(w_{gj}\): SNP weight for gene \(g\) in a given tissue model
- \(Z_j\): GWAS SNP z-score
- \(\Sigma\): SNP covariance (LD/covariance used by the model)

Interpretation:

- Significant gene p-values suggest genetically regulated expression of that gene is associated with the phenotype in that tissue context.

### 2) S-MultiXcan / cross-tissue aggregation

This pipeline supports two modes:

1. **Official SMulTiXcan** (`USE_OFFICIAL_SMULTIXCAN=1`) with MetaXcan `SMulTiXcan.py`.
2. **Fallback Stouffer aggregation** (default), combining per-tissue z-scores:

\[
Z_{\text{Stouffer}} = \frac{\sum_{t=1}^{k} Z_t}{\sqrt{k}}
\]

This gives a robust cross-tissue summary when full official setup is unavailable.

### 3) Integrated prioritization score

`scripts/prioritize_genes.py` merges MAGMA, S-PrediXcan, S-MultiXcan, and optional coloc evidence into a combined ranking with:

- per-method significance flags,
- convergent evidence count,
- weighted score based on `-log10(p)` and coloc PP.H4.

## Scripts and workflow

Run from `5_Functional_Genomics/`:

```bash
bash 3_S-PrediXcan/scripts/run_spredixcan.sh
bash 3_S-PrediXcan/scripts/run_smultixcan.sh
python 3_S-PrediXcan/scripts/prioritize_genes.py
```

### `scripts/run_spredixcan.sh`

- Input GWAS tables: `0_munge/output/spredixcan/*.spredixcan.tsv.gz`
- Tissue models: `reference_data/predixcan/models/*.db`
- Covariances: `reference_data/predixcan/covariances/*`
- Output: `output/spredixcan/<gwas>__<tissue>.csv`

### `scripts/run_smultixcan.sh`

- Input: single-tissue S-PrediXcan outputs
- Output:
  - official mode: `output/smultixcan/smultixcan_official.tsv`
  - fallback mode: `output/smultixcan/*__smultixcan_fallback.tsv`

### `scripts/prioritize_genes.py`

- Input:
  - MAGMA gene results (`2_MAGMA/output/02_gene_analysis`)
  - S-PrediXcan results (`output/spredixcan`)
  - S-MultiXcan results (`output/smultixcan`)
  - optional coloc summary (`5_Colocalization/output/summary/coloc_prioritized.tsv`)
- Output: `output/prioritized_genes.tsv`

## Input requirements

- Harmonized GWAS columns from `0_munge`:
  - `rsid`, `chromosome`, `position`, `effect_allele`, `non_effect_allele`, `beta`, `se`, `pvalue`
- Compatible genome build between GWAS files and prediction models (this pipeline assumes b37/hg19 conventions)
- Matching covariance files for model DBs

## Interpretation tips

- Treat TWAS associations as gene-prioritization evidence, not definitive causal proof.
- Cross-tissue consistency can increase confidence, but tissue-specific hits may also be biologically meaningful.
- Strongest candidates are genes supported by multiple modules (e.g., MAGMA + TWAS + coloc).

## References and documentation

- MetaXcan repository: https://github.com/hakyimlab/MetaXcan
- PredictDB resources: https://predictdb.org/
- S-PrediXcan paper (Barbeira et al., 2018): https://www.nature.com/articles/s41467-018-03621-1
- S-MultiXcan method (Barbeira et al., 2019): https://doi.org/10.1371/journal.pgen.1007889
