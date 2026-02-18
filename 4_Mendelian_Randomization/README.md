# Mendelian Randomization (Scaffold)

This directory is a starter scaffold for TwoSampleMR analyses that use GWAS/MTAG summary statistics produced in this repository.

## Intended workflow

1. Prepare exposure and outcome summary statistics from `2_GWAS/` or `3_MTAG/`.
2. Harmonize alleles between exposure and outcome datasets.
3. Run core MR estimators (IVW, Egger, weighted median).
4. Run sensitivity checks (heterogeneity, pleiotropy, leave-one-out).

## Files

- `scripts/run_twosample_mr.R`: template end-to-end TwoSampleMR pipeline for binary-coded traits
- `scripts/run_twosample_mr_continuous.R`: template pipeline for continuous-coded traits
- `requirements.txt`: required R packages

## Input expectations

Template code assumes tabular summary-stat files with at least:

- SNP identifier (`SNP` or `snpid`)
- effect allele / other allele
- effect estimate (`beta`)
- standard error (`se`)
- p-value (`pval`)
- sample size (`n`)

Update the paths and column mappings in `scripts/run_twosample_mr.R` before running.

## Related modules

This module runs in parallel with `5_Functional_Genomics/`, which provides complementary post-GWAS interpretation:

- LDSC heritability and genetic correlations
- MAGMA gene-based and pathway analyses
- S-PrediXcan transcriptome-wide association
- SuSiE fine-mapping and GWAS-eQTL colocalization
- Integrated gene prioritization

See `../5_Functional_Genomics/README.md` for details.
