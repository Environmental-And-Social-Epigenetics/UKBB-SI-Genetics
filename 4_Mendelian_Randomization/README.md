# Mendelian Randomization (Scaffold)

This directory is a starter scaffold for TwoSampleMR analyses that use GWAS/MTAG summary statistics produced in this repository.

## Intended workflow

1. Prepare exposure and outcome summary statistics from `2_GWAS/` or `3_MTAG/`.
2. Harmonize alleles between exposure and outcome datasets.
3. Run core MR estimators (IVW, Egger, weighted median).
4. Run sensitivity checks (heterogeneity, pleiotropy, leave-one-out).

## Files

- `scripts/run_twosample_mr.R`: template end-to-end TwoSampleMR pipeline
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
