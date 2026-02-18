# MTAG Analysis

This directory contains scripts and analysis assets for Multi-Trait Analysis of GWAS (MTAG).

## Layout

- `SI/`: MTAG run scripts for Social Isolation traits (`Loneliness`, `FreqSoc`, `AbilityToConfide`)
  - `SI/EUR_MM/`, `SI/EUR_Male_MM/`, `SI/EUR_Female_MM/` -- population-stratified runs
- `SI_continuous/`: Continuous-coded Social Isolation MTAG (same structure as `SI/`)
- `MRI/`: MTAG run scripts for MRI white-matter traits (`FA`, `MD`, `MO`, `OD`)
- `Analysis/SI/`: plotting and summary files for downstream interpretation
- `run_all_mtag_local.sh`: convenience script to run all 6 SI MTAG analyses locally

## Output Files

Each MTAG run produces (inside `results/` subdirectory):

| File | Description |
|---|---|
| `*_trait_N.txt` | Per-trait MTAG summary statistics (columns: `SNP CHR BP A1 A2 Z N FRQ mtag_beta mtag_se mtag_z mtag_pval`) |
| `*_omega_hat.txt` | Estimated genetic covariance matrix |
| `*_sigma_hat.txt` | Estimated residual covariance matrix |
| `*.log` | Execution log with chi-squared statistics, weight factors, and run time |

Trait mapping for SI analyses: Trait 1 = AbilityToConfide, Trait 2 = FreqSoc, Trait 3 = Loneliness.

## Downstream Analysis

MTAG trait output files serve as primary inputs for two downstream modules:

- **`4_Mendelian_Randomization/`** -- Two-sample MR testing causal relationships between SI traits and neuroimaging outcomes.
- **`5_Functional_Genomics/`** -- Comprehensive post-GWAS functional interpretation pipeline including LDSC heritability/genetic correlations, MAGMA gene/pathway/tissue enrichment, S-PrediXcan TWAS, SuSiE fine-mapping, GWAS-eQTL colocalization, and integrated gene prioritization. See `5_Functional_Genomics/README.md` for full documentation and run instructions.

## Notes

- Raw MTAG sumstats and result outputs are intentionally excluded from version control via `.gitignore`.
- Use `run_mtag.sh` scripts in each population-specific subdirectory to launch MTAG with the expected local environment paths.
- All analyses assume GRCh37/hg19 genome build (matching BOLT-LMM output).
