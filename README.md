# UKBB-SI-Genetics

Organized UK Biobank genetics workflow for social isolation and loneliness analyses, from phenotype construction through GWAS, MTAG, and Mendelian Randomization (MR).

## Repository Organization

### 1) `1_Phenotype_Creation/`

- Source notebook for phenotype formatting:
  - `Pheno_Formatting.ipynb`
- Builds cleaned trait files and recoded variables used in downstream analyses.

### 2) `2_GWAS/`

- BOLT-LMM pipeline scripts for three social traits across three population sets:
  - Traits: `Loneliness`, `FreqSoc`, `AbilityToConfide`
  - Populations: `EUR_MM`, `EUR_Male`, `EUR_Female`
- Includes SLURM jobs, population filtering, BOLT execution, and conversion to MTAG input format.

### 3) `3_MTAG/`

- MTAG run scripts for:
  - Social Isolation traits (`SI/`)
  - MRI traits (`MRI/`)
- Includes plotting and analysis scripts under `Analysis/SI/`.

### 4) `4_Mendelian_Randomization/`

- Starter scaffold for TwoSampleMR-based analyses.
- Contains:
  - `scripts/run_twosample_mr.R`
  - `requirements.txt`

## Pipeline Flow

1. Create and format phenotypes in `1_Phenotype_Creation/`
2. Run BOLT-LMM GWAS in `2_GWAS/`
3. Run MTAG in `3_MTAG/`
4. Run MR analyses in `4_Mendelian_Randomization/`

## Data and Results Tracking

Large data files and generated outputs are excluded from git via `.gitignore`, including:

- compressed GWAS outputs (`*.tsv.gz`, `*.stats.gz`, `*.log.gz`)
- MTAG summary statistics (`*.sumstats.txt`)
- MTAG matrix outputs (`*_omega_hat.txt`, `*_sigma_hat.txt`)
- generated figures (`*.png`)

This repository is focused on reproducible code and workflow structure.
