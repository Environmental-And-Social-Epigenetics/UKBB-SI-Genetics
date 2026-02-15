# UKBB-SI-Genetics

Reproducible genetics workflow for social isolation and loneliness analyses in the UK Biobank, spanning phenotype construction, GWAS, multi-trait analysis (MTAG), and Mendelian Randomization (MR).

## Background

This project investigates the genetic architecture of social isolation using UK Biobank data. The phenotype definitions and binary coding follow:

> Day, F.R., Ong, K.K. & Perry, J.R.B. **Elucidating the genetic basis of social interaction and isolation.** *Nature Communications* 9, 2457 (2018). https://doi.org/10.1038/s41467-018-04930-1

In addition to replicating the binary GWAS from Day et al., we extend the analysis with a **continuous coding** of the same traits and apply MTAG to jointly analyze social isolation and neuroimaging phenotypes.

## Traits and Populations

### Social Isolation Traits

Three derived traits are used, each oriented so that higher values indicate greater isolation:

| Trait | UKBB Fields | Binary Coding (Day et al.) | Continuous Coding (this study) |
|---|---|---|---|
| **Loneliness** | 2020 | Yes = case, No = control | Raw 0/1 |
| **AbilityToConfide** | 2110 | Never/almost never = case | 5 - raw value |
| **FreqSoc** | 1031 + 709 | Lives alone & rarely visited = case | Mean of z-scored components |

### Population Stratifications

Each GWAS is run across three European-ancestry subsets:

| Population | Description |
|---|---|
| `EUR_MM` | All European-ancestry participants (mixed male/female) |
| `EUR_Male` | European-ancestry males only |
| `EUR_Female` | European-ancestry females only |

## Pipeline Overview

```
1_Phenotype_Creation          2_GWAS                          3_MTAG                  4_MR
─────────────────────         ──────────────────────────       ─────────────────       ──────────────
                              Filter populations
UKBB basket file              (EUR_MM, EUR_Male, EUR_Female)
       │                              │
       ├── Binary coding ──────► BOLT-LMM binary GWAS         MTAG SI binary ───┐
       │                              │                              │            │
       │                         Convert to MTAG format ────►  Analysis/plots    │
       │                                                                         ├──► TwoSampleMR
       ├── Continuous coding ──► BOLT-LMM continuous GWAS      MTAG SI cont. ───┤
                                      │                                          │
                                 Convert to MTAG format ────►                    │
                                                               MTAG MRI ────────┘
```

Each GWAS step runs **3 traits x 3 populations = 9 jobs** per coding scheme (18 total).

## Repository Organization

### `1_Phenotype_Creation/`

Extracts raw phenotype columns from the UKBB basket file and formats them into PLINK-compatible phenotype files for downstream GWAS.

- **`extract_phenotype.sh`** / **`extract_all_phenotypes.sh`** -- Shell scripts to pull fields 2020, 2110, 1031, and 709 from the basket file
- **`Binary_Pheno_Formatting.ipynb`** -- Applies Day et al. binary case/control coding
- **`Continuous_Pheno_Formatting.ipynb`** -- Applies continuous coding (this study)
- Outputs: `isolation_run_binary.tsv.gz`, `isolation_run_continuous.tsv.gz`

### `2_GWAS/`

BOLT-LMM mixed-model GWAS pipeline running 3 traits x 3 populations = 9 jobs per coding scheme (18 total).

- **Step 0**: `0a_filter_populations.sbatch.sh` -- Filters phenotype/covariate files to each population subset
- **Step 0b**: `0b_test_run.sbatch.sh` / `0b_test_run_continuous.sbatch.sh` -- Single-job validation runs
- **Step 1**: `1_run_bolt_lmm.sbatch.sh` / `1b_run_bolt_lmm_continuous.sbatch.sh` -- Full SLURM array jobs (9 tasks each)
- **Step 2**: `2_convert_to_MTAG.sbatch.sh` / `2b_convert_to_MTAG_continuous.sbatch.sh` -- Convert BOLT-LMM output to MTAG input format
- Core execution: `run_single_phenotype.sh` -- Runs BOLT-LMM for one phenotype-population combination
- Key parameters: ~1.3M autosomal variants, ~444K model SNPs for GRM, covariates = age + sex + array

### `3_MTAG/`

Multi-Trait Analysis of GWAS summary statistics to boost power by leveraging genetic correlations between related traits.

- **`SI/`** -- Binary-coded social isolation MTAG (3 populations x 3 traits)
- **`SI_continuous/`** -- Continuous-coded social isolation MTAG (3 populations x 3 traits)
- **`MRI/`** -- MRI white-matter trait MTAG (FA, MD, MO, OD across 3 populations)
- **`Analysis/SI/`** -- Downstream analysis including Manhattan plot generation and result summaries

### `4_Mendelian_Randomization/`

Template scripts for two-sample MR analyses testing causal relationships between social isolation traits and neuroimaging outcomes.

- **`scripts/run_twosample_mr.R`** -- MR pipeline for binary-coded traits
- **`scripts/run_twosample_mr_continuous.R`** -- MR pipeline for continuous-coded traits
- Methods: IVW, MR-Egger, weighted median, with heterogeneity and pleiotropy sensitivity checks
- **`requirements.txt`** -- R package dependencies (TwoSampleMR, MendelianRandomization, MRPRESSO, etc.)

## Compute Environment

All GWAS and MTAG jobs are designed for the MIT Luria SLURM cluster:

- BOLT-LMM runs use the `bolt_lmm` conda environment
- Python conversion scripts use the `Python_Analysis` conda environment
- Typical GWAS job: 100 GB RAM, 32 CPUs, 1--2 hours per phenotype-population combination

## Data and Results Tracking

Large data files and generated outputs are excluded from git via `.gitignore`:

- Compressed GWAS outputs (`*.tsv.gz`, `*.stats.gz`, `*.log.gz`)
- MTAG summary statistics (`*.sumstats.txt`)
- MTAG matrix outputs (`*_omega_hat.txt`, `*_sigma_hat.txt`)
- Generated figures (`*.png`)
- SLURM log files (`*.out`, `*.err`)

This repository tracks reproducible code and workflow structure only.

## References

- Day, F.R., Ong, K.K. & Perry, J.R.B. Elucidating the genetic basis of social interaction and isolation. *Nat Commun* 9, 2457 (2018). https://doi.org/10.1038/s41467-018-04930-1
