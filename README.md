# UKBB-SI-Genetics

Reproducible genetics workflow for social isolation and loneliness analyses in the UK Biobank, spanning phenotype construction, GWAS, multi-trait analysis (MTAG), functional genomics, and Mendelian Randomization (MR).

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
1_Phenotype         2_GWAS                     3_MTAG               4_MR              5_Functional_Genomics
────────────        ─────────────────          ─────────────        ──────────────     ─────────────────────
                    Filter populations
UKBB basket         (EUR_MM, EUR_Male,
       │             EUR_Female)
       │                    │
       ├── Binary ────► BOLT-LMM binary        MTAG SI binary ───┐                  ┌── LDSC (h2, rg)
       │                    │                        │             │                  ├── MAGMA (gene, pathway, tissue)
       │               Convert to MTAG ────►   Analysis/plots     ├──► TwoSampleMR   ├── S-PrediXcan / S-MultiXcan
       │                                                          │                  ├── Fine-mapping (SuSiE)
       ├── Continuous ► BOLT-LMM cont.         MTAG SI cont. ────┤                  ├── Colocalization (coloc)
                            │                                     │                  └── Visualization + gene table
                       Convert to MTAG ────►                      │
                                               MTAG MRI ─────────┘
                                                     │
                                                     └────────────────────────────────► 5_Functional_Genomics
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

### `5_Functional_Genomics/`

Comprehensive post-GWAS functional interpretation pipeline operating entirely on summary statistics from `3_MTAG/`.

- **`reference_data/`** -- Download script for LD scores, 1000G EUR, gene sets, GTEx models, eQTL data, and external GWAS sumstats
- **`0_munge/`** -- Reformats MTAG outputs for each downstream tool (LDSC, MAGMA, S-PrediXcan)
- **`1_LDSC/`** -- SNP heritability, genetic correlations (internal and external), partitioned h2 (baselineLD + cell-type)
- **`2_MAGMA/`** -- Gene-based association, gene-set (GO, KEGG, MSigDB) enrichment, and tissue expression analysis
- **`3_S-PrediXcan/`** -- Per-tissue and cross-tissue transcriptome-wide association, gene prioritization
- **`4_Fine_Mapping/`** -- Locus definition, LD matrix computation, SuSiE credible sets
- **`5_Colocalization/`** -- GWAS-eQTL colocalization with coloc, colocalized gene summary
- **`6_Visualization/`** -- Integrated plots (rg heatmap, partitioned h2, tissue enrichment, Miami plot, regional loci) and gene prioritization table
- **`environment.yml`** -- Conda environment covering all Python, R, and CLI dependencies

## Compute Environment

GWAS and MTAG jobs are designed for the MIT Luria SLURM cluster. Post-GWAS functional genomics runs locally on a laptop/workstation.

- BOLT-LMM runs use the `bolt_lmm` conda environment
- Python conversion scripts use the `Python_Analysis` conda environment
- Typical GWAS job: 100 GB RAM, 32 CPUs, 1--2 hours per phenotype-population combination
- **Functional genomics** (`5_Functional_Genomics/`) uses the `functional_genomics` conda environment (see `5_Functional_Genomics/environment.yml`). All analyses run locally on summary statistics only.

## Data and Results Tracking

Large data files and generated outputs are excluded from git via `.gitignore`:

- Compressed GWAS outputs (`*.tsv.gz`, `*.stats.gz`, `*.log.gz`)
- MTAG summary statistics (`*.sumstats.txt`)
- MTAG matrix outputs (`*_omega_hat.txt`, `*_sigma_hat.txt`)
- Generated figures (`*.png`)
- SLURM log files (`*.out`, `*.err`)
- Functional genomics reference data (`5_Functional_Genomics/reference_data/`)
- Munged summary statistics and analysis outputs under `5_Functional_Genomics/*/output/`

This repository tracks reproducible code and workflow structure only.

## References

- Day, F.R., Ong, K.K. & Perry, J.R.B. Elucidating the genetic basis of social interaction and isolation. *Nat Commun* 9, 2457 (2018). https://doi.org/10.1038/s41467-018-04930-1
