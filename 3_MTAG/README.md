# MTAG Analysis

This directory contains scripts and analysis assets for Multi-Trait Analysis of GWAS (MTAG).

## Layout

- `SI/`: MTAG run scripts for Social Isolation traits (`Loneliness`, `FreqSoc`, `AbilityToConfide`)
- `MRI/`: MTAG run scripts for MRI traits (`FA`, `MD`, `MO`, `OD`)
- `Analysis/SI/`: plotting and summary files for downstream interpretation

## Notes

- Raw MTAG sumstats and result outputs are intentionally excluded from version control via `.gitignore`.
- Use `run_mtag.sh` scripts in each population-specific subdirectory to launch MTAG with the expected local environment paths.
