# Phenotype File Creation

This directory contains the phenotype formatting workflow used to prepare UK Biobank traits for downstream GWAS and MR analyses.

## Source

- Notebook: `Pheno_Formatting.ipynb`
- Original location: `Run_1/Pheno_Formatting.ipynb`

## What this notebook does

- Loads raw phenotype files and fixes column formatting (`FID`, `IID`, phenotype columns)
- Creates social isolation traits including `Loneliness`, `FreqSoc`, and `AbilityToConfide`
- Applies recoding and missing-value handling
- Produces sex-stratified phenotype datasets
- Prepares phenotype outputs used by GWAS (BOLT-LMM) and MR workflows

## Output usage

Outputs from this stage are consumed by:

- `2_GWAS/` for BOLT-LMM association analyses
- `4_Mendelian_Randomization/` for exposure/outcome preparation
