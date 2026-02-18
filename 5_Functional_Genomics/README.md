# Functional Genomics (Post-GWAS)

This module runs local, summary-statistics-based functional interpretation analyses on MTAG/GWAS outputs from `3_MTAG/`.

It is organized as a reproducible pipeline:

1. Reference download and setup
2. Summary-stat munging
3. LDSC heritability/genetic correlation analysis
4. MAGMA gene/pathway/tissue analysis
5. S-PrediXcan / S-MultiXcan (or fallback cross-tissue aggregation)
6. Fine-mapping (locus definition, LD matrices, SuSiE)
7. Colocalization (GWAS vs eQTL)
8. Visualization and integrated gene prioritization

## Directory Layout

```text
5_Functional_Genomics/
├── reference_data/
│   └── download_references.sh
├── 0_munge/
│   └── scripts/munge_sumstats.py
├── 1_LDSC/
│   ├── config/
│   │   ├── external_traits.tsv
│   │   └── cell_type_groups.ldcts
│   └── scripts/
├── 2_MAGMA/
│   └── scripts/
├── 3_S-PrediXcan/
│   └── scripts/
├── 4_Fine_Mapping/
│   └── scripts/
├── 5_Colocalization/
│   ├── config/eqtl_manifest.tsv
│   └── scripts/
├── 6_Visualization/
│   └── scripts/
└── environment.yml
```

## Inputs

- Primary input files are MTAG trait outputs, for example:
  - `3_MTAG/SI/EUR_MM/results/SI_EUR_MM_Output_trait_3.txt`
  - `3_MTAG/SI_continuous/EUR_MM/results/SI_EUR_MM_Output_continuous_trait_3.txt`
- Expected columns:
  - `SNP CHR BP A1 A2 Z N FRQ mtag_beta mtag_se mtag_z mtag_pval`

## Environment Setup

```bash
cd UKBB-SI-Genetics/5_Functional_Genomics
conda env create -f environment.yml
conda activate functional_genomics
```

### External tools not bundled by conda env

- **MAGMA**: standalone binary (`magma`) must be on `PATH`.
- **LDSC**: expected at `reference_data/ldsc/ldsc/ldsc.py` (download script can clone repository).
- **MetaXcan**: clone repository and set `SPREDIXCAN_SCRIPT` / `SMULTIXCAN_SCRIPT` to corresponding scripts.

## Step 0: Download Reference Data

```bash
bash reference_data/download_references.sh
```

Notes:
- Some resources require manual URL setup due licensing/login (MSigDB, some GTEx/PredictDB mirrors, some external GWAS datasets).
- Follow script messages and set required environment variables before re-running.

## Step 1: Munge Summary Statistics

```bash
python 0_munge/scripts/munge_sumstats.py \
  --run-ldsc-munge \
  --hm3-snplist reference_data/ldsc/w_hm3.snplist \
  --overwrite
```

Outputs:
- `0_munge/output/ldsc_input/*.ldsc.tsv.gz`
- `0_munge/output/ldsc_munged/*.sumstats.gz`
- `0_munge/output/magma/*.genes.raw`
- `0_munge/output/spredixcan/*.spredixcan.tsv.gz`
- `0_munge/output/manifest.tsv`

## Step 2: LDSC

```bash
bash 1_LDSC/scripts/run_h2.sh
bash 1_LDSC/scripts/run_rg_internal.sh
bash 1_LDSC/scripts/run_rg_external.sh
bash 1_LDSC/scripts/run_partitioned_h2.sh
bash 1_LDSC/scripts/run_celltype_h2.sh
```

Before `run_rg_external.sh`:
- Edit `1_LDSC/config/external_traits.tsv` to provide local paths for external trait `.sumstats.gz` files.

Before `run_celltype_h2.sh`:
- Edit `1_LDSC/config/cell_type_groups.ldcts` with real annotation prefixes.

## Step 3: MAGMA

```bash
bash 2_MAGMA/scripts/01_annotate.sh
bash 2_MAGMA/scripts/02_gene_analysis.sh
bash 2_MAGMA/scripts/03_geneset_analysis.sh
bash 2_MAGMA/scripts/04_tissue_expression.sh
```

Optional overrides:
- `MAGMA_BIN` (path/name of binary)
- `BFILE_PREFIX` for 1000G EUR PLINK reference
- `GENE_COVAR_FILE` for tissue expression covariates

## Step 4: S-PrediXcan / S-MultiXcan

```bash
bash 3_S-PrediXcan/scripts/run_spredixcan.sh
bash 3_S-PrediXcan/scripts/run_smultixcan.sh
python 3_S-PrediXcan/scripts/prioritize_genes.py
```

If official `SMulTiXcan.py` is not configured, `run_smultixcan.sh` performs a fallback cross-tissue Stouffer aggregation.

## Step 5: Fine-Mapping

Example for one MTAG trait file:

```bash
python 4_Fine_Mapping/scripts/01_define_loci.py \
  --sumstats 3_MTAG/SI/EUR_MM/results/SI_EUR_MM_Output_trait_3.txt \
  --plink-bfile 5_Functional_Genomics/reference_data/1000G_EUR/g1000_eur \
  --allow-distance-fallback
```

Then:

```bash
bash 4_Fine_Mapping/scripts/02_compute_ld.sh
Rscript 4_Fine_Mapping/scripts/03_run_susie.R \
  --sumstats 3_MTAG/SI/EUR_MM/results/SI_EUR_MM_Output_trait_3.txt \
  --loci 4_Fine_Mapping/output/01_defined_loci/loci.loci.tsv \
  --ld-dir 4_Fine_Mapping/output/02_ld_matrices \
  --out-dir 4_Fine_Mapping/output/03_susie
```

## Step 6: Colocalization

1. Update `5_Colocalization/config/eqtl_manifest.tsv` with real eQTL file paths and column mappings.
2. Run:

```bash
Rscript 5_Colocalization/scripts/run_coloc.R \
  --gwas 3_MTAG/SI/EUR_MM/results/SI_EUR_MM_Output_trait_3.txt \
  --loci 4_Fine_Mapping/output/01_defined_loci/loci.loci.tsv \
  --eqtl-manifest 5_Colocalization/config/eqtl_manifest.tsv \
  --out-dir 5_Colocalization/output

Rscript 5_Colocalization/scripts/summarize_coloc.R
```

## Step 7: Visualization and Integration

```bash
python 6_Visualization/scripts/plot_rg_heatmap.py
python 6_Visualization/scripts/plot_partitioned_h2.py
python 6_Visualization/scripts/plot_celltype_enrichment.py
python 6_Visualization/scripts/plot_magma_tissue_enrichment.py
python 6_Visualization/scripts/build_gene_prioritization_table.py
python 6_Visualization/scripts/plot_regional_loci.py --sumstats 3_MTAG/SI/EUR_MM/results/SI_EUR_MM_Output_trait_3.txt
python 6_Visualization/scripts/plot_miami.py \
  --binary 3_MTAG/SI/EUR_MM/results/SI_EUR_MM_Output_trait_3.txt \
  --continuous 3_MTAG/SI_continuous/EUR_MM/results/SI_EUR_MM_Output_continuous_trait_3.txt
```

## Recommended Run Order

Start with:

1. `EUR_MM` + `Loneliness` (trait 3) for complete pipeline shakeout.
2. Then expand to trait 1 and trait 2.
3. Then replicate in `EUR_Male_MM` and `EUR_Female_MM`.

## Common Pitfalls

- Genome build mismatch (`hg19` vs `hg38`) between summary stats and references.
- Missing login-protected downloads (MSigDB, some summary-stat sources).
- Running LDSC scripts before generating `.sumstats.gz` via `--run-ldsc-munge`.
- Missing MetaXcan scripts/models on local system.
