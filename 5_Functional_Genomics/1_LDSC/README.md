# 1_LDSC: Heritability and Genetic Correlation Architecture

This module runs LD Score Regression (LDSC) analyses on munged GWAS/MTAG summary statistics to answer four major questions:

1. How polygenic is each trait (SNP heritability, `h2`)?
2. How much genetic signal is shared across traits (`rg`)?
3. Which functional annotations carry disproportionate heritability (partitioned `h2`)?
4. Which cell-type-relevant annotations are enriched (LDSC-CTS)?

## Biological significance

LDSC connects GWAS summary statistics to biological interpretation at the architecture level:

- **Trait-level architecture**: whether trait signal is diffuse (polygenic) or sparse.
- **Cross-trait overlap**: whether two traits likely share upstream genetic mechanisms.
- **Annotation enrichment**: whether specific regulatory genomic categories contribute more than expected by SNP count.
- **Cell-type relevance**: whether disease or behavioral genetics point to particular tissue/cell programs.

For social isolation traits, this helps separate:

- global polygenicity from confounding,
- shared psychiatric/neurological liability from trait-specific effects,
- broad vs cell-type-specific biological pathways.

## Statistical framework (high level)

### 1) Univariate LDSC (`--h2`)

For SNP \(j\), LDSC regresses GWAS test statistics on LD score \(l_j\):

\[
E[\chi_j^2] \approx 1 + \frac{N h^2}{M} l_j + c
\]

- \(N\): sample size
- \(M\): number of SNPs
- \(h^2\): SNP heritability
- \(c\): intercept capturing inflation (population stratification, cryptic relatedness, etc.)

Key interpretation:

- Higher slope -> larger SNP heritability.
- Intercept separates true polygenicity from confounding.

### 2) Bivariate LDSC (`--rg`)

Uses cross-trait summary statistics to estimate genetic covariance/correlation:

\[
rg = \frac{\mathrm{cov}_g(T_1, T_2)}{\sqrt{h_1^2 h_2^2}}
\]

Interpretation:

- `rg ~ 1`: nearly same additive genetic architecture
- `rg ~ 0`: little shared additive signal
- `rg < 0`: opposite signed genetic effects on average

### 3) Partitioned heritability

Uses annotation-stratified LD scores:

\[
E[\chi_j^2] \approx 1 + N \sum_c \tau_c \, l_{j,c}
\]

where \(\tau_c\) is annotation-specific contribution. Enrichment asks whether annotation \(c\) explains more heritability than expected from its SNP proportion.

### 4) Cell-type-specific LDSC (LDSC-CTS)

Tests whether trait heritability is enriched in annotations linked to cell-type expression programs (while conditioning on baseline LD annotations).

## Scripts and what each one does

- `scripts/run_h2.sh`
  - Runs univariate LDSC heritability for each `.sumstats.gz`.
  - Output: `output/h2/*.log`

- `scripts/run_rg_internal.sh`
  - Creates pairwise trait combinations within each analysis group/population.
  - Runs LDSC `--rg` for each pair.
  - Output: `output/rg_internal/*.log` and `_internal_pairs.tsv`

- `scripts/run_rg_external.sh`
  - Computes `rg` between local SI traits and external GWAS listed in config.
  - Output: `output/rg_external/*.log`

- `scripts/run_partitioned_h2.sh`
  - Runs baselineLD v2.2 partitioned heritability.
  - Output: `output/partitioned_h2/*.results`, `*.log`

- `scripts/run_celltype_h2.sh`
  - Runs LDSC-CTS using annotation groups from `config/cell_type_groups.ldcts`.
  - Output: `output/celltype_h2/*.results`, `*.log`

## Configuration files

- `config/external_traits.tsv`
  - Manifest of external traits for cross-trait `rg`.
  - You must fill `sumstats_path` with local paths to pre-munged `.sumstats.gz`.

- `config/cell_type_groups.ldcts`
  - LDSC-CTS annotation list (`<cell_type_name><TAB><annotation_prefix>`).
  - Replace placeholder prefixes with real annotation paths.

## Inputs and outputs

### Required inputs

- Munged files from `0_munge/output/ldsc_munged/*.sumstats.gz`
- LDSC software: `reference_data/ldsc/ldsc/ldsc.py`
- LD references (EUR):
  - `reference_data/ldsc/eur_w_ld_chr/`
  - `reference_data/ldsc/weights_hm3_no_hla/weights.hm3_noMHC.`
- For partitioned `h2`: baselineLD v2.2 files
- For CTS: `.ldcts` annotation mapping

### Main outputs

- `output/h2/`: per-trait heritability logs
- `output/rg_internal/`: within-cohort trait-pair genetic correlations
- `output/rg_external/`: local-vs-external trait genetic correlations
- `output/partitioned_h2/`: annotation enrichment statistics
- `output/celltype_h2/`: cell-type enrichment statistics

## Example usage

From `5_Functional_Genomics/`:

```bash
bash 1_LDSC/scripts/run_h2.sh
bash 1_LDSC/scripts/run_rg_internal.sh
bash 1_LDSC/scripts/run_rg_external.sh
bash 1_LDSC/scripts/run_partitioned_h2.sh
bash 1_LDSC/scripts/run_celltype_h2.sh
```

## Interpretation tips

- Always check LDSC intercept and ratio in logs before over-interpreting `h2`.
- For many trait pairs, correct for multiple testing when interpreting `rg`.
- Partitioned/cell-type enrichments are annotation-dependent; interpret in the context of baseline model and annotation quality.
- Compare binary and continuous SI codings to test robustness of biological signals.

## References and documentation

- LDSC repository/documentation: https://github.com/bulik/ldsc
- Bulik-Sullivan et al. (2015), LD Score regression: https://www.nature.com/articles/ng.3211
- Bulik-Sullivan et al. (2015), atlas of genetic correlations: https://www.nature.com/articles/ng.3406
- Finucane et al. (2015), partitioned heritability: https://www.nature.com/articles/ng.3404
- Finucane et al. (2018), specifically expressed genes / CTS framework: https://www.nature.com/articles/s41588-018-0081-4
