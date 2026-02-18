# 4_Fine_Mapping: From GWAS Loci to Credible Causal Variant Sets

This module narrows broad GWAS signals into plausible causal SNP subsets by combining locus definition, local LD, and Bayesian fine-mapping (SuSiE).

## Biological significance

Genome-wide significant loci usually contain many correlated SNPs due to linkage disequilibrium (LD). Fine-mapping improves biological interpretability by estimating which variants are most likely causal, which supports:

- variant-to-gene linking,
- functional follow-up experiments,
- integration with colocalization and molecular QTL data.

Without fine-mapping, lead SNPs are often just proxies.

## Statistical framework

## 1) Locus definition

`scripts/01_define_loci.py` identifies independent lead loci from significant SNPs (`p < 5e-8` by default) using:

- **Preferred**: PLINK clumping (`--clump-r2`, `--clump-kb`, `--clump-p1/p2`)
- **Fallback**: distance-based pruning within `window_kb`

Each lead SNP defines a region:

\[
[\text{lead BP} - \text{window}, \ \text{lead BP} + \text{window}]
\]

default `window_kb = 500`.

## 2) Locus-specific LD matrices

`scripts/02_compute_ld.sh` uses a reference panel (`1000G_EUR`) to compute per-locus LD matrices (`R`) via PLINK (`--r square gz`). These matrices encode correlation structure needed by summary-statistic fine-mapping.

## 3) SuSiE fine-mapping (`susie_rss`)

`scripts/03_run_susie.R` runs `susieR::susie_rss` on each locus:

- Inputs: z-scores (`beta/se`), LD matrix `R`, effective sample size `n`
- Model: sum of single-effect components (`L`, default 10)
- Outputs:
  - **PIP** (posterior inclusion probability) per SNP
  - **Credible sets** (default 95% coverage)

Interpretation:

- Higher PIP -> higher posterior probability SNP is causal (within model assumptions).
- A credible set is a compact SNP set expected to contain the causal variant with specified posterior coverage.

## Scripts and expected outputs

### `scripts/01_define_loci.py`

Input:

- GWAS/MTAG summary stats (`--sumstats`)

Output prefix `<out-prefix>` creates:

- `<out-prefix>.sig_snps.tsv`
- `<out-prefix>.lead_snps.tsv`
- `<out-prefix>.loci.tsv`

### `scripts/02_compute_ld.sh`

Input:

- Locus table (default `output/01_defined_loci/loci.loci.tsv`)
- PLINK reference panel (`reference_data/1000G_EUR/g1000_eur`)

Output per locus:

- `<locus_id>.snplist`
- `<locus_id>.ld.gz` (square LD matrix)

### `scripts/03_run_susie.R`

Input:

- `--sumstats`, `--loci`, `--ld-dir`

Output:

- `susie_pip.tsv`
- `susie_credible_sets.tsv`
- `susie_locus_summary.tsv`

## Typical run sequence

```bash
python 4_Fine_Mapping/scripts/01_define_loci.py \
  --sumstats 3_MTAG/SI/EUR_MM/results/SI_EUR_MM_Output_trait_3.txt \
  --plink-bfile 5_Functional_Genomics/reference_data/1000G_EUR/g1000_eur \
  --allow-distance-fallback

bash 4_Fine_Mapping/scripts/02_compute_ld.sh

Rscript 4_Fine_Mapping/scripts/03_run_susie.R \
  --sumstats 3_MTAG/SI/EUR_MM/results/SI_EUR_MM_Output_trait_3.txt \
  --loci 4_Fine_Mapping/output/01_defined_loci/loci.loci.tsv \
  --ld-dir 4_Fine_Mapping/output/02_ld_matrices \
  --out-dir 4_Fine_Mapping/output/03_susie
```

## Practical interpretation notes

- PIP values are model-based probabilities, not frequentist p-values.
- Credible sets depend strongly on LD reference ancestry match.
- If loci are very large/high-LD, credible sets may remain broad.
- Integrate fine-mapping with colocalization and TWAS evidence before selecting top biological targets.

## References and documentation

- SuSiE R package: https://stephenslab.github.io/susieR/
- Wang et al. (2020), SuSiE method: https://doi.org/10.1111/rssb.12388
- PLINK documentation: https://www.cog-genomics.org/plink/
- Fine-mapping review (Schaid et al., 2018): https://doi.org/10.1038/s41576-018-0018-x
