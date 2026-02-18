# 7_Sex_Differences: Sex-Stratified Genetic Architecture of SI Traits

This module performs a multi-level comparison of male and female genetic architecture for social isolation traits (binary and continuous codings), from genome-wide architecture down to genes and pathways.

## Biological significance

Sex differences in genetic architecture can arise from:

- sex-specific regulatory environments,
- hormone-sensitive pathways,
- sex-dependent developmental trajectories,
- differential interaction with environmental exposures.

This pipeline is designed to detect where those differences appear:

- globally (`h2`, cross-sex `rg`),
- at SNP/locus level (interaction and locus class),
- at gene/pathway level (MAGMA differential analysis and pathway enrichment).

## End-to-end analysis layers

1. **Cross-sex genetic correlation** (`rg`)
2. **Sex comparison of SNP heritability** (`h2_male` vs `h2_female`)
3. **SNP-level sex interaction**
4. **Locus-level sex-specificity classification**
5. **Sex-stratified MAGMA gene analysis**
6. **Sex-differential genes**
7. **Sex-differential pathway enrichment**
8. **Integrated visualization outputs**

---

## Statistical details by step

### 1) Cross-sex LDSC `rg`

`scripts/01_cross_sex_rg.sh` runs LDSC `--rg` for male vs female summary stats for each trait/coding.

Interpretation:

- `rg ~ 1`: highly similar additive architecture across sexes
- `rg < 1`: partial divergence
- lower/negative values: stronger sex heterogeneity

### 2) Heritability difference test

`scripts/02_h2_sex_comparison.py` parses LDSC logs and computes:

\[
Z_{\Delta h^2} = \frac{h^2_{\text{female}} - h^2_{\text{male}}}{\sqrt{SE_{\text{female}}^2 + SE_{\text{male}}^2}}
\]

Two-sided p-values are then used to label significance.

### 3) SNP-level sex interaction

`scripts/03_snp_sex_interaction.py` harmonizes male/female alleles and computes:

\[
Z_{\text{int}} = \frac{\beta_{\text{male}} - \beta_{\text{female}}}{\sqrt{SE_{\text{male}}^2 + SE_{\text{female}}^2}}
\]

\[
p_{\text{int}} = 2\Phi(-|Z_{\text{int}}|)
\]

Outputs include genome-wide interaction tables and top hits per trait.

### 4) Locus-level sex classification

`scripts/04_locus_sex_classification.py`:

- defines male/female/combined loci via `4_Fine_Mapping/scripts/01_define_loci.py`,
- merges nearby leads,
- classifies each locus using male/female p-value patterns and direction concordance.

Categories include:

- `shared`
- `female-specific` / `male-specific`
- `female-amplified` / `male-amplified`
- `discordant`
- `suggestive-other`

### 5) Sex-stratified MAGMA

`scripts/05_magma_sex_stratified.sh` runs MAGMA annotation + gene analysis separately for male and female datasets (both codings, all traits).

### 6) Sex-differential genes

`scripts/06_sex_differential_genes.py` merges male/female MAGMA gene outputs and computes:

\[
Z_{\text{diff}} = \frac{Z_{\text{male}} - Z_{\text{female}}}{\sqrt{2}}
\]

with p-values, BH FDR, and category labels:

- `shared_significant`
- `male_specific`, `female_specific`
- `male_enriched`, `female_enriched`
- `sex_differential`

### 7) Sex-differential pathways

`scripts/07_sex_differential_pathways.py` performs rank-based pathway testing on `|Z_diff|` using Mann-Whitney U tests (per coding/trait), then BH correction.

Curated sets include:

- hormone-related genes
- HPA-axis stress genes
- serotonin / dopamine signaling
- oxytocin-vasopressin
- GABA/glutamate signaling
- microglia/immune
- X-linked genes

### 8) Integrated figures

`scripts/08_sex_difference_plots.py` generates:

- cross-sex effect scatter plots
- sex-interaction Manhattan plots
- male/female Miami plots
- h2 comparison barplots
- cross-sex `rg` forest plot
- locus category stacked bars
- sex-differential gene volcano
- pathway sex-difference heatmap

---

## Inputs and outputs

### Inputs

- Sex-stratified MTAG outputs in:
  - `3_MTAG/SI/EUR_Male_MM/results/*`
  - `3_MTAG/SI/EUR_Female_MM/results/*`
  - `3_MTAG/SI_continuous/EUR_Male_MM/results/*`
  - `3_MTAG/SI_continuous/EUR_Female_MM/results/*`
- Munged LDSC and MAGMA files from upstream modules
- Reference resources used by LDSC/MAGMA/fine-mapping modules

### Key output folders

- `output/01_cross_sex_rg/`
- `output/02_h2_sex_comparison/`
- `output/03_snp_sex_interaction/`
- `output/04_locus_classification/`
- `output/05_magma_sex_stratified/`
- `output/06_sex_differential_genes/`
- `output/07_sex_differential_pathways/`
- `output/08_sex_difference_plots/`

## Recommended run order

From `5_Functional_Genomics/`:

```bash
bash 7_Sex_Differences/scripts/01_cross_sex_rg.sh
python 7_Sex_Differences/scripts/02_h2_sex_comparison.py
python 7_Sex_Differences/scripts/03_snp_sex_interaction.py
python 7_Sex_Differences/scripts/04_locus_sex_classification.py
bash 7_Sex_Differences/scripts/05_magma_sex_stratified.sh
python 7_Sex_Differences/scripts/06_sex_differential_genes.py
python 7_Sex_Differences/scripts/07_sex_differential_pathways.py
python 7_Sex_Differences/scripts/08_sex_difference_plots.py
```

## Interpretation caveats

- Sex-stratified analyses have reduced effective sample size per stratum, lowering power.
- Null results can reflect power limits, not absence of biological differences.
- Direction/sign differences should be confirmed with strict allele harmonization and external replication when possible.
- Pathway differences are best treated as hypothesis-generating unless convergent with independent evidence.

## References and supporting resources

- Rawlik et al. (2016), sex-specific genetic architecture concepts: https://doi.org/10.1038/ncomms11483
- Khramtsova et al. (2019), sex differences in psychiatric genetics review: https://doi.org/10.1038/s41380-018-0135-6
- Ngun et al. (2011), sex differences in brain/gene expression context: https://doi.org/10.1177/1073858410386494
- LDSC documentation: https://github.com/bulik/ldsc
- MAGMA documentation: https://ctg.cncr.nl/software/magma
