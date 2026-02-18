# 6_Visualization: Integrated Figures and Evidence Synthesis

This module turns the outputs of LDSC, MAGMA, fine-mapping, colocalization, and transcriptome analyses into interpretable visual summaries and a consolidated gene-prioritization table.

## Biological significance

Each upstream method captures a different biological dimension:

- genome-wide architecture (LDSC),
- gene/pathway/tissue enrichment (MAGMA),
- expression-mediated signal (S-PrediXcan),
- variant-level causality support (SuSiE),
- shared causal mechanisms with eQTL (coloc).

Visualization is where these dimensions are integrated into coherent biological narratives (for example, "trait X shows high genetic overlap with mood phenotypes, neuronal cell-type enrichment, and convergent gene-level evidence at locus Y").

## What this module produces

- Correlation structure heatmaps
- Top annotation/cell-type/tissue enrichment barplots
- Regional locus plots with optional PIP overlays
- Miami plots comparing binary vs continuous coding
- A multi-source gene-prioritization table

## Scripts and interpretation

### 1) `scripts/plot_rg_heatmap.py`

- Parses LDSC `rg` result files (`*.results` or `*.log`)
- Builds trait-by-trait matrix with `rg` values
- Output: heatmap (`-1` to `+1` color scale)

Interpretation:

- Red (positive) and blue (negative) correlation patterns reveal shared vs opposed genetic architecture.

### 2) `scripts/plot_partitioned_h2.py`

- Reads partitioned heritability results
- Ranks categories by `-log10(p)` (or enrichment fallback)
- Plots top categories

Interpretation:

- Highlights annotation classes disproportionately contributing to trait heritability.

### 3) `scripts/plot_celltype_enrichment.py`

- Parses LDSC-CTS outputs
- Plots top cell-type enrichments by `-log10(p)`

Interpretation:

- Suggests tissues/cell types where trait-associated regulatory variation is concentrated.

### 4) `scripts/plot_magma_tissue_enrichment.py`

- Parses MAGMA tissue-expression (`gene-property`) outputs
- Plots strongest tissue associations

Interpretation:

- Provides expression-context clues for implicated genes.

### 5) `scripts/build_gene_prioritization_table.py`

Merges evidence from:

- S-PrediXcan prioritized genes
- coloc summary (`coloc_max_pp_h4`)
- MAGMA best gene p-values
- optional fine-mapping support (`fine_map_max_pip`)

Builds binary evidence flags and an integrated score:

\[
\text{score} \approx \sum w_m \cdot [-\log_{10}(p_m)] + 5\cdot PP.H4 + 4\cdot PIP
\]

with method-specific weights in script defaults.

Output:

- `output/gene_prioritization_table.tsv`

### 6) `scripts/plot_regional_loci.py`

- Plots per-locus regional association (`-log10(p)` vs position)
- Optional SuSiE PIP overlay (point color/size)

Interpretation:

- Visualizes local signal shape and whether high-PIP variants coincide with strongest association peaks.

### 7) `scripts/plot_miami.py`

- Compares binary and continuous SI summary statistics in a single Miami plot
- Binary (top) vs continuous (bottom), with genome-wide significance lines

Interpretation:

- Quickly identifies shared vs coding-specific loci across phenotype definitions.

## Typical usage

From `5_Functional_Genomics/`:

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

## Inputs and output location

- Main input directories:
  - `1_LDSC/output/*`
  - `2_MAGMA/output/*`
  - `3_S-PrediXcan/output/*`
  - `4_Fine_Mapping/output/*`
  - `5_Colocalization/output/*`
- Main outputs:
  - `6_Visualization/output/*.png`
  - `6_Visualization/output/gene_prioritization_table.tsv`

## Interpretation caveats

- Visual prominence does not replace statistical correction.
- Ranking scores are heuristic prioritization tools, not formal posterior probabilities.
- Cross-method convergence is usually more informative than any single method.

## References and supporting resources

- Matplotlib: https://matplotlib.org/stable/
- Seaborn: https://seaborn.pydata.org/
- LocusZoom-style regional plotting concept: https://my.locuszoom.org/
- Miami/manhattan plot conventions overview: https://www.cog-genomics.org/plink/1.9/manhattan
