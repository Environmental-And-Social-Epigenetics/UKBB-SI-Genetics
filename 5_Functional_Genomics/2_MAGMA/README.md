# 2_MAGMA: Gene, Pathway, and Tissue Enrichment Analysis

This module translates SNP-level GWAS/MTAG signal into gene-level and pathway-level biology using MAGMA. It asks: which genes, biological pathways, and tissues are most strongly implicated by the association architecture?

## Biological significance

Single SNP hits are often hard to interpret biologically. MAGMA addresses this by:

- aggregating nearby SNP evidence into **gene-level association scores**,
- testing whether known **gene sets/pathways** are enriched,
- testing whether genes highly expressed in particular **tissues** show stronger association.

For social isolation phenotypes, this helps move from loci to putative mechanisms (for example, neurobiological, immune, or endocrine pathways).

## Statistical framework

MAGMA implements a multiple regression / SNP-wise mean framework that combines SNP association statistics while accounting for LD between SNPs.

Conceptually:

1. Map SNPs to genes (here using symmetric \(\pm 10\) kb window).
2. Build a gene-level test statistic from SNP p-values and LD structure.
3. Test whether genes in a pathway (or genes with high expression in a tissue) have stronger gene-level statistics than background genes.

### Gene-level test

MAGMA's gene analysis aggregates SNP evidence while correcting for LD, avoiding naive independence assumptions.

### Gene-set analysis

This pipeline uses competitive testing: are genes in set \(S\) more strongly associated than genes outside \(S\), conditional on gene-level signal distribution?

### Tissue-expression analysis

`--gene-covar` and `--model direction=pos condition-hide=Average` test whether increased tissue-specific expression predicts stronger association signal (while conditioning on average expression).

## Scripts and pipeline order

Run from `5_Functional_Genomics/`:

```bash
bash 2_MAGMA/scripts/01_annotate.sh
bash 2_MAGMA/scripts/02_gene_analysis.sh
bash 2_MAGMA/scripts/03_geneset_analysis.sh
bash 2_MAGMA/scripts/04_tissue_expression.sh
```

### `scripts/01_annotate.sh`

- Input: `0_munge/output/magma/*.genes.raw`
- Action: SNP-to-gene annotation using `NCBI37.3.gene.loc` with `window=10,10` kb
- Output: `output/01_annotation/*.genes.annot`

### `scripts/02_gene_analysis.sh`

- Input: `.genes.annot` plus SNP p-value files
- Reference LD: `reference_data/1000G_EUR/g1000_eur`
- Output: `output/02_gene_analysis/*.genes.out` (+ MAGMA logs/raw intermediates)

### `scripts/03_geneset_analysis.sh`

- Input: gene-level results and gene-set files (`reference_data/magma/msigdb/*.gmt`)
- Output: `output/03_geneset_analysis/*`

### `scripts/04_tissue_expression.sh`

- Input: gene-level results + GTEx tissue covariates
- Covariate file default: `reference_data/magma/gtex_v8/GTEx_v8_tissue_expression.gene_covar.txt`
- Output: `output/04_tissue_expression/*`

## Input and output expectations

### Inputs

- SNP-level MAGMA table from munging (`SNP CHR BP P N`)
- Gene location annotation file (`NCBI37.3.gene.loc`)
- EUR LD reference PLINK files (`.bed/.bim/.fam`)
- Optional/extended:
  - MSigDB `.gmt` pathway files
  - GTEx v8 tissue expression covariate table

### Outputs

- **Gene annotation**: `.genes.annot`
- **Gene associations**: `.genes.out` (main gene p-values and z-statistics)
- **Pathway enrichment**: gene-set outputs per trait and set collection
- **Tissue enrichment**: gene-property outputs per trait

## Interpretation guidance

- Gene significance should be interpreted with multiple-testing correction at the gene level.
- Pathway enrichment depends on gene-set definition quality and overlap among sets.
- Tissue enrichment is suggestive of regulatory context, not proof of causal tissue.
- Compare concordance with S-PrediXcan and coloc outputs for stronger mechanistic confidence.

## References and documentation

- MAGMA software/documentation: https://ctg.cncr.nl/software/magma
- MAGMA paper (de Leeuw et al., 2015): https://doi.org/10.1371/journal.pcbi.1004219
- MSigDB database: https://www.gsea-msigdb.org/gsea/msigdb
- GTEx portal: https://gtexportal.org/home/
