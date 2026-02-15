# Loneliness (Trait 3) MTAG Analysis Summary

## Overview
This document summarizes the analysis of loneliness GWAS results from MTAG across different stratifications of the Social Isolation (SI) phenotype.

**Analysis Date:** November 4, 2025  
**Significance Threshold:** p < 5e-8

---

## Significant GWAS Hits Count

| Stratification | Sample Description | Significant SNPs | Total SNPs |
|----------------|-------------------|------------------|------------|
| **EUR_Female_MM** | European Females | **7** | 1,205,618 |
| **EUR_Male_MM** | European Males | **1** | 1,205,404 |
| **EUR_MM** | Combined European | **152** | 1,206,657 |

---

## Key Findings

### 1. Sex-Specific Differences
- **Females show more significant associations** than males (7 vs 1 significant SNPs)
- This suggests potential sex-specific genetic architecture for loneliness
- The male population shows notably fewer genome-wide significant hits

### 2. Combined Analysis Power
- The combined EUR_MM analysis identified **152 significant SNPs**
- This dramatic increase (compared to sex-stratified analyses) demonstrates:
  - Increased statistical power from larger sample size
  - Shared genetic architecture across sexes
  - Potentially sex-averaged effects that are consistent across populations

### 3. Statistical Considerations
- The MTAG analysis benefited from leveraging correlations between:
  - Ability to Confide (Trait 1)
  - Frequency of Social Contact (Trait 2)
  - Loneliness (Trait 3)
- All analyses included the `--force` flag due to some traits having mean chi² < 1.02

---

## Generated Visualizations

### Individual Manhattan Plots
1. **EUR_Female_MM**: `EUR_Female_MM/manhattan_loneliness_EUR_Female_MM.png`
2. **EUR_Male_MM**: `EUR_Male_MM/manhattan_loneliness_EUR_Male_MM.png`
3. **EUR_MM**: `EUR_MM/manhattan_loneliness_EUR_MM.png`

### Comparison Plot
- **Male vs Female Comparison**: `loneliness_male_vs_female_comparison.png`
  - Side-by-side Manhattan plots for direct comparison
  - Highlights sex-specific signal differences

---

## Technical Details

### MTAG Input Files
- Trait 3 corresponds to **Loneliness** in the social isolation phenotype battery
- Input summary statistics: `*.Day_NoPCs.mtag.sumstats.txt`
- MTAG leveraged genetic correlations between all three SI traits

### Output Files Analyzed
- `SI_EUR_Female_MM_Output_trait_3.txt`
- `SI_EUR_Male_MM_Output_trait_3.txt`
- `SI_EUR_MM_Output_trait_3.txt`

### Plot Specifications
- **Resolution**: 300 DPI
- **Significance threshold**: Horizontal line at -log10(5e-8)
- **Chromosome coloring**: Alternating colors for visual clarity
- **Statistics overlay**: SNP counts and total variants displayed

---

## Interpretation

### Biological Implications
1. **Sex-Specific Genetic Effects**: The disparity between male and female significant hits suggests:
   - Different genetic pathways may contribute to loneliness by sex
   - Hormonal or social factors may interact differently with genetic variants
   - Power differences due to sample size or phenotype measurement differences

2. **Shared Genetic Architecture**: The large number of hits in the combined analysis indicates:
   - Substantial overlap in genetic susceptibility across sexes
   - Common biological pathways underlying loneliness
   - Value of increasing sample size for discovery

3. **Phenotypic Correlation**: MTAG's effectiveness demonstrates that:
   - Loneliness shares genetic etiology with related social traits
   - Multi-trait analysis enhances discovery power
   - Social phenotypes form a genetically-related cluster

### Future Directions
- Investigate the 7 female-specific loci for sex-differential effects
- Fine-mapping of the 152 significant loci from combined analysis
- Gene-based and pathway analyses
- Genetic correlation with psychiatric and behavioral traits
- Polygenic risk score development and validation

---

## Files Generated

```
MTAG_Analysis/SI/
├── EUR_Female_MM/
│   └── manhattan_loneliness_EUR_Female_MM.png (1.5M)
├── EUR_Male_MM/
│   └── manhattan_loneliness_EUR_Male_MM.png (1.5M)
├── EUR_MM/
│   └── manhattan_loneliness_EUR_MM.png (974K)
├── loneliness_male_vs_female_comparison.png (1.9M)
├── generate_manhattan_plots.py
└── LONELINESS_ANALYSIS_SUMMARY.md
```

---

## Methods Summary

**MTAG Version**: 1.0.8  
**LD Reference Panel**: EUR (European) from 1000 Genomes  
**Analysis Parameters**:
- Minimum N: 0.0
- Force flag: Enabled (to handle low mean chi² in some traits)
- Stream stdout: Enabled

**Visualization**:
- Python libraries: pandas, matplotlib, seaborn, numpy
- Plot type: Manhattan plot (-log10 p-value vs genomic position)
- Color scheme: Chromosome-alternating for clarity

---

*Analysis conducted as part of the TsaiKellis Lab UKBB GWAS project*



