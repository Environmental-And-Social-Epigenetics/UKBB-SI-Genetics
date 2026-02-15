# BOLT-LMM GWAS Analysis: Social Isolation and Loneliness Phenotypes

**BOLT-LMM analysis for 3 social isolation/loneliness phenotypes across 3 population stratifications**

---

## Overview

This analysis performs genome-wide association studies (GWAS) for social isolation and loneliness phenotypes using BOLT-LMM v2.5 with mixed linear models. The analysis is stratified by sex to detect sex-specific genetic effects on social behavior traits.

### Phenotypes (3 total)

| Phenotype | Description | Type | Prevalence |
|-----------|-------------|------|------------|
| **Loneliness** | Self-reported feelings of loneliness | Binary | ~5-10% |
| **FreqSoc** | Frequency of social contact | Binary | Variable |
| **AbilityToConfide** | Having someone to confide in | Binary | ~90% yes |

**Source**: Based on Day et al. (2018) phenotype definitions for UK Biobank social isolation traits.

### Population Stratifications (3 total)

| Population | Description | Expected N |
|------------|-------------|------------|
| **EUR_MM** | European ancestry (includes related individuals) | ~426,000 |
| **EUR_Male** | European ancestry males only | ~200,000 |
| **EUR_Female** | European ancestry females only | ~226,000 |

### Covariate Model

**Day_NoPCs** (Primary model without principal components):
- Age (quantitative)
- Sex (categorical) - except in sex-stratified analyses
- Genotyping array (categorical)
- Population structure controlled via genetic relationship matrix (GRM)

**Note**: No principal components included because:
1. GRM captures population structure efficiently
2. Maximizes degrees of freedom for better power
3. Consistent with Day et al. (2018) primary analysis model

---

## Analysis Configuration

### Total Jobs

```
3 phenotypes × 3 populations × 1 covariate set = 9 jobs
```

### Job Mapping (SLURM array tasks 1-9)

| Task | Phenotype | Population |
|------|-----------|------------|
| 1-3 | Loneliness | EUR_MM, EUR_Male, EUR_Female |
| 4-6 | FreqSoc | EUR_MM, EUR_Male, EUR_Female |
| 7-9 | AbilityToConfide | EUR_MM, EUR_Male, EUR_Female |

---

## Quick Start

### Prerequisites

On the HPC, ensure:
- [x] Genotype files converted to BED format (`ukb_genoHM3_bed.*`)
- [x] Model SNPs prepared (`ukb_genoHM3_modelSNPs.txt`)
- [x] Population keep files exist in `sqc/population.20220316/`
- [x] Phenotype file: `pheno/isolation_run_control.tsv.gz`
- [x] Covariate file: `sqc/sqc.20220316.tsv.gz`

### Step 1: Filter to Populations

```bash
# On HPC, in this directory
cd /home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness

# Submit batch job to filter phenotype and covariate files for all 3 populations
sbatch 0a_filter_populations.sbatch.sh

# Monitor
tail -f 0a_filter.out

# This creates:
# - isolation_run_control.EUR_MM.tsv.gz
# - isolation_run_control.EUR_Male.tsv.gz
# - isolation_run_control.EUR_Female.tsv.gz
# - sqc.EUR_MM.tsv.gz
# - sqc.EUR_Male.tsv.gz
# - sqc.EUR_Female.tsv.gz
```

### Step 2: Test Run

```bash
# Test with one phenotype and population
sbatch 0b_test_run.sbatch.sh

# Monitor
tail -f 0b_test.out

# Check for "TEST PASSED" message
# Review output: results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.stats.gz
```

### Step 3: Full Analysis

```bash
# If test passes, submit all 9 jobs
sbatch 1_run_bolt_lmm.sbatch.sh

# Monitor progress
squeue -u mabdel03

# Check individual job outputs
tail -f 1_*.out
```

---

## Output Structure

```
BOLT-LMM_SI-Loneliness/
├── results/
│   └── Day_NoPCs/
│       ├── EUR_MM/
│       │   ├── bolt_Loneliness.Day_NoPCs.stats.gz
│       │   ├── bolt_Loneliness.Day_NoPCs.log.gz
│       │   ├── bolt_FreqSoc.Day_NoPCs.stats.gz
│       │   ├── bolt_FreqSoc.Day_NoPCs.log.gz
│       │   ├── bolt_AbilityToConfide.Day_NoPCs.stats.gz
│       │   └── bolt_AbilityToConfide.Day_NoPCs.log.gz
│       ├── EUR_Male/
│       │   └── [same 6 files]
│       └── EUR_Female/
│           └── [same 6 files]
├── isolation_run_control.EUR_MM.tsv.gz
├── isolation_run_control.EUR_Male.tsv.gz
├── isolation_run_control.EUR_Female.tsv.gz
├── sqc.EUR_MM.tsv.gz
├── sqc.EUR_Male.tsv.gz
└── sqc.EUR_Female.tsv.gz
```

---

## Expected Runtime

Per job (150GB RAM, 100 CPUs):
- Sample size: ~200K-426K individuals
- Variants: ~1.3M autosomal (HapMap3)
- Expected time: **1-2 hours per job**

Total wall-clock time (if run concurrently): **~2 hours**

---

## Quality Control

### From Log Files

For each phenotype-population combination, check:

```bash
# Sample size
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log.gz | grep "Analyzing"

# Heritability estimate
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log.gz | grep "h2:"

# Genomic inflation
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log.gz | grep -i "lambda\|inflation"

# Warnings
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log.gz | grep -i "warning"
```

**Expected QC metrics:**
- λ_GC: 1.00-1.05 (well-calibrated)
- h²: 3-15% for social isolation traits (Day et al. 2018)
- Sample size: matches expected for population

### From Summary Statistics

```bash
# Count variants
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.stats.gz | wc -l
# Expected: ~1.3M

# Check for genome-wide significant hits (p < 5×10⁻⁸)
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.stats.gz | \
    awk 'NR>1 && $NF < 5e-8' | wc -l

# Preview top associations
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.stats.gz | \
    awk 'NR>1' | sort -k12,12g | head -20
```

---

## Understanding Binary Trait Results

### Liability Threshold Model

BOLT-LMM uses a liability threshold model for binary traits:

- **Effect size (BETA)**: On the liability scale (continuous latent variable)
- **Interpretation**: Each copy of effect allele changes liability by BETA standard deviations
- **Conversion to odds ratio**: OR ≈ exp(BETA) for small effects
- **Example**: BETA=0.05 → OR ≈ 1.051 (5.1% increased odds per allele)

### Phenotype-Specific Considerations

**Loneliness** (Case prevalence ~5-10%):
- Lower prevalence = higher effective sample size penalty
- Look for associations in psychological/neurological pathways

**AbilityToConfide** (Case prevalence ~90%):
- High prevalence trait (most people have someone to confide in)
- Cases are actually "0" (no one to confide in) - minority group
- Consider inverting for interpretation

**FreqSoc** (Frequency of social contact):
- More balanced prevalence
- May show associations with social behavior and personality genes

---

## Sex-Stratified Analysis Interpretation

### Comparing Results Across Populations

**1. Identify sex-specific associations:**

Compare male vs female results to find variants that show:
- Genome-wide significance in one sex but not the other
- Substantially different effect sizes between sexes
- Potential sex-by-SNP interactions

**2. Known sex differences in social behavior:**

- Females: Generally report more social connectedness
- Males: Higher rates of social isolation in some studies
- Social support networks differ by sex
- Expect some sex-specific genetic effects

### Biological Interpretation

**Potential genetic pathways:**
- **Neurotransmitter systems**: Serotonin, dopamine, oxytocin
- **Personality traits**: Extraversion, neuroticism
- **Mental health**: Depression, anxiety (genetic correlation)
- **Social cognition**: Empathy, theory of mind
- **Stress response**: HPA axis function

---

## Comparison to Day et al. (2018)

### Methodological Alignment

This analysis replicates the Day et al. (2018) UK Biobank analysis with additional sex stratification:

| Aspect | Day et al. (2018) | This Analysis |
|--------|-------------------|---------------|
| Software | BOLT-LMM v2.3 | BOLT-LMM v2.5 |
| Phenotypes | Same definitions | Same |
| Sample size | ~456,000 EUR | ~426,000 EUR |
| Covariates | Age, sex, array, center, PCs | Age, sex, array (GRM) |
| **NEW** | Combined sexes | **Sex-stratified** |

### Expected Findings

Based on Day et al. (2018), we expect:

**Loneliness**:
- h² ≈ 3-5% on liability scale
- Multiple genome-wide significant loci
- Genetic correlation with depression, neuroticism

**Frequency of Social Contact**:
- h² ≈ 4-6%
- Associations with personality and behavior genes
- Genetic correlation with educational attainment

**Ability to Confide**:
- h² ≈ 3-5%
- Associations with social support and mental health
- Potential overlap with loneliness loci

---

## Downstream Analyses

### 1. LD Score Regression

Estimate heritability and genetic correlations:

```bash
# Heritability for each phenotype-population
for pop in EUR_MM EUR_Male EUR_Female; do
    for pheno in Loneliness FreqSoc AbilityToConfide; do
        ldsc.py \
            --h2 results/Day_NoPCs/${pop}/bolt_${pheno}.Day_NoPCs.stats.gz \
            --ref-ld-chr eur_w_ld_chr/ \
            --w-ld-chr eur_w_ld_chr/ \
            --out ${pheno}_${pop}.h2
    done
done

# Genetic correlation between males and females
ldsc.py \
    --rg results/Day_NoPCs/EUR_Male/bolt_Loneliness.Day_NoPCs.stats.gz,results/Day_NoPCs/EUR_Female/bolt_Loneliness.Day_NoPCs.stats.gz \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out Loneliness_Male_vs_Female.rg
```

### 2. Cross-Trait Genetic Correlations

Compare with other traits:
- Mental health: Depression, anxiety, bipolar disorder
- Personality: Big Five traits (especially extraversion, neuroticism)
- Cognitive: Educational attainment, intelligence
- Physical health: BMI, cardiovascular disease

### 3. Polygenic Risk Scores

```bash
# Build PRS using PRSice-2 or LDpred2
# Can predict loneliness in independent cohorts
# Assess sex-specific PRS performance
```

### 4. Meta-Analysis

Combine male and female results:
```bash
# Use METAL or similar
# Weight by sample size or inverse variance
# Test for heterogeneity (sex differences)
```

---

## Troubleshooting

### Common Issues

**Issue: Population-filtered files not found**
```bash
# Solution: Run filter script
sbatch 0a_filter_populations.sbatch.sh
```

**Issue: Low case/control counts**
```bash
# Check phenotype distributions
zcat isolation_run_control.EUR_MM.tsv.gz | awk -F'\t' '{print $3}' | sort | uniq -c

# For binary traits, need sufficient cases and controls (>100 each ideally)
```

**Issue: High λ_GC (>1.10)**
```bash
# Check: Population stratification control adequate?
# Check: Any batch effects in phenotypes?
# Consider: Using Day_10PCs covariate set instead
```

---

## File Descriptions

| File | Purpose |
|------|---------|
| `run_single_phenotype.sh` | Core BOLT-LMM execution for one phenotype-population combo |
| `1_run_bolt_lmm.sbatch.sh` | SLURM array job (9 tasks) |
| `0a_filter_populations.sbatch.sh` | Filter all 3 populations |
| `0b_test_run.sbatch.sh` | Test script for validation |
| `filter_to_population.sh` | Population-specific filtering (called by 0a) |
| `paths.sh` | Configuration for file paths and BOLT-LMM settings |

---

## Resources

- **BOLT-LMM Manual**: https://alkesgroup.broadinstitute.org/BOLT-LMM/
- **Day et al. (2018)**: Nature Communications (methodology reference)
- **UK Biobank**: https://biobank.ndph.ox.ac.uk/

---

## Citation

If you use this analysis pipeline, please cite:

1. **BOLT-LMM**:
   - Loh, P.-R., et al. (2015). Nature Genetics, 47(3), 284-290.
   - Loh, P.-R., et al. (2018). Nature Genetics, 50(7), 906-908.

2. **UK Biobank**:
   - Bycroft, C., et al. (2018). Nature, 562(7726), 203-209.

3. **Phenotype definitions and methodology**:
   - Day, F. R., et al. (2018). Nature Communications, 9(1), 2457.

---

*Analysis Directory: `/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness`*

*Last Updated: October 30, 2025*
