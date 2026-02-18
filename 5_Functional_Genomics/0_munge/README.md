# 0_munge: Summary Statistics Harmonization for Functional Genomics

This step converts raw MTAG outputs into clean, tool-specific inputs for LDSC, MAGMA, and S-PrediXcan. Its main purpose is to make all downstream analyses comparable and biologically interpretable by enforcing consistent variant identifiers, allele orientation, sample size fields, and valid p-value/frequency ranges.

## Why this analysis matters biologically

Post-GWAS interpretation methods answer different biological questions (heritability architecture, gene-level signal, tissue expression mediation, etc.), but they all rely on the same core SNP-level summary statistics. If those statistics are not harmonized first, downstream biological conclusions can be biased by technical artifacts (for example, allele mislabeling, duplicated SNP IDs, or invalid frequencies), not true biology.

In practice, this step is quality control plus format translation:

- Quality control protects against spurious enrichment or false-positive gene hits.
- Format translation lets one GWAS result be reused consistently across LDSC, MAGMA, and transcriptome-based methods.

## Core statistical / data-processing logic

`scripts/munge_sumstats.py` performs:

1. **Column validation**: requires `SNP`, `CHR`, `BP`, `A1`, `A2`, `N`, `FRQ`, `mtag_beta`, `mtag_se`, `mtag_z`, `mtag_pval`.
2. **Numeric coercion and validity filters**:
   - `0 < p <= 1` (with clipping to floating-point epsilon for numerical stability)
   - `0 < FRQ < 1`
   - optional MAF filter: `min(FRQ, 1-FRQ) >= --min-maf`
3. **Duplicate handling**: if a SNP appears multiple times, keep the row with smallest p-value.
4. **Stable ordering**: sort by chromosome and base-pair position.
5. **Derived quantities and format mapping**:
   - Uses `mtag_z` directly (equivalent to `beta / se` under standard assumptions).
   - Builds `panel_variant_id = chr<CHR>_<BP>_<A1>_<A2>_b37` for PrediXcan harmonization.

## What gets generated

From each MTAG file under `3_MTAG/**/results/*_trait_*.txt`, the script creates:

- `output/ldsc_input/*.ldsc.tsv.gz`: LDSC-friendly table
- `output/magma/*.genes.raw`: MAGMA SNP-level p-value table
- `output/spredixcan/*.spredixcan.tsv.gz`: S-PrediXcan GWAS table
- `output/manifest.tsv`: provenance table linking every output back to source file and trait metadata

Optional:

- `output/ldsc_munged/*.sumstats.gz` if `--run-ldsc-munge` is used

## Expected input and output columns

### Input MTAG columns

- Variant identity: `SNP`, `CHR`, `BP`, `A1`, `A2`
- Study/statistics: `N`, `FRQ`, `mtag_beta`, `mtag_se`, `mtag_z`, `mtag_pval`

### Output schemas

- **LDSC input**: `SNP A1 A2 N Z P FRQ CHR BP BETA SE`
- **MAGMA input**: `SNP CHR BP P N`
- **S-PrediXcan input**: `rsid chromosome position effect_allele non_effect_allele beta se pvalue zscore sample_size effect_allele_frequency panel_variant_id`

## Running the step

From `5_Functional_Genomics/`:

```bash
python 0_munge/scripts/munge_sumstats.py \
  --run-ldsc-munge \
  --hm3-snplist reference_data/ldsc/w_hm3.snplist \
  --overwrite
```

Useful flags:

- `--min-maf`: require a minimum minor allele frequency
- `--run-ldsc-munge`: additionally call LDSC's `munge_sumstats.py`
- `--input-glob`: customize MTAG file discovery pattern

## Notes on trait naming and provenance

- Trait and cohort metadata are inferred from directory/file names and encoded into an `analysis_id`.
- `output/manifest.tsv` is the canonical record for reproducibility (source file, row counts before/after QC, and all generated paths).

## References and documentation

- LDSC repository and docs: https://github.com/bulik/ldsc
- LDSC `munge_sumstats.py` wiki: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
- MTAG paper (Turley et al., 2018): https://www.nature.com/articles/s41588-017-0009-4
- PrediXcan / MetaXcan repository: https://github.com/hakyimlab/MetaXcan
