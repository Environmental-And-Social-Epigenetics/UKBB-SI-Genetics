# Phenotype Creation

This directory contains the full workflow for extracting and formatting UK Biobank phenotypes used in the social isolation GWAS pipeline.

## Background

The phenotypes and binary coding scheme used in this analysis are based on:

> Day, F.R., Ong, K.K. & Perry, J.R.B. **Elucidating the genetic basis of social interaction and isolation.** *Nature Communications* 9, 2457 (2018). https://doi.org/10.1038/s41467-018-04930-1

Day et al. defined three social isolation traits -- Loneliness, Ability to Confide, and Frequency of Social Contact (FreqSoc) -- from UK Biobank touchscreen questionnaire items and applied binary case/control coding for GWAS. We adopt their phenotype definitions and binary coding directly.

In addition to the binary analysis, we introduce a **continuous coding** of the same traits to capture finer-grained variation in social isolation. The continuous coding approach was conceived independently for this study and is not part of the Day et al. methodology.

## Phenotypes

Four raw phenotypes are extracted from the UKBB basket file. All use **Instance 0** (initial assessment visit).

| Trait Name | UKBB Field ID | Description | Showcase Link |
|---|---|---|---|
| Loneliness | 2020 | "Do you often feel lonely?" (No/Yes) | [Field 2020](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2020) |
| Ability to Confide | 2110 | "How often are you able to confide in someone close to you?" | [Field 2110](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2110) |
| Frequency of Friend/Family Visits | 1031 | "How often do you visit friends or family or have them visit you?" | [Field 1031](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1031) |
| Number in Household | 709 | "Including yourself, how many people are living together in your household?" | [Field 709](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=709) |

These four raw phenotypes produce **three derived GWAS traits**:

| Derived Trait | Raw Fields | Description |
|---|---|---|
| **Loneliness** | Field 2020 | Direct mapping from the loneliness question |
| **AbilityToConfide** | Field 2110 | Recoded measure of social support |
| **FreqSoc** | Fields 1031 + 709 | Composite of visit frequency and household size |

## Workflow

### Step 0: Decompress the Basket File (if needed)

The UKBB basket file on Luria is stored compressed. Decompress it before extraction:

```bash
# Basket file location: /home/mabdel03/data/files/Isolation_Genetics/GWAS/Basket_File/
zstd -d ukb675079.tab.zst -o ukb675079.tab
```

### Step 1: Extract Raw Phenotypes from Basket File

Use the extraction scripts to pull phenotype columns from the decompressed basket file:

```bash
# Extract all four phenotypes at once (uses default paths on Luria)
bash extract_all_phenotypes.sh

# Or with custom paths:
bash extract_all_phenotypes.sh /path/to/ukb675079.tab /path/to/output_dir
```

This calls `extract_phenotype.sh` for each field, which:
1. Finds the column matching `f.FIELD_ID.0.0` in the basket header
2. Extracts the participant ID and phenotype value using `awk`
3. Renames the header to PLINK-compatible format (`FID`, `IID`, `Phenotype`)
4. Removes rows with missing values

Output files: `phenotype_2020.tsv`, `phenotype_2110.tsv`, `phenotype_1031.tsv`, `phenotype_709.tsv`

### Step 2: Format Binary Phenotype File (Day et al. coding)

Run `Binary_Pheno_Formatting.ipynb` to produce PLINK case/control coding (1 = control, 2 = case), following the definitions in Day et al. (2018). All traits are oriented so that **higher values indicate greater social isolation**:

- **Loneliness**: 0 (No) -> 1 (control), 1 (Yes) -> 2 (case)
- **AbilityToConfide**: 0 (never/almost never) -> 2 (case, isolated); any other response -> 1 (control)
- **FreqSoc**: Lives alone (NumHousehold = 1) AND rarely visited (FreqVisit in {6, 7}) -> 2 (case); otherwise -> 1 (control)

Output: `isolation_run_binary.tsv.gz`

### Step 3: Format Continuous Phenotype File (novel coding)

Run `Continuous_Pheno_Formatting.ipynb` to produce continuous-valued coding. This coding was developed for this study to complement the binary analysis. All traits are oriented so that **higher values indicate greater social isolation**:

- **Loneliness**: Raw 0/1 value
- **AbilityToConfide**: 5 - raw value (reverses the UKBB scale)
- **FreqSoc**: Average of z-scored FreqVisit and z-scored negated NumHousehold

Output: `isolation_run_continuous.tsv.gz`

## Output Files

Both output files share the same column structure:

| Column | Description |
|---|---|
| `FID` | UKBB participant ID (family ID) |
| `IID` | UKBB participant ID (individual ID, same as FID) |
| `Loneliness` | Loneliness trait value |
| `AbilityToConfide` | Ability to confide trait value |
| `FreqSoc` | Frequency of socialising composite trait value |

Output location on Luria: `/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/pheno/`

These files are consumed by the BOLT-LMM GWAS scripts in `../2_GWAS/`.

## Files in This Directory

| File | Description |
|---|---|
| `extract_phenotype.sh` | Reusable script to extract a single field from the UKBB basket file |
| `extract_all_phenotypes.sh` | Master script that extracts all four SI phenotypes |
| `Binary_Pheno_Formatting.ipynb` | Notebook to create binary-coded (case/control) phenotype file |
| `Continuous_Pheno_Formatting.ipynb` | Notebook to create continuous-coded phenotype file |
| `Pheno_Formatting.ipynb` | Original combined notebook (retained as reference) |

## Data Locations on Luria

| Resource | Path |
|---|---|
| Basket file (compressed) | `/home/mabdel03/data/files/Isolation_Genetics/GWAS/Basket_File/ukb675079.tab.zst` |
| Basket file (decompressed) | `/home/mabdel03/data/files/Isolation_Genetics/GWAS/Basket_File/ukb675079.tab` |
| Extracted phenotype TSVs | `/home/mabdel03/data/files/Isolation_Genetics/GWAS/Basket_File/extracted_phenotypes/` |
| Formatted output files | `/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/pheno/` |

## References

- Day, F.R., Ong, K.K. & Perry, J.R.B. Elucidating the genetic basis of social interaction and isolation. *Nat Commun* 9, 2457 (2018). https://doi.org/10.1038/s41467-018-04930-1
