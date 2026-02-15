#!/bin/bash
set -beEo pipefail

# Filter phenotype and covariate files to EUR_Male ancestry samples
# This is simpler than using --remove in BOLT-LMM

ukb21942_d='/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942'
SRCDIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness"

echo "========================================"
echo "Filter Phenotype & Covariate Files to EUR_Male"
echo "========================================"
echo ""

# Input files
keep_file="${ukb21942_d}/sqc/population.20220316/EUR_Male.keep"
pheno_file="${ukb21942_d}/pheno/isolation_run_control.tsv.gz"
covar_file="${ukb21942_d}/sqc/sqc.20220316.tsv.gz"

# Output files (in analysis directory for easy access)
pheno_eur_male="${SRCDIR}/isolation_run_control.EUR_Male.tsv.gz"
covar_eur_male="${SRCDIR}/sqc.EUR_Male.tsv.gz"

echo "Input files:"
echo "  EUR_Male samples: ${keep_file}"
echo "  Phenotypes: ${pheno_file}"
echo "  Covariates: ${covar_file}"
echo ""
echo "Output files:"
echo "  EUR_Male phenotypes: ${pheno_eur_male}"
echo "  EUR_Male covariates: ${covar_eur_male}"
echo ""

# Check if EUR_Male.keep exists
if [ ! -f "${keep_file}" ]; then
    echo "ERROR: EUR_Male.keep file not found: ${keep_file}" >&2
    exit 1
fi

# Count EUR_Male samples
n_eur_male=$(wc -l < "${keep_file}")
echo "EUR_Male samples to keep: ${n_eur_male}"
echo ""

# Filter phenotype file
echo "Filtering phenotype file..."
echo "This may take a minute..."

# Create temporary ID lookup file (just IIDs from keep file)
awk '{print $2}' "${keep_file}" > /tmp/eur_male_iids.txt
n_ids=$(wc -l < /tmp/eur_male_iids.txt)
echo "  EUR_Male IDs to match: ${n_ids}"

# Method: Use grep with file of patterns (much faster and more reliable)
{
    # Extract and write header
    zcat "${pheno_file}" | head -1
    
    # Extract data rows and filter using grep
    # -F: fixed strings (not regex)
    # -f: patterns from file
    zcat "${pheno_file}" | tail -n +2 | grep -F -f /tmp/eur_male_iids.txt
    
} | gzip > "${pheno_eur_male}"

n_pheno_out=$(zcat "${pheno_eur_male}" | tail -n +2 | wc -l)
echo "✓ EUR_Male phenotype file created"
echo "  Input samples: $(zcat "${pheno_file}" | tail -n +2 | wc -l)"
echo "  Output samples: ${n_pheno_out}"
echo ""

# Filter covariate file
echo "Filtering covariate file..."
echo "This may take a minute..."

{
    # Header
    zcat "${covar_file}" | head -1
    
    # Data rows filtered to EUR_Male
    zcat "${covar_file}" | tail -n +2 | grep -F -f /tmp/eur_male_iids.txt
    
} | gzip > "${covar_eur_male}"

n_covar_out=$(zcat "${covar_eur_male}" | tail -n +2 | wc -l)
echo "✓ EUR_Male covariate file created"
echo "  Input samples: $(zcat "${covar_file}" | tail -n +2 | wc -l)"
echo "  Output samples: ${n_covar_out}"
echo ""

# Clean up
rm -f /tmp/eur_male_iids.txt

echo ""
echo "========================================"
echo "Filtering Complete!"
echo "========================================"
echo ""
echo "Summary:"
echo "  EUR_Male samples requested: ${n_eur_male}"
echo "  Phenotype file output: ${n_pheno_out}"
echo "  Covariate file output: ${n_covar_out}"
echo ""

if [ ${n_pheno_out} -ne ${n_eur_male} ]; then
    echo "⚠️  WARNING: Phenotype sample count doesn't match EUR_Male count"
    echo "   This is normal if some EUR_Male samples have missing phenotype data"
fi

if [ ${n_covar_out} -ne ${n_eur_male} ]; then
    echo "⚠️  WARNING: Covariate sample count doesn't match EUR_Male count"
    echo "   This is normal if some EUR_Male samples have missing covariate data"
fi

echo ""
echo "Next steps:"
echo "1. Run test: sbatch 0b_test_run.sbatch.sh"
echo "2. If test passes: sbatch 1_run_bolt_lmm.sbatch.sh"

