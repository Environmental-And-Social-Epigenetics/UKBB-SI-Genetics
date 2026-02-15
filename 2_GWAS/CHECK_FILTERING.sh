#!/bin/bash
# Diagnostic script to check filtering status and re-run if needed

echo "========================================"
echo "BOLT-LMM_SI-Loneliness Filtering Diagnostics"
echo "========================================"
echo ""

# Check what was created
echo "1. Files that exist:"
ls -lh *.EUR*.tsv.gz 2>/dev/null || echo "  No filtered files found"
echo ""

# Check the filter output and errors
echo "2. Last filtering job output (0a_filter.out):"
if [ -f "0a_filter.out" ]; then
    tail -20 0a_filter.out
else
    echo "  File not found"
fi
echo ""

echo "3. Last filtering job errors (0a_filter.err):"
if [ -f "0a_filter.err" ]; then
    cat 0a_filter.err
else
    echo "  File not found"
fi
echo ""

# Check if source files exist
echo "4. Checking source files:"
ukb21942_d='/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942'

files_to_check=(
    "${ukb21942_d}/sqc/population.20220316/EUR_MM.keep"
    "${ukb21942_d}/sqc/population.20220316/EUR_Male.keep"
    "${ukb21942_d}/sqc/population.20220316/EUR_Female.keep"
    "${ukb21942_d}/pheno/isolation_run_control.tsv.gz"
    "${ukb21942_d}/sqc/sqc.20220316.tsv.gz"
)

for file in "${files_to_check[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ MISSING: $file"
    fi
done
echo ""

# Expected output files
echo "5. Expected filtered files (6 total):"
expected_files=(
    "isolation_run_control.EUR_MM.tsv.gz"
    "isolation_run_control.EUR_Male.tsv.gz"
    "isolation_run_control.EUR_Female.tsv.gz"
    "sqc.EUR_MM.tsv.gz"
    "sqc.EUR_Male.tsv.gz"
    "sqc.EUR_Female.tsv.gz"
)

for file in "${expected_files[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file ($(du -h $file | cut -f1))"
    else
        echo "  ✗ MISSING: $file"
    fi
done
echo ""

echo "6. Next steps:"
echo "  - If source files are missing, they need to be created first"
echo "  - If filtering failed, check the error messages above"
echo "  - To re-run filtering: sbatch 0a_filter_populations.sbatch.sh"
echo "  - To manually filter one population: bash filter_to_population.sh EUR_MM"
echo ""

