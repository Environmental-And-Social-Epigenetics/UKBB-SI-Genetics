#!/bin/bash
# extract_all_phenotypes.sh
#
# Extract all four social-isolation phenotype columns from the UKBB basket file.
# This is the master script that calls extract_phenotype.sh for each field.
#
# Prerequisites:
#   - The basket file must be decompressed first if stored as .tab.zst:
#       zstd -d ukb675079.tab.zst -o ukb675079.tab
#
# Usage:
#   bash extract_all_phenotypes.sh [BASKET_FILE] [OUTPUT_DIR]
#
# Defaults (on Luria):
#   BASKET_FILE = /home/mabdel03/data/files/Isolation_Genetics/GWAS/Basket_File/ukb675079.tab
#   OUTPUT_DIR  = /home/mabdel03/data/files/Isolation_Genetics/GWAS/Basket_File/extracted_phenotypes

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BASKET_FILE="${1:-/home/mabdel03/data/files/Isolation_Genetics/GWAS/Basket_File/ukb675079.tab}"
OUTPUT_DIR="${2:-/home/mabdel03/data/files/Isolation_Genetics/GWAS/Basket_File/extracted_phenotypes}"

INSTANCE=0

echo "============================================================"
echo "Social Isolation Phenotype Extraction"
echo "============================================================"
echo "Basket file : ${BASKET_FILE}"
echo "Output dir  : ${OUTPUT_DIR}"
echo "Instance    : ${INSTANCE} (initial assessment)"
echo ""
echo "Phenotypes to extract:"
echo "  - Loneliness                      (Field 2020)"
echo "  - Ability to confide              (Field 2110)"
echo "  - Frequency of friend/family visits (Field 1031)"
echo "  - Number in household             (Field 709)"
echo "============================================================"
echo ""

# Verify basket file exists
if [[ ! -f "${BASKET_FILE}" ]]; then
    echo "ERROR: Basket file not found: ${BASKET_FILE}"
    echo ""
    echo "If the file is compressed, decompress it first:"
    echo "  zstd -d ukb675079.tab.zst -o ukb675079.tab"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

FAILED=0

# --- Loneliness (Field 2020) ---
echo "[1/4] Extracting Loneliness (Field 2020)..."
if bash "${SCRIPT_DIR}/extract_phenotype.sh" 2020 "${INSTANCE}" "${BASKET_FILE}" "${OUTPUT_DIR}"; then
    echo "[1/4] OK"
else
    echo "[1/4] FAILED"
    FAILED=$((FAILED + 1))
fi
echo ""

# --- Ability to Confide (Field 2110) ---
echo "[2/4] Extracting Ability to Confide (Field 2110)..."
if bash "${SCRIPT_DIR}/extract_phenotype.sh" 2110 "${INSTANCE}" "${BASKET_FILE}" "${OUTPUT_DIR}"; then
    echo "[2/4] OK"
else
    echo "[2/4] FAILED"
    FAILED=$((FAILED + 1))
fi
echo ""

# --- Frequency of Friend/Family Visits (Field 1031) ---
echo "[3/4] Extracting Frequency of Friend/Family Visits (Field 1031)..."
if bash "${SCRIPT_DIR}/extract_phenotype.sh" 1031 "${INSTANCE}" "${BASKET_FILE}" "${OUTPUT_DIR}"; then
    echo "[3/4] OK"
else
    echo "[3/4] FAILED"
    FAILED=$((FAILED + 1))
fi
echo ""

# --- Number in Household (Field 709) ---
echo "[4/4] Extracting Number in Household (Field 709)..."
if bash "${SCRIPT_DIR}/extract_phenotype.sh" 709 "${INSTANCE}" "${BASKET_FILE}" "${OUTPUT_DIR}"; then
    echo "[4/4] OK"
else
    echo "[4/4] FAILED"
    FAILED=$((FAILED + 1))
fi
echo ""

# --- Summary ---
echo "============================================================"
echo "Extraction complete."
echo ""
if [[ ${FAILED} -eq 0 ]]; then
    echo "All 4 phenotypes extracted successfully."
    echo ""
    echo "Output files:"
    ls -lh "${OUTPUT_DIR}"/phenotype_*.tsv 2>/dev/null || echo "  (no files found)"
else
    echo "WARNING: ${FAILED} extraction(s) failed. Check output above."
fi
echo ""
echo "Next steps:"
echo "  1. Review extracted files with: head <file>"
echo "  2. Run Binary_Pheno_Formatting.ipynb for binary phenotype coding"
echo "  3. Run Continuous_Pheno_Formatting.ipynb for continuous phenotype coding"
echo "============================================================"

exit ${FAILED}
