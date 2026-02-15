#!/bin/bash
# extract_phenotype.sh
#
# Extract a single phenotype column from a UKBB basket (.tab) file
# and produce a PLINK-compatible TSV (FID, IID, Phenotype).
#
# Usage:
#   bash extract_phenotype.sh FIELD_ID INSTANCE BASKET_FILE OUTPUT_DIR
#
# Arguments:
#   FIELD_ID    - UKBB data-field ID (e.g. 2020)
#   INSTANCE    - Assessment instance (0 = initial visit, 2 = first repeat, etc.)
#   BASKET_FILE - Path to the decompressed .tab basket file
#   OUTPUT_DIR  - Directory to write the output TSV into
#
# Output:
#   OUTPUT_DIR/phenotype_FIELDID.tsv
#
# Example:
#   bash extract_phenotype.sh 2020 0 /path/to/ukb675079.tab /path/to/output

set -euo pipefail

FIELD_ID="${1:?Usage: extract_phenotype.sh FIELD_ID INSTANCE BASKET_FILE OUTPUT_DIR}"
INSTANCE="${2:?Missing INSTANCE argument}"
BASKET_FILE="${3:?Missing BASKET_FILE argument}"
OUTPUT_DIR="${4:?Missing OUTPUT_DIR argument}"

COLUMN_PATTERN="f.${FIELD_ID}.${INSTANCE}.0"
OUTPUT_FILE="${OUTPUT_DIR}/phenotype_${FIELD_ID}.tsv"

echo "============================================"
echo "Extracting UKBB Field ${FIELD_ID} (instance ${INSTANCE})"
echo "  Basket file : ${BASKET_FILE}"
echo "  Column match: ${COLUMN_PATTERN}"
echo "  Output      : ${OUTPUT_FILE}"
echo "============================================"

# Verify basket file exists
if [[ ! -f "${BASKET_FILE}" ]]; then
    echo "ERROR: Basket file not found: ${BASKET_FILE}"
    exit 1
fi

# Create output directory if needed
mkdir -p "${OUTPUT_DIR}"

# Step 1: Find the column number for the target field
echo "Finding column number for ${COLUMN_PATTERN}..."
COL_NUM=$(head -n 1 "${BASKET_FILE}" | tr '\t' '\n' | nl | grep "${COLUMN_PATTERN}" | awk '{print $1}')

if [[ -z "${COL_NUM}" ]]; then
    echo "ERROR: Column '${COLUMN_PATTERN}' not found in basket file header."
    echo "Available columns matching field ${FIELD_ID}:"
    head -n 1 "${BASKET_FILE}" | tr '\t' '\n' | nl | grep "f.${FIELD_ID}" || echo "  (none found)"
    exit 1
fi

echo "  Found at column ${COL_NUM}"

# Step 2: Extract f.eid (column 1) and the phenotype column
echo "Extracting columns..."
awk -F'\t' -v col="${COL_NUM}" '{print $1"\t"$1"\t"$col}' "${BASKET_FILE}" > "${OUTPUT_FILE}.tmp"

# Step 3: Replace header with PLINK-compatible names
sed -i '1s/.*/FID\tIID\tPhenotype/' "${OUTPUT_FILE}.tmp"

# Step 4: Remove rows where phenotype is NA or empty
echo "Removing missing values..."
awk -F'\t' 'NR==1 || ($3 != "NA" && $3 != "" && $3 != "nan")' "${OUTPUT_FILE}.tmp" > "${OUTPUT_FILE}"
rm -f "${OUTPUT_FILE}.tmp"

# Step 5: Report summary
TOTAL=$(wc -l < "${OUTPUT_FILE}")
echo ""
echo "Done. Output: ${OUTPUT_FILE}"
echo "  Total rows (including header): ${TOTAL}"
echo "  Preview:"
head -n 5 "${OUTPUT_FILE}"
echo "============================================"
