#!/usr/bin/env bash
set -euo pipefail

# MAGMA step 4: Tissue-expression enrichment using gene-property analysis.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

MAGMA_BIN="${MAGMA_BIN:-magma}"
GENE_RESULTS_DIR="${GENE_RESULTS_DIR:-${FG_ROOT}/2_MAGMA/output/02_gene_analysis}"
GENE_COVAR_FILE="${GENE_COVAR_FILE:-${FG_ROOT}/reference_data/magma/gtex_v8/GTEx_v8_tissue_expression.gene_covar.txt}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/2_MAGMA/output/04_tissue_expression}"

mkdir -p "${OUT_DIR}"

if ! command -v "${MAGMA_BIN}" >/dev/null 2>&1; then
    echo "ERROR: MAGMA binary not found on PATH: ${MAGMA_BIN}" >&2
    exit 1
fi

if [[ ! -f "${GENE_COVAR_FILE}" ]]; then
    echo "ERROR: Tissue covariate file not found: ${GENE_COVAR_FILE}" >&2
    echo "Set GENE_COVAR_FILE or download GTEx covariates into reference_data/magma/gtex_v8." >&2
    exit 1
fi

shopt -s nullglob
gene_results=("${GENE_RESULTS_DIR}"/*.genes.raw)
shopt -u nullglob

if [[ "${#gene_results[@]}" -eq 0 ]]; then
    echo "ERROR: No gene results (*.genes.raw) found in ${GENE_RESULTS_DIR}" >&2
    exit 1
fi

for gr in "${gene_results[@]}"; do
    trait_base="$(basename "${gr}" .genes.raw)"
    out_prefix="${OUT_DIR}/${trait_base}"
    echo "[magma-tissue] ${trait_base}"
    "${MAGMA_BIN}" \
        --gene-results "${gr}" \
        --gene-covar "${GENE_COVAR_FILE}" \
        --model direction=pos condition-hide=Average \
        --out "${out_prefix}"
done

echo "Completed MAGMA tissue-expression analysis."
