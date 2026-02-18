#!/usr/bin/env bash
set -euo pipefail

# MAGMA step 3: Gene-set enrichment analysis.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

MAGMA_BIN="${MAGMA_BIN:-magma}"
GENE_RESULTS_DIR="${GENE_RESULTS_DIR:-${FG_ROOT}/2_MAGMA/output/02_gene_analysis}"
GENESET_DIR="${GENESET_DIR:-${FG_ROOT}/reference_data/magma/msigdb}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/2_MAGMA/output/03_geneset_analysis}"

mkdir -p "${OUT_DIR}"

if ! command -v "${MAGMA_BIN}" >/dev/null 2>&1; then
    echo "ERROR: MAGMA binary not found on PATH: ${MAGMA_BIN}" >&2
    exit 1
fi

shopt -s nullglob
gene_results=("${GENE_RESULTS_DIR}"/*.genes.raw)
set_files=("${GENESET_DIR}"/*.gmt)
shopt -u nullglob

if [[ "${#gene_results[@]}" -eq 0 ]]; then
    echo "ERROR: No gene results (*.genes.raw) found in ${GENE_RESULTS_DIR}" >&2
    echo "Run 2_MAGMA/scripts/02_gene_analysis.sh first." >&2
    exit 1
fi

if [[ "${#set_files[@]}" -eq 0 ]]; then
    echo "ERROR: No gene set files (*.gmt) found in ${GENESET_DIR}" >&2
    exit 1
fi

for gr in "${gene_results[@]}"; do
    trait_base="$(basename "${gr}" .genes.raw)"
    for set_file in "${set_files[@]}"; do
        set_base="$(basename "${set_file}" .gmt)"
        out_prefix="${OUT_DIR}/${trait_base}__${set_base}"
        echo "[magma-geneset] ${trait_base} vs ${set_base}"
        "${MAGMA_BIN}" \
            --gene-results "${gr}" \
            --set-annot "${set_file}" \
            --out "${out_prefix}"
    done
done

echo "Completed MAGMA gene-set analysis."
