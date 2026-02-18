#!/usr/bin/env bash
set -euo pipefail

# MAGMA step 1: SNP-to-gene annotation.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

MAGMA_BIN="${MAGMA_BIN:-magma}"
MAGMA_INPUT_DIR="${MAGMA_INPUT_DIR:-${FG_ROOT}/0_munge/output/magma}"
GENE_LOC_FILE="${GENE_LOC_FILE:-${FG_ROOT}/reference_data/magma/NCBI37.3.gene.loc}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/2_MAGMA/output/01_annotation}"
WINDOW_KB="${WINDOW_KB:-10}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${GENE_LOC_FILE}" ]]; then
    echo "ERROR: GENE_LOC_FILE not found: ${GENE_LOC_FILE}" >&2
    exit 1
fi

if ! command -v "${MAGMA_BIN}" >/dev/null 2>&1; then
    echo "ERROR: MAGMA binary not found on PATH: ${MAGMA_BIN}" >&2
    exit 1
fi

shopt -s nullglob
inputs=("${MAGMA_INPUT_DIR}"/*.genes.raw)
shopt -u nullglob

if [[ "${#inputs[@]}" -eq 0 ]]; then
    echo "ERROR: No *.genes.raw files found in ${MAGMA_INPUT_DIR}" >&2
    exit 1
fi

for snp_file in "${inputs[@]}"; do
    base="$(basename "${snp_file}" .genes.raw)"
    out_prefix="${OUT_DIR}/${base}"
    echo "[magma-annotate] ${base}"
    "${MAGMA_BIN}" \
        --annotate window="${WINDOW_KB},${WINDOW_KB}" \
        --snp-loc "${snp_file}" \
        --gene-loc "${GENE_LOC_FILE}" \
        --out "${out_prefix}"
done

echo "Completed MAGMA annotation for ${#inputs[@]} files."
