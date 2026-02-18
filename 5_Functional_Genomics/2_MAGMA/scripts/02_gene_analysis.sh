#!/usr/bin/env bash
set -euo pipefail

# MAGMA step 2: Gene-level association analysis.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

MAGMA_BIN="${MAGMA_BIN:-magma}"
MAGMA_INPUT_DIR="${MAGMA_INPUT_DIR:-${FG_ROOT}/0_munge/output/magma}"
ANNOT_DIR="${ANNOT_DIR:-${FG_ROOT}/2_MAGMA/output/01_annotation}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/2_MAGMA/output/02_gene_analysis}"
BFILE_PREFIX="${BFILE_PREFIX:-${FG_ROOT}/reference_data/1000G_EUR/g1000_eur}"

mkdir -p "${OUT_DIR}"

if ! command -v "${MAGMA_BIN}" >/dev/null 2>&1; then
    echo "ERROR: MAGMA binary not found on PATH: ${MAGMA_BIN}" >&2
    exit 1
fi

if [[ ! -f "${BFILE_PREFIX}.bed" || ! -f "${BFILE_PREFIX}.bim" || ! -f "${BFILE_PREFIX}.fam" ]]; then
    echo "ERROR: 1000G EUR reference PLINK files not found for prefix ${BFILE_PREFIX}" >&2
    exit 1
fi

shopt -s nullglob
annot_files=("${ANNOT_DIR}"/*.genes.annot)
shopt -u nullglob

if [[ "${#annot_files[@]}" -eq 0 ]]; then
    echo "ERROR: No *.genes.annot files found in ${ANNOT_DIR}" >&2
    echo "Run 2_MAGMA/scripts/01_annotate.sh first." >&2
    exit 1
fi

for annot_file in "${annot_files[@]}"; do
    base="$(basename "${annot_file}" .genes.annot)"
    pval_file="${MAGMA_INPUT_DIR}/${base}.genes.raw"
    if [[ ! -f "${pval_file}" ]]; then
        echo "[skip] Missing p-value file for ${base}: ${pval_file}" >&2
        continue
    fi

    out_prefix="${OUT_DIR}/${base}"
    echo "[magma-gene] ${base}"
    "${MAGMA_BIN}" \
        --bfile "${BFILE_PREFIX}" \
        --pval "${pval_file}" use=SNP,P ncol=N \
        --gene-annot "${annot_file}" \
        --out "${out_prefix}"
done

echo "Completed MAGMA gene-level analysis."
