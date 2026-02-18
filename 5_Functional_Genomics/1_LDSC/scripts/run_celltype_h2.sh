#!/usr/bin/env bash
set -euo pipefail

# Run cell-type-specific partitioned heritability (LDSC-CTS).

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

SUMSTATS_DIR="${SUMSTATS_DIR:-${FG_ROOT}/0_munge/output/ldsc_munged}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/1_LDSC/output/celltype_h2}"
PYTHON_BIN="${PYTHON_BIN:-python}"
LDSC_SCRIPT="${LDSC_SCRIPT:-${FG_ROOT}/reference_data/ldsc/ldsc/ldsc.py}"
BASELINE_LD_PREFIX="${BASELINE_LD_PREFIX:-${FG_ROOT}/reference_data/ldsc/baselineLD_v2.2/baselineLD.}"
W_LD_CHR="${W_LD_CHR:-${FG_ROOT}/reference_data/ldsc/weights_hm3_no_hla/weights.hm3_noMHC.}"
LDCTS_FILE="${LDCTS_FILE:-${FG_ROOT}/1_LDSC/config/cell_type_groups.ldcts}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${LDSC_SCRIPT}" ]]; then
    echo "ERROR: Could not find ldsc.py at ${LDSC_SCRIPT}" >&2
    exit 1
fi

if [[ ! -f "${LDCTS_FILE}" ]]; then
    echo "ERROR: LDCTS file not found: ${LDCTS_FILE}" >&2
    echo "Create this file with lines: <cell_type_name><TAB><annotation_prefix>" >&2
    exit 1
fi

shopt -s nullglob
sumstats_files=("${SUMSTATS_DIR}"/*.sumstats.gz)
shopt -u nullglob

if [[ "${#sumstats_files[@]}" -eq 0 ]]; then
    echo "ERROR: No *.sumstats.gz files found in ${SUMSTATS_DIR}" >&2
    exit 1
fi

for sumstats in "${sumstats_files[@]}"; do
    base="$(basename "${sumstats}" .sumstats.gz)"
    out_prefix="${OUT_DIR}/${base}"
    echo "[celltype-h2] ${base}"
    "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
        --h2-cts "${sumstats}" \
        --ref-ld-chr "${BASELINE_LD_PREFIX}" \
        --w-ld-chr "${W_LD_CHR}" \
        --ref-ld-chr-cts "${LDCTS_FILE}" \
        --out "${out_prefix}"
done

echo "Completed cell-type heritability runs."
