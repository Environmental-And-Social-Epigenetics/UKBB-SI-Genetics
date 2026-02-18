#!/usr/bin/env bash
set -euo pipefail

# Estimate SNP heritability for each munged LDSC summary statistics file.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

SUMSTATS_DIR="${SUMSTATS_DIR:-${FG_ROOT}/0_munge/output/ldsc_munged}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/1_LDSC/output/h2}"
PYTHON_BIN="${PYTHON_BIN:-python}"
LDSC_SCRIPT="${LDSC_SCRIPT:-${FG_ROOT}/reference_data/ldsc/ldsc/ldsc.py}"
REF_LD_CHR="${REF_LD_CHR:-${FG_ROOT}/reference_data/ldsc/eur_w_ld_chr/}"
W_LD_CHR="${W_LD_CHR:-${FG_ROOT}/reference_data/ldsc/weights_hm3_no_hla/weights.hm3_noMHC.}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${LDSC_SCRIPT}" ]]; then
    echo "ERROR: Could not find ldsc.py at ${LDSC_SCRIPT}" >&2
    echo "Set LDSC_SCRIPT or run reference_data/download_references.sh first." >&2
    exit 1
fi

if [[ ! -d "${SUMSTATS_DIR}" ]]; then
    echo "ERROR: SUMSTATS_DIR not found: ${SUMSTATS_DIR}" >&2
    exit 1
fi

shopt -s nullglob
sumstats_files=("${SUMSTATS_DIR}"/*.sumstats.gz)
shopt -u nullglob

if [[ "${#sumstats_files[@]}" -eq 0 ]]; then
    echo "ERROR: No *.sumstats.gz files found in ${SUMSTATS_DIR}" >&2
    echo "Run 0_munge/scripts/munge_sumstats.py with --run-ldsc-munge first." >&2
    exit 1
fi

for sumstats in "${sumstats_files[@]}"; do
    base="$(basename "${sumstats}" .sumstats.gz)"
    out_prefix="${OUT_DIR}/${base}"
    echo "[h2] ${base}"

    "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
        --h2 "${sumstats}" \
        --ref-ld-chr "${REF_LD_CHR}" \
        --w-ld-chr "${W_LD_CHR}" \
        --out "${out_prefix}"
done

echo "Completed h2 estimation for ${#sumstats_files[@]} traits."
