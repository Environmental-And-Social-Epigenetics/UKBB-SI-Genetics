#!/usr/bin/env bash
set -euo pipefail

# Compute genetic correlations between local SI/MRI traits and external GWAS.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

SUMSTATS_DIR="${SUMSTATS_DIR:-${FG_ROOT}/0_munge/output/ldsc_munged}"
EXTERNAL_CONFIG="${EXTERNAL_CONFIG:-${FG_ROOT}/1_LDSC/config/external_traits.tsv}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/1_LDSC/output/rg_external}"
PYTHON_BIN="${PYTHON_BIN:-python}"
LDSC_SCRIPT="${LDSC_SCRIPT:-${FG_ROOT}/reference_data/ldsc/ldsc/ldsc.py}"
REF_LD_CHR="${REF_LD_CHR:-${FG_ROOT}/reference_data/ldsc/eur_w_ld_chr/}"
W_LD_CHR="${W_LD_CHR:-${FG_ROOT}/reference_data/ldsc/weights_hm3_no_hla/weights.hm3_noMHC.}"

# Optional regex to focus local traits, e.g. "Loneliness|AbilityToConfide"
LOCAL_TRAIT_REGEX="${LOCAL_TRAIT_REGEX:-.*}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${LDSC_SCRIPT}" ]]; then
    echo "ERROR: Could not find ldsc.py at ${LDSC_SCRIPT}" >&2
    exit 1
fi

if [[ ! -d "${SUMSTATS_DIR}" ]]; then
    echo "ERROR: SUMSTATS_DIR not found: ${SUMSTATS_DIR}" >&2
    exit 1
fi

if [[ ! -f "${EXTERNAL_CONFIG}" ]]; then
    echo "ERROR: EXTERNAL_CONFIG not found: ${EXTERNAL_CONFIG}" >&2
    exit 1
fi

shopt -s nullglob
local_files=("${SUMSTATS_DIR}"/*.sumstats.gz)
shopt -u nullglob

if [[ "${#local_files[@]}" -eq 0 ]]; then
    echo "ERROR: No local *.sumstats.gz files in ${SUMSTATS_DIR}" >&2
    exit 1
fi

run_count=0
while IFS=$'\t' read -r trait_name sumstats_path sample_size source notes; do
    [[ "${trait_name}" == "trait_name" ]] && continue
    [[ -z "${trait_name}" ]] && continue

    if [[ -z "${sumstats_path}" ]]; then
        echo "[skip] ${trait_name}: empty sumstats_path in ${EXTERNAL_CONFIG}"
        continue
    fi
    if [[ ! -f "${sumstats_path}" ]]; then
        echo "[skip] ${trait_name}: file not found -> ${sumstats_path}"
        continue
    fi

    ext_token="$(echo "${trait_name}" | tr '[:space:]/' '__')"
    for local_sumstats in "${local_files[@]}"; do
        local_base="$(basename "${local_sumstats}" .sumstats.gz)"
        if [[ ! "${local_base}" =~ ${LOCAL_TRAIT_REGEX} ]]; then
            continue
        fi

        out_prefix="${OUT_DIR}/${local_base}__vs__${ext_token}"
        echo "[rg-external] ${local_base} vs ${trait_name}"
        "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
            --rg "${local_sumstats},${sumstats_path}" \
            --ref-ld-chr "${REF_LD_CHR}" \
            --w-ld-chr "${W_LD_CHR}" \
            --out "${out_prefix}"
        run_count=$((run_count + 1))
    done
done < "${EXTERNAL_CONFIG}"

echo "Completed ${run_count} external rg runs."
