#!/usr/bin/env bash
set -euo pipefail

# Run S-PrediXcan for each local GWAS file against all available tissue models.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python}"
GWAS_DIR="${GWAS_DIR:-${FG_ROOT}/0_munge/output/spredixcan}"
MODELS_DIR="${MODELS_DIR:-${FG_ROOT}/reference_data/predixcan/models}"
COVARIANCE_DIR="${COVARIANCE_DIR:-${FG_ROOT}/reference_data/predixcan/covariances}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/3_S-PrediXcan/output/spredixcan}"
SPREDIXCAN_SCRIPT="${SPREDIXCAN_SCRIPT:-${FG_ROOT}/reference_data/MetaXcan/software/SPrediXcan.py}"
MODEL_SNP_KEY="${MODEL_SNP_KEY:-varID}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${SPREDIXCAN_SCRIPT}" ]]; then
    echo "ERROR: SPrediXcan.py not found at ${SPREDIXCAN_SCRIPT}" >&2
    echo "Clone MetaXcan and point SPREDIXCAN_SCRIPT to software/SPrediXcan.py." >&2
    exit 1
fi

shopt -s nullglob
gwas_files=("${GWAS_DIR}"/*.spredixcan.tsv.gz)
model_files=("${MODELS_DIR}"/*.db)
shopt -u nullglob

if [[ "${#gwas_files[@]}" -eq 0 ]]; then
    echo "ERROR: No SPrediXcan-formatted GWAS files in ${GWAS_DIR}" >&2
    exit 1
fi
if [[ "${#model_files[@]}" -eq 0 ]]; then
    echo "ERROR: No tissue model DB files (*.db) found in ${MODELS_DIR}" >&2
    exit 1
fi

find_covariance_file() {
    local model_basename="$1"
    local candidates=(
        "${COVARIANCE_DIR}/${model_basename}.txt.gz"
        "${COVARIANCE_DIR}/${model_basename}.covariances.txt.gz"
        "${COVARIANCE_DIR}/${model_basename}.txt"
    )
    local c
    for c in "${candidates[@]}"; do
        if [[ -f "${c}" ]]; then
            echo "${c}"
            return 0
        fi
    done
    return 1
}

run_count=0
for gwas_file in "${gwas_files[@]}"; do
    gwas_base="$(basename "${gwas_file}" .spredixcan.tsv.gz)"
    for model_db in "${model_files[@]}"; do
        model_base="$(basename "${model_db}" .db)"
        if ! cov_file="$(find_covariance_file "${model_base}")"; then
            echo "[skip] Covariance file not found for ${model_base}"
            continue
        fi

        out_file="${OUT_DIR}/${gwas_base}__${model_base}.csv"
        echo "[spredixcan] ${gwas_base} x ${model_base}"
        "${PYTHON_BIN}" "${SPREDIXCAN_SCRIPT}" \
            --gwas_file "${gwas_file}" \
            --model_db_path "${model_db}" \
            --covariance "${cov_file}" \
            --snp_column rsid \
            --chromosome_column chromosome \
            --position_column position \
            --effect_allele_column effect_allele \
            --non_effect_allele_column non_effect_allele \
            --beta_column beta \
            --se_column se \
            --pvalue_column pvalue \
            --keep_non_rsid \
            --model_db_snp_key "${MODEL_SNP_KEY}" \
            --throw \
            --output_file "${out_file}"
        run_count=$((run_count + 1))
    done
done

echo "Completed ${run_count} S-PrediXcan runs."
