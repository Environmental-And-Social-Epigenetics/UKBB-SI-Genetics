#!/usr/bin/env bash
set -euo pipefail

# Run partitioned heritability using baselineLD annotations.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

SUMSTATS_DIR="${SUMSTATS_DIR:-${FG_ROOT}/0_munge/output/ldsc_munged}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/1_LDSC/output/partitioned_h2}"
PYTHON_BIN="${PYTHON_BIN:-python}"
LDSC_SCRIPT="${LDSC_SCRIPT:-${FG_ROOT}/reference_data/ldsc/ldsc/ldsc.py}"
BASELINE_LD_PREFIX="${BASELINE_LD_PREFIX:-${FG_ROOT}/reference_data/ldsc/baselineLD_v2.2/baselineLD.}"
W_LD_CHR="${W_LD_CHR:-${FG_ROOT}/reference_data/ldsc/weights_hm3_no_hla/weights.hm3_noMHC.}"
FRQFILE_CHR="${FRQFILE_CHR:-${FG_ROOT}/reference_data/ldsc/1000G_Phase3_frq/1000G.EUR.QC.}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${LDSC_SCRIPT}" ]]; then
    echo "ERROR: Could not find ldsc.py at ${LDSC_SCRIPT}" >&2
    exit 1
fi

shopt -s nullglob
sumstats_files=("${SUMSTATS_DIR}"/*.sumstats.gz)
shopt -u nullglob

if [[ "${#sumstats_files[@]}" -eq 0 ]]; then
    echo "ERROR: No *.sumstats.gz files found in ${SUMSTATS_DIR}" >&2
    exit 1
fi

extra_args=(--overlap-annot --print-coefficients)
if [[ -f "${FRQFILE_CHR}1.frq.gz" ]]; then
    extra_args+=(--frqfile-chr "${FRQFILE_CHR}")
else
    echo "[warn] FRQ files not found at ${FRQFILE_CHR}*. Proceeding without --frqfile-chr."
fi

for sumstats in "${sumstats_files[@]}"; do
    base="$(basename "${sumstats}" .sumstats.gz)"
    out_prefix="${OUT_DIR}/${base}"
    echo "[partitioned-h2] ${base}"
    "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
        --h2 "${sumstats}" \
        --ref-ld-chr "${BASELINE_LD_PREFIX}" \
        --w-ld-chr "${W_LD_CHR}" \
        "${extra_args[@]}" \
        --out "${out_prefix}"
done

echo "Completed partitioned heritability runs."
