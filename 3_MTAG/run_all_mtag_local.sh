#!/usr/bin/env bash
set -euo pipefail

# Run SI MTAG locally for all 6 analyses:
#   - Binary: EUR_MM, EUR_Male, EUR_Female
#   - Continuous: EUR_MM, EUR_Male, EUR_Female
#
# Optional env overrides:
#   CONDA_SH       (default: "${HOME}/opt/anaconda3/etc/profile.d/conda.sh")
#   MTAG_ENV_PATH  (default: "${HOME}/Desktop/MIT/Software/Research_Software/conda_envs/MTAG")
#   MTAG_DIR       (default: "${HOME}/Desktop/MIT/Software/Research_Software/mtag")

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPODIR="$(cd "${SCRIPTDIR}/.." && pwd)"

CONDA_SH="${CONDA_SH:-${HOME}/opt/anaconda3/etc/profile.d/conda.sh}"
MTAG_ENV_PATH="${MTAG_ENV_PATH:-${HOME}/Desktop/MIT/Software/Research_Software/conda_envs/MTAG}"
MTAG_DIR="${MTAG_DIR:-${HOME}/Desktop/MIT/Software/Research_Software/mtag}"

if [[ ! -f "${CONDA_SH}" ]]; then
    echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
    exit 1
fi

if [[ ! -f "${MTAG_DIR}/mtag.py" ]]; then
    echo "ERROR: mtag.py not found at: ${MTAG_DIR}/mtag.py" >&2
    exit 1
fi

source "${CONDA_SH}"
conda activate "${MTAG_ENV_PATH}"

require_file() {
    local fpath="$1"
    if [[ ! -f "${fpath}" ]]; then
        echo "ERROR: Required file missing: ${fpath}" >&2
        exit 1
    fi
}

run_si_mtag() {
    local sumstats_root="$1"
    local pop_dir="$2"
    local output_dir="$3"
    local output_prefix="$4"
    local label="$5"

    local s1="${sumstats_root}/${pop_dir}/AbilityToConfide.Day_NoPCs.mtag.sumstats.txt"
    local s2="${sumstats_root}/${pop_dir}/FreqSoc.Day_NoPCs.mtag.sumstats.txt"
    local s3="${sumstats_root}/${pop_dir}/Loneliness.Day_NoPCs.mtag.sumstats.txt"

    require_file "${s1}"
    require_file "${s2}"
    require_file "${s3}"

    mkdir -p "${output_dir}"

    local sumstats_csv="${s1},${s2},${s3}"
    local out_prefix="${output_dir}/${output_prefix}"

    echo
    echo "============================================================"
    echo "Running MTAG: ${label}"
    echo "  Sumstats: ${sumstats_root}/${pop_dir}"
    echo "  Output:   ${out_prefix}"
    echo "============================================================"

    python "${MTAG_DIR}/mtag.py" \
        --sumstats "${sumstats_csv}" \
        --out "${out_prefix}" \
        --ld_ref_panel "${MTAG_DIR}/ld_ref_panel/eur_w_ld_chr/" \
        --snp_name snpid \
        --a1_name a1 \
        --a2_name a2 \
        --z_name z \
        --p_name pval \
        --n_name n \
        --chr_name chr \
        --bpos_name bpos \
        --eaf_name freq \
        --n_min 0.0 \
        --force \
        --stream_stdout
}

BINARY_SUMSTATS_ROOT="${REPODIR}/2_GWAS/mtag_results"
CONT_SUMSTATS_ROOT="${REPODIR}/2_GWAS/mtag_results_continuous"

# Binary analyses
run_si_mtag "${BINARY_SUMSTATS_ROOT}" "EUR_MM" \
    "${REPODIR}/3_MTAG/SI/EUR_MM/results" \
    "SI_EUR_MM_Output" \
    "SI binary EUR_MM"

run_si_mtag "${BINARY_SUMSTATS_ROOT}" "EUR_Male" \
    "${REPODIR}/3_MTAG/SI/EUR_Male_MM/results" \
    "SI_EUR_Male_MM_Output" \
    "SI binary EUR_Male"

run_si_mtag "${BINARY_SUMSTATS_ROOT}" "EUR_Female" \
    "${REPODIR}/3_MTAG/SI/EUR_Female_MM/results" \
    "SI_EUR_Female_MM_Output" \
    "SI binary EUR_Female"

# Continuous analyses
run_si_mtag "${CONT_SUMSTATS_ROOT}" "EUR_MM" \
    "${REPODIR}/3_MTAG/SI_continuous/EUR_MM/results" \
    "SI_EUR_MM_Output_continuous" \
    "SI continuous EUR_MM"

run_si_mtag "${CONT_SUMSTATS_ROOT}" "EUR_Male" \
    "${REPODIR}/3_MTAG/SI_continuous/EUR_Male_MM/results" \
    "SI_EUR_Male_MM_Output_continuous" \
    "SI continuous EUR_Male"

run_si_mtag "${CONT_SUMSTATS_ROOT}" "EUR_Female" \
    "${REPODIR}/3_MTAG/SI_continuous/EUR_Female_MM/results" \
    "SI_EUR_Female_MM_Output_continuous" \
    "SI continuous EUR_Female"

echo
echo "All SI MTAG analyses completed."
