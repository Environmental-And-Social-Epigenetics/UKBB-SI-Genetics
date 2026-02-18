#!/usr/bin/env bash
set -euo pipefail

# Run MAGMA annotation and gene analysis for all sex-stratified SI files.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

MAGMA_BIN="${MAGMA_BIN:-${FG_ROOT}/reference_data/magma/bin/magma}"
GENE_LOC_FILE="${GENE_LOC_FILE:-${FG_ROOT}/reference_data/magma/NCBI37.3.gene.loc}"
BFILE_PREFIX="${BFILE_PREFIX:-${FG_ROOT}/reference_data/1000G_EUR/g1000_eur}"
MAGMA_INPUT_DIR="${MAGMA_INPUT_DIR:-${FG_ROOT}/0_munge/output/magma}"
OUT_BASE="${OUT_BASE:-${FG_ROOT}/7_Sex_Differences/output/05_magma_sex_stratified}"

ANNOT_DIR="${OUT_BASE}/annotation"
GENE_DIR="${OUT_BASE}/gene"
mkdir -p "${ANNOT_DIR}" "${GENE_DIR}"

if [[ ! -x "${MAGMA_BIN}" ]]; then
    echo "ERROR: MAGMA binary not executable: ${MAGMA_BIN}" >&2
    exit 1
fi
if [[ ! -f "${GENE_LOC_FILE}" ]]; then
    echo "ERROR: Gene location file missing: ${GENE_LOC_FILE}" >&2
    exit 1
fi
if [[ ! -f "${BFILE_PREFIX}.bed" || ! -f "${BFILE_PREFIX}.bim" || ! -f "${BFILE_PREFIX}.fam" ]]; then
    echo "ERROR: Missing merged 1000G PLINK files for prefix ${BFILE_PREFIX}" >&2
    exit 1
fi

trait_name_from_idx() {
    local trait_idx="$1"
    if [[ "${trait_idx}" == "1" ]]; then
        echo "AbilityToConfide"
    elif [[ "${trait_idx}" == "2" ]]; then
        echo "FreqSoc"
    elif [[ "${trait_idx}" == "3" ]]; then
        echo "Loneliness"
    else
        echo "UnknownTrait"
    fi
}

run_one() {
    local coding="$1"
    local sex="$2"
    local trait_idx="$3"
    local trait_name="$4"
    local input_file="${MAGMA_INPUT_DIR}/${coding}__EUR_${sex}_MM__trait${trait_idx}__${trait_name}.genes.raw"
    local base="${coding}__EUR_${sex}_MM__trait${trait_idx}__${trait_name}"
    local annot_prefix="${ANNOT_DIR}/${base}"
    local gene_prefix="${GENE_DIR}/${base}"

    if [[ ! -f "${input_file}" ]]; then
        echo "[skip] Missing input: ${input_file}"
        return 0
    fi
    if [[ -f "${gene_prefix}.genes.out" ]]; then
        echo "[skip] Already complete: ${base}"
        return 0
    fi

    echo "[magma] ${base} -- annotate"
    "${MAGMA_BIN}" \
        --annotate window=10,10 \
        --snp-loc "${input_file}" \
        --gene-loc "${GENE_LOC_FILE}" \
        --out "${annot_prefix}"

    echo "[magma] ${base} -- gene analysis"
    "${MAGMA_BIN}" \
        --bfile "${BFILE_PREFIX}" \
        --pval "${input_file}" use=SNP,P ncol=N \
        --gene-annot "${annot_prefix}.genes.annot" \
        --out "${gene_prefix}"
}

for coding in SI SI_continuous; do
    for sex in Male Female; do
        for trait_idx in 1 2 3; do
            run_one "${coding}" "${sex}" "${trait_idx}" "$(trait_name_from_idx "${trait_idx}")"
        done
    done
done

echo "MAGMA sex-stratified analyses complete."
