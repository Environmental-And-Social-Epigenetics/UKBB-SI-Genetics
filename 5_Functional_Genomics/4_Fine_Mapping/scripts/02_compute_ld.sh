#!/usr/bin/env bash
set -euo pipefail

# Compute locus-specific LD matrices from a reference panel using PLINK.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

PLINK_BIN="${PLINK_BIN:-plink}"
BFILE_PREFIX="${BFILE_PREFIX:-${FG_ROOT}/reference_data/1000G_EUR/g1000_eur}"
LOCI_FILE="${LOCI_FILE:-${FG_ROOT}/4_Fine_Mapping/output/01_defined_loci/loci.loci.tsv}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/4_Fine_Mapping/output/02_ld_matrices}"
MAX_SNPS_PER_LOCUS="${MAX_SNPS_PER_LOCUS:-4000}"

mkdir -p "${OUT_DIR}"

if ! command -v "${PLINK_BIN}" >/dev/null 2>&1; then
    echo "ERROR: PLINK not found on PATH: ${PLINK_BIN}" >&2
    exit 1
fi

if [[ ! -f "${BFILE_PREFIX}.bed" || ! -f "${BFILE_PREFIX}.bim" || ! -f "${BFILE_PREFIX}.fam" ]]; then
    echo "ERROR: Reference PLINK files not found for prefix: ${BFILE_PREFIX}" >&2
    exit 1
fi

if [[ ! -f "${LOCI_FILE}" ]]; then
    echo "ERROR: Loci file not found: ${LOCI_FILE}" >&2
    exit 1
fi

run_count=0
skip_count=0

while IFS=$'\t' read -r locus_id chr start end lead_snp lead_bp lead_p method; do
    [[ "${locus_id}" == "locus_id" ]] && continue
    [[ -z "${locus_id}" ]] && continue

    out_prefix="${OUT_DIR}/${locus_id}"
    snplist="${out_prefix}.snplist"

    echo "[ld] ${locus_id} chr${chr}:${start}-${end}"
    "${PLINK_BIN}" \
        --bfile "${BFILE_PREFIX}" \
        --chr "${chr}" \
        --from-bp "${start}" \
        --to-bp "${end}" \
        --write-snplist \
        --out "${out_prefix}"

    if [[ ! -s "${snplist}" ]]; then
        echo "[skip] ${locus_id}: no SNPs found in reference panel."
        skip_count=$((skip_count + 1))
        continue
    fi

    snp_count="$(wc -l < "${snplist}" | tr -d ' ')"
    if (( snp_count < 2 )); then
        echo "[skip] ${locus_id}: fewer than 2 SNPs (${snp_count})."
        skip_count=$((skip_count + 1))
        continue
    fi
    if (( snp_count > MAX_SNPS_PER_LOCUS )); then
        echo "[skip] ${locus_id}: ${snp_count} SNPs exceeds MAX_SNPS_PER_LOCUS=${MAX_SNPS_PER_LOCUS}."
        skip_count=$((skip_count + 1))
        continue
    fi

    "${PLINK_BIN}" \
        --bfile "${BFILE_PREFIX}" \
        --extract "${snplist}" \
        --keep-allele-order \
        --r square gz yes-really \
        --out "${out_prefix}"

    run_count=$((run_count + 1))
done < "${LOCI_FILE}"

echo "LD matrix computation complete. Ran ${run_count} loci, skipped ${skip_count} loci."
