#!/usr/bin/env bash
set -euo pipefail

# Compute pairwise genetic correlations within each analysis group/population.
# Expected filename pattern in SUMSTATS_DIR:
#   <analysis>__<population>__trait<k>__<trait>.sumstats.gz

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

SUMSTATS_DIR="${SUMSTATS_DIR:-${FG_ROOT}/0_munge/output/ldsc_munged}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/1_LDSC/output/rg_internal}"
PYTHON_BIN="${PYTHON_BIN:-python}"
LDSC_SCRIPT="${LDSC_SCRIPT:-${FG_ROOT}/reference_data/ldsc/ldsc/ldsc.py}"
REF_LD_CHR="${REF_LD_CHR:-${FG_ROOT}/reference_data/ldsc/eur_w_ld_chr/}"
W_LD_CHR="${W_LD_CHR:-${FG_ROOT}/reference_data/ldsc/weights_hm3_no_hla/weights.hm3_noMHC.}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${LDSC_SCRIPT}" ]]; then
    echo "ERROR: Could not find ldsc.py at ${LDSC_SCRIPT}" >&2
    exit 1
fi

if [[ ! -d "${SUMSTATS_DIR}" ]]; then
    echo "ERROR: SUMSTATS_DIR not found: ${SUMSTATS_DIR}" >&2
    exit 1
fi

pair_tsv="${OUT_DIR}/_internal_pairs.tsv"
"${PYTHON_BIN}" - "${SUMSTATS_DIR}" "${pair_tsv}" <<'PY'
import itertools
import sys
from pathlib import Path

sumstats_dir = Path(sys.argv[1])
pair_path = Path(sys.argv[2])

files = sorted(sumstats_dir.glob("*.sumstats.gz"))
groups = {}
for f in files:
    stem = f.name[:-len(".sumstats.gz")]
    parts = stem.split("__")
    if len(parts) < 4:
        continue
    group_key = "__".join(parts[:2])
    groups.setdefault(group_key, []).append((stem, f))

rows = []
for group_key, group_files in groups.items():
    if len(group_files) < 2:
        continue
    for (name_a, path_a), (name_b, path_b) in itertools.combinations(group_files, 2):
        label = f"{group_key}__{name_a.split('__')[-1]}_vs_{name_b.split('__')[-1]}"
        rows.append((str(path_a), str(path_b), label))

with pair_path.open("w") as handle:
    for row in rows:
        handle.write("\t".join(row) + "\n")

print(len(rows))
PY

if [[ ! -s "${pair_tsv}" ]]; then
    echo "ERROR: No internal trait pairs could be built from ${SUMSTATS_DIR}" >&2
    exit 1
fi

pair_count=0
while IFS=$'\t' read -r sumstats_a sumstats_b label; do
    out_prefix="${OUT_DIR}/${label}"
    echo "[rg-internal] ${label}"
    "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
        --rg "${sumstats_a},${sumstats_b}" \
        --ref-ld-chr "${REF_LD_CHR}" \
        --w-ld-chr "${W_LD_CHR}" \
        --out "${out_prefix}"
    pair_count=$((pair_count + 1))
done < "${pair_tsv}"

echo "Completed ${pair_count} internal rg runs."
