#!/usr/bin/env bash
set -euo pipefail

# Cross-sex LDSC genetic correlation for all SI traits and codings.

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

LDSC_SCRIPT="${LDSC_SCRIPT:-${FG_ROOT}/reference_data/ldsc/ldsc/ldsc.py}"
SUMSTATS_DIR="${SUMSTATS_DIR:-${FG_ROOT}/0_munge/output/ldsc_munged}"
REF_LD_CHR="${REF_LD_CHR:-${FG_ROOT}/reference_data/ldsc/LDscore/LDscore.}"
W_LD_CHR="${W_LD_CHR:-${FG_ROOT}/reference_data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.}"
PYTHON_BIN="${PYTHON_BIN:-python}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/7_Sex_Differences/output/01_cross_sex_rg}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${LDSC_SCRIPT}" ]]; then
    echo "ERROR: LDSC script not found: ${LDSC_SCRIPT}" >&2
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

run_rg() {
    local coding="$1"
    local trait_idx="$2"
    local trait_name="$3"

    local male_file="${SUMSTATS_DIR}/${coding}__EUR_Male_MM__trait${trait_idx}__${trait_name}.sumstats.gz"
    local female_file="${SUMSTATS_DIR}/${coding}__EUR_Female_MM__trait${trait_idx}__${trait_name}.sumstats.gz"
    local out_prefix="${OUT_DIR}/${coding}__trait${trait_idx}__${trait_name}__Male_vs_Female"

    if [[ ! -f "${male_file}" ]]; then
        echo "ERROR: Missing male file: ${male_file}" >&2
        return 1
    fi
    if [[ ! -f "${female_file}" ]]; then
        echo "ERROR: Missing female file: ${female_file}" >&2
        return 1
    fi

    echo "[cross-sex-rg] ${coding} trait${trait_idx} ${trait_name}"
    "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
        --rg "${male_file},${female_file}" \
        --ref-ld-chr "${REF_LD_CHR}" \
        --w-ld-chr "${W_LD_CHR}" \
        --out "${out_prefix}"
}

for coding in SI SI_continuous; do
    for trait_idx in 1 2 3; do
        run_rg "${coding}" "${trait_idx}" "$(trait_name_from_idx "${trait_idx}")"
    done
done

summary_tsv="${OUT_DIR}/cross_sex_rg_summary.tsv"
"${PYTHON_BIN}" - "${OUT_DIR}" "${summary_tsv}" <<'PY'
import re
import sys
from pathlib import Path

out_dir = Path(sys.argv[1])
summary_path = Path(sys.argv[2])
log_files = sorted(out_dir.glob("*.log"))

rows = []
for log_file in log_files:
    text = log_file.read_text()
    rg = se = z = p = None
    rg_match = re.search(r"Genetic Correlation:\s+([-\d.eE]+)\s+\(([-\d.eE]+)\)", text)
    z_match = re.search(r"Z-score:\s+([-\d.eE]+)", text)
    p_match = re.search(r"P:\s+([-\d.eE]+)", text)
    if rg_match:
        rg = float(rg_match.group(1))
        se = float(rg_match.group(2))
    if z_match:
        z = float(z_match.group(1))
    if p_match:
        p = float(p_match.group(1))

    stem = log_file.stem
    # Example: SI__trait3__Loneliness__Male_vs_Female
    parts = stem.split("__")
    coding = parts[0] if len(parts) > 0 else "unknown"
    trait = parts[2] if len(parts) > 2 else "unknown"
    trait_idx = parts[1].replace("trait", "") if len(parts) > 1 else "NA"
    rows.append((coding, trait_idx, trait, rg, se, z, p, str(log_file)))

with summary_path.open("w") as handle:
    handle.write("coding\ttrait_index\ttrait\trg\tse\tz\tp\tlog_file\n")
    for row in rows:
        handle.write("\t".join("" if x is None else str(x) for x in row) + "\n")

print(f"Wrote {len(rows)} rows: {summary_path}")
PY

echo "Cross-sex rg complete."
