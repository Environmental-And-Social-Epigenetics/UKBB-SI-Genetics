#!/usr/bin/env bash
set -euo pipefail

# Cross-tissue aggregation for S-PrediXcan outputs.
# Two modes:
#   1) Official S-MultiXcan (set USE_OFFICIAL_SMULTIXCAN=1 and configure args)
#   2) Fallback Stouffer aggregation across tissues (default).

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FG_ROOT="$(cd "${SCRIPTDIR}/../.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python}"
SPREDIXCAN_RESULTS_DIR="${SPREDIXCAN_RESULTS_DIR:-${FG_ROOT}/3_S-PrediXcan/output/spredixcan}"
OUT_DIR="${OUT_DIR:-${FG_ROOT}/3_S-PrediXcan/output/smultixcan}"

USE_OFFICIAL_SMULTIXCAN="${USE_OFFICIAL_SMULTIXCAN:-0}"
SMULTIXCAN_SCRIPT="${SMULTIXCAN_SCRIPT:-${FG_ROOT}/reference_data/MetaXcan/software/SMulTiXcan.py}"
MODELS_DIR="${MODELS_DIR:-${FG_ROOT}/reference_data/predixcan/models}"
MODELS_NAME_PATTERN="${MODELS_NAME_PATTERN:-^(.*)\\.db$}"
SNP_COVARIANCE="${SNP_COVARIANCE:-}"

mkdir -p "${OUT_DIR}"

if [[ ! -d "${SPREDIXCAN_RESULTS_DIR}" ]]; then
    echo "ERROR: SPrediXcan result directory not found: ${SPREDIXCAN_RESULTS_DIR}" >&2
    exit 1
fi

if [[ "${USE_OFFICIAL_SMULTIXCAN}" == "1" ]]; then
    if [[ ! -f "${SMULTIXCAN_SCRIPT}" ]]; then
        echo "ERROR: SMulTiXcan.py not found at ${SMULTIXCAN_SCRIPT}" >&2
        exit 1
    fi
    if [[ -z "${SNP_COVARIANCE}" ]]; then
        echo "ERROR: SNP_COVARIANCE must be set for official S-MultiXcan mode." >&2
        exit 1
    fi
    echo "[info] Running official S-MultiXcan mode."
    "${PYTHON_BIN}" "${SMULTIXCAN_SCRIPT}" \
        --models_folder "${MODELS_DIR}" \
        --models_name_pattern "${MODELS_NAME_PATTERN}" \
        --metaxcan_folder "${SPREDIXCAN_RESULTS_DIR}" \
        --metaxcan_filter ".*\\.csv" \
        --metaxcan_file_name_parse_pattern "^(.*)__.*\\.csv$" \
        --snp_covariance "${SNP_COVARIANCE}" \
        --cutoff_condition_number 30 \
        --verbosity 9 \
        --throw \
        --output "${OUT_DIR}/smultixcan_official.tsv"
    echo "Completed official S-MultiXcan."
    exit 0
fi

echo "[info] Running fallback cross-tissue Stouffer aggregation."
"${PYTHON_BIN}" - "${SPREDIXCAN_RESULTS_DIR}" "${OUT_DIR}" <<'PY'
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

in_dir = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)

files = sorted(in_dir.glob("*.csv"))
if not files:
    raise SystemExit(f"No CSV files found in {in_dir}")

grouped = {}
for f in files:
    stem = f.stem
    if "__" not in stem:
        continue
    gwas_id, _tissue = stem.rsplit("__", 1)
    grouped.setdefault(gwas_id, []).append(f)

if not grouped:
    raise SystemExit("No grouped S-PrediXcan files found with '<gwas>__<tissue>.csv' pattern.")

def pick_column(columns, candidates):
    for c in candidates:
        if c in columns:
            return c
    return None

for gwas_id, gwas_files in grouped.items():
    gene_to_z = {}
    gene_to_absz = {}
    gene_to_p = {}
    for f in gwas_files:
        try:
            df = pd.read_csv(f)
        except Exception:
            # Some runs write tab-delimited text.
            df = pd.read_csv(f, sep="\t")

        gene_col = pick_column(df.columns, ["gene", "gene_name", "ensg"])
        z_col = pick_column(df.columns, ["zscore", "z"])
        p_col = pick_column(df.columns, ["pvalue", "p"])
        if gene_col is None or z_col is None:
            continue

        sub = df[[gene_col, z_col] + ([p_col] if p_col else [])].copy()
        sub = sub.dropna(subset=[gene_col, z_col])
        sub[gene_col] = sub[gene_col].astype(str).str.upper()
        sub[z_col] = pd.to_numeric(sub[z_col], errors="coerce")
        if p_col:
            sub[p_col] = pd.to_numeric(sub[p_col], errors="coerce")
        sub = sub[np.isfinite(sub[z_col])]

        for row in sub.itertuples(index=False):
            gene = getattr(row, gene_col)
            z = float(getattr(row, z_col))
            p = float(getattr(row, p_col)) if p_col else np.nan
            gene_to_z.setdefault(gene, []).append(z)
            gene_to_absz.setdefault(gene, []).append(abs(z))
            if np.isfinite(p):
                gene_to_p.setdefault(gene, []).append(p)

    rows = []
    for gene, zvals in gene_to_z.items():
        k = len(zvals)
        if k == 0:
            continue
        stouffer_z = float(np.sum(zvals) / math.sqrt(k))
        stouffer_p = math.erfc(abs(stouffer_z) / math.sqrt(2.0))
        rows.append(
            {
                "gene": gene,
                "n_tissues": k,
                "stouffer_z": stouffer_z,
                "stouffer_p": stouffer_p,
                "mean_abs_z": float(np.mean(gene_to_absz.get(gene, [np.nan]))),
                "min_tissue_p": float(np.min(gene_to_p[gene])) if gene in gene_to_p else np.nan,
            }
        )

    out_file = out_dir / f"{gwas_id}__smultixcan_fallback.tsv"
    out_df = pd.DataFrame(rows)
    if out_df.empty:
        out_df = pd.DataFrame(columns=["gene", "n_tissues", "stouffer_z", "stouffer_p", "mean_abs_z", "min_tissue_p"])
    else:
        out_df = out_df.sort_values("stouffer_p", ascending=True)
    out_df.to_csv(out_file, sep="\t", index=False)
    print(f"[write] {out_file} ({len(out_df)} genes)")
PY

echo "Completed fallback cross-tissue aggregation."
