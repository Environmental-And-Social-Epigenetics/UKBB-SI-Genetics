#!/usr/bin/env python3
"""
Integrate MAGMA, S-PrediXcan/S-MultiXcan, and coloc evidence into
a single ranked gene-prioritization table.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description="Prioritize genes across post-GWAS analyses.")
    parser.add_argument(
        "--magma-dir",
        type=Path,
        default=fg_root / "2_MAGMA" / "output" / "02_gene_analysis",
        help="Directory containing MAGMA gene outputs (*.genes.out).",
    )
    parser.add_argument(
        "--spredixcan-dir",
        type=Path,
        default=fg_root / "3_S-PrediXcan" / "output" / "spredixcan",
        help="Directory containing S-PrediXcan tissue outputs (*.csv).",
    )
    parser.add_argument(
        "--smultixcan-dir",
        type=Path,
        default=fg_root / "3_S-PrediXcan" / "output" / "smultixcan",
        help="Directory containing S-MultiXcan outputs (*.tsv).",
    )
    parser.add_argument(
        "--coloc-summary",
        type=Path,
        default=fg_root / "5_Colocalization" / "output" / "summary" / "coloc_prioritized.tsv",
        help="Optional coloc summary table.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=fg_root / "3_S-PrediXcan" / "output" / "prioritized_genes.tsv",
        help="Output TSV path.",
    )
    parser.add_argument("--magma-alpha", type=float, default=2.7e-6)
    parser.add_argument("--spredixcan-alpha", type=float, default=2.5e-6)
    parser.add_argument("--smultixcan-alpha", type=float, default=2.5e-6)
    parser.add_argument("--coloc-threshold", type=float, default=0.8)
    return parser.parse_args()


def _best_column(columns: list[str], candidates: list[str]) -> str | None:
    for c in candidates:
        if c in columns:
            return c
    return None


def _safe_read(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)


def load_magma(magma_dir: Path) -> pd.DataFrame:
    rows = []
    files = sorted(magma_dir.glob("*.genes.out"))
    for f in files:
        df = pd.read_csv(f, delim_whitespace=True, low_memory=False)
        gene_col = _best_column(list(df.columns), ["GENE", "Gene", "gene"])
        p_col = _best_column(list(df.columns), ["P", "p"])
        if gene_col is None or p_col is None:
            continue

        trait = f.name.replace(".genes.out", "")
        sub = df[[gene_col, p_col]].copy()
        sub[gene_col] = sub[gene_col].astype(str).str.upper()
        sub[p_col] = pd.to_numeric(sub[p_col], errors="coerce")
        sub = sub.dropna()
        sub["trait"] = trait
        sub = sub.sort_values(p_col, ascending=True).drop_duplicates(gene_col, keep="first")

        for row in sub.itertuples(index=False):
            rows.append(
                {
                    "gene": getattr(row, gene_col),
                    "magma_p": float(getattr(row, p_col)),
                    "magma_trait": getattr(row, "trait"),
                }
            )

    if not rows:
        return pd.DataFrame(columns=["gene", "magma_p", "magma_trait"])
    out = pd.DataFrame(rows).sort_values("magma_p", ascending=True).drop_duplicates("gene", keep="first")
    return out


def load_spredixcan(sp_dir: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(sp_dir.glob("*.csv")):
        try:
            df = pd.read_csv(f)
        except Exception:
            df = pd.read_csv(f, sep="\t")

        gene_col = _best_column(list(df.columns), ["gene", "gene_name", "ensg"])
        p_col = _best_column(list(df.columns), ["pvalue", "p"])
        z_col = _best_column(list(df.columns), ["zscore", "z"])
        if gene_col is None or p_col is None:
            continue

        stem = f.stem
        tissue = stem.rsplit("__", 1)[-1] if "__" in stem else "unknown_tissue"
        sub = df[[gene_col, p_col] + ([z_col] if z_col else [])].copy()
        sub[gene_col] = sub[gene_col].astype(str).str.upper()
        sub[p_col] = pd.to_numeric(sub[p_col], errors="coerce")
        if z_col:
            sub[z_col] = pd.to_numeric(sub[z_col], errors="coerce")
        sub = sub.dropna(subset=[gene_col, p_col])

        for row in sub.itertuples(index=False):
            rows.append(
                {
                    "gene": getattr(row, gene_col),
                    "spredixcan_p": float(getattr(row, p_col)),
                    "spredixcan_z": float(getattr(row, z_col)) if z_col else np.nan,
                    "best_tissue_candidate": tissue,
                }
            )

    if not rows:
        return pd.DataFrame(
            columns=["gene", "spredixcan_min_p", "spredixcan_best_tissue", "spredixcan_max_abs_z", "spredixcan_n_tissues"]
        )

    df = pd.DataFrame(rows)
    df["abs_z"] = df["spredixcan_z"].abs()
    agg = (
        df.sort_values("spredixcan_p", ascending=True)
        .groupby("gene", as_index=False)
        .agg(
            spredixcan_min_p=("spredixcan_p", "min"),
            spredixcan_best_tissue=("best_tissue_candidate", "first"),
            spredixcan_max_abs_z=("abs_z", "max"),
            spredixcan_n_tissues=("best_tissue_candidate", "nunique"),
        )
    )
    return agg


def load_smultixcan(smx_dir: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(smx_dir.glob("*.tsv")):
        df = _safe_read(f)
        gene_col = _best_column(list(df.columns), ["gene", "gene_name", "ensg"])
        p_col = _best_column(list(df.columns), ["pvalue", "stouffer_p", "p"])
        z_col = _best_column(list(df.columns), ["zscore", "stouffer_z", "z"])
        if gene_col is None or p_col is None:
            continue

        sub = df[[gene_col, p_col] + ([z_col] if z_col else [])].copy()
        sub[gene_col] = sub[gene_col].astype(str).str.upper()
        sub[p_col] = pd.to_numeric(sub[p_col], errors="coerce")
        if z_col:
            sub[z_col] = pd.to_numeric(sub[z_col], errors="coerce")
        sub = sub.dropna(subset=[gene_col, p_col])

        for row in sub.itertuples(index=False):
            rows.append(
                {
                    "gene": getattr(row, gene_col),
                    "smultixcan_p": float(getattr(row, p_col)),
                    "smultixcan_z": float(getattr(row, z_col)) if z_col else np.nan,
                }
            )

    if not rows:
        return pd.DataFrame(columns=["gene", "smultixcan_p", "smultixcan_abs_z"])

    out = pd.DataFrame(rows)
    out["smultixcan_abs_z"] = out["smultixcan_z"].abs()
    out = (
        out.sort_values("smultixcan_p", ascending=True)
        .drop_duplicates("gene", keep="first")
        [["gene", "smultixcan_p", "smultixcan_abs_z"]]
    )
    return out


def load_coloc(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["gene", "coloc_pp_h4"])

    df = _safe_read(path)
    gene_col = _best_column(list(df.columns), ["gene", "gene_name", "ensg"])
    h4_col = _best_column(list(df.columns), ["PP.H4", "PP_H4", "pp_h4", "coloc_pp_h4"])
    if gene_col is None or h4_col is None:
        return pd.DataFrame(columns=["gene", "coloc_pp_h4"])

    out = df[[gene_col, h4_col]].copy()
    out[gene_col] = out[gene_col].astype(str).str.upper()
    out[h4_col] = pd.to_numeric(out[h4_col], errors="coerce")
    out = out.dropna().sort_values(h4_col, ascending=False).drop_duplicates(gene_col, keep="first")
    out = out.rename(columns={gene_col: "gene", h4_col: "coloc_pp_h4"})
    return out


def score_row(row: pd.Series) -> float:
    score = 0.0
    for pcol, weight in [
        ("magma_p", 1.0),
        ("spredixcan_min_p", 1.0),
        ("smultixcan_p", 1.25),
    ]:
        pval = row.get(pcol, np.nan)
        if pd.notna(pval) and pval > 0:
            score += weight * (-math.log10(float(pval)))

    h4 = row.get("coloc_pp_h4", np.nan)
    if pd.notna(h4):
        score += 5.0 * float(h4)
    return score


def main() -> None:
    args = parse_args()

    magma = load_magma(args.magma_dir)
    spredixcan = load_spredixcan(args.spredixcan_dir)
    smultixcan = load_smultixcan(args.smultixcan_dir)
    coloc = load_coloc(args.coloc_summary)

    merged = pd.DataFrame({"gene": sorted(set(magma["gene"]) | set(spredixcan["gene"]) | set(smultixcan["gene"]) | set(coloc["gene"]))})
    merged = merged.merge(magma, on="gene", how="left")
    merged = merged.merge(spredixcan, on="gene", how="left")
    merged = merged.merge(smultixcan, on="gene", how="left")
    merged = merged.merge(coloc, on="gene", how="left")

    merged["magma_sig"] = merged["magma_p"] < args.magma_alpha
    merged["spredixcan_sig"] = merged["spredixcan_min_p"] < args.spredixcan_alpha
    merged["smultixcan_sig"] = merged["smultixcan_p"] < args.smultixcan_alpha
    merged["coloc_sig"] = merged["coloc_pp_h4"] >= args.coloc_threshold
    merged["convergent_hits"] = (
        merged["magma_sig"].fillna(False).astype(int)
        + merged["spredixcan_sig"].fillna(False).astype(int)
        + merged["smultixcan_sig"].fillna(False).astype(int)
        + merged["coloc_sig"].fillna(False).astype(int)
    )
    merged["prioritization_score"] = merged.apply(score_row, axis=1)
    merged = merged.sort_values(
        ["convergent_hits", "prioritization_score"],
        ascending=[False, False],
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote prioritized genes: {args.out}")


if __name__ == "__main__":
    main()
