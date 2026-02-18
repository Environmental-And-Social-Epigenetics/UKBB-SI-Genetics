#!/usr/bin/env python3
"""
Plot MAGMA tissue-expression enrichment results from gene-property analysis.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description="Plot MAGMA tissue enrichments.")
    parser.add_argument(
        "--magma-tissue-dir",
        type=Path,
        default=fg_root / "2_MAGMA" / "output" / "04_tissue_expression",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=fg_root / "6_Visualization" / "output" / "magma_tissue_enrichment.png",
    )
    parser.add_argument("--top-n", type=int, default=20)
    return parser.parse_args()


def parse_file(path: Path) -> pd.DataFrame:
    # MAGMA outputs can be whitespace-delimited with varying headers.
    try:
        df = pd.read_csv(path, sep=r"\s+")
    except Exception:
        return pd.DataFrame()

    var_col = None
    p_col = None
    for c in ["VARIABLE", "Variable", "variable"]:
        if c in df.columns:
            var_col = c
            break
    for c in ["P", "P_VALUE", "p", "pval"]:
        if c in df.columns:
            p_col = c
            break
    if var_col is None or p_col is None:
        return pd.DataFrame()

    out = df[[var_col, p_col]].copy()
    out.columns = ["tissue", "pvalue"]
    out["pvalue"] = pd.to_numeric(out["pvalue"], errors="coerce")
    out = out.dropna()
    out = out[(out["pvalue"] > 0) & (out["pvalue"] <= 1)]
    out["neglog10p"] = -np.log10(out["pvalue"])
    out["trait"] = path.name.replace(".gsa.out", "").replace(".out", "")
    return out


def main() -> None:
    args = parse_args()
    files = sorted(args.magma_tissue_dir.glob("*.gsa.out"))
    if not files:
        files = sorted(args.magma_tissue_dir.glob("*.out"))
    if not files:
        raise FileNotFoundError(f"No MAGMA tissue output files found in {args.magma_tissue_dir}")

    parsed = [parse_file(f) for f in files]
    parsed = [p for p in parsed if not p.empty]
    if not parsed:
        raise RuntimeError("No parseable MAGMA tissue output tables found.")

    merged = pd.concat(parsed, ignore_index=True)
    top = merged.sort_values("neglog10p", ascending=False).head(args.top_n).copy()
    top = top.sort_values("neglog10p", ascending=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(11, max(6, 0.35 * len(top))))
    sns.barplot(data=top, x="neglog10p", y="tissue", hue="trait", dodge=False)
    plt.xlabel("-log10(p)")
    plt.ylabel("Tissue")
    plt.title("MAGMA Tissue Enrichment (Top Signals)")
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close()
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
