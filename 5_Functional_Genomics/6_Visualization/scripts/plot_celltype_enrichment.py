#!/usr/bin/env python3
"""
Plot cell-type-specific LDSC enrichment signals from *.cell_type_results.txt files.
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
    parser = argparse.ArgumentParser(description="Plot LDSC cell-type enrichments.")
    parser.add_argument(
        "--celltype-dir",
        type=Path,
        default=fg_root / "1_LDSC" / "output" / "celltype_h2",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=fg_root / "6_Visualization" / "output" / "celltype_enrichment.png",
    )
    parser.add_argument("--top-n", type=int, default=25)
    return parser.parse_args()


def parse_one(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=r"\s+")
    except Exception:
        return pd.DataFrame()

    name_col = None
    p_col = None
    for c in ["Name", "name", "Cell_type", "CellType"]:
        if c in df.columns:
            name_col = c
            break
    for c in ["Coefficient_P_value", "P", "P_value", "p"]:
        if c in df.columns:
            p_col = c
            break
    if name_col is None or p_col is None:
        return pd.DataFrame()

    out = df[[name_col, p_col]].copy()
    out.columns = ["cell_type", "pvalue"]
    out["pvalue"] = pd.to_numeric(out["pvalue"], errors="coerce")
    out = out.dropna()
    out = out[(out["pvalue"] > 0) & (out["pvalue"] <= 1)]
    out["neglog10p"] = -np.log10(out["pvalue"])
    out["trait"] = path.name.replace(".cell_type_results.txt", "").replace(".results", "")
    return out


def main() -> None:
    args = parse_args()
    files = sorted(args.celltype_dir.glob("*.cell_type_results.txt"))
    if not files:
        files = sorted(args.celltype_dir.glob("*.results"))
    if not files:
        raise FileNotFoundError(f"No LDSC CTS result files found in {args.celltype_dir}")

    all_df = [parse_one(f) for f in files]
    all_df = [df for df in all_df if not df.empty]
    if not all_df:
        raise RuntimeError("No parseable cell-type result tables found.")

    merged = pd.concat(all_df, ignore_index=True)
    top = merged.sort_values("neglog10p", ascending=False).head(args.top_n).copy()
    top = top.sort_values("neglog10p", ascending=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(11, max(6, 0.35 * len(top))))
    sns.barplot(data=top, x="neglog10p", y="cell_type", hue="trait", dodge=False)
    plt.xlabel("-log10(p)")
    plt.ylabel("Cell type / annotation")
    plt.title("LDSC Cell-Type Enrichment")
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close()
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
