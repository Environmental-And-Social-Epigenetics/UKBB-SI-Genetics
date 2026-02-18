#!/usr/bin/env python3
"""
Plot top partitioned heritability enrichments from LDSC output.
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
    parser = argparse.ArgumentParser(description="Plot partitioned h2 enrichments.")
    parser.add_argument(
        "--partitioned-dir",
        type=Path,
        default=fg_root / "1_LDSC" / "output" / "partitioned_h2",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=fg_root / "6_Visualization" / "output" / "partitioned_h2_top_categories.png",
    )
    parser.add_argument("--top-n", type=int, default=15)
    return parser.parse_args()


def parse_partitioned_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+")
    required_any = {"Category", "Enrichment"}
    if not required_any.issubset(df.columns):
        return pd.DataFrame()

    p_col = "Enrichment_p"
    if p_col not in df.columns:
        p_col = "Coefficient_z-score" if "Coefficient_z-score" in df.columns else None

    out = df.copy()
    out["trait"] = path.stem.replace(".results", "")
    if p_col is not None:
        out["plot_p"] = pd.to_numeric(out[p_col], errors="coerce")
    else:
        out["plot_p"] = np.nan
    out["Enrichment"] = pd.to_numeric(out["Enrichment"], errors="coerce")
    out = out.dropna(subset=["Category", "Enrichment"])
    return out


def main() -> None:
    args = parse_args()
    files = sorted(args.partitioned_dir.glob("*.results"))
    if not files:
        raise FileNotFoundError(f"No *.results files found in {args.partitioned_dir}")

    all_rows = []
    for f in files:
        df = parse_partitioned_file(f)
        if not df.empty:
            all_rows.append(df)
    if not all_rows:
        raise RuntimeError("No parseable partitioned heritability results found.")

    full = pd.concat(all_rows, ignore_index=True)
    full["score"] = np.where(
        full["plot_p"].notna() & (full["plot_p"] > 0),
        -np.log10(full["plot_p"]),
        full["Enrichment"],
    )
    full = full.sort_values("score", ascending=False)
    top = full.head(args.top_n).copy()
    top = top.sort_values("score", ascending=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(11, max(6, 0.35 * len(top))))
    sns.barplot(data=top, x="score", y="Category", hue="trait", dodge=False)
    plt.xlabel("-log10(p) (or enrichment if p unavailable)")
    plt.ylabel("Annotation category")
    plt.title("Top Partitioned Heritability Signals")
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close()
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
