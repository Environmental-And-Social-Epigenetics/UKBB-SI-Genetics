#!/usr/bin/env python3
"""
Generate a Miami plot comparing binary vs continuous GWAS/MTAG results.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description="Create Miami plot from two summary statistics files.")
    parser.add_argument("--binary", type=Path, required=True, help="Binary-coded summary stats file.")
    parser.add_argument("--continuous", type=Path, required=True, help="Continuous-coded summary stats file.")
    parser.add_argument(
        "--out",
        type=Path,
        default=fg_root / "6_Visualization" / "output" / "miami_plot.png",
    )
    parser.add_argument("--title", default="Binary vs Continuous GWAS/MTAG Comparison")
    return parser.parse_args()


def best_col(cols: list[str], candidates: list[str]) -> str | None:
    for c in candidates:
        if c in cols:
            return c
    return None


def load_sumstats(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        df = pd.read_csv(path, sep=r"\s+")

    chr_col = best_col(list(df.columns), ["CHR", "chr"])
    bp_col = best_col(list(df.columns), ["BP", "bpos"])
    p_col = best_col(list(df.columns), ["mtag_pval", "pval", "P"])
    if None in (chr_col, bp_col, p_col):
        raise ValueError(f"Could not infer CHR/BP/P columns in {path}")

    out = pd.DataFrame(
        {
            "CHR": pd.to_numeric(df[chr_col], errors="coerce"),
            "BP": pd.to_numeric(df[bp_col], errors="coerce"),
            "P": pd.to_numeric(df[p_col], errors="coerce"),
        }
    ).dropna()
    out = out[(out["P"] > 0) & (out["P"] <= 1)]
    out["neglog10p"] = -np.log10(out["P"])
    out["CHR"] = out["CHR"].astype(int)
    out["BP"] = out["BP"].astype(int)
    return out


def add_cumulative_pos(df: pd.DataFrame, chr_offsets: dict[int, int]) -> pd.DataFrame:
    out = df.copy()
    out["x"] = out.apply(lambda r: r["BP"] + chr_offsets.get(int(r["CHR"]), 0), axis=1)
    return out


def main() -> None:
    args = parse_args()
    bin_df = load_sumstats(args.binary)
    cont_df = load_sumstats(args.continuous)

    chroms = sorted(set(bin_df["CHR"]).union(set(cont_df["CHR"])))
    chr_max = {}
    for c in chroms:
        max_bp = max(
            bin_df.loc[bin_df["CHR"] == c, "BP"].max() if (bin_df["CHR"] == c).any() else 0,
            cont_df.loc[cont_df["CHR"] == c, "BP"].max() if (cont_df["CHR"] == c).any() else 0,
        )
        chr_max[c] = int(max_bp)

    offsets = {}
    current = 0
    tick_pos = []
    tick_labels = []
    for c in chroms:
        offsets[c] = current
        tick_pos.append(current + chr_max[c] / 2)
        tick_labels.append(str(c))
        current += chr_max[c]

    bin_df = add_cumulative_pos(bin_df, offsets)
    cont_df = add_cumulative_pos(cont_df, offsets)
    cont_df["neglog10p"] = -cont_df["neglog10p"]  # bottom half

    fig, ax = plt.subplots(figsize=(15, 6))
    colors = ["#4C78A8", "#F58518"]

    for idx, c in enumerate(chroms):
        b = bin_df[bin_df["CHR"] == c]
        d = cont_df[cont_df["CHR"] == c]
        ax.scatter(b["x"], b["neglog10p"], s=4, color=colors[idx % 2], alpha=0.7)
        ax.scatter(d["x"], d["neglog10p"], s=4, color=colors[idx % 2], alpha=0.7)

    sig_line = -np.log10(5e-8)
    ax.axhline(sig_line, color="red", linestyle="--", linewidth=1)
    ax.axhline(-sig_line, color="red", linestyle="--", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.8)

    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(p)  [top: binary, bottom: continuous]")
    ax.set_title(args.title)
    ax.grid(alpha=0.2)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close()
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
