#!/usr/bin/env python3
"""
Create regional association plots for top loci.
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
    parser = argparse.ArgumentParser(description="Plot regional association signals for lead loci.")
    parser.add_argument("--sumstats", type=Path, required=True, help="GWAS/MTAG summary stats file.")
    parser.add_argument(
        "--loci",
        type=Path,
        default=fg_root / "4_Fine_Mapping" / "output" / "01_defined_loci" / "loci.loci.tsv",
    )
    parser.add_argument(
        "--pip",
        type=Path,
        default=fg_root / "4_Fine_Mapping" / "output" / "03_susie" / "susie_pip.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=fg_root / "6_Visualization" / "output" / "regional_loci.png",
    )
    parser.add_argument("--top-n", type=int, default=6)
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

    snp_col = best_col(list(df.columns), ["SNP", "snpid"])
    chr_col = best_col(list(df.columns), ["CHR", "chr"])
    bp_col = best_col(list(df.columns), ["BP", "bpos"])
    p_col = best_col(list(df.columns), ["mtag_pval", "pval", "P"])
    if None in (snp_col, chr_col, bp_col, p_col):
        raise ValueError("Could not infer SNP/CHR/BP/P columns from sumstats.")

    out = pd.DataFrame(
        {
            "SNP": df[snp_col].astype(str),
            "CHR": pd.to_numeric(df[chr_col], errors="coerce"),
            "BP": pd.to_numeric(df[bp_col], errors="coerce"),
            "P": pd.to_numeric(df[p_col], errors="coerce"),
        }
    ).dropna()
    out = out[(out["P"] > 0) & (out["P"] <= 1)]
    out["neglog10p"] = -np.log10(out["P"])
    return out


def main() -> None:
    args = parse_args()
    sumstats = load_sumstats(args.sumstats)
    loci = pd.read_csv(args.loci, sep="\t")
    if loci.empty:
        raise RuntimeError(f"No loci found in {args.loci}")
    loci = loci.sort_values("LEAD_P" if "LEAD_P" in loci.columns else "locus_id").head(args.top_n)

    pip_df = None
    if args.pip.exists():
        pip_df = pd.read_csv(args.pip, sep="\t")
        if "SNP" in pip_df.columns and "PIP" in pip_df.columns:
            pip_df = pip_df[["SNP", "PIP"]].copy()
        else:
            pip_df = None

    n = len(loci)
    fig, axes = plt.subplots(n, 1, figsize=(11, max(3 * n, 4)), sharex=False)
    if n == 1:
        axes = [axes]

    for ax, (_, locus) in zip(axes, loci.iterrows()):
        locus_id = locus["locus_id"]
        chr_val = int(locus["CHR"])
        start = int(locus["START"])
        end = int(locus["END"])
        lead_snp = locus.get("LEAD_SNP", "NA")

        sub = sumstats[(sumstats["CHR"] == chr_val) & (sumstats["BP"] >= start) & (sumstats["BP"] <= end)].copy()
        if sub.empty:
            ax.text(0.5, 0.5, f"{locus_id}: no SNPs", ha="center", va="center")
            ax.axis("off")
            continue

        if pip_df is not None:
            sub = sub.merge(pip_df, on="SNP", how="left")
            sc = ax.scatter(
                sub["BP"],
                sub["neglog10p"],
                c=sub["PIP"].fillna(0.0),
                s=12 + 40 * sub["PIP"].fillna(0.0),
                cmap="viridis",
                alpha=0.8,
                edgecolor="none",
            )
            cbar = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01)
            cbar.set_label("PIP")
        else:
            ax.scatter(sub["BP"], sub["neglog10p"], s=10, alpha=0.7, color="#1f77b4")

        ax.axhline(-np.log10(5e-8), color="red", linestyle="--", linewidth=1)
        ax.set_title(f"{locus_id} | chr{chr_val}:{start:,}-{end:,} | lead={lead_snp}")
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("-log10(p)")
        ax.grid(alpha=0.2)

    plt.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, dpi=300)
    plt.close()
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
