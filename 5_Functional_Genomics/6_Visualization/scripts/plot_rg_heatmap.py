#!/usr/bin/env python3
"""
Plot an LDSC genetic-correlation heatmap from *.results files.
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
    parser = argparse.ArgumentParser(description="Generate heatmap from LDSC rg results.")
    parser.add_argument(
        "--rg-dir",
        type=Path,
        default=fg_root / "1_LDSC" / "output" / "rg_external",
        help="Directory containing LDSC rg *.results files.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=fg_root / "6_Visualization" / "output" / "rg_heatmap.png",
    )
    parser.add_argument("--min-abs-rg", type=float, default=0.0)
    return parser.parse_args()


def tidy_trait_name(path_or_name: str) -> str:
    p = Path(str(path_or_name))
    stem = p.stem
    stem = stem.replace(".sumstats", "")
    return stem


def read_results(path: Path) -> tuple[str, str, float] | None:
    # Try tabular .results format first
    try:
        df = pd.read_csv(path, sep=r"\s+")
        if not df.empty and {"p1", "p2", "rg"}.issubset(df.columns):
            row = df.iloc[0]
            p1 = tidy_trait_name(row["p1"])
            p2 = tidy_trait_name(row["p2"])
            rg = pd.to_numeric(row["rg"], errors="coerce")
            if pd.notna(rg):
                return p1, p2, float(rg)
    except Exception:
        pass

    # Fall back to parsing LDSC .log files
    import re
    try:
        text = path.read_text()
    except Exception:
        return None

    rg_match = re.search(r"Genetic Correlation:\s+([-.e\d]+)", text)
    if rg_match is None:
        return None
    rg = float(rg_match.group(1))

    # Extract trait names from filename
    stem = path.stem
    if "__vs__" in stem:
        parts = stem.split("__vs__")
        p1, p2 = tidy_trait_name(parts[0]), tidy_trait_name(parts[1])
    elif "_vs_" in stem:
        idx = stem.index("_vs_")
        group_and_t1 = stem[:idx]
        t2 = stem[idx + 4:]
        last_sep = group_and_t1.rfind("__")
        if last_sep > 0:
            p1 = group_and_t1[last_sep + 2:]
            p2 = t2
            group_prefix = group_and_t1[:last_sep]
            p1 = f"{group_prefix}__{p1}"
            p2 = f"{group_prefix}__{p2}"
        else:
            p1 = group_and_t1
            p2 = t2
    else:
        p1 = stem
        p2 = stem
    return p1, p2, rg


def main() -> None:
    args = parse_args()
    files = sorted(args.rg_dir.glob("*.results"))
    if not files:
        files = sorted(args.rg_dir.glob("*.log"))
    if not files:
        raise FileNotFoundError(f"No *.results or *.log files found in {args.rg_dir}")

    rows = []
    for f in files:
        val = read_results(f)
        if val is None:
            continue
        p1, p2, rg = val
        if abs(rg) < args.min_abs_rg:
            continue
        rows.append((p1, p2, rg))
        rows.append((p2, p1, rg))

    if not rows:
        raise RuntimeError("No valid rg records parsed from LDSC result files.")

    rg_df = pd.DataFrame(rows, columns=["trait1", "trait2", "rg"])
    traits = sorted(set(rg_df["trait1"]).union(set(rg_df["trait2"])))
    mat = pd.DataFrame(np.nan, index=traits, columns=traits)
    np.fill_diagonal(mat.values, 1.0)

    for row in rg_df.itertuples(index=False):
        mat.loc[row.trait1, row.trait2] = row.rg

    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(max(8, len(traits) * 0.5), max(6, len(traits) * 0.45)))
    sns.heatmap(
        mat,
        cmap="coolwarm",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.5,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Genetic correlation (rg)"},
    )
    plt.title("LDSC Genetic Correlation Heatmap")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close()
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
