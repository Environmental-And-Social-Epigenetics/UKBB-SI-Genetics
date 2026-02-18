#!/usr/bin/env python3
"""
Formal comparison of male vs female SNP heritability using LDSC logs.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import norm


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description="Compare male vs female h2 from LDSC output logs.")
    parser.add_argument(
        "--h2-dir",
        type=Path,
        default=fg_root / "1_LDSC" / "output" / "h2",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=fg_root / "7_Sex_Differences" / "output" / "02_h2_sex_comparison",
    )
    return parser.parse_args()


def parse_h2_from_log(log_path: Path) -> tuple[float, float] | None:
    text = log_path.read_text()
    match = re.search(r"Total Observed scale h2:\s+([-\d.eE]+)\s+\(([-\d.eE]+)\)", text)
    if not match:
        return None
    return float(match.group(1)), float(match.group(2))


def parse_meta_from_name(stem: str) -> tuple[str, str, str, str]:
    # Example: SI_continuous__EUR_Female_MM__trait1__AbilityToConfide
    parts = stem.split("__")
    coding = parts[0]
    population = parts[1]
    trait_index = parts[2].replace("trait", "")
    trait = parts[3]
    coding_label = "continuous" if coding == "SI_continuous" else "binary"
    sex = "female" if "Female" in population else "male" if "Male" in population else "combined"
    return coding_label, sex, trait_index, trait


def add_sig_label(p: float) -> str:
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for log_path in sorted(args.h2_dir.glob("*.log")):
        parsed = parse_h2_from_log(log_path)
        if parsed is None:
            continue
        h2, se = parsed
        coding, sex, trait_index, trait = parse_meta_from_name(log_path.stem)
        if sex not in {"male", "female"}:
            continue
        rows.append(
            {
                "coding": coding,
                "sex": sex,
                "trait_index": int(trait_index),
                "trait": trait,
                "h2": h2,
                "se": se,
                "log_file": str(log_path),
            }
        )

    if not rows:
        raise RuntimeError(f"No male/female h2 logs parsed from {args.h2_dir}")

    df = pd.DataFrame(rows)
    pivot = (
        df.pivot_table(
            index=["coding", "trait_index", "trait"],
            columns="sex",
            values=["h2", "se"],
            aggfunc="first",
        )
        .reset_index()
    )
    # Flatten columns.
    pivot.columns = [
        c if isinstance(c, str) else f"{c[0]}_{c[1]}".rstrip("_")
        for c in pivot.columns.to_flat_index()
    ]

    pivot["h2_diff_female_minus_male"] = pivot["h2_female"] - pivot["h2_male"]
    pivot["se_diff"] = (pivot["se_female"] ** 2 + pivot["se_male"] ** 2).pow(0.5)
    pivot["z_diff"] = pivot["h2_diff_female_minus_male"] / pivot["se_diff"]
    pivot["p_diff"] = 2.0 * norm.sf(pivot["z_diff"].abs())
    pivot["sig_label"] = pivot["p_diff"].map(add_sig_label)

    out_tsv = args.out_dir / "h2_sex_comparison.tsv"
    pivot.sort_values(["coding", "trait_index"]).to_csv(out_tsv, sep="\t", index=False)

    long_df = []
    for _, row in pivot.iterrows():
        for sex in ["male", "female"]:
            long_df.append(
                {
                    "coding": row["coding"],
                    "trait_index": row["trait_index"],
                    "trait": row["trait"],
                    "sex": sex,
                    "h2": row[f"h2_{sex}"],
                    "se": row[f"se_{sex}"],
                    "p_diff": row["p_diff"],
                    "sig_label": row["sig_label"],
                }
            )
    long_df = pd.DataFrame(long_df)
    long_df["group"] = long_df["coding"] + "__" + long_df["trait"]

    sns.set_style("whitegrid")

    # Plot 1: side-by-side h2 bars with error bars.
    fig, ax = plt.subplots(figsize=(12, 6))
    x_order = (
        long_df[["group", "coding", "trait_index", "trait"]]
        .drop_duplicates()
        .sort_values(["coding", "trait_index"])
    )
    order = x_order["group"].tolist()
    sns.barplot(
        data=long_df,
        x="group",
        y="h2",
        hue="sex",
        order=order,
        ax=ax,
    )
    # Add error bars manually.
    x_positions = {g: i for i, g in enumerate(order)}
    offsets = {"male": -0.2, "female": 0.2}
    for _, row in long_df.iterrows():
        ax.errorbar(
            x=x_positions[row["group"]] + offsets[row["sex"]],
            y=row["h2"],
            yerr=row["se"],
            fmt="none",
            ecolor="black",
            elinewidth=1,
            capsize=3,
        )

    # Add significance labels per group.
    top_by_group = long_df.groupby("group")["h2"].max().to_dict()
    se_by_group = long_df.groupby("group")["se"].max().to_dict()
    pmap = pivot.set_index(pivot["coding"] + "__" + pivot["trait"])["sig_label"].to_dict()
    for g in order:
        y = top_by_group[g] + se_by_group[g] + 0.002
        ax.text(x_positions[g], y, pmap.get(g, "ns"), ha="center", va="bottom", fontsize=10)

    ax.set_xlabel("Trait (coding)")
    ax.set_ylabel("Observed-scale SNP heritability (h2)")
    ax.set_title("Male vs Female h2 by Trait and Coding")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=25, ha="right")
    ax.legend(title="Sex")
    plt.tight_layout()
    h2_plot = args.out_dir / "h2_male_vs_female.png"
    plt.savefig(h2_plot, dpi=300)
    plt.close(fig)

    # Plot 2: female-minus-male z-differences.
    zdf = pivot.sort_values(["coding", "trait_index"]).copy()
    zdf["group"] = zdf["coding"] + "__" + zdf["trait"]
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(data=zdf, x="group", y="z_diff", color="#4C78A8", ax=ax)
    ax.axhline(0, color="black", linewidth=1)
    ax.axhline(1.96, color="red", linestyle="--", linewidth=1)
    ax.axhline(-1.96, color="red", linestyle="--", linewidth=1)
    for i, row in zdf.reset_index(drop=True).iterrows():
        ax.text(i, row["z_diff"] + (0.08 if row["z_diff"] >= 0 else -0.08), row["sig_label"], ha="center")
    ax.set_xlabel("Trait (coding)")
    ax.set_ylabel("Z difference (female - male)")
    ax.set_title("Sex Difference in h2 (female - male)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=25, ha="right")
    plt.tight_layout()
    z_plot = args.out_dir / "h2_sex_difference_zscores.png"
    plt.savefig(z_plot, dpi=300)
    plt.close(fig)

    print(f"Wrote: {out_tsv}")
    print(f"Wrote: {h2_plot}")
    print(f"Wrote: {z_plot}")


if __name__ == "__main__":
    main()
