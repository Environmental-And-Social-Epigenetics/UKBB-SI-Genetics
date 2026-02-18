#!/usr/bin/env python3
"""
Generate integrated sex-difference visualizations.
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
    out_root = fg_root / "7_Sex_Differences" / "output"
    parser = argparse.ArgumentParser(description="Generate sex-difference plots.")
    parser.add_argument("--cross-sex-rg", type=Path, default=out_root / "01_cross_sex_rg" / "cross_sex_rg_summary.tsv")
    parser.add_argument("--h2-comparison", type=Path, default=out_root / "02_h2_sex_comparison" / "h2_sex_comparison.tsv")
    parser.add_argument("--snp-interaction-dir", type=Path, default=out_root / "03_snp_sex_interaction")
    parser.add_argument("--locus-summary", type=Path, default=out_root / "04_locus_classification" / "locus_classification_summary.tsv")
    parser.add_argument("--diff-genes", type=Path, default=out_root / "06_sex_differential_genes" / "sex_differential_genes_all.tsv")
    parser.add_argument("--pathway-table", type=Path, default=out_root / "07_sex_differential_pathways" / "sex_differential_pathways.tsv")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=out_root / "08_sex_difference_plots",
    )
    parser.add_argument("--scatter-max-points", type=int, default=120000)
    parser.add_argument("--manhattan-max-points", type=int, default=450000)
    return parser.parse_args()


def add_cumulative_position(df: pd.DataFrame, chr_col: str = "CHR", bp_col: str = "BP") -> tuple[pd.DataFrame, list[float], list[str]]:
    out = df.copy()
    chroms = sorted(out[chr_col].astype(int).unique().tolist())
    chr_max = {c: int(out.loc[out[chr_col] == c, bp_col].max()) for c in chroms}
    offsets = {}
    tick_pos = []
    tick_labels = []
    current = 0
    for c in chroms:
        offsets[c] = current
        tick_pos.append(current + chr_max[c] / 2)
        tick_labels.append(str(c))
        current += chr_max[c]
    out["x"] = out.apply(lambda r: int(r[bp_col]) + offsets[int(r[chr_col])], axis=1)
    return out, tick_pos, tick_labels


def parse_interaction_name(path: Path) -> tuple[str, str]:
    # Example: SI__trait1__AbilityToConfide.sex_interaction.tsv.gz
    stem = path.name.replace(".sex_interaction.tsv.gz", "")
    parts = stem.split("__")
    coding = "continuous" if parts[0] == "SI_continuous" else "binary"
    trait = parts[2] if len(parts) > 2 else "trait"
    return coding, trait


def plot_scatter(df: pd.DataFrame, title: str, out_file: Path, max_points: int) -> None:
    if len(df) > max_points:
        df = df.sample(n=max_points, random_state=42)
    color = -np.log10(df["p_interaction"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(7, 7))
    sc = ax.scatter(df["beta_male"], df["beta_female"], c=color, s=4, alpha=0.5, cmap="viridis")
    lim = np.nanpercentile(np.abs(np.concatenate([df["beta_male"].values, df["beta_female"].values])), 99.5)
    lim = max(lim, 0.02)
    ax.plot([-lim, lim], [-lim, lim], linestyle="--", linewidth=1, color="red")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("beta male")
    ax.set_ylabel("beta female")
    ax.set_title(title)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("-log10(p interaction)")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close(fig)


def plot_manhattan(df: pd.DataFrame, title: str, out_file: Path, max_points: int) -> None:
    work = df.copy()
    work["CHR"] = work["CHR"].astype(int)
    work["BP"] = work["BP"].astype(int)
    work["neglog10p"] = -np.log10(work["p_interaction"].clip(lower=1e-300))
    if len(work) > max_points:
        # Keep the strongest points and sample the remainder for clarity/performance.
        top = work.nsmallest(100000, "p_interaction")
        rest = work.drop(top.index)
        n_remaining = max(0, max_points - len(top))
        if n_remaining > 0 and len(rest) > n_remaining:
            rest = rest.sample(n=n_remaining, random_state=42)
        work = pd.concat([top, rest], ignore_index=True)

    work, ticks, labels = add_cumulative_position(work)
    fig, ax = plt.subplots(figsize=(14, 5))
    colors = ["#4C78A8", "#F58518"]
    for idx, chrom in enumerate(sorted(work["CHR"].unique())):
        sub = work[work["CHR"] == chrom]
        ax.scatter(sub["x"], sub["neglog10p"], s=3, alpha=0.65, color=colors[idx % 2])
    ax.axhline(-np.log10(5e-8), color="red", linestyle="--", linewidth=1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(p interaction)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close(fig)


def plot_miami(df: pd.DataFrame, title: str, out_file: Path, max_points: int) -> None:
    work = df.copy()
    work["CHR"] = work["CHR"].astype(int)
    work["BP"] = work["BP"].astype(int)
    if len(work) > max_points:
        work = work.sample(n=max_points, random_state=42)
    work, ticks, labels = add_cumulative_position(work)
    work["male_y"] = -np.log10(work["p_male"].clip(lower=1e-300))
    work["female_y"] = -np.log10(work["p_female"].clip(lower=1e-300)) * -1
    fig, ax = plt.subplots(figsize=(14, 6))
    colors = ["#4C78A8", "#F58518"]
    for idx, chrom in enumerate(sorted(work["CHR"].unique())):
        sub = work[work["CHR"] == chrom]
        ax.scatter(sub["x"], sub["male_y"], s=3, alpha=0.65, color=colors[idx % 2])
        ax.scatter(sub["x"], sub["female_y"], s=3, alpha=0.65, color=colors[idx % 2])
    sig = -np.log10(5e-8)
    ax.axhline(sig, color="red", linestyle="--", linewidth=1)
    ax.axhline(-sig, color="red", linestyle="--", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(p) male (top) / female (bottom)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sns.set_style("whitegrid")

    # 1/2/3: per-combination scatter, interaction Manhattan, male-vs-female Miami.
    interaction_files = sorted(args.snp_interaction_dir.glob("*.sex_interaction.tsv.gz"))
    if not interaction_files:
        raise FileNotFoundError(f"No interaction files found in {args.snp_interaction_dir}")
    for f in interaction_files:
        coding, trait = parse_interaction_name(f)
        df = pd.read_csv(
            f,
            sep="\t",
            usecols=["CHR", "BP", "beta_male", "beta_female", "p_male", "p_female", "p_interaction"],
        )
        base = f"{coding}__{trait}"
        plot_scatter(
            df=df,
            title=f"Cross-sex effect scatter: {trait} ({coding})",
            out_file=args.out_dir / f"{base}.cross_sex_scatter.png",
            max_points=args.scatter_max_points,
        )
        plot_manhattan(
            df=df,
            title=f"Sex-interaction Manhattan: {trait} ({coding})",
            out_file=args.out_dir / f"{base}.sex_interaction_manhattan.png",
            max_points=args.manhattan_max_points,
        )
        plot_miami(
            df=df,
            title=f"Male vs Female Miami: {trait} ({coding})",
            out_file=args.out_dir / f"{base}.male_female_miami.png",
            max_points=args.manhattan_max_points,
        )

    # 4: h2 comparison bar chart.
    h2 = pd.read_csv(args.h2_comparison, sep="\t")
    h2["group"] = h2["coding"] + "__" + h2["trait"]
    long_h2 = pd.concat(
        [
            h2[["group", "coding", "trait", "h2_male", "se_male", "p_diff"]].rename(
                columns={"h2_male": "h2", "se_male": "se"}
            ).assign(sex="male"),
            h2[["group", "coding", "trait", "h2_female", "se_female", "p_diff"]].rename(
                columns={"h2_female": "h2", "se_female": "se"}
            ).assign(sex="female"),
        ],
        ignore_index=True,
    )
    order = h2.sort_values(["coding", "trait_index"])["group"].tolist()
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.barplot(data=long_h2, x="group", y="h2", hue="sex", order=order, ax=ax)
    offsets = {"male": -0.2, "female": 0.2}
    x_map = {g: i for i, g in enumerate(order)}
    for _, r in long_h2.iterrows():
        ax.errorbar(x=x_map[r["group"]] + offsets[r["sex"]], y=r["h2"], yerr=r["se"], fmt="none", ecolor="black", capsize=3)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=25, ha="right")
    ax.set_xlabel("Trait (coding)")
    ax.set_ylabel("h2")
    ax.set_title("Heritability comparison: male vs female")
    plt.tight_layout()
    plt.savefig(args.out_dir / "h2_sex_comparison_bar.png", dpi=300)
    plt.close(fig)

    # 5: cross-sex rg forest plot.
    rg = pd.read_csv(args.cross_sex_rg, sep="\t")
    rg["coding"] = np.where(rg["coding"] == "SI_continuous", "continuous", "binary")
    rg["label"] = rg["coding"] + "__" + rg["trait"]
    rg = rg.sort_values(["coding", "trait_index"])
    fig, ax = plt.subplots(figsize=(9, 5))
    y = np.arange(len(rg))
    ax.errorbar(rg["rg"], y, xerr=1.96 * rg["se"], fmt="o", color="#4C78A8", ecolor="gray", capsize=3)
    ax.axvline(1.0, color="red", linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(rg["label"])
    ax.set_xlabel("Cross-sex rg")
    ax.set_title("Cross-sex LDSC rg (male vs female)")
    plt.tight_layout()
    plt.savefig(args.out_dir / "cross_sex_rg_forest.png", dpi=300)
    plt.close(fig)

    # 6: locus classification stacked bar chart.
    locus = pd.read_csv(args.locus_summary, sep="\t")
    count_cols = [c for c in locus.columns if c.startswith("n_") and c not in {"n_loci_total"}]
    locus_long = locus.melt(
        id_vars=["coding", "trait"],
        value_vars=count_cols,
        var_name="category",
        value_name="count",
    )
    locus_long["group"] = locus_long["coding"] + "__" + locus_long["trait"]
    pivot = locus_long.pivot_table(index="group", columns="category", values="count", aggfunc="sum").fillna(0)
    pivot = pivot.loc[sorted(pivot.index)]
    fig, ax = plt.subplots(figsize=(12, 6))
    bottom = np.zeros(len(pivot))
    for col in pivot.columns:
        ax.bar(pivot.index, pivot[col].values, bottom=bottom, label=col)
        bottom += pivot[col].values
    ax.set_xticklabels(pivot.index, rotation=25, ha="right")
    ax.set_ylabel("Locus count")
    ax.set_title("Locus-level sex classification")
    ax.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(args.out_dir / "locus_classification_stacked.png", dpi=300)
    plt.close(fig)

    # 7: sex-differential gene volcano plot (global).
    genes = pd.read_csv(args.diff_genes, sep="\t")
    genes["neglog10_pdiff"] = -np.log10(genes["p_diff"].clip(lower=1e-300))
    genes["coding_trait"] = genes["coding"] + "__" + genes["trait"]
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.scatterplot(
        data=genes.sample(n=min(150000, len(genes)), random_state=42),
        x="z_diff_male_minus_female",
        y="neglog10_pdiff",
        hue="coding",
        alpha=0.35,
        s=10,
        linewidth=0,
        ax=ax,
    )
    ax.axvline(0, color="black", linewidth=1)
    ax.axhline(-np.log10(1e-5), color="red", linestyle="--", linewidth=1)
    ax.set_xlabel("Z diff (male - female)")
    ax.set_ylabel("-log10(p diff)")
    ax.set_title("Sex-differential genes volcano")
    plt.tight_layout()
    plt.savefig(args.out_dir / "sex_differential_genes_volcano.png", dpi=300)
    plt.close(fig)

    # 8: pathway sex-difference heatmap.
    path = pd.read_csv(args.pathway_table, sep="\t")
    path["column"] = path["coding"] + "__" + path["trait"]
    path["signed"] = -np.log10(path["q_enrichment_bh"].fillna(1.0).clip(lower=1e-300))
    path["signed"] *= np.where(path["stronger_sex"] == "male", 1.0, -1.0)
    mat = path.pivot_table(index="gene_set", columns="column", values="signed", aggfunc="first").fillna(0.0)
    fig, ax = plt.subplots(figsize=(12, max(4, 0.5 * len(mat))))
    sns.heatmap(mat, cmap="coolwarm", center=0, ax=ax, cbar_kws={"label": "Signed -log10(q): male(+) female(-)"})
    ax.set_title("Pathway sex-difference heatmap")
    plt.tight_layout()
    plt.savefig(args.out_dir / "pathway_sex_difference_heatmap.png", dpi=300)
    plt.close(fig)

    print(f"Wrote plots to: {args.out_dir}")


if __name__ == "__main__":
    main()
