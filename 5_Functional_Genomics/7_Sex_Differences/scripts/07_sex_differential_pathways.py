#!/usr/bin/env python3
"""
Rank-based pathway enrichment on gene-level sex-differential Z scores.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description="Pathway enrichment from sex-differential gene scores.")
    parser.add_argument(
        "--sex-diff-genes",
        type=Path,
        default=fg_root / "7_Sex_Differences" / "output" / "06_sex_differential_genes" / "sex_differential_genes_all.tsv",
    )
    parser.add_argument(
        "--gene-info",
        type=Path,
        default=fg_root / "reference_data" / "magma" / "gene_info_human.tsv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=fg_root / "7_Sex_Differences" / "output" / "07_sex_differential_pathways",
    )
    parser.add_argument("--min-genes-in-set", type=int, default=5)
    return parser.parse_args()


def bh_fdr(p: pd.Series) -> pd.Series:
    p = p.fillna(1.0).astype(float)
    n = len(p)
    order = np.argsort(p.values)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    q = p.values * n / ranks
    q_sorted = q[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_final = np.empty(n, dtype=float)
    q_final[order] = np.clip(q_sorted, 0, 1)
    return pd.Series(q_final, index=p.index)


def curated_gene_sets(x_linked_symbols: set[str]) -> dict[str, set[str]]:
    return {
        "hormone_related": {
            "ESR1",
            "ESR2",
            "AR",
            "PGR",
            "CYP19A1",
            "HSD17B1",
            "NR3C1",
            "NR3C2",
            "CRH",
            "CRHR1",
            "OXTR",
            "AVPR1A",
        },
        "hpa_axis_stress": {
            "CRH",
            "CRHR1",
            "CRHR2",
            "NR3C1",
            "FKBP5",
            "MC2R",
            "POMC",
            "AVP",
            "ADCYAP1",
            "RORA",
        },
        "serotonin_signaling": {
            "SLC6A4",
            "TPH1",
            "TPH2",
            "HTR1A",
            "HTR1B",
            "HTR2A",
            "HTR2C",
            "MAOA",
            "MAOB",
            "DDC",
        },
        "dopamine_signaling": {
            "TH",
            "DDC",
            "SLC6A3",
            "DRD1",
            "DRD2",
            "DRD3",
            "DRD4",
            "DRD5",
            "COMT",
            "MAOA",
        },
        "oxytocin_vasopressin": {
            "OXTR",
            "OXT",
            "AVP",
            "AVPR1A",
            "AVPR1B",
            "AVPR2",
            "CD38",
            "LNPEP",
        },
        "gaba_glutamate": {
            "GAD1",
            "GAD2",
            "GABRA1",
            "GABRB2",
            "GABRG2",
            "SLC32A1",
            "GRIN1",
            "GRIN2A",
            "GRM5",
            "SLC17A7",
        },
        "microglia_immune": {
            "CX3CR1",
            "TREM2",
            "TYROBP",
            "AIF1",
            "CSF1R",
            "C1QA",
            "C1QB",
            "C1QC",
            "IL6",
            "TNF",
        },
        "x_linked_genes": set(x_linked_symbols),
    }


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(args.sex_diff_genes, sep="\t")
    if "gene_symbol" not in genes.columns or "z_diff_male_minus_female" not in genes.columns:
        raise ValueError("Input sex differential genes file missing required columns.")
    genes["gene_symbol"] = genes["gene_symbol"].astype(str).str.upper()

    x_linked = set()
    if args.gene_info.exists():
        info = pd.read_csv(args.gene_info, sep="\t")
        if {"chr", "symbol"}.issubset(info.columns):
            x_linked = set(
                info.loc[info["chr"].astype(str).isin(["23", "X"]), "symbol"]
                .astype(str)
                .str.upper()
                .tolist()
            )

    gene_sets = curated_gene_sets(x_linked_symbols=x_linked)
    rows = []

    for (coding, trait), sub in genes.groupby(["coding", "trait"]):
        z = sub["z_diff_male_minus_female"].astype(float)
        abs_z = z.abs()
        symbol_series = sub["gene_symbol"]
        universe_n = len(sub)
        if universe_n < 100:
            continue

        for set_name, gene_set in gene_sets.items():
            in_set_mask = symbol_series.isin(gene_set)
            n_in = int(in_set_mask.sum())
            n_out = int((~in_set_mask).sum())
            if n_in < args.min_genes_in_set or n_out < args.min_genes_in_set:
                rows.append(
                    {
                        "coding": coding,
                        "trait": trait,
                        "gene_set": set_name,
                        "n_in_set": n_in,
                        "n_out_set": n_out,
                        "mean_z_in_set": np.nan,
                        "median_abs_z_in_set": np.nan,
                        "median_abs_z_out_set": np.nan,
                        "u_stat": np.nan,
                        "p_enrichment": np.nan,
                    }
                )
                continue

            in_abs = abs_z[in_set_mask]
            out_abs = abs_z[~in_set_mask]
            stat = mannwhitneyu(in_abs, out_abs, alternative="greater")
            rows.append(
                {
                    "coding": coding,
                    "trait": trait,
                    "gene_set": set_name,
                    "n_in_set": n_in,
                    "n_out_set": n_out,
                    "mean_z_in_set": float(z[in_set_mask].mean()),
                    "median_abs_z_in_set": float(in_abs.median()),
                    "median_abs_z_out_set": float(out_abs.median()),
                    "u_stat": float(stat.statistic),
                    "p_enrichment": float(stat.pvalue),
                }
            )

    if not rows:
        raise RuntimeError("No pathway tests were performed. Check input files and min genes threshold.")

    out = pd.DataFrame(rows)
    out["q_enrichment_bh"] = np.nan
    for (coding, trait), idx in out.groupby(["coding", "trait"]).groups.items():
        pvals = out.loc[idx, "p_enrichment"]
        out.loc[idx, "q_enrichment_bh"] = bh_fdr(pvals)
    out["stronger_sex"] = np.where(out["mean_z_in_set"] > 0, "male", "female")

    out_tsv = args.out_dir / "sex_differential_pathways.tsv"
    out.sort_values(["coding", "trait", "q_enrichment_bh", "p_enrichment"]).to_csv(out_tsv, sep="\t", index=False)

    # Build heatmap-ready matrix: signed -log10(q)
    heat = out.copy()
    eps = 1e-300
    heat["signed_log10_q"] = -np.log10(heat["q_enrichment_bh"].fillna(1.0).clip(lower=eps))
    heat["signed_log10_q"] *= np.where(heat["stronger_sex"] == "male", 1.0, -1.0)
    heat["column"] = heat["coding"] + "__" + heat["trait"]
    heat_wide = heat.pivot_table(index="gene_set", columns="column", values="signed_log10_q", aggfunc="first")
    heat_wide.to_csv(args.out_dir / "sex_differential_pathways_heatmap_matrix.tsv", sep="\t")

    print(f"Wrote: {out_tsv}")
    print(f"Wrote: {args.out_dir / 'sex_differential_pathways_heatmap_matrix.tsv'}")


if __name__ == "__main__":
    main()
