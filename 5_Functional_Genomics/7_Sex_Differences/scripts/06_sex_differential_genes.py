#!/usr/bin/env python3
"""
Compare male vs female MAGMA gene associations and identify sex-differential genes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

TRAITS = {
    1: "AbilityToConfide",
    2: "FreqSoc",
    3: "Loneliness",
}


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description="Sex-differential MAGMA gene comparison.")
    parser.add_argument(
        "--magma-gene-dir",
        type=Path,
        default=fg_root / "7_Sex_Differences" / "output" / "05_magma_sex_stratified" / "gene",
    )
    parser.add_argument(
        "--gene-info",
        type=Path,
        default=fg_root / "reference_data" / "magma" / "gene_info_human.tsv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=fg_root / "7_Sex_Differences" / "output" / "06_sex_differential_genes",
    )
    parser.add_argument("--magma-alpha", type=float, default=2.7e-6)
    parser.add_argument("--p-null", type=float, default=0.05)
    parser.add_argument("--p-diff-threshold", type=float, default=1e-5)
    return parser.parse_args()


def bh_fdr(p: pd.Series) -> pd.Series:
    p = p.fillna(1.0).astype(float)
    n = len(p)
    order = np.argsort(p.values)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    q = p.values * n / ranks
    # Enforce monotonicity.
    q_sorted = q[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_final = np.empty(n, dtype=float)
    q_final[order] = np.clip(q_sorted, 0, 1)
    return pd.Series(q_final, index=p.index)


def read_magma_gene_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", low_memory=False)
    required = {"GENE", "ZSTAT", "P"}
    if not required.issubset(df.columns):
        raise ValueError(f"Missing required columns in {path}: {required - set(df.columns)}")
    return df[["GENE", "ZSTAT", "P"]].rename(columns={"GENE": "gene_id", "ZSTAT": "z", "P": "p"})


def load_gene_symbol_map(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return {}
    if not {"gene_id", "symbol"}.issubset(df.columns):
        return {}
    return {str(g): str(s) for g, s in zip(df["gene_id"], df["symbol"])}


def classify_row(
    p_male: float,
    p_female: float,
    z_diff: float,
    p_diff: float,
    magma_alpha: float,
    p_null: float,
    p_diff_threshold: float,
) -> str:
    male_sig = p_male < magma_alpha
    female_sig = p_female < magma_alpha
    if male_sig and female_sig:
        return "shared_significant"
    if female_sig and p_male > p_null:
        return "female_specific"
    if male_sig and p_female > p_null:
        return "male_specific"
    if female_sig and not male_sig and z_diff < 0:
        return "female_enriched"
    if male_sig and not female_sig and z_diff > 0:
        return "male_enriched"
    if p_diff < p_diff_threshold:
        return "sex_differential"
    return "not_significant"


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    symbol_map = load_gene_symbol_map(args.gene_info)
    all_rows: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []

    for coding in ["SI", "SI_continuous"]:
        coding_label = "binary" if coding == "SI" else "continuous"
        for trait_idx, trait_name in TRAITS.items():
            male_path = args.magma_gene_dir / f"{coding}__EUR_Male_MM__trait{trait_idx}__{trait_name}.genes.out"
            female_path = args.magma_gene_dir / f"{coding}__EUR_Female_MM__trait{trait_idx}__{trait_name}.genes.out"
            if not male_path.exists() or not female_path.exists():
                print(f"[skip] Missing MAGMA pair: {male_path.name} / {female_path.name}")
                continue

            male = read_magma_gene_file(male_path).rename(columns={"z": "z_male", "p": "p_male"})
            female = read_magma_gene_file(female_path).rename(columns={"z": "z_female", "p": "p_female"})
            merged = male.merge(female, on="gene_id", how="inner")
            merged["coding"] = coding_label
            merged["coding_raw"] = coding
            merged["trait_index"] = trait_idx
            merged["trait"] = trait_name
            merged["z_diff_male_minus_female"] = (merged["z_male"] - merged["z_female"]) / np.sqrt(2.0)
            merged["p_diff"] = 2.0 * norm.sf(np.abs(merged["z_diff_male_minus_female"]))
            merged["q_diff_bh"] = bh_fdr(merged["p_diff"])
            merged["male_sig"] = merged["p_male"] < args.magma_alpha
            merged["female_sig"] = merged["p_female"] < args.magma_alpha
            merged["stronger_sex"] = np.where(
                merged["z_diff_male_minus_female"] > 0, "male", "female"
            )
            merged["gene_symbol"] = merged["gene_id"].astype(str).map(symbol_map).fillna(merged["gene_id"].astype(str))
            merged["category"] = merged.apply(
                lambda r: classify_row(
                    p_male=float(r["p_male"]),
                    p_female=float(r["p_female"]),
                    z_diff=float(r["z_diff_male_minus_female"]),
                    p_diff=float(r["p_diff"]),
                    magma_alpha=args.magma_alpha,
                    p_null=args.p_null,
                    p_diff_threshold=args.p_diff_threshold,
                ),
                axis=1,
            )
            merged = merged.sort_values("p_diff", ascending=True)

            per_combo = args.out_dir / f"{coding}__trait{trait_idx}__{trait_name}.sex_differential_genes.tsv"
            merged.to_csv(per_combo, sep="\t", index=False)
            all_rows.append(merged)

            category_counts = merged["category"].value_counts().to_dict()
            summary_row = {
                "coding": coding_label,
                "coding_raw": coding,
                "trait_index": trait_idx,
                "trait": trait_name,
                "n_genes_tested": int(len(merged)),
                "n_male_sig": int(merged["male_sig"].sum()),
                "n_female_sig": int(merged["female_sig"].sum()),
                "n_sex_differential_pdiff": int((merged["p_diff"] < args.p_diff_threshold).sum()),
                "n_sex_differential_qdiff": int((merged["q_diff_bh"] < 0.05).sum()),
            }
            for k in [
                "shared_significant",
                "female_specific",
                "male_specific",
                "female_enriched",
                "male_enriched",
                "sex_differential",
                "not_significant",
            ]:
                summary_row[f"n_{k}"] = int(category_counts.get(k, 0))
            summary_rows.append(summary_row)

    if not all_rows:
        raise RuntimeError(f"No male/female MAGMA gene result pairs found in {args.magma_gene_dir}")

    all_df = pd.concat(all_rows, ignore_index=True)
    all_df.to_csv(args.out_dir / "sex_differential_genes_all.tsv", sep="\t", index=False)
    summary_df = pd.DataFrame(summary_rows).sort_values(["coding", "trait_index"])
    summary_df.to_csv(args.out_dir / "sex_differential_genes_summary.tsv", sep="\t", index=False)

    top = all_df.sort_values("p_diff", ascending=True).head(500)
    top.to_csv(args.out_dir / "sex_differential_genes_top500.tsv", sep="\t", index=False)

    print(f"Wrote: {args.out_dir / 'sex_differential_genes_all.tsv'}")
    print(f"Wrote: {args.out_dir / 'sex_differential_genes_summary.tsv'}")
    print(f"Wrote: {args.out_dir / 'sex_differential_genes_top500.tsv'}")


if __name__ == "__main__":
    main()
