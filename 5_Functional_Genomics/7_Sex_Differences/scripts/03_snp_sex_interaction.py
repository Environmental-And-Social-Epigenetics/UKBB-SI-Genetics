#!/usr/bin/env python3
"""
Compute SNP-level sex interaction statistics from male vs female MTAG outputs.
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
    repo_root = fg_root.parent
    parser = argparse.ArgumentParser(description="SNP-level male vs female sex interaction from MTAG.")
    parser.add_argument("--repo-root", type=Path, default=repo_root)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=fg_root / "7_Sex_Differences" / "output" / "03_snp_sex_interaction",
    )
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--top-n", type=int, default=200)
    return parser.parse_args()


def read_mtag(path: Path) -> pd.DataFrame:
    usecols = [
        "SNP",
        "CHR",
        "BP",
        "A1",
        "A2",
        "N",
        "FRQ",
        "mtag_beta",
        "mtag_se",
        "mtag_z",
        "mtag_pval",
    ]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    df = df.rename(
        columns={
            "N": "n",
            "FRQ": "frq",
            "mtag_beta": "beta",
            "mtag_se": "se",
            "mtag_z": "z",
            "mtag_pval": "p",
        }
    )
    df["SNP"] = df["SNP"].astype(str)
    for col in ["CHR", "BP", "n", "frq", "beta", "se", "z", "p"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["SNP", "CHR", "BP", "A1", "A2", "beta", "se", "p"])
    df = df[(df["se"] > 0) & (df["p"] > 0) & (df["p"] <= 1)]
    df = df.sort_values("p").drop_duplicates("SNP", keep="first").reset_index(drop=True)
    return df


def harmonize_male_female(male: pd.DataFrame, female: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    merged = male.merge(
        female,
        on="SNP",
        suffixes=("_male", "_female"),
        how="inner",
    )
    merged = merged[(merged["CHR_male"] == merged["CHR_female"]) & (merged["BP_male"] == merged["BP_female"])].copy()
    merged = merged.rename(columns={"CHR_male": "CHR", "BP_male": "BP"})

    same = (merged["A1_male"] == merged["A1_female"]) & (merged["A2_male"] == merged["A2_female"])
    flipped = (merged["A1_male"] == merged["A2_female"]) & (merged["A2_male"] == merged["A1_female"])

    # Align flipped female alleles to male orientation.
    merged.loc[flipped, "beta_female"] = -merged.loc[flipped, "beta_female"]
    merged.loc[flipped, "z_female"] = -merged.loc[flipped, "z_female"]
    merged.loc[flipped, "frq_female"] = 1.0 - merged.loc[flipped, "frq_female"]
    merged.loc[flipped, "A1_female"] = merged.loc[flipped, "A1_male"]
    merged.loc[flipped, "A2_female"] = merged.loc[flipped, "A2_male"]

    aligned = merged[same | flipped].copy()
    aligned["allele_alignment"] = np.where(same[same | flipped], "same", "flipped")

    stats = {
        "n_overlap_snps": int(len(merged)),
        "n_aligned_same": int(same.sum()),
        "n_aligned_flipped": int(flipped.sum()),
        "n_discarded_allele_mismatch": int((~(same | flipped)).sum()),
        "n_aligned_total": int(len(aligned)),
    }
    return aligned, stats


def compute_interaction(df: pd.DataFrame) -> pd.DataFrame:
    se_int = np.sqrt(df["se_male"] ** 2 + df["se_female"] ** 2)
    z_int = (df["beta_male"] - df["beta_female"]) / se_int
    p_int = 2.0 * norm.sf(np.abs(z_int))

    out = pd.DataFrame(
        {
            "SNP": df["SNP"],
            "CHR": df["CHR"].astype(int),
            "BP": df["BP"].astype(int),
            "A1": df["A1_male"],
            "A2": df["A2_male"],
            "allele_alignment": df["allele_alignment"],
            "beta_male": df["beta_male"],
            "se_male": df["se_male"],
            "p_male": df["p_male"],
            "beta_female": df["beta_female"],
            "se_female": df["se_female"],
            "p_female": df["p_female"],
            "beta_diff_male_minus_female": df["beta_male"] - df["beta_female"],
            "se_interaction": se_int,
            "z_interaction": z_int,
            "p_interaction": p_int,
            "same_direction": (np.sign(df["beta_male"]) == np.sign(df["beta_female"])),
            "n_male": df["n_male"],
            "n_female": df["n_female"],
        }
    )
    out = out.sort_values("p_interaction", ascending=True).reset_index(drop=True)
    return out


def get_mtag_paths(repo_root: Path, coding: str, trait_index: int) -> tuple[Path, Path]:
    if coding == "SI":
        male = repo_root / "3_MTAG" / "SI" / "EUR_Male_MM" / "results" / f"SI_EUR_Male_MM_Output_trait_{trait_index}.txt"
        female = repo_root / "3_MTAG" / "SI" / "EUR_Female_MM" / "results" / f"SI_EUR_Female_MM_Output_trait_{trait_index}.txt"
    elif coding == "SI_continuous":
        male = (
            repo_root
            / "3_MTAG"
            / "SI_continuous"
            / "EUR_Male_MM"
            / "results"
            / f"SI_EUR_Male_MM_Output_continuous_trait_{trait_index}.txt"
        )
        female = (
            repo_root
            / "3_MTAG"
            / "SI_continuous"
            / "EUR_Female_MM"
            / "results"
            / f"SI_EUR_Female_MM_Output_continuous_trait_{trait_index}.txt"
        )
    else:
        raise ValueError(f"Unsupported coding: {coding}")
    return male, female


def main() -> None:
    args = parse_args()
    args.repo_root = args.repo_root.resolve()
    args.out_dir = args.out_dir.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []

    for coding in ["SI", "SI_continuous"]:
        coding_label = "binary" if coding == "SI" else "continuous"
        for trait_index, trait_name in TRAITS.items():
            male_path, female_path = get_mtag_paths(args.repo_root, coding, trait_index)
            if not male_path.exists() or not female_path.exists():
                raise FileNotFoundError(f"Missing files for {coding} trait{trait_index}: {male_path} / {female_path}")

            print(f"[sex-interaction] {coding_label} trait{trait_index} {trait_name}")
            male = read_mtag(male_path)
            female = read_mtag(female_path)
            aligned, align_stats = harmonize_male_female(male, female)
            interaction = compute_interaction(aligned)

            prefix = f"{coding}__trait{trait_index}__{trait_name}"
            out_all = args.out_dir / f"{prefix}.sex_interaction.tsv.gz"
            out_top = args.out_dir / f"{prefix}.top_hits.tsv"
            interaction.to_csv(out_all, sep="\t", index=False, compression="gzip")
            interaction.head(args.top_n).to_csv(out_top, sep="\t", index=False)

            n_sig = int((interaction["p_interaction"] < args.p_threshold).sum())
            n_concordant = int(interaction["same_direction"].sum())
            n_discordant = int((~interaction["same_direction"]).sum())

            summary_rows.append(
                {
                    "coding": coding_label,
                    "coding_raw": coding,
                    "trait_index": trait_index,
                    "trait": trait_name,
                    "n_overlap_snps": align_stats["n_overlap_snps"],
                    "n_aligned_total": align_stats["n_aligned_total"],
                    "n_aligned_same": align_stats["n_aligned_same"],
                    "n_aligned_flipped": align_stats["n_aligned_flipped"],
                    "n_discarded_allele_mismatch": align_stats["n_discarded_allele_mismatch"],
                    "n_direction_concordant": n_concordant,
                    "n_direction_discordant": n_discordant,
                    "frac_direction_concordant": n_concordant / max(1, len(interaction)),
                    "n_sig_interaction_p5e8": n_sig,
                    "min_p_interaction": float(interaction["p_interaction"].min()),
                    "median_abs_beta_diff": float((interaction["beta_diff_male_minus_female"].abs()).median()),
                    "output_file": str(out_all),
                    "top_hits_file": str(out_top),
                }
            )

    summary = pd.DataFrame(summary_rows).sort_values(["coding", "trait_index"])
    summary_tsv = args.out_dir / "sex_interaction_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)
    print(f"Wrote: {summary_tsv}")


if __name__ == "__main__":
    main()
