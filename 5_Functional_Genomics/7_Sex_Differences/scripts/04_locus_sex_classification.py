#!/usr/bin/env python3
"""
Define male/female/combined loci and classify each locus by sex-specificity.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TRAITS = {
    1: "AbilityToConfide",
    2: "FreqSoc",
    3: "Loneliness",
}


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    repo_root = fg_root.parent
    parser = argparse.ArgumentParser(description="Classify sex-specific loci from MTAG outputs.")
    parser.add_argument("--repo-root", type=Path, default=repo_root)
    parser.add_argument(
        "--define-loci-script",
        type=Path,
        default=fg_root / "4_Fine_Mapping" / "scripts" / "01_define_loci.py",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=fg_root / "7_Sex_Differences" / "output" / "04_locus_classification",
    )
    parser.add_argument("--window-kb", type=int, default=500)
    parser.add_argument("--merge-distance-kb", type=int, default=250)
    parser.add_argument("--p-gws", type=float, default=5e-8)
    parser.add_argument("--p-suggestive", type=float, default=1e-5)
    parser.add_argument("--p-null", type=float, default=0.05)
    parser.add_argument("--force-rerun-define", action="store_true")
    return parser.parse_args()


def mtag_paths(repo_root: Path, coding: str, trait_index: int) -> dict[str, Path]:
    if coding == "SI":
        return {
            "male": repo_root / "3_MTAG" / "SI" / "EUR_Male_MM" / "results" / f"SI_EUR_Male_MM_Output_trait_{trait_index}.txt",
            "female": repo_root / "3_MTAG" / "SI" / "EUR_Female_MM" / "results" / f"SI_EUR_Female_MM_Output_trait_{trait_index}.txt",
            "combined": repo_root / "3_MTAG" / "SI" / "EUR_MM" / "results" / f"SI_EUR_MM_Output_trait_{trait_index}.txt",
        }
    if coding == "SI_continuous":
        return {
            "male": repo_root
            / "3_MTAG"
            / "SI_continuous"
            / "EUR_Male_MM"
            / "results"
            / f"SI_EUR_Male_MM_Output_continuous_trait_{trait_index}.txt",
            "female": repo_root
            / "3_MTAG"
            / "SI_continuous"
            / "EUR_Female_MM"
            / "results"
            / f"SI_EUR_Female_MM_Output_continuous_trait_{trait_index}.txt",
            "combined": repo_root
            / "3_MTAG"
            / "SI_continuous"
            / "EUR_MM"
            / "results"
            / f"SI_EUR_MM_Output_continuous_trait_{trait_index}.txt",
        }
    raise ValueError(coding)


def run_define_loci(
    define_script: Path,
    sumstats: Path,
    out_prefix: Path,
    force_rerun: bool,
) -> Path:
    loci_path = out_prefix.with_suffix(".loci.tsv")
    if loci_path.exists() and not force_rerun:
        return loci_path

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(define_script),
        "--sumstats",
        str(sumstats),
        "--out-prefix",
        str(out_prefix),
        "--allow-distance-fallback",
    ]
    subprocess.run(cmd, check=True)
    return loci_path


def read_sumstats(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["SNP", "CHR", "BP", "A1", "A2", "mtag_beta", "mtag_se", "mtag_pval"],
        low_memory=False,
    )
    df = df.rename(columns={"mtag_beta": "beta", "mtag_se": "se", "mtag_pval": "p"})
    df["SNP"] = df["SNP"].astype(str)
    for c in ["CHR", "BP", "beta", "se", "p"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["SNP", "CHR", "BP", "A1", "A2", "beta", "se", "p"])
    df = df[(df["p"] > 0) & (df["p"] <= 1) & (df["se"] > 0)]
    return df


def dedupe_leads(leads: pd.DataFrame, merge_distance_kb: int) -> pd.DataFrame:
    if leads.empty:
        return leads
    merge_bp = merge_distance_kb * 1000
    leads = leads.sort_values(["CHR", "LEAD_P", "LEAD_BP"]).reset_index(drop=True)
    kept = []
    for row in leads.itertuples(index=False):
        conflict = False
        for prev in kept:
            if row.CHR == prev["CHR"] and abs(row.LEAD_BP - prev["LEAD_BP"]) <= merge_bp:
                conflict = True
                break
        if not conflict:
            kept.append(
                {
                    "CHR": int(row.CHR),
                    "LEAD_BP": int(row.LEAD_BP),
                    "LEAD_SNP": row.LEAD_SNP,
                    "LEAD_P": float(row.LEAD_P),
                    "source": row.source,
                }
            )
    return pd.DataFrame(kept)


def get_beta_aligned(male_row: pd.Series | None, female_row: pd.Series | None) -> tuple[float | None, float | None, bool | None]:
    if male_row is None or female_row is None:
        return None, None, None
    beta_m = float(male_row["beta"])
    beta_f = float(female_row["beta"])
    if male_row["A1"] == female_row["A1"] and male_row["A2"] == female_row["A2"]:
        aligned_f = beta_f
    elif male_row["A1"] == female_row["A2"] and male_row["A2"] == female_row["A1"]:
        aligned_f = -beta_f
    else:
        return beta_m, None, None
    discordant = np.sign(beta_m) != np.sign(aligned_f)
    return beta_m, aligned_f, bool(discordant)


def classify_locus(
    p_male: float,
    p_female: float,
    sign_discordant: bool | None,
    p_gws: float,
    p_suggestive: float,
    p_null: float,
) -> str:
    if sign_discordant is True:
        return "discordant"
    if p_male < p_gws and p_female < p_gws:
        return "shared"
    if p_female < p_gws and p_male < p_suggestive:
        return "female-amplified"
    if p_male < p_gws and p_female < p_suggestive:
        return "male-amplified"
    if p_female < p_gws and p_male > p_null:
        return "female-specific"
    if p_male < p_gws and p_female > p_null:
        return "male-specific"
    return "suggestive-other"


def main() -> None:
    args = parse_args()
    args.repo_root = args.repo_root.resolve()
    args.define_loci_script = args.define_loci_script.resolve()
    args.out_dir = args.out_dir.resolve()
    loci_dir = args.out_dir / "loci_definitions"
    loci_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    all_rows: list[pd.DataFrame] = []

    for coding in ["SI", "SI_continuous"]:
        coding_label = "binary" if coding == "SI" else "continuous"
        for trait_idx, trait_name in TRAITS.items():
            paths = mtag_paths(args.repo_root, coding, trait_idx)
            male_df = read_sumstats(paths["male"])
            female_df = read_sumstats(paths["female"])

            # Run/collect loci for male, female, combined.
            lead_tables = []
            for sex_key in ["male", "female", "combined"]:
                out_prefix = loci_dir / f"{coding}__{sex_key}__trait{trait_idx}__{trait_name}"
                run_define_loci(args.define_loci_script, paths[sex_key], out_prefix, args.force_rerun_define)
                lead_path = out_prefix.with_suffix(".lead_snps.tsv")
                if not lead_path.exists():
                    continue
                lead_df = pd.read_csv(lead_path, sep="\t")
                if lead_df.empty:
                    continue
                lead_df["source"] = sex_key
                lead_tables.append(lead_df[["CHR", "BP", "SNP", "P", "source"]].rename(columns={"BP": "LEAD_BP", "SNP": "LEAD_SNP", "P": "LEAD_P"}))

            if not lead_tables:
                continue
            combined_leads = pd.concat(lead_tables, ignore_index=True).drop_duplicates(["LEAD_SNP"])
            loci = dedupe_leads(combined_leads, merge_distance_kb=args.merge_distance_kb)

            records = []
            window_bp = args.window_kb * 1000
            male_idx = male_df.set_index("SNP", drop=False)
            female_idx = female_df.set_index("SNP", drop=False)

            for row in loci.itertuples(index=False):
                chr_ = int(row.CHR)
                start = max(1, int(row.LEAD_BP) - window_bp)
                end = int(row.LEAD_BP) + window_bp
                male_window = male_df[(male_df["CHR"] == chr_) & (male_df["BP"] >= start) & (male_df["BP"] <= end)]
                female_window = female_df[(female_df["CHR"] == chr_) & (female_df["BP"] >= start) & (female_df["BP"] <= end)]

                p_male = float(male_window["p"].min()) if not male_window.empty else 1.0
                p_female = float(female_window["p"].min()) if not female_window.empty else 1.0

                male_row = male_idx.loc[row.LEAD_SNP] if row.LEAD_SNP in male_idx.index else None
                female_row = female_idx.loc[row.LEAD_SNP] if row.LEAD_SNP in female_idx.index else None
                if isinstance(male_row, pd.DataFrame):
                    male_row = male_row.iloc[0]
                if isinstance(female_row, pd.DataFrame):
                    female_row = female_row.iloc[0]

                beta_m, beta_f_aligned, sign_discordant = get_beta_aligned(male_row, female_row)
                category = classify_locus(
                    p_male=p_male,
                    p_female=p_female,
                    sign_discordant=sign_discordant,
                    p_gws=args.p_gws,
                    p_suggestive=args.p_suggestive,
                    p_null=args.p_null,
                )
                records.append(
                    {
                        "coding": coding_label,
                        "coding_raw": coding,
                        "trait_index": trait_idx,
                        "trait": trait_name,
                        "lead_snp": row.LEAD_SNP,
                        "chr": chr_,
                        "lead_bp": int(row.LEAD_BP),
                        "window_start": start,
                        "window_end": end,
                        "source_lead_set": row.source,
                        "p_male_min_window": p_male,
                        "p_female_min_window": p_female,
                        "beta_male_lead": beta_m,
                        "beta_female_lead_aligned": beta_f_aligned,
                        "sign_discordant": sign_discordant,
                        "locus_category": category,
                    }
                )

            out_df = pd.DataFrame(records).sort_values(["chr", "lead_bp"])
            per_combo = args.out_dir / f"{coding}__trait{trait_idx}__{trait_name}.locus_classification.tsv"
            out_df.to_csv(per_combo, sep="\t", index=False)
            all_rows.append(out_df)

            counts = out_df["locus_category"].value_counts().to_dict()
            summary_row = {
                "coding": coding_label,
                "coding_raw": coding,
                "trait_index": trait_idx,
                "trait": trait_name,
                "n_loci_total": int(len(out_df)),
            }
            for key in [
                "shared",
                "female-amplified",
                "male-amplified",
                "female-specific",
                "male-specific",
                "discordant",
                "suggestive-other",
            ]:
                summary_row[f"n_{key.replace('-', '_')}"] = int(counts.get(key, 0))
            summary_rows.append(summary_row)

    summary = pd.DataFrame(summary_rows).sort_values(["coding", "trait_index"])
    summary_tsv = args.out_dir / "locus_classification_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)

    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True)
        combined.to_csv(args.out_dir / "locus_classification_all.tsv", sep="\t", index=False)

    print(f"Wrote: {summary_tsv}")


if __name__ == "__main__":
    main()
