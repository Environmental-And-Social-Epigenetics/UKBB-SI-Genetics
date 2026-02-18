#!/usr/bin/env python3
"""
Build an integrated gene prioritization table across post-GWAS modules.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description="Build integrated gene prioritization table.")
    parser.add_argument(
        "--spredixcan-prioritized",
        type=Path,
        default=fg_root / "3_S-PrediXcan" / "output" / "prioritized_genes.tsv",
    )
    parser.add_argument(
        "--coloc-summary",
        type=Path,
        default=fg_root / "5_Colocalization" / "output" / "summary" / "coloc_prioritized.tsv",
    )
    parser.add_argument(
        "--magma-dir",
        type=Path,
        default=fg_root / "2_MAGMA" / "output" / "02_gene_analysis",
    )
    parser.add_argument(
        "--susie-pip",
        type=Path,
        default=fg_root / "4_Fine_Mapping" / "output" / "03_susie" / "susie_pip.tsv",
    )
    parser.add_argument(
        "--snp-gene-map",
        type=Path,
        default=None,
        help="Optional SNP-to-gene mapping table with columns SNP and gene.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=fg_root / "6_Visualization" / "output" / "gene_prioritization_table.tsv",
    )
    parser.add_argument("--magma-alpha", type=float, default=2.7e-6)
    parser.add_argument("--transcriptome-alpha", type=float, default=2.5e-6)
    parser.add_argument("--coloc-threshold", type=float, default=0.8)
    parser.add_argument("--pip-threshold", type=float, default=0.1)
    return parser.parse_args()


def safe_read(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)


def best_col(columns: list[str], candidates: list[str]) -> str | None:
    for c in candidates:
        if c in columns:
            return c
    return None


def load_spredixcan_prior(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["gene"])
    df = safe_read(path)
    if "gene" not in df.columns:
        return pd.DataFrame(columns=["gene"])
    df["gene"] = df["gene"].astype(str).str.upper()
    return df


def load_coloc(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["gene"])
    df = safe_read(path)
    if "gene" not in df.columns:
        return pd.DataFrame(columns=["gene"])
    df["gene"] = df["gene"].astype(str).str.upper()
    rename_map = {}
    if "coloc_max_pp_h4" not in df.columns:
        if "PP.H4" in df.columns:
            rename_map["PP.H4"] = "coloc_max_pp_h4"
        elif "coloc_pp_h4" in df.columns:
            rename_map["coloc_pp_h4"] = "coloc_max_pp_h4"
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


def load_magma_gene_best(magma_dir: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(magma_dir.glob("*.genes.out")):
        try:
            df = pd.read_csv(f, sep=r"\s+")
        except Exception:
            continue
        if "GENE" not in df.columns or "P" not in df.columns:
            continue
        trait = f.name.replace(".genes.out", "")
        sub = df[["GENE", "P"]].copy()
        sub["GENE"] = sub["GENE"].astype(str).str.upper()
        sub["P"] = pd.to_numeric(sub["P"], errors="coerce")
        sub = sub.dropna()
        sub["magma_trait"] = trait
        rows.append(sub.rename(columns={"GENE": "gene", "P": "magma_p"}))
    if not rows:
        return pd.DataFrame(columns=["gene", "magma_p", "magma_trait"])
    out = pd.concat(rows, ignore_index=True)
    out = out.sort_values("magma_p", ascending=True).drop_duplicates("gene", keep="first")
    return out


def load_finemap_gene_support(pip_file: Path, snp_gene_map: Path | None) -> pd.DataFrame:
    if snp_gene_map is None or not pip_file.exists() or not snp_gene_map.exists():
        return pd.DataFrame(columns=["gene", "fine_map_max_pip", "fine_map_n_snps_high_pip"])

    pip = safe_read(pip_file)
    if "SNP" not in pip.columns or "PIP" not in pip.columns:
        return pd.DataFrame(columns=["gene", "fine_map_max_pip", "fine_map_n_snps_high_pip"])
    mapping = safe_read(snp_gene_map)
    if "SNP" not in mapping.columns or "gene" not in mapping.columns:
        return pd.DataFrame(columns=["gene", "fine_map_max_pip", "fine_map_n_snps_high_pip"])

    pip["SNP"] = pip["SNP"].astype(str)
    pip["PIP"] = pd.to_numeric(pip["PIP"], errors="coerce")
    mapping["SNP"] = mapping["SNP"].astype(str)
    mapping["gene"] = mapping["gene"].astype(str).str.upper()
    merged = pip.merge(mapping[["SNP", "gene"]], on="SNP", how="inner")
    if merged.empty:
        return pd.DataFrame(columns=["gene", "fine_map_max_pip", "fine_map_n_snps_high_pip"])

    merged["high_pip"] = merged["PIP"] >= 0.1
    out = (
        merged.groupby("gene", as_index=False)
        .agg(
            fine_map_max_pip=("PIP", "max"),
            fine_map_n_snps_high_pip=("high_pip", "sum"),
        )
    )
    return out


def score_row(row: pd.Series) -> float:
    score = 0.0
    for col, weight in [("magma_p", 1.0), ("spredixcan_min_p", 1.0), ("smultixcan_p", 1.2)]:
        val = row.get(col, np.nan)
        if pd.notna(val) and val > 0:
            score += weight * (-math.log10(float(val)))
    pp4 = row.get("coloc_max_pp_h4", np.nan)
    if pd.notna(pp4):
        score += 5.0 * float(pp4)
    pip = row.get("fine_map_max_pip", np.nan)
    if pd.notna(pip):
        score += 4.0 * float(pip)
    return score


def main() -> None:
    args = parse_args()

    sp = load_spredixcan_prior(args.spredixcan_prioritized)
    coloc = load_coloc(args.coloc_summary)
    magma = load_magma_gene_best(args.magma_dir)
    fm = load_finemap_gene_support(args.susie_pip, args.snp_gene_map)

    genes = set(sp.get("gene", [])) | set(coloc.get("gene", [])) | set(magma.get("gene", [])) | set(fm.get("gene", []))
    if not genes:
        raise RuntimeError("No genes found in supplied inputs.")

    merged = pd.DataFrame({"gene": sorted(genes)})
    if not sp.empty:
        merged = merged.merge(sp, on="gene", how="left")
    if not coloc.empty:
        merged = merged.merge(coloc, on="gene", how="left")
    if not magma.empty:
        merged = merged.merge(magma, on="gene", how="left")
    if not fm.empty:
        merged = merged.merge(fm, on="gene", how="left")

    merged["magma_sig"] = merged["magma_p"] < args.magma_alpha
    merged["transcriptome_sig"] = (
        (merged.get("spredixcan_min_p", np.nan) < args.transcriptome_alpha)
        | (merged.get("smultixcan_p", np.nan) < args.transcriptome_alpha)
    )
    merged["coloc_sig"] = merged.get("coloc_max_pp_h4", np.nan) >= args.coloc_threshold
    merged["fine_map_sig"] = merged.get("fine_map_max_pip", np.nan) >= args.pip_threshold
    merged["evidence_count"] = (
        merged["magma_sig"].fillna(False).astype(int)
        + merged["transcriptome_sig"].fillna(False).astype(int)
        + merged["coloc_sig"].fillna(False).astype(int)
        + merged["fine_map_sig"].fillna(False).astype(int)
    )
    merged["integrated_score"] = merged.apply(score_row, axis=1)
    merged = merged.sort_values(["evidence_count", "integrated_score"], ascending=[False, False])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
