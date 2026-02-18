#!/usr/bin/env python3
"""
Convert MTAG output files into downstream tool-specific input formats.

Inputs:
  3_MTAG/**/results/*_trait_*.txt

Outputs:
  - LDSC input tables (+ optional ldsc munge_sumstats output)
  - MAGMA SNP location/p-value tables
  - S-PrediXcan harmonized GWAS tables
  - A manifest describing all generated files
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

SI_TRAITS = {
    1: "AbilityToConfide",
    2: "FreqSoc",
    3: "Loneliness",
}

MRI_TRAITS = {
    1: "FA",
    2: "MD",
    3: "MO",
    4: "OD",
}

REQUIRED_COLUMNS = {
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
}


@dataclass
class FileMeta:
    source_file: Path
    analysis_group: str
    population: str
    coding: str
    trait_index: int
    trait_name: str
    analysis_id: str


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    repo_root_default = script_path.parents[3]
    out_dir_default = repo_root_default / "5_Functional_Genomics" / "0_munge" / "output"
    mtag_root_default = repo_root_default / "3_MTAG"
    ldsc_munge_default = (
        repo_root_default
        / "5_Functional_Genomics"
        / "reference_data"
        / "ldsc"
        / "ldsc"
        / "munge_sumstats.py"
    )

    parser = argparse.ArgumentParser(
        description="Munge MTAG outputs for LDSC, MAGMA, and S-PrediXcan."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=repo_root_default,
        help="Path to UKBB-SI-Genetics repository root.",
    )
    parser.add_argument(
        "--mtag-root",
        type=Path,
        default=mtag_root_default,
        help="Path to the 3_MTAG directory.",
    )
    parser.add_argument(
        "--input-glob",
        default="**/results/*_trait_*.txt",
        help="Glob pattern under --mtag-root for MTAG output files.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=out_dir_default,
        help="Output directory for munged files.",
    )
    parser.add_argument(
        "--min-maf",
        type=float,
        default=0.0,
        help="Minimum minor allele frequency filter. Use 0 to disable.",
    )
    parser.add_argument(
        "--run-ldsc-munge",
        action="store_true",
        help="Run LDSC munge_sumstats.py to generate *.sumstats.gz files.",
    )
    parser.add_argument(
        "--ldsc-munge-script",
        type=Path,
        default=ldsc_munge_default,
        help="Path to LDSC munge_sumstats.py script.",
    )
    parser.add_argument(
        "--hm3-snplist",
        type=Path,
        default=None,
        help="Optional HapMap3 snplist path passed to LDSC munge with --merge-alleles.",
    )
    parser.add_argument(
        "--python-bin",
        default=sys.executable,
        help="Python binary to use when running LDSC munge script.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def infer_trait_name(analysis_group: str, trait_index: int) -> str:
    if analysis_group.startswith("SI"):
        return SI_TRAITS.get(trait_index, f"Trait{trait_index}")
    if analysis_group.startswith("MRI"):
        return MRI_TRAITS.get(trait_index, f"Trait{trait_index}")
    return f"Trait{trait_index}"


def sanitize_token(token: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", token)


def parse_file_meta(path: Path, mtag_root: Path) -> FileMeta:
    rel = path.relative_to(mtag_root)
    # Expected shape: ANALYSIS/POPULATION/results/file_trait_N.txt
    if len(rel.parts) < 4:
        raise ValueError(f"Unexpected MTAG path shape: {path}")

    analysis_group = rel.parts[0]
    population = rel.parts[1]
    coding = "continuous" if "continuous" in analysis_group.lower() else "binary"

    match = re.search(r"_trait_(\d+)\.txt$", path.name)
    if not match:
        raise ValueError(f"Could not parse trait index from filename: {path.name}")
    trait_index = int(match.group(1))
    trait_name = infer_trait_name(analysis_group, trait_index)

    analysis_id = "__".join(
        [
            sanitize_token(analysis_group),
            sanitize_token(population),
            f"trait{trait_index}",
            sanitize_token(trait_name),
        ]
    )

    return FileMeta(
        source_file=path,
        analysis_group=analysis_group,
        population=population,
        coding=coding,
        trait_index=trait_index,
        trait_name=trait_name,
        analysis_id=analysis_id,
    )


def discover_input_files(mtag_root: Path, pattern: str) -> list[Path]:
    return sorted([p for p in mtag_root.glob(pattern) if p.is_file()])


def _safe_pvalues(series: pd.Series) -> pd.Series:
    eps = np.nextafter(0, 1, dtype=np.float64)
    pvals = pd.to_numeric(series, errors="coerce").astype(float)
    pvals = pvals.clip(lower=eps, upper=1.0)
    return pvals


def _safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def read_and_filter_input(path: Path, min_maf: float) -> tuple[pd.DataFrame, int]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")

    before_n = len(df)

    # Normalize and validate key columns.
    df["SNP"] = df["SNP"].astype(str)
    df["CHR"] = _safe_numeric(df["CHR"])
    df["BP"] = _safe_numeric(df["BP"])
    df["N"] = _safe_numeric(df["N"])
    df["FRQ"] = _safe_numeric(df["FRQ"])
    df["mtag_beta"] = _safe_numeric(df["mtag_beta"])
    df["mtag_se"] = _safe_numeric(df["mtag_se"])
    df["mtag_z"] = _safe_numeric(df["mtag_z"])
    df["mtag_pval"] = _safe_pvalues(df["mtag_pval"])

    valid = (
        df["SNP"].notna()
        & df["SNP"].ne("nan")
        & df["CHR"].notna()
        & df["BP"].notna()
        & df["N"].notna()
        & df["FRQ"].notna()
        & df["mtag_beta"].notna()
        & df["mtag_se"].notna()
        & df["mtag_z"].notna()
        & df["mtag_pval"].notna()
        & df["A1"].notna()
        & df["A2"].notna()
    )
    df = df.loc[valid].copy()

    # Standard frequency/MAF QC.
    df = df[(df["FRQ"] > 0.0) & (df["FRQ"] < 1.0)]
    if min_maf > 0:
        maf = np.minimum(df["FRQ"], 1.0 - df["FRQ"])
        df = df[maf >= min_maf]

    # Remove duplicate SNP IDs, keeping the smallest p-value row.
    df = df.sort_values("mtag_pval", ascending=True).drop_duplicates("SNP", keep="first")
    df = df.sort_values(["CHR", "BP"]).reset_index(drop=True)
    return df, before_n


def build_ldsc_table(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "SNP": df["SNP"],
            "A1": df["A1"],
            "A2": df["A2"],
            "N": df["N"].astype(int),
            "Z": df["mtag_z"],
            "P": df["mtag_pval"],
            "FRQ": df["FRQ"],
            "CHR": df["CHR"].astype(int),
            "BP": df["BP"].astype(int),
            "BETA": df["mtag_beta"],
            "SE": df["mtag_se"],
        }
    )


def build_magma_table(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "SNP": df["SNP"],
            "CHR": df["CHR"].astype(int),
            "BP": df["BP"].astype(int),
            "P": df["mtag_pval"],
            "N": df["N"].astype(int),
        }
    )


def build_spredixcan_table(df: pd.DataFrame) -> pd.DataFrame:
    chr_int = df["CHR"].astype(int).astype(str)
    bp_int = df["BP"].astype(int).astype(str)
    a1 = df["A1"].astype(str)
    a2 = df["A2"].astype(str)
    panel_variant_id = "chr" + chr_int + "_" + bp_int + "_" + a1 + "_" + a2 + "_b37"

    return pd.DataFrame(
        {
            "rsid": df["SNP"],
            "chromosome": df["CHR"].astype(int),
            "position": df["BP"].astype(int),
            "effect_allele": df["A1"],
            "non_effect_allele": df["A2"],
            "beta": df["mtag_beta"],
            "se": df["mtag_se"],
            "pvalue": df["mtag_pval"],
            "zscore": df["mtag_z"],
            "sample_size": df["N"].astype(int),
            "effect_allele_frequency": df["FRQ"],
            "panel_variant_id": panel_variant_id,
        }
    )


def maybe_run_ldsc_munge(
    args: argparse.Namespace,
    ldsc_input_path: Path,
    out_prefix: Path,
) -> Path | None:
    if not args.run_ldsc_munge:
        return None

    if not args.ldsc_munge_script.exists():
        raise FileNotFoundError(
            f"--run-ldsc-munge was requested, but script is missing: {args.ldsc_munge_script}"
        )

    cmd = [
        args.python_bin,
        str(args.ldsc_munge_script),
        "--sumstats",
        str(ldsc_input_path),
        "--snp",
        "SNP",
        "--a1",
        "A1",
        "--a2",
        "A2",
        "--N-col",
        "N",
        "--signed-sumstats",
        "Z,0",
        "--p",
        "P",
        "--frq",
        "FRQ",
        "--chunksize",
        "500000",
        "--out",
        str(out_prefix),
    ]
    if args.hm3_snplist is not None:
        cmd.extend(["--merge-alleles", str(args.hm3_snplist)])

    subprocess.run(cmd, check=True)
    out_path = out_prefix.with_suffix(".sumstats.gz")
    return out_path


def ensure_dirs(out_dir: Path) -> dict[str, Path]:
    paths = {
        "ldsc_input": out_dir / "ldsc_input",
        "ldsc_munged": out_dir / "ldsc_munged",
        "magma": out_dir / "magma",
        "spredixcan": out_dir / "spredixcan",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths


def write_manifest(rows: Iterable[dict[str, str]], manifest_path: Path) -> None:
    fieldnames = [
        "analysis_id",
        "analysis_group",
        "coding",
        "population",
        "trait_index",
        "trait_name",
        "source_file",
        "n_input_rows",
        "n_output_rows",
        "ldsc_input",
        "ldsc_sumstats",
        "magma_input",
        "spredixcan_input",
    ]
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    args.repo_root = args.repo_root.resolve()
    args.mtag_root = args.mtag_root.resolve()
    args.out_dir = args.out_dir.resolve()
    if args.hm3_snplist is not None:
        args.hm3_snplist = args.hm3_snplist.resolve()
    args.ldsc_munge_script = args.ldsc_munge_script.resolve()

    if not args.mtag_root.exists():
        raise FileNotFoundError(f"MTAG root does not exist: {args.mtag_root}")

    files = discover_input_files(args.mtag_root, args.input_glob)
    if not files:
        raise FileNotFoundError(
            f"No files matched '{args.input_glob}' under {args.mtag_root}"
        )

    out_paths = ensure_dirs(args.out_dir)
    manifest_rows: list[dict[str, str]] = []

    for fpath in files:
        meta = parse_file_meta(fpath, args.mtag_root)
        print(f"[process] {fpath}")

        df, input_rows = read_and_filter_input(fpath, min_maf=args.min_maf)
        output_rows = len(df)
        if output_rows == 0:
            print(f"[warn] No rows passed QC for {fpath}; skipping.")
            continue

        ldsc_table = build_ldsc_table(df)
        magma_table = build_magma_table(df)
        spredixcan_table = build_spredixcan_table(df)

        ldsc_input_path = out_paths["ldsc_input"] / f"{meta.analysis_id}.ldsc.tsv.gz"
        magma_path = out_paths["magma"] / f"{meta.analysis_id}.genes.raw"
        spredixcan_path = (
            out_paths["spredixcan"] / f"{meta.analysis_id}.spredixcan.tsv.gz"
        )

        if not args.overwrite:
            for output_path in [ldsc_input_path, magma_path, spredixcan_path]:
                if output_path.exists():
                    raise FileExistsError(
                        f"Output exists (use --overwrite to replace): {output_path}"
                    )

        ldsc_table.to_csv(ldsc_input_path, sep="\t", index=False, compression="gzip")
        magma_table.to_csv(magma_path, sep="\t", index=False)
        spredixcan_table.to_csv(
            spredixcan_path, sep="\t", index=False, compression="gzip"
        )

        ldsc_sumstats_path: Path | None = None
        if args.run_ldsc_munge:
            ldsc_out_prefix = out_paths["ldsc_munged"] / meta.analysis_id
            ldsc_sumstats_path = maybe_run_ldsc_munge(
                args=args,
                ldsc_input_path=ldsc_input_path,
                out_prefix=ldsc_out_prefix,
            )

        manifest_rows.append(
            {
                "analysis_id": meta.analysis_id,
                "analysis_group": meta.analysis_group,
                "coding": meta.coding,
                "population": meta.population,
                "trait_index": str(meta.trait_index),
                "trait_name": meta.trait_name,
                "source_file": str(meta.source_file),
                "n_input_rows": str(input_rows),
                "n_output_rows": str(output_rows),
                "ldsc_input": str(ldsc_input_path),
                "ldsc_sumstats": "" if ldsc_sumstats_path is None else str(ldsc_sumstats_path),
                "magma_input": str(magma_path),
                "spredixcan_input": str(spredixcan_path),
            }
        )

    manifest_path = args.out_dir / "manifest.tsv"
    write_manifest(manifest_rows, manifest_path)
    print(f"[done] Wrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()
