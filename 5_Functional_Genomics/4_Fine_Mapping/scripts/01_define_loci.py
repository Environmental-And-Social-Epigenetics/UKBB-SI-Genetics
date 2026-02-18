#!/usr/bin/env python3
"""
Define genome-wide significant loci from MTAG/GWAS summary statistics.

Preferred mode:
  - PLINK clumping (r2 + distance)

Fallback mode:
  - Distance-based lead SNP pruning within +/- window_kb
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class ColMap:
    snp: str
    chr: str
    bp: str
    p: str


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    fg_root = script_path.parents[2]
    out_default = fg_root / "4_Fine_Mapping" / "output" / "01_defined_loci" / "loci"

    parser = argparse.ArgumentParser(description="Define independent significant loci.")
    parser.add_argument("--sumstats", type=Path, required=True, help="Input summary statistics file.")
    parser.add_argument("--out-prefix", type=Path, default=out_default)
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--window-kb", type=int, default=500)
    parser.add_argument("--plink-bin", default="plink")
    parser.add_argument("--plink-bfile", type=Path, default=None, help="PLINK reference prefix for clumping.")
    parser.add_argument("--clump-r2", type=float, default=0.1)
    parser.add_argument("--clump-p2", type=float, default=1e-4)
    parser.add_argument(
        "--allow-distance-fallback",
        action="store_true",
        help="Use distance-based pruning if PLINK clumping is unavailable.",
    )
    return parser.parse_args()


def detect_columns(df: pd.DataFrame) -> ColMap:
    snp = "SNP" if "SNP" in df.columns else "snpid" if "snpid" in df.columns else None
    chr_col = "CHR" if "CHR" in df.columns else "chr" if "chr" in df.columns else None
    bp_col = "BP" if "BP" in df.columns else "bpos" if "bpos" in df.columns else None
    p_col = (
        "mtag_pval"
        if "mtag_pval" in df.columns
        else "pval"
        if "pval" in df.columns
        else "P"
        if "P" in df.columns
        else None
    )

    if snp is None or chr_col is None or bp_col is None or p_col is None:
        raise ValueError(
            f"Could not infer required columns from: {list(df.columns)}. "
            "Need SNP/snpid, CHR/chr, BP/bpos, and p-value column."
        )
    return ColMap(snp=snp, chr=chr_col, bp=bp_col, p=p_col)


def read_sumstats(path: Path) -> tuple[pd.DataFrame, ColMap]:
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception:
        df = pd.read_csv(path, delim_whitespace=True, low_memory=False)

    colmap = detect_columns(df)
    df = df[[colmap.snp, colmap.chr, colmap.bp, colmap.p]].copy()
    df.columns = ["SNP", "CHR", "BP", "P"]
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna()
    df = df[(df["P"] > 0) & (df["P"] <= 1)]
    df = df.sort_values("P", ascending=True).drop_duplicates("SNP", keep="first")
    return df, colmap


def distance_prune(sig_df: pd.DataFrame, window_kb: int) -> pd.DataFrame:
    selected_rows = []
    window_bp = window_kb * 1000

    for row in sig_df.sort_values("P", ascending=True).itertuples(index=False):
        chr_sel = row.CHR
        bp_sel = row.BP

        conflict = False
        for s in selected_rows:
            if s["CHR"] == chr_sel and abs(s["BP"] - bp_sel) <= window_bp:
                conflict = True
                break
        if not conflict:
            selected_rows.append({"SNP": row.SNP, "CHR": chr_sel, "BP": bp_sel, "P": row.P})

    return pd.DataFrame(selected_rows)


def parse_plink_clumped(clumped_path: Path) -> pd.DataFrame:
    if not clumped_path.exists() or clumped_path.stat().st_size == 0:
        return pd.DataFrame(columns=["SNP", "CHR", "BP", "P"])

    df = pd.read_csv(clumped_path, delim_whitespace=True, comment="#", low_memory=False)
    required = {"SNP", "CHR", "BP", "P"}
    if not required.issubset(df.columns):
        return pd.DataFrame(columns=["SNP", "CHR", "BP", "P"])
    out = df[["SNP", "CHR", "BP", "P"]].copy()
    out["CHR"] = pd.to_numeric(out["CHR"], errors="coerce")
    out["BP"] = pd.to_numeric(out["BP"], errors="coerce")
    out["P"] = pd.to_numeric(out["P"], errors="coerce")
    out = out.dropna()
    return out


def run_plink_clump(
    sig_df: pd.DataFrame,
    plink_bin: str,
    plink_bfile: Path,
    p1: float,
    p2: float,
    r2: float,
    kb: int,
) -> pd.DataFrame:
    if shutil.which(plink_bin) is None:
        raise FileNotFoundError(f"PLINK binary not found: {plink_bin}")
    if not Path(f"{plink_bfile}.bed").exists():
        raise FileNotFoundError(f"PLINK reference missing: {plink_bfile}.bed")

    with tempfile.TemporaryDirectory(prefix="define_loci_") as tmpd:
        tmpdir = Path(tmpd)
        clump_input = tmpdir / "clump_input.tsv"
        clump_out_prefix = tmpdir / "plink_clump"
        sig_df[["SNP", "P"]].to_csv(clump_input, sep="\t", index=False)

        cmd = [
            plink_bin,
            "--bfile",
            str(plink_bfile),
            "--clump",
            str(clump_input),
            "--clump-snp-field",
            "SNP",
            "--clump-field",
            "P",
            "--clump-p1",
            str(p1),
            "--clump-p2",
            str(p2),
            "--clump-r2",
            str(r2),
            "--clump-kb",
            str(kb),
            "--out",
            str(clump_out_prefix),
        ]
        subprocess.run(cmd, check=True)
        return parse_plink_clumped(clump_out_prefix.with_suffix(".clumped"))


def to_loci_table(leads: pd.DataFrame, window_kb: int, method: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if leads.empty:
        empty_loci = pd.DataFrame(
            columns=["locus_id", "CHR", "START", "END", "LEAD_SNP", "LEAD_BP", "LEAD_P", "method"]
        )
        empty_leads = pd.DataFrame(
            columns=["locus_id", "SNP", "CHR", "BP", "P", "method"]
        )
        return empty_loci, empty_leads

    leads = leads.sort_values(["CHR", "BP"]).reset_index(drop=True)
    lead_rows = []
    locus_rows = []
    window_bp = window_kb * 1000

    for idx, row in enumerate(leads.itertuples(index=False), start=1):
        locus_id = f"locus_{idx:04d}_chr{int(row.CHR)}_{int(row.BP)}"
        start = max(1, int(row.BP) - window_bp)
        end = int(row.BP) + window_bp
        lead_rows.append(
            {
                "locus_id": locus_id,
                "SNP": row.SNP,
                "CHR": int(row.CHR),
                "BP": int(row.BP),
                "P": float(row.P),
                "method": method,
            }
        )
        locus_rows.append(
            {
                "locus_id": locus_id,
                "CHR": int(row.CHR),
                "START": start,
                "END": end,
                "LEAD_SNP": row.SNP,
                "LEAD_BP": int(row.BP),
                "LEAD_P": float(row.P),
                "method": method,
            }
        )

    return pd.DataFrame(locus_rows), pd.DataFrame(lead_rows)


def main() -> None:
    args = parse_args()
    args.sumstats = args.sumstats.resolve()
    args.out_prefix = args.out_prefix.resolve()
    if args.plink_bfile is not None:
        args.plink_bfile = args.plink_bfile.resolve()

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    df, _ = read_sumstats(args.sumstats)
    sig_df = df[df["P"] < args.p_threshold].copy()

    sig_out = args.out_prefix.with_suffix(".sig_snps.tsv")
    sig_df.to_csv(sig_out, sep="\t", index=False)

    method = "distance"
    leads = pd.DataFrame(columns=["SNP", "CHR", "BP", "P"])

    if len(sig_df) > 0 and args.plink_bfile is not None:
        try:
            leads = run_plink_clump(
                sig_df=sig_df,
                plink_bin=args.plink_bin,
                plink_bfile=args.plink_bfile,
                p1=args.p_threshold,
                p2=args.clump_p2,
                r2=args.clump_r2,
                kb=args.window_kb,
            )
            method = "plink_clump"
        except Exception as exc:
            if not args.allow_distance_fallback:
                raise RuntimeError(
                    f"PLINK clumping failed ({exc}). Re-run with --allow-distance-fallback."
                ) from exc
            print(f"[warn] PLINK clumping failed, falling back to distance pruning: {exc}")
            leads = distance_prune(sig_df, window_kb=args.window_kb)
            method = "distance_fallback"
    elif len(sig_df) > 0:
        leads = distance_prune(sig_df, window_kb=args.window_kb)
        method = "distance"

    loci_df, lead_df = to_loci_table(leads=leads, window_kb=args.window_kb, method=method)
    loci_out = args.out_prefix.with_suffix(".loci.tsv")
    lead_out = args.out_prefix.with_suffix(".lead_snps.tsv")

    loci_df.to_csv(loci_out, sep="\t", index=False)
    lead_df.to_csv(lead_out, sep="\t", index=False)

    print(f"Input summary stats: {args.sumstats}")
    print(f"Genome-wide significant SNPs: {len(sig_df)}")
    print(f"Lead loci: {len(lead_df)} (method={method})")
    print(f"Wrote: {sig_out}")
    print(f"Wrote: {lead_out}")
    print(f"Wrote: {loci_out}")


if __name__ == "__main__":
    main()
