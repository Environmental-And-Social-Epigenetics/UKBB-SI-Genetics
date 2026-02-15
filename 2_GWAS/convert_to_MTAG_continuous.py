#!/usr/bin/env python3
"""
Convert continuous-coded BOLT-LMM output files to MTAG format with rsID mapping.
Processes all population stratifications (EUR_MM, EUR_Male, EUR_Female).
"""

import os
import pandas as pd

from convert_to_MTAG import load_rsid_mapping, convert_bolt_to_mtag


def main():
    # Paths
    srcdir = "/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness"

    # rsID mapping file
    annot_file = "/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/geno/ukb_genoHM3.annot.genesymbol_mapped.pvar.gz"

    # Analysis configuration
    phenotypes = ["Loneliness", "FreqSoc", "AbilityToConfide"]
    populations = ["EUR_MM", "EUR_Male", "EUR_Female"]
    covar_set = "Day_NoPCs"
    pheno_prefix = "isolation_run_continuous"
    results_subdir = "results_continuous"
    output_root = "mtag_results_continuous"

    print("=" * 70)
    print("BOLT-LMM to MTAG Format Conversion (Continuous)")
    print("Social Isolation Phenotypes - Population Stratified")
    print("=" * 70)
    print()

    # Load rsID mapping (once for all)
    rsid_lookup = load_rsid_mapping(annot_file)
    print()

    # Process each population
    for population in populations:
        print("=" * 70)
        print(f"Processing {population}")
        print("=" * 70)

        # Population-specific phenotype file for sample sizes
        pheno_file = f"{srcdir}/{pheno_prefix}.{population}.tsv.gz"

        # Check if phenotype file exists
        if not os.path.exists(pheno_file):
            print(f"Phenotype file not found: {pheno_file}")
            print(f"Skipping {population}")
            print()
            continue

        # Get sample sizes for each phenotype
        print(f"Determining sample sizes from: {os.path.basename(pheno_file)}")
        sample_sizes = {}
        df_pheno = pd.read_csv(pheno_file, sep="\t", compression="gzip")
        for pheno in phenotypes:
            n = df_pheno[pheno].notna().sum()
            sample_sizes[pheno] = n
            print(f"  {pheno}: {n:,} samples")

        # Check if results directory exists
        results_dir = f"{srcdir}/{results_subdir}/{covar_set}/{population}"
        if not os.path.exists(results_dir):
            print(f"Results directory not found: {results_dir}")
            print(f"Skipping {population}")
            print()
            continue

        # Create output directory for this population
        output_dir = f"{srcdir}/{output_root}/{population}"
        os.makedirs(output_dir, exist_ok=True)
        print(f"\nOutput directory: {output_dir}")

        # Process each phenotype
        converted_count = 0
        for pheno in phenotypes:
            bolt_file = f"{results_dir}/bolt_{pheno}.{covar_set}.stats.gz"

            if not os.path.exists(bolt_file):
                print(f"\n  File not found: {bolt_file}")
                print(f"  Skipping {pheno}")
                continue

            output_file = f"{output_dir}/{pheno}.{covar_set}.mtag.sumstats.txt"

            try:
                convert_bolt_to_mtag(
                    bolt_file=bolt_file,
                    rsid_lookup=rsid_lookup,
                    trait_name=f"{pheno} ({population}, continuous)",
                    sample_size=sample_sizes[pheno],
                    output_file=output_file,
                )
                converted_count += 1
            except Exception as exc:
                print(f"\n  Error processing {pheno}: {exc}")
                continue

        print(f"\nConverted {converted_count}/{len(phenotypes)} phenotypes for {population}")
        print()

    print("=" * 70)
    print("Continuous conversion complete")
    print("=" * 70)
    print(f"Output directory structure: {srcdir}/{output_root}/")
    print()


if __name__ == "__main__":
    main()
