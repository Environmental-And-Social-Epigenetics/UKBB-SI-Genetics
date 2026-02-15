#!/usr/bin/env python3
"""
Filter phenotype and covariate files to EUR_MM, EUR_Male, EUR_Female populations.
More robust than bash/awk approach - matches working repo's filter_to_EUR_python.py
"""

import gzip
import sys
import os

def read_keep_ids(keep_file):
    """Read IIDs from keep file"""
    keep_ids = set()
    with open(keep_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                iid = parts[1]  # Second column is IID
                keep_ids.add(iid)
    return keep_ids

def filter_file(input_file, output_file, keep_ids, ensure_fid_iid_header=False):
    """Filter a gzipped TSV file to specified samples only"""
    n_in = 0
    n_out = 0
    
    with gzip.open(input_file, 'rt') as fin, gzip.open(output_file, 'wt') as fout:
        # Read header
        header = fin.readline()
        header_parts = header.strip().split('\t')
        
        # BOLT-LMM requires "FID IID" as first two columns
        if ensure_fid_iid_header:
            if header_parts[0] != 'FID' or header_parts[1] != 'IID':
                # Fix header if needed
                print(f"  NOTE: Adjusting header to start with 'FID IID'", file=sys.stderr)
                print(f"    Original: {header_parts[0]} {header_parts[1]}", file=sys.stderr)
                header_parts[0] = 'FID'
                header_parts[1] = 'IID'
                header = '\t'.join(header_parts) + '\n'
        
        fout.write(header)
        
        # Process data rows
        for line in fin:
            n_in += 1
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                iid = parts[1]  # IID is second column
                if iid in keep_ids:
                    fout.write(line)
                    n_out += 1
            
            # Progress indicator
            if n_in % 100000 == 0:
                print(f"  Processed {n_in} samples, kept {n_out}...", file=sys.stderr)
    
    return n_in, n_out

def filter_population(population_name):
    """Filter phenotype and covariate files for one population"""
    # Paths
    ukb21942_d = '/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942'
    srcdir = '/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness'
    
    keep_file = f'{ukb21942_d}/sqc/population.20220316/{population_name}.keep'
    pheno_in = f'{ukb21942_d}/pheno/isolation_run_control.tsv.gz'
    covar_in = f'{ukb21942_d}/sqc/sqc.20220316.tsv.gz'
    
    pheno_out = f'{srcdir}/isolation_run_control.{population_name}.tsv.gz'
    covar_out = f'{srcdir}/sqc.{population_name}.tsv.gz'
    
    print("=" * 60)
    print(f"Filter to {population_name}")
    print("=" * 60)
    print()
    
    # Check if keep file exists
    if not os.path.exists(keep_file):
        print(f"ERROR: Keep file not found: {keep_file}", file=sys.stderr)
        return False
    
    # Read sample IDs
    print(f"Reading sample IDs from: {keep_file}")
    keep_ids = read_keep_ids(keep_file)
    print(f"  Samples to keep: {len(keep_ids)}")
    print()
    
    # Filter phenotype file
    print(f"Filtering phenotype file...")
    print(f"  Input:  {pheno_in}")
    print(f"  Output: {pheno_out}")
    n_pheno_in, n_pheno_out = filter_file(pheno_in, pheno_out, keep_ids, ensure_fid_iid_header=True)
    print(f"  ✓ Complete: {n_pheno_in} input → {n_pheno_out} {population_name} samples")
    print()
    
    # Filter covariate file (BOLT requires "FID IID" header)
    print(f"Filtering covariate file...")
    print(f"  Input:  {covar_in}")
    print(f"  Output: {covar_out}")
    n_covar_in, n_covar_out = filter_file(covar_in, covar_out, keep_ids, ensure_fid_iid_header=True)
    print(f"  ✓ Complete: {n_covar_in} input → {n_covar_out} {population_name} samples")
    print()
    
    print("=" * 60)
    print(f"Filtering Complete for {population_name}!")
    print("=" * 60)
    print()
    print("Summary:")
    print(f"  IDs in keep file:         {len(keep_ids)}")
    print(f"  Phenotype samples:        {n_pheno_out}")
    print(f"  Covariate samples:        {n_covar_out}")
    print()
    
    if n_pheno_out != len(keep_ids):
        print(f"  Note: {len(keep_ids) - n_pheno_out} samples missing from phenotype file")
        print("        (normal if some samples have missing data)")
    
    if n_covar_out != len(keep_ids):
        print(f"  Note: {len(keep_ids) - n_covar_out} samples missing from covariate file")
        print("        (normal if some samples have missing data)")
    print()
    
    return True

def main():
    """Filter all three populations"""
    populations = ['EUR_MM', 'EUR_Male', 'EUR_Female']
    
    print()
    print("=" * 60)
    print("BOLT-LMM Loneliness: Filter All Populations")
    print("=" * 60)
    print()
    
    success_count = 0
    for pop in populations:
        if filter_population(pop):
            success_count += 1
        print()
    
    print("=" * 60)
    print(f"ALL FILTERING COMPLETED: {success_count}/{len(populations)} populations")
    print("=" * 60)
    print()
    
    if success_count == len(populations):
        print("Next steps:")
        print("  1. Verify files: ls -lh *.EUR*.tsv.gz")
        print("  2. Run test: sbatch 0b_test_run.sbatch.sh")
        print("  3. If test passes: sbatch 1_run_bolt_lmm.sbatch.sh")
    else:
        print("ERROR: Some populations failed to filter")
        return 1
    
    return 0

if __name__ == '__main__':
    sys.exit(main())

