#!/bin/bash
set -beEo pipefail

# Simplified BOLT-LMM script: runs one phenotype with one covariate set for one population
# Processes the FULL GENOME (no variant splitting)

if [ $# -ne 3 ]; then
    echo "Usage: $0 <phenotype> <covar_str> <keep_set>" >&2
    echo "Example: $0 Loneliness Day_NoPCs EUR_MM" >&2
    exit 1
fi

phenotype=$1
covar_str=$2
keep_set=$3

# Directories
REPODIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942"
SRCDIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness"
ukb21942_d='/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942'

echo "========================================"
echo "Running BOLT-LMM for ${phenotype}"
echo "Covariate model: ${covar_str}"
echo "Population: ${keep_set}"
echo "========================================"

# Output directory
out_dir="${SRCDIR}/results/${covar_str}/${keep_set}"
mkdir -p ${out_dir}

out_file="${out_dir}/bolt_${phenotype}.${covar_str}"

# Remove existing output files to ensure clean run
echo "Checking for existing output files..."
for ext in stats stats.gz log log.gz; do
    if [ -f "${out_file}.${ext}" ]; then
        echo "  Removing old file: ${out_file}.${ext}"
        rm -f "${out_file}.${ext}"
    fi
done
echo "✓ Ready for clean run"
echo ""

# Input files
genotype_bfile="${ukb21942_d}/geno/ukb_genoHM3/ukb_genoHM3_bed"
model_snps_file="${ukb21942_d}/geno/ukb_genoHM3/ukb_genoHM3_modelSNPs.txt"
ld_scores_file="/home/mabdel03/data/software/BOLT-LMM_v2.5/tables/LDSCORE.1000G_EUR.GRCh38.tab.gz"
genetic_map_file="/home/mabdel03/data/software/BOLT-LMM_v2.5/tables/genetic_map_hg19_withX.txt.gz"

# Use population-filtered files (created by filter_to_population.sh)
pheno_file_pop="${SRCDIR}/isolation_run_control.${keep_set}.tsv.gz"
covar_file_pop="${SRCDIR}/sqc.${keep_set}.tsv.gz"

# Set up covariates based on covar_str
# Note: BOLT needs separate arguments for each covariate (both categorical and quantitative)
if [ "${covar_str}" == "Day_NoPCs" ]; then
    # Quantitative: age only
    qcovar_col_args="--qCovarCol=age"
    # Categorical: sex and array
    covar_col_args="--covarCol=sex --covarCol=array"
else
    echo "ERROR: Unknown covar_str: ${covar_str}" >&2
    exit 1
fi

echo "Configuration:"
echo "  Genotype: ${genotype_bfile}"
echo "  Phenotype: ${phenotype}"
echo "  Phenotype file: ${pheno_file_pop}"
echo "  Covariate file: ${covar_file_pop}"
echo "  Quantitative covariates: age"
echo "  Categorical covariates: sex, array"
echo "  Model SNPs: ${model_snps_file} (~444K SNPs)"
echo "  Output: ${out_file}.stats"
echo ""

# Verify all required files exist
echo "Verifying input files..."
for file in "${genotype_bfile}.bed" "${genotype_bfile}.bim" "${genotype_bfile}.fam" \
            "${pheno_file_pop}" \
            "${covar_file_pop}" \
            "${model_snps_file}" \
            "${ld_scores_file}" \
            "${genetic_map_file}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file" >&2
        if [[ "$file" == *".${keep_set}.tsv.gz" ]]; then
            echo "Create population-filtered files by running: sbatch 0a_filter_populations.sbatch.sh" >&2
        fi
        exit 1
    fi
done
echo "✓ All input files verified"
echo ""

# Run BOLT-LMM
echo "Starting BOLT-LMM analysis..."
echo "This will analyze ~1.3M autosomal variants (full genome)"
echo ""

    # BOLT-LMM command
    # Using population-filtered phenotype and covariate files
    # No need for --remove since files only contain desired population samples
    
    bolt \
        --bfile=${genotype_bfile} \
        --phenoFile=${pheno_file_pop} \
        --phenoCol=${phenotype} \
        --covarFile=${covar_file_pop} \
        ${qcovar_col_args} \
        ${covar_col_args} \
        --covarMaxLevels=30 \
        --modelSnps=${model_snps_file} \
        --LDscoresFile=${ld_scores_file} \
        --geneticMapFile=${genetic_map_file} \
        --lmm \
        --LDscoresMatchBp \
        --numThreads=32 \
        --statsFile=${out_file}.stats \
        --verboseStats \
        2>&1 | tee ${out_file}.log

bolt_exit_code=$?

echo ""
echo "BOLT-LMM exit code: ${bolt_exit_code}"

if [ ${bolt_exit_code} -ne 0 ]; then
    echo "ERROR: BOLT-LMM failed" >&2
    exit 1
fi

# Compress output files
echo "Compressing output files..."
if [ -s "${out_file}.stats" ]; then
    gzip -f ${out_file}.stats
    echo "✓ Created: ${out_file}.stats.gz"
fi

if [ -s "${out_file}.log" ]; then
    gzip -f ${out_file}.log
    echo "✓ Created: ${out_file}.log.gz"
fi

# Report success
echo ""
echo "========================================"
echo "✅ COMPLETED: ${phenotype} with ${covar_str} for ${keep_set}"
echo "Output files:"
ls -lh ${out_file}.stats.gz
ls -lh ${out_file}.log.gz
echo "========================================"

