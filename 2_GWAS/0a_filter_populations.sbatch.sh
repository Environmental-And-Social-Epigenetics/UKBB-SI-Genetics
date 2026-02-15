#!/bin/bash
#SBATCH --job-name=filter_populations
#SBATCH --partition=kellis
#SBATCH --mem=16G
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --output=2_GWAS/logs/0a_filter.out
#SBATCH --error=2_GWAS/logs/0a_filter.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -beEo pipefail

echo "========================================"
echo "Filter Phenotype & Covariate Files to Populations"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Start time: $(date)"
echo "========================================"
echo ""

# Activate conda environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate /home/mabdel03/data/conda_envs/Python_Analysis

# Script directory (git repo) and data directory (working area)
SCRIPTDIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/UKBB-SI-Genetics/2_GWAS"
SRCDIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness"
cd ${SRCDIR}

echo "Python version:"
python3 --version
echo ""

# ---- Binary phenotype filtering ----
echo "========================================"
echo "Filtering BINARY phenotype files"
echo "========================================"
echo ""

python3 ${SCRIPTDIR}/filter_populations.py --pheno-prefix isolation_run_binary

binary_exit=$?

if [ ${binary_exit} -ne 0 ]; then
    echo "ERROR: Binary phenotype filtering failed (exit code ${binary_exit})" >&2
    exit ${binary_exit}
fi

echo ""

# ---- Continuous phenotype filtering ----
echo "========================================"
echo "Filtering CONTINUOUS phenotype files"
echo "========================================"
echo ""

python3 ${SCRIPTDIR}/filter_populations.py --pheno-prefix isolation_run_continuous

continuous_exit=$?

if [ ${continuous_exit} -ne 0 ]; then
    echo "ERROR: Continuous phenotype filtering failed (exit code ${continuous_exit})" >&2
    exit ${continuous_exit}
fi

echo ""
echo "========================================"
echo "ALL FILTERING COMPLETED SUCCESSFULLY"
echo "========================================"
echo ""
echo "Created population-filtered files:"
echo ""
echo "  Binary:"
for pop in EUR_MM EUR_Male EUR_Female; do
    f="isolation_run_binary.${pop}.tsv.gz"
    if [ -f "$f" ]; then
        echo "    $f ($(du -h "$f" | cut -f1))"
    else
        echo "    MISSING: $f"
    fi
done
echo ""
echo "  Continuous:"
for pop in EUR_MM EUR_Male EUR_Female; do
    f="isolation_run_continuous.${pop}.tsv.gz"
    if [ -f "$f" ]; then
        echo "    $f ($(du -h "$f" | cut -f1))"
    else
        echo "    MISSING: $f"
    fi
done
echo ""
echo "  Covariates:"
for pop in EUR_MM EUR_Male EUR_Female; do
    f="sqc.${pop}.tsv.gz"
    if [ -f "$f" ]; then
        echo "    $f ($(du -h "$f" | cut -f1))"
    else
        echo "    MISSING: $f"
    fi
done
echo ""
echo "Next steps:"
echo "1. Test binary:      sbatch 0b_test_run.sbatch.sh"
echo "2. Test continuous:   sbatch 0b_test_run_continuous.sbatch.sh"
echo "3. Full binary GWAS:  sbatch 1_run_bolt_lmm.sbatch.sh"
echo "4. Full continuous:   sbatch 1b_run_bolt_lmm_continuous.sbatch.sh"
echo ""
echo "End time: $(date)"
echo "========================================"
