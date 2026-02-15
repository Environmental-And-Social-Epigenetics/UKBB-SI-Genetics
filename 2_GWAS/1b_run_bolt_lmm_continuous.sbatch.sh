#!/bin/bash
#SBATCH --job-name=bolt_loneliness_cont
#SBATCH --partition=kellis
#SBATCH --mem=100G
#SBATCH -n 32
#SBATCH --time=47:00:00
#SBATCH --output=2_GWAS/logs/1b_%a.out
#SBATCH --error=2_GWAS/logs/1b_%a.err
#SBATCH --array=1-9
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

set -beEo pipefail

# Continuous-coded BOLT-LMM GWAS: 9 jobs total
# 3 social isolation phenotypes × 3 population stratifications × 1 covariate set = 9 jobs
# No variant splitting - each job processes the full genome

echo "========================================"
echo "BOLT-LMM Social Isolation GWAS Analysis (Continuous)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Start time: $(date)"
echo "========================================"

# Activate conda environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate /home/mabdel03/data/conda_envs/bolt_lmm

# Script directory (git repo) and data directory (working area)
SCRIPTDIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/UKBB-SI-Genetics/2_GWAS"
SRCDIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness"
cd ${SRCDIR}

# Define phenotypes and population stratifications
phenotypes=(Loneliness FreqSoc AbilityToConfide)  # 3 continuous-coded social isolation phenotypes
keep_sets=(EUR_MM EUR_Male EUR_Female)  # 3 population stratifications
covar_str="Day_NoPCs"  # Single covariate set

# Map array task ID to phenotype and population combination
# Task 1-3: Loneliness with EUR_MM, EUR_Male, EUR_Female
# Task 4-6: FreqSoc with EUR_MM, EUR_Male, EUR_Female
# Task 7-9: AbilityToConfide with EUR_MM, EUR_Male, EUR_Female

n_keeps=${#keep_sets[@]}
pheno_idx=$(( (SLURM_ARRAY_TASK_ID - 1) / n_keeps ))
keep_idx=$(( (SLURM_ARRAY_TASK_ID - 1) % n_keeps ))

phenotype=${phenotypes[$pheno_idx]}
keep_set=${keep_sets[$keep_idx]}

echo "Processing:"
echo "  Phenotype: ${phenotype}"
echo "  Population: ${keep_set}"
echo "  Covariate set: ${covar_str}"
echo "  Phenotype prefix: isolation_run_continuous"
echo "  Results subdir: results_continuous"
echo ""

# Run BOLT-LMM using continuous phenotype files and results directory
bash ${SCRIPTDIR}/run_single_phenotype.sh \
    ${phenotype} \
    ${covar_str} \
    ${keep_set} \
    isolation_run_continuous \
    results_continuous

# Check if successful
exit_code=$?

echo ""
echo "========================================"
if [ ${exit_code} -eq 0 ]; then
    echo "✅ SUCCESS: ${phenotype} with ${covar_str} for ${keep_set} (continuous)"
else
    echo "❌ FAILED: ${phenotype} with ${covar_str} for ${keep_set} (continuous)"
    echo "Exit code: ${exit_code}"
fi
echo "End time: $(date)"
echo "========================================"

exit ${exit_code}
