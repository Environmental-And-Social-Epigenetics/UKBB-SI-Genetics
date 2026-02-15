#!/bin/bash
#SBATCH --job-name=bolt_loneliness
#SBATCH --partition=kellis
#SBATCH --mem=100G
#SBATCH -n 32
#SBATCH --time=47:00:00
#SBATCH --output=1_%a.out
#SBATCH --error=1_%a.err
#SBATCH --array=1-9
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

set -beEo pipefail

# Simplified BOLT-LMM GWAS: 9 jobs total
# 3 social isolation phenotypes × 3 population stratifications × 1 covariate set = 9 jobs
# No variant splitting - each job processes the full genome

echo "========================================"
echo "BOLT-LMM Social Isolation GWAS Analysis"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Start time: $(date)"
echo "========================================"

# Activate conda environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate /home/mabdel03/data/conda_envs/bolt_lmm

# Navigate to analysis directory
SRCDIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness"
cd ${SRCDIR}

# Define phenotypes and population stratifications
phenotypes=(Loneliness FreqSoc AbilityToConfide)  # 3 binary social isolation phenotypes
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
echo ""

# Run BOLT-LMM using the simplified workflow script
bash ${SRCDIR}/run_single_phenotype.sh ${phenotype} ${covar_str} ${keep_set}

# Check if successful
exit_code=$?

echo ""
echo "========================================"
if [ ${exit_code} -eq 0 ]; then
    echo "✅ SUCCESS: ${phenotype} with ${covar_str} for ${keep_set}"
else
    echo "❌ FAILED: ${phenotype} with ${covar_str} for ${keep_set}"
    echo "Exit code: ${exit_code}"
fi
echo "End time: $(date)"
echo "========================================"

exit ${exit_code}

