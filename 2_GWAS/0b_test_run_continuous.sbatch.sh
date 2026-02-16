#!/bin/bash
#SBATCH --job-name=bolt_loneliness_cont_test
#SBATCH --partition=kellis
#SBATCH --mem=100G
#SBATCH -n 32
#SBATCH --time=47:00:00
#SBATCH --output=2_GWAS/logs/0b_test_continuous.out
#SBATCH --error=2_GWAS/logs/0b_test_continuous.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -beEo pipefail

# Test run for continuous-coded phenotypes: One phenotype, one population, full genome

echo "========================================"
echo "BOLT-LMM Social Isolation Test Run (Continuous)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Resources: 100GB RAM, 32 CPUs"
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

echo ""
echo "Testing with Loneliness phenotype (continuous), Day_NoPCs covariate set, EUR_MM population"
echo "This tests the full pipeline on the complete genome (~1.3M variants)"
echo ""

# Clean up any previous test outputs
echo "Removing any previous continuous test outputs..."
rm -f ${SCRIPTDIR}/results_continuous/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.stats*
rm -f ${SCRIPTDIR}/results_continuous/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log*
echo "Ready for clean test run"
echo ""

# Run test with continuous phenotype prefix and results directory
bash ${SCRIPTDIR}/run_single_phenotype.sh Loneliness Day_NoPCs EUR_MM isolation_run_continuous results_continuous

test_exit=$?

echo ""
echo "========================================"
if [ ${test_exit} -eq 0 ]; then
    echo "TEST PASSED! (Continuous)"
    echo ""
    echo "Verification:"
    ls -lh ${SCRIPTDIR}/results_continuous/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.stats.gz 2>/dev/null || echo "Stats file not found"
    ls -lh ${SCRIPTDIR}/results_continuous/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log.gz 2>/dev/null || echo "Log file not found"
    echo ""
    echo "Next steps:"
    echo "1. Review the output files and log"
    echo "2. Check for any warnings or issues"
    echo "3. If everything looks good, submit full continuous analysis:"
    echo "   sbatch 1b_run_bolt_lmm_continuous.sbatch.sh"
    echo ""
else
    echo "TEST FAILED (Continuous)"
    echo "Check error messages above and in 0b_test_continuous.err"
    echo "Do NOT proceed to full analysis"
    exit 1
fi

echo "End time: $(date)"
echo "========================================"
