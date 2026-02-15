#!/bin/bash
#SBATCH --job-name=bolt_loneliness_test
#SBATCH --partition=kellis
#SBATCH --mem=100G
#SBATCH -n 32
#SBATCH --time=47:00:00
#SBATCH --output=0b_test.out
#SBATCH --error=0b_test.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -beEo pipefail

# Simplified test run: One phenotype, one population, full genome

echo "========================================"
echo "BOLT-LMM Social Isolation Test Run"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Resources: 100GB RAM, 32 CPUs"
echo "Start time: $(date)"
echo "========================================"

# Activate conda environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate /home/mabdel03/data/conda_envs/bolt_lmm

# Navigate to analysis directory
SRCDIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness"
cd ${SRCDIR}

echo ""
echo "Testing with Loneliness phenotype, Day_NoPCs covariate set, EUR_MM population"
echo "This tests the full pipeline on the complete genome (~1.3M variants)"
echo ""

# Clean up any previous test outputs
echo "Removing any previous test outputs..."
rm -f ${SRCDIR}/results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.stats*
rm -f ${SRCDIR}/results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log*
echo "✓ Ready for clean test run"
echo ""

# Run test
bash run_single_phenotype.sh Loneliness Day_NoPCs EUR_MM

test_exit=$?

echo ""
echo "========================================"
if [ ${test_exit} -eq 0 ]; then
    echo "🎉 TEST PASSED!"
    echo ""
    echo "Verification:"
    ls -lh results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.stats.gz 2>/dev/null || echo "Stats file not found"
    ls -lh results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log.gz 2>/dev/null || echo "Log file not found"
    echo ""
    echo "Next steps:"
    echo "1. Review the output files and log"
    echo "2. Check for any warnings or issues"
    echo "3. If everything looks good, submit full analysis:"
    echo "   sbatch 1_run_bolt_lmm.sbatch.sh"
    echo ""
else
    echo "❌ TEST FAILED"
    echo "Check error messages above and in 0b_test.err"
    echo "Do NOT proceed to full analysis"
    exit 1
fi

echo "End time: $(date)"
echo "========================================"

