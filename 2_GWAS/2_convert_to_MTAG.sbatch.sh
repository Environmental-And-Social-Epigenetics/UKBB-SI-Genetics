#!/bin/bash
#SBATCH --job-name=mtag_convert_loneliness
#SBATCH --partition=kellis
#SBATCH --mem=32G
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --output=2_GWAS/logs/2_mtag_convert.out
#SBATCH --error=2_GWAS/logs/2_mtag_convert.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -beEo pipefail

echo "========================================"
echo "BOLT-LMM to MTAG Format Conversion"
echo "Social Isolation Phenotypes"
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

echo "Pandas version:"
python3 -c "import pandas; print(pandas.__version__)"
echo ""

# Run conversion script
echo "Running MTAG conversion..."
echo ""

python3 ${SCRIPTDIR}/convert_to_MTAG.py

convert_exit=$?

echo ""
echo "========================================"
if [ ${convert_exit} -eq 0 ]; then
    echo "✅ MTAG CONVERSION COMPLETED SUCCESSFULLY"
    echo ""
    echo "Output directory: ${SCRIPTDIR}/mtag_results/"
    echo ""
    echo "Created files by population:"
    if [ -d "${SCRIPTDIR}/mtag_results" ]; then
        for pop in EUR_MM EUR_Male EUR_Female; do
            if [ -d "${SCRIPTDIR}/mtag_results/${pop}" ]; then
                echo ""
                echo "  ${pop}:"
                ls -lh "${SCRIPTDIR}/mtag_results/${pop}/"*.mtag.sumstats.txt 2>/dev/null | \
                    awk '{printf "    %s (%s)\n", $9, $5}'
            fi
        done
    fi
    echo ""
    echo "Next steps:"
    echo "1. Verify output files in mtag_results/<population>/"
    echo "2. Run MTAG analysis on these formatted files"
else
    echo "❌ MTAG CONVERSION FAILED"
    echo "Exit code: ${convert_exit}"
    echo ""
    echo "Check error log: 2_mtag_convert.err"
fi
echo ""
echo "End time: $(date)"
echo "========================================"

exit ${convert_exit}











