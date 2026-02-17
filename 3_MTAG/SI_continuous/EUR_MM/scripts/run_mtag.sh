#!/usr/bin/env bash

module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate /home/mabdel03/data/conda_envs/mtag

MTAG_DIR="/home/mabdel03/data/software/mtag"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPODIR="$(cd "${SCRIPTDIR}/../../../.." && pwd)"
SUMSTATS_DIR="${REPODIR}/2_GWAS/mtag_results_continuous/EUR_MM"
RESULTS_DIR="${SCRIPTDIR}/../results"
mkdir -p ${RESULTS_DIR}

# SI continuous EUR_MM - Social Isolation traits
python ${MTAG_DIR}/mtag.py \
        --sumstats ${SUMSTATS_DIR}/AbilityToConfide.Day_NoPCs.mtag.sumstats.txt,${SUMSTATS_DIR}/FreqSoc.Day_NoPCs.mtag.sumstats.txt,${SUMSTATS_DIR}/Loneliness.Day_NoPCs.mtag.sumstats.txt \
        --out ${RESULTS_DIR}/SI_EUR_MM_Output_continuous \
        --ld_ref_panel ${MTAG_DIR}/ld_ref_panel/eur_w_ld_chr/ \
        --snp_name snpid \
        --a1_name a1 \
        --a2_name a2 \
        --z_name z \
        --p_name pval \
        --n_name n \
        --chr_name chr \
        --bpos_name bpos \
        --eaf_name freq \
        --n_min 0.0 \
        --force \
        --stream_stdout
