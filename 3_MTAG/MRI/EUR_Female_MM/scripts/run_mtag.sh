#!/usr/bin/env bash

module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate /home/mabdel03/data/conda_envs/mtag

MTAG_DIR="/home/mabdel03/data/software/mtag"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUMSTATS_DIR="${SCRIPTDIR}/../sumstats"
RESULTS_DIR="${SCRIPTDIR}/../results"
mkdir -p ${RESULTS_DIR}

# MRI EUR_Female_MM - MRI traits (Female)
python ${MTAG_DIR}/mtag.py \
        --sumstats ${SUMSTATS_DIR}/FA.Day_NoPCs.mtag.sumstats.txt,${SUMSTATS_DIR}/MD.Day_NoPCs.mtag.sumstats.txt,${SUMSTATS_DIR}/MO.Day_NoPCs.mtag.sumstats.txt,${SUMSTATS_DIR}/OD.Day_NoPCs.mtag.sumstats.txt \
        --out ${RESULTS_DIR}/MRI_EUR_Female_MM_Output \
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
