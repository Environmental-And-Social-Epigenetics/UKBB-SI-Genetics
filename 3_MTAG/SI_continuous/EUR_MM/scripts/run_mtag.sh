#!/usr/bin/env bash

source /Users/mahmoudabdelmoneum/opt/anaconda3/etc/profile.d/conda.sh
conda activate /Users/mahmoudabdelmoneum/Desktop/MIT/Software/Research_Software/conda_envs/MTAG

BOLT_DIR="/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness"

# SI continuous EUR_MM - Social Isolation traits
python /Users/mahmoudabdelmoneum/Desktop/MIT/Software/Research_Software/mtag/mtag.py \
        --sumstats ${BOLT_DIR}/mtag_results_continuous/EUR_MM/AbilityToConfide.Day_NoPCs.mtag.sumstats.txt,${BOLT_DIR}/mtag_results_continuous/EUR_MM/FreqSoc.Day_NoPCs.mtag.sumstats.txt,${BOLT_DIR}/mtag_results_continuous/EUR_MM/Loneliness.Day_NoPCs.mtag.sumstats.txt \
        --out ../results/SI_EUR_MM_Output_continuous \
        --ld_ref_panel /Users/mahmoudabdelmoneum/Desktop/MIT/Software/Research_Software/mtag/ld_ref_panel/eur_w_ld_chr/ \
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
