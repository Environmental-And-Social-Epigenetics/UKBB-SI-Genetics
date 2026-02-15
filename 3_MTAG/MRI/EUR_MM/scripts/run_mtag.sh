#!/usr/bin/env bash

source /Users/mahmoudabdelmoneum/opt/anaconda3/etc/profile.d/conda.sh
conda activate /Users/mahmoudabdelmoneum/Desktop/MIT/Software/Research_Software/conda_envs/MTAG

# MRI EUR_MM - MRI traits
python /Users/mahmoudabdelmoneum/Desktop/MIT/Software/Research_Software/mtag/mtag.py \
        --sumstats ../sumstats/FA.Day_NoPCs.mtag.sumstats.txt,../sumstats/MD.Day_NoPCs.mtag.sumstats.txt,../sumstats/MO.Day_NoPCs.mtag.sumstats.txt,../sumstats/OD.Day_NoPCs.mtag.sumstats.txt \
        --out ../results/MRI_EUR_MM_Output \
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


