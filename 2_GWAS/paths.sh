# Path configuration for BOLT-LMM Social Isolation/Loneliness analysis

ukb21942_d='/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942'
tmp_ukb21942_d='/tmp/kellis/ukb21942'

# BOLT-LMM specific paths
# BOLT-LMM v2.5 installation
BOLT_LMM_DIR="/home/mabdel03/data/software/BOLT-LMM_v2.5"
BOLT_TABLES_DIR="${BOLT_LMM_DIR}/tables"

# LD scores for BOLT-LMM calibration
# For European ancestry (GRCh38 coordinates, but matched by rsID/position)
LD_SCORES_FILE="${BOLT_TABLES_DIR}/LDSCORE.1000G_EUR.GRCh38.tab.gz"

# Genetic map for interpolation (hg19/GRCh37 - matches UK Biobank coordinates)
GENETIC_MAP_FILE="${BOLT_TABLES_DIR}/genetic_map_hg19_withX.txt.gz"

# Model SNPs for genetic relationship matrix
MODEL_SNPS_FILE="${ukb21942_d}/geno/ukb_genoHM3/ukb_genoHM3_modelSNPs.txt"

