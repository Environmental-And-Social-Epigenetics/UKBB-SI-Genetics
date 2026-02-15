#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
})

# Update these paths for your continuous-coded study
exposure_file <- "path/to/3_MTAG/SI_continuous/EUR_MM/results/SI_EUR_MM_Output_continuous_trait_1.txt"
outcome_file  <- "path/to/3_MTAG/MRI/EUR_MM/results/MRI_EUR_MM_Output_trait_1.txt"
out_prefix    <- "mr_results_continuous/loneliness_to_mri"

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

# Example column mapping. Adjust to your file schema.
read_exposure_data_local <- function(path) {
  TwoSampleMR::read_exposure_data(
    filename = path,
    sep = "\t",
    snp_col = "snpid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "a1",
    other_allele_col = "a2",
    eaf_col = "freq",
    pval_col = "pval",
    samplesize_col = "n",
    chr_col = "chr",
    pos_col = "bpos"
  )
}

read_outcome_data_local <- function(path, snps) {
  TwoSampleMR::read_outcome_data(
    filename = path,
    snps = snps,
    sep = "\t",
    snp_col = "snpid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "a1",
    other_allele_col = "a2",
    eaf_col = "freq",
    pval_col = "pval",
    samplesize_col = "n",
    chr_col = "chr",
    pos_col = "bpos"
  )
}

message("Loading exposure data...")
exposure_dat <- read_exposure_data_local(exposure_file)

message("Clumping exposure instruments...")
exposure_dat <- clump_data(exposure_dat)

message("Loading outcome data...")
outcome_dat <- read_outcome_data_local(outcome_file, exposure_dat$SNP)

message("Harmonizing datasets...")
harmonized <- harmonise_data(exposure_dat, outcome_dat)

message("Running MR methods...")
mr_results <- mr(harmonized)
het_results <- mr_heterogeneity(harmonized)
pleio_results <- mr_pleiotropy_test(harmonized)
loo_results <- mr_leaveoneout(harmonized)

write_tsv(mr_results, paste0(out_prefix, "_main.tsv"))
write_tsv(het_results, paste0(out_prefix, "_heterogeneity.tsv"))
write_tsv(pleio_results, paste0(out_prefix, "_pleiotropy.tsv"))
write_tsv(loo_results, paste0(out_prefix, "_leaveoneout.tsv"))

message("Generating basic plots...")
pdf(paste0(out_prefix, "_plots.pdf"), width = 10, height = 8)
print(mr_scatter_plot(mr_results, harmonized)[[1]])
print(mr_forest_plot(loo_results)[[1]])
print(mr_funnel_plot(loo_results)[[1]])
dev.off()

message("Continuous MR scaffold run complete.")
message("Outputs written with prefix: ", out_prefix)
