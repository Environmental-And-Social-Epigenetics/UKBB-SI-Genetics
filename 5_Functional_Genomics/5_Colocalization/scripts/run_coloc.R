#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[length(idx)] == length(args)) stop(sprintf("Missing value for %s", flag))
  args[idx[length(idx)] + 1]
}

gwas_file <- get_arg("--gwas")
loci_file <- get_arg("--loci")
eqtl_manifest <- get_arg("--eqtl-manifest")
out_dir <- get_arg("--out-dir")
pp4_threshold <- as.numeric(get_arg("--pp4-threshold", "0.8"))
min_overlap <- as.integer(get_arg("--min-overlap-snps", "50"))
p1 <- as.numeric(get_arg("--p1", "1e-4"))
p2 <- as.numeric(get_arg("--p2", "1e-4"))
p12 <- as.numeric(get_arg("--p12", "1e-5"))

if (is.null(gwas_file) || is.null(loci_file) || is.null(eqtl_manifest) || is.null(out_dir)) {
  stop(
    "Usage: run_coloc.R --gwas <file> --loci <file> --eqtl-manifest <file> --out-dir <dir> ",
    "[--pp4-threshold 0.8 --min-overlap-snps 50 --p1 1e-4 --p2 1e-4 --p12 1e-5]"
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pick_col <- function(cols, candidates) {
  for (c in candidates) {
    if (c %in% cols) return(c)
  }
  return(NA_character_)
}

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

read_table_flexible <- function(path) {
  tryCatch(
    suppressMessages(fread(path, sep = "\t", data.table = FALSE)),
    error = function(e) suppressMessages(fread(path, data.table = FALSE))
  )
}

normalize_gwas <- function(path) {
  df <- read_table_flexible(path)
  snp_col <- pick_col(names(df), c("SNP", "snpid"))
  chr_col <- pick_col(names(df), c("CHR", "chr"))
  bp_col <- pick_col(names(df), c("BP", "bpos"))
  beta_col <- pick_col(names(df), c("mtag_beta", "beta", "BETA"))
  se_col <- pick_col(names(df), c("mtag_se", "se", "SE"))
  p_col <- pick_col(names(df), c("mtag_pval", "pval", "P"))
  n_col <- pick_col(names(df), c("N", "n"))
  maf_col <- pick_col(names(df), c("FRQ", "freq", "maf", "MAF"))

  required <- c(snp_col, chr_col, bp_col, beta_col, se_col, p_col)
  if (any(is.na(required))) {
    stop("GWAS file is missing required columns. Found columns: ", paste(names(df), collapse = ", "))
  }

  out <- df %>%
    transmute(
      SNP = as.character(.data[[snp_col]]),
      CHR = to_numeric_safe(.data[[chr_col]]),
      BP = to_numeric_safe(.data[[bp_col]]),
      BETA = to_numeric_safe(.data[[beta_col]]),
      SE = to_numeric_safe(.data[[se_col]]),
      P = to_numeric_safe(.data[[p_col]]),
      N = if (!is.na(n_col)) to_numeric_safe(.data[[n_col]]) else NA_real_,
      MAF = if (!is.na(maf_col)) to_numeric_safe(.data[[maf_col]]) else NA_real_
    ) %>%
    filter(!is.na(SNP), !is.na(CHR), !is.na(BP), !is.na(BETA), !is.na(SE), !is.na(P), P > 0, P <= 1, SE > 0) %>%
    distinct(SNP, .keep_all = TRUE)

  out
}

normalize_eqtl <- function(df, config_row) {
  get_cfg <- function(name, default_val) {
    v <- as.character(config_row[[name]])
    if (length(v) == 0 || is.na(v) || v == "") default_val else v
  }

  snp_col <- get_cfg("snp_col", "SNP")
  chr_col <- get_cfg("chr_col", "CHR")
  bp_col <- get_cfg("bp_col", "BP")
  beta_col <- get_cfg("beta_col", "beta")
  se_col <- get_cfg("se_col", "se")
  p_col <- get_cfg("p_col", "pval")
  n_col <- get_cfg("n_col", "n")
  maf_col <- get_cfg("maf_col", "maf")
  gene_col <- get_cfg("gene_col", "gene")
  tissue_col <- get_cfg("tissue_col", "tissue")

  required <- c(snp_col, chr_col, bp_col, beta_col, se_col, p_col)
  miss <- required[!required %in% names(df)]
  if (length(miss) > 0) {
    stop(sprintf("Missing required eQTL columns: %s", paste(miss, collapse = ", ")))
  }

  out <- df %>%
    transmute(
      SNP = as.character(.data[[snp_col]]),
      CHR = to_numeric_safe(.data[[chr_col]]),
      BP = to_numeric_safe(.data[[bp_col]]),
      BETA = to_numeric_safe(.data[[beta_col]]),
      SE = to_numeric_safe(.data[[se_col]]),
      P = to_numeric_safe(.data[[p_col]]),
      N = if (n_col %in% names(df)) to_numeric_safe(.data[[n_col]]) else NA_real_,
      MAF = if (maf_col %in% names(df)) to_numeric_safe(.data[[maf_col]]) else NA_real_,
      gene = if (gene_col %in% names(df)) as.character(.data[[gene_col]]) else as.character(config_row$gene),
      tissue = if (tissue_col %in% names(df)) as.character(.data[[tissue_col]]) else as.character(config_row$tissue)
    ) %>%
    filter(!is.na(SNP), !is.na(CHR), !is.na(BP), !is.na(BETA), !is.na(SE), !is.na(P), P > 0, P <= 1, SE > 0) %>%
    distinct(SNP, .keep_all = TRUE)

  out
}

build_coloc_dataset <- function(df, type = "quant", fallback_n = 100000) {
  n_val <- suppressWarnings(median(df$N, na.rm = TRUE))
  if (!is.finite(n_val) || n_val <= 0) n_val <- fallback_n

  ds <- list(
    beta = df$BETA,
    varbeta = df$SE^2,
    snp = df$SNP,
    N = n_val,
    type = type
  )
  if (all(is.finite(df$MAF))) ds$MAF <- df$MAF
  ds
}

gwas <- normalize_gwas(gwas_file)
loci <- read_table_flexible(loci_file)
manifest <- read_table_flexible(eqtl_manifest)

if (!all(c("locus_id", "CHR", "START", "END") %in% names(loci))) {
  stop("Loci file must have columns: locus_id, CHR, START, END")
}
if (!all(c("dataset_id", "eqtl_file", "gene", "tissue") %in% names(manifest))) {
  stop("Manifest must include columns: dataset_id, eqtl_file, gene, tissue")
}

eqtl_cache <- list()
results <- list()

for (i in seq_len(nrow(manifest))) {
  cfg <- manifest[i, ]
  eqtl_path <- as.character(cfg$eqtl_file)
  if (!file.exists(eqtl_path)) {
    message(sprintf("[skip] missing eQTL file: %s", eqtl_path))
    next
  }

  cache_key <- paste0(eqtl_path, "::", cfg$dataset_id)
  if (!(cache_key %in% names(eqtl_cache))) {
    raw_eqtl <- read_table_flexible(eqtl_path)
    eqtl_cache[[cache_key]] <- tryCatch(
      normalize_eqtl(raw_eqtl, cfg),
      error = function(e) {
        message(sprintf("[skip] could not normalize eQTL %s (%s)", eqtl_path, e$message))
        NULL
      }
    )
  }
  eqtl_df <- eqtl_cache[[cache_key]]
  if (is.null(eqtl_df)) next

  for (j in seq_len(nrow(loci))) {
    locus <- loci[j, ]
    locus_id <- as.character(locus$locus_id)
    chr <- as.numeric(locus$CHR)
    start <- as.numeric(locus$START)
    end <- as.numeric(locus$END)

    gwas_locus <- gwas %>%
      filter(CHR == chr, BP >= start, BP <= end)
    if (nrow(gwas_locus) < min_overlap) next

    eqtl_locus <- eqtl_df %>%
      filter(CHR == chr, BP >= start, BP <= end)
    if (nrow(eqtl_locus) < min_overlap) next

    merged <- inner_join(
      gwas_locus %>% select(SNP, BETA_g = BETA, SE_g = SE, P_g = P, N_g = N, MAF_g = MAF),
      eqtl_locus %>% select(SNP, BETA_e = BETA, SE_e = SE, P_e = P, N_e = N, MAF_e = MAF),
      by = "SNP"
    )
    if (nrow(merged) < min_overlap) next

    gwas_ds <- data.frame(
      SNP = merged$SNP,
      BETA = merged$BETA_g,
      SE = merged$SE_g,
      P = merged$P_g,
      N = merged$N_g,
      MAF = merged$MAF_g
    )
    eqtl_ds <- data.frame(
      SNP = merged$SNP,
      BETA = merged$BETA_e,
      SE = merged$SE_e,
      P = merged$P_e,
      N = merged$N_e,
      MAF = merged$MAF_e
    )

    fit <- tryCatch(
      coloc.abf(
        dataset1 = build_coloc_dataset(gwas_ds, type = "quant"),
        dataset2 = build_coloc_dataset(eqtl_ds, type = "quant"),
        p1 = p1,
        p2 = p2,
        p12 = p12
      ),
      error = function(e) {
        message(sprintf("[skip] coloc failed (%s, %s): %s", locus_id, cfg$dataset_id, e$message))
        NULL
      }
    )
    if (is.null(fit)) next

    top_snp <- NA_character_
    top_snp_pp <- NA_real_
    if ("results" %in% names(fit) && "SNP.PP.H4" %in% names(fit$results)) {
      ord <- order(fit$results$SNP.PP.H4, decreasing = TRUE)
      if (length(ord) > 0) {
        top_snp <- fit$results$snp[ord[1]]
        top_snp_pp <- fit$results$SNP.PP.H4[ord[1]]
      }
    }

    summary_row <- fit$summary
    results[[length(results) + 1]] <- data.frame(
      locus_id = locus_id,
      dataset_id = as.character(cfg$dataset_id),
      gene = as.character(cfg$gene),
      tissue = as.character(cfg$tissue),
      eqtl_file = eqtl_path,
      n_overlap_snps = nrow(merged),
      PP.H0 = as.numeric(summary_row["PP.H0.abf"]),
      PP.H1 = as.numeric(summary_row["PP.H1.abf"]),
      PP.H2 = as.numeric(summary_row["PP.H2.abf"]),
      PP.H3 = as.numeric(summary_row["PP.H3.abf"]),
      PP.H4 = as.numeric(summary_row["PP.H4.abf"]),
      top_snp = top_snp,
      top_snp_pp_h4 = top_snp_pp,
      stringsAsFactors = FALSE
    )
  }
}

all_results <- bind_rows(results)
results_path <- file.path(out_dir, "coloc_results.tsv")
strong_path <- file.path(out_dir, "coloc_strong_hits.tsv")

if (nrow(all_results) == 0) {
  write_tsv(data.frame(), results_path)
  write_tsv(data.frame(), strong_path)
  message("No coloc results produced.")
} else {
  all_results <- all_results %>% arrange(desc(PP.H4), desc(n_overlap_snps))
  strong_hits <- all_results %>% filter(PP.H4 >= pp4_threshold)
  write_tsv(all_results, results_path)
  write_tsv(strong_hits, strong_path)
  message("Wrote coloc results: ", results_path)
  message("Wrote strong coloc hits: ", strong_path)
}
