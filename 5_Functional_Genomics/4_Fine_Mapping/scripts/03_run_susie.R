#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(susieR)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[length(idx)] == length(args)) stop(sprintf("Missing value for %s", flag))
  args[idx[length(idx)] + 1]
}

sumstats_file <- get_arg("--sumstats")
loci_file <- get_arg("--loci")
ld_dir <- get_arg("--ld-dir")
out_dir <- get_arg("--out-dir")
L_components <- as.integer(get_arg("--L", "10"))
coverage <- as.numeric(get_arg("--coverage", "0.95"))
min_snps <- as.integer(get_arg("--min-snps", "20"))

if (is.null(sumstats_file) || is.null(loci_file) || is.null(ld_dir) || is.null(out_dir)) {
  stop("Usage: 03_run_susie.R --sumstats <file> --loci <file> --ld-dir <dir> --out-dir <dir> [--L 10] [--coverage 0.95]")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

detect_col <- function(cols, candidates) {
  for (c in candidates) {
    if (c %in% cols) return(c)
  }
  return(NA_character_)
}

sumstats <- suppressMessages(fread(sumstats_file, sep = "\t", data.table = FALSE))
snp_col <- detect_col(names(sumstats), c("SNP", "snpid"))
chr_col <- detect_col(names(sumstats), c("CHR", "chr"))
bp_col <- detect_col(names(sumstats), c("BP", "bpos"))
beta_col <- detect_col(names(sumstats), c("mtag_beta", "beta", "BETA"))
se_col <- detect_col(names(sumstats), c("mtag_se", "se", "SE"))
p_col <- detect_col(names(sumstats), c("mtag_pval", "pval", "P"))
n_col <- detect_col(names(sumstats), c("N", "n"))

required <- c(snp_col, chr_col, bp_col, beta_col, se_col, p_col)
if (any(is.na(required))) {
  stop(sprintf("Could not infer required columns from sumstats: %s", paste(names(sumstats), collapse = ", ")))
}

sumstats <- sumstats %>%
  transmute(
    SNP = .data[[snp_col]],
    CHR = as.integer(.data[[chr_col]]),
    BP = as.integer(.data[[bp_col]]),
    BETA = as.numeric(.data[[beta_col]]),
    SE = as.numeric(.data[[se_col]]),
    P = as.numeric(.data[[p_col]]),
    N = if (!is.na(n_col)) as.numeric(.data[[n_col]]) else NA_real_
  ) %>%
  filter(!is.na(SNP), !is.na(CHR), !is.na(BP), !is.na(BETA), !is.na(SE), !is.na(P), SE > 0)

loci <- suppressMessages(fread(loci_file, sep = "\t", data.table = FALSE))
if (!all(c("locus_id", "CHR", "START", "END") %in% names(loci))) {
  stop("Loci file must contain columns: locus_id, CHR, START, END")
}

all_pip <- list()
all_cs <- list()
all_meta <- list()

for (i in seq_len(nrow(loci))) {
  locus <- loci[i, ]
  locus_id <- as.character(locus$locus_id)
  chr <- as.integer(locus$CHR)
  start <- as.integer(locus$START)
  end <- as.integer(locus$END)

  ld_file <- file.path(ld_dir, paste0(locus_id, ".ld.gz"))
  snplist_file <- file.path(ld_dir, paste0(locus_id, ".snplist"))

  if (!file.exists(ld_file) || !file.exists(snplist_file)) {
    message(sprintf("[skip] %s: missing LD files.", locus_id))
    next
  }

  locus_sumstats <- sumstats %>%
    filter(CHR == chr, BP >= start, BP <= end) %>%
    distinct(SNP, .keep_all = TRUE)

  if (nrow(locus_sumstats) < min_snps) {
    message(sprintf("[skip] %s: only %d SNPs in locus summary stats.", locus_id, nrow(locus_sumstats)))
    next
  }

  snplist <- suppressMessages(fread(snplist_file, header = FALSE, data.table = FALSE))
  colnames(snplist) <- c("SNP")

  locus_sumstats <- snplist %>%
    inner_join(locus_sumstats, by = "SNP")

  if (nrow(locus_sumstats) < min_snps) {
    message(sprintf("[skip] %s: %d overlapping SNPs after LD/SNP intersection.", locus_id, nrow(locus_sumstats)))
    next
  }

  R <- as.matrix(suppressMessages(fread(ld_file, header = FALSE, data.table = FALSE)))
  if (nrow(R) != ncol(R)) {
    message(sprintf("[skip] %s: LD matrix is not square.", locus_id))
    next
  }

  if (nrow(R) != nrow(snplist)) {
    # Align to the minimum common length if PLINK output and SNP list disagree.
    m <- min(nrow(R), nrow(snplist), nrow(locus_sumstats))
    R <- R[seq_len(m), seq_len(m), drop = FALSE]
    locus_sumstats <- locus_sumstats[seq_len(m), , drop = FALSE]
  } else {
    # Ensure same row count as LD dimension.
    m <- min(nrow(R), nrow(locus_sumstats))
    R <- R[seq_len(m), seq_len(m), drop = FALSE]
    locus_sumstats <- locus_sumstats[seq_len(m), , drop = FALSE]
  }

  if (nrow(locus_sumstats) < min_snps) {
    message(sprintf("[skip] %s: insufficient SNPs after matrix alignment.", locus_id))
    next
  }

  z <- locus_sumstats$BETA / locus_sumstats$SE
  n_eff <- suppressWarnings(median(locus_sumstats$N, na.rm = TRUE))
  if (!is.finite(n_eff) || n_eff <= 0) n_eff <- 100000

  fit <- tryCatch(
    susie_rss(
      z = z,
      R = R,
      n = n_eff,
      L = L_components,
      coverage = coverage,
      estimate_residual_variance = TRUE
    ),
    error = function(e) {
      message(sprintf("[skip] %s: SuSiE failed (%s)", locus_id, e$message))
      NULL
    }
  )
  if (is.null(fit)) next

  pip_df <- locus_sumstats %>%
    mutate(
      locus_id = locus_id,
      Z = z,
      PIP = fit$pip
    ) %>%
    select(locus_id, SNP, CHR, BP, BETA, SE, Z, P, N, PIP)
  all_pip[[length(all_pip) + 1]] <- pip_df

  cs_list <- fit$sets$cs
  if (!is.null(cs_list) && length(cs_list) > 0) {
    cs_rows <- lapply(seq_along(cs_list), function(k) {
      idxs <- cs_list[[k]]
      snps <- locus_sumstats$SNP[idxs]
      data.frame(
        locus_id = locus_id,
        cs_id = paste0(locus_id, "_CS", k),
        n_snps = length(snps),
        snps = paste(snps, collapse = ","),
        min_p = min(locus_sumstats$P[idxs], na.rm = TRUE),
        max_pip = max(fit$pip[idxs], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
    all_cs[[length(all_cs) + 1]] <- bind_rows(cs_rows)
  }

  all_meta[[length(all_meta) + 1]] <- data.frame(
    locus_id = locus_id,
    n_snps = nrow(locus_sumstats),
    n_cs = ifelse(is.null(cs_list), 0, length(cs_list)),
    converged = fit$converged,
    elbo = tail(fit$elbo, 1),
    stringsAsFactors = FALSE
  )
}

pip_out <- file.path(out_dir, "susie_pip.tsv")
cs_out <- file.path(out_dir, "susie_credible_sets.tsv")
meta_out <- file.path(out_dir, "susie_locus_summary.tsv")

if (length(all_pip) > 0) {
  write_tsv(bind_rows(all_pip), pip_out)
} else {
  write_tsv(data.frame(), pip_out)
}

if (length(all_cs) > 0) {
  write_tsv(bind_rows(all_cs), cs_out)
} else {
  write_tsv(data.frame(), cs_out)
}

if (length(all_meta) > 0) {
  write_tsv(bind_rows(all_meta), meta_out)
} else {
  write_tsv(data.frame(), meta_out)
}

message("SuSiE fine-mapping complete.")
message("Wrote: ", pip_out)
message("Wrote: ", cs_out)
message("Wrote: ", meta_out)
