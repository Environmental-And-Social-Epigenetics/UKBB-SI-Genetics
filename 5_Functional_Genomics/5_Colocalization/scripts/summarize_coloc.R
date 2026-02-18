#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[length(idx)] == length(args)) stop(sprintf("Missing value for %s", flag))
  args[idx[length(idx)] + 1]
}

script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
fg_root <- normalizePath(file.path(script_dir, "../.."))

coloc_results <- get_arg("--coloc-results", file.path(fg_root, "5_Colocalization/output/coloc_results.tsv"))
spredixcan_prioritized <- get_arg("--spredixcan-prioritized", file.path(fg_root, "3_S-PrediXcan/output/prioritized_genes.tsv"))
magma_dir <- get_arg("--magma-dir", file.path(fg_root, "2_MAGMA/output/02_gene_analysis"))
out_file <- get_arg("--out", file.path(fg_root, "5_Colocalization/output/summary/coloc_prioritized.tsv"))
pp4_threshold <- as.numeric(get_arg("--pp4-threshold", "0.8"))

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(coloc_results)) {
  stop("coloc results file not found: ", coloc_results)
}

coloc <- suppressMessages(fread(coloc_results, sep = "\t", data.table = FALSE))
if (nrow(coloc) == 0) {
  write_tsv(data.frame(), out_file)
  message("No coloc rows found. Wrote empty file: ", out_file)
  quit(save = "no", status = 0)
}

required <- c("gene", "locus_id", "tissue", "PP.H4")
if (!all(required %in% names(coloc))) {
  stop("coloc results must contain columns: ", paste(required, collapse = ", "))
}

coloc_summary <- coloc %>%
  mutate(
    gene = toupper(as.character(gene)),
    PP.H4 = as.numeric(PP.H4)
  ) %>%
  filter(!is.na(gene), !is.na(PP.H4)) %>%
  group_by(gene) %>%
  summarise(
    coloc_max_pp_h4 = max(PP.H4, na.rm = TRUE),
    coloc_n_loci = n_distinct(locus_id),
    coloc_n_strong_loci = sum(PP.H4 >= pp4_threshold, na.rm = TRUE),
    best_locus = locus_id[which.max(PP.H4)],
    best_tissue = tissue[which.max(PP.H4)],
    .groups = "drop"
  )

if (file.exists(spredixcan_prioritized)) {
  sp <- suppressMessages(fread(spredixcan_prioritized, sep = "\t", data.table = FALSE))
  if ("gene" %in% names(sp)) {
    sp <- sp %>%
      mutate(gene = toupper(as.character(gene))) %>%
      select(any_of(c("gene", "prioritization_score", "convergent_hits", "magma_p", "spredixcan_min_p", "smultixcan_p", "coloc_pp_h4")))
    coloc_summary <- left_join(coloc_summary, sp, by = "gene")
  }
}

if (dir.exists(magma_dir)) {
  magma_files <- list.files(magma_dir, pattern = "\\.genes\\.out$", full.names = TRUE)
  magma_rows <- list()
  for (f in magma_files) {
    df <- tryCatch(
      suppressMessages(fread(f, data.table = FALSE)),
      error = function(e) NULL
    )
    if (is.null(df)) next
    if (!all(c("GENE", "P") %in% names(df))) next
    trait <- sub("\\.genes\\.out$", "", basename(f))
    magma_rows[[length(magma_rows) + 1]] <- df %>%
      transmute(gene = toupper(as.character(GENE)), magma_gene_p = as.numeric(P), magma_trait = trait)
  }
  if (length(magma_rows) > 0) {
    magma_df <- bind_rows(magma_rows) %>%
      filter(!is.na(gene), !is.na(magma_gene_p)) %>%
      arrange(magma_gene_p) %>%
      distinct(gene, .keep_all = TRUE)
    coloc_summary <- left_join(coloc_summary, magma_df, by = "gene")
  }
}

coloc_summary <- coloc_summary %>%
  mutate(
    strong_coloc = coloc_max_pp_h4 >= pp4_threshold,
    coloc_rank_score = coloc_max_pp_h4 * 5 + log1p(coloc_n_strong_loci) + log1p(coloc_n_loci)
  ) %>%
  arrange(desc(strong_coloc), desc(coloc_rank_score), desc(coloc_max_pp_h4))

write_tsv(coloc_summary, out_file)
message("Wrote coloc-prioritized summary: ", out_file)
