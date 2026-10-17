#!/usr/bin/env Rscript
# Aggregate per-tumour CUPPA tissue-of-origin predictions into a cohort CSV.
# Input:  <src_dir>/cuppa/<tumor_lib>.cuppa.pred_summ.tsv
# Output: data/wgs/oncoanalyser_aggregate/cuppa_predictions.csv
#         (one row per tumour: the combined DNA classifier's top-3 calls)
#
# Run from project root:
#   Rscript modules/oncoanalyser_aggregate/scripts/agg_cuppa.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

KEEP <- c("clf_group", "clf_name",
          "pred_class_1", "pred_class_2", "pred_class_3",
          "pred_prob_1", "pred_prob_2", "pred_prob_3")

read_cuppa <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "cuppa", paste0(tumor_lib, ".cuppa.pred_summ.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0 || !"clf_name" %in% names(df)) return(tibble())
  # the combined DNA classifier is the headline tissue-of-origin call
  df <- df %>% filter(clf_name == "dna_combined")
  if (nrow(df) == 0) return(tibble())
  df[, intersect(KEEP, names(df)), drop = FALSE]
}

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)), read_cuppa(src_dir, tumor_lib))

write_cohort(rows, "cuppa_predictions.csv")
