#!/usr/bin/env Rscript
# Aggregate per-tumour CHORD HRD predictions into a cohort CSV.
# Input:  <src_dir>/chord/<tumor_lib>.chord.prediction.tsv
# Output: data/wgs/oncoanalyser_aggregate/chord_hrd.csv
#         (one row per tumour)
#
# Run from project root:
#   Rscript modules/oncoanalyser_aggregate/scripts/agg_chord.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

KEEP <- c("p_BRCA1", "p_BRCA2", "p_hrd", "hr_status", "hrd_type",
          "remarks_hr_status", "remarks_hrd_type")

read_chord <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "chord", paste0(tumor_lib, ".chord.prediction.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())
  df[, intersect(KEEP, names(df)), drop = FALSE]
}

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)), read_chord(src_dir, tumor_lib))

write_cohort(rows, "chord_hrd.csv")
