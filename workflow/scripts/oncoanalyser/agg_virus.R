#!/usr/bin/env Rscript
# Aggregate per-tumour Virus Interpreter calls into a cohort CSV.
# Input:  <src_dir>/virusinterpreter/<tumor_lib>.virus.annotated.tsv
# Output: data/wgs/oncoanalyser_aggregate/virus_annotated.csv
#         (one row per tumour x detected virus; tumours with no virus contribute
#          no rows)
#
# Run from project root:
#   Rscript modules/oncoanalyser_aggregate/scripts/agg_virus.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

KEEP <- c("taxid", "name", "qcStatus", "integrations", "interpretation",
          "percentageCovered", "meanCoverage", "reported", "driverLikelihood")

read_virus <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "virusinterpreter",
            paste0(tumor_lib, ".virus.annotated.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())
  df[, intersect(KEEP, names(df)), drop = FALSE]
}

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)), read_virus(src_dir, tumor_lib))

write_cohort(rows, "virus_annotated.csv")
