#!/usr/bin/env Rscript
# Aggregate per-tumour LINX reported gene fusions into a cohort CSV.
# Input:  <src_dir>/linx/somatic_annotations/<tumor_lib>.linx.fusion.tsv
# Output: data/wgs/oncoanalyser_aggregate/linx_fusions.csv
#         (one row per tumour x REPORTED fusion; tumours with no reported fusion
#          contribute no rows). The full candidate-fusion table — every plausible
#          fusion, not just reported ones — stays in the per-tumour linx output.
#
# Run from project root:
#   Rscript modules/oncoanalyser_aggregate/scripts/agg_fusions.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

KEEP <- c("name", "reportedType", "reportableReasons", "phased", "likelihood",
          "geneStart", "geneEnd", "junctionCopyNumber", "chainLength")

read_fusions <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "linx", "somatic_annotations",
            paste0(tumor_lib, ".linx.fusion.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0 || !"reported" %in% names(df)) return(tibble())
  df <- df %>% filter(reported == "true")
  if (nrow(df) == 0) return(tibble())
  df[, intersect(KEEP, names(df)), drop = FALSE]
}

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)), read_fusions(src_dir, tumor_lib))

write_cohort(rows, "linx_fusions.csv")
