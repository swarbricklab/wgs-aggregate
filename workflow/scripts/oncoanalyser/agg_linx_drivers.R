#!/usr/bin/env Rscript
# Aggregate per-tumour LINX structural-variant drivers into a cohort CSV.
# Input:  <src_dir>/linx/somatic_annotations/<tumor_lib>.linx.driver.catalog.tsv
# Output: data/wgs/oncoanalyser_aggregate/linx_drivers.csv
#         (one row per tumour x SV-mediated driver)
#
# Scope: SV-mediated gene-disruption drivers only — the events PURPLE's point /
# copy-number driver catalog does NOT capture. The point-mutation, AMP and DEL
# drivers (which LINX imports verbatim from PURPLE) are already aggregated by
# sash_aggregate (purple_drivers_*), so we deliberately exclude them here to
# avoid duplication. Fusions-as-drivers are covered by agg_fusions.R.
#
# Run from project root:
#   Rscript modules/oncoanalyser_aggregate/scripts/agg_linx_drivers.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

# LINX driver classes that arise from structural variants (not in PURPLE's
# point/CN driver catalog).
SV_DRIVERS <- c("DISRUPTION", "HOM_DUP_DISRUPTION", "HOM_DEL_DISRUPTION",
                "GERMLINE_DISRUPTION")

KEEP <- c("gene", "driver", "category", "likelihoodMethod", "driverLikelihood",
          "biallelic", "minCopyNumber", "maxCopyNumber",
          "chromosome", "chromosomeBand")

read_drivers <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "linx", "somatic_annotations",
            paste0(tumor_lib, ".linx.driver.catalog.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0 || !"driver" %in% names(df)) return(tibble())
  df <- df %>% filter(driver %in% SV_DRIVERS)
  if (nrow(df) == 0) return(tibble())
  df[, intersect(KEEP, names(df)), drop = FALSE]
}

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)), read_drivers(src_dir, tumor_lib))

write_cohort(rows, "linx_drivers.csv")
