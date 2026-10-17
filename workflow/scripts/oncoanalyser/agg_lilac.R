#!/usr/bin/env Rscript
# Aggregate per-tumour LILAC HLA typing into a cohort CSV.
# Input:  <src_dir>/lilac/<tumor_lib>.lilac.tsv
# Output: data/wgs/oncoanalyser_aggregate/lilac_alleles.csv
#         (one row per tumour x HLA allele)
#
# Run from project root:
#   Rscript modules/oncoanalyser_aggregate/scripts/agg_lilac.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

COL_MAP <- c(
  Allele                       = "allele",
  TumorCopyNumber              = "tumor_copy_number",
  SomaticMissense              = "somatic_missense",
  SomaticNonsenseOrFrameshift  = "somatic_nonsense_or_frameshift",
  SomaticSplice                = "somatic_splice",
  SomaticInframeIndel          = "somatic_inframe_indel",
  SomaticSynonymous            = "somatic_synonymous"
)

read_lilac <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "lilac", paste0(tumor_lib, ".lilac.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())
  present <- intersect(names(COL_MAP), names(df))
  df <- df[, present, drop = FALSE]
  names(df) <- unname(COL_MAP[present])
  df
}

ss <- load_samplesheet()

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)), read_lilac(src_dir, tumor_lib))

write_cohort(rows, "lilac_alleles.csv")
