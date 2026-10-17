#!/usr/bin/env Rscript
# Aggregate per-tumour ORANGE immune-escape flags into a cohort CSV.
# Input:  <src_dir>/orange/<tumor_lib>.orange.json   ($immuneEscape)
# Output: data/wgs/oncoanalyser_aggregate/immune_escape.csv
#         (one row per tumour; six boolean immune-evasion mechanism flags)
#
# These flags are an ORANGE-integrated call (drawing on LILAC, PURPLE drivers,
# expression-independent genomic signals) and are not surfaced by any other
# cohort aggregate.
#
# Run from project root:
#   Rscript modules/oncoanalyser_aggregate/scripts/agg_immune_escape.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))
suppressPackageStartupMessages(library(jsonlite))

ss <- load_samplesheet()

# ORANGE immuneEscape boolean fields -> output column names.
FLAGS <- c(
  hasHlaEscape                       = "hla_escape",
  hasAntigenPresentationPathwayEscape = "antigen_presentation_escape",
  hasIFNGammaPathwayEscape           = "ifn_gamma_escape",
  hasPDL1OverexpressionEscape        = "pdl1_overexpression_escape",
  hasCD58InactivationEscape          = "cd58_inactivation_escape",
  hasEpigeneticSETDB1Escape          = "setdb1_epigenetic_escape"
)

read_immune <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "orange", paste0(tumor_lib, ".orange.json"))
  if (!file.exists(p)) return(tibble())
  # NB: do NOT fromJSON() the whole file. Two problems with the full document:
  #   1. ORANGE writes non-standard `NaN` literals (e.g. "...":NaN) in the
  #      SV/fusion sections, which jsonlite rejects (silently losing those runs).
  #   2. these files are 18-88 MB, so a full parse / PCRE scan is far too slow.
  # Instead: read once, locate the flat (brace-free) immuneEscape object with a
  # cheap fixed-string search, and parse only that small substring (~0.3s/file).
  txt <- readChar(p, file.info(p)$size, useBytes = TRUE)
  k <- regexpr('"immuneEscape":', txt, fixed = TRUE)
  if (k < 0) return(tibble())
  start <- k + attr(k, "match.length")          # position of the opening '{'
  sub <- substring(txt, start, start + 500)      # object is < 500 chars
  endrel <- regexpr("}", sub, fixed = TRUE)
  if (endrel < 0) return(tibble())
  ie <- tryCatch(fromJSON(substring(sub, 1, endrel)), error = function(e) NULL)
  if (is.null(ie)) return(tibble())
  vals <- lapply(names(FLAGS), function(k) {
    v <- ie[[k]]
    if (is.null(v)) NA else as.logical(v)
  })
  names(vals) <- unname(FLAGS)
  as_tibble(vals)
}

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)), read_immune(src_dir, tumor_lib))

write_cohort(rows, "immune_escape.csv")
