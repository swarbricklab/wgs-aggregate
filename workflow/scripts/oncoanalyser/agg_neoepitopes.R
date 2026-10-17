#!/usr/bin/env Rscript
# Aggregate per-tumour LINX fusion-derived neoepitope candidates into a cohort CSV.
# Input:  <src_dir>/linx/somatic_annotations/<tumor_lib>.linx.neoepitope.tsv
# Output: data/wgs/oncoanalyser_aggregate/linx_neoepitopes.csv
#         (one row per tumour x candidate fusion neoepitope; tumours with no
#          candidate contribute no rows).
#
# LINX emits this table when run with `write_neo_epitopes`: every plausible
# fusion junction (5' "Up" breakend x 3' "Down" breakend) that could yield a
# neoepitope, with the participating transcripts. It is a *candidate* input to
# the HMF Neo pipeline (NeoEpitopeAnnotator -> NeoScorer), NOT a scored result:
# Neo was not run in this DNA-only track (no *.neo.* outputs, no ISOFOX RNA), so
# these rows are the raw fusion neoepitope candidates only. Scored neoepitopes
# (peptides, expression, HLA context) from the DNA+RNA track are in
# oncoanalyser_dna_rna_aggregate/neo_neoepitopes.csv.
#
# Canonical source: oncoanalyser. sash re-runs LINX on its own PURPLE fit and
# produces a near-identical neoepitope table (same LINX 2.0, ~95% shared
# junctions, sub-0.01 copy-number drift), but all LINX-derived cohort tables
# live here in oncoanalyser_aggregate (see agg_fusions.R / agg_linx_drivers.R),
# and the Neo pipeline — if run — consumes the oncoanalyser LINX output.
#
# Run from project root:
#   Rscript modules/oncoanalyser_aggregate/scripts/agg_neoepitopes.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

# Full neoepitope schema, source order. Kept verbatim (no thresholding) — this
# is a candidate dump; downstream filters by copy number / transcript as needed.
KEEP <- c("geneIdUp", "geneNameUp", "chromosomeUp", "positionUp",
          "orientationUp", "svIdUp",
          "geneIdDown", "geneNameDown", "chromosomeDown", "positionDown",
          "orientationDown", "svIdDown",
          "junctionCopyNumber", "copyNumber", "insertSeq", "chainLength",
          "transcriptsUp", "transcriptsDown")

read_neoepitopes <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "linx", "somatic_annotations",
            paste0(tumor_lib, ".linx.neoepitope.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())
  df[, intersect(KEEP, names(df)), drop = FALSE]
}

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)), read_neoepitopes(src_dir, tumor_lib))

write_cohort(rows, "linx_neoepitopes.csv")
