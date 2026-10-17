#!/usr/bin/env Rscript
# Aggregate per-tumour HMF Neo (NeoScorer) scored neoepitopes into a cohort CSV.
# Input:  <src_dir>/neo/scorer/<tumor_lib>.neo.neoepitope.tsv
# Output: data/wgs/oncoanalyser_dna_rna_aggregate/neo_neoepitopes.csv
#         (one row per tumour x neoepitope (NeId); tumours whose scorer file is
#          header-only contribute no rows).
#
# This is the RNA-integrated *scored* summary from the HMF Neo pipeline
# (NeoFinder -> NeoScorer): one row per candidate neoepitope with its
# amino-acid context, RNA-informed expression (TPM), NMD, copy-number and
# subclonal context. It is the scored counterpart to the LINX fusion neoepitope
# *candidates* aggregated by oncoanalyser_aggregate/agg_neoepitopes.R
# (data/wgs/oncoanalyser_aggregate/linx_neoepitopes.csv) — those are the Neo input;
# these are the Neo output.
#
# Kept verbatim (all 32 columns, no thresholding / reshaping), mirroring that
# LINX aggregator. Reading all-character preserves scientific-notation TPM
# values and the -1 NMD sentinel exactly. The finer per-peptide x HLA-allele
# binding table (<tumor_lib>.neo.peptide_scores.tsv) is intentionally NOT
# aggregated here (see README) - enormous, finest grain; PeptideCount already
# summarises it per neoepitope.
#
# Run from project root:
#   Rscript modules/oncoanalyser_dna_rna_aggregate/scripts/agg_neo_neoepitopes.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

# Full NeoScorer neoepitope schema, source order. Kept verbatim (no thresholding).
KEEP <- c("NeId", "VariantType", "VariantInfo", "GeneName",
          "UpAminoAcids", "NovelAminoAcids", "DownAminoAcids", "PeptideCount",
          "TpmSource", "RnaFrags", "RnaDepth", "TpmUp", "TpmDown", "ExpectedTpm",
          "RawEffectiveTpm", "EffectiveTpm", "TpmCancerUp", "TpmCancerDown",
          "TpmPanCancerUp", "TpmPanCancerDown", "NmdMin", "NmdMax",
          "CodingBasesLengthMin", "CodingBasesLengthMax", "FusedIntronLength",
          "SkippedDonors", "SkippedAcceptors", "TranscriptsUp", "TranscriptsDown",
          "VariantCopyNumber", "CopyNumber", "SubclonalLikelihood")

read_neo_neoepitopes <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "neo", "scorer",
            paste0(tumor_lib, ".neo.neoepitope.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())          # header-only scorer file
  df[, intersect(KEEP, names(df)), drop = FALSE]
}

rows <- ss %>%
  rowwise() %>%
  reframe(pick(any_of(IDENTITY_COLS)),
          read_neo_neoepitopes(src_dir, tumor_lib))

write_cohort(rows, "neo_neoepitopes.csv")
