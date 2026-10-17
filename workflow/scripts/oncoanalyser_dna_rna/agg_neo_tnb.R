#!/usr/bin/env Rscript
# Aggregate per-tumour HMF Neo (NeoScorer) peptide-binding scores into the
# cohort tumour neoantigen burden (TNB) tables.
#
# Input:  <src_dir>/neo/scorer/<tumor_lib>.neo.peptide_scores.tsv
#         <src_dir>/neo/scorer/<tumor_lib>.neo.neoepitope.tsv   (clonality join)
# Output: data/wgs/oncoanalyser_dna_rna_aggregate/neo_presented_peptides.csv
#         data/wgs/oncoanalyser_dna_rna_aggregate/neo_tnb_summary.csv
#
# TNB is the neoantigen counterpart of TMB: instead of counting somatic
# mutations, count the mutation-derived peptides predicted to be presented on
# the patient's own HLA alleles. NeoScorer emits one row per
# (peptide x HLA allele) candidate — ~10^5 rows per tumour, ~4x10^7 cohort-wide,
# which is why the raw table is not aggregated (see agg_neo_neoepitopes.R and
# the module README). Restricting to *presented* peptides
# (LikelihoodRank < 0.02, the HMF top-2% presentation-likelihood convention)
# reduces that ~200x to a tractable cohort table, and that presented set is
# exactly the grain the three TNB metrics count.
#
# Definitions (all counts of presented peptide x HLA-allele rows):
#   tnb_total_dna   LikelihoodRank < 0.02
#   tnb_expressed   ... and ExpLikelihoodRank < 0.02
#                   (RNA-informed presentation likelihood — needs the DNA+RNA track)
#   tnb_clonal_core ... and the source variant is a point mutation / indel
#                   (not a fusion) and is clonal (SubclonalLikelihood < 0.85)
#
# These reproduce the TNB_Total_DNA / TNB_Expressed / TNB_Clonal_Core columns in
# the externally supplied cohort summary exactly for 194 of 195 tumours. The one
# deviation (L2500941, donor 4270) is expected: that sample's TNB was computed
# upstream against a cross-tissue RNA pairing (rna_lib L2400848 = tissue 4270-2)
# that commit 40d6b53 corrected to the tissue-matched L2401182. Our value is the
# corrected one. Thresholds are declared as constants below rather than inlined
# so this provenance stays auditable.
#
# SubclonalLikelihood and VariantType live on the neoepitope, not the peptide, so
# they are joined in via NeIds. ~19% of presented peptides carry several NeIds
# (one peptide reachable from several neoepitope contexts of the same variant);
# the first is used. Verified immaterial: across 7,174 multi-NeId presented rows
# no group disagreed on either fusion status or clonality, and the assertion
# below re-checks that on every run rather than trusting the spot check.
#
# Run from project root:
#   Rscript modules/oncoanalyser_dna_rna_aggregate/scripts/agg_neo_tnb.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

suppressPackageStartupMessages({
  library(purrr)
})

# --- TNB definition constants (see header) -----------------------------------
PRESENTED_RANK_MAX   <- 0.02   # LikelihoodRank / ExpLikelihoodRank cut
CLONAL_SUBCLONAL_MAX <- 0.85   # SubclonalLikelihood cut for "clonal"

# Full NeoScorer peptide_scores schema, source order. Kept verbatim for the
# retained (presented) rows; reading all-character preserves the scientific
# notation NeoScorer writes for EffectiveTpm.
PEPTIDE_KEEP <- c("NeIds", "VariantType", "VariantInfo", "Gene", "Allele",
                  "Peptide", "FlankUp", "FlankDown", "Score", "Rank",
                  "Likelihood", "LikelihoodRank", "ExpLikelihood",
                  "ExpLikelihoodRank", "EffectiveTpm", "RecogSim",
                  "OtherAlleleRecogSim", "AlleleCN", "AlleleDisrupted")

# Neoepitope-level context joined onto each presented peptide so consumers can
# re-derive tnb_clonal_core without re-reading the per-sample scorer files.
NEO_JOIN <- c("VariantCopyNumber", "CopyNumber", "SubclonalLikelihood")

num <- function(x) suppressWarnings(as.numeric(x))
# NA-safe strict less-than: an unparseable / missing score never passes a cut.
lt  <- function(x, t) { v <- num(x); !is.na(v) & v < t }

first_neid <- function(ne_ids) sub(";.*$", "", ne_ids)
is_fusion  <- function(variant_type) grepl("FUSION", variant_type, fixed = TRUE)

scorer_path <- function(src_dir, tumor_lib, suffix) {
  path(src_dir, "neo", "scorer", paste0(tumor_lib, ".neo.", suffix, ".tsv"))
}

# Guard the first-NeId shortcut: every NeId on a multi-NeId peptide must agree on
# fusion status and clonality, else the join is ambiguous and the counts drift.
check_neid_agreement <- function(pep, neo, tumor_lib) {
  multi <- pep[grepl(";", pep$NeIds, fixed = TRUE), , drop = FALSE]
  if (nrow(multi) == 0) return(invisible(NULL))
  ctx <- neo %>%
    mutate(.fus = is_fusion(VariantType),
           .clo = lt(SubclonalLikelihood, CLONAL_SUBCLONAL_MAX))
  lookup <- setNames(seq_len(nrow(ctx)), ctx$NeId)
  bad <- keep(strsplit(multi$NeIds, ";", fixed = TRUE), function(ids) {
    idx <- lookup[ids]
    idx <- idx[!is.na(idx)]
    if (length(idx) < 2) return(FALSE)
    length(unique(ctx$.fus[idx])) > 1 || length(unique(ctx$.clo[idx])) > 1
  })
  if (length(bad) > 0) {
    stop(sprintf(paste0("%s: %d presented peptide(s) map to NeIds that disagree ",
                        "on fusion status or clonality — the first-NeId join is ",
                        "no longer safe (e.g. NeIds '%s')"),
                 tumor_lib, length(bad), paste(bad[[1]], collapse = ";")))
  }
  invisible(NULL)
}

# Read one tumour's presented peptides + the summary counts. Returns NULL when
# the scorer file is absent so the tumour is skipped rather than NA-filled.
read_tumour_tnb <- function(src_dir, tumor_lib) {
  p_pep <- scorer_path(src_dir, tumor_lib, "peptide_scores")
  p_neo <- scorer_path(src_dir, tumor_lib, "neoepitope")
  if (!file.exists(p_pep) || !file.exists(p_neo)) return(NULL)

  neo <- read_tsv_quiet(p_neo, col_types = cols(.default = "c"))
  pep <- read_tsv_quiet(p_pep, col_types = cols(.default = "c"))

  n_scored      <- nrow(pep)
  n_neoepitopes <- nrow(neo)

  # Filter to presented rows first: this is the ~200x reduction that makes a
  # cohort-wide peptide table feasible at all.
  pep <- pep[lt(pep$LikelihoodRank, PRESENTED_RANK_MAX), , drop = FALSE]

  if (nrow(pep) == 0 || n_neoepitopes == 0) {
    return(list(
      peptides = tibble(),
      summary  = tibble(n_neoepitopes = n_neoepitopes,
                        n_peptide_allele_scored = n_scored,
                        n_hla_alleles = 0L,
                        n_neoepitopes_presented = 0L,
                        tnb_total_dna = 0L, tnb_expressed = 0L,
                        tnb_clonal_core = 0L,
                        tnb_total_dna_point = 0L, tnb_total_dna_fusion = 0L)
    ))
  }

  check_neid_agreement(pep, neo, tumor_lib)

  ctx <- neo[, intersect(c("NeId", NEO_JOIN), names(neo)), drop = FALSE]
  pep <- pep[, intersect(PEPTIDE_KEEP, names(pep)), drop = FALSE] %>%
    mutate(.neid = first_neid(NeIds)) %>%
    left_join(ctx, by = c(".neid" = "NeId")) %>%
    select(-.neid)

  # Nested, not independent: expressed is a subset of total, clonal_core of
  # expressed. Materialised per row so consumers filter instead of re-deriving.
  pep <- pep %>%
    mutate(
      tnb_expressed   = lt(ExpLikelihoodRank, PRESENTED_RANK_MAX),
      tnb_clonal_core = tnb_expressed &
                        !is_fusion(VariantType) &
                        lt(SubclonalLikelihood, CLONAL_SUBCLONAL_MAX)
    )

  list(
    peptides = pep,
    summary  = tibble(
      n_neoepitopes           = n_neoepitopes,
      n_peptide_allele_scored = n_scored,
      n_hla_alleles           = length(unique(pep$Allele)),
      n_neoepitopes_presented = length(unique(unlist(strsplit(pep$NeIds, ";", fixed = TRUE)))),
      tnb_total_dna           = nrow(pep),
      tnb_expressed           = sum(pep$tnb_expressed),
      tnb_clonal_core         = sum(pep$tnb_clonal_core),
      tnb_total_dna_point     = sum(!is_fusion(pep$VariantType)),
      tnb_total_dna_fusion    = sum(is_fusion(pep$VariantType))
    )
  )
}

ss <- load_samplesheet()

peptide_rows <- vector("list", nrow(ss))
summary_rows <- vector("list", nrow(ss))

for (i in seq_len(nrow(ss))) {
  rec <- ss[i, ]
  res <- read_tumour_tnb(rec$src_dir, rec$tumor_lib)
  if (is.null(res)) next
  ident <- rec %>% select(any_of(IDENTITY_COLS))
  if (nrow(res$peptides) > 0) {
    peptide_rows[[i]] <- bind_cols(ident[rep(1, nrow(res$peptides)), ], res$peptides)
  }
  summary_rows[[i]] <- bind_cols(ident, res$summary)
  if (i %% 25 == 0) message(sprintf("  ... %d/%d tumours", i, nrow(ss)))
}

write_cohort(bind_rows(peptide_rows), "neo_presented_peptides.csv")
write_cohort(bind_rows(summary_rows), "neo_tnb_summary.csv")
