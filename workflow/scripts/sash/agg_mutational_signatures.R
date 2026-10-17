#!/usr/bin/env Rscript
# Aggregate per-tumour PCGR mutational signature fits into a cohort CSV.
# Input:  <src_dir>/smlv_somatic/report/pcgr/<tumor_lib>.pcgr_acmg.grch38.mutational_signatures.tsv
# Output: data/wgs/sash_aggregate/mutational_signatures.csv
#
# Schema differs from umccrise_aggregate/signatures.csv: PCGR fits to a single
# SBS reference panel (no separate sig_type ∈ snv_2020/snv_2015/dbs/indel) and
# emits only relative proportions, not absolute counts. The verbose
# `reference_signatures` column (long comma-joined list of every panel
# signature, repeated identically across rows) is dropped — the
# `reference_collection` column (e.g. COSMIC_v32) is the load-bearing pointer.
# Likewise `sample_id` (= tumor_lib) is dropped as redundant.
#
# Run from project root:
#   Rscript modules/sash_aggregate/scripts/agg_mutational_signatures.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

NUMERIC_COLS <- c("prop_signature", "fitting_accuracy")
LOGICAL_COLS <- c("all_reference_signatures")

# Snake_case identity already; the source TSV uses snake_case throughout.
KEEP_COLS <- c("signature_id", "prop_signature", "group",
               "all_reference_signatures", "tumor_type",
               "reference_collection", "fitting_accuracy")

read_signatures_tsv <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "smlv_somatic", "report", "pcgr",
            paste0(tumor_lib, ".pcgr_acmg.grch38.mutational_signatures.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())

  present <- intersect(KEEP_COLS, names(df))
  df[, present, drop = FALSE]
}

rows <- ss %>%
  rowwise() %>%
  reframe(
    pick(any_of(IDENTITY_COLS)),
    read_signatures_tsv(src_dir, tumor_lib)
  )

if (nrow(rows) > 0) {
  num_present <- intersect(NUMERIC_COLS, names(rows))
  rows <- rows %>%
    mutate(across(all_of(num_present), ~ suppressWarnings(as.numeric(.))))

  log_present <- intersect(LOGICAL_COLS, names(rows))
  rows <- rows %>%
    mutate(across(all_of(log_present), ~ tolower(.) == "true"))
}

dir_create(OUT_DIR)
out_csv <- file.path(OUT_DIR, "mutational_signatures.csv")
write_csv(rows, out_csv, na = "NA")
cat(sprintf("Wrote %d signature row(s) from %d tumour(s) to %s\n",
            nrow(rows), length(unique(rows$tumor_gussid)), out_csv))
