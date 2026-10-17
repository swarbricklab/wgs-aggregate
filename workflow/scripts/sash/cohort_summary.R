#!/usr/bin/env Rscript
# Build a single wide cohort summary CSV from per-tumour sash outputs.
# Sash equivalent of modules/umccrise_aggregate/scripts/cohort_summary.R; the
# portal's `wgs_summary` table is migrating from the umccrise summary to this
# one (sash is the canonical somatic pipeline now).
#
# Inputs (per tumour, all under <src_dir>):
#   purple/<tumor_lib>.purple.purity.tsv
#   purple/<tumor_lib>.purple.qc                                    (key-value)
#   cancer_report/cancer_report_tables/hrd/*-hrdetect.tsv.gz        (globbed)
# Plus the already-aggregated:
#   data/wgs/sash_aggregate/mutational_signatures.csv                   (top-3 SBS)
#
# Output: data/wgs/sash_aggregate/sash_cohort_summary.csv  (one row per tumor_gussid)
#
# Schema gaps vs umccrise_cohort_summary.csv (documented in module README):
#   - No CHORD: sash does not emit CHORD HRD; HRDetect is the only HRD score
#     available here. CHORD columns are simply absent (not NA-filled).
#   - No Conpair: sash does not surface Conpair concordance in cancer_report;
#     concordance columns are absent.
#   - PCGR signatures expose proportions only (no absolute mutation counts),
#     so the top3_sbs string is "SBSx(prop);..." rather than
#     "SBSx(rel_freq);..." with separate contribution counts.
#
# Note on msi_indels_per_mb: PURPLE writes msIndelsPerMb = 0.0000 for every
# tumour in this cohort (checked across all 208 purity TSVs; msStatus is MSS
# throughout), so the column is genuinely zero rather than mis-read. An actual
# MSI indel rate is available from umccrise_aggregate/purple_purity.csv
# (ms_indels_per_mb), which comes from a different PURPLE execution.
#
# Run from project root (usually via modules/sash_aggregate/run.sh):
#   Rscript modules/sash_aggregate/scripts/cohort_summary.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
})

ss <- load_samplesheet()

# Re-read the samplesheet to recover the broader identity set
# (load_samplesheet() returns only IDENTITY_COLS + src_dir; the cohort
# summary wants normal_gussid / normal_tissue_id / batch as well).
identity_extra <- read_csv(SAMPLESHEET, show_col_types = FALSE,
                           comment = "#",
                           col_types = cols(.default = "c")) %>%
  select(any_of(c("tumor_gussid", "normal_gussid",
                  "normal_tissue_id", "batch"))) %>%
  distinct(tumor_gussid, .keep_all = TRUE)

# Post-rename (output) names — read_purity() applies PURPLE_PURITY_COL_MAP
# before this is used. Four of these previously carried PURPLE's source-side
# spelling (polyclonal_proportion / diploid_proportion / ms_indels_per_mb /
# sv_tumor_mutational_burden), so they silently escaped coercion and were
# written as raw PURPLE strings.
PURPLE_PURITY_NUMERIC <- c("purity", "ploidy", "polyclonal_prop",
                           "diploid_prop", "msi_indels_per_mb",
                           "tmb_per_mb", "tml", "sv_tmb")

# Numeric columns sourced from purple.qc (same post-rename convention).
PURPLE_QC_NUMERIC <- c("amber_mean_depth", "loh_percent")

PURPLE_PURITY_COL_MAP <- c(
  status                    = "fit_method",
  purity                    = "purity",
  ploidy                    = "ploidy",
  gender                    = "gender",
  wholeGenomeDuplication    = "wgd",
  polyclonalProportion      = "polyclonal_prop",
  diploidProportion         = "diploid_prop",
  msIndelsPerMb             = "msi_indels_per_mb",
  msStatus                  = "msi_status",
  tmbPerMb                  = "tmb_per_mb",
  tmbStatus                 = "tmb_status",
  tml                       = "tml",
  tmlStatus                 = "tml_status",
  svTumorMutationalBurden   = "sv_tmb"
)

# PURPLE QC is a 2-column key/value file with these keys of interest.
PURPLE_QC_COL_MAP <- c(
  QCStatus                          = "qc_status",
  Contamination                     = "contamination",
  CopyNumberSegments                = "cn_segments",
  UnsupportedCopyNumberSegments     = "cn_segments_unsupp",
  DeletedGenes                      = "deleted_genes",
  GermlineAberrations               = "germline_aberrations",
  AmberMeanDepth                    = "amber_mean_depth",
  LohPercent                        = "loh_percent"
)
# Not mapped: Method (duplicates purity.tsv `status` -> fit_method), AmberGender
# and CobaltGender (purity.tsv `gender` already covers it), and QC's own Purity
# (same value as purity.tsv, fewer decimals).

read_purity <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "purple", paste0(tumor_lib, ".purple.purity.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())
  present <- intersect(names(PURPLE_PURITY_COL_MAP), names(df))
  df <- df[, present, drop = FALSE]
  names(df) <- unname(PURPLE_PURITY_COL_MAP[present])
  df %>% slice(1)
}

read_qc <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "purple", paste0(tumor_lib, ".purple.qc"))
  if (!file.exists(p)) return(tibble())
  kv <- read_tsv_quiet(p, col_names = c("key", "value"),
                       col_types = cols(.default = "c"))
  if (nrow(kv) == 0) return(tibble())
  pres <- kv %>% filter(key %in% names(PURPLE_QC_COL_MAP))
  if (nrow(pres) == 0) return(tibble())
  pres %>%
    pivot_wider(names_from = "key", values_from = "value") %>%
    rename_with(~ PURPLE_QC_COL_MAP[.x])
}

read_hrdetect <- function(src_dir) {
  hits <- tryCatch(dir_ls(path(src_dir, "cancer_report", "cancer_report_tables", "hrd"),
                          regexp = "-hrdetect\\.tsv\\.gz$"),
                   error = function(e) character(0))
  if (length(hits) == 0) return(tibble(hrdetect_probability = NA_character_))
  df <- read_tsv_quiet(hits[[1]], col_types = cols(.default = "c"))
  if (nrow(df) == 0 || !"Probability" %in% names(df)) {
    return(tibble(hrdetect_probability = NA_character_))
  }
  tibble(hrdetect_probability = df$Probability[[1]])
}

# Per-tumour bind: identity + purity + qc + hrdetect (purity is the "usable" key)
per_sample <- vector("list", nrow(ss))
for (i in seq_len(nrow(ss))) {
  rec <- ss[i, ]
  pur <- read_purity(rec$src_dir, rec$tumor_lib)
  if (nrow(pur) == 0) next
  qc <- read_qc(rec$src_dir, rec$tumor_lib)
  hrd <- read_hrdetect(rec$src_dir)
  per_sample[[i]] <- bind_cols(
    rec %>% select(any_of(IDENTITY_COLS)),
    pur, qc, hrd
  )
}

base <- bind_rows(per_sample) %>%
  left_join(identity_extra, by = "tumor_gussid") %>%
  select(any_of(c("tumor_gussid", "normal_gussid", "donor_id",
                  "tumor_tissue_id", "normal_tissue_id",
                  "tumor_lib", "normal_lib", "batch", "sash_run_id")),
         everything())

# Coerce numerics
num_present <- intersect(c(PURPLE_PURITY_NUMERIC, PURPLE_QC_NUMERIC), names(base))
base <- base %>%
  mutate(across(all_of(num_present), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(hrdetect_probability = suppressWarnings(as.numeric(hrdetect_probability)))

# ---- top-3 SBS signatures from data/wgs/sash_aggregate/mutational_signatures.csv ----
# PCGR proportions only - no separate contribution counts (unlike umccrise).

sigs_path <- file.path(OUT_DIR, "mutational_signatures.csv")
if (file.exists(sigs_path)) {
  sigs <- read_csv(sigs_path, show_col_types = FALSE,
                   col_types = cols(.default = "c"))
  top3_sbs <- sigs %>%
    mutate(prop_signature = suppressWarnings(as.numeric(prop_signature))) %>%
    arrange(tumor_gussid, desc(prop_signature)) %>%
    group_by(tumor_gussid) %>%
    slice_head(n = 3) %>%
    summarise(
      top3_sbs = paste0(signature_id, "(", round(prop_signature, 3), ")",
                       collapse = ";"),
      .groups = "drop"
    )
  base <- base %>% left_join(top3_sbs, by = "tumor_gussid")
} else {
  message("INFO: data/wgs/sash_aggregate/mutational_signatures.csv not present yet — top3_sbs omitted. Run agg_mutational_signatures.R first.")
}

# Render numerics in fixed notation. The columns above are coerced with
# as.numeric() so the schema is typed, but readr writes small doubles in
# scientific notation (polyclonal_prop 0.0002 -> "2e-04"), which is a change
# from what PURPLE emits and from what this file published before the coercion
# was fixed. Format per element (not vector-wide, which would pad every value
# to a common width) so each keeps its own precision.
fixed_notation <- function(x) {
  ifelse(is.na(x), NA_character_,
         vapply(x, function(v) format(v, scientific = FALSE, trim = TRUE,
                                      digits = 15),
                character(1)))
}
fmt_cols <- intersect(c(num_present, "hrdetect_probability"), names(base))
base <- base %>% mutate(across(all_of(fmt_cols), fixed_notation))

dir_create(OUT_DIR)
out_csv <- file.path(OUT_DIR, "sash_cohort_summary.csv")
write_csv(base, out_csv, na = "NA")
message(sprintf("Wrote %d rows x %d cols to %s",
                nrow(base), ncol(base), out_csv))
