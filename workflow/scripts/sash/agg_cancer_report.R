#!/usr/bin/env Rscript
# Aggregate the per-tumour sash cancer-report tables into two cohort CSVs.
#
# Inputs (per pair, under <src_dir>/cancer_report/cancer_report_tables/):
#   <donor>_<tumor_lib>-qc_summary.tsv.gz    variable/value/details dashboard
#   sigs/<donor>_<tumor_lib>-<sig_type>.tsv.gz   Rank/Signature/Contribution/RelFreq
# Outputs:
#   data/wgs/sash_aggregate/cancer_report_qc_summary.csv   (one row per tumour)
#   data/wgs/sash_aggregate/cancer_report_signatures.csv   (tumour x sig_type x signature)
#
# Why this exists
# ---------------
# sash's cancer report already carries several quantities that have no other
# home in this repo: the DRAGEN HRD component, the per-process somatic SNV
# counts, and the SV type breakdown. Downstream consumers were reconstructing
# them by hand from a spreadsheet because nothing here exposed them. Verified
# 2026-08-27 against that spreadsheet: 24 of its columns reproduce byte-exactly
# from these two tables (see modules/sash_aggregate/README.md).
#
# Relationship to the other outputs in this module:
#   - sash_cohort_summary.csv parses PURPLE's own TSVs and is the place to go
#     for atomic purity/ploidy/TMB/TML. It does NOT carry snvs_somatic,
#     svs_melted, hypermutated, mutsigs_*, cnvs_somatic or the DRAGEN HRD score.
#     Those are only in the cancer report, hence this stage.
#   - mutational_signatures.csv is PCGR's fit (proportions only, one panel).
#     cancer_report_signatures.csv is sigrap's fit (absolute contributions +
#     relative frequencies, four panels). They are DIFFERENT FITS and disagree;
#     do not mix them. The composite strings in the qc summary and in the
#     reference spreadsheet come from the sigrap fit, i.e. from this table.
#
# Composite cells are emitted twice: verbatim in the column named after the
# report variable, and parsed into atomic columns alongside. The verbatim form
# is what makes the byte-comparison against the reference spreadsheet possible;
# the parsed form is so consumers stop writing their own regexes.
#
# Coverage note: a tumour with no somatic calls has no cancer report at all
# (e.g. L2501198, tissue 4338-2, FAIL_NO_TUMOR), so it is absent here rather
# than present with NA. The `dbs` signature panel is fitted for fewer tumours
# than the other three; absence there is not an error.
#
# Run from project root (usually via modules/sash_aggregate/run.sh):
#   Rscript modules/sash_aggregate/scripts/agg_cancer_report.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
})

ss <- load_samplesheet()

REPORT_SUBDIR <- file.path("cancer_report", "cancer_report_tables")

# The report filenames are prefixed with the donor id as sash saw it, which is
# not always the samplesheet's donor_id (tissue suffixes differ). Glob on the
# tumour library, which is unambiguous, rather than reconstructing the name.
find_report_file <- function(src_dir, tumor_lib, subdir = NULL, suffix) {
  dir <- if (is.null(subdir)) path(src_dir, REPORT_SUBDIR) else path(src_dir, REPORT_SUBDIR, subdir)
  if (!dir_exists(dir)) return(NA_character_)
  hits <- dir_ls(dir, glob = sprintf("*_%s-%s.tsv.gz", tumor_lib, suffix), fail = FALSE)
  if (length(hits) == 0) return(NA_character_)
  if (length(hits) > 1) {
    stop(sprintf("Ambiguous %s report for %s: %d matches in %s",
                 suffix, tumor_lib, length(hits), dir))
  }
  as.character(hits[[1]])
}

# ---------------------------------------------------------------- qc summary

# Report variable -> output column. Kept explicit rather than mechanically
# snake_cased so that a new variable appearing upstream is surfaced (below)
# instead of silently landing under a generated name.
QC_VARS <- c(
  "QC_Status"           = "qc_status",
  "Purity"              = "purity",
  "Ploidy"              = "ploidy",
  "Hypermutated"        = "hypermutated",
  "Gender"              = "gender",
  "MutSigs (Old)"       = "mutsigs_old",
  "MutSigs (New)"       = "mutsigs_new",
  "HRD"                 = "hrd",
  "WGD"                 = "wgd",
  "PolyclonalProp"      = "polyclonal_prop",
  "DiploidyProp"        = "diploid_prop",
  "TMB"                 = "tmb",
  "TML"                 = "tml",
  "TMB-SV"              = "tmb_sv",
  "Method"              = "fit_method",
  "CopyNumberSegments"  = "cn_segments",
  "DeletedGenes"        = "deleted_genes",
  "Contamination"       = "contamination",
  "GermlineAberrations" = "germline_aberrations",
  "SNVs (Somatic)"      = "snvs_somatic",
  "SVs (unmelted)"      = "svs_unmelted",
  "SVs (melted)"        = "svs_melted",
  "AmberMeanDepth"      = "amber_mean_depth",
  "CNVs (Somatic)"      = "cnvs_somatic",
  "Genome"              = "genome"
)

read_qc_summary <- function(src_dir, tumor_lib) {
  p <- find_report_file(src_dir, tumor_lib, suffix = "qc_summary")
  if (is.na(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0 || !all(c("variable", "value") %in% names(df))) return(tibble())

  unknown <- setdiff(df$variable, names(QC_VARS))
  if (length(unknown) > 0) {
    message(sprintf("WARN: %s: unmapped qc_summary variable(s): %s",
                    tumor_lib, paste(unknown, collapse = ", ")))
  }

  vals <- setNames(as.list(df$value), df$variable)
  present <- intersect(names(QC_VARS), names(vals))
  out <- as_tibble(vals[present])
  names(out) <- unname(QC_VARS[present])
  out
}

# --- parsers for the composite cells ---------------------------------------
# Each returns a one-row tibble; a cell that does not match yields NA rather
# than a partial parse, so a format change upstream shows up as an empty column
# rather than as silently wrong numbers.

num_or_na <- function(x) suppressWarnings(as.numeric(str_remove_all(x, ",")))

# "dragen: 50,251; sage: 50,251; annotated: 50,251; filter_pass: 7,246"
parse_snvs <- function(s) {
  get <- function(k) num_or_na(str_match(s, paste0(k, ":\\s*([\\d,]+)"))[, 2])
  tibble(snv_dragen = get("dragen"), snv_sage = get("sage"),
         snv_annotated = get("annotated"), snv_filter_pass = get("filter_pass"))
}

# "Sum: 1554, SGL: 19, PURPLE_inf: 111, INV: 271, INS: 2, DUP: 93, DEL: 832, BND: 226"
#
# sash OMITS a type whose count is zero rather than writing `BND: 0`, so an
# absent field means zero, not unknown. Verified across all 208 tumours: the
# type fields sum exactly to `Sum:` in every case, for both the melted and
# unmelted strings. Zero-filling is therefore safe and NOT zero-filling would
# silently turn "no BNDs" into a missing value in 6 tumours.
# `Sum` itself is never zero-filled: if it is absent the whole cell failed to
# parse, and every field stays NA so the failure is visible.
SV_TYPES  <- c("SGL", "PURPLE_inf", "INV", "INS", "DUP", "DEL", "BND")
parse_svs <- function(s, prefix) {
  total <- num_or_na(str_match(s, "\\bSum:\\s*(\\d+)")[, 2])
  parsed <- !is.na(total)
  vals <- lapply(SV_TYPES, function(k) {
    v <- num_or_na(str_match(s, paste0("\\b", k, ":\\s*(\\d+)"))[, 2])
    ifelse(parsed & is.na(v), 0, v)
  })
  vals <- c(list(total), vals)
  names(vals) <- paste0(prefix, "_", tolower(c("Sum", SV_TYPES)))
  as_tibble(vals)
}

# "DRAGEN: 61; CHORD: 0.954; HRDetect: 0.876" -- sash emits DRAGEN + HRDetect
# only; a CHORD term, if ever present, is read but CHORD's own table
# (oncoanalyser_aggregate/chord_hrd.csv) remains the source of record for it.
parse_hrd <- function(s) {
  get <- function(k) num_or_na(str_match(s, paste0(k, ":\\s*([\\d.eE+-]+)"))[, 2])
  tibble(hrd_dragen_score = get("DRAGEN"), hrd_hrdetect_prob = get("HRDetect"))
}

# "Min: -0.9; Max: 12.3; N: 1455"
parse_cnvs <- function(s) {
  get <- function(k) num_or_na(str_match(s, paste0(k, ":\\s*(-?[\\d.]+)"))[, 2])
  tibble(cnv_min = get("Min"), cnv_max = get("Max"), cnv_n = get("N"))
}

qc <- ss %>%
  rowwise() %>%
  reframe(
    pick(any_of(IDENTITY_COLS)),
    read_qc_summary(src_dir, tumor_lib)
  )

if (nrow(qc) > 0) {
  blank <- function(col) if (col %in% names(qc)) qc[[col]] else NA_character_
  qc <- bind_cols(
    qc,
    parse_snvs(blank("snvs_somatic")),
    parse_svs(blank("svs_melted"), "sv_melted"),
    parse_svs(blank("svs_unmelted"), "sv_unmelted"),
    parse_hrd(blank("hrd")),
    parse_cnvs(blank("cnvs_somatic"))
  )
}

dir_create(OUT_DIR)
qc_csv <- file.path(OUT_DIR, "cancer_report_qc_summary.csv")
write_csv(qc, qc_csv, na = "NA")
cat(sprintf("Wrote %d tumour row(s) x %d column(s) to %s\n",
            nrow(qc), ncol(qc), qc_csv))

missing_report <- setdiff(ss$tumor_lib, qc$tumor_lib)
if (length(missing_report) > 0) {
  cat(sprintf("  no cancer report for %d tumour(s): %s\n",
              length(missing_report), paste(missing_report, collapse = ", ")))
}

# ---------------------------------------------------------------- signatures

# sigrap fits four reference panels; `dbs` is fitted for fewer tumours than the
# rest. Names are the file suffixes, carried through as `sig_type` so the
# schema matches umccrise_aggregate/signatures.csv.
SIG_TYPES <- c("snv_2020", "snv_2015", "dbs", "indel")

read_sigs <- function(src_dir, tumor_lib) {
  out <- lapply(SIG_TYPES, function(st) {
    p <- find_report_file(src_dir, tumor_lib, subdir = "sigs", suffix = st)
    if (is.na(p)) return(NULL)
    df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
    if (nrow(df) == 0) return(NULL)
    tibble(
      sig_type     = st,
      rank         = suppressWarnings(as.integer(df$Rank)),
      signature    = df$Signature,
      contribution = suppressWarnings(as.numeric(df$Contribution)),
      rel_freq     = suppressWarnings(as.numeric(df$RelFreq))
    )
  })
  bind_rows(out)
}

sigs <- ss %>%
  rowwise() %>%
  reframe(
    pick(any_of(IDENTITY_COLS)),
    read_sigs(src_dir, tumor_lib)
  )

sig_csv <- file.path(OUT_DIR, "cancer_report_signatures.csv")
write_csv(sigs, sig_csv, na = "NA")
cat(sprintf("Wrote %d signature row(s) from %d tumour(s) to %s\n",
            nrow(sigs), length(unique(sigs$tumor_gussid)), sig_csv))
if (nrow(sigs) > 0) {
  per_type <- sigs %>% count(sig_type, name = "rows") %>%
    left_join(sigs %>% distinct(sig_type, tumor_lib) %>% count(sig_type, name = "tumours"),
              by = "sig_type")
  for (i in seq_len(nrow(per_type))) {
    cat(sprintf("  %-9s %6d row(s) from %d tumour(s)\n",
                per_type$sig_type[i], per_type$rows[i], per_type$tumours[i]))
  }
}
