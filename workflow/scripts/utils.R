# Consolidated utilities for wgs-aggregate (replaces the 4 per-module utils.R).
#
# Config comes from env vars set by the calling Snakemake rule (defaults match the
# atlas layout so scripts also run standalone from a dataset root):
#   WGS_KIND        sash | oncoanalyser | oncoanalyser_dna_rna   (selects the loader)
#   WGS_PIPELINE_DIR per-sample pipeline outputs  (e.g. data/wgs/sash)
#   WGS_OUT_DIR     cohort-table output dir       (e.g. data/wgs/sash_aggregate)
#   WGS_SAMPLESHEET consolidated samplesheet for this pipeline
#
# No index.csv / config/projects: all identity comes from the samplesheet.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(fs)
})

.cfg <- function(k, d = "") { v <- Sys.getenv(k); if (nzchar(v)) v else d }
KIND        <- .cfg("WGS_KIND", "sash")
DATA_DIR    <- .cfg("WGS_PIPELINE_DIR", file.path("data/wgs", KIND))
OUT_DIR     <- .cfg("WGS_OUT_DIR", file.path("data/wgs", paste0(KIND, "_aggregate")))
SAMPLESHEET <- .cfg("WGS_SAMPLESHEET", file.path("config", KIND, "samplesheet.csv"))

# Identity columns preserved in every per-source output, per pipeline kind.
.IDENTITY <- list(
  sash                 = c("tumor_gussid","donor_id","tumor_tissue_id","tumor_lib","normal_lib","sash_run_id"),
  oncoanalyser         = c("tumor_gussid","donor_id","tumor_tissue_id","tumor_lib","normal_lib","onco_run_id"),
  oncoanalyser_dna_rna = c("tumor_gussid","donor_id","tumor_tissue_id","tumor_lib","normal_lib","rna_lib","dna_rna_run_id")
)
IDENTITY_COLS <- .IDENTITY[[KIND]]

read_tsv_quiet <- function(path, ...) {
  suppressMessages(suppressWarnings(read_tsv(path, show_col_types = FALSE, ...)))
}

# --- per-row pipeline-output dir resolvers ------------------------------------
.resolve_sash <- function(ss) ss %>% rename(src_dir = sash_dir)

.resolve_oncoanalyser <- function(ss) {
  one <- function(run, t, n) {
    if (is.na(run) || run == "") return(NA_character_)
    for (sep in c("__", "_")) { d <- path(DATA_DIR, run, paste0(t, sep, n)); if (dir_exists(d)) return(as.character(d)) }
    NA_character_
  }
  ss$src_dir <- mapply(one, ss$onco_run_id, ss$tumor_lib, ss$normal_lib); ss
}

.resolve_dna_rna <- function(ss) {
  one <- function(run, lib_dir) {
    if (is.na(run) || run == "" || is.na(lib_dir) || lib_dir == "") return(NA_character_)
    d <- path(DATA_DIR, run, lib_dir); if (dir_exists(d)) as.character(d) else NA_character_
  }
  ss$src_dir <- mapply(one, ss$dna_rna_run_id, ss$lib_dir); ss
}

.REQUIRED <- list(
  sash                 = c("tumor_gussid","tumor_lib","normal_lib","sash_dir"),
  oncoanalyser         = c("tumor_gussid","tumor_lib","normal_lib","onco_run_id"),
  oncoanalyser_dna_rna = c("tumor_gussid","tumor_lib","normal_lib","rna_lib","dna_rna_run_id","lib_dir")
)

# Load usable rows + a resolved `src_dir` column (the per-sample pipeline output dir).
load_samplesheet <- function(path = SAMPLESHEET, kind = KIND) {
  ss <- read_csv(path, show_col_types = FALSE, comment = "#", col_types = cols(.default = "c"))
  missing <- setdiff(.REQUIRED[[kind]], names(ss))
  if (length(missing)) stop("Samplesheet missing required columns: ", paste(missing, collapse = ", "))

  if ("skip_reason" %in% names(ss)) {
    n <- sum(!is.na(ss$skip_reason) & ss$skip_reason != "")
    ss <- ss %>% filter(is.na(skip_reason) | skip_reason == "")
    if (n) message(sprintf("INFO: %d row(s) skipped via skip_reason", n))
  }

  ss <- switch(kind,
    sash                 = .resolve_sash(ss),
    oncoanalyser         = .resolve_oncoanalyser(ss),
    oncoanalyser_dna_rna = .resolve_dna_rna(ss),
    stop("unknown WGS_KIND: ", kind))

  n_bad <- sum(is.na(ss$src_dir) | !dir_exists(ss$src_dir))
  if (n_bad) message(sprintf("WARN: skipping %d row(s) whose pipeline dir is absent", n_bad))
  ss <- ss %>% filter(!is.na(src_dir) & dir_exists(src_dir))

  message(sprintf("Loaded %d usable row(s) from %s", nrow(ss), path))
  ss %>% select(any_of(IDENTITY_COLS), src_dir)
}

# Write a cohort table with identity cols leading.
write_cohort <- function(df, filename, out_dir = OUT_DIR) {
  dir_create(out_dir)
  out <- path(out_dir, filename)
  lead <- intersect(IDENTITY_COLS, names(df))
  df <- df %>% select(all_of(lead), everything())
  write_csv(df, out, na = "")
  message(sprintf("Wrote %d row(s) -> %s", nrow(df), out))
}

# --- DRAGEN metrics helpers (decoupled from index.csv) ------------------------
# tissue_id now comes from the sash samplesheet (tumor_lib -> tumor_tissue_id);
# the registry-only wgs_sbj/wgs_prj identifiers are intentionally dropped.
load_run_metadata <- function(runs_csv = .cfg("WGS_DRAGEN_RUNS", "config/tumor-normal/dragen_runs.csv"),
                              sash_ss  = .cfg("WGS_SASH_SAMPLESHEET", "config/sash/samplesheet.csv")) {
  runs <- read_csv(runs_csv, col_types = cols(.default = "c"))
  ss   <- read_csv(sash_ss,  col_types = cols(.default = "c"))
  tumor_meta <- ss %>% transmute(tumor_library = tumor_lib, tissue_id = tumor_tissue_id) %>% distinct()
  runs %>% left_join(tumor_meta, by = "tumor_library")
}

read_dragen_metrics <- function(path) {
  as_tibble(read.csv(path, header = FALSE,
    col.names = c("section","sample","metric","value","pct"),
    colClasses = "character", fill = TRUE, strip.white = TRUE))
}

build_metric_paths <- function(meta, pipeline = "somatic", suffix,
                               data_dir = .cfg("WGS_DRAGEN_DIR", "data/wgs/tumor-normal")) {
  if (pipeline == "somatic") { lib_col <- "tumor_library"; dir_suffix <- "_dragen_somatic" }
  else                       { lib_col <- "germline_library"; dir_suffix <- "_dragen_germline" }
  meta %>% mutate(metric_path = file.path(data_dir, run_id,
      paste0(.data[[lib_col]], dir_suffix), paste0(.data[[lib_col]], ".", suffix)))
}
