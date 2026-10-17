#!/usr/bin/env Rscript
# Aggregate per-tumour PURPLE somatic copy number segment tables into a cohort CSV.
# Input:  <src_dir>/purple/<tumor_lib>.purple.cnv.somatic.tsv
# Output: data/wgs/sash_aggregate/purple_cnv_somatic.csv
#
# Full segment track: one row per (tumour, CN segment) — no thresholding.
# ~65-4800 segments per tumour x ~208 tumours. This is the segment-level
# companion to purple_cnv_gene.csv and the natural comparator to the (coarser)
# SNP-array segmented copy number. Downstream consumers compare segment
# copyNumber against per-tumour ploidy as needed.
#
# Run from project root:
#   Rscript modules/sash_aggregate/scripts/agg_purple_cnv_somatic.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

NUMERIC_COLS <- c("start", "end", "copy_number", "baf_count", "observed_baf",
                  "baf", "depth_window_count", "gc_content",
                  "min_start", "max_start",
                  "minor_allele_copy_number", "major_allele_copy_number")

# Mapping from PURPLE camelCase -> snake_case. Explicit so we'd notice if a
# future sash bump adds/drops columns.
COL_MAP <- c(
  chromosome            = "chromosome",
  start                 = "start",
  end                   = "end",
  copyNumber            = "copy_number",
  bafCount              = "baf_count",
  observedBAF           = "observed_baf",
  baf                   = "baf",
  segmentStartSupport   = "segment_start_support",
  segmentEndSupport     = "segment_end_support",
  method                = "method",
  depthWindowCount      = "depth_window_count",
  gcContent             = "gc_content",
  minStart              = "min_start",
  maxStart              = "max_start",
  minorAlleleCopyNumber = "minor_allele_copy_number",
  majorAlleleCopyNumber = "major_allele_copy_number"
)

read_cnv_somatic_tsv <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "purple",
            paste0(tumor_lib, ".purple.cnv.somatic.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())

  present <- intersect(names(COL_MAP), names(df))
  df <- df[, present, drop = FALSE]
  names(df) <- unname(COL_MAP[present])
  df
}

# Stream per-sample reads + bind once, matching agg_purple_cnv.R.
per_sample <- vector("list", nrow(ss))
for (i in seq_len(nrow(ss))) {
  rec <- ss[i, ]
  cnv <- read_cnv_somatic_tsv(rec$src_dir, rec$tumor_lib)
  if (nrow(cnv) == 0) next
  per_sample[[i]] <- cnv %>%
    mutate(across(everything(), as.character)) %>%
    bind_cols(rec %>% select(any_of(IDENTITY_COLS)) %>% slice(rep(1, nrow(cnv))))
}

rows <- bind_rows(per_sample) %>%
  select(any_of(IDENTITY_COLS), everything(), -any_of("src_dir"))

if (nrow(rows) > 0) {
  num_present <- intersect(NUMERIC_COLS, names(rows))
  rows <- rows %>%
    mutate(across(all_of(num_present), ~ suppressWarnings(as.numeric(.))))
}

dir_create(OUT_DIR)
out_csv <- file.path(OUT_DIR, "purple_cnv_somatic.csv")
write_csv(rows, out_csv, na = "NA")
cat(sprintf("Wrote %d row(s) from %d tumour(s) to %s\n",
            nrow(rows), length(unique(rows$tumor_gussid)), out_csv))
