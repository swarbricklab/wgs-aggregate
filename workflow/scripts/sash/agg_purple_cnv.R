#!/usr/bin/env Rscript
# Aggregate per-tumour PURPLE per-gene copy number tables into a cohort CSV.
# Input:  <src_dir>/purple/<tumor_lib>.purple.cnv.gene.tsv
# Output: data/wgs/sash_aggregate/purple_cnv_gene.csv
#
# Full matrix: one row per (tumour, gene) — no thresholding. ~39k genes per
# tumour x ~209 tumours = ~8.2M rows. Downstream consumers apply amp/del
# thresholds against ploidy as needed.
#
# Run from project root:
#   Rscript modules/sash_aggregate/scripts/agg_purple_cnv.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

NUMERIC_COLS <- c("start", "end", "min_copy_number", "max_copy_number",
                  "somatic_regions", "min_regions",
                  "min_region_start", "min_region_end",
                  "min_minor_allele_copy_number", "depth_window_count")
LOGICAL_COLS <- c("is_canonical")

# Mapping from PURPLE camelCase -> snake_case. Explicit so we'd notice if a
# future sash bump adds/drops columns (e.g. gcContent appeared in an interim
# sash version on now-superseded Nov-2025 runs; the canonical post-cleanup
# cohort omits it).
COL_MAP <- c(
  chromosome              = "chromosome",
  start                   = "start",
  end                     = "end",
  gene                    = "gene",
  minCopyNumber           = "min_copy_number",
  maxCopyNumber           = "max_copy_number",
  somaticRegions          = "somatic_regions",
  transcriptId            = "transcript",
  isCanonical             = "is_canonical",
  chromosomeBand          = "chromosome_band",
  minRegions              = "min_regions",
  minRegionStart          = "min_region_start",
  minRegionEnd            = "min_region_end",
  minRegionStartSupport   = "min_region_start_support",
  minRegionEndSupport     = "min_region_end_support",
  minRegionMethod         = "min_region_method",
  minMinorAlleleCopyNumber = "min_minor_allele_copy_number",
  depthWindowCount        = "depth_window_count"
)

read_cnv_tsv <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "purple",
            paste0(tumor_lib, ".purple.cnv.gene.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())

  present <- intersect(names(COL_MAP), names(df))
  df <- df[, present, drop = FALSE]
  names(df) <- unname(COL_MAP[present])
  df
}

# Stream per-sample reads + bind once. rowwise() + reframe() like agg_drivers.R
# works but is slow on ~8M rows; an explicit lapply + bind_rows is a few times
# faster and clearer for this scale.
per_sample <- vector("list", nrow(ss))
for (i in seq_len(nrow(ss))) {
  rec <- ss[i, ]
  cnv <- read_cnv_tsv(rec$src_dir, rec$tumor_lib)
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

  log_present <- intersect(LOGICAL_COLS, names(rows))
  rows <- rows %>%
    mutate(across(all_of(log_present), ~ tolower(.) == "true"))
}

dir_create(OUT_DIR)
out_csv <- file.path(OUT_DIR, "purple_cnv_gene.csv")
write_csv(rows, out_csv, na = "NA")
cat(sprintf("Wrote %d row(s) from %d tumour(s) to %s\n",
            nrow(rows), length(unique(rows$tumor_gussid)), out_csv))
