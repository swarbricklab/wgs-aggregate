#!/usr/bin/env Rscript
# Aggregate germline ROH and het/hom ratio metrics
# Input:  {germline}.roh_metrics.csv            (2 rows)
#         {germline}.vc_hethom_ratio_metrics.csv (144 rows: per-chr het/hom)
# Output: data/wgs/dragen_aggregate/germline_roh_hethom_metrics.csv
#
# For hethom: extracts PREFILTER Het/Hom ratio for autosomes (chr1-22) + X + Y.
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_roh_hethom_metrics.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

clean_metric <- function(x) {
  x |> tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("_+$", "", x = _)
}

meta <- load_run_metadata()

# ROH: 2 rows
roh <- meta |>
  build_metric_paths(pipeline = "germline", suffix = "roh_metrics.csv") |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_dragen_metrics(metric_path) |>
      mutate(col = clean_metric(metric), value = as.numeric(value)) |>
      select(col, value)
  ) |>
  pivot_wider(names_from = col, values_from = value)

# Het/Hom: extract PREFILTER Het/Hom ratio rows only (rows like "1 Het/Hom ratio", "2 Het/Hom ratio", ...)
# The metric field looks like "1 Het/Hom ratio", "2 Het/Hom ratio", ..., "X Het/Hom ratio", "Y Het/Hom ratio"
hethom <- meta |>
  build_metric_paths(pipeline = "germline", suffix = "vc_hethom_ratio_metrics.csv") |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_dragen_metrics(metric_path) |>
      filter(
        section == "VARIANT CALLER PREFILTER",
        grepl("Het/Hom ratio", metric)
      ) |>
      mutate(
        chr   = sub("^(\\S+)\\s+Het/Hom ratio$", "\\1", metric),
        col   = paste0("chr", chr, "_hethom_ratio"),
        value = as.numeric(value)
      ) |>
      select(col, value)
  ) |>
  pivot_wider(names_from = col, values_from = value)

result <- roh |>
  left_join(hethom, by = c("run_id", "donor_id", "tissue_id", "tumor_library", "germline_library"))

write_csv(result, file.path(OUT_DIR, "germline_roh_hethom_metrics.csv"))
cat(sprintf("Wrote %d rows to germline_roh_hethom_metrics.csv\n", nrow(result)))
