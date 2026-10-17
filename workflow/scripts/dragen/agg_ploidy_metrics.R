#!/usr/bin/env Rscript
# Aggregate ploidy estimation metrics across all runs (germline file)
# Input:  {germline}.ploidy_estimation_metrics.csv
# Output: data/wgs/dragen_aggregate/ploidy_estimation_metrics.csv
#
# Uses the germline run file (unambiguous constitutional ploidy).
# Key metrics: autosomal/X/Y coverage, per-chromosome ratios, sex call.
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_ploidy_metrics.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

clean_metric <- function(x) {
  x |> tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("_+$", "", x = _)
}

ploidy <- load_run_metadata() |>
  build_metric_paths(pipeline = "germline", suffix = "ploidy_estimation_metrics.csv") |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_dragen_metrics(metric_path) |>
      mutate(
        col   = clean_metric(metric),
        value = value
      ) |>
      select(col, value)
  ) |>
  # ploidy_estimation is a string (e.g. "XX"); keep others as numeric
  mutate(value = if_else(col == "ploidy_estimation", value, suppressWarnings(as.character(as.numeric(value))))) |>
  pivot_wider(names_from = col, values_from = value) |>
  mutate(across(-c(run_id, donor_id, tissue_id, tumor_library, germline_library, ploidy_estimation),
                as.numeric))

write_csv(ploidy, file.path(OUT_DIR, "ploidy_estimation_metrics.csv"))
cat(sprintf("Wrote %d rows to ploidy_estimation_metrics.csv\n", nrow(ploidy)))
