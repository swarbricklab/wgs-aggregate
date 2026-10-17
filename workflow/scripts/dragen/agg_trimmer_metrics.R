#!/usr/bin/env Rscript
# Aggregate trimmer metrics across all runs
# Input:  {tumor}.trimmer_metrics.csv   (somatic run)
#         {germline}.trimmer_metrics.csv (germline run)
# Output: data/wgs/dragen_aggregate/somatic_trimmer_metrics.csv
#         data/wgs/dragen_aggregate/germline_trimmer_metrics.csv
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_trimmer_metrics.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

clean_metric <- function(x) {
  x |> tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("_+$", "", x = _)
}

agg_trimmer <- function(meta) {
  meta |>
    rowwise() |>
    reframe(
      pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
      read_dragen_metrics(metric_path) |>
        mutate(col = clean_metric(metric), value = as.numeric(value)) |>
        select(col, value)
    ) |>
    pivot_wider(names_from = col, values_from = value)
}

somatic_trim <- load_run_metadata() |>
  build_metric_paths(pipeline = "somatic", suffix = "trimmer_metrics.csv") |>
  agg_trimmer()
write_csv(somatic_trim, file.path(OUT_DIR, "somatic_trimmer_metrics.csv"))
cat(sprintf("Wrote %d rows to somatic_trimmer_metrics.csv\n", nrow(somatic_trim)))

germline_trim <- load_run_metadata() |>
  build_metric_paths(pipeline = "germline", suffix = "trimmer_metrics.csv") |>
  agg_trimmer()
write_csv(germline_trim, file.path(OUT_DIR, "germline_trimmer_metrics.csv"))
cat(sprintf("Wrote %d rows to germline_trimmer_metrics.csv\n", nrow(germline_trim)))
