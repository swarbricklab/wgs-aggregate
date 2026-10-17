#!/usr/bin/env Rscript
# Aggregate WGS coverage metrics across all runs
# Input:  {tumor}.wgs_coverage_metrics_tumor.csv   (somatic run, tumor)
#         {tumor}.wgs_coverage_metrics_normal.csv  (somatic run, normal)
#         {germline}.wgs_coverage_metrics.csv       (germline run)
# Output: data/wgs/dragen_aggregate/somatic_coverage_metrics.csv
#         data/wgs/dragen_aggregate/germline_coverage_metrics.csv
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_coverage_metrics.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

clean_metric <- function(x) {
  x |> tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("_+$", "", x = _)
}

agg_coverage <- function(meta_tumor, meta_normal) {
  bind_rows(
    meta_tumor |>
      rowwise() |>
      reframe(
        pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
        read_dragen_metrics(metric_path) |>
          mutate(col = paste0("tumor_", clean_metric(metric)), value = as.numeric(value)) |>
          select(col, value)
      ),
    meta_normal |>
      rowwise() |>
      reframe(
        pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
        read_dragen_metrics(metric_path) |>
          mutate(col = paste0("normal_", clean_metric(metric)), value = as.numeric(value)) |>
          select(col, value)
      )
  ) |>
  pivot_wider(names_from = col, values_from = value)
}

meta <- load_run_metadata()

somatic_cov <- agg_coverage(
  meta |> build_metric_paths(pipeline = "somatic", suffix = "wgs_coverage_metrics_tumor.csv"),
  meta |> build_metric_paths(pipeline = "somatic", suffix = "wgs_coverage_metrics_normal.csv")
)
write_csv(somatic_cov, file.path(OUT_DIR, "somatic_coverage_metrics.csv"))
cat(sprintf("Wrote %d rows to somatic_coverage_metrics.csv\n", nrow(somatic_cov)))

germline_cov <- load_run_metadata() |>
  build_metric_paths(pipeline = "germline", suffix = "wgs_coverage_metrics.csv") |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_dragen_metrics(metric_path) |>
      mutate(col = clean_metric(metric), value = as.numeric(value)) |>
      select(col, value)
  ) |>
  pivot_wider(names_from = col, values_from = value)
write_csv(germline_cov, file.path(OUT_DIR, "germline_coverage_metrics.csv"))
cat(sprintf("Wrote %d rows to germline_coverage_metrics.csv\n", nrow(germline_cov)))
