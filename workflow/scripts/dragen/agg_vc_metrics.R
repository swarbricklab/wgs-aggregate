#!/usr/bin/env Rscript
# Aggregate variant calling metrics across all runs
# Input:  {tumor}.vc_metrics.csv    (somatic: PREFILTER + POSTFILTER)
#         {germline}.vc_metrics.csv (germline: PREFILTER + POSTFILTER)
# Output: data/wgs/dragen_aggregate/somatic_vc_metrics.csv
#         data/wgs/dragen_aggregate/germline_vc_metrics.csv
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_vc_metrics.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

clean_metric <- function(x) {
  x |> tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("_+$", "", x = _)
}

agg_vc <- function(meta) {
  meta |>
    rowwise() |>
    reframe(
      pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
      read_dragen_metrics(metric_path) |>
        filter(section %in% c(
          "VARIANT CALLER SUMMARY",
          "VARIANT CALLER PREFILTER",
          "VARIANT CALLER POSTFILTER"
        )) |>
        mutate(
          filter_level = section |>
            sub("VARIANT CALLER ", "", x = _) |>
            tolower() |>
            sub("summary", "summary", x = _),
          col   = paste0(filter_level, "_", clean_metric(metric)),
          value = suppressWarnings(as.numeric(value))
        ) |>
        select(col, value)
    ) |>
    pivot_wider(names_from = col, values_from = value)
}

somatic_vc <- load_run_metadata() |>
  build_metric_paths(pipeline = "somatic", suffix = "vc_metrics.csv") |>
  agg_vc()
write_csv(somatic_vc, file.path(OUT_DIR, "somatic_vc_metrics.csv"))
cat(sprintf("Wrote %d rows to somatic_vc_metrics.csv\n", nrow(somatic_vc)))

germline_vc <- load_run_metadata() |>
  build_metric_paths(pipeline = "germline", suffix = "vc_metrics.csv") |>
  agg_vc()
write_csv(germline_vc, file.path(OUT_DIR, "germline_vc_metrics.csv"))
cat(sprintf("Wrote %d rows to germline_vc_metrics.csv\n", nrow(germline_vc)))
