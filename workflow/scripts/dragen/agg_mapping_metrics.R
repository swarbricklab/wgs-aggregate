#!/usr/bin/env Rscript
# Aggregate mapping metrics across all runs
# Input:  {tumor}.mapping_metrics.csv  (somatic run: TUMOR + NORMAL SUMMARY sections)
#         {germline}.mapping_metrics.csv (germline run: MAPPING/ALIGNING SUMMARY)
# Output: data/wgs/dragen_aggregate/somatic_mapping_metrics.csv
#         data/wgs/dragen_aggregate/germline_mapping_metrics.csv
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_mapping_metrics.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

clean_metric <- function(x) {
  x |> tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("_+$", "", x = _)
}

# --- Somatic: extract TUMOR and NORMAL SUMMARY rows, output wide with prefix ---
somatic_meta <- load_run_metadata() |>
  build_metric_paths(pipeline = "somatic", suffix = "mapping_metrics.csv")

somatic_map <- somatic_meta |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_dragen_metrics(metric_path) |>
      filter(section %in% c(
        "TUMOR MAPPING/ALIGNING SUMMARY",
        "NORMAL MAPPING/ALIGNING SUMMARY"
      )) |>
      mutate(
        prefix  = if_else(grepl("^TUMOR", section), "tumor", "normal"),
        col     = paste0(prefix, "_", clean_metric(metric)),
        value   = as.numeric(value)
      ) |>
      select(col, value)
  ) |>
  pivot_wider(names_from = col, values_from = value)

write_csv(somatic_map, file.path(OUT_DIR, "somatic_mapping_metrics.csv"))
cat(sprintf("Wrote %d rows to somatic_mapping_metrics.csv\n", nrow(somatic_map)))

# --- Germline: extract MAPPING/ALIGNING SUMMARY rows ---
germline_meta <- load_run_metadata() |>
  build_metric_paths(pipeline = "germline", suffix = "mapping_metrics.csv")

germline_map <- germline_meta |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_dragen_metrics(metric_path) |>
      filter(section == "MAPPING/ALIGNING SUMMARY") |>
      mutate(
        col   = clean_metric(metric),
        value = as.numeric(value)
      ) |>
      select(col, value)
  ) |>
  pivot_wider(names_from = col, values_from = value)

write_csv(germline_map, file.path(OUT_DIR, "germline_mapping_metrics.csv"))
cat(sprintf("Wrote %d rows to germline_mapping_metrics.csv\n", nrow(germline_map)))
