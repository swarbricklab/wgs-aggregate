#!/usr/bin/env Rscript
# Aggregate CNV metrics across all somatic runs
# Input:  {tumor}.cnv_metrics.csv (DRAGEN format, 16 rows)
# Output: data/wgs/dragen_aggregate/somatic_cnv_metrics.csv
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_cnv_metrics.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

meta <- load_run_metadata() |>
  build_metric_paths(pipeline = "somatic", suffix = "cnv_metrics.csv")

cnv <- meta |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_dragen_metrics(metric_path) |>
      select(metric, value)
  ) |>
  mutate(
    metric = metric |>
      tolower() |>
      gsub("[\\s()]+", "_", x = _, perl = TRUE) |>
      gsub("_+$", "", x = _) |>
      gsub("_+", "_", x = _),
    value = as.numeric(value)
  ) |>
  pivot_wider(names_from = metric, values_from = value)

write_csv(cnv, file.path(OUT_DIR, "somatic_cnv_metrics.csv"))
cat(sprintf("Wrote %d rows to somatic_cnv_metrics.csv\n", nrow(cnv)))
