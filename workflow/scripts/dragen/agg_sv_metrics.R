#!/usr/bin/env Rscript
# Aggregate SV metrics across all somatic runs
# Input:  {tumor}.sv_metrics.csv (DRAGEN format, 5 rows)
# Output: data/wgs/dragen_aggregate/somatic_sv_metrics.csv
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_sv_metrics.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

meta <- load_run_metadata() |>
  build_metric_paths(pipeline = "somatic", suffix = "sv_metrics.csv")

sv <- meta |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_dragen_metrics(metric_path) |>
      select(metric, value)
  ) |>
  mutate(
    metric = metric |>
      tolower() |>
      gsub("\\s*\\(pass\\)", "", x = _) |>
      gsub(" ", "_", x = _),
    value = as.integer(value)
  ) |>
  pivot_wider(names_from = metric, values_from = value)

write_csv(sv, file.path(OUT_DIR, "somatic_sv_metrics.csv"))
cat(sprintf("Wrote %d rows to somatic_sv_metrics.csv\n", nrow(sv)))
