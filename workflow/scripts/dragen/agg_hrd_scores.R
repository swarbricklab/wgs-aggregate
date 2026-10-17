#!/usr/bin/env Rscript
# Aggregate HRD scores across all somatic runs
# Input:  {tumor}.hrdscore.csv (standard CSV with header)
# Output: data/wgs/dragen_aggregate/somatic_hrd_scores.csv
#
# Run from project root:
#   Rscript modules/dragen_aggregate/scripts/agg_hrd_scores.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

meta <- load_run_metadata() |>
  build_metric_paths(pipeline = "somatic", suffix = "hrdscore.csv")

hrd <- meta |>
  rowwise() |>
  reframe(
    pick(run_id, donor_id, tissue_id, tumor_library, germline_library),
    read_csv(metric_path, col_types = "cccc") |>
      rename(
        loh_score = LOH_Score,
        tai_score = TAI_Score,
        lst_score = LST_Score,
        hrd_score = HRD_Score
      ) |>
      select(-Sample)
  ) |>
  mutate(across(ends_with("_score"), as.integer))

write_csv(hrd, file.path(OUT_DIR, "somatic_hrd_scores.csv"))
cat(sprintf("Wrote %d rows to somatic_hrd_scores.csv\n", nrow(hrd)))
