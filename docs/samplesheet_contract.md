# Samplesheet contract

The workflow consumes the dataset's **consolidated** samplesheets (one row per
tumour/run). No registry files (`index.csv`, `config/projects/*`) are read — all
sample identity comes from these columns.

## Common identity columns (all pipelines)
`tumor_gussid, donor_id, tumor_tissue_id, tumor_lib, normal_lib`

## Per-pipeline keys
- **sash** (`samplesheets.sash`): `sash_run_id`, `sash_dir` (per-sample output dir), optional `skip_reason`.
- **oncoanalyser** (`samplesheets.oncoanalyser`): `onco_run_id` (dir resolved as `<data_dir>/oncoanalyser/<onco_run_id>/<tumor_lib><sep><normal_lib>`, sep `_`|`__`).
- **oncoanalyser_dna_rna** (`samplesheets.oncoanalyser_dna_rna`): `rna_lib, dna_rna_run_id, lib_dir`.
- **dragen** (`samplesheets.dragen_runs`): `donor_id, tumor_library, germline_library, run_id`; `tissue_id` is joined from the sash samplesheet (not a registry index).

## Germline cohort (`cohort_samplesheet`)
`donor_id, assay, usage` — the germline cohort is `assay == wgs` and
`usage ∈ germline_usage_include`. Replaces the old `config/projects/wgs.csv`.

## Dropped vs the registry
`wgs_sbj` / `wgs_prj` (Garvan↔UMCCR id links) are **not** carried — a registry
concern, not a dataset one.
