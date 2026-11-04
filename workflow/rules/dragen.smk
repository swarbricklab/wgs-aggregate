# DRAGEN metrics aggregation -- one rule per agg script (each emits 1-2 tables).
# Decoupled from index.csv: load_run_metadata() (in utils.R) joins dragen_runs.csv
# with the sash samplesheet for tissue_id. Dispatched via the shared r_cmd() helper
# (no `kind`; dragen reads WGS_DRAGEN_* + WGS_OUT_DIR, which r_cmd supplies).

DRAGEN_OUT = OUT + "/dragen_aggregate"

# agg script -> output table(s)
DRAGEN_SCRIPT = {
    "agg_cnv_metrics.R":        ["somatic_cnv_metrics"],
    "agg_coverage_metrics.R":   ["somatic_coverage_metrics", "germline_coverage_metrics"],
    "agg_hrd_scores.R":         ["somatic_hrd_scores"],
    "agg_mapping_metrics.R":    ["somatic_mapping_metrics", "germline_mapping_metrics"],
    "agg_ploidy_metrics.R":     ["ploidy_estimation_metrics"],
    "agg_roh_hethom_metrics.R": ["germline_roh_hethom_metrics"],
    "agg_sv_metrics.R":         ["somatic_sv_metrics"],
    "agg_trimmer_metrics.R":    ["somatic_trimmer_metrics", "germline_trimmer_metrics"],
    "agg_vc_metrics.R":         ["somatic_vc_metrics", "germline_vc_metrics"],
}

for _script, _tables in DRAGEN_SCRIPT.items():
    rule:
        name:   "dragen_" + _script.replace("agg_", "").replace(".R", "")
        input:  runs = SS["dragen_runs"], sash_ss = SS["sash"]
        output: [DRAGEN_OUT + "/" + t + ".csv" for t in _tables]
        params: cmd = r_cmd("dragen/" + _script, out_dir=DRAGEN_OUT)
        conda:  ENV_R
        shell:  "{params.cmd}"
