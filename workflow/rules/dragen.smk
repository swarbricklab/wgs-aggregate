# DRAGEN metrics aggregation — one rule per agg script (each emits 1-2 tables).
# Decoupled from index.csv: load_run_metadata() (in utils.R) joins dragen_runs.csv
# with the sash samplesheet for tissue_id and DROPS the registry-only wgs_sbj/wgs_prj.

_DG_OUT = OUT + "/dragen_aggregate"
_DG_DIR = SCRIPTS + "/dragen"

def _dg_env():
    return (
        'WGS_UTILS="' + UTILS + '" WGS_OUT_DIR="' + _DG_OUT + '" '
        'WGS_DRAGEN_RUNS="' + SS["dragen_runs"] + '" '
        'WGS_SASH_SAMPLESHEET="' + SS["sash"] + '" '
        'WGS_DRAGEN_DIR="' + DATA + '/tumor-normal"'
    )

# agg script -> output table(s)
_DG = {
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

for _script, _tables in _DG.items():
    rule:
        name:   "dragen_" + _script.replace("agg_", "").replace(".R", "")
        input:  runs = SS["dragen_runs"], sash_ss = SS["sash"]
        output: [_DG_OUT + "/" + t + ".csv" for t in _tables]
        params: cmd = _dg_env() + " Rscript " + _DG_DIR + "/" + _script
        conda:  "../envs/r-aggregate.yaml"
        shell:  "{params.cmd}"
