# sash aggregation — per-table rules over config/sash/samplesheet.csv.
# Each rule runs one ported agg_*.R against the consolidated utils.R (WGS_UTILS).
# The per-sample sash output dirs (resolved from the samplesheet's sash_dir col)
# are read at runtime, so only the samplesheet is a declared input.

_SASH_SS  = SS["sash"]
_SASH_OUT = OUT + "/sash_aggregate"
_SASH_DIR = SCRIPTS + "/sash"

def _sash_env():
    return (
        'WGS_KIND=sash WGS_UTILS="' + UTILS + '" WGS_SAMPLESHEET="' + _SASH_SS + '" '
        'WGS_PIPELINE_DIR="' + DATA + '/sash" WGS_OUT_DIR="' + _SASH_OUT + '"'
    )

# table -> agg script (1:1 except cancer_report -> 2 tables, handled separately)
_SASH_SCRIPT = {
    "purple_drivers_somatic":  "agg_drivers.R",
    "purple_drivers_germline": "agg_drivers_germline.R",
    "purple_cnv_gene":         "agg_purple_cnv.R",
    "purple_cnv_somatic":      "agg_purple_cnv_somatic.R",
    "mutational_signatures":   "agg_mutational_signatures.R",
    "cpsr_germline_variants":  "agg_cpsr_germline.R",
    "sash_cohort_summary":     "cohort_summary.R",
}

rule sash_table:
    """One PURPLE/PCGR/CPSR/cohort table per ported agg script (1:1)."""
    input:  ss = _SASH_SS
    output: _SASH_OUT + "/{table}.csv"
    params: cmd = lambda wc: _sash_env() + " Rscript " + _SASH_DIR + "/" + _SASH_SCRIPT[wc.table]
    conda:  "../envs/r-aggregate.yaml"
    wildcard_constraints: table = "|".join(_SASH_SCRIPT)
    shell:  "{params.cmd}"

rule sash_smlv:
    """Somatic small variants — R agg that shells out to bcftools (needs htslib)."""
    input:  ss = _SASH_SS
    output: _SASH_OUT + "/smlv_somatic_variants.csv"
    params: cmd = _sash_env() + " Rscript " + _SASH_DIR + "/agg_smlv_somatic.R"
    conda:  "../envs/r-aggregate.yaml"
    shell:  "{params.cmd}"

rule sash_cancer_report:
    """Cancer-report QC summary + signatures (one script, two tables)."""
    input:  ss = _SASH_SS
    output: qc  = _SASH_OUT + "/cancer_report_qc_summary.csv",
            sig = _SASH_OUT + "/cancer_report_signatures.csv"
    params: cmd = _sash_env() + " Rscript " + _SASH_DIR + "/agg_cancer_report.R"
    conda:  "../envs/r-aggregate.yaml"
    shell:  "{params.cmd}"
