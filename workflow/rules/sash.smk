# sash aggregation -- per-table rules over the sash samplesheet. Each rule runs one
# ported agg_*.R via the shared r_cmd() helper (kind="sash"); the per-sample sash
# output dirs are resolved from the samplesheet at runtime by utils.R, so only the
# samplesheet is a declared input. Paths/env come from rules/common.smk.

# table -> agg script (1:1 except cancer_report -> 2 tables, handled separately).
SASH_SCRIPT = {
    "purple_drivers_somatic":  "agg_drivers.R",
    "purple_drivers_germline": "agg_drivers_germline.R",
    "purple_cnv_gene":         "agg_purple_cnv.R",
    "purple_cnv_somatic":      "agg_purple_cnv_somatic.R",
    "mutational_signatures":   "agg_mutational_signatures.R",
    "cpsr_germline_variants":  "agg_cpsr_germline.R",
}

rule sash_table:
    """One PURPLE/PCGR/CPSR/cohort table per ported agg script (1:1)."""
    input:  ss = SS["sash"]
    output: AGG("sash", "{table}")
    params: cmd = lambda wc: r_cmd("sash/" + SASH_SCRIPT[wc.table], kind="sash")
    conda:  ENV_R
    wildcard_constraints: table = "|".join(SASH_SCRIPT)
    shell:  "{params.cmd}"

rule sash_smlv:
    """Somatic small variants -- R agg that shells out to bcftools (needs htslib)."""
    input:  ss = SS["sash"]
    output: AGG("sash", "smlv_somatic_variants")
    params: cmd = r_cmd("sash/agg_smlv_somatic.R", kind="sash")
    conda:  ENV_R
    shell:  "{params.cmd}"

rule sash_cancer_report:
    """Cancer-report QC summary + signatures (one script, two tables)."""
    input:  ss = SS["sash"]
    output: qc  = AGG("sash", "cancer_report_qc_summary"),
            sig = AGG("sash", "cancer_report_signatures")
    params: cmd = r_cmd("sash/agg_cancer_report.R", kind="sash")
    conda:  ENV_R
    shell:  "{params.cmd}"

rule sash_cohort_summary:
    """Per-tumour cohort summary. Reads the mutational_signatures table (top-3 SBS
    column) from WGS_OUT_DIR at runtime, so declare it as an input to order this
    after sash_table(mutational_signatures) -- otherwise the column is dropped."""
    input:  ss   = SS["sash"],
            sigs = AGG("sash", "mutational_signatures")
    output: AGG("sash", "sash_cohort_summary")
    params: cmd = r_cmd("sash/cohort_summary.R", kind="sash")
    conda:  ENV_R
    shell:  "{params.cmd}"
