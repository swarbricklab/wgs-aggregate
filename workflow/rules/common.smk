# Shared wiring for every rule file: config-derived paths, the single R-dispatch
# helper, and the per-step target lists. This replaces the five bespoke per-file
# env-builders (_sash_env / _onco_env / _dr_env / _dg_env / germline's inline env)
# with one `r_cmd()` so every R aggregation rule is dispatched the same way.
#
# NOTE: no f-strings anywhere in Snakemake-parsed files -- this snakemake + py3.12
# combo mangles f-strings in rule directives (spaces appear at every replacement
# field). Use "+" concatenation / str.format instead. (See the git history.)

import os

# --- config-derived paths (deps / outs / params already pulled in the Snakefile) --
DATA    = deps["data_dir"]
OUT     = outs["out_dir"]
REPORTS = outs["reports_dir"]
SS      = deps["samplesheets"]
STEPS   = set(params.get("steps", []))

SCRIPTS = os.path.join(workflow.basedir, "scripts")   # workflow/scripts
UTILS   = os.path.join(SCRIPTS, "utils.R")            # consolidated utils.R (WGS_UTILS)
# conda env path is resolved relative to the rule FILE (workflow/rules/), so hop
# up twice to repo-root env/ -- not once like a rule defined in workflow/Snakefile.
ENV_R   = "../../env/environment.yaml"                 # single conda env for all rules


# --- path helpers -------------------------------------------------------------
def AGG(step, table):
    """Cohort-table output path: <out_dir>/<step>_aggregate/<table>.csv."""
    return OUT + "/" + step + "_aggregate/" + table + ".csv"


# --- the single R dispatch helper ---------------------------------------------
# Every WGS_* variable the R scripts read is derivable from config. We export the
# full set from one place; `kind` adds the per-pipeline bits (sash/oncoanalyser/
# oncoanalyser_dna_rna), `out_dir` overrides the output dir (dragen). Scripts read
# only the vars they need and ignore the rest, so one env builder serves them all.
def wgs_env(kind=None, out_dir=None):
    e = {
        "WGS_UTILS":            UTILS,
        "WGS_SASH_SAMPLESHEET": SS["sash"],
        "WGS_DRAGEN_RUNS":      SS["dragen_runs"],
        "WGS_DRAGEN_DIR":       DATA + "/tumor-normal",
    }
    if kind:
        e["WGS_KIND"]         = kind
        e["WGS_SAMPLESHEET"]  = SS[kind]
        e["WGS_PIPELINE_DIR"] = DATA + "/" + kind
        e["WGS_OUT_DIR"]      = out_dir or (OUT + "/" + kind + "_aggregate")
    elif out_dir:
        e["WGS_OUT_DIR"] = out_dir
    return " ".join('{}="{}"'.format(k, v) for k, v in e.items())


def r_cmd(script_rel, kind=None, out_dir=None):
    """Shell command to run one aggregation R script under the resolved WGS_* env.

    `script_rel` is relative to workflow/scripts (e.g. "sash/agg_drivers.R").
    """
    return wgs_env(kind, out_dir) + " Rscript " + os.path.join(SCRIPTS, script_rel)


# --- per-step output tables (the rules themselves live in rules/<step>.smk) ----
SASH_TABLES = [
    "purple_drivers_somatic", "purple_drivers_germline", "purple_cnv_gene",
    "purple_cnv_somatic", "mutational_signatures", "cpsr_germline_variants",
    "smlv_somatic_variants", "cancer_report_qc_summary", "cancer_report_signatures",
    "sash_cohort_summary",
]
ONCO_TABLES = [
    "cuppa_predictions", "lilac_alleles", "chord_hrd", "virus_annotated",
    "linx_fusions", "linx_drivers", "linx_neoepitopes", "immune_escape",
]
DNARNA_TABLES = ["neo_neoepitopes", "neo_presented_peptides", "neo_tnb_summary"]
DRAGEN_TABLES = [
    "somatic_cnv_metrics", "somatic_coverage_metrics", "germline_coverage_metrics",
    "somatic_hrd_scores", "somatic_mapping_metrics", "germline_mapping_metrics",
    "ploidy_estimation_metrics", "germline_roh_hethom_metrics",
    "somatic_sv_metrics", "somatic_trimmer_metrics", "germline_trimmer_metrics",
    "somatic_vc_metrics", "germline_vc_metrics",
]
REPORT_INDICES = ["sash", "tumor-normal", "oncoanalyser", "oncoanalyser_dna_rna"]
GERMLINE_VCF   = OUT + "/germline_cohort_vcf/cohort.germline.common_snps.vcf.gz"


def targets():
    """Final outputs for `rule all`, gated by params.steps."""
    t = []
    if "sash" in STEPS:                 t += [AGG("sash", x)                 for x in SASH_TABLES]
    if "oncoanalyser" in STEPS:         t += [AGG("oncoanalyser", x)         for x in ONCO_TABLES]
    if "oncoanalyser_dna_rna" in STEPS: t += [AGG("oncoanalyser_dna_rna", x) for x in DNARNA_TABLES]
    if "dragen" in STEPS:               t += [AGG("dragen", x)               for x in DRAGEN_TABLES]
    if "report_index" in STEPS:         t += [REPORTS + "/" + x + ".csv"     for x in REPORT_INDICES]
    if "germline" in STEPS:             t += [GERMLINE_VCF]
    return t
