# oncoanalyser (DNA-only) aggregation -- per-table rules over the oncoanalyser
# samplesheet. onco_run_id -> per-sample dir is resolved in utils.R. One agg_*.R
# per table, dispatched via the shared r_cmd() helper (kind="oncoanalyser").

ONCO_SCRIPT = {
    "cuppa_predictions": "agg_cuppa.R",
    "lilac_alleles":     "agg_lilac.R",
    "chord_hrd":         "agg_chord.R",
    "virus_annotated":   "agg_virus.R",
    "linx_fusions":      "agg_fusions.R",
    "linx_drivers":      "agg_linx_drivers.R",
    "linx_neoepitopes":  "agg_neoepitopes.R",
    "immune_escape":     "agg_immune_escape.R",
}

rule oncoanalyser_table:
    input:  ss = SS["oncoanalyser"]
    output: AGG("oncoanalyser", "{table}")
    params: cmd = lambda wc: r_cmd("oncoanalyser/" + ONCO_SCRIPT[wc.table], kind="oncoanalyser")
    conda:  ENV_R
    wildcard_constraints: table = "|".join(ONCO_SCRIPT)
    shell:  "{params.cmd}"
