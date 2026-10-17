# oncoanalyser (DNA-only) aggregation — per-table rules over
# config/oncoanalyser/samplesheet.csv. onco_run_id -> per-sample dir resolved in utils.R.

_ONCO_SS  = SS["oncoanalyser"]
_ONCO_OUT = OUT + "/oncoanalyser_aggregate"
_ONCO_DIR = SCRIPTS + "/oncoanalyser"

def _onco_env():
    return (
        'WGS_KIND=oncoanalyser WGS_UTILS="' + UTILS + '" WGS_SAMPLESHEET="' + _ONCO_SS + '" '
        'WGS_PIPELINE_DIR="' + DATA + '/oncoanalyser" WGS_OUT_DIR="' + _ONCO_OUT + '"'
    )

_ONCO_SCRIPT = {
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
    input:  ss = _ONCO_SS
    output: _ONCO_OUT + "/{table}.csv"
    params: cmd = lambda wc: _onco_env() + " Rscript " + _ONCO_DIR + "/" + _ONCO_SCRIPT[wc.table]
    conda:  "../envs/r-aggregate.yaml"
    wildcard_constraints: table = "|".join(_ONCO_SCRIPT)
    shell:  "{params.cmd}"
