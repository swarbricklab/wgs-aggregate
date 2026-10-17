# oncoanalyser DNA+RNA aggregation — per-table rules over
# config/oncoanalyser_dna_rna/samplesheet.csv (dna_rna_run_id + lib_dir -> per-sample dir).

_DR_SS  = SS["oncoanalyser_dna_rna"]
_DR_OUT = OUT + "/oncoanalyser_dna_rna_aggregate"
_DR_DIR = SCRIPTS + "/oncoanalyser_dna_rna"

def _dr_env():
    return (
        'WGS_KIND=oncoanalyser_dna_rna WGS_UTILS="' + UTILS + '" WGS_SAMPLESHEET="' + _DR_SS + '" '
        'WGS_PIPELINE_DIR="' + DATA + '/oncoanalyser_dna_rna" WGS_OUT_DIR="' + _DR_OUT + '"'
    )

rule dna_rna_neoepitopes:
    """HMF Neo scored neoepitopes."""
    input:  ss = _DR_SS
    output: _DR_OUT + "/neo_neoepitopes.csv"
    params: cmd = _dr_env() + " Rscript " + _DR_DIR + "/agg_neo_neoepitopes.R"
    conda:  "../envs/r-aggregate.yaml"
    shell:  "{params.cmd}"

rule dna_rna_tnb:
    """Presented peptides + TNB summary (one script, two tables)."""
    input:  ss = _DR_SS
    output: pep = _DR_OUT + "/neo_presented_peptides.csv",
            tnb = _DR_OUT + "/neo_tnb_summary.csv"
    params: cmd = _dr_env() + " Rscript " + _DR_DIR + "/agg_neo_tnb.R"
    conda:  "../envs/r-aggregate.yaml"
    shell:  "{params.cmd}"
