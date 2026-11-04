# oncoanalyser DNA+RNA aggregation -- per-table rules over the oncoanalyser_dna_rna
# samplesheet (dna_rna_run_id + lib_dir -> per-sample dir, resolved in utils.R).
# Dispatched via the shared r_cmd() helper (kind="oncoanalyser_dna_rna").

rule dna_rna_neoepitopes:
    """HMF Neo scored neoepitopes."""
    input:  ss = SS["oncoanalyser_dna_rna"]
    output: AGG("oncoanalyser_dna_rna", "neo_neoepitopes")
    params: cmd = r_cmd("oncoanalyser_dna_rna/agg_neo_neoepitopes.R", kind="oncoanalyser_dna_rna")
    conda:  ENV_R
    shell:  "{params.cmd}"

rule dna_rna_tnb:
    """Presented peptides + TNB summary (one script, two tables)."""
    input:  ss = SS["oncoanalyser_dna_rna"]
    output: pep = AGG("oncoanalyser_dna_rna", "neo_presented_peptides"),
            tnb = AGG("oncoanalyser_dna_rna", "neo_tnb_summary")
    params: cmd = r_cmd("oncoanalyser_dna_rna/agg_neo_tnb.R", kind="oncoanalyser_dna_rna")
    conda:  ENV_R
    shell:  "{params.cmd}"
