# report_index — searchable HTML/CSV indices of per-sample reports.
# sash + tumor_normal get identity from the sash samplesheet (WGS_SASH_SAMPLESHEET);
# oncoanalyser + dna_rna read their own samplesheets (already registry-free).
# Outputs live under reports/wgs/ (fixed location, not out_dir).

_RI_DIR = SCRIPTS + "/report_index"
_SASH_ENV = 'WGS_SASH_SAMPLESHEET="' + SS["sash"] + '"'

rule report_index_sash:
    input:  ss = SS["sash"]
    output: csv = "reports/wgs/sash.csv", html = "reports/wgs/sash.html"
    params: cmd = _SASH_ENV + " Rscript " + _RI_DIR + "/build_sash_report_index.R"
    conda:  "../envs/r-aggregate.yaml"
    shell:  "{params.cmd}"

rule report_index_tumor_normal:
    input:  ss = SS["sash"], runs = SS["dragen_runs"]
    output: csv = "reports/wgs/tumor-normal.csv", html = "reports/wgs/tumor-normal.html"
    params: cmd = _SASH_ENV + " Rscript " + _RI_DIR + "/build_tumor_normal_report_index.R"
    conda:  "../envs/r-aggregate.yaml"
    shell:  "{params.cmd}"

rule report_index_oncoanalyser:
    input:  ss = SS["oncoanalyser"]
    output: csv = "reports/wgs/oncoanalyser.csv", html = "reports/wgs/oncoanalyser.html"
    params: cmd = "Rscript " + _RI_DIR + "/build_oncoanalyser_report_index.R"
    conda:  "../envs/r-aggregate.yaml"
    shell:  "{params.cmd}"

rule report_index_oncoanalyser_dna_rna:
    input:  ss = SS["oncoanalyser_dna_rna"]
    output: csv = "reports/wgs/oncoanalyser_dna_rna.csv", html = "reports/wgs/oncoanalyser_dna_rna.html"
    params: cmd = "Rscript " + _RI_DIR + "/build_oncoanalyser_dna_rna_report_index.R"
    conda:  "../envs/r-aggregate.yaml"
    shell:  "{params.cmd}"
