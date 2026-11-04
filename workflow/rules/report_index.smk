# report_index -- searchable HTML/CSV indices of per-sample reports. Dispatched via
# the shared r_cmd() helper; sash + tumor_normal read identity from the sash
# samplesheet (WGS_SASH_SAMPLESHEET, always exported by r_cmd), the oncoanalyser
# scripts read their own samplesheet. Outputs live under reports_dir (fixed in the
# report_index R scripts; see config/config.yaml note).
# (No f-strings in Snakemake-parsed files -- see common.smk header.)

rule report_index_sash:
    input:  ss = SS["sash"]
    output: csv = REPORTS + "/sash.csv", html = REPORTS + "/sash.html"
    params: cmd = r_cmd("report_index/build_sash_report_index.R")
    conda:  ENV_R
    shell:  "{params.cmd}"

rule report_index_tumor_normal:
    input:  ss = SS["sash"], runs = SS["dragen_runs"]
    output: csv = REPORTS + "/tumor-normal.csv", html = REPORTS + "/tumor-normal.html"
    params: cmd = r_cmd("report_index/build_tumor_normal_report_index.R")
    conda:  ENV_R
    shell:  "{params.cmd}"

rule report_index_oncoanalyser:
    input:  ss = SS["oncoanalyser"]
    output: csv = REPORTS + "/oncoanalyser.csv", html = REPORTS + "/oncoanalyser.html"
    params: cmd = r_cmd("report_index/build_oncoanalyser_report_index.R")
    conda:  ENV_R
    shell:  "{params.cmd}"

rule report_index_oncoanalyser_dna_rna:
    input:  ss = SS["oncoanalyser_dna_rna"]
    output: csv = REPORTS + "/oncoanalyser_dna_rna.csv", html = REPORTS + "/oncoanalyser_dna_rna.html"
    params: cmd = r_cmd("report_index/build_oncoanalyser_dna_rna_report_index.R")
    conda:  ENV_R
    shell:  "{params.cmd}"
