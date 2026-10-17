# germline cohort VCF — per-donor DRAGEN germline VCFs -> cohort common-SNP panel
# + atlas-cohort subset. Decoupled from index.csv (germline libs from the sash
# samplesheet) and config/projects/wgs.csv (cohort donors from the atlas samplesheet).
# Heavy (bcftools merge across ~all donors); run on a compute node.

_GL_OUT = OUT + "/germline_cohort_vcf"
_GL_DIR = SCRIPTS + "/germline"
_USAGE  = ",".join(config.get("germline_usage_include", ["analysed", "warn_qc"]))

rule germline_cohort_vcf:
    input:
        sash_ss   = SS["sash"],
        cohort_ss = config["cohort_samplesheet"],
    output:
        vcf      = _GL_OUT + "/cohort.germline.common_snps.vcf.gz",
        manifest = _GL_OUT + "/donor_manifest.csv",
        projects = _GL_OUT + "/projects/project_donors.csv",
    params:
        env = ('REPO="$PWD" GERMLINE_SCRIPT_DIR="' + _GL_DIR + '" '
               'WGS_DRAGEN_RUNS="' + SS["dragen_runs"] + '" '
               'WGS_SASH_SAMPLESHEET="' + SS["sash"] + '" '
               'WGS_COHORT_SAMPLESHEET="' + config["cohort_samplesheet"] + '" '
               'WGS_GERMLINE_USAGE="' + _USAGE + '"')
    threads: 16
    conda: "../envs/bcftools.yaml"
    shell: "{params.env} THREADS={threads} bash " + _GL_DIR + "/run.sh"
