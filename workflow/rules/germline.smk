# germline cohort VCF -- per-donor DRAGEN germline VCFs -> cohort common-SNP panel
# + per-project subsets. Decoupled from index.csv (germline libs from the sash
# samplesheet) and config/projects/wgs.csv (cohort donors from the cohort
# samplesheet). Heavy (bcftools merge across ~all donors); run on a compute node.
#
# Unlike the R aggregations this is a bash pipeline (scripts/germline/run.sh), so it
# keeps its own env block rather than r_cmd(); all values still come from config,
# including the panel knobs (maf_min / autosomes) that are now honoured.
# (No f-strings in Snakemake-parsed files -- see common.smk header.)

GL_OUT = OUT + "/germline_cohort_vcf"
GL_DIR = os.path.join(SCRIPTS, "germline")
_gl    = params["germline"]

rule germline_cohort_vcf:
    input:
        sash_ss   = SS["sash"],
        cohort_ss = deps["cohort_samplesheet"],
    output:
        vcf      = GL_OUT + "/cohort.germline.common_snps.vcf.gz",
        manifest = GL_OUT + "/donor_manifest.csv",
        projects = GL_OUT + "/projects/project_donors.csv",
    params:
        env = " ".join([
            'REPO="$PWD"',
            'GERMLINE_SCRIPT_DIR="' + GL_DIR + '"',
            'WGS_DRAGEN_RUNS="' + SS["dragen_runs"] + '"',
            'WGS_SASH_SAMPLESHEET="' + SS["sash"] + '"',
            'WGS_COHORT_SAMPLESHEET="' + deps["cohort_samplesheet"] + '"',
            'WGS_GERMLINE_USAGE="' + ",".join(_gl["usage_include"]) + '"',
            'WGS_MAF_MIN="' + str(_gl["maf_min"]) + '"',
            'WGS_AUTOSOMES="' + str(_gl["autosomes"]) + '"',
        ])
    threads: 16
    conda: ENV_R
    shell: "{params.env} THREADS={threads} bash " + GL_DIR + "/run.sh"
