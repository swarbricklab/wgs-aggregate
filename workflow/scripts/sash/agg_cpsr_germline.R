#!/usr/bin/env Rscript
# Aggregate per-tumour CPSR germline variant classifications into a cohort CSV.
# Input:  <src_dir>/smlv_germline/report/cpsr/<normal_lib>.cpsr.grch38.snvs_indels.tiers.tsv
# Output: data/wgs/sash_aggregate/cpsr_germline_variants.csv
#
# Unlike the other aggregators, the source file is prefixed with *normal_lib*
# (CPSR classifies the germline sample). Rows are still keyed per samplesheet
# row (tumor_gussid), so donors with >1 tumour repeat their germline variants
# once per tumour — dedupe on donor_id (or normal_lib) for donor-level
# carrier counts.
#
# All variants are kept (Benign through Pathogenic; ~7.5k rows per sample) —
# no thresholding, downstream filters on final_classification. Dropped as
# redundant/verbose: VEP_ALL_CSQ (full VEP dump; ~40% of file width),
# CPSR_CLASSIFICATION_DOC (boilerplate criterion text), the 31 per-criterion
# ACMG_* booleans (cpsr_classification_code carries the fired criteria,
# pipe-joined), GENOME_VERSION (constant grch38), VCF_SAMPLE_ID (= normal_lib).
#
# Run from project root:
#   Rscript modules/sash_aggregate/scripts/agg_cpsr_germline.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

NUMERIC_COLS <- c("exon_position", "intron_position", "cancer_phenotype",
                  "clinvar_review_status_stars", "gerp_score",
                  "n_insilico_called", "n_insilico_damaging",
                  "n_insilico_tolerated", "n_insilico_splicing_neutral",
                  "n_insilico_splicing_affected", "global_af_gnomad",
                  "non_cancer_af_global", "cpsr_pathogenicity_score")
LOGICAL_COLS <- c("oncogene", "tumor_suppressor", "last_exon",
                  "loss_of_function", "clinvar_conflicted")

# UPPER_SNAKE source -> snake_case output. entrez_id / clinvar_msid stay
# character: they are identifiers, not quantities.
COL_MAP <- c(
  GENOMIC_CHANGE               = "genomic_change",
  VAR_ID                       = "var_id",
  GENOTYPE                     = "genotype",
  CPSR_CLASSIFICATION_SOURCE   = "cpsr_classification_source",
  VARIANT_CLASS                = "variant_class",
  CODING_STATUS                = "coding_status",
  SYMBOL                       = "symbol",
  GENE_NAME                    = "gene_name",
  CCDS                         = "ccds",
  ENTREZ_ID                    = "entrez_id",
  UNIPROT_ID                   = "uniprot_id",
  ENSEMBL_GENE_ID              = "ensembl_gene_id",
  ENSEMBL_TRANSCRIPT_ID        = "ensembl_transcript_id",
  REFSEQ_MRNA                  = "refseq_mrna",
  ONCOGENE                     = "oncogene",
  TUMOR_SUPPRESSOR             = "tumor_suppressor",
  CONSEQUENCE                  = "consequence",
  REGULATORY_ANNOTATION        = "regulatory_annotation",
  PROTEIN_CHANGE               = "protein_change",
  PROTEIN_DOMAIN               = "protein_domain",
  DBSNP                        = "dbsnp",
  HGVSp                        = "hgvsp",
  HGVSc                        = "hgvsc",
  LAST_EXON                    = "last_exon",
  EXON_POSITION                = "exon_position",
  INTRON_POSITION              = "intron_position",
  CDS_CHANGE                   = "cds_change",
  MUTATION_HOTSPOT             = "mutation_hotspot",
  RMSK_HIT                     = "rmsk_hit",
  PROTEIN_FEATURE              = "protein_feature",
  EFFECT_PREDICTIONS           = "effect_predictions",
  LOSS_OF_FUNCTION             = "loss_of_function",
  CANCER_PHENOTYPE             = "cancer_phenotype",
  CLINVAR_CLASSIFICATION       = "clinvar_classification",
  CLINVAR_MSID                 = "clinvar_msid",
  CLINVAR_VARIANT_ORIGIN       = "clinvar_variant_origin",
  CLINVAR_CONFLICTED           = "clinvar_conflicted",
  CLINVAR_PHENOTYPE            = "clinvar_phenotype",
  CLINVAR_REVIEW_STATUS_STARS  = "clinvar_review_status_stars",
  DBMTS                        = "dbmts",
  miRNA_TARGET_HIT             = "mirna_target_hit",
  miRNA_TARGET_HIT_PREDICTION  = "mirna_target_hit_prediction",
  TF_BINDING_SITE_VARIANT      = "tf_binding_site_variant",
  TF_BINDING_SITE_VARIANT_INFO = "tf_binding_site_variant_info",
  GERP_SCORE                   = "gerp_score",
  N_INSILICO_CALLED            = "n_insilico_called",
  N_INSILICO_DAMAGING          = "n_insilico_damaging",
  N_INSILICO_TOLERATED         = "n_insilico_tolerated",
  N_INSILICO_SPLICING_NEUTRAL  = "n_insilico_splicing_neutral",
  N_INSILICO_SPLICING_AFFECTED = "n_insilico_splicing_affected",
  GLOBAL_AF_GNOMAD             = "global_af_gnomad",
  NON_CANCER_AF_GLOBAL         = "non_cancer_af_global",
  FINAL_CLASSIFICATION         = "final_classification",
  CPSR_CLASSIFICATION          = "cpsr_classification",
  CPSR_PATHOGENICITY_SCORE     = "cpsr_pathogenicity_score",
  CPSR_CLASSIFICATION_CODE     = "cpsr_classification_code"
)

cpsr_tiers_path <- function(src_dir, normal_lib) {
  path(src_dir, "smlv_germline", "report", "cpsr",
       paste0(normal_lib, ".cpsr.grch38.snvs_indels.tiers.tsv"))
}

read_cpsr_tsv <- function(src_dir, normal_lib) {
  p <- cpsr_tiers_path(src_dir, normal_lib)
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())

  present <- intersect(names(COL_MAP), names(df))
  df <- df[, present, drop = FALSE]
  names(df) <- unname(COL_MAP[present])
  df
}

missing <- ss %>% filter(!file.exists(cpsr_tiers_path(src_dir, normal_lib)))
if (nrow(missing) > 0) {
  message(sprintf("WARN: %d tumour(s) with no CPSR tiers file (0 rows contributed):",
                  nrow(missing)))
  print(missing %>% select(any_of(c("tumor_gussid", "sash_run_id"))))
}

rows <- ss %>%
  rowwise() %>%
  reframe(
    pick(any_of(IDENTITY_COLS)),
    read_cpsr_tsv(src_dir, normal_lib)
  )

if (nrow(rows) > 0) {
  num_present <- intersect(NUMERIC_COLS, names(rows))
  rows <- rows %>%
    mutate(across(all_of(num_present), ~ suppressWarnings(as.numeric(.))))

  log_present <- intersect(LOGICAL_COLS, names(rows))
  rows <- rows %>%
    mutate(across(all_of(log_present), ~ tolower(.) == "true"))
}

dir_create(OUT_DIR)
out_csv <- file.path(OUT_DIR, "cpsr_germline_variants.csv")
write_csv(rows, out_csv, na = "NA")
cat(sprintf("Wrote %d variant row(s) from %d tumour(s) to %s\n",
            nrow(rows), length(unique(rows$tumor_gussid)), out_csv))
