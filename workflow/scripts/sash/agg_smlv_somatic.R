#!/usr/bin/env Rscript
# Aggregate per-tumour sash somatic small variants (SNV/indel) into a cohort CSV.
#
# Input:  <src_dir>/smlv_somatic/filter/<tumor_lib>.pass.vcf.gz          (row set)
#         <src_dir>/smlv_somatic/report/pcgr/<tumor_lib>.pcgr_acmg.grch38.snvs_indels.tiers.tsv
# Output: data/wgs/sash_aggregate/smlv_somatic_variants.csv
#
# Long form — one row per (tumour, PASS somatic small variant), ~4.3k rows per
# tumour, ~895k rows for the cohort.
#
# INCLUSION CRITERIA: every record in the sash PASS VCF, nothing else. This is
# the same row set as results/BrCa_combine_sash_filter.maf, which was produced
# by running vcf2maf over these same per-sample VCFs and concatenating; per-
# sample record counts agree exactly for all 208 tumours. No thresholding is
# applied here either — filter on pcgr_tier / vep_gnomad_af / t_af downstream.
#
# The VCF supplies the row set, the VEP annotation (INFO/PCGR_CSQ, picked
# transcript) and the read support; the PCGR tiers TSV is left-joined on
# chrom/pos/ref/alt for the clinical layer (tier, hotspots, COSMIC, TCGA,
# ChEMBL, OpenTargets). The two files are 1:1 across the cohort.
#
# A variant usually overlaps several transcripts. The unprefixed annotation
# columns report the one chosen by extract_smlv_somatic.sh's canonical-first
# order (see its header for the derivation); PCGR's tiers TSV picks its own,
# differing on ~46% of variants, so every transcript-dependent column taken
# from that file is prefixed `pcgr_` and the two are never silently mixed.
#
# Requires bcftools on PATH:
#   module load bcftools/1.21
#
# Run from project root:
#   Rscript modules/sash_aggregate/scripts/agg_smlv_somatic.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

# Co-located with this script (workflow/scripts/sash/); resolve it relative to
# WGS_UTILS (workflow/scripts/utils.R) rather than a hardcoded inline-module path.
EXTRACT_SH <- file.path(dirname(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R")),
                        "sash", "extract_smlv_somatic.sh")
OUT_CSV    <- file.path(OUT_DIR, "smlv_somatic_variants.csv")

if (Sys.which("bcftools") == "") {
  stop("bcftools not on PATH — run 'module load bcftools/1.21' before this script")
}

# ---------------------------------------------------------------- extractor --
# 117 fields, in the order extract_smlv_somatic.sh emits them.
VCF_COLS <- c(
  # 1-27: locus + scalar/flag INFO
  "chrom", "pos", "vcf_id", "ref", "alt", "filter", "vcf_qual",
  "info_dp", "mq", "mq_rank_sum", "read_pos_rank_sum",
  "fraction_informative_reads", "soft_clip_ratio", "pon_count", "gnomad_af",
  "hmf_hotspot", "sage_hotspot", "sage_novel", "sage_rescue",
  "sage_hotspot_rescue", "pcgr_tier_rescue", "clinical_potential_rescue",
  "hotspot", "giab_conf", "encode",
  "rescued_filters_existing", "rescued_filters_pending",
  # 28-29
  "difficult_regions", "n_csq",
  # 30-89: the picked PCGR_CSQ entry, one column per VEP subfield
  "csq_allele", "consequence", "impact", "symbol", "gene", "feature_type",
  "feature", "biotype", "exon", "intron", "hgvsc", "hgvsp", "cdna_position",
  "cds_position", "protein_position", "amino_acids", "codons",
  "existing_variation", "allele_num", "distance", "strand", "flags", "pick",
  "variant_class", "symbol_source", "hgnc_id", "canonical", "mane_select",
  "mane_plus_clinical", "tsl", "appris", "ccds", "ensp", "swissprot", "trembl",
  "uniparc", "uniprot_isoform", "refseq", "domains", "hgvs_offset",
  "af_1kg", "afr_af_1kg", "amr_af_1kg", "eas_af_1kg", "eur_af_1kg", "sas_af_1kg",
  "vep_gnomad_af", "vep_gnomad_afr_af", "vep_gnomad_amr_af",
  "vep_gnomad_asj_af", "vep_gnomad_eas_af", "vep_gnomad_fin_af",
  "vep_gnomad_nfe_af", "vep_gnomad_oth_af", "vep_gnomad_sas_af",
  "clin_sig", "somatic", "pheno", "check_ref", "nearest_exon_jb",
  # 90-117: FORMAT, tumour then normal
  "t_gt", "t_ad", "t_af", "t_dp", "t_f1r2", "t_f2r1", "t_sb", "t_mb", "t_sq",
  "t_ps", "t_sage_ad", "t_sage_af", "t_sage_dp", "t_sage_sb",
  "n_gt", "n_ad", "n_af", "n_dp", "n_f1r2", "n_f2r1", "n_sb", "n_mb", "n_sq",
  "n_ps", "n_sage_ad", "n_sage_af", "n_sage_dp", "n_sage_sb"
)

# -------------------------------------------------------------- tiers TSV ----
# UPPER_SNAKE source -> output. Columns that describe *PCGR's* picked
# transcript get a `pcgr_` prefix; the locus-level annotation does not.
# Dropped as redundant against the CSQ columns or as noise: GENOME_VERSION
# (constant grch38), VCF_SAMPLE_ID (= tumor_lib), VARIANT_CLASS, CANONICAL,
# ENSEMBL_GENE_ID, HGVSp, HGVSc, CCDS, REFSEQ_MRNA, REFSEQ_PEPTIDE,
# UNIPROT_ID, TRANSCRIPT_MANE_SELECT, GLOBAL_AF_GNOMAD, GLOBAL_AF_1KG,
# DP_TUMOR/AF_TUMOR/DP_CONTROL/AF_CONTROL (= t_dp/t_af/n_dp/n_af),
# ONCOGENE_EVIDENCE/TUMOR_SUPPRESSOR_EVIDENCE (citation dumps), and
# VEP_ALL_CSQ (full per-transcript dump; n_csq carries the count).
TIERS_MAP <- c(
  GENOMIC_CHANGE              = "genomic_change",
  DBSNPRSID                   = "dbsnp_rs",
  COSMIC_MUTATION_ID          = "cosmic_mutation_id",
  TCGA_PANCANCER_COUNT        = "tcga_pancancer_count",
  TCGA_FREQUENCY              = "tcga_frequency",
  ICGC_PCAWG_OCCURRENCE       = "icgc_pcawg_occurrence",
  CLINVAR                     = "clinvar",
  CLINVAR_CLNSIG              = "clinvar_clnsig",
  NCER_PERCENTILE             = "ncer_percentile",
  SIMPLEREPEATS_HIT           = "simplerepeats_hit",
  WINMASKER_HIT               = "winmasker_hit",
  OPENTARGETS_RANK            = "opentargets_rank",
  CHEMBL_COMPOUND_ID          = "chembl_compound_id",
  CHEMBL_COMPOUND_TERMS       = "chembl_compound_terms",
  PUTATIVE_DRIVER_MUTATION    = "putative_driver_mutation",
  CHASMPLUS_DRIVER            = "chasmplus_driver",
  CHASMPLUS_TTYPE             = "chasmplus_ttype",
  MUTATION_HOTSPOT            = "mutation_hotspot",
  MUTATION_HOTSPOT_TRANSCRIPT = "mutation_hotspot_transcript",
  MUTATION_HOTSPOT_CANCERTYPE = "mutation_hotspot_cancertype",
  CALL_CONFIDENCE             = "call_confidence",
  TIER                        = "pcgr_tier",
  TIER_DESCRIPTION            = "pcgr_tier_description",
  SYMBOL                      = "pcgr_symbol",
  GENE_NAME                   = "pcgr_gene_name",
  ENTREZ_ID                   = "pcgr_entrez_id",
  ENSEMBL_TRANSCRIPT_ID       = "pcgr_transcript",
  CONSEQUENCE                 = "pcgr_consequence",
  PROTEIN_CHANGE              = "pcgr_protein_change",
  PROTEIN_DOMAIN              = "pcgr_protein_domain",
  CDS_CHANGE                  = "pcgr_cds_change",
  CODING_STATUS               = "pcgr_coding_status",
  EXONIC_STATUS               = "pcgr_exonic_status",
  EFFECT_PREDICTIONS          = "pcgr_effect_predictions",
  ONCOGENE                    = "pcgr_oncogene",
  TUMOR_SUPPRESSOR            = "pcgr_tumor_suppressor"
)

# ---------------------------------------------------------------- typing -----
NUMERIC_COLS <- c(
  "pos", "start_position", "end_position", "vcf_qual", "info_dp", "mq",
  "mq_rank_sum", "read_pos_rank_sum", "fraction_informative_reads",
  "soft_clip_ratio", "pon_count", "gnomad_af", "n_csq", "distance", "strand",
  "hgvs_offset",
  "af_1kg", "afr_af_1kg", "amr_af_1kg", "eas_af_1kg", "eur_af_1kg", "sas_af_1kg",
  "vep_gnomad_af", "vep_gnomad_afr_af", "vep_gnomad_amr_af",
  "vep_gnomad_asj_af", "vep_gnomad_eas_af", "vep_gnomad_fin_af",
  "vep_gnomad_nfe_af", "vep_gnomad_oth_af", "vep_gnomad_sas_af",
  "t_af", "t_dp", "t_sq", "t_ps", "t_sage_af", "t_sage_dp",
  "n_af", "n_dp", "n_sq", "n_ps", "n_sage_af", "n_sage_dp",
  "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count",
  "tcga_pancancer_count", "ncer_percentile", "opentargets_rank"
)

# NOT in NUMERIC_COLS, despite looking scalar: t_sage_sb / n_sage_sb are declared
# SAGE_SB,Number=1,Type=Float in every VCF header but actually carry a
# comma-joined pair (e.g. "0.53,0.593"), so as.numeric() would silently NA the
# lot. Left as character alongside the other comma-joined FORMAT fields
# (t_ad, t_sb, t_mb, t_f1r2, t_f2r1).

# bcftools renders INFO Flags as "1" when set and "." when not.
FLAG_COLS <- c("hmf_hotspot", "sage_hotspot", "sage_novel", "sage_rescue",
               "sage_hotspot_rescue", "pcgr_tier_rescue",
               "clinical_potential_rescue", "hotspot", "giab_conf", "encode")

# TRUE/FALSE strings out of the PCGR tiers TSV.
TF_COLS <- c("simplerepeats_hit", "winmasker_hit", "putative_driver_mutation",
             "pcgr_oncogene", "pcgr_tumor_suppressor")

# Consequence -> MAF Variant_Classification. vcf2maf resolves the &-joined
# consequence list by taking its first (most severe) term and looking it up
# here, defaulting to Targeted_Region; this table is the one observed in
# results/BrCa_combine_sash_filter.maf, extended with the SO terms vcf2maf
# handles that this cohort happens not to contain.
SO_TO_MAF <- c(
  splice_acceptor_variant                    = "Splice_Site",
  splice_donor_variant                       = "Splice_Site",
  transcript_ablation                        = "Splice_Site",
  exon_loss_variant                          = "Splice_Site",
  stop_gained                                = "Nonsense_Mutation",
  stop_lost                                  = "Nonstop_Mutation",
  start_lost                                 = "Translation_Start_Site",
  initiator_codon_variant                    = "Translation_Start_Site",
  inframe_insertion                          = "In_Frame_Ins",
  disruptive_inframe_insertion               = "In_Frame_Ins",
  inframe_deletion                           = "In_Frame_Del",
  disruptive_inframe_deletion                = "In_Frame_Del",
  missense_variant                           = "Missense_Mutation",
  conservative_missense_variant              = "Missense_Mutation",
  rare_amino_acid_variant                    = "Missense_Mutation",
  coding_sequence_variant                    = "Missense_Mutation",
  splice_region_variant                      = "Splice_Region",
  synonymous_variant                         = "Silent",
  stop_retained_variant                      = "Silent",
  start_retained_variant                     = "Silent",
  incomplete_terminal_codon_variant          = "Silent",
  intron_variant                             = "Intron",
  transcript_amplification                   = "Intron",
  INTRAGENIC                                 = "Intron",
  mature_miRNA_variant                       = "RNA",
  exon_variant                               = "RNA",
  non_coding_exon_variant                    = "RNA",
  non_coding_transcript_exon_variant         = "RNA",
  non_coding_transcript_variant              = "RNA",
  nc_transcript_variant                      = "RNA",
  NMD_transcript_variant                     = "Silent",
  `5_prime_UTR_variant`                      = "5'UTR",
  `5_prime_UTR_premature_start_codon_gain_variant` = "5'UTR",
  `3_prime_UTR_variant`                      = "3'UTR",
  upstream_gene_variant                      = "5'Flank",
  downstream_gene_variant                    = "3'Flank",
  intergenic_variant                         = "IGR",
  intergenic_region                          = "IGR",
  regulatory_region_variant                  = "IGR",
  regulatory_region                          = "IGR",
  TF_binding_site_variant                    = "IGR"
)

# MAF Variant_Type + the MAF's trimmed coordinates/alleles. vcf2maf strips the
# shared anchor base from indels: an insertion becomes ref "-" spanning the two
# flanking bases, a deletion becomes alt "-" spanning the deleted bases.
maf_coords <- function(pos, ref, alt) {
  pos <- as.numeric(pos)
  nr  <- nchar(ref)
  na_ <- nchar(alt)

  vtype <- ifelse(nr == na_,
                  c("SNP", "DNP", "TNP")[pmin(nr, 3)],
                  ifelse(nr < na_, "INS", "DEL"))
  vtype[nr == na_ & nr > 3] <- "ONP"

  is_ins <- nr < na_ & substr(alt, 1, nr) == ref
  is_del <- nr > na_ & substr(ref, 1, na_) == alt

  start <- ifelse(is_ins, pos + nr - 1, ifelse(is_del, pos + na_, pos))
  end   <- ifelse(is_ins, pos + nr,     pos + nr - 1)
  r_all <- ifelse(is_ins, "-", ifelse(is_del, substring(ref, na_ + 1), ref))
  a_all <- ifelse(is_ins, substring(alt, nr + 1), ifelse(is_del, "-", alt))

  # Complex indels (neither allele a prefix of the other) keep untrimmed
  # coordinates; vcf2maf still labels them INS/DEL by length.
  tibble(variant_type = vtype, start_position = start, end_position = end,
         ref_allele = r_all, tumor_seq_allele2 = a_all)
}

classify <- function(consequence, variant_type) {
  first <- sub("[&,].*$", "", consequence)
  cls   <- unname(SO_TO_MAF[first])
  cls[first == "frameshift_variant"] <-
    ifelse(variant_type[first == "frameshift_variant"] == "INS",
           "Frame_Shift_Ins", "Frame_Shift_Del")
  cls[first == "protein_altering_variant"] <-
    ifelse(variant_type[first == "protein_altering_variant"] == "INS",
           "In_Frame_Ins", "In_Frame_Del")
  ifelse(is.na(cls), "Targeted_Region", cls)
}

# readr writes doubles at full 17-digit precision (0.1078 -> 0.10780000000000001),
# which diverges from both the VCF's own text and the reference MAF. as.character()
# gives R's shortest 15-significant-digit representation, which round-trips and
# reproduces the source string exactly; only the handful of values it would render
# in scientific notation need the slower per-element format(). Same problem, and
# same fix, as fixed_notation() in cohort_summary.R -- vectorised format() is not
# an option because it pads every value to a common width.
fixed_notation <- function(x) {
  out <- as.character(x)
  sci <- !is.na(out) & grepl("e", out, fixed = TRUE)
  if (any(sci)) {
    out[sci] <- vapply(x[sci],
                       function(v) format(v, scientific = FALSE, trim = TRUE, digits = 15),
                       character(1))
  }
  out
}

# vcf2maf strips the versioned accession that VEP prefixes onto HGVS strings
# ("ENST00000379370.7:c.2182G>A" -> "c.2182G>A"); this table matches that, since
# `feature` / `ensp` already name the transcript and the version is a constant
# property of the GENCODE 39 cache, not per-row information.
strip_hgvs_prefix <- function(x) sub("^[^:]*:", "", x)

# MAF HGVSp_Short: VEP's three-letter HGVSp with amino acids in one-letter form
# and Ter as "*" (p.Gln307Ter -> p.Q307*, p.Ter350SerextTer51 -> p.*350Sext*51).
AA3TO1 <- c(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E",
            Gly="G", His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F",
            Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V", Ter="*",
            Sec="U", Pyl="O", Asx="B", Glx="Z", Xaa="X")
AA3_RE <- paste0("(", paste(names(AA3TO1), collapse = "|"), ")")
hgvsp_to_short <- function(x) {
  out <- x
  ok  <- !is.na(out)
  if (!any(ok)) return(out)
  m <- gregexpr(AA3_RE, out[ok])
  regmatches(out[ok], m) <- lapply(regmatches(out[ok], m),
                                   function(v) unname(AA3TO1[v]))
  out
}

# Comma-joined AD -> ref/alt counts. "." (no call) yields NA.
ad_part <- function(ad, i) {
  suppressWarnings(as.numeric(vapply(strsplit(ad, ",", fixed = TRUE),
                                     function(x) if (length(x) >= i) x[i] else NA_character_,
                                     character(1))))
}

read_tiers <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "smlv_somatic", "report", "pcgr",
            paste0(tumor_lib, ".pcgr_acmg.grch38.snvs_indels.tiers.tsv"))
  if (!file.exists(p)) return(NULL)
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(NULL)

  present <- intersect(names(TIERS_MAP), names(df))
  keys <- df[, c("CHROM", "POS", "REF", "ALT"), drop = FALSE]
  out  <- df[, present, drop = FALSE]
  names(out) <- unname(TIERS_MAP[present])
  out$.chrom <- sub("^chr", "", keys$CHROM)
  out$.pos   <- keys$POS
  out$.ref   <- keys$REF
  out$.alt   <- keys$ALT
  # Explicit join indicator: every value column here can legitimately be NA,
  # so an unmatched row is not distinguishable from a matched-but-empty one.
  out$.matched <- TRUE
  out
}

extract_vcf <- function(vcf, tumor_lib, normal_lib) {
  con <- pipe(sprintf("bash %s %s %s %s", shQuote(EXTRACT_SH), shQuote(vcf),
                      shQuote(tumor_lib), shQuote(normal_lib)))
  df <- read_tsv_quiet(con, col_names = VCF_COLS, col_types = cols(.default = "c"),
                       na = character(0), quote = "")
  # A PASS VCF always has records, and pipe() will not propagate the
  # extractor's exit status — so an empty or ragged read means bcftools or the
  # CSQ-layout assertion failed, not that the tumour has no variants.
  if (ncol(df) != length(VCF_COLS)) {
    stop(sprintf("extractor returned %d columns for %s, expected %d",
                 ncol(df), tumor_lib, length(VCF_COLS)))
  }
  if (nrow(df) == 0) {
    stop(sprintf("extractor returned no rows for %s (%s)", tumor_lib, vcf))
  }
  if (anyNA(df$chrom) || anyNA(df$pos)) {
    stop(sprintf("extractor returned short row(s) for %s (%s)", tumor_lib, vcf))
  }
  df
}

# ------------------------------------------------------------------ main -----
ss <- load_samplesheet()

dir_create(OUT_DIR)
if (file.exists(OUT_CSV)) invisible(file.remove(OUT_CSV))

n_rows <- 0L
n_tumours <- 0L
missing_vcf <- character(0)
missing_tiers <- character(0)
unjoined <- character(0)

for (i in seq_len(nrow(ss))) {
  rec <- ss[i, ]
  vcf <- path(rec$src_dir, "smlv_somatic", "filter",
              paste0(rec$tumor_lib, ".pass.vcf.gz"))
  if (!file.exists(vcf)) {
    missing_vcf <- c(missing_vcf, rec$tumor_gussid)
    next
  }

  df <- extract_vcf(vcf, rec$tumor_lib, rec$normal_lib)

  # bcftools writes "." for absent scalars; flags become logical below.
  df <- df %>% mutate(across(-all_of(FLAG_COLS), ~ ifelse(.x == "." | .x == "", NA_character_, .x)))

  coords <- maf_coords(df$pos, df$ref, df$alt)
  df <- bind_cols(df, coords) %>%
    mutate(
      chrom                  = sub("^chr", "", chrom),
      variant_classification = classify(ifelse(is.na(consequence), "", consequence),
                                        variant_type),
      t_ref_count = ad_part(t_ad, 1), t_alt_count = ad_part(t_ad, 2),
      n_ref_count = ad_part(n_ad, 1), n_alt_count = ad_part(n_ad, 2),
      hgvsc       = strip_hgvs_prefix(hgvsc),
      hgvsp       = strip_hgvs_prefix(hgvsp),
      hgvsp_short = hgvsp_to_short(hgvsp),
      # The MAF's Hugo_Symbol is SYMBOL, falling back to the transcript ID and
      # then to the literal "Unknown" -- reproduced so the table is drop-in for
      # MAF tooling. `symbol` keeps VEP's value (NA when absent).
      hugo_symbol = ifelse(!is.na(symbol), symbol,
                    ifelse(!is.na(feature), feature, "Unknown"))
    )

  tiers <- read_tiers(rec$src_dir, rec$tumor_lib)
  if (is.null(tiers)) {
    missing_tiers <- c(missing_tiers, rec$tumor_gussid)
    for (nm in unname(TIERS_MAP)) df[[nm]] <- NA_character_
  } else {
    before <- nrow(df)
    df <- df %>%
      left_join(tiers, by = c(chrom = ".chrom", pos = ".pos",
                              ref = ".ref", alt = ".alt"))
    if (nrow(df) != before) {
      stop(sprintf("tiers join changed row count for %s (%d -> %d): duplicate keys",
                   rec$tumor_gussid, before, nrow(df)))
    }
    n_unmatched <- sum(is.na(df$.matched))
    if (n_unmatched > 0) {
      unjoined <- c(unjoined, sprintf("%s (%d)", rec$tumor_gussid, n_unmatched))
    }
    for (nm in setdiff(unname(TIERS_MAP), names(df))) df[[nm]] <- NA_character_
  }

  df <- df %>%
    mutate(
      across(any_of(FLAG_COLS), ~ .x == "1"),
      across(any_of(TF_COLS),   ~ toupper(.x) == "TRUE"),
      pick      = !is.na(pick) & pick == "1",
      canonical = !is.na(canonical) & toupper(canonical) == "YES",
      across(any_of(NUMERIC_COLS), ~ fixed_notation(suppressWarnings(as.numeric(.x))))
    ) %>%
    bind_cols(rec %>% select(any_of(IDENTITY_COLS)) %>% slice(rep(1, nrow(df)))) %>%
    select(
      any_of(IDENTITY_COLS),
      # locus, in both VCF and MAF conventions
      chrom, pos, ref, alt, filter, variant_type, variant_classification,
      start_position, end_position, ref_allele, tumor_seq_allele2,
      genomic_change, dbsnp_rs, vcf_id, vcf_qual,
      # VEP annotation for the picked transcript
      hugo_symbol, symbol, gene, feature, feature_type, transcript_biotype = biotype,
      consequence, impact, exon, intron, hgvsc, hgvsp, hgvsp_short, cdna_position,
      cds_position, protein_position, amino_acids, codons, existing_variation,
      distance, strand, flags, pick, canonical, variant_class, symbol_source,
      hgnc_id, mane_select, mane_plus_clinical, tsl, appris, ccds, ensp,
      swissprot, trembl, uniparc, uniprot_isoform, refseq, domains,
      hgvs_offset, clin_sig, somatic, pheno, check_ref, nearest_exon_jb, n_csq,
      # population frequencies
      af_1kg, afr_af_1kg, amr_af_1kg, eas_af_1kg, eur_af_1kg, sas_af_1kg,
      starts_with("vep_gnomad_"), gnomad_af,
      # read support
      t_dp, t_ref_count, t_alt_count, t_af, t_ad, t_gt, t_f1r2, t_f2r1, t_sb,
      t_mb, t_sq, t_ps, t_sage_ad, t_sage_af, t_sage_dp, t_sage_sb,
      n_dp, n_ref_count, n_alt_count, n_af, n_ad, n_gt, n_f1r2, n_f2r1, n_sb,
      n_mb, n_sq, n_ps, n_sage_ad, n_sage_af, n_sage_dp, n_sage_sb,
      # caller / filtering context
      info_dp, mq, mq_rank_sum, read_pos_rank_sum, fraction_informative_reads,
      soft_clip_ratio, pon_count, difficult_regions, giab_conf, encode,
      hotspot, hmf_hotspot, sage_hotspot, sage_novel, sage_rescue,
      sage_hotspot_rescue, pcgr_tier_rescue, clinical_potential_rescue,
      rescued_filters_existing, rescued_filters_pending,
      # PCGR clinical layer (pcgr_* = PCGR's own picked transcript)
      pcgr_tier, pcgr_tier_description, call_confidence, pcgr_symbol,
      pcgr_gene_name, pcgr_entrez_id, pcgr_transcript, pcgr_consequence,
      pcgr_protein_change, pcgr_protein_domain, pcgr_cds_change,
      pcgr_coding_status, pcgr_exonic_status, pcgr_effect_predictions,
      pcgr_oncogene, pcgr_tumor_suppressor, putative_driver_mutation,
      chasmplus_driver, chasmplus_ttype, mutation_hotspot,
      mutation_hotspot_transcript, mutation_hotspot_cancertype,
      cosmic_mutation_id, tcga_pancancer_count, tcga_frequency,
      icgc_pcawg_occurrence, clinvar, clinvar_clnsig, ncer_percentile,
      opentargets_rank, chembl_compound_id, chembl_compound_terms,
      simplerepeats_hit, winmasker_hit
    )

  write_csv(df, OUT_CSV, na = "NA", append = n_tumours > 0)
  n_rows    <- n_rows + nrow(df)
  n_tumours <- n_tumours + 1L
}

if (length(missing_vcf) > 0) {
  message(sprintf("WARN: %d tumour(s) have no smlv_somatic PASS VCF and contribute 0 rows: %s",
                  length(missing_vcf), paste(missing_vcf, collapse = ", ")))
}
if (length(missing_tiers) > 0) {
  message(sprintf("WARN: %d tumour(s) have no PCGR tiers TSV; pcgr_* columns are NA: %s",
                  length(missing_tiers), paste(missing_tiers, collapse = ", ")))
}
if (length(unjoined) > 0) {
  message(sprintf("WARN: variants with no PCGR tiers match in %d tumour(s): %s",
                  length(unjoined), paste(unjoined, collapse = ", ")))
}

cat(sprintf("Wrote %d variant row(s) from %d tumour(s) to %s\n",
            n_rows, n_tumours, OUT_CSV))
