#!/usr/bin/env bash
# Flatten one sash somatic PASS VCF into a headerless TSV on stdout.
#
# Usage: extract_smlv_somatic.sh <pass.vcf.gz> <tumor_lib> <normal_lib>
#
# Emits 117 fields per record, in the order agg_smlv_somatic.R names them:
#   1-27    locus + scalar/flag INFO fields
#   28      difficult_regions  (the 12 DIFFICULT_* flags, pipe-joined)
#   29      n_csq              (number of PCGR_CSQ transcript entries)
#   30-89   the picked PCGR_CSQ entry, split into its 60 VEP subfields
#   90-103  FORMAT fields for the tumour sample
#   104-117 FORMAT fields for the normal sample
#
# Sample order is forced with `-s <tumor>,<normal>` — the VCF itself stores the
# normal first, so relying on file order would silently swap t_/n_ columns.
#
# Transcript pick rule: lowest (canonical, MANE Select, consequence rank,
# protein_coding) tuple, ties broken by CSQ order. This reproduces the pick in
# results/BrCa_combine_sash_filter.maf on 95.0% of variants that have one.
#
# It is NOT VEP's own PICK=1 flag. PCGR runs VEP with an ordering that favours
# consequence severity and biotype, so its flagged pick is canonical only ~25%
# of the time where it disagrees with the MAF, and reproduces the MAF on just
# 78.6% of variants. vcf2maf orders canonical-first, so this does too.
#
# The order was derived empirically, not copied from a vcf2maf version: greedy
# forward selection over the pick keys available in PCGR_CSQ, scored against
# the MAF's own Transcript_ID across 33,265 annotated variants in 12 tumours.
# Named vcf2maf/VEP orders scored 92.3-93.7%; this one scores 95.0%. The
# remaining 5.0% is essentially all ties: 5.69% of variants have >=2
# transcripts equal on every key here, which VEP breaks on transcript length --
# a field PCGR_CSQ does not carry, so it is a genuine ceiling.
#
# Both VEP runs saw the same transcript set (VEP v105 / GENCODE 39; every MAF
# transcript is present in PCGR_CSQ for its variant), so this is purely a
# question of ordering, not of annotation source.
set -euo pipefail

VCF="${1:?usage: extract_smlv_somatic.sh <pass.vcf.gz> <tumor_lib> <normal_lib>}"
TLIB="${2:?missing tumor_lib}"
NLIB="${3:?missing normal_lib}"

command -v bcftools >/dev/null 2>&1 || {
  echo "ERROR: bcftools not on PATH — run 'module load bcftools/1.21' first." >&2
  exit 1
}

# Field order of PCGR_CSQ. Identical across all 208 cohort VCFs; assert rather
# than trust, because every index below is hardcoded against it.
EXPECTED_CSQ='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|RefSeq|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|CHECK_REF|NearestExonJB'

ACTUAL_CSQ="$(bcftools view -h "$VCF" \
  | sed -nE 's/^##INFO=<ID=PCGR_CSQ.*Format: ([^"]*)">$/\1/p' | head -1)"

if [ "$ACTUAL_CSQ" != "$EXPECTED_CSQ" ]; then
  echo "ERROR: PCGR_CSQ field order differs from the expected layout in $VCF" >&2
  echo "  expected: $EXPECTED_CSQ" >&2
  echo "  actual:   $ACTUAL_CSQ" >&2
  exit 1
fi

# The 12 DIFFICULT_* flags, collapsed downstream into one pipe-joined column.
DIFFICULT_TAGS='bad_promoter gc15 gc70to75 gc75to80 gc80 gc80to85
                low_complexity_ditr low_complexity_quadtr
                low_complexity_tandemrepeats low_complexity_tritr
                mappability_nonunique segdup'

FMT='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL'
FMT+='\t%INFO/DP\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum'
FMT+='\t%INFO/FractionInformativeReads\t%INFO/SoftClipRatio'
FMT+='\t%INFO/PON_COUNT\t%INFO/gnomAD_AF'
FMT+='\t%INFO/HMF_HOTSPOT\t%INFO/SAGE_HOTSPOT\t%INFO/SAGE_NOVEL'
FMT+='\t%INFO/SAGE_RESCUE\t%INFO/SAGE_HOTSPOT_RESCUE\t%INFO/PCGR_TIER_RESCUE'
FMT+='\t%INFO/CLINICAL_POTENTIAL_RESCUE\t%INFO/hotspot\t%INFO/GIAB_CONF'
FMT+='\t%INFO/ENCODE'
FMT+='\t%INFO/RESCUED_FILTERS_EXISTING\t%INFO/RESCUED_FILTERS_PENDING'
for tag in $DIFFICULT_TAGS; do FMT+="\t%INFO/DIFFICULT_${tag}"; done
FMT+='\t%INFO/PCGR_CSQ'
FMT+='[\t%GT\t%AD\t%AF\t%DP\t%F1R2\t%F2R1\t%SB\t%MB\t%SQ\t%PS'
FMT+='\t%SAGE_AD\t%SAGE_AF\t%SAGE_DP\t%SAGE_SB]\n'

# shellcheck disable=SC2086  # DIFFICULT_TAGS is a deliberate word list
bcftools query -s "${TLIB},${NLIB}" -f "$FMT" "$VCF" \
  | awk -F'\t' -v OFS='\t' -v tags="$(echo $DIFFICULT_TAGS | tr ' ' ',')" '
    # Ensembl consequence severity, most severe first. An entry carrying several
    # &-joined terms ranks by its most severe one; unknown terms sort last.
    function csq_rank(cons,   m, parts, j, r, best) {
      best = 999
      m = split(cons, parts, "&")
      for (j = 1; j <= m; j++) { r = SORANK[parts[j]]; if (r > 0 && r < best) best = r }
      return best
    }
    BEGIN {
      ntags = split(tags, tagname, ",")
      split("transcript_ablation splice_acceptor_variant splice_donor_variant " \
            "stop_gained frameshift_variant stop_lost start_lost " \
            "transcript_amplification feature_elongation feature_truncation " \
            "inframe_insertion inframe_deletion missense_variant " \
            "protein_altering_variant splice_donor_5th_base_variant " \
            "splice_region_variant splice_donor_region_variant " \
            "splice_polypyrimidine_tract_variant incomplete_terminal_codon_variant " \
            "start_retained_variant stop_retained_variant synonymous_variant " \
            "coding_sequence_variant mature_miRNA_variant 5_prime_UTR_variant " \
            "3_prime_UTR_variant non_coding_transcript_exon_variant intron_variant " \
            "NMD_transcript_variant non_coding_transcript_variant " \
            "coding_transcript_variant upstream_gene_variant downstream_gene_variant " \
            "TFBS_ablation TFBS_amplification TF_binding_site_variant " \
            "regulatory_region_ablation regulatory_region_amplification " \
            "regulatory_region_variant intergenic_variant sequence_variant", so, " ")
      for (j in so) SORANK[so[j]] = j
    }
    {
      # 28: collapse the DIFFICULT_* flag block (fields 28..39) into one string
      difficult = ""
      for (i = 1; i <= ntags; i++) {
        if ($(27 + i) == "1") difficult = (difficult == "" ? tagname[i] : difficult "|" tagname[i])
      }

      # 29-89: pick one CSQ entry and split it. Keys ascending, first wins ties.
      n = split($40, csq, ",")
      chosen = csq[1]; have = 0
      for (i = 1; i <= n; i++) {
        split(csq[i], f, "|")
        k1 = (f[27] == "YES")            ? 0 : 1   # CANONICAL
        k2 = (f[28] != "")               ? 0 : 1   # MANE_SELECT
        k3 = csq_rank(f[2])                        # consequence severity
        k4 = (f[8] == "protein_coding")  ? 0 : 1   # BIOTYPE
        if (!have || k1 < b1 || (k1 == b1 && (k2 < b2 || (k2 == b2 && \
            (k3 < b3 || (k3 == b3 && k4 < b4)))))) {
          chosen = csq[i]; b1 = k1; b2 = k2; b3 = k3; b4 = k4; have = 1
        }
      }
      nf = split(chosen, vep, "|")
      for (i = nf + 1; i <= 60; i++) vep[i] = ""

      out = ""
      for (i = 1; i <= 27; i++) out = out $i OFS
      out = out difficult OFS n
      for (i = 1; i <= 60; i++) out = out OFS vep[i]
      for (i = 41; i <= 68; i++) out = out OFS $i
      print out
    }'
