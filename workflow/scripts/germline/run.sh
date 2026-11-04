#!/usr/bin/env bash
# Build the cohort-wide multi-sample germline genotype VCF + per-project subsets
# for single-cell genotype demux (issue #56).
#
# Source: per-donor DRAGEN genome-wide germline hard-filtered VCFs
#   data/wgs/tumor-normal/<run>/<lib>_dragen_germline/<lib>.hard-filtered.vcf.gz
# Outputs (data/wgs/germline_cohort_vcf/):
#   cohort.germline.common_snps.vcf.gz(.csi)   multi-sample, donor-ID columns
#   donor_manifest.csv                         donor -> source VCF provenance
#   projects/<project>.germline.common_snps.vcf.gz(.csi)
#   projects/project_donors.csv                per-project donor sets
#
# Panel definition (all decisions documented in README.md):
#   autosomal (chr1-22), PASS, biallelic SNPs, GT only, missing genotype -> 0/0
#   (a donor absent at a site is hom-ref, not missing), cohort MAF >= 0.05.
#
# Requires bcftools/bgzip/tabix on PATH (conda env `match_vcf`).
# Usage: conda activate match_vcf && bash modules/germline_cohort_vcf/run.sh
set -euo pipefail

# REPO is the dataset root (where config/ + data/ live). Honour an existing REPO env
# (set by the Snakemake rule); otherwise fall back to the legacy in-tree layout.
REPO="${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}"; export REPO
cd "$REPO"
SCRIPT_DIR="${GERMLINE_SCRIPT_DIR:-$REPO/modules/germline_cohort_vcf/scripts}"
OUT="$REPO/data/wgs/germline_cohort_vcf"
COHORT="$OUT/cohort.germline.common_snps.vcf.gz"
# Panel knobs: honour config (params.germline.*), falling back to the historical
# defaults (chr1..22, MAF >= 0.05) so an unset env reproduces the prior output.
AUTOSOMES="$(printf "chr%d," $(seq 1 "${WGS_AUTOSOMES:-22}") | sed 's/,$//')"
MAF_MIN="${WGS_MAF_MIN:-0.05}"
THREADS="${THREADS:-$(nproc)}"

# bcftools does all the VCF work; pull it from the environment modules if the
# active env (e.g. the `meta` conda env this stage runs under) lacks it.
if ! command -v bcftools >/dev/null 2>&1; then
  module load bcftools/1.21 2>/dev/null || true
fi
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH — 'module load bcftools/1.21' (or 'conda activate match_vcf') first." >&2; exit 1; }

# Scratch work dir (NOT jobfs — the ~195 per-donor + per-chr VCFs are several GB,
# far over a default 100 MB jobfs). Override with GERMLINE_TMPDIR if needed.
WORK="$(mktemp -d -p "${GERMLINE_TMPDIR:-/scratch/a56/jr9959/tmp}" germline_cohort.XXXXXX)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$OUT/projects"

echo "=== [1/4] Building donor -> germline VCF manifest ==="
python "$SCRIPT_DIR/build_manifest.py"

# --- [2/4] Per-donor prep: normalise -> autosomal PASS biallelic SNPs -> GT only
#           -> rename sample to donor_id. Parallel across donors. ------------------
prep_one() {
  local donor="$1" src="$2"
  local tmp="$WORK/$donor.tmp.vcf.gz" final="$WORK/$donor.vcf.gz"
  bcftools norm -m-any "$src" -Ou \
    | bcftools view -f PASS -v snps -m2 -M2 -t "$AUTOSOMES" -Ou \
    | bcftools annotate -x 'INFO,^FORMAT/GT' -Oz -o "$tmp"
  printf '%s\n' "$donor" > "$WORK/$donor.name"
  bcftools reheader -s "$WORK/$donor.name" "$tmp" -o "$final"
  bcftools index -f "$final"
  rm -f "$tmp" "$WORK/$donor.name"
}

echo "=== [2/4] Per-donor prep ($THREADS-way parallel) ==="
# Bounded-parallel over donors: donor_id is field 1, source_vcf the last field.
while IFS=, read -r donor src; do
  [ -n "$donor" ] || continue
  prep_one "$donor" "$src" &
  while (( $(jobs -rp | wc -l) >= THREADS )); do wait -n; done
done < <(tail -n +2 "$OUT/donor_manifest.csv" | awk -F, '{print $1","$NF}')
wait
n_prepped=$(find "$WORK" -maxdepth 1 -name '*.vcf.gz' | wc -l)
n_expected=$(( $(wc -l < "$OUT/donor_manifest.csv") - 1 ))
echo "prepared $n_prepped / $n_expected per-donor VCFs"
[ "$n_prepped" -eq "$n_expected" ] || { echo "ERROR: prepared $n_prepped but expected $n_expected donors — aborting rather than merging an incomplete cohort"; exit 1; }

# --- [3/4] Merge -> cohort common-SNP panel ------------------------------------
# Merge per chromosome in parallel (each donor VCF is indexed, so `merge -r` does
# random access), then concat. A single whole-genome merge of ~195 samples is a
# multi-hour serial bottleneck; per-chromosome fans it out ~22x.
echo "=== [3/4] Per-chromosome merge + common-SNP filter (MAF >= $MAF_MIN) ==="
find "$WORK" -maxdepth 1 -name '*.vcf.gz' | sort > "$WORK/filelist.txt"

merge_chr() {
  local chr="$1"
  # merge -0: genotypes missing from a donor's VCF are hom-ref (0/0), not ./.
  bcftools merge -0 -r "$chr" -l "$WORK/filelist.txt" -Ou \
    | bcftools norm -m-any -Ou \
    | bcftools view -v snps -m2 -M2 -Ou \
    | bcftools +fill-tags -Ou -- -t AN,AC,AF,MAF \
    | bcftools view -i "INFO/MAF>=$MAF_MIN" -Oz -o "$WORK/merged.$chr.vcf.gz"
  bcftools index -f "$WORK/merged.$chr.vcf.gz"
}
export -f merge_chr; export WORK MAF_MIN

for n in $(seq 1 22); do
  merge_chr "chr$n" &
  while (( $(jobs -rp | wc -l) >= THREADS )); do wait -n; done
done
wait

# concat autosomes in order into the cohort VCF
chr_files=(); for n in $(seq 1 22); do chr_files+=("$WORK/merged.chr$n.vcf.gz"); done
bcftools concat "${chr_files[@]}" --threads "$THREADS" -Oz -o "$COHORT"
bcftools index -f "$COHORT"
n_sites=$(bcftools index -n "$COHORT")
n_samp=$(bcftools query -l "$COHORT" | wc -l)
echo "cohort VCF: $n_samp donors x $n_sites common biallelic SNPs"

# --- [4/4] Per-project subsets --------------------------------------------------
echo "=== [4/4] Per-project subsets ==="
python "$SCRIPT_DIR/project_donors.py"
while IFS= read -r proj; do
  donors="$(cat "$OUT/projects/$proj.donors.txt")"
  [ -n "$donors" ] || { echo "  $proj: no donors, skipping"; continue; }
  pv="$OUT/projects/$proj.germline.common_snps.vcf.gz"
  # -c 1:minor: keep only sites polymorphic among this project's donors
  bcftools view -s "$donors" -c 1:minor "$COHORT" --threads "$THREADS" -Oz -o "$pv"
  bcftools index -f "$pv"
  echo "  $proj: $(bcftools index -n "$pv") polymorphic sites"
done < <(tail -n +2 "$OUT/projects/project_donors.csv" | cut -d, -f1)

echo "=== Done -> $OUT ==="
