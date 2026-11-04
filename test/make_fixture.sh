#!/bin/bash
# Regenerate the committed synthetic test fixture under test/fixtures/.
#
# A tiny 2-donor cohort exercising the two representative steps the golden test
# covers: `sash` PURPLE somatic drivers (agg_drivers.R) and `dragen` WGS coverage
# metrics (agg_coverage_metrics.R). Inputs are synthetic and data-independent.
# The golden outputs (test/golden/) are produced from these by THIS code -- see
# run_test.sh. Re-run this only to change the fixture; then regenerate golden.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
F="$HERE/fixtures"
rm -rf "$F"
mkdir -p "$F"

# --- sash samplesheet + per-sample PURPLE driver catalogs ---------------------
cat > "$F/sash.samplesheet.csv" <<'CSV'
tumor_gussid,donor_id,tumor_tissue_id,tumor_lib,normal_lib,sash_run_id,sash_dir
G1,D1,TIS1,T1,N1,RUN1,test/fixtures/sash/RUN1/T1
G2,D2,TIS2,T2,N2,RUN2,test/fixtures/sash/RUN2/T2
CSV

write_drivers() {  # <dir> <gene> <chrom> <band> <transcript> <likelihood>
  mkdir -p "$1/purple"
  { printf 'chromosome\tchromosomeBand\tgene\ttranscript\tisCanonical\tdriver\tcategory\tlikelihoodMethod\tdriverLikelihood\tmissense\tnonsense\tsplice\tinframe\tframeshift\tbiallelic\tminCopyNumber\tmaxCopyNumber\n'
    printf '%s\t%s\t%s\t%s\ttrue\tMUTATION\tTSG\tDNDS\t%s\t1\t0\t0\t0\t0\ttrue\t1.20\t2.30\n' "$3" "$4" "$2" "$5" "$6"
  } > "$1/purple/$(basename "$1").purple.driver.catalog.somatic.tsv"
}
write_drivers "$F/sash/RUN1/T1" BRCA1 chr17 q21.31 ENST00000357654 0.95
write_drivers "$F/sash/RUN2/T2" BRCA2 chr13 q13.1  ENST00000380152 0.88

# --- dragen runs + per-run coverage metric CSVs -------------------------------
cat > "$F/dragen_runs.csv" <<'CSV'
run_id,donor_id,tumor_library,germline_library
RUN1,D1,T1,N1
RUN2,D2,T2,N2
CSV

cov() {  # <path> <avg> <uniformity>
  mkdir -p "$(dirname "$1")"
  { printf 'COVERAGE SUMMARY,,Average alignment coverage over genome,%s\n' "$2"
    printf 'COVERAGE SUMMARY,,Uniformity of coverage (PCT > 0.2*mean) over genome,%s\n' "$3"
  } > "$1"
}
for r in RUN1:T1:N1:35.2:95.1 RUN2:T2:N2:40.8:96.3; do
  IFS=: read -r run t n avg unif <<<"$r"
  d="$F/tumor-normal/$run"
  cov "$d/${t}_dragen_somatic/${t}.wgs_coverage_metrics_tumor.csv"  "$avg" "$unif"
  cov "$d/${t}_dragen_somatic/${t}.wgs_coverage_metrics_normal.csv" "30.1" "94.0"
  cov "$d/${n}_dragen_germline/${n}.wgs_coverage_metrics.csv"       "31.5" "93.2"
done

echo "fixture written -> $F"
