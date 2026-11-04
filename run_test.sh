#! /bin/bash
# Build the two representative aggregation tables from the committed synthetic
# fixture (test/fixtures/) and, with --check, assert they match the golden outputs
# produced by the SAME code (test/golden/). Fixture + golden are committed, so the
# test is fully self-contained and data-independent.
#
# Usage:
#   ./run_test.sh            # build into test/out/
#   ./run_test.sh --check    # build, then compare test/out to test/golden
#   ./run_test.sh -n         # dry run
#
# Launches via Snakemake with --use-conda (per-rule env = env/environment.yaml),
# matching run.sh. To regenerate the fixture: bash test/make_fixture.sh.
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate "${WGS_AGGREGATE_ENV:-snakemake_7.32.4}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TARGETS=(
    test/out/sash_aggregate/purple_drivers_somatic.csv
    test/out/dragen_aggregate/somatic_coverage_metrics.csv
    test/out/dragen_aggregate/germline_coverage_metrics.csv
)

check=false
sm_args=()
for a in "$@"; do
    if [[ "$a" == "--check" ]]; then check=true; else sm_args+=("$a"); fi
done

snakemake "${TARGETS[@]}" \
    --snakefile "$HERE/workflow/Snakefile" \
    --configfile "$HERE/test/config.yaml" \
    --workflow-profile "$HERE/profiles/workflow" \
    --cores "${CORES:-2}" \
    ${sm_args[@]+"${sm_args[@]}"}

if [[ "$check" == true ]]; then
    echo "=== comparing test/out to test/golden ==="
    for t in "${TARGETS[@]}"; do
        rel="${t#test/out/}"
        python3 "$HERE/test/compare_csv.py" "$HERE/test/golden/$rel" "$HERE/$t"
    done
    echo "OK: test output matches golden"
fi
