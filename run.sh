#! /bin/bash
# Run the wgs-aggregate workflow.
#
# Standalone:   ./run.sh [step|target ...] [extra snakemake args]
#   e.g.  ./run.sh sash
#         ./run.sh germline -j16
#         ./run.sh data/wgs/sash_aggregate/purple_drivers_somatic.csv   # explicit target
#
# As a submodule (installed under <parent>/modules/wgs-aggregate), set
# WGS_AGGREGATE_CONFIG to the parent repo's config and invoke this script
# (mirrors bulk-rna-preprocessing's run.sh).
#
# Uses --use-conda (the per-rule env is env/environment.yaml). On a compute node
# without internet, build that env once first (on a data-mover node) so the run is
# offline-safe. The global NCI profile (PBS submission) is used when present as a
# submodule under profiles/global/nci.
set -e
eval "$(conda shell.bash hook)"
conda activate "${WGS_AGGREGATE_ENV:-snakemake_7.32.4}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

global_profile=""
if [[ "$(hostname)" == *"nci"* && -f "$HERE/profiles/global/nci/config.yaml" ]]; then
    global_profile="--profile $HERE/profiles/global/nci"
    mkdir -p logs/joblogs
fi

# Arg order matters. `--cores N` first gives a default core count (snakemake 7
# errors without one off-profile); then "$@" (targets / extra flags) -- placed
# where the preceding option is bounded, so targets bind as positionals and a
# user-supplied -j/--cores later wins. --configfile/--profile (nargs="+") go LAST,
# with nothing after them to swallow. Local run (no PBS); --rerun-triggers mtime
# matches the dataset's incremental-rebuild expectation.
snakemake \
    --cores "${CORES:-${THREADS:-4}}" \
    "$@" \
    $global_profile \
    --rerun-triggers mtime \
    --workflow-profile "$HERE/profiles/workflow" \
    --snakefile "$HERE/workflow/Snakefile" \
    --configfile "${WGS_AGGREGATE_CONFIG:-$HERE/config/config.yaml}"
