#! /bin/bash
# Run the wgs-aggregate workflow.
#
# Two toolchains:
#   Default (NCI / offline): run under the WGS_AGGREGATE_ENV conda env (meta) with
#     bcftools from the environment module and snakemake invoked by full path -- no
#     --use-conda, so it works on internet-less compute nodes. This is how the
#     datasets drive it (and the behaviour of the former run_mod.sh).
#   Portable (WGS_AGGREGATE_USE_CONDA=1): --use-conda with env/environment.yaml +
#     the workflow profile (the bulk-rna-preprocessing standalone standard). Builds
#     the per-rule env on first run, so it needs internet / a prebuilt env.
#
# Standalone:   ./run.sh [step|target ...] [extra snakemake args]
#   e.g.  ./run.sh sash
#         THREADS=16 ./run.sh germline
#         ./run.sh data/wgs/sash_aggregate/purple_drivers_somatic.csv   # explicit target
#
# As a submodule (installed under <parent>/modules/wgs-aggregate), set
# WGS_AGGREGATE_CONFIG to the parent repo's config and invoke this script.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${WGS_AGGREGATE_CONFIG:-$HERE/config/config.yaml}"

source /g/data/a56/conda/miniconda/etc/profile.d/conda.sh 2>/dev/null || eval "$(conda shell.bash hook)"

extra=()
if [[ "${WGS_AGGREGATE_USE_CONDA:-0}" == "1" ]]; then
    conda activate "${WGS_AGGREGATE_LAUNCHER_ENV:-snakemake_7.32.4}"
    SNAKEMAKE="${SNAKEMAKE:-snakemake}"
    extra=(--use-conda --conda-frontend mamba --workflow-profile "$HERE/profiles/workflow")
    if [[ "$(hostname)" == *"nci"* && -f "$HERE/profiles/global/nci/config.yaml" ]]; then
        extra+=(--profile "$HERE/profiles/global/nci"); mkdir -p logs/joblogs
    fi
else
    conda activate "${WGS_AGGREGATE_ENV:-meta}"
    module load bcftools/1.21 2>/dev/null || true
    SNAKEMAKE="${SNAKEMAKE:-/g/data/a56/conda/envs/snakemake_7.32.4/bin/snakemake}"
fi

# Arg order: `--cores N` gives a default core count (snakemake 7 requires one when
# no profile supplies it); then "$@" (targets / extra flags) bind as positionals
# because the preceding option is bounded, and a user -j/--cores later wins. The
# nargs="+" options (--configfile, --rerun-triggers) go LAST, with nothing after
# them to swallow a target as an extra value.
exec "$SNAKEMAKE" \
    --cores "${CORES:-${THREADS:-4}}" \
    "$@" \
    ${extra[@]+"${extra[@]}"} \
    --rerun-triggers mtime \
    --snakefile "$HERE/workflow/Snakefile" \
    --configfile "$CONFIG"
