#!/usr/bin/env bash
# Run the wgs-aggregate workflow for the given step(s)/target(s), invoked from a
# dataset root (where config/ + data/ live):
#   bash modules/wgs-aggregate/run_mod.sh <step|target> ...
#   e.g.  bash modules/wgs-aggregate/run_mod.sh sash
#         THREADS=16 bash modules/wgs-aggregate/run_mod.sh germline
#
# Resolves its own Snakefile; the dataset supplies the config (default
# config/wgs-aggregate.yaml, override with WGS_CONFIG=).
#
# Env (NCI defaults; override for another cluster): the R aggregators run under the
# `meta` conda env, bcftools (sash smlv + germline) from the environment module, and
# snakemake is invoked by full path from its own env — so no --use-conda is needed
# and this works on internet-less compute nodes. For full portability, switch to
# `--use-conda` with workflow/envs/ instead of the env activation below.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${WGS_CONFIG:-config/wgs-aggregate.yaml}"

source /g/data/a56/conda/miniconda/etc/profile.d/conda.sh
conda activate "${WGS_CONDA_ENV:-meta}"
module load bcftools/1.21 2>/dev/null || true

SNAKEMAKE="${SNAKEMAKE:-/g/data/a56/conda/envs/snakemake_7.32.4/bin/snakemake}"
# --rerun-triggers is greedy (nargs=+); keep it before -j so an option (not the
# target args in "$@") terminates its value list.
exec "$SNAKEMAKE" \
  --snakefile "$HERE/workflow/Snakefile" \
  --configfile "$CONFIG" \
  --rerun-triggers mtime \
  -j "${THREADS:-4}" \
  "$@"
