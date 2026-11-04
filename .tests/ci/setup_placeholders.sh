#!/bin/bash
# Touch empty placeholder inputs so `snakemake -n` can build the full DAG in CI
# without any real data (the R scripts never run in a dry-run, so empty files are
# enough). Pair with .tests/ci/config.yaml.
set -eu
base=.tests/ci
mkdir -p "$base/data"
touch "$base/sash.samplesheet.csv" \
      "$base/oncoanalyser.samplesheet.csv" \
      "$base/oncoanalyser_dna_rna.samplesheet.csv" \
      "$base/dragen_runs.csv" \
      "$base/cohort.samplesheet.csv"
echo "placeholders touched"
