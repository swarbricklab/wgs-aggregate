# wgs-aggregate

A standalone [Snakemake](https://snakemake.github.io) workflow that aggregates
per-sample WGS pipeline outputs (sash, oncoanalyser, oncoanalyser DNA+RNA, DRAGEN)
into cohort-level tables, builds the per-pipeline report indices, and assembles
the cohort germline genotype VCF. It plugs into **any dataset** — the atlas
(`brca_mega_atlas-data`), the consortium registry, or a future cohort — via a
small `deps`/`outs`/`params` config contract.

It follows the lab "a stage is a product" pattern (config-driven, schema-validated,
tested, containerisable; runs standalone or as a git submodule inside a DVC-tracked
dataset) — the same shape as [`bulk-rna-preprocessing`](https://github.com/swarbricklab/bulk-rna-preprocessing).

## Why standalone

These aggregation steps were previously duplicated as inline mini-modules inside
each dataset repo. Consolidating them into one versioned workflow gives:

- **One source of truth** for the aggregation logic (same tables everywhere;
  verified byte-identical to the prior per-module outputs).
- **No registry coupling.** Sample identity comes entirely from the dataset's
  consolidated samplesheets — there is no `index.csv` and no `config/projects/*`.
  (Garvan↔UMCCR id links like `wgs_sbj`/`wgs_prj` are a registry concern, dropped.)
- **Managed environments.** Every rule declares its conda env, so consumers don't
  hand-manage the R/bcftools toolchain.

## Steps & outputs

| step | input pipeline dir | key outputs (`<out_dir>/<step>_aggregate/…`) |
| --- | --- | --- |
| `sash` | `data/wgs/sash` | PURPLE drivers/CNV, PCGR/CPSR, signatures, smlv, cohort summary, cancer report |
| `oncoanalyser` | `data/wgs/oncoanalyser` | CUPPA, LILAC, CHORD, virus, LINX fusions/drivers/neoepitopes, immune escape |
| `oncoanalyser_dna_rna` | `data/wgs/oncoanalyser_dna_rna` | Neo scored neoepitopes, presented peptides, TNB summary |
| `dragen` | `data/wgs/tumor-normal` | DRAGEN CNV/coverage/HRD/mapping/ploidy/ROH/SV/trimmer/VC metrics |
| `report_index` | the above | per-pipeline HTML+CSV report indices (under `reports_dir`) |
| `germline` | `data/wgs/tumor-normal` | cohort germline common-SNP VCF + per-project subsets |

## Layout

```
config/config.yaml                 # the deps/outs/params config contract (template)
env/environment.yaml               # single conda env (R stack + bcftools + snakemake)
profiles/workflow/config.yaml      # per-rule resources, use-conda
workflow/
  Snakefile                        # validates config, pulls deps/outs/params, includes rules
  schemas/config.schema.yaml       # config schema, validated at workflow start
  rules/
    common.smk                     # the ONE dispatch helper (r_cmd) + path helpers + target lists
    {sash,oncoanalyser,oncoanalyser_dna_rna,dragen,report_index,germline}.smk
  scripts/
    utils.R                        # consolidated R helpers (samplesheet load, write_cohort, …)
    <step>/*.R, germline/*.{sh,py} # the aggregation logic (unchanged on extraction)
test/                              # self-contained synthetic fixture + golden + comparator
.tests/ci/                         # data-independent config + placeholder setup (lint / schema / dry-run)
run.sh  run_test.sh
```

**One dispatch convention.** Every R aggregation rule is a thin wrapper that calls
`r_cmd()` in `rules/common.smk`; that single helper builds the `WGS_*` environment
the R scripts read. (This replaces the five bespoke per-file env-builders the
module used to carry.) The rule *name* and its `output` say what it produces; the
heavy logic stays in `workflow/scripts/`. Note: no f-strings appear in the
Snakefile/`.smk` files — this snakemake + py3.12 combo mangles them in rule
directives; use `+` concatenation.

## Configure

Everything dataset-specific is in `config/config.yaml`, validated against
`workflow/schemas/config.schema.yaml` at workflow start (a bad/misnamed key fails
fast). Three namespaces:

| key | meaning |
| --- | --- |
| `deps.data_dir` | root of per-sample pipeline outputs (`<data_dir>/<pipeline>/…`) |
| `deps.samplesheets.*` | consolidated samplesheet per pipeline (identity + per-row dir key) |
| `deps.cohort_samplesheet` | cohort samplesheet (`donor_id,assay,usage`) — germline cohort set |
| `outs.out_dir` | root for cohort tables (`<out_dir>/<step>_aggregate/…`) |
| `outs.reports_dir` | report-index output root (kept at `reports/wgs`; see note below) |
| `params.steps` | which aggregation steps to run |
| `params.germline.*` | germline panel knobs: `usage_include`, `maf_min`, `autosomes` |

Samplesheet column contract: [`docs/samplesheet_contract.md`](docs/samplesheet_contract.md).

> The `report_index` R scripts currently write to `reports/wgs/` hardcoded, so
> `outs.reports_dir` should stay `reports/wgs` unless those scripts are updated.

## Run

```bash
# standalone (builds config/config.yaml's outputs; target first, then flags)
./run.sh                                   # all steps
./run.sh sash                              # one step
./run.sh germline -j16
./run.sh data/wgs/sash_aggregate/purple_drivers_somatic.csv   # explicit target

# as a submodule, point at the parent repo's config
WGS_AGGREGATE_CONFIG=/path/to/parent/config.yaml ./run.sh sash
```

`run.sh` activates the Snakemake launcher env (`snakemake_7.32.4`), uses
`--use-conda` (per-rule env `env/environment.yaml`) and the workflow profile, and
on NCI adds the global PBS profile when present (`profiles/global/nci`, a
submodule). The aggregation reads every per-sample file into memory and the
germline step merges per-donor VCFs across the whole cohort — run it via the
cluster scheduler, not a login node.

> **Offline compute nodes.** `--use-conda` must build the per-rule env on first
> run; build it once on a data-mover node so compute runs are offline-safe. Or
> drive snakemake directly without conda if the toolchain is already on `PATH`.

## Test

Two tiers, both data-independent of the atlas:

```bash
# fast: schema validation + full-DAG dry-run (no real data)
bash .tests/ci/setup_placeholders.sh
snakemake -n --snakefile workflow/Snakefile --configfile .tests/ci/config.yaml

# end-to-end golden: build from test/fixtures/, compare to test/golden/
./run_test.sh --check
```

`run_test.sh` builds the two representative tables (sash PURPLE drivers + DRAGEN
coverage) from a tiny committed synthetic cohort (`test/fixtures/`, regenerated by
`test/make_fixture.sh`) and asserts they match `test/golden/` via the NaN-aware
`test/compare_csv.py`. The golden was produced by this same code, so the test is
self-contained and proves the aggregation logic is unchanged through the wiring.

## Provenance

Extracted from `brca_mega_atlas-data` (the inline `modules/*_aggregate`,
`modules/germline_cohort_vcf`, and `scripts/build_*_report_index.R`). The
extraction drops the registry-only inputs and is validated to reproduce the prior
outputs byte-for-byte on the atlas cohort. See `CITATION.cff`.
