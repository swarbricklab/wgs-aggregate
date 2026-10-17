# wgs-aggregate

A standalone **Snakemake** workflow that aggregates per-sample WGS pipeline
outputs (sash, oncoanalyser, oncoanalyser DNA+RNA, DRAGEN) into cohort-level
tables, builds the per-pipeline report indices, and assembles the cohort germline
genotype VCF. It plugs into **any dataset** — the atlas (`brca_mega_atlas-data`),
the consortium registry, or a future cohort — via a small config contract.

## Why standalone

These aggregation steps were previously duplicated as inline mini-modules inside
each dataset repo. Consolidating them into one versioned Snakemake workflow:

- **One source of truth** for the aggregation logic (same tables everywhere;
  verified byte-identical to the prior per-module outputs).
- **No registry coupling.** Sample identity comes entirely from the dataset's
  consolidated samplesheets — there is no `index.csv` and no
  `config/projects/*.csv`. (Garvan↔UMCCR id links like `wgs_sbj`/`wgs_prj` are a
  registry concern, not a dataset one, and are dropped.)
- **Managed environments.** Each rule declares its conda env, so consumers don't
  hand-manage R/bcftools toolchains.

## Usage

Add as a submodule and run against the dataset's config:

```bash
snakemake --use-conda --configfile config/config.yaml -j<N>
# or a subset:
snakemake --use-conda sash oncoanalyser -j<N>
```

## The contract (`config/config.yaml`)

| key | meaning |
| --- | --- |
| `data_dir` | root of per-sample pipeline outputs (`<data_dir>/<pipeline>/…`) |
| `out_dir` | root for cohort tables (`<out_dir>/<step>_aggregate/…`) |
| `samplesheets.*` | consolidated samplesheet per pipeline (identity + per-row dir key) |
| `cohort_samplesheet` | atlas/cohort samplesheet (`donor_id,assay,usage`) — germline cohort set |
| `germline_usage_include` | usage verdicts admitted to the germline cohort |
| `steps` | which aggregation steps to run |

Samplesheet column contract: `docs/samplesheet_contract.md`.

## Steps

| step | input pipeline dir | key outputs |
| --- | --- | --- |
| `sash` | `data/wgs/sash` | PURPLE drivers/CNV, PCGR/CPSR, signatures, smlv, cohort summary, cancer report |
| `oncoanalyser` | `data/wgs/oncoanalyser` | CUPPA, LILAC, CHORD, virus, LINX fusions/drivers/neoepitopes, immune escape |
| `oncoanalyser_dna_rna` | `data/wgs/oncoanalyser_dna_rna` | Neo scored neoepitopes, presented peptides, TNB summary |
| `dragen` | `data/wgs/tumor-normal` | DRAGEN CNV/coverage/HRD/mapping/ploidy/ROH/SV/trimmer/VC metrics |
| `report_index` | the above | per-pipeline HTML report indices (identity from samplesheets) |
| `germline` | `data/wgs/tumor-normal` | cohort germline common-SNP VCF + per-cohort subset |

## Provenance

Extracted from `brca_mega_atlas-data` (the inline `modules/*_aggregate`,
`modules/germline_cohort_vcf`, and `scripts/build_*_report_index.R`). The extraction
drops the registry-only inputs and is validated to reproduce the prior outputs
byte-for-byte on the atlas cohort.
