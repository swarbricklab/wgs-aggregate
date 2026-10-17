#!/usr/bin/env python3
"""Resolve the dataset's donor set for the germline VCF subset.

This is a single-dataset workflow, so there is one project: the atlas WGS cohort.
Its donors come from the cohort samplesheet (assay == wgs, usage in the configured
include set) — replacing the old config/projects/wgs.csv x index.csv lookup. They are
intersected with the donors actually present in the cohort germline VCF (the
donor_manifest); a donor assigned but absent from the cohort VCF (tumour-only,
germline QC fail, or not sequenced) is reported, not silently dropped.

Env:
  WGS_COHORT_SAMPLESHEET   cohort samplesheet (default config/samplesheet.csv)
  WGS_GERMLINE_USAGE       comma-separated usage verdicts to include
                           (default "analysed,warn_qc")
  WGS_PROJECT              project label for the output files (default brca_mega_atlas)

Writes under data/wgs/germline_cohort_vcf/projects/:
  project_donors.csv          project, n_donors, n_assigned_absent_from_cohort_vcf, donors
  <project>.donors.txt        comma-separated donor ids for `bcftools view -s`
"""
import csv, os, sys

REPO = os.environ.get("REPO", os.getcwd())
SAMPLESHEET = os.environ.get("WGS_COHORT_SAMPLESHEET", os.path.join(REPO, "config/samplesheet.csv"))
MANIFEST = os.path.join(REPO, "data/wgs/germline_cohort_vcf/donor_manifest.csv")
OUTDIR = os.path.join(REPO, "data/wgs/germline_cohort_vcf/projects")
PROJECT = os.environ.get("WGS_PROJECT", "brca_mega_atlas")
USAGE_INCLUDE = {u.strip() for u in os.environ.get("WGS_GERMLINE_USAGE", "analysed,warn_qc").split(",") if u.strip()}


def main():
    atlas_donors = {r["donor_id"] for r in csv.DictReader(open(SAMPLESHEET))
                    if r.get("assay") == "wgs" and r.get("usage") in USAGE_INCLUDE}
    cohort_donors = {r["donor_id"] for r in csv.DictReader(open(MANIFEST))}

    present = sorted(atlas_donors & cohort_donors, key=lambda d: (len(d), d))
    absent = sorted(atlas_donors - cohort_donors)

    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, f"{PROJECT}.donors.txt"), "w") as fh:
        fh.write(",".join(present) + "\n")
    with open(os.path.join(OUTDIR, "project_donors.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["project", "n_donors",
                                           "n_assigned_absent_from_cohort_vcf", "donors"],
                           lineterminator="\n")
        w.writeheader()
        w.writerow({"project": PROJECT, "n_donors": len(present),
                    "n_assigned_absent_from_cohort_vcf": len(absent),
                    "donors": ";".join(present)})

    if absent:
        print(f"{PROJECT}: {len(absent)} atlas WGS donor(s) not in cohort VCF "
              f"(tumour-only / germline QC fail / not sequenced): {absent}", file=sys.stderr)
    print(f"  {PROJECT}: {len(present)} donors", file=sys.stderr)


if __name__ == "__main__":
    main()
