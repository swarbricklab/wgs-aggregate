#!/usr/bin/env python3
"""Build the donor -> DRAGEN germline VCF manifest for the cohort germline panel.

One row per donor that has a usable genome-wide germline VCF. The germline library
for each donor is the germline_library from dragen_runs.csv (no registry index.csv) —
dragen_runs enumerates the DRAGEN runs, so it is the source of truth for the VCF
locations; scoping the germline cohort to the atlas therefore flows through the
dragen_runs atlas subset, consistent with dragen_aggregate. Its DRAGEN hard-filtered
VCF is located on disk under
  data/wgs/tumor-normal/<run>/<lib>_dragen_germline/<lib>.hard-filtered.vcf.gz
A handful of libraries have more than one run on disk (re-runs); the run with the
most recent DATE prefix is kept. When several runs share that latest date (the
trailing run id is opaque, not a recency signal) the choice is made deterministic
(lexical-max run id) and flagged via ambiguous_same_date + a stderr warning, because
such same-date re-runs can differ in content (see bcarc_wgs#6).

Emits data/wgs/germline_cohort_vcf/donor_manifest.csv with provenance columns so the
cohort VCF's sample<->source mapping is auditable. Donors whose germline library has
no VCF (the 'Fail' placeholder lib, or a not-yet-fetched run) are reported and
skipped, not silently dropped.

Registry-only columns (wgs_sbj, compartment) are intentionally not carried — the
dataset samplesheet is the sole identity source.
"""
import csv, glob, os, sys, collections

REPO = os.environ.get("REPO", os.getcwd())
DRAGEN_RUNS = os.environ.get("WGS_DRAGEN_RUNS", os.path.join(REPO, "config/tumor-normal/dragen_runs.csv"))
SASH_SS = os.environ.get("WGS_SASH_SAMPLESHEET", os.path.join(REPO, "config/sash/samplesheet.csv"))
TN_GLOB = os.path.join(REPO, "data/wgs/tumor-normal/*/{lib}_dragen_germline/{lib}.hard-filtered.vcf.gz")
OUT = os.path.join(REPO, "data/wgs/germline_cohort_vcf/donor_manifest.csv")

INVALID_LIBS = {"", "Fail", "NA", "na"}


def _run_of(p):
    # run dir is the parent-of-parent: .../tumor-normal/<run>/<lib>_dragen_germline/<file>
    return os.path.basename(os.path.dirname(os.path.dirname(p)))


def _date_of(run):
    # run dirs are "<YYYYMMDD><opaque-hash>"; recency lives in the date prefix.
    return run[:8] if run[:8].isdigit() else ""


def find_vcfs(lib):
    hits = glob.glob(TN_GLOB.format(lib=lib))
    return sorted(hits, key=lambda p: (_date_of(_run_of(p)), _run_of(p)))


def _tissue_by_lib():
    """normal_lib -> normal_tissue_id from the sash samplesheet (best-effort enrichment)."""
    try:
        return {(r.get("normal_lib") or "").strip(): (r.get("normal_tissue_id") or "")
                for r in csv.DictReader(open(SASH_SS))}
    except FileNotFoundError:
        return {}


def germline_libs():
    """Distinct (donor_id, germline_library) from dragen_runs.csv; tissue_id enriched
    from the sash samplesheet where available."""
    tissue = _tissue_by_lib()
    seen = {}
    for r in csv.DictReader(open(DRAGEN_RUNS)):
        lib = (r.get("germline_library") or "").strip()
        if lib in INVALID_LIBS:
            continue
        key = (r["donor_id"], lib)
        if key not in seen:
            seen[key] = {"donor_id": r["donor_id"], "wgs_lib": lib,
                         "tissue_id": tissue.get(lib, "")}
    return list(seen.values())


def main():
    germ = germline_libs()
    by_donor = collections.Counter(r["donor_id"] for r in germ)
    multi = {d: n for d, n in by_donor.items() if n > 1}
    if multi:
        print(f"WARNING: {len(multi)} donor(s) have >1 germline library; keeping all "
              f"as separate rows: {multi}", file=sys.stderr)

    rows, skipped, ambiguous = [], [], []
    for r in germ:
        lib = r["wgs_lib"]
        vcfs = find_vcfs(lib)
        if not vcfs:
            skipped.append((r["donor_id"], lib))
            continue
        latest_date = _date_of(_run_of(vcfs[-1]))
        same_date = [p for p in vcfs if _date_of(_run_of(p)) == latest_date]
        chosen = max(same_date, key=_run_of)
        run = _run_of(chosen)
        is_ambiguous = len(same_date) > 1
        if is_ambiguous:
            ambiguous.append((lib, latest_date, sorted(_run_of(p) for p in same_date)))
        rows.append({
            "donor_id": r["donor_id"],
            "wgs_lib": lib,
            "tissue_id": r["tissue_id"],
            "source_run": run,
            "n_runs_on_disk": len(vcfs),
            "ambiguous_same_date": int(is_ambiguous),
            "source_vcf": os.path.relpath(chosen, REPO),
        })

    rows.sort(key=lambda x: x["donor_id"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"donors with a germline VCF: {len(rows)}", file=sys.stderr)
    if skipped:
        print(f"donors skipped (germline lib has no VCF): {len(skipped)} -> {skipped}",
              file=sys.stderr)
    if ambiguous:
        print(f"WARNING: {len(ambiguous)} donor(s) have >1 germline run sharing the latest "
              f"date; kept the lexical-max run id deterministically "
              f"(ambiguous_same_date=1; may differ in content, see bcarc_wgs#6):",
              file=sys.stderr)
        for lib, d, runs in ambiguous:
            print(f"  {lib} {d}: {runs}", file=sys.stderr)
    print(OUT)


if __name__ == "__main__":
    main()
