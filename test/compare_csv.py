#!/usr/bin/env python3
"""NaN-aware CSV comparator for the golden test.

Compares two CSVs for equality of shape, columns (as a set, order-independent),
and cell values -- treating empty / NA / NaN as equal, and comparing numeric
columns with a small tolerance so platform float formatting doesn't cause
spurious diffs. Exits non-zero with a readable message on the first mismatch.

Usage: compare_csv.py GOLDEN ACTUAL
Stdlib only (no pandas) so it runs under any python.
"""
from __future__ import annotations

import csv
import sys


def _load(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def _is_blank(v):
    return v is None or str(v).strip().lower() in ("", "na", "nan")


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _cell_eq(a, b, tol=1e-9):
    if _is_blank(a) and _is_blank(b):
        return True
    na, nb = _num(a), _num(b)
    if na is not None and nb is not None:
        return abs(na - nb) <= tol * max(1.0, abs(na), abs(nb))
    return str(a).strip() == str(b).strip()


def main(golden, actual):
    gh, gr = _load(golden)
    ah, ar = _load(actual)

    if set(gh) != set(ah):
        sys.exit(f"FAIL {actual}: columns differ\n  golden: {gh}\n  actual: {ah}")
    if len(gr) != len(ar):
        sys.exit(f"FAIL {actual}: row count {len(ar)} != golden {len(gr)}")

    gi = {c: i for i, c in enumerate(gh)}
    ai = {c: i for i, c in enumerate(ah)}
    for r, (grow, arow) in enumerate(zip(gr, ar), start=2):
        for col in gh:
            a, b = grow[gi[col]], arow[ai[col]]
            if not _cell_eq(a, b):
                sys.exit(f"FAIL {actual}: row {r} col '{col}': golden={a!r} actual={b!r}")
    print(f"OK {actual}  ({len(ar)} rows x {len(ah)} cols)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: compare_csv.py GOLDEN ACTUAL")
    main(sys.argv[1], sys.argv[2])
