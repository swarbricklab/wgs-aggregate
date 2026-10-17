#!/usr/bin/env Rscript

# Build a searchable HTML index of oncoanalyser DNA+RNA per-sample reports.
#
# Outputs:
#   reports/wgs/oncoanalyser_dna_rna.csv   - one row per dna-rna run, identity +
#                                        the orange (RNA-integrated) PDF path.
#   reports/wgs/oncoanalyser_dna_rna.html  - DataTables table fetching the CSV via
#                                        PapaParse (serve over HTTP).
#
# These RNA-append runs publish only the ORANGE report (no LINX html); the
# DNA-only LINX/cancer reports are in reports/wgs/oncoanalyser/. Identity is joined
# from config/oncoanalyser_dna_rna/samplesheet.csv by dna_rna_run_id.
#
# Usage:
#   Rscript scripts/build_oncoanalyser_dna_rna_report_index.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fs)
})

MANIFEST    <- "config/oncoanalyser_dna_rna/samplesheet.csv"
REPORTS_DIR <- "reports/wgs/oncoanalyser_dna_rna"
OUT_CSV     <- "reports/wgs/oncoanalyser_dna_rna.csv"
OUT_HTML    <- "reports/wgs/oncoanalyser_dna_rna.html"

REPORT_SPECS <- list(
  list(suffix = ".orange.pdf", col = "orange", lib = "tumor")
)

# ---- discover lib_dir dirs under reports/wgs/oncoanalyser_dna_rna/<run>/ --------

run_dirs <- dir_ls(REPORTS_DIR, type = "directory")
pairs <- lapply(run_dirs, function(rd) {
  ld <- dir_ls(rd, type = "directory")
  if (length(ld) == 0) return(NULL)
  tibble(dna_rna_run_id = basename(rd), lib_dir = basename(ld))
}) |> bind_rows()

cat(sprintf("Discovered %d lib_dir entries across %d dna-rna runs\n",
            nrow(pairs), length(unique(pairs$dna_rna_run_id))))

# ---- join identity from the master manifest ---------------------------------

manifest <- read_csv(MANIFEST, show_col_types = FALSE, col_types = cols(.default = "c"))

cohort <- pairs |>
  left_join(
    manifest |> select(dna_rna_run_id, sash_run_id, donor_id,
                        tumor_gussid, tumor_tissue_id, tumor_lib,
                        normal_gussid, normal_tissue_id, normal_lib, rna_lib),
    by = "dna_rna_run_id"
  )

# ---- build per-report relative paths ----------------------------------------

make_path <- function(run, lib_dir, lib, suffix) {
  filename  <- paste0(lib, suffix)
  rel_path  <- file.path("oncoanalyser_dna_rna", run, lib_dir, filename)
  full_path <- file.path(REPORTS_DIR, run, lib_dir, filename)
  if_else(file_exists(full_path), rel_path, NA_character_)
}

for (spec in REPORT_SPECS) {
  lib_col <- if (spec$lib == "tumor") "tumor_lib" else "normal_lib"
  cohort[[spec$col]] <- mapply(
    function(r, ld, l) make_path(r, ld, l, spec$suffix),
    cohort$dna_rna_run_id, cohort$lib_dir, cohort[[lib_col]],
    USE.NAMES = FALSE
  )
}

# ---- emit CSV ---------------------------------------------------------------

out <- cohort |>
  select(
    dna_rna_run_id, sash_run_id, donor_id,
    tumor_gussid, tumor_tissue_id, tumor_lib,
    normal_gussid, normal_tissue_id, normal_lib, rna_lib,
    all_of(vapply(REPORT_SPECS, function(s) s$col, character(1)))
  )

dir_create(path_dir(OUT_CSV))
write_csv(out, OUT_CSV, na = "")
cat(sprintf("Wrote %d rows to %s\n", nrow(out), OUT_CSV))

# ---- HTML (CSV-driven via fetch + PapaParse + DataTables) -------------------

html <- '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>oncoanalyser DNA+RNA cohort report index</title>
<link rel="stylesheet" href="https://cdn.datatables.net/2.1.8/css/dataTables.dataTables.min.css">
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 24px; color: #1f2328; }
  h1 { margin: 0 0 4px 0; font-size: 1.5em; }
  .subtitle { color: #57606a; margin: 0 0 16px 0; font-size: 0.95em; }
  #status { color: #57606a; font-size: 0.9em; margin-top: 10px; }
  #status.error { color: #b30000; }
  table.dataTable td { vertical-align: middle; }
  table.dataTable td a {
    display: inline-block; padding: 2px 8px; margin: 1px;
    border-radius: 3px; background: #eef3fa; color: #0969da;
    text-decoration: none; font-size: 0.85em; white-space: nowrap;
  }
  table.dataTable td a:hover { background: #d8e6f5; }
  table.dataTable td .missing {
    display: inline-block; padding: 2px 8px; margin: 1px;
    color: #ccc; font-size: 0.85em; white-space: nowrap;
  }
  table.dataTable thead th { background: #f6f8fa; font-weight: 600; }
  .dt-search input { width: 320px; }
</style>
</head>
<body>
<h1>oncoanalyser DNA+RNA cohort report index</h1>
<p class="subtitle">RNA-integrated ORANGE reports. Search any column (donor, tissue, gussid, library, rna_lib, dna_rna_run_id).</p>
<p id="status">Loading…</p>
<table id="cohort" class="display compact stripe" style="width:100%"></table>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/2.1.8/js/dataTables.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/papaparse@5.4.1/papaparse.min.js"></script>
<script>
const CSV_URL = "oncoanalyser_dna_rna.csv";
const ACRONYMS = { id: "ID", orange: "ORANGE", rna: "RNA", gussid: "gussid" };

function looksLikePath(v) {
  if (typeof v !== "string" || v.indexOf("/") < 0) return false;
  const lower = v.toLowerCase();
  return lower.endsWith(".html") || lower.endsWith(".pdf");
}

function isPathColumn(rows, field) {
  return rows.some(function (r) { return looksLikePath(r[field]); });
}

function prettyLabel(field) {
  return field.split(/[_\\s]+/).map(function (w) {
    const lower = w.toLowerCase();
    return ACRONYMS[lower] || (w.charAt(0).toUpperCase() + w.slice(1));
  }).join(" ");
}

function renderCell(val, field, isPath) {
  if (val == null || val === "") {
    return isPath ? "<span class=\\"missing\\">" + prettyLabel(field) + "</span>" : "";
  }
  if (isPath) {
    return "<a href=\\"" + val + "\\" target=\\"_blank\\" rel=\\"noopener\\">" + prettyLabel(field) + "</a>";
  }
  return $("<div>").text(val).html();
}

fetch(CSV_URL)
  .then(function (resp) {
    if (!resp.ok) throw new Error("HTTP " + resp.status + " fetching " + CSV_URL);
    return resp.text();
  })
  .then(function (text) {
    const parsed = Papa.parse(text, { header: true, skipEmptyLines: true });
    if (parsed.errors && parsed.errors.length) {
      console.warn("CSV parse warnings:", parsed.errors);
    }
    const fields = parsed.meta.fields;
    const rows   = parsed.data;
    const pathCols = {};
    fields.forEach(function (f) { pathCols[f] = isPathColumn(rows, f); });

    const headerHtml = "<thead><tr>" +
      fields.map(function (f) { return "<th>" + prettyLabel(f) + "</th>"; }).join("") +
      "</tr></thead>";
    const bodyHtml = "<tbody>" + rows.map(function (r) {
      return "<tr>" + fields.map(function (f) {
        return "<td>" + renderCell(r[f], f, pathCols[f]) + "</td>";
      }).join("") + "</tr>";
    }).join("") + "</tbody>";

    $("#cohort").html(headerHtml + bodyHtml);
    const pathColIdx = fields.map(function (f, i) { return pathCols[f] ? i : -1; })
                              .filter(function (i) { return i >= 0; });
    $("#cohort").DataTable({
      pageLength: 25,
      lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
      order: [[0, "asc"]],
      columnDefs: [{ targets: pathColIdx, orderable: false, searchable: false }]
    });
    $("#status").text(rows.length + " rows loaded.");
  })
  .catch(function (err) {
    console.error(err);
    $("#status").addClass("error").text(
      "Failed to load " + CSV_URL + ": " + err.message +
      " (note: this page must be served over HTTP, not opened as file://)"
    );
  });
</script>
</body>
</html>
'

writeLines(html, OUT_HTML)
cat(sprintf("Wrote %s (%.1f KB)\n", OUT_HTML, file.size(OUT_HTML) / 1024))
