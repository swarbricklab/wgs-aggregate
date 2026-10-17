#!/usr/bin/env Rscript

# Build a searchable HTML index of sash per-sample reports.
#
# Outputs two companion files:
#   reports/wgs/sash.csv   - one row per lib_pair (within a sash run), with
#                        identity columns + 6 report path columns. Standalone-
#                        useful (grep/awk OK).
#   reports/wgs/sash.html  - DataTables-powered table that fetches the CSV at
#                        load time via PapaParse. Must be served over HTTP
#                        (LiveServer or `python -m http.server`).
#
# Identity comes from parsing the lib_pair dir name (<tumor_lib>__<normal_lib>)
# and joining with index.csv via wgs_lib to surface tumor/normal gussid +
# tissue_id, donor_id.
#
# Usage:
#   Rscript scripts/build_sash_report_index.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(fs)
})

SAMPLESHEET <- Sys.getenv("WGS_SASH_SAMPLESHEET", "config/sash/samplesheet.csv")
REPORTS_DIR <- "reports/wgs/sash"
OUT_CSV     <- "reports/wgs/sash.csv"
OUT_HTML    <- "reports/wgs/sash.html"

# (filename suffix, CSV column name). Note the linx uses _underscore, others .dot.
REPORT_SPECS <- list(
  list(suffix = ".cancer_report.html",          col = "cancer_report", lib = "tumor"),
  list(suffix = ".multiqc.html",                col = "multiqc",       lib = "tumor"),
  list(suffix = ".pcgr.html",                   col = "pcgr",          lib = "tumor"),
  list(suffix = ".pcgr_acmg.grch38.flexdb.html",col = "pcgr_flexdb",   lib = "tumor"),
  list(suffix = "_linx.html",                   col = "linx",          lib = "tumor"),
  list(suffix = ".cpsr.html",                   col = "cpsr",          lib = "normal")
)

# ---- discover lib_pair dirs under reports/wgs/sash/<run>/ -----------------------

run_dirs <- dir_ls(REPORTS_DIR, type = "directory")
pairs <- lapply(run_dirs, function(rd) {
  lp <- dir_ls(rd, type = "directory")
  if (length(lp) == 0) return(NULL)
  tibble(
    sash_run_id = basename(rd),
    lib_pair    = basename(lp)
  )
}) |> bind_rows()

# lib_pair = "<tumor_lib>_<normal_lib>" or "<tumor_lib>__<normal_lib>"; both
# conventions appear in the cohort. Library IDs are L<digits> with no internal
# underscores, so collapsing 1+ underscores is safe.
pairs <- pairs |>
  separate(lib_pair, into = c("tumor_lib", "normal_lib"), sep = "_+", remove = FALSE)

cat(sprintf("Discovered %d lib_pair entries across %d sash runs\n",
            nrow(pairs), length(unique(pairs$sash_run_id))))

# ---- join identity from the sash samplesheet (no registry index.csv) --------

ss <- read_csv(SAMPLESHEET, show_col_types = FALSE, col_types = cols(.default = "c"))

tumour_meta <- ss |>
  transmute(tumor_lib, tumor_gussid, donor_id, tumor_tissue_id) |>
  distinct()

normal_meta <- ss |>
  transmute(normal_lib, normal_gussid, normal_tissue_id) |>
  distinct()

cohort <- pairs |>
  left_join(tumour_meta, by = "tumor_lib") |>
  left_join(normal_meta, by = "normal_lib")

# ---- build per-report relative paths ----------------------------------------

make_path <- function(sash_run_id, lib_pair, lib, suffix) {
  filename  <- paste0(lib, suffix)
  rel_path  <- file.path("sash", sash_run_id, lib_pair, filename)
  full_path <- file.path(REPORTS_DIR, sash_run_id, lib_pair, filename)
  if_else(file_exists(full_path), rel_path, NA_character_)
}

for (spec in REPORT_SPECS) {
  lib_col <- if (spec$lib == "tumor") "tumor_lib" else "normal_lib"
  cohort[[spec$col]] <- mapply(
    function(r, lp, l) make_path(r, lp, l, spec$suffix),
    cohort$sash_run_id, cohort$lib_pair, cohort[[lib_col]],
    USE.NAMES = FALSE
  )
}

# ---- emit CSV ---------------------------------------------------------------

out <- cohort |>
  select(
    sash_run_id, donor_id,
    tumor_gussid, tumor_tissue_id, tumor_lib,
    normal_gussid, normal_tissue_id, normal_lib,
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
<title>sash cohort report index</title>
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
<h1>sash cohort report index</h1>
<p class="subtitle">Search any column (donor, tissue, gussid, library, sash_run_id).</p>
<p id="status">Loading…</p>
<table id="cohort" class="display compact stripe" style="width:100%"></table>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/2.1.8/js/dataTables.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/papaparse@5.4.1/papaparse.min.js"></script>
<script>
const CSV_URL = "sash.csv";
const ACRONYMS = { id: "ID", multiqc: "MultiQC", pcgr: "PCGR", cpsr: "CPSR", cacao: "CACAO", linx: "LINX", flexdb: "FlexDB" };

function looksLikePath(v) {
  return typeof v === "string" && v.indexOf("/") >= 0 && v.toLowerCase().endsWith(".html");
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
