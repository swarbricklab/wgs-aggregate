#!/usr/bin/env Rscript

# Build a searchable HTML index of DRAGEN tumor-normal multiqc reports.
#
# Outputs two companion files:
#   reports/wgs/tumor-normal.csv  - one row per dragen run with identity + report
#                               path columns. Standalone-useful (grep/awk OK).
#   reports/wgs/tumor-normal.html - DataTables-powered table that fetches the CSV
#                               at load time. Requires being served over HTTP
#                               (LiveServer or `python -m http.server`); the
#                               browser blocks fetch() from file:// URLs.
#
# Identity comes from joining config/tumor-normal/dragen_runs.csv with
# index.csv (via wgs_lib) to surface tumor/normal gussid + tissue_id.
#
# Usage:
#   Rscript scripts/build_tumor_normal_report_index.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fs)
})

DRAGEN_RUNS <- "config/tumor-normal/dragen_runs.csv"
SAMPLESHEET <- Sys.getenv("WGS_SASH_SAMPLESHEET", "config/sash/samplesheet.csv")
REPORTS_DIR <- "reports/wgs/tumor-normal"
OUT_CSV     <- "reports/wgs/tumor-normal.csv"
OUT_HTML    <- "reports/wgs/tumor-normal.html"

runs <- read_csv(DRAGEN_RUNS, show_col_types = FALSE, col_types = cols(.default = "c"))
ss   <- read_csv(SAMPLESHEET,  show_col_types = FALSE, col_types = cols(.default = "c"))

# identity from the sash samplesheet (no registry index.csv)
tumour_meta <- ss %>%
  transmute(tumor_lib, tumor_gussid, tumor_tissue_id) %>%
  distinct()

normal_meta <- ss %>%
  transmute(normal_lib, normal_gussid, normal_tissue_id) %>%
  distinct()

# multiqc filename: {tumor_lib}__{normal_lib}_dragen_somatic_and_germline_multiqc.html
multiqc_filename <- function(tumor_lib, normal_lib) {
  sprintf("%s__%s_dragen_somatic_and_germline_multiqc.html", tumor_lib, normal_lib)
}

cohort <- runs %>%
  rename(tumor_lib = tumor_library, normal_lib = germline_library) %>%
  left_join(tumour_meta, by = "tumor_lib") %>%
  left_join(normal_meta, by = "normal_lib") %>%
  mutate(
    multiqc_basename = multiqc_filename(tumor_lib, normal_lib),
    multiqc          = file.path("tumor-normal", run_id, multiqc_basename),
    multiqc_exists   = file_exists(file.path(REPORTS_DIR, run_id, multiqc_basename))
  )

n_missing <- sum(!cohort$multiqc_exists)
if (n_missing > 0) {
  message(sprintf("WARN: %d run(s) have no multiqc HTML in %s — leaving blank:",
                  n_missing, REPORTS_DIR))
  print(cohort %>% filter(!multiqc_exists) %>% select(run_id, tumor_lib, normal_lib))
}

cohort <- cohort %>%
  mutate(multiqc = if_else(multiqc_exists, multiqc, NA_character_)) %>%
  select(
    run_id, donor_id,
    tumor_gussid, tumor_tissue_id, tumor_lib,
    normal_gussid, normal_tissue_id, normal_lib,
    multiqc
  )

dir_create(path_dir(OUT_CSV))
write_csv(cohort, OUT_CSV, na = "")
cat(sprintf("Wrote %d rows to %s\n", nrow(cohort), OUT_CSV))


# ---- HTML (CSV-driven via fetch + PapaParse + DataTables) -------------------

html <- '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>DRAGEN tumor-normal report index</title>
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
<h1>DRAGEN tumor-normal report index</h1>
<p class="subtitle">Search any column (donor, tissue, gussid, library, run_id).</p>
<p id="status">Loading…</p>
<table id="cohort" class="display compact stripe" style="width:100%"></table>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/2.1.8/js/dataTables.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/papaparse@5.4.1/papaparse.min.js"></script>
<script>
const CSV_URL = "tumor-normal.csv";
const ACRONYMS = { id: "ID", multiqc: "MultiQC", pcgr: "PCGR", cpsr: "CPSR", cacao: "CACAO" };

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
