#!/usr/bin/env Rscript
# Aggregate per-tumour PURPLE somatic driver catalogs into a cohort CSV.
# Input:  <src_dir>/purple/<tumor_lib>.purple.driver.catalog.somatic.tsv
# Output: data/wgs/sash_aggregate/purple_drivers_somatic.csv
#
# Run from project root:
#   Rscript modules/sash_aggregate/scripts/agg_drivers.R

source(Sys.getenv("WGS_UTILS", "workflow/scripts/utils.R"))

ss <- load_samplesheet()

# Numeric columns from PURPLE driver catalog. Counts (missense/nonsense/...) are
# integers; likelihood and CN are doubles. Coerce via as.numeric and let R store
# integers as doubles in the output - the loader's schema will type them.
NUMERIC_COLS <- c("driver_likelihood", "missense", "nonsense", "splice",
                  "inframe", "frameshift", "min_copy_number", "max_copy_number")
LOGICAL_COLS <- c("is_canonical", "biallelic")

# Mapping from PURPLE camelCase -> snake_case (explicit so columns survive
# upstream schema additions cleanly).
COL_MAP <- c(
  chromosome        = "chromosome",
  chromosomeBand    = "chromosome_band",
  gene              = "gene",
  transcript        = "transcript",
  isCanonical       = "is_canonical",
  driver            = "driver",
  category          = "category",
  likelihoodMethod  = "likelihood_method",
  driverLikelihood  = "driver_likelihood",
  missense          = "missense",
  nonsense          = "nonsense",
  splice            = "splice",
  inframe           = "inframe",
  frameshift        = "frameshift",
  biallelic         = "biallelic",
  minCopyNumber     = "min_copy_number",
  maxCopyNumber     = "max_copy_number"
)

read_driver_tsv <- function(src_dir, tumor_lib) {
  p <- path(src_dir, "purple",
            paste0(tumor_lib, ".purple.driver.catalog.somatic.tsv"))
  if (!file.exists(p)) return(tibble())
  df <- read_tsv_quiet(p, col_types = cols(.default = "c"))
  if (nrow(df) == 0) return(tibble())

  present <- intersect(names(COL_MAP), names(df))
  df <- df[, present, drop = FALSE]
  names(df) <- unname(COL_MAP[present])
  df
}

rows <- ss %>%
  rowwise() %>%
  reframe(
    pick(any_of(IDENTITY_COLS)),
    read_driver_tsv(src_dir, tumor_lib)
  )

if (nrow(rows) > 0) {
  num_present <- intersect(NUMERIC_COLS, names(rows))
  rows <- rows %>%
    mutate(across(all_of(num_present), ~ suppressWarnings(as.numeric(.))))

  log_present <- intersect(LOGICAL_COLS, names(rows))
  rows <- rows %>%
    mutate(across(all_of(log_present), ~ tolower(.) == "true"))
}

dir_create(OUT_DIR)
out_csv <- file.path(OUT_DIR, "purple_drivers_somatic.csv")
write_csv(rows, out_csv, na = "NA")
cat(sprintf("Wrote %d driver row(s) from %d tumour(s) to %s\n",
            nrow(rows), length(unique(rows$tumor_gussid)), out_csv))
