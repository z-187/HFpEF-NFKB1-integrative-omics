############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 1 FIXED v2: Complete input-data audit and manifest
##
## FIXED v2 change:
##   - Renames the GEO metadata output column from `key` to `metadata_key`.
##     In data.table(), `key` is a reserved constructor argument; using it as
##     a column name caused setkeyv() to interpret GEO header text as column names.
##
## Project root:
##   <HFPEF_PROJECT_DIR>
##
## Read-only input directory:
##   <HFPEF_PROJECT_DIR>/0.GEO
##
## Output directory:
##   <HFPEF_PROJECT_DIR>/01_stage1_data_audit_FIXED_v2
##
## Review package:
##   <HFPEF_PROJECT_DIR>/01_stage1_data_audit_FIXED_v2_CHECK.zip
##
## Purpose:
##   1) Inventory every local input file.
##   2) Calculate file hashes for reproducibility.
##   3) Parse GEO series/family metadata without analysing biology.
##   4) Inspect archive contents, 10x matrix dimensions, H5AD structure,
##      and NicheNet RDS object structure.
##   5) Build a dataset-role manifest and required-input checklist.
##   6) Produce a compact CHECK.zip for external review.
##
## Scientific boundary:
##   - This stage performs data auditing only.
##   - It does not perform normalization, differential expression,
##     clustering, cell annotation, pathway analysis, or candidate ranking.
##   - No biological conclusion should be drawn from Stage 1 outputs.
##
## Recommended run:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage1_Data_Audit_FIXED_v2.R",
##     encoding = "UTF-8"
##   )
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(warn = 1)
options(encoding = "UTF-8")

############################################################
## 0. User settings
############################################################

PROJECT_DIR <- Sys.getenv("HFPEF_PROJECT_DIR", unset = "")
if (!nzchar(PROJECT_DIR)) {
  stop(
    "HFPEF_PROJECT_DIR is not set. Define it as the local project root ",
    "containing 0.GEO and the stage output folders before running this script."
  )
}
PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)
DATA_DIR    <- file.path(PROJECT_DIR, "0.GEO")

STAGE_NAME  <- "01_stage1_data_audit_FIXED_v2"
OUT_DIR     <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP   <- file.path(PROJECT_DIR, paste0(STAGE_NAME, "_CHECK.zip"))

## Protect existing results. Set TRUE only when intentionally rerunning
## this exact version after manually backing up the old output.
ALLOW_OVERWRITE <- FALSE

## Reproducibility settings.
COMPUTE_SHA256       <- TRUE
HASH_MAX_FILE_GB     <- 20
DEEP_INSPECT_MAX_GB  <- 5
COUNT_LINES_MAX_GB   <- 0.75
MAX_GEO_HEADER_LINES <- 250000L

## Optional deep inspection of H5AD and RDS resources.
INSPECT_H5AD <- TRUE
INSPECT_RDS  <- TRUE

############################################################
## 1. Output directories and package setup
############################################################

if (!dir.exists(PROJECT_DIR)) {
  stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
}
if (!dir.exists(DATA_DIR)) {
  stop("DATA_DIR does not exist: ", DATA_DIR)
}
if (dir.exists(OUT_DIR) && !ALLOW_OVERWRITE) {
  stop(
    "Output folder already exists and overwrite protection is active:\n",
    OUT_DIR,
    "\nUse a new version name or set ALLOW_OVERWRITE <- TRUE after backup."
  )
}
if (dir.exists(OUT_DIR) && ALLOW_OVERWRITE) {
  unlink(OUT_DIR, recursive = TRUE, force = TRUE)
}
if (file.exists(CHECK_ZIP) && !ALLOW_OVERWRITE) {
  stop(
    "CHECK zip already exists and overwrite protection is active:\n",
    CHECK_ZIP
  )
}
if (file.exists(CHECK_ZIP) && ALLOW_OVERWRITE) {
  unlink(CHECK_ZIP, force = TRUE)
}

DIRS <- list(
  logs       = file.path(OUT_DIR, "00_logs"),
  tables     = file.path(OUT_DIR, "01_tables"),
  objects    = file.path(OUT_DIR, "02_objects"),
  figures    = file.path(OUT_DIR, "03_figures"),
  source     = file.path(OUT_DIR, "04_source_data"),
  methods    = file.path(OUT_DIR, "05_methods"),
  check      = file.path(OUT_DIR, "06_review_check")
)

for (d in c(OUT_DIR, unlist(DIRS, use.names = FALSE))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(DIRS$logs, "stage1_data_audit.log")
WARN_FILE <- file.path(DIRS$logs, "stage1_warnings.log")

log_msg <- function(..., level = "INFO") {
  txt <- paste0(..., collapse = "")
  line <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, txt)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
  invisible(line)
}

warning_records <- list()
add_warning <- function(category, item, message) {
  rec <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    category = as.character(category),
    item = as.character(item),
    message = as.character(message),
    stringsAsFactors = FALSE
  )
  warning_records[[length(warning_records) + 1L]] <<- rec
  cat(sprintf("[%s] [%s] %s: %s\n", rec$timestamp, category, item, message),
      file = WARN_FILE, append = TRUE)
  log_msg(category, " | ", item, " | ", message, level = "WARN")
  invisible(rec)
}

ensure_cran <- function(pkgs, required = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    log_msg("Installing missing CRAN package(s): ", paste(missing, collapse = ", "))
    try(
      install.packages(
        missing,
        repos = "https://cloud.r-project.org",
        dependencies = TRUE
      ),
      silent = TRUE
    )
  }
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0L) {
    msg <- paste("Package(s) unavailable:", paste(still_missing, collapse = ", "))
    if (required) stop(msg) else add_warning("PACKAGE", paste(still_missing, collapse = ";"), msg)
  }
  invisible(setdiff(pkgs, still_missing))
}

ensure_cran(c("data.table", "openxlsx", "digest", "ggplot2", "zip"), required = TRUE)
if (INSPECT_H5AD) ensure_cran("hdf5r", required = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(digest)
  library(ggplot2)
})

log_msg("Stage 1 data audit started.")
log_msg("PROJECT_DIR: ", PROJECT_DIR)
log_msg("DATA_DIR: ", DATA_DIR)
log_msg("OUT_DIR: ", OUT_DIR)

############################################################
## 2. Utility functions
############################################################

normalize_slash <- function(x) gsub("\\\\", "/", x)

relative_to <- function(path, root) {
  path_n <- normalize_slash(normalizePath(path, winslash = "/", mustWork = FALSE))
  root_n <- normalize_slash(normalizePath(root, winslash = "/", mustWork = FALSE))
  prefix <- paste0(sub("/+$", "", root_n), "/")
  if (startsWith(path_n, prefix)) substring(path_n, nchar(prefix) + 1L) else path_n
}

extract_gse <- function(x) {
  m <- regexpr("GSE[0-9]+", x, ignore.case = TRUE)
  if (length(m) == 0L || is.na(m[1L]) || m[1L] < 0L) return(NA_character_)
  toupper(regmatches(x, m)[1L])
}

classify_dataset <- function(path) {
  b <- basename(path)
  gse <- extract_gse(b)
  if (!is.na(gse)) return(gse)
  if (grepl("HFpEF_snRNAseq.*\\.h5ad$", b, ignore.case = TRUE) ||
      grepl("^file_supplemental_info\\.tsv$", b, ignore.case = TRUE)) return("SCP3342")
  if (grepl("ligand_target_matrix.*mouse.*\\.rds$", b, ignore.case = TRUE) ||
      grepl("lr_network_mouse.*\\.rds$", b, ignore.case = TRUE)) return("NicheNet_mouse")
  if (grepl("HFpEF.*方案|执行总方案|reanalysis.*plan", b, ignore.case = TRUE)) return("Protocol_document")
  "Miscellaneous"
}

classify_file_type <- function(path) {
  b <- tolower(basename(path))
  if (grepl("series_matrix\\.txt\\.gz$", b)) return("GEO_series_matrix")
  if (grepl("family\\.soft\\.gz$", b)) return("GEO_family_SOFT")
  if (grepl("raw\\.tar$", b)) return("GEO_RAW_tar")
  if (grepl("cell_ranger_outs\\.tar\\.gz$", b)) return("CellRanger_tar_gz")
  if (grepl("matrix\\.mtx\\.gz$", b)) return("10x_matrix_mtx_gz")
  if (grepl("barcodes\\.tsv\\.gz$", b)) return("10x_barcodes_tsv_gz")
  if (grepl("features\\.tsv\\.gz$", b)) return("10x_features_tsv_gz")
  if (grepl("counts?\\.txt\\.gz$", b)) return("bulk_counts_txt_gz")
  if (grepl("tpm\\.txt\\.gz$", b)) return("bulk_TPM_txt_gz")
  if (grepl("length\\.txt\\.gz$", b)) return("gene_length_txt_gz")
  if (grepl("deg.*\\.txt\\.gz$", b)) return("precomputed_DEG_txt_gz")
  if (grepl("bulk.*rna.*seq.*\\.txt\\.gz$", b)) return("bulk_expression_txt_gz")
  if (grepl("integrated_counts.*\\.csv\\.gz$", b)) return("single_cell_counts_csv_gz")
  if (grepl("\\.h5ad$", b)) return("AnnData_h5ad")
  if (grepl("\\.rds$", b)) return("RDS_resource")
  if (grepl("\\.docx$", b)) return("Word_document")
  if (grepl("\\.xlsx$", b)) return("Excel_workbook")
  if (grepl("\\.csv\\.gz$", b)) return("CSV_gz")
  if (grepl("\\.txt\\.gz$", b)) return("TXT_gz")
  if (grepl("\\.tar\\.gz$|\\.tgz$", b)) return("tar_gz_archive")
  if (grepl("\\.tar$", b)) return("tar_archive")
  tools::file_ext(b)
}

safe_file_size_gb <- function(bytes) bytes / (1024^3)

open_text_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt", encoding = "UTF-8")
  } else {
    file(path, open = "rt", encoding = "UTF-8")
  }
}

strip_outer_quotes <- function(x) {
  x <- trimws(x)
  x <- sub('^"', '', x)
  x <- sub('"$', '', x)
  x
}

parse_geo_header_line <- function(line) {
  if (!grepl("=", line, fixed = TRUE)) {
    return(list(key = trimws(line), values = character()))
  }
  key <- trimws(sub("=.*$", "", line))
  rhs <- sub("^[^=]*=", "", line)
  vals <- strsplit(rhs, "\t", fixed = TRUE)[[1]]
  vals <- strip_outer_quotes(vals)
  list(key = key, values = vals)
}

read_geo_header <- function(path, max_lines = MAX_GEO_HEADER_LINES) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  kept <- character()
  total <- 0L
  truncated <- FALSE
  repeat {
    chunk <- readLines(con, n = 5000L, warn = FALSE)
    if (length(chunk) == 0L) break
    total <- total + length(chunk)
    kept <- c(kept, chunk[grepl("^!", chunk)])
    if (any(grepl("^!series_matrix_table_begin", chunk, ignore.case = TRUE))) break
    if (total >= max_lines) {
      truncated <- TRUE
      break
    }
  }
  list(lines = kept, lines_scanned = total, truncated = truncated)
}

parse_series_sample_metadata <- function(path, dataset_id) {
  hdr <- read_geo_header(path)
  parsed <- lapply(hdr$lines, parse_geo_header_line)
  keys <- vapply(parsed, `[[`, character(1), "key")
  acc_idx <- which(tolower(keys) == "!sample_geo_accession")
  if (length(acc_idx) == 0L) return(data.table())
  accessions <- parsed[[acc_idx[1L]]]$values
  n_samples <- length(accessions)
  if (n_samples == 0L) return(data.table())

  sample_idx <- which(grepl("^!Sample_", keys, ignore.case = TRUE))
  out <- list()
  field_counter <- list()
  for (i in sample_idx) {
    vals <- parsed[[i]]$values
    if (length(vals) != n_samples) next
    field <- sub("^!Sample_", "", keys[i], ignore.case = TRUE)
    field_counter[[field]] <- (field_counter[[field]] %||% 0L) + 1L
    field_instance <- field_counter[[field]]
    field_out <- if (field_instance > 1L) paste0(field, "_", field_instance) else field
    out[[length(out) + 1L]] <- data.table(
      dataset_id = dataset_id,
      source_file = basename(path),
      sample_accession = accessions,
      field = field_out,
      value = vals
    )
  }
  rbindlist(out, fill = TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

count_lines_stream <- function(path, chunk_n = 100000L) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  n <- 0
  repeat {
    x <- readLines(con, n = chunk_n, warn = FALSE)
    if (length(x) == 0L) break
    n <- n + length(x)
  }
  n
}

read_first_lines <- function(path, n = 5L) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

read_mtx_dimensions <- function(path) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  lines <- readLines(con, n = 200L, warn = FALSE)
  non_comment <- lines[!grepl("^%", lines)]
  non_comment <- non_comment[nzchar(trimws(non_comment))]
  if (length(non_comment) == 0L) return(c(NA_integer_, NA_integer_, NA_real_))
  vals <- suppressWarnings(as.numeric(strsplit(trimws(non_comment[1L]), "[[:space:]]+")[[1]]))
  if (length(vals) < 3L) return(c(NA_integer_, NA_integer_, NA_real_))
  vals[1:3]
}

inspect_delimited_file <- function(path, file_size_gb) {
  first <- tryCatch(read_first_lines(path, n = 3L), error = function(e) character())
  if (length(first) == 0L) {
    return(list(delimiter = NA_character_, n_columns_header = NA_integer_,
                n_lines = NA_real_, header_preview = NA_character_, status = "READ_FAILED"))
  }
  header <- first[1L]
  tab_n <- lengths(regmatches(header, gregexpr("\t", header, fixed = TRUE)))
  comma_n <- lengths(regmatches(header, gregexpr(",", header, fixed = TRUE)))
  delimiter <- if (tab_n >= comma_n) "TAB" else "COMMA"
  delim_char <- if (delimiter == "TAB") "\t" else ","
  n_cols <- length(strsplit(header, delim_char, fixed = TRUE)[[1]])
  n_lines <- NA_real_
  status <- "HEADER_ONLY"
  if (is.finite(file_size_gb) && file_size_gb <= COUNT_LINES_MAX_GB) {
    n_lines <- tryCatch(count_lines_stream(path), error = function(e) NA_real_)
    status <- ifelse(is.na(n_lines), "HEADER_ONLY", "HEADER_AND_LINE_COUNT")
  }
  list(
    delimiter = delimiter,
    n_columns_header = n_cols,
    n_lines = n_lines,
    header_preview = substr(header, 1L, 1000L),
    status = status
  )
}

safe_sha256 <- function(path, size_gb) {
  if (!COMPUTE_SHA256) return(NA_character_)
  if (!is.finite(size_gb) || size_gb > HASH_MAX_FILE_GB) return(NA_character_)
  tryCatch(
    digest::digest(file = path, algo = "sha256", serialize = FALSE),
    error = function(e) {
      add_warning("HASH", basename(path), conditionMessage(e))
      NA_character_
    }
  )
}

sanitize_sheet_name <- function(x) {
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
  substr(x, 1L, 31L)
}

write_sheet_safe <- function(wb, sheet, x) {
  sheet <- sanitize_sheet_name(sheet)
  addWorksheet(wb, sheet)
  if (is.null(x) || nrow(x) == 0L) {
    writeData(wb, sheet, data.frame(note = "No records generated."))
    return(invisible(NULL))
  }
  x_out <- as.data.frame(x, stringsAsFactors = FALSE)
  char_cols <- vapply(x_out, is.character, logical(1))
  if (any(char_cols)) {
    x_out[char_cols] <- lapply(x_out[char_cols], function(z) substr(z, 1L, 30000L))
  }
  max_excel_rows <- 1048575L
  if (nrow(x_out) > max_excel_rows) {
    writeData(wb, sheet, x_out[seq_len(max_excel_rows), , drop = FALSE])
    add_warning("EXCEL", sheet, paste0("Table truncated in XLSX at ", max_excel_rows, " rows; full CSV retained."))
  } else {
    writeData(wb, sheet, x_out)
  }
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = seq_len(min(ncol(x_out), 30L)), widths = "auto")
  invisible(NULL)
}

safe_rbindlist <- function(x) {
  if (length(x) == 0L) data.table() else rbindlist(x, fill = TRUE, use.names = TRUE)
}

############################################################
## 3. Expected dataset-role reference
############################################################

role_manifest <- data.table(
  dataset_id = c(
    "GSE237156", "GSE208425", "GSE209548", "GSE236585", "GSE236584",
    "GSE223527", "GSE183852", "GSE245034", "GSE249412", "GSE270896",
    "GSE275031", "SCP3342", "NicheNet_mouse"
  ),
  species = c(
    "Mus musculus", "Mus musculus", "Mus musculus", "Mus musculus", "Mus musculus",
    "Homo sapiens", "Homo sapiens", "Mus musculus", "Mus musculus", "Mus musculus",
    "Mus musculus", "Homo sapiens", "Mus musculus knowledge resource"
  ),
  tissue_or_source = c(
    "Sorted cardiac macrophage populations",
    "Cardiac immune single-cell context",
    "Cardiomyocyte bulk transcriptome",
    "Whole-heart/cardiac single-cell transcriptome",
    "Matched cardiac bulk transcriptome",
    "Peripheral blood mononuclear cells",
    "Human non-failing and dilated-cardiomyopathy heart atlas",
    "HFpEF heart bulk transcriptome with empagliflozin/TYA-018 arms",
    "HFpEF heart single-nucleus transcriptome with empagliflozin/TYA-018 arms",
    "Mouse heart single-nucleus model comparison",
    "Independent mouse HFpEF cardiac single-cell dataset",
    "Human HFpEF myocardium single-nucleus transcriptome",
    "Mouse ligand-receptor and ligand-target prior knowledge"
  ),
  modality = c(
    "Bulk RNA-seq", "scRNA-seq", "Bulk RNA-seq", "scRNA-seq", "Bulk RNA-seq",
    "10x GEX/CITE-seq-like public matrices", "Single-cell count matrix", "Bulk RNA-seq",
    "snRNA-seq", "snRNA-seq", "scRNA-seq", "snRNA-seq AnnData", "NicheNet RDS resources"
  ),
  planned_role = c(
    "Primary drug-opposed macrophage-program discovery",
    "Internal immune-cell contextualization; same project family as drug anchor",
    "Auxiliary cardiomyocyte response; not independent validation",
    "Primary cardiac single-cell discovery",
    "Matched orthogonal bulk support; not independent validation",
    "Auxiliary human circulating-immune support",
    "Human disease comparator; not HFpEF validation",
    "External sample-level SGLT2-inhibitor response validation",
    "Cell-type-resolved SGLT2-inhibitor validation; paired study with GSE245034",
    "HFpEF/model specificity and broader heart-failure comparison",
    "Independent external mouse HFpEF validation",
    "Primary independent human HFpEF myocardial validation",
    "Prior-knowledge resource for mouse NicheNet only"
  ),
  discovery_or_validation = c(
    "Discovery", "Internal context", "Internal context", "Discovery", "Internal support",
    "Human support", "Disease comparator", "External validation", "External validation",
    "Specificity validation", "External validation", "External human validation", "Resource"
  ),
  expected_sample_structure = c(
    "Diet/disease and dapagliflozin contrasts in Ccr2-positive and Ccr2-negative macrophages; verify exact sample metadata",
    "Cardiac immune cells; verify sample identifiers and relation to GSE237156",
    "Cardiomyocyte samples; verify treatment and diet groups",
    "HFpEF versus control biological samples; preserve sample identity for pseudobulk",
    "HFpEF versus control bulk biological samples",
    "Public control and HFpEF matrices; donor-level metadata may be limited",
    "Non-failing donor and DCM cells; use only as a non-HFpEF comparator",
    "Control, HFpEF vehicle, HFpEF empagliflozin, and HFpEF TYA-018 biological samples",
    "Control, HFpEF vehicle, HFpEF empagliflozin, and HFpEF TYA-018 nuclei",
    "Multiple mouse model groups; verify exact HFpEF/HFrEF phenotype assignment",
    "Independent HFpEF and control heart samples; small sample size expected",
    "Expected 43 donors (19 HFpEF and 24 controls) and 48,866 nuclei; verify adata.obs fields",
    "Not a biological cohort"
  ),
  claim_boundary = c(
    "Defines drug-opposed candidates but does not prove direct molecular causality",
    "Cannot be called independent replication of the drug-anchor project",
    "Cannot establish macrophage-to-cardiomyocyte causality",
    "Discovery dataset; no validation claim by itself",
    "Same-study support only",
    "Do not treat individual cells as independent patients; confirm donor IDs",
    "Cannot be described as human HFpEF validation",
    "Supports cross-SGLT2-inhibitor reversal at sample level",
    "Same study as GSE245034; provides orthogonal cell-type context, not a second independent cohort",
    "Model-specificity evidence depends on correct phenotype mapping",
    "Report effect sizes and direction consistency; avoid overreliance on P values",
    "Use donor-level pseudobulk; do not treat nuclei as independent replicates",
    "Resource-based prioritization is not physical binding or causal proof"
  )
)

############################################################
## 4. File inventory and SHA256 hashes
############################################################

all_files <- list.files(
  DATA_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)
all_files <- all_files[file.exists(all_files)]
if (length(all_files) == 0L) stop("No files found in DATA_DIR: ", DATA_DIR)

fi <- file.info(all_files)
inventory <- data.table(
  full_path = normalize_slash(all_files),
  relative_path = vapply(all_files, relative_to, character(1), root = DATA_DIR),
  filename = basename(all_files),
  dataset_id = vapply(all_files, classify_dataset, character(1)),
  file_type = vapply(all_files, classify_file_type, character(1)),
  size_bytes = as.numeric(fi$size),
  size_mb = as.numeric(fi$size) / (1024^2),
  size_gb = as.numeric(fi$size) / (1024^3),
  modified_time = format(fi$mtime, "%Y-%m-%d %H:%M:%S"),
  readable = file.access(all_files, 4L) == 0L
)
setorder(inventory, dataset_id, filename)

log_msg("Found ", nrow(inventory), " input files across ", uniqueN(inventory$dataset_id), " dataset/resource groups.")

inventory[, sha256 := NA_character_]
if (COMPUTE_SHA256) {
  for (i in seq_len(nrow(inventory))) {
    log_msg("SHA256 ", i, "/", nrow(inventory), ": ", inventory$filename[i])
    inventory$sha256[i] <- safe_sha256(inventory$full_path[i], inventory$size_gb[i])
  }
}

fwrite(inventory, file.path(DIRS$tables, "01_file_inventory_with_SHA256.csv"))
saveRDS(inventory, file.path(DIRS$objects, "stage1_file_inventory.rds"))

############################################################
## 5. Dataset manifest joined to local-file status
############################################################

local_summary <- inventory[, .(
  local_file_count = .N,
  total_size_gb = sum(size_gb, na.rm = TRUE),
  local_files = paste(filename, collapse = "; ")
), by = dataset_id]

dataset_manifest <- merge(role_manifest, local_summary, by = "dataset_id", all = TRUE)
dataset_manifest[is.na(local_file_count), local_file_count := 0L]
dataset_manifest[is.na(total_size_gb), total_size_gb := 0]
dataset_manifest[, local_status := fifelse(local_file_count > 0L, "PRESENT", "NOT_FOUND")]
role_levels <- c(
  "Discovery", "Internal context", "Internal support", "External validation",
  "Specificity validation", "External human validation", "Human support",
  "Disease comparator", "Resource"
)
dataset_manifest[, role_order := match(discovery_or_validation, role_levels)]
dataset_manifest[is.na(role_order), role_order := length(role_levels) + 1L]
setorder(dataset_manifest, role_order, dataset_id)
dataset_manifest[, role_order := NULL]

fwrite(dataset_manifest, file.path(DIRS$tables, "02_dataset_manifest.csv"))

############################################################
## 6. Required-input completeness checklist
############################################################

requirements <- data.table::rbindlist(list(
  data.table(dataset_id = "GSE237156", requirement = c("RAW archive", "Series matrix"),
             regex = c("^GSE237156_RAW\\.tar$", "^GSE237156_series_matrix\\.txt\\.gz$"), required = TRUE),
  data.table(dataset_id = "GSE208425", requirement = c("RAW archive", "Series matrix"),
             regex = c("^GSE208425_RAW\\.tar$", "^GSE208425_series_matrix\\.txt\\.gz$"), required = TRUE),
  data.table(dataset_id = "GSE209548", requirement = c("RAW archive", "Series matrix"),
             regex = c("^GSE209548_RAW\\.tar$", "^GSE209548_series_matrix\\.txt\\.gz$"), required = FALSE),
  data.table(dataset_id = "GSE236585", requirement = c("RAW archive", "Series matrix"),
             regex = c("^GSE236585_RAW\\.tar$", "^GSE236585_series_matrix\\.txt\\.gz$"), required = TRUE),
  data.table(dataset_id = "GSE236584", requirement = c("Bulk expression/count table", "Series matrix"),
             regex = c("^GSE236584_.*(bulk|count|RNA-seq).*\\.txt\\.gz$", "^GSE236584_series_matrix\\.txt\\.gz$"), required = TRUE),
  data.table(dataset_id = "GSE223527", requirement = c(
    "Control barcodes", "Control features", "Control matrix",
    "HFpEF barcodes", "HFpEF features", "HFpEF matrix", "Series matrix"),
    regex = c(
      "^GSE223527_ctrl_barcodes\\.tsv\\.gz$", "^GSE223527_ctrl_features\\.tsv\\.gz$", "^GSE223527_ctrl_matrix\\.mtx\\.gz$",
      "^GSE223527_hfpef_barcodes\\.tsv\\.gz$", "^GSE223527_hfpef_features\\.tsv\\.gz$", "^GSE223527_hfpef_matrix\\.mtx\\.gz$",
      "^GSE223527_series_matrix\\.txt\\.gz$"), required = TRUE),
  data.table(dataset_id = "GSE183852", requirement = c("Integrated count matrix", "Series matrix"),
             regex = c("^GSE183852_Integrated_Counts\\.csv\\.gz$", "^GSE183852_series_matrix\\.txt\\.gz$"), required = FALSE),
  data.table(dataset_id = "GSE245034", requirement = c("Raw count table", "TPM table", "Gene-length table", "Series matrix", "Family SOFT"),
             regex = c("^GSE245034_.*counts\\.txt\\.gz$", "^GSE245034_.*TPM\\.txt\\.gz$", "^GSE245034_.*length\\.txt\\.gz$",
                       "^GSE245034_series_matrix\\.txt\\.gz$", "^GSE245034_family\\.soft\\.gz$"),
             required = c(TRUE, FALSE, FALSE, TRUE, FALSE)),
  data.table(dataset_id = "GSE249412", requirement = c("RAW archive", "Series matrix", "Family SOFT"),
             regex = c("^GSE249412_RAW\\.tar$", "^GSE249412_series_matrix\\.txt\\.gz$", "^GSE249412_family\\.soft\\.gz$"),
             required = c(TRUE, TRUE, FALSE)),
  data.table(dataset_id = "GSE270896", requirement = c("RAW archive", "Series matrix", "Family SOFT"),
             regex = c("^GSE270896_RAW\\.tar$", "^GSE270896_series_matrix\\.txt\\.gz$", "^GSE270896_family\\.soft\\.gz$"),
             required = c(TRUE, TRUE, FALSE)),
  data.table(dataset_id = "GSE275031", requirement = c("Cell Ranger outputs archive", "Series matrix", "Family SOFT"),
             regex = c("^GSE275031_cell_ranger_outs\\.tar\\.gz$", "^GSE275031_series_matrix\\.txt\\.gz$", "^GSE275031_family\\.soft\\.gz$"),
             required = c(TRUE, TRUE, FALSE)),
  data.table(dataset_id = "SCP3342", requirement = c("Human HFpEF AnnData", "Supplemental file manifest"),
             regex = c("^HFpEF_snRNAseq_single_cell_portal_10\\.14\\.2025\\.h5ad$", "^file_supplemental_info\\.tsv$"),
             required = c(TRUE, FALSE)),
  data.table(dataset_id = "NicheNet_mouse", requirement = c("Mouse ligand-target matrix", "Mouse ligand-receptor network"),
             regex = c("^ligand_target_matrix_nsga2r_final_mouse\\.rds$", "^lr_network_mouse_21122021\\.rds$"), required = TRUE)
), fill = TRUE)

requirements[, `:=`(found_n = 0L, matched_files = NA_character_, status = NA_character_)]
for (i in seq_len(nrow(requirements))) {
  hits <- inventory[dataset_id == requirements$dataset_id[i] &
                      grepl(requirements$regex[i], filename, ignore.case = TRUE)]
  requirements$found_n[i] <- nrow(hits)
  requirements$matched_files[i] <- if (nrow(hits) > 0L) paste(hits$filename, collapse = "; ") else ""
  requirements$status[i] <- if (nrow(hits) > 0L) "PRESENT" else if (requirements$required[i]) "MISSING_REQUIRED" else "MISSING_OPTIONAL"
}

fwrite(requirements, file.path(DIRS$tables, "03_required_input_completeness.csv"))

############################################################
## 7. GEO metadata parsing
############################################################

geo_files <- inventory[file_type %in% c("GEO_series_matrix", "GEO_family_SOFT")]
geo_header_records <- list()
series_sample_records <- list()

if (nrow(geo_files) > 0L) {
  for (i in seq_len(nrow(geo_files))) {
    f <- geo_files$full_path[i]
    ds <- geo_files$dataset_id[i]
    log_msg("Parsing GEO metadata: ", basename(f))
    hdr <- tryCatch(read_geo_header(f), error = function(e) {
      add_warning("GEO_METADATA", basename(f), conditionMessage(e))
      NULL
    })
    if (is.null(hdr)) next

    if (length(hdr$lines) > 0L) {
      parsed <- lapply(hdr$lines, parse_geo_header_line)
      for (j in seq_along(parsed)) {
        vals <- parsed[[j]]$values
        geo_header_records[[length(geo_header_records) + 1L]] <- data.table(
          dataset_id = ds,
          source_file = basename(f),
          file_type = geo_files$file_type[i],
          line_index = j,
          metadata_key = parsed[[j]]$key,
          n_values = length(vals),
          values_joined = paste(vals, collapse = " | "),
          lines_scanned = hdr$lines_scanned,
          scan_truncated = hdr$truncated
        )
      }
    }

    if (geo_files$file_type[i] == "GEO_series_matrix") {
      sm <- tryCatch(parse_series_sample_metadata(f, ds), error = function(e) {
        add_warning("SERIES_SAMPLE_METADATA", basename(f), conditionMessage(e))
        data.table()
      })
      if (nrow(sm) > 0L) series_sample_records[[length(series_sample_records) + 1L]] <- sm
    }
  }
}

geo_headers <- safe_rbindlist(geo_header_records)
series_sample_metadata <- safe_rbindlist(series_sample_records)

fwrite(geo_headers, file.path(DIRS$tables, "04_GEO_header_metadata_long.csv"))
fwrite(series_sample_metadata, file.path(DIRS$tables, "05_GEO_series_sample_metadata_long.csv"))

############################################################
## 8. Archive-content audit
############################################################

archive_files <- inventory[file_type %in% c("GEO_RAW_tar", "CellRanger_tar_gz", "tar_gz_archive", "tar_archive")]
archive_records <- list()
archive_summary_records <- list()

if (nrow(archive_files) > 0L) {
  for (i in seq_len(nrow(archive_files))) {
    f <- archive_files$full_path[i]
    log_msg("Listing archive contents: ", basename(f))
    members <- tryCatch(
      utils::untar(f, list = TRUE),
      error = function(e) {
        add_warning("ARCHIVE", basename(f), conditionMessage(e))
        character()
      }
    )
    if (length(members) == 0L) {
      archive_summary_records[[length(archive_summary_records) + 1L]] <- data.table(
        dataset_id = archive_files$dataset_id[i],
        archive_file = basename(f),
        member_count = 0L,
        matrix_files = 0L,
        barcode_files = 0L,
        feature_files = 0L,
        h5_files = 0L,
        rds_files = 0L,
        listing_status = "EMPTY_OR_FAILED"
      )
      next
    }
    member_type <- fifelse(grepl("matrix\\.mtx(\\.gz)?$", members, ignore.case = TRUE), "matrix",
                    fifelse(grepl("barcodes\\.tsv(\\.gz)?$", members, ignore.case = TRUE), "barcodes",
                    fifelse(grepl("(features|genes)\\.tsv(\\.gz)?$", members, ignore.case = TRUE), "features",
                    fifelse(grepl("\\.h5$", members, ignore.case = TRUE), "h5",
                    fifelse(grepl("\\.rds$", members, ignore.case = TRUE), "rds", "other")))))
    archive_records[[length(archive_records) + 1L]] <- data.table(
      dataset_id = archive_files$dataset_id[i],
      archive_file = basename(f),
      member_index = seq_along(members),
      member_path = members,
      member_type = member_type
    )
    archive_summary_records[[length(archive_summary_records) + 1L]] <- data.table(
      dataset_id = archive_files$dataset_id[i],
      archive_file = basename(f),
      member_count = length(members),
      matrix_files = sum(member_type == "matrix"),
      barcode_files = sum(member_type == "barcodes"),
      feature_files = sum(member_type == "features"),
      h5_files = sum(member_type == "h5"),
      rds_files = sum(member_type == "rds"),
      listing_status = "OK"
    )
  }
}

archive_contents <- safe_rbindlist(archive_records)
archive_summary <- safe_rbindlist(archive_summary_records)
fwrite(archive_contents, file.path(DIRS$tables, "06_archive_contents.csv"))
fwrite(archive_summary, file.path(DIRS$tables, "07_archive_summary.csv"))

############################################################
## 9. Matrix and text-file structure audit
############################################################

matrix_structure <- list()

## 9a. Standalone 10x Matrix Market files
mtx_files <- inventory[file_type == "10x_matrix_mtx_gz"]
if (nrow(mtx_files) > 0L) {
  for (i in seq_len(nrow(mtx_files))) {
    dims <- tryCatch(read_mtx_dimensions(mtx_files$full_path[i]), error = function(e) {
      add_warning("MTX", mtx_files$filename[i], conditionMessage(e))
      c(NA_real_, NA_real_, NA_real_)
    })
    matrix_structure[[length(matrix_structure) + 1L]] <- data.table(
      dataset_id = mtx_files$dataset_id[i],
      source_file = mtx_files$filename[i],
      structure_type = "MatrixMarket",
      n_features_or_rows = dims[1],
      n_cells_or_columns = dims[2],
      nonzero_entries = dims[3],
      n_lines = NA_real_,
      delimiter = NA_character_,
      header_preview = NA_character_,
      inspection_status = ifelse(all(!is.na(dims)), "OK", "FAILED")
    )
  }
}

## 9b. Standalone barcode and feature files
line_count_files <- inventory[file_type %in% c("10x_barcodes_tsv_gz", "10x_features_tsv_gz")]
if (nrow(line_count_files) > 0L) {
  for (i in seq_len(nrow(line_count_files))) {
    n_lines <- tryCatch(count_lines_stream(line_count_files$full_path[i]), error = function(e) {
      add_warning("LINE_COUNT", line_count_files$filename[i], conditionMessage(e))
      NA_real_
    })
    matrix_structure[[length(matrix_structure) + 1L]] <- data.table(
      dataset_id = line_count_files$dataset_id[i],
      source_file = line_count_files$filename[i],
      structure_type = line_count_files$file_type[i],
      n_features_or_rows = ifelse(line_count_files$file_type[i] == "10x_features_tsv_gz", n_lines, NA_real_),
      n_cells_or_columns = ifelse(line_count_files$file_type[i] == "10x_barcodes_tsv_gz", n_lines, NA_real_),
      nonzero_entries = NA_real_,
      n_lines = n_lines,
      delimiter = "TAB",
      header_preview = NA_character_,
      inspection_status = ifelse(is.na(n_lines), "FAILED", "OK")
    )
  }
}

## 9c. Generic compressed tables
text_types <- c(
  "bulk_counts_txt_gz", "bulk_TPM_txt_gz", "gene_length_txt_gz",
  "precomputed_DEG_txt_gz", "bulk_expression_txt_gz",
  "single_cell_counts_csv_gz", "CSV_gz", "TXT_gz"
)
text_files <- inventory[file_type %in% text_types]
if (nrow(text_files) > 0L) {
  for (i in seq_len(nrow(text_files))) {
    log_msg("Inspecting delimited file header: ", text_files$filename[i])
    ins <- tryCatch(
      inspect_delimited_file(text_files$full_path[i], text_files$size_gb[i]),
      error = function(e) {
        add_warning("DELIMITED", text_files$filename[i], conditionMessage(e))
        list(delimiter = NA_character_, n_columns_header = NA_integer_, n_lines = NA_real_,
             header_preview = NA_character_, status = "FAILED")
      }
    )
    matrix_structure[[length(matrix_structure) + 1L]] <- data.table(
      dataset_id = text_files$dataset_id[i],
      source_file = text_files$filename[i],
      structure_type = text_files$file_type[i],
      n_features_or_rows = ifelse(is.na(ins$n_lines), NA_real_, max(ins$n_lines - 1, 0)),
      n_cells_or_columns = ifelse(is.na(ins$n_columns_header), NA_real_, max(ins$n_columns_header - 1, 0)),
      nonzero_entries = NA_real_,
      n_lines = ins$n_lines,
      delimiter = ins$delimiter,
      header_preview = ins$header_preview,
      inspection_status = ins$status
    )
  }
}

matrix_structure_dt <- safe_rbindlist(matrix_structure)
fwrite(matrix_structure_dt, file.path(DIRS$tables, "08_matrix_and_table_structure.csv"))

############################################################
## 10. H5AD structural audit without loading the expression matrix
############################################################

h5ad_structure <- data.table()
h5ad_summary <- data.table()

inspect_h5ad_file <- function(path, dataset_id) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  ls_all <- as.data.table(h5$ls(recursive = TRUE))
  if (nrow(ls_all) > 0L) {
    ls_all[, `:=`(
      dataset_id = dataset_id,
      source_file = basename(path)
    )]
  }

  get_index_length <- function(group_name) {
    grp <- tryCatch(h5[[group_name]], error = function(e) NULL)
    if (is.null(grp)) return(NA_real_)

    index_name <- tryCatch(hdf5r::h5attr(grp, "_index"), error = function(e) NA_character_)
    candidates <- unique(c(as.character(index_name), "_index"))
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

    for (nm in candidates) {
      idx <- tryCatch(grp[[nm]], error = function(e) NULL)
      if (!is.null(idx)) {
        d <- tryCatch(idx$dims, error = function(e) NA)
        if (length(d) > 0L && !is.na(d[1L])) return(as.numeric(d[1L]))
      }
    }

    grp_ls <- tryCatch(grp$ls(recursive = FALSE), error = function(e) NULL)
    if (!is.null(grp_ls) && nrow(grp_ls) > 0L) {
      for (nm in grp_ls$name) {
        obj <- tryCatch(grp[[nm]], error = function(e) NULL)
        d <- if (!is.null(obj)) tryCatch(obj$dims, error = function(e) NA) else NA
        if (length(d) > 0L && !is.na(d[1L])) return(as.numeric(d[1L]))
      }
    }
    NA_real_
  }

  obs_n <- get_index_length("obs")
  var_n <- get_index_length("var")
  obs_grp <- tryCatch(h5[["obs"]], error = function(e) NULL)
  obs_fields <- if (!is.null(obs_grp)) {
    tryCatch(obs_grp$ls(recursive = FALSE)$name, error = function(e) character())
  } else character()
  var_grp <- tryCatch(h5[["var"]], error = function(e) NULL)
  var_fields <- if (!is.null(var_grp)) {
    tryCatch(var_grp$ls(recursive = FALSE)$name, error = function(e) character())
  } else character()

  top_level <- tryCatch(h5$ls(recursive = FALSE)$name, error = function(e) character())
  summary <- data.table(
    dataset_id = dataset_id,
    source_file = basename(path),
    obs_count = obs_n,
    var_count = var_n,
    obs_field_count = length(obs_fields),
    obs_fields = paste(obs_fields, collapse = "; "),
    var_field_count = length(var_fields),
    var_fields_preview = paste(head(var_fields, 100L), collapse = "; "),
    top_level_groups = paste(top_level, collapse = "; "),
    has_X = "X" %in% top_level,
    has_raw = "raw" %in% top_level,
    has_layers = "layers" %in% top_level,
    has_obsm = "obsm" %in% top_level,
    inspection_status = "OK"
  )
  list(structure = ls_all, summary = summary)
}

h5ad_files <- inventory[file_type == "AnnData_h5ad"]
if (INSPECT_H5AD && nrow(h5ad_files) > 0L) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    add_warning("H5AD", "hdf5r", "hdf5r is unavailable; H5AD internal structure was not inspected.")
  } else {
    for (i in seq_len(nrow(h5ad_files))) {
      f <- h5ad_files$full_path[i]
      log_msg("Inspecting H5AD structure: ", basename(f))
      result <- tryCatch(
        inspect_h5ad_file(f, h5ad_files$dataset_id[i]),
        error = function(e) {
          add_warning("H5AD", basename(f), conditionMessage(e))
          list(
            structure = data.table(),
            summary = data.table(
              dataset_id = h5ad_files$dataset_id[i],
              source_file = basename(f),
              obs_count = NA_real_, var_count = NA_real_,
              obs_field_count = NA_integer_, obs_fields = NA_character_,
              var_field_count = NA_integer_, var_fields_preview = NA_character_,
              top_level_groups = NA_character_,
              has_X = NA, has_raw = NA, has_layers = NA, has_obsm = NA,
              inspection_status = "FAILED"
            )
          )
        }
      )
      h5ad_structure <- rbind(h5ad_structure, result$structure, fill = TRUE)
      h5ad_summary <- rbind(h5ad_summary, result$summary, fill = TRUE)
    }
  }
}

fwrite(h5ad_structure, file.path(DIRS$tables, "09_H5AD_structure.csv"))
fwrite(h5ad_summary, file.path(DIRS$tables, "10_H5AD_summary.csv"))

############################################################
## 11. NicheNet RDS resource audit
############################################################

rds_structure <- data.table()
rds_files <- inventory[file_type == "RDS_resource"]
if (INSPECT_RDS && nrow(rds_files) > 0L) {
  for (i in seq_len(nrow(rds_files))) {
    f <- rds_files$full_path[i]
    if (rds_files$size_gb[i] > DEEP_INSPECT_MAX_GB) {
      add_warning("RDS", basename(f), "Skipped deep inspection because the file exceeds DEEP_INSPECT_MAX_GB.")
      next
    }
    log_msg("Inspecting RDS resource: ", basename(f))
    rec <- tryCatch({
      obj <- readRDS(f)
      d <- dim(obj)
      nms <- names(obj)
      data.table(
        dataset_id = rds_files$dataset_id[i],
        source_file = basename(f),
        object_class = paste(class(obj), collapse = "; "),
        dim_1 = ifelse(length(d) >= 1L, d[1L], NA_real_),
        dim_2 = ifelse(length(d) >= 2L, d[2L], NA_real_),
        object_length = length(obj),
        names_preview = paste(head(nms, 50L), collapse = "; "),
        memory_mb_after_load = as.numeric(object.size(obj)) / (1024^2),
        inspection_status = "OK"
      )
    }, error = function(e) {
      add_warning("RDS", basename(f), conditionMessage(e))
      data.table(
        dataset_id = rds_files$dataset_id[i],
        source_file = basename(f),
        object_class = NA_character_, dim_1 = NA_real_, dim_2 = NA_real_,
        object_length = NA_real_, names_preview = NA_character_,
        memory_mb_after_load = NA_real_, inspection_status = "FAILED"
      )
    })
    rds_structure <- rbind(rds_structure, rec, fill = TRUE)
    rm(rec)
    gc()
  }
}

fwrite(rds_structure, file.path(DIRS$tables, "11_RDS_resource_structure.csv"))

############################################################
## 12. Cross-file consistency checks
############################################################

consistency_checks <- list()

## Match standalone 10x matrix dimensions to barcode/feature line counts.
for (prefix in c("GSE223527_ctrl", "GSE223527_hfpef")) {
  mtx_row <- matrix_structure_dt[grepl(paste0("^", prefix, "_matrix\\.mtx\\.gz$"), source_file, ignore.case = TRUE)]
  bc_row  <- matrix_structure_dt[grepl(paste0("^", prefix, "_barcodes\\.tsv\\.gz$"), source_file, ignore.case = TRUE)]
  ft_row  <- matrix_structure_dt[grepl(paste0("^", prefix, "_features\\.tsv\\.gz$"), source_file, ignore.case = TRUE)]

  mtx_cells <- if (nrow(mtx_row) > 0L) mtx_row$n_cells_or_columns[1L] else NA_real_
  mtx_genes <- if (nrow(mtx_row) > 0L) mtx_row$n_features_or_rows[1L] else NA_real_
  bc_n <- if (nrow(bc_row) > 0L) bc_row$n_lines[1L] else NA_real_
  ft_n <- if (nrow(ft_row) > 0L) ft_row$n_lines[1L] else NA_real_

  consistency_checks[[length(consistency_checks) + 1L]] <- data.table(
    check_id = paste0(prefix, "_10x_dimension_match"),
    dataset_id = "GSE223527",
    expected = paste0("matrix rows=features lines; matrix columns=barcode lines"),
    observed = paste0("matrix=", mtx_genes, "x", mtx_cells, "; features=", ft_n, "; barcodes=", bc_n),
    status = ifelse(!is.na(mtx_cells) && !is.na(mtx_genes) && !is.na(bc_n) && !is.na(ft_n) &&
                      mtx_cells == bc_n && mtx_genes == ft_n, "PASS", "CHECK")
  )
}

## H5AD expected cell and gene dimensions from the portal description.
if (nrow(h5ad_summary) > 0L) {
  for (i in seq_len(nrow(h5ad_summary))) {
    consistency_checks[[length(consistency_checks) + 1L]] <- data.table(
      check_id = "SCP3342_expected_dimensions",
      dataset_id = h5ad_summary$dataset_id[i],
      expected = "approximately 48,866 nuclei and 36,601 genes",
      observed = paste0("obs=", h5ad_summary$obs_count[i], "; var=", h5ad_summary$var_count[i]),
      status = ifelse(
        isTRUE(h5ad_summary$obs_count[i] == 48866) && isTRUE(h5ad_summary$var_count[i] == 36601),
        "PASS", "CHECK"
      )
    )
  }
}

consistency_dt <- safe_rbindlist(consistency_checks)
fwrite(consistency_dt, file.path(DIRS$tables, "12_cross_file_consistency_checks.csv"))

############################################################
## 13. Summary figures
############################################################

file_size_summary <- inventory[, .(
  file_count = .N,
  total_size_gb = sum(size_gb, na.rm = TRUE)
), by = dataset_id]
setorder(file_size_summary, -total_size_gb)
fwrite(file_size_summary, file.path(DIRS$source, "FigS1_file_size_by_dataset_source.csv"))

p1 <- ggplot(file_size_summary, aes(x = reorder(dataset_id, total_size_gb), y = total_size_gb)) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Total local file size (GB)",
    title = "Stage 1 input-data footprint by dataset/resource"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(DIRS$figures, "FigS1_input_file_size_by_dataset.png"), p1, width = 10, height = 7, dpi = 300)
ggsave(file.path(DIRS$figures, "FigS1_input_file_size_by_dataset.pdf"), p1, width = 10, height = 7)

file_type_summary <- inventory[, .N, by = .(dataset_id, file_type)]
fwrite(file_type_summary, file.path(DIRS$source, "FigS2_file_type_counts_source.csv"))

p2 <- ggplot(file_type_summary, aes(x = reorder(dataset_id, N), y = N, fill = file_type)) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of local files",
    fill = "File type",
    title = "Stage 1 file-type composition"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right")

ggsave(file.path(DIRS$figures, "FigS2_file_type_composition.png"), p2, width = 12, height = 8, dpi = 300)
ggsave(file.path(DIRS$figures, "FigS2_file_type_composition.pdf"), p2, width = 12, height = 8)

############################################################
## 14. Warnings, run status, methods, and README
############################################################

warnings_dt <- if (length(warning_records) > 0L) rbindlist(warning_records, fill = TRUE) else data.table(
  timestamp = character(), category = character(), item = character(), message = character()
)
fwrite(warnings_dt, file.path(DIRS$tables, "13_warnings_and_nonfatal_issues.csv"))

required_missing <- requirements[required == TRUE & status == "MISSING_REQUIRED"]
nonreadable <- inventory[readable == FALSE]
critical_failures <- nrow(required_missing) + nrow(nonreadable)

run_status <- data.table(
  stage = STAGE_NAME,
  start_time = if (file.exists(LOG_FILE)) format(file.info(LOG_FILE)$ctime, "%Y-%m-%d %H:%M:%S") else NA_character_,
  finish_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input_directory = DATA_DIR,
  output_directory = OUT_DIR,
  input_file_count = nrow(inventory),
  dataset_group_count = uniqueN(inventory$dataset_id),
  required_item_count = requirements[required == TRUE, .N],
  missing_required_count = nrow(required_missing),
  nonreadable_file_count = nrow(nonreadable),
  warning_count = nrow(warnings_dt),
  overall_status = ifelse(
    critical_failures > 0L,
    "FAIL_REQUIRED_INPUT_CHECK",
    ifelse(nrow(warnings_dt) > 0L, "PASS_WITH_WARNINGS", "PASS")
  )
)
fwrite(run_status, file.path(DIRS$tables, "00_stage1_run_status.csv"))

package_versions <- data.table(
  package = c("R", "data.table", "openxlsx", "digest", "ggplot2", "zip", "hdf5r"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(utils::packageVersion("data.table")),
    as.character(utils::packageVersion("openxlsx")),
    as.character(utils::packageVersion("digest")),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("zip")),
    if (requireNamespace("hdf5r", quietly = TRUE)) as.character(utils::packageVersion("hdf5r")) else NA_character_
  )
)
fwrite(package_versions, file.path(DIRS$tables, "14_software_versions.csv"))

methods_text <- c(
  "HFpEF Reanalysis Project — Stage 1 Data Audit",
  "",
  paste0("Input directory: ", DATA_DIR),
  paste0("Output directory: ", OUT_DIR),
  "",
  "Stage objective:",
  "Create a reproducible inventory and metadata audit before any biological analysis.",
  "",
  "Core operations:",
  "1. Recursive local-file inventory with file type, size, modification time, readability, and SHA256 hash.",
  "2. Dataset-role manifest separating discovery, internal support, external validation, human validation, disease comparison, and knowledge resources.",
  "3. Required-input completeness checks using explicit filename patterns.",
  "4. GEO series/family metadata parsing without reading expression matrices into memory.",
  "5. TAR archive listing to verify the presence of matrix, barcode, feature, H5, and RDS members.",
  "6. Structural checks for standalone 10x Matrix Market files and compressed tabular matrices.",
  "7. H5AD structural inspection through HDF5 metadata only; the full expression matrix is not loaded.",
  "8. NicheNet RDS class and dimension inspection.",
  "9. Cross-file consistency checks for 10x dimensions and expected SCP3342 object dimensions.",
  "",
  "Analysis boundary:",
  "No normalization, clustering, cell annotation, differential expression, enrichment, communication inference, or candidate prioritization is performed in this stage.",
  "",
  paste0("SHA256 enabled: ", COMPUTE_SHA256),
  paste0("Maximum file size for hashing: ", HASH_MAX_FILE_GB, " GB"),
  paste0("Maximum file size for full line counting: ", COUNT_LINES_MAX_GB, " GB"),
  paste0("Maximum GEO metadata lines scanned per file: ", MAX_GEO_HEADER_LINES)
)
writeLines(methods_text, file.path(DIRS$methods, "Stage1_methods_and_boundaries.txt"), useBytes = TRUE)

readme_text <- c(
  "HFpEF Reanalysis Project",
  "Stage 1: Complete input-data audit",
  "",
  paste0("Status: ", run_status$overall_status),
  paste0("Input files: ", run_status$input_file_count),
  paste0("Dataset/resource groups: ", run_status$dataset_group_count),
  paste0("Missing required items: ", run_status$missing_required_count),
  paste0("Warnings/nonfatal issues: ", run_status$warning_count),
  "",
  "Main outputs:",
  "- 01_tables/01_file_inventory_with_SHA256.csv",
  "- 01_tables/02_dataset_manifest.csv",
  "- 01_tables/03_required_input_completeness.csv",
  "- 01_tables/04_GEO_header_metadata_long.csv",
  "- 01_tables/05_GEO_series_sample_metadata_long.csv",
  "- 01_tables/06_archive_contents.csv",
  "- 01_tables/08_matrix_and_table_structure.csv",
  "- 01_tables/10_H5AD_summary.csv",
  "- 01_tables/11_RDS_resource_structure.csv",
  "- 01_tables/12_cross_file_consistency_checks.csv",
  "- 01_tables/Stage1_data_audit_report.xlsx",
  "",
  "Interpretation:",
  "This stage only determines whether the downloaded inputs are present, readable, structurally interpretable, and correctly assigned to discovery or validation roles.",
  "Biological analysis begins only after this CHECK package has been reviewed."
)
writeLines(readme_text, file.path(OUT_DIR, "README_stage1.txt"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(DIRS$logs, "sessionInfo.txt"))

############################################################
## 15. Consolidated Excel audit workbook
############################################################

wb <- createWorkbook(creator = "HFpEF Stage 1 data audit")
write_sheet_safe(wb, "RunStatus", run_status)
write_sheet_safe(wb, "FileInventory", inventory)
write_sheet_safe(wb, "DatasetManifest", dataset_manifest)
write_sheet_safe(wb, "Completeness", requirements)
write_sheet_safe(wb, "GEOHeaders", geo_headers)
write_sheet_safe(wb, "SampleMetadata", series_sample_metadata)
write_sheet_safe(wb, "ArchiveSummary", archive_summary)
write_sheet_safe(wb, "ArchiveContents", archive_contents)
write_sheet_safe(wb, "MatrixStructure", matrix_structure_dt)
write_sheet_safe(wb, "H5ADSummary", h5ad_summary)
write_sheet_safe(wb, "H5ADStructure", h5ad_structure)
write_sheet_safe(wb, "RDSStructure", rds_structure)
write_sheet_safe(wb, "Consistency", consistency_dt)
write_sheet_safe(wb, "Warnings", warnings_dt)
write_sheet_safe(wb, "SoftwareVersions", package_versions)

saveWorkbook(
  wb,
  file.path(DIRS$tables, "Stage1_data_audit_report.xlsx"),
  overwrite = TRUE
)

############################################################
## 16. Build compact review folder and CHECK.zip
############################################################

copy_into_check <- function(path, subdir = "") {
  if (!file.exists(path)) return(FALSE)
  target_dir <- file.path(DIRS$check, subdir)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(path, file.path(target_dir, basename(path)), overwrite = TRUE)
}

## Copy the current script when it can be resolved.
cmd_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(script_arg) > 0L) sub("^--file=", "", script_arg[1L]) else file.path(PROJECT_DIR, "HFpEF_Stage1_Data_Audit_FIXED_v2.R")
if (file.exists(script_path)) {
  copy_into_check(script_path, "code")
} else {
  writeLines(
    c(
      "The running script could not be automatically copied.",
      paste0("Expected path: ", script_path),
      "Please place HFpEF_Stage1_Data_Audit_FIXED_v2.R in the CHECK package manually if needed."
    ),
    file.path(DIRS$check, "SCRIPT_COPY_NOTICE.txt")
  )
}

copy_into_check(file.path(OUT_DIR, "README_stage1.txt"), "")
copy_into_check(LOG_FILE, "logs")
copy_into_check(WARN_FILE, "logs")
copy_into_check(file.path(DIRS$logs, "sessionInfo.txt"), "logs")
copy_into_check(file.path(DIRS$methods, "Stage1_methods_and_boundaries.txt"), "methods")

## All audit tables are compact enough for review; large raw data are never copied.
for (f in list.files(DIRS$tables, full.names = TRUE)) copy_into_check(f, "tables")
for (f in list.files(DIRS$source, full.names = TRUE)) copy_into_check(f, "source_data")
for (f in list.files(DIRS$figures, full.names = TRUE, pattern = "\\.(png|pdf)$", ignore.case = TRUE)) copy_into_check(f, "figures")

old_wd <- getwd()
tryCatch({
  setwd(DIRS$check)
  check_files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
  zip::zipr(CHECK_ZIP, files = check_files)
}, finally = {
  setwd(old_wd)
})

check_zip_size_mb <- file.info(CHECK_ZIP)$size / (1024^2)
CHECK_MAX_MB <- 45
if (is.finite(check_zip_size_mb) && check_zip_size_mb > CHECK_MAX_MB) {
  add_warning("CHECK_ZIP", basename(CHECK_ZIP), paste0(
    "Initial CHECK package was ", sprintf("%.2f", check_zip_size_mb),
    " MB. PDF previews and the largest long-form CSV files were removed to keep the package uploadable."
  ))
  unlink(list.files(file.path(DIRS$check, "figures"), pattern = "\\.pdf$", full.names = TRUE), force = TRUE)
  optional_large <- c(
    file.path(DIRS$check, "tables", "04_GEO_header_metadata_long.csv"),
    file.path(DIRS$check, "tables", "06_archive_contents.csv"),
    file.path(DIRS$check, "tables", "09_H5AD_structure.csv")
  )
  unlink(optional_large[file.exists(optional_large)], force = TRUE)
  unlink(CHECK_ZIP, force = TRUE)
  tryCatch({
    setwd(DIRS$check)
    check_files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
    zip::zipr(CHECK_ZIP, files = check_files)
  }, finally = {
    setwd(old_wd)
  })
  check_zip_size_mb <- file.info(CHECK_ZIP)$size / (1024^2)
}
log_msg("CHECK package created: ", CHECK_ZIP)
log_msg("CHECK package size: ", sprintf("%.2f MB", check_zip_size_mb))

############################################################
## 17. Final status and controlled stop if required inputs fail
############################################################

log_msg("Stage 1 completed with status: ", run_status$overall_status)
log_msg("Main report: ", file.path(DIRS$tables, "Stage1_data_audit_report.xlsx"))
log_msg("Review package: ", CHECK_ZIP)

cat("\n============================================================\n")
cat("HFpEF Stage 1 data audit completed\n")
cat("Status: ", run_status$overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat("CHECK size: ", sprintf("%.2f MB", check_zip_size_mb), "\n", sep = "")
cat("============================================================\n\n")

if (critical_failures > 0L) {
  stop(
    "Stage 1 finished its audit outputs, but one or more required inputs are missing or unreadable. ",
    "Review 03_required_input_completeness.csv and upload the CHECK.zip before proceeding."
  )
}

invisible(list(
  status = run_status,
  output_directory = OUT_DIR,
  check_zip = CHECK_ZIP,
  inventory = inventory,
  dataset_manifest = dataset_manifest,
  completeness = requirements
))
