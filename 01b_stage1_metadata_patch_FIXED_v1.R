############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 1 Metadata Patch FIXED v1
##
## Project root:
##   <HFPEF_PROJECT_DIR>
##
## Read-only input directory:
##   <HFPEF_PROJECT_DIR>/0.GEO
##
## Required previous output:
##   <HFPEF_PROJECT_DIR>/01_stage1_data_audit_FIXED_v2
##
## Output directory:
##   <HFPEF_PROJECT_DIR>/01_stage1_metadata_patch_FIXED_v1
##
## Review package:
##   <HFPEF_PROJECT_DIR>/01_stage1_metadata_patch_FIXED_v1_CHECK.zip
##
## Purpose:
##   1) Correctly parse GEO sample-level metadata from series-matrix and SOFT files.
##   2) Produce one-row-per-sample metadata tables.
##   3) Expand repeated characteristics fields into structured key-value tables.
##   4) Generate conservative, reviewable grouping candidates without forcing labels.
##   5) Inspect SCP3342 donor-level metadata from AnnData obs without loading X.
##   6) Archive the exact executed R script in the output and CHECK package.
##
## Scientific boundary:
##   - This patch only repairs and locks metadata.
##   - It does not perform normalization, differential expression,
##     clustering, annotation, pathway analysis, or candidate ranking.
##   - Inferred grouping fields are candidates for manual confirmation.
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
PREV_STAGE  <- file.path(PROJECT_DIR, "01_stage1_data_audit_FIXED_v2")

STAGE_NAME  <- "01_stage1_metadata_patch_FIXED_v1"
OUT_DIR     <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP   <- file.path(PROJECT_DIR, paste0(STAGE_NAME, "_CHECK.zip"))

ALLOW_OVERWRITE <- FALSE
INSPECT_SCP3342 <- TRUE
MAX_HEADER_LINES <- 300000L

EXPECTED_SAMPLE_COUNTS <- data.frame(
  dataset_id = c(
    "GSE223527", "GSE236584", "GSE236585", "GSE237156",
    "GSE245034", "GSE249412", "GSE270896", "GSE275031"
  ),
  expected_n = c(4L, 12L, 6L, 16L, 38L, 8L, 8L, 4L),
  expectation_basis = c(
    "Four public group-level GEX/ADT sample entries; donor-level structure may be unavailable",
    "Matched HFpEF versus control bulk cohort expected from study design",
    "Three HFpEF and three control cardiac single-cell samples",
    "Diet/disease x dapagliflozin x Ccr2-status macrophage design",
    "Control, HFpEF vehicle, HFpEF empagliflozin, and HFpEF TYA-018 sample-level bulk design",
    "Two samples per four treatment/model groups",
    "Eight mouse-heart snRNA-seq samples across model groups",
    "Two HFpEF and two control cardiac single-cell samples"
  ),
  stringsAsFactors = FALSE
)

############################################################
## 1. Output setup and packages
############################################################

if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
if (!dir.exists(DATA_DIR)) stop("DATA_DIR does not exist: ", DATA_DIR)
if (!dir.exists(PREV_STAGE)) warning("Previous Stage 1 folder was not found: ", PREV_STAGE)

if (dir.exists(OUT_DIR) && !ALLOW_OVERWRITE) {
  stop(
    "Output folder already exists and overwrite protection is active:\n",
    OUT_DIR,
    "\nUse a new version name or back up the old output before rerunning."
  )
}
if (dir.exists(OUT_DIR) && ALLOW_OVERWRITE) unlink(OUT_DIR, recursive = TRUE, force = TRUE)
if (file.exists(CHECK_ZIP) && !ALLOW_OVERWRITE) {
  stop("CHECK zip already exists and overwrite protection is active:\n", CHECK_ZIP)
}
if (file.exists(CHECK_ZIP) && ALLOW_OVERWRITE) unlink(CHECK_ZIP, force = TRUE)

DIRS <- list(
  logs    = file.path(OUT_DIR, "00_logs"),
  tables  = file.path(OUT_DIR, "01_tables"),
  objects = file.path(OUT_DIR, "02_objects"),
  figures = file.path(OUT_DIR, "03_figures"),
  source  = file.path(OUT_DIR, "04_source_data"),
  methods = file.path(OUT_DIR, "05_methods"),
  check   = file.path(OUT_DIR, "06_review_check")
)
for (d in c(OUT_DIR, unlist(DIRS, use.names = FALSE))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE  <- file.path(DIRS$logs, "stage1_metadata_patch.log")
WARN_FILE <- file.path(DIRS$logs, "stage1_metadata_patch_warnings.log")
START_TIME <- Sys.time()

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
      install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE),
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

ensure_cran(c("data.table", "openxlsx", "zip"), required = TRUE)
if (INSPECT_SCP3342) ensure_cran("hdf5r", required = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

log_msg("Stage 1 metadata patch started.")
log_msg("PROJECT_DIR: ", PROJECT_DIR)
log_msg("DATA_DIR: ", DATA_DIR)
log_msg("OUT_DIR: ", OUT_DIR)

############################################################
## 2. General utility functions
############################################################

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

normalize_slash <- function(x) gsub("\\\\", "/", x)

extract_gse <- function(x) {
  m <- regexpr("GSE[0-9]+", x, ignore.case = TRUE)
  if (length(m) == 0L || is.na(m[1L]) || m[1L] < 0L) return(NA_character_)
  toupper(regmatches(x, m)[1L])
}

open_text_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt", encoding = "UTF-8")
  } else {
    file(path, open = "rt", encoding = "UTF-8")
  }
}

strip_quotes <- function(x) {
  x <- trimws(as.character(x))
  x <- sub('^"', '', x)
  x <- sub('"$', '', x)
  x
}

clean_field_name <- function(x) {
  x <- sub("^!Sample_", "", x, ignore.case = TRUE)
  x <- sub("^!Series_", "", x, ignore.case = TRUE)
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

sanitize_value <- function(x) {
  x <- strip_quotes(x)
  x <- gsub("\\r|\\n", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

get_script_path <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    of <- tryCatch(frames[[i]]$ofile, error = function(e) NULL)
    if (!is.null(of) && length(of) == 1L && nzchar(of)) {
      return(normalizePath(of, winslash = "/", mustWork = FALSE))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

safe_rbind <- function(x) {
  if (length(x) == 0L) data.table() else rbindlist(x, fill = TRUE, use.names = TRUE)
}

write_csv_safe <- function(x, path) {
  if (is.null(x) || ncol(x) == 0L) {
    fwrite(data.table(note = "No records generated."), path)
  } else {
    fwrite(x, path)
  }
}

sanitize_sheet_name <- function(x) {
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
  substr(x, 1L, 31L)
}

write_sheet_safe <- function(wb, sheet, x) {
  sheet <- sanitize_sheet_name(sheet)
  addWorksheet(wb, sheet)
  if (is.null(x) || nrow(x) == 0L || ncol(x) == 0L) {
    writeData(wb, sheet, data.frame(note = "No records generated."))
    return(invisible(NULL))
  }
  y <- as.data.frame(x, stringsAsFactors = FALSE)
  char_cols <- vapply(y, is.character, logical(1))
  if (any(char_cols)) y[char_cols] <- lapply(y[char_cols], function(z) substr(z, 1L, 30000L))
  writeData(wb, sheet, y)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = seq_len(min(ncol(y), 30L)), widths = "auto")
  invisible(NULL)
}

############################################################
## 3. GEO series-matrix parser
############################################################

## GEO series-matrix metadata lines normally use:
##   !Sample_title<TAB>"sample 1"<TAB>"sample 2"...
## Some files may instead use an equals sign. This parser supports both.
parse_series_metadata_line <- function(line) {
  if (!grepl("^!Sample_", line, ignore.case = TRUE)) return(NULL)

  tab_pos <- regexpr("\\t", line)
  eq_pos  <- regexpr("=", line, fixed = TRUE)

  split_pos <- NA_integer_
  split_type <- NA_character_
  if (tab_pos[1L] > 0L && (eq_pos[1L] < 0L || tab_pos[1L] < eq_pos[1L])) {
    split_pos <- tab_pos[1L]
    split_type <- "TAB"
  } else if (eq_pos[1L] > 0L) {
    split_pos <- eq_pos[1L]
    split_type <- "EQUALS"
  }

  if (is.na(split_pos)) {
    return(list(raw_key = trimws(line), field = clean_field_name(line), values = character(), split_type = "NONE"))
  }

  lhs <- trimws(substr(line, 1L, split_pos - 1L))
  rhs <- substr(line, split_pos + 1L, nchar(line))
  vals <- if (split_type == "TAB") {
    strsplit(rhs, "\\t", fixed = FALSE)[[1L]]
  } else {
    ## Equals-form sample metadata may still contain tab-separated vectors.
    strsplit(trimws(rhs), "\\t", fixed = FALSE)[[1L]]
  }

  vals <- vapply(vals, sanitize_value, character(1))
  list(raw_key = lhs, field = clean_field_name(lhs), values = vals, split_type = split_type)
}

read_series_sample_metadata <- function(path, dataset_id, max_lines = MAX_HEADER_LINES) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  sample_lines <- character()
  lines_scanned <- 0L
  reached_table <- FALSE

  repeat {
    chunk <- readLines(con, n = 5000L, warn = FALSE)
    if (length(chunk) == 0L) break
    lines_scanned <- lines_scanned + length(chunk)
    sample_lines <- c(sample_lines, chunk[grepl("^!Sample_", chunk, ignore.case = TRUE)])
    if (any(grepl("^!series_matrix_table_begin", chunk, ignore.case = TRUE))) {
      reached_table <- TRUE
      break
    }
    if (lines_scanned >= max_lines) break
  }

  if (length(sample_lines) == 0L) {
    add_warning("SERIES_MATRIX", basename(path), "No !Sample_ metadata lines were detected.")
    return(list(long = data.table(), wide = data.table(), parse_summary = data.table(
      dataset_id = dataset_id,
      source_file = basename(path),
      lines_scanned = lines_scanned,
      reached_table_begin = reached_table,
      sample_metadata_lines = 0L,
      sample_count = 0L,
      parse_status = "NO_SAMPLE_LINES"
    )))
  }

  parsed <- lapply(sample_lines, parse_series_metadata_line)
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  fields <- vapply(parsed, function(z) z$field, character(1))

  accession_idx <- which(fields == "geo_accession")
  if (length(accession_idx) == 0L) {
    add_warning("SERIES_MATRIX", basename(path), "!Sample_geo_accession was not found.")
    return(list(long = data.table(), wide = data.table(), parse_summary = data.table(
      dataset_id = dataset_id,
      source_file = basename(path),
      lines_scanned = lines_scanned,
      reached_table_begin = reached_table,
      sample_metadata_lines = length(parsed),
      sample_count = 0L,
      parse_status = "NO_ACCESSION_VECTOR"
    )))
  }

  accessions <- parsed[[accession_idx[1L]]]$values
  n_samples <- length(accessions)
  if (n_samples == 0L) {
    add_warning("SERIES_MATRIX", basename(path), "The accession vector was empty.")
    return(list(long = data.table(), wide = data.table(), parse_summary = data.table(
      dataset_id = dataset_id,
      source_file = basename(path),
      lines_scanned = lines_scanned,
      reached_table_begin = reached_table,
      sample_metadata_lines = length(parsed),
      sample_count = 0L,
      parse_status = "EMPTY_ACCESSION_VECTOR"
    )))
  }

  field_counter <- list()
  long_records <- list()

  for (z in parsed) {
    vals <- z$values
    if (length(vals) == 0L) next

    field_counter[[z$field]] <- (field_counter[[z$field]] %||% 0L) + 1L
    instance <- field_counter[[z$field]]
    field_instance <- if (instance > 1L) paste0(z$field, "_", instance) else z$field

    if (length(vals) != n_samples) {
      add_warning(
        "SERIES_VECTOR_LENGTH",
        paste0(basename(path), "::", field_instance),
        paste0("Expected ", n_samples, " values but found ", length(vals), "; line retained in parse summary but omitted from sample table.")
      )
      next
    }

    long_records[[length(long_records) + 1L]] <- data.table(
      dataset_id = dataset_id,
      source_file = basename(path),
      sample_accession = accessions,
      field = field_instance,
      base_field = z$field,
      field_instance = instance,
      value = vals,
      parser = "series_matrix"
    )
  }

  long_dt <- safe_rbind(long_records)
  if (nrow(long_dt) == 0L) {
    wide_dt <- data.table()
  } else {
    wide_dt <- dcast(
      long_dt,
      dataset_id + source_file + sample_accession ~ field,
      value.var = "value",
      fun.aggregate = function(x) paste(unique(x[nzchar(x)]), collapse = " | ")
    )
  }

  list(
    long = long_dt,
    wide = wide_dt,
    parse_summary = data.table(
      dataset_id = dataset_id,
      source_file = basename(path),
      lines_scanned = lines_scanned,
      reached_table_begin = reached_table,
      sample_metadata_lines = length(parsed),
      sample_count = n_samples,
      parse_status = "OK"
    )
  )
}

############################################################
## 4. GEO family-SOFT parser
############################################################

parse_soft_assignment <- function(line) {
  pos <- regexpr("=", line, fixed = TRUE)
  if (pos[1L] < 0L) return(list(lhs = trimws(line), rhs = ""))
  lhs <- trimws(substr(line, 1L, pos[1L] - 1L))
  rhs <- sanitize_value(substr(line, pos[1L] + 1L, nchar(line)))
  list(lhs = lhs, rhs = rhs)
}

read_soft_sample_metadata <- function(path, dataset_id) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  records <- list()
  sample_counter <- 0L
  current_accession <- NA_character_
  field_counter <- list()
  lines_scanned <- 0L

  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (length(chunk) == 0L) break
    lines_scanned <- lines_scanned + length(chunk)

    for (line in chunk) {
      if (grepl("^\\^SAMPLE", line, ignore.case = TRUE)) {
        asg <- parse_soft_assignment(line)
        current_accession <- asg$rhs
        sample_counter <- sample_counter + 1L
        field_counter <- list()
        records[[length(records) + 1L]] <- data.table(
          dataset_id = dataset_id,
          source_file = basename(path),
          sample_accession = current_accession,
          field = "geo_accession",
          base_field = "geo_accession",
          field_instance = 1L,
          value = current_accession,
          parser = "family_soft"
        )
      } else if (!is.na(current_accession) && grepl("^!Sample_", line, ignore.case = TRUE)) {
        asg <- parse_soft_assignment(line)
        field <- clean_field_name(asg$lhs)
        field_counter[[field]] <- (field_counter[[field]] %||% 0L) + 1L
        instance <- field_counter[[field]]
        field_instance <- if (instance > 1L) paste0(field, "_", instance) else field
        records[[length(records) + 1L]] <- data.table(
          dataset_id = dataset_id,
          source_file = basename(path),
          sample_accession = current_accession,
          field = field_instance,
          base_field = field,
          field_instance = instance,
          value = asg$rhs,
          parser = "family_soft"
        )
      }
    }
  }

  long_dt <- safe_rbind(records)
  if (nrow(long_dt) == 0L) {
    wide_dt <- data.table()
    status <- "NO_SAMPLE_RECORDS"
  } else {
    wide_dt <- dcast(
      long_dt,
      dataset_id + source_file + sample_accession ~ field,
      value.var = "value",
      fun.aggregate = function(x) paste(unique(x[nzchar(x)]), collapse = " | ")
    )
    status <- "OK"
  }

  list(
    long = long_dt,
    wide = wide_dt,
    parse_summary = data.table(
      dataset_id = dataset_id,
      source_file = basename(path),
      lines_scanned = lines_scanned,
      sample_metadata_lines = nrow(long_dt),
      sample_count = uniqueN(long_dt$sample_accession),
      parse_status = status
    )
  )
}

############################################################
## 5. Parse all GEO metadata files
############################################################

series_files <- list.files(DATA_DIR, pattern = "series_matrix\\.txt\\.gz$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
soft_files   <- list.files(DATA_DIR, pattern = "family\\.soft\\.gz$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

series_long_list <- list()
series_wide_list <- list()
soft_long_list   <- list()
soft_wide_list   <- list()
parse_summary_list <- list()

for (f in series_files) {
  ds <- extract_gse(basename(f))
  log_msg("Parsing series-matrix sample metadata: ", basename(f))
  res <- tryCatch(
    read_series_sample_metadata(f, ds),
    error = function(e) {
      add_warning("SERIES_MATRIX_FATAL", basename(f), conditionMessage(e))
      list(long = data.table(), wide = data.table(), parse_summary = data.table(
        dataset_id = ds,
        source_file = basename(f),
        lines_scanned = NA_integer_,
        reached_table_begin = NA,
        sample_metadata_lines = NA_integer_,
        sample_count = 0L,
        parse_status = "FAILED"
      ))
    }
  )
  if (nrow(res$long) > 0L) series_long_list[[length(series_long_list) + 1L]] <- res$long
  if (nrow(res$wide) > 0L) series_wide_list[[length(series_wide_list) + 1L]] <- res$wide
  parse_summary_list[[length(parse_summary_list) + 1L]] <- res$parse_summary
}

for (f in soft_files) {
  ds <- extract_gse(basename(f))
  log_msg("Parsing family-SOFT sample metadata: ", basename(f))
  res <- tryCatch(
    read_soft_sample_metadata(f, ds),
    error = function(e) {
      add_warning("FAMILY_SOFT_FATAL", basename(f), conditionMessage(e))
      list(long = data.table(), wide = data.table(), parse_summary = data.table(
        dataset_id = ds,
        source_file = basename(f),
        lines_scanned = NA_integer_,
        sample_metadata_lines = NA_integer_,
        sample_count = 0L,
        parse_status = "FAILED"
      ))
    }
  )
  if (nrow(res$long) > 0L) soft_long_list[[length(soft_long_list) + 1L]] <- res$long
  if (nrow(res$wide) > 0L) soft_wide_list[[length(soft_wide_list) + 1L]] <- res$wide
  parse_summary_list[[length(parse_summary_list) + 1L]] <- res$parse_summary
}

series_long <- safe_rbind(series_long_list)
series_wide <- safe_rbind(series_wide_list)
soft_long   <- safe_rbind(soft_long_list)
soft_wide   <- safe_rbind(soft_wide_list)
parse_summary <- safe_rbind(parse_summary_list)

## Prefer series-matrix values for a dataset when available, while retaining
## SOFT-only fields and preserving the source parser for auditability.
all_long <- rbindlist(list(series_long, soft_long), fill = TRUE, use.names = TRUE)
if (nrow(all_long) > 0L) {
  all_long[, parser_priority := fifelse(parser == "series_matrix", 1L, 2L)]
  setorder(all_long, dataset_id, sample_accession, field, parser_priority)

  dedup_long <- all_long[, .SD[1L], by = .(dataset_id, sample_accession, field)]
  dedup_long[, parser_priority := NULL]

  combined_wide <- dcast(
    dedup_long,
    dataset_id + sample_accession ~ field,
    value.var = "value",
    fun.aggregate = function(x) paste(unique(x[nzchar(x)]), collapse = " | ")
  )
} else {
  dedup_long <- data.table()
  combined_wide <- data.table()
}

write_csv_safe(series_long, file.path(DIRS$tables, "01_GEO_series_sample_metadata_long.csv"))
write_csv_safe(series_wide, file.path(DIRS$tables, "02_GEO_series_sample_metadata_wide.csv"))
write_csv_safe(soft_long, file.path(DIRS$tables, "03_GEO_SOFT_sample_metadata_long.csv"))
write_csv_safe(soft_wide, file.path(DIRS$tables, "04_GEO_SOFT_sample_metadata_wide.csv"))
write_csv_safe(dedup_long, file.path(DIRS$tables, "05_GEO_combined_sample_metadata_long.csv"))
write_csv_safe(combined_wide, file.path(DIRS$tables, "06_GEO_combined_sample_metadata_wide.csv"))
write_csv_safe(parse_summary, file.path(DIRS$tables, "07_GEO_metadata_parse_summary.csv"))
saveRDS(combined_wide, file.path(DIRS$objects, "GEO_combined_sample_metadata_wide.rds"))

############################################################
## 6. Expand characteristics fields into key-value metadata
############################################################

characteristics_long <- data.table()
characteristics_wide <- data.table()

if (nrow(dedup_long) > 0L) {
  char_rows <- dedup_long[grepl("^characteristics_ch", base_field, ignore.case = TRUE)]

  if (nrow(char_rows) > 0L) {
    char_records <- vector("list", nrow(char_rows))
    for (i in seq_len(nrow(char_rows))) {
      value <- sanitize_value(char_rows$value[i])
      colon_pos <- regexpr(":", value, fixed = TRUE)
      if (colon_pos[1L] > 0L) {
        key <- trimws(substr(value, 1L, colon_pos[1L] - 1L))
        val <- trimws(substr(value, colon_pos[1L] + 1L, nchar(value)))
        parse_status <- "KEY_VALUE"
      } else {
        key <- paste0("unparsed_characteristic_", char_rows$field_instance[i])
        val <- value
        parse_status <- "NO_COLON"
      }
      key_clean <- clean_field_name(key)
      if (!nzchar(key_clean)) key_clean <- paste0("characteristic_", char_rows$field_instance[i])

      char_records[[i]] <- data.table(
        dataset_id = char_rows$dataset_id[i],
        sample_accession = char_rows$sample_accession[i],
        source_field = char_rows$field[i],
        characteristic_key_raw = key,
        characteristic_key = key_clean,
        characteristic_value = val,
        parse_status = parse_status
      )
    }
    characteristics_long <- safe_rbind(char_records)

    characteristics_wide <- dcast(
      characteristics_long,
      dataset_id + sample_accession ~ characteristic_key,
      value.var = "characteristic_value",
      fun.aggregate = function(x) paste(unique(x[nzchar(x)]), collapse = " | ")
    )
  }
}

write_csv_safe(characteristics_long, file.path(DIRS$tables, "08_GEO_characteristics_key_value_long.csv"))
write_csv_safe(characteristics_wide, file.path(DIRS$tables, "09_GEO_characteristics_key_value_wide.csv"))

############################################################
## 7. Build conservative grouping candidates
############################################################

collapse_nonmissing <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
  paste(unique(vals), collapse = " | ")
}

first_match <- function(text, patterns, labels) {
  text_low <- tolower(text %||% "")
  hits <- which(vapply(patterns, function(p) grepl(p, text_low, perl = TRUE), logical(1)))
  if (length(hits) == 0L) NA_character_ else labels[hits[1L]]
}

if (nrow(combined_wide) > 0L) {
  grouping_base <- copy(combined_wide)
  if (nrow(characteristics_wide) > 0L) {
    grouping_base <- merge(grouping_base, characteristics_wide,
                           by = c("dataset_id", "sample_accession"), all.x = TRUE)
  }

  grouping_cols <- setdiff(names(grouping_base), c("dataset_id", "sample_accession"))
  grouping_base[, metadata_text := apply(.SD, 1L, function(z) collapse_nonmissing(z)), .SDcols = grouping_cols]

  grouping_base[, condition_candidate := vapply(metadata_text, first_match, character(1),
    patterns = c(
      "hfpef", "heart failure with preserved", "control|non[- ]?failing|normal|healthy|sham",
      "dcm|dilated cardiomyopathy", "hfr?ef|reduced ejection", "tac"
    ),
    labels = c("HFpEF", "HFpEF", "Control_or_nonfailing", "DCM", "HFrEF_or_reduced_EF", "TAC_model")
  )]

  grouping_base[, drug_candidate := vapply(metadata_text, first_match, character(1),
    patterns = c("dapagliflozin|dapa", "empagliflozin|empa", "tya[-_ ]?018", "vehicle|veh", "untreated|no treatment"),
    labels = c("Dapagliflozin", "Empagliflozin", "TYA-018", "Vehicle", "Untreated")
  )]

  grouping_base[, diet_candidate := vapply(metadata_text, first_match, character(1),
    patterns = c("high[- ]?fat|hfd", "chow|control diet|normal diet|cd\\b", "ob/ob|obob"),
    labels = c("HFD", "Control_diet", "ob_ob")
  )]

  grouping_base[, macrophage_subset_candidate := vapply(metadata_text, first_match, character(1),
    patterns = c("ccr2\\+|ccr2 positive|ccr2pos", "ccr2[-−]|ccr2 negative|ccr2neg"),
    labels = c("Ccr2_positive", "Ccr2_negative")
  )]

  grouping_base[, sample_type_candidate := vapply(metadata_text, first_match, character(1),
    patterns = c("single nucleus|snrna|nuclei", "single cell|scrna", "bulk rna", "pbmc", "cardiomyocyte", "macrophage"),
    labels = c("Single_nucleus", "Single_cell", "Bulk_RNA", "PBMC", "Cardiomyocyte", "Macrophage")
  )]

  grouping_base[, manual_review_required :=
                  is.na(condition_candidate) & is.na(drug_candidate) &
                  is.na(diet_candidate) & is.na(macrophage_subset_candidate)]

  grouping_base[, grouping_note := fifelse(
    manual_review_required,
    "No reliable grouping token was inferred; inspect title/characteristics manually.",
    "Candidate labels inferred from public metadata text; confirm against the study design before analysis."
  )]

  grouping_candidates <- grouping_base[, .(
    dataset_id,
    sample_accession,
    title = if ("title" %in% names(grouping_base)) title else NA_character_,
    source_name_ch1 = if ("source_name_ch1" %in% names(grouping_base)) source_name_ch1 else NA_character_,
    condition_candidate,
    drug_candidate,
    diet_candidate,
    macrophage_subset_candidate,
    sample_type_candidate,
    manual_review_required,
    grouping_note,
    metadata_text
  )]
} else {
  grouping_candidates <- data.table()
}

write_csv_safe(grouping_candidates, file.path(DIRS$tables, "10_GEO_grouping_candidates_FOR_MANUAL_REVIEW.csv"))

############################################################
## 8. Sample-count checks
############################################################

sample_counts <- if (nrow(combined_wide) > 0L) {
  combined_wide[, .(observed_n = uniqueN(sample_accession)), by = dataset_id]
} else data.table(dataset_id = character(), observed_n = integer())

sample_count_check <- merge(
  as.data.table(EXPECTED_SAMPLE_COUNTS),
  sample_counts,
  by = "dataset_id",
  all = TRUE
)
sample_count_check[, status := fifelse(
  is.na(expected_n), "NO_PRESET_EXPECTATION",
  fifelse(is.na(observed_n), "NO_PARSED_SAMPLES",
          fifelse(expected_n == observed_n, "PASS", "CHECK"))
)]

write_csv_safe(sample_count_check, file.path(DIRS$tables, "11_GEO_sample_count_checks.csv"))

############################################################
## 9. SCP3342 donor-level AnnData obs inspection
############################################################

scp3342_obs_fields <- data.table()
scp3342_donor_summary <- data.table()
scp3342_donor_celltype <- data.table()
scp3342_status <- data.table(
  source_file = NA_character_,
  obs_count = NA_real_,
  donor_field = NA_character_,
  disease_field = NA_character_,
  cell_type_field = NA_character_,
  donor_count = NA_integer_,
  hfpef_donor_count = NA_integer_,
  control_donor_count = NA_integer_,
  inspection_status = "NOT_RUN"
)

read_h5_dataset_safe <- function(obj) {
  tryCatch(obj[], error = function(e) {
    tryCatch(obj$read(), error = function(e2) NULL)
  })
}

read_anndata_obs_column <- function(obs_group, field_name) {
  obj <- tryCatch(obs_group[[field_name]], error = function(e) NULL)
  if (is.null(obj)) return(NULL)

  ## Plain HDF5 dataset.
  if (inherits(obj, "H5D")) {
    vals <- read_h5_dataset_safe(obj)
    if (is.raw(vals)) vals <- rawToChar(vals)
    return(as.vector(vals))
  }

  ## AnnData categorical group: codes + categories.
  if (inherits(obj, "H5Group")) {
    names_here <- tryCatch(obj$ls(recursive = FALSE)$name, error = function(e) character())
    if (all(c("codes", "categories") %in% names_here)) {
      codes <- read_h5_dataset_safe(obj[["codes"]])
      cats  <- read_h5_dataset_safe(obj[["categories"]])
      if (is.null(codes) || is.null(cats)) return(NULL)
      codes <- as.integer(codes)
      cats <- as.character(cats)
      out <- rep(NA_character_, length(codes))
      valid <- !is.na(codes) & codes >= 0L & (codes + 1L) <= length(cats)
      out[valid] <- cats[codes[valid] + 1L]
      return(out)
    }

    ## Some encodings store values in a subgroup named "values".
    if ("values" %in% names_here) {
      vals <- read_h5_dataset_safe(obj[["values"]])
      return(as.vector(vals))
    }
  }

  NULL
}

h5ad_files <- list.files(DATA_DIR, pattern = "\\.h5ad$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if (INSPECT_SCP3342 && length(h5ad_files) > 0L) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    add_warning("SCP3342", "hdf5r", "hdf5r is unavailable; donor-level AnnData metadata was not inspected.")
    scp3342_status$inspection_status <- "HDF5R_UNAVAILABLE"
  } else {
    h5ad <- h5ad_files[grepl("HFpEF_snRNAseq", basename(h5ad_files), ignore.case = TRUE)][1L]
    if (is.na(h5ad) || !file.exists(h5ad)) h5ad <- h5ad_files[1L]
    log_msg("Inspecting SCP3342 donor-level obs metadata: ", basename(h5ad))

    tryCatch({
      h5 <- hdf5r::H5File$new(h5ad, mode = "r")
      obs <- h5[["obs"]]
      obs_fields <- obs$ls(recursive = FALSE)$name

      scp3342_obs_fields <- data.table(
        source_file = basename(h5ad),
        obs_field = obs_fields
      )

      donor_candidates <- c("donor_id", "donor", "patient_id", "patient", "individual_id", "subject_id", "sample_id")
      disease_candidates <- c("disease", "condition", "group", "diagnosis", "phenotype")
      celltype_candidates <- c("cell_type", "celltype", "cell_type_final", "annotation", "cell_annotation")

      pick_field <- function(candidates) {
        exact <- candidates[candidates %in% obs_fields]
        if (length(exact) > 0L) return(exact[1L])
        low_fields <- tolower(obs_fields)
        idx <- match(tolower(candidates), low_fields, nomatch = 0L)
        idx <- idx[idx > 0L]
        if (length(idx) > 0L) obs_fields[idx[1L]] else NA_character_
      }

      donor_field <- pick_field(donor_candidates)
      disease_field <- pick_field(disease_candidates)
      celltype_field <- pick_field(celltype_candidates)

      if (is.na(donor_field)) add_warning("SCP3342", "donor_field", "No donor identifier field was found in adata.obs.")
      if (is.na(disease_field)) add_warning("SCP3342", "disease_field", "No disease/condition field was found in adata.obs.")
      if (is.na(celltype_field)) add_warning("SCP3342", "cell_type_field", "No cell-type field was found in adata.obs.")

      donor <- if (!is.na(donor_field)) read_anndata_obs_column(obs, donor_field) else NULL
      disease <- if (!is.na(disease_field)) read_anndata_obs_column(obs, disease_field) else NULL
      celltype <- if (!is.na(celltype_field)) read_anndata_obs_column(obs, celltype_field) else NULL

      lengths_found <- c(length(donor), length(disease), length(celltype))
      obs_n <- max(lengths_found, na.rm = TRUE)
      if (!is.finite(obs_n)) obs_n <- NA_real_

      if (!is.null(donor) && !is.null(disease) && length(donor) == length(disease)) {
        obs_dt <- data.table(
          donor_id = as.character(donor),
          disease = as.character(disease)
        )
        if (!is.null(celltype) && length(celltype) == nrow(obs_dt)) {
          obs_dt[, cell_type := as.character(celltype)]
        } else {
          obs_dt[, cell_type := NA_character_]
        }

        obs_dt <- obs_dt[!is.na(donor_id) & nzchar(donor_id)]
        scp3342_donor_summary <- obs_dt[, .(
          n_nuclei = .N,
          n_cell_types = uniqueN(cell_type[!is.na(cell_type) & nzchar(cell_type)])
        ), by = .(donor_id, disease)]
        setorder(scp3342_donor_summary, disease, donor_id)

        scp3342_donor_celltype <- obs_dt[!is.na(cell_type) & nzchar(cell_type), .N,
                                        by = .(donor_id, disease, cell_type)]
        setnames(scp3342_donor_celltype, "N", "n_nuclei")
        setorder(scp3342_donor_celltype, disease, donor_id, -n_nuclei)

        disease_low <- tolower(scp3342_donor_summary$disease)
        hfpef_n <- uniqueN(scp3342_donor_summary$donor_id[grepl("hfpef|preserved", disease_low)])
        control_n <- uniqueN(scp3342_donor_summary$donor_id[grepl("control|non[- ]?failing|normal|healthy", disease_low)])

        scp3342_status <- data.table(
          source_file = basename(h5ad),
          obs_count = obs_n,
          donor_field = donor_field,
          disease_field = disease_field,
          cell_type_field = celltype_field,
          donor_count = uniqueN(scp3342_donor_summary$donor_id),
          hfpef_donor_count = hfpef_n,
          control_donor_count = control_n,
          inspection_status = "OK"
        )
      } else {
        add_warning("SCP3342", "obs_lengths", "Donor and disease columns were missing or had incompatible lengths.")
        scp3342_status <- data.table(
          source_file = basename(h5ad),
          obs_count = obs_n,
          donor_field = donor_field,
          disease_field = disease_field,
          cell_type_field = celltype_field,
          donor_count = NA_integer_,
          hfpef_donor_count = NA_integer_,
          control_donor_count = NA_integer_,
          inspection_status = "INCOMPLETE_OBS_COLUMNS"
        )
      }
      try(h5$close_all(), silent = TRUE)
    }, error = function(e) {
      if (exists("h5", inherits = FALSE)) try(h5$close_all(), silent = TRUE)
      add_warning("SCP3342_FATAL", basename(h5ad), conditionMessage(e))
      scp3342_status$inspection_status <- "FAILED"
    })
  }
}

write_csv_safe(scp3342_obs_fields, file.path(DIRS$tables, "12_SCP3342_obs_fields.csv"))
write_csv_safe(scp3342_donor_summary, file.path(DIRS$tables, "13_SCP3342_donor_summary.csv"))
write_csv_safe(scp3342_donor_celltype, file.path(DIRS$tables, "14_SCP3342_donor_celltype_counts.csv"))
write_csv_safe(scp3342_status, file.path(DIRS$tables, "15_SCP3342_metadata_status.csv"))

############################################################
## 10. Metadata lock summary and manual-review flags
############################################################

manual_review_summary <- data.table()
if (nrow(grouping_candidates) > 0L) {
  manual_review_summary <- grouping_candidates[, .(
    parsed_samples = .N,
    samples_requiring_manual_review = sum(manual_review_required, na.rm = TRUE),
    proportion_requiring_manual_review = mean(manual_review_required, na.rm = TRUE),
    inferred_condition_n = sum(!is.na(condition_candidate)),
    inferred_drug_n = sum(!is.na(drug_candidate)),
    inferred_diet_n = sum(!is.na(diet_candidate)),
    inferred_macrophage_subset_n = sum(!is.na(macrophage_subset_candidate))
  ), by = dataset_id]
}
write_csv_safe(manual_review_summary, file.path(DIRS$tables, "16_manual_review_summary_by_dataset.csv"))

metadata_lock <- merge(
  sample_count_check,
  manual_review_summary,
  by = "dataset_id",
  all = TRUE
)
metadata_lock[, metadata_lock_status := fifelse(
  status == "PASS" & (is.na(samples_requiring_manual_review) | samples_requiring_manual_review == 0L),
  "READY_AFTER_MANUAL_SPOT_CHECK",
  fifelse(status == "PASS", "SAMPLE_COUNT_OK_GROUPS_NEED_REVIEW",
          fifelse(status == "NO_PRESET_EXPECTATION", "PARSED_NO_PRESET_COUNT", "CHECK_REQUIRED"))
)]
write_csv_safe(metadata_lock, file.path(DIRS$tables, "17_metadata_lock_status.csv"))

############################################################
## 11. Workbook, methods, README, and exact script archival
############################################################

warnings_dt <- if (length(warning_records) > 0L) {
  rbindlist(warning_records, fill = TRUE)
} else {
  data.table(timestamp = character(), category = character(), item = character(), message = character())
}
write_csv_safe(warnings_dt, file.path(DIRS$tables, "18_warnings_and_nonfatal_issues.csv"))

wb <- createWorkbook()
write_sheet_safe(wb, "Parse_summary", parse_summary)
write_sheet_safe(wb, "Combined_metadata", combined_wide)
write_sheet_safe(wb, "Characteristics", characteristics_wide)
write_sheet_safe(wb, "Grouping_candidates", grouping_candidates)
write_sheet_safe(wb, "Sample_count_checks", sample_count_check)
write_sheet_safe(wb, "SCP3342_status", scp3342_status)
write_sheet_safe(wb, "SCP3342_donors", scp3342_donor_summary)
write_sheet_safe(wb, "Metadata_lock", metadata_lock)
write_sheet_safe(wb, "Warnings", warnings_dt)
saveWorkbook(wb, file.path(DIRS$tables, "Stage1_Metadata_Patch_Summary.xlsx"), overwrite = TRUE)

script_path <- get_script_path()
script_copy_status <- "NOT_FOUND"
if (!is.na(script_path) && file.exists(script_path)) {
  file.copy(script_path, file.path(DIRS$methods, basename(script_path)), overwrite = TRUE)
  file.copy(script_path, file.path(DIRS$check, basename(script_path)), overwrite = TRUE)
  script_copy_status <- "COPIED"
} else {
  add_warning(
    "SCRIPT_ARCHIVE",
    "executed_script",
    "The script path could not be resolved automatically. Manually copy the executed .R file into 05_methods and 06_review_check before sharing the CHECK package."
  )
  writeLines(
    c(
      "The exact executed R script path could not be resolved automatically.",
      "Before external review, manually copy the executed .R file into this folder."
    ),
    file.path(DIRS$methods, "SCRIPT_COPY_REQUIRED.txt")
  )
}

methods_text <- c(
  "Stage 1 Metadata Patch FIXED v1",
  "",
  "Purpose:",
  "- Correct GEO sample-level metadata parsing from series-matrix and family-SOFT files.",
  "- Build one-row-per-sample metadata and structured characteristics tables.",
  "- Generate conservative grouping candidates requiring manual confirmation.",
  "- Inspect SCP3342 donor, disease, and cell-type metadata without loading the expression matrix.",
  "",
  "Important boundaries:",
  "- Grouping candidates are inferred from metadata text and must be manually checked.",
  "- No biological analysis was performed.",
  "- Individual cells/nuclei are not treated as independent biological replicates.",
  "",
  paste0("Executed script archive status: ", script_copy_status),
  paste0("Run started: ", format(START_TIME, "%Y-%m-%d %H:%M:%S")),
  paste0("Run completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(methods_text, file.path(DIRS$methods, "metadata_patch_methods_and_boundaries.txt"), useBytes = TRUE)

readme_text <- c(
  "HFpEF Stage 1 Metadata Patch FIXED v1",
  "",
  paste0("Input directory: ", DATA_DIR),
  paste0("Previous Stage 1 directory: ", PREV_STAGE),
  paste0("Output directory: ", OUT_DIR),
  "",
  "Core outputs:",
  "01-06: GEO sample metadata in long and wide formats",
  "08-09: structured characteristics key-value tables",
  "10: conservative grouping candidates for manual review",
  "11: expected versus observed GEO sample-count checks",
  "13-15: SCP3342 donor-level metadata summaries",
  "17: metadata-lock status by dataset",
  "Stage1_Metadata_Patch_Summary.xlsx: compact review workbook",
  "",
  "Do not start differential expression until grouping candidates and metadata-lock flags have been manually reviewed."
)
writeLines(readme_text, file.path(OUT_DIR, "README_stage1_metadata_patch.txt"), useBytes = TRUE)

############################################################
## 12. Run status and CHECK package
############################################################

critical_issues <- 0L
if (nrow(parse_summary) == 0L || all(parse_summary$sample_count == 0L, na.rm = TRUE)) critical_issues <- critical_issues + 1L
if (nrow(sample_count_check[status == "NO_PARSED_SAMPLES"]) > 0L) critical_issues <- critical_issues + 1L
if (INSPECT_SCP3342 && nrow(scp3342_status) > 0L && scp3342_status$inspection_status[1L] == "FAILED") critical_issues <- critical_issues + 1L

END_TIME <- Sys.time()
run_status <- data.table(
  stage = STAGE_NAME,
  start_time = format(START_TIME, "%Y-%m-%d %H:%M:%S"),
  end_time = format(END_TIME, "%Y-%m-%d %H:%M:%S"),
  elapsed_minutes = round(as.numeric(difftime(END_TIME, START_TIME, units = "mins")), 2),
  series_files_parsed = length(series_files),
  soft_files_parsed = length(soft_files),
  geo_samples_parsed = if (nrow(combined_wide) > 0L) uniqueN(combined_wide$sample_accession) else 0L,
  datasets_with_samples = if (nrow(combined_wide) > 0L) uniqueN(combined_wide$dataset_id) else 0L,
  scp3342_donors = scp3342_status$donor_count[1L] %||% NA_integer_,
  warning_count = nrow(warnings_dt),
  critical_issue_count = critical_issues,
  script_copy_status = script_copy_status,
  overall_status = ifelse(critical_issues == 0L, "COMPLETED_REVIEW_REQUIRED", "COMPLETED_WITH_CRITICAL_ISSUES")
)
write_csv_safe(run_status, file.path(DIRS$tables, "19_run_status.csv"))

## Refresh warnings after any late script-copy warning.
warnings_dt <- if (length(warning_records) > 0L) {
  rbindlist(warning_records, fill = TRUE)
} else {
  data.table(timestamp = character(), category = character(), item = character(), message = character())
}
write_csv_safe(warnings_dt, file.path(DIRS$tables, "18_warnings_and_nonfatal_issues.csv"))

## Copy compact review files.
check_files <- c(
  file.path(DIRS$tables, "06_GEO_combined_sample_metadata_wide.csv"),
  file.path(DIRS$tables, "09_GEO_characteristics_key_value_wide.csv"),
  file.path(DIRS$tables, "10_GEO_grouping_candidates_FOR_MANUAL_REVIEW.csv"),
  file.path(DIRS$tables, "11_GEO_sample_count_checks.csv"),
  file.path(DIRS$tables, "13_SCP3342_donor_summary.csv"),
  file.path(DIRS$tables, "15_SCP3342_metadata_status.csv"),
  file.path(DIRS$tables, "17_metadata_lock_status.csv"),
  file.path(DIRS$tables, "18_warnings_and_nonfatal_issues.csv"),
  file.path(DIRS$tables, "19_run_status.csv"),
  file.path(DIRS$tables, "Stage1_Metadata_Patch_Summary.xlsx"),
  file.path(DIRS$logs, "stage1_metadata_patch.log"),
  file.path(DIRS$methods, "metadata_patch_methods_and_boundaries.txt"),
  file.path(OUT_DIR, "README_stage1_metadata_patch.txt")
)
check_files <- check_files[file.exists(check_files)]
for (f in check_files) file.copy(f, file.path(DIRS$check, basename(f)), overwrite = TRUE)

## Zip the contents of 06_review_check, not the parent directory.
zip_inputs <- list.files(DIRS$check, full.names = TRUE, recursive = TRUE)
if (length(zip_inputs) == 0L) stop("No files were available for the CHECK package.")

old_wd <- getwd()
setwd(DIRS$check)
zip_result <- tryCatch(
  zip::zipr(
    zipfile = CHECK_ZIP,
    files = list.files(".", recursive = TRUE, all.files = FALSE),
    recurse = TRUE,
    include_directories = FALSE
  ),
  error = function(e) {
    setwd(old_wd)
    stop("Failed to create CHECK package: ", conditionMessage(e))
  }
)
setwd(old_wd)

log_msg("Stage 1 metadata patch completed.")
log_msg("Elapsed minutes: ", run_status$elapsed_minutes)
log_msg("GEO samples parsed: ", run_status$geo_samples_parsed)
log_msg("Datasets with parsed samples: ", run_status$datasets_with_samples)
log_msg("SCP3342 donors detected: ", run_status$scp3342_donors)
log_msg("Warnings: ", run_status$warning_count)
log_msg("Overall status: ", run_status$overall_status)
log_msg("CHECK package: ", CHECK_ZIP)

cat("\n============================================================\n")
cat("Stage 1 Metadata Patch FIXED v1 completed\n")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK:  ", CHECK_ZIP, "\n", sep = "")
cat("Status: ", run_status$overall_status, "\n", sep = "")
cat("Review 10_GEO_grouping_candidates_FOR_MANUAL_REVIEW.csv before Stage 2.\n")
cat("============================================================\n")
