############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 1 Metadata Lock FIXED v3
##
## FIXED v3 change:
##   - Detects the actual script path automatically when run with source().
##   - The script no longer has to be stored under one exact local pathname.
##
## Project root:
##   <HFPEF_PROJECT_DIR>
##
## Read-only raw-data directory:
##   <HFPEF_PROJECT_DIR>/0.GEO
##
## Required previous output:
##   <HFPEF_PROJECT_DIR>/01_stage1_metadata_patch_FIXED_v1
##
## Output directory:
##   <HFPEF_PROJECT_DIR>/01_stage1_metadata_lock_FIXED_v3
##
## Review package:
##   <HFPEF_PROJECT_DIR>/01_stage1_metadata_lock_FIXED_v3_CHECK.zip
##
## Purpose:
##   1) Replace broad keyword-based grouping with dataset-specific rules.
##   2) Lock one definitive sample manifest for all local GEO datasets.
##   3) Validate every expected group count before biological analysis.
##   4) Map SCP3342 ontology codes to 19 HFpEF and 24 control donors.
##   5) Define each dataset's analytical role, unit, and planned contrasts.
##   6) Archive the exact executed script in the output and CHECK package.
##
## Scientific boundary:
##   - This stage locks metadata only.
##   - It does not normalize expression data or perform biological analysis.
##   - GSE270896 model groups are locked as genotype/background x surgery;
##     no HFpEF/HFrEF phenotype is assigned at this stage.
##   - GSE223527 is group-level public data; cells must not be treated as
##     independent patients.
##
## Required run command:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage1_Metadata_Lock_FIXED_v3.R",
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
PREV_PATCH  <- file.path(PROJECT_DIR, "01_stage1_metadata_patch_FIXED_v1")

## Automatically identify the exact script being executed.
## This works when the file is run with source(), regardless of whether it is
## stored in the project folder or the Downloads folder.
detect_script_file <- function() {
  frames <- sys.frames()
  candidates <- character()

  for (i in rev(seq_along(frames))) {
    ofile <- tryCatch(frames[[i]]$ofile, error = function(e) NULL)
    if (!is.null(ofile) && length(ofile) == 1L && nzchar(ofile)) {
      candidates <- c(candidates, ofile)
    }
  }

  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0L) {
    candidates <- c(candidates, sub("^--file=", "", file_arg[1L]))
  }

  candidates <- unique(candidates)
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0L) {
    return(NA_character_)
  }

  normalizePath(
    candidates[1L],
    winslash = "/",
    mustWork = TRUE
  )
}

SCRIPT_FILE <- detect_script_file()

STAGE_NAME <- "01_stage1_metadata_lock_FIXED_v3"
OUT_DIR   <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP <- file.path(PROJECT_DIR, paste0(STAGE_NAME, "_CHECK.zip"))

ALLOW_OVERWRITE <- FALSE
INSPECT_SCP3342 <- TRUE

PREV_METADATA_FILE <- file.path(
  PREV_PATCH,
  "01_tables",
  "06_GEO_combined_sample_metadata_wide.csv"
)

H5AD_PATTERN <- "^HFpEF_snRNAseq_single_cell_portal_10\\.14\\.2025\\.h5ad$"

############################################################
## 1. Output setup and packages
############################################################

if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
if (!dir.exists(DATA_DIR)) stop("DATA_DIR does not exist: ", DATA_DIR)
if (!file.exists(PREV_METADATA_FILE)) {
  stop(
    "Required parsed metadata file does not exist:\n",
    PREV_METADATA_FILE,
    "\nRun Stage 1 Metadata Patch FIXED v1 first."
  )
}
SCRIPT_ARCHIVE_AVAILABLE <- (
  length(SCRIPT_FILE) == 1L &&
  !is.na(SCRIPT_FILE) &&
  nzchar(SCRIPT_FILE) &&
  file.exists(SCRIPT_FILE)
)

if (!SCRIPT_ARCHIVE_AVAILABLE) {
  warning(
    paste0(
      "The running script path could not be detected automatically. ",
      "The metadata lock will continue, but the script-copy check will be ",
      "reported as unavailable. Run the file with source() rather than ",
      "pasting all lines directly into the console."
    ),
    call. = FALSE
  )
}

if (dir.exists(OUT_DIR) && !ALLOW_OVERWRITE) {
  stop(
    "Output folder already exists and overwrite protection is active:\n",
    OUT_DIR,
    "\nUse a new version name or back up the old output."
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

LOG_FILE  <- file.path(DIRS$logs, "stage1_metadata_lock.log")
WARN_FILE <- file.path(DIRS$logs, "stage1_metadata_lock_warnings.log")
START_TIME <- Sys.time()

log_msg <- function(..., level = "INFO") {
  txt <- paste0(..., collapse = "")
  line <- sprintf(
    "[%s] [%s] %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    level,
    txt
  )
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
  cat(
    sprintf(
      "[%s] [%s] %s: %s\n",
      rec$timestamp,
      category,
      item,
      message
    ),
    file = WARN_FILE,
    append = TRUE
  )
  log_msg(category, " | ", item, " | ", message, level = "WARN")
  invisible(rec)
}

ensure_cran <- function(pkgs, required = TRUE) {
  missing <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    log_msg(
      "Installing missing CRAN package(s): ",
      paste(missing, collapse = ", ")
    )
    try(
      install.packages(
        missing,
        repos = "https://cloud.r-project.org",
        dependencies = TRUE
      ),
      silent = TRUE
    )
  }

  still_missing <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(still_missing) > 0L) {
    msg <- paste(
      "Package(s) unavailable:",
      paste(still_missing, collapse = ", ")
    )
    if (required) {
      stop(msg)
    } else {
      add_warning(
        "PACKAGE",
        paste(still_missing, collapse = ";"),
        msg
      )
    }
  }
  invisible(setdiff(pkgs, still_missing))
}

ensure_cran(c("data.table", "openxlsx", "zip"), required = TRUE)
if (INSPECT_SCP3342) ensure_cran("hdf5r", required = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

log_msg("Stage 1 metadata lock started.")
log_msg("PROJECT_DIR: ", PROJECT_DIR)
log_msg("DATA_DIR: ", DATA_DIR)
log_msg("PREV_METADATA_FILE: ", PREV_METADATA_FILE)
log_msg("OUT_DIR: ", OUT_DIR)
log_msg(
  "Detected SCRIPT_FILE: ",
  ifelse(SCRIPT_ARCHIVE_AVAILABLE, SCRIPT_FILE, "NOT_DETECTED")
)

############################################################
## 2. Utility functions
############################################################

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\r|\\n", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

lower_text <- function(x) tolower(normalize_text(x))

extract_first_regex <- function(x, pattern, group = 1L) {
  x <- as.character(x)
  m <- regexec(pattern, x, perl = TRUE, ignore.case = TRUE)
  hit <- regmatches(x, m)
  vapply(
    hit,
    function(z) {
      if (length(z) > group) z[group + 1L] else NA_character_
    },
    character(1)
  )
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
    writeData(
      wb,
      sheet,
      data.frame(note = "No records generated.")
    )
    return(invisible(NULL))
  }

  y <- as.data.frame(x, stringsAsFactors = FALSE)
  char_cols <- vapply(y, is.character, logical(1))
  if (any(char_cols)) {
    y[char_cols] <- lapply(
      y[char_cols],
      function(z) substr(z, 1L, 30000L)
    )
  }

  writeData(wb, sheet, y)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(
    wb,
    sheet,
    cols = seq_len(min(ncol(y), 40L)),
    widths = "auto"
  )
  invisible(NULL)
}

read_h5_dataset_safe <- function(obj) {
  tryCatch(
    obj[],
    error = function(e) {
      tryCatch(obj$read(), error = function(e2) NULL)
    }
  )
}

read_anndata_obs_column <- function(obs_group, field_name) {
  obj <- tryCatch(obs_group[[field_name]], error = function(e) NULL)
  if (is.null(obj)) return(NULL)

  if (inherits(obj, "H5D")) {
    vals <- read_h5_dataset_safe(obj)
    if (is.raw(vals)) vals <- rawToChar(vals)
    return(as.vector(vals))
  }

  if (inherits(obj, "H5Group")) {
    names_here <- tryCatch(
      obj$ls(recursive = FALSE)$name,
      error = function(e) character()
    )

    if (all(c("codes", "categories") %in% names_here)) {
      codes <- read_h5_dataset_safe(obj[["codes"]])
      cats  <- read_h5_dataset_safe(obj[["categories"]])

      if (is.null(codes) || is.null(cats)) return(NULL)

      codes <- as.integer(codes)
      cats  <- as.character(cats)

      out <- rep(NA_character_, length(codes))
      valid <- !is.na(codes) &
        codes >= 0L &
        (codes + 1L) <= length(cats)
      out[valid] <- cats[codes[valid] + 1L]
      return(out)
    }

    if ("values" %in% names_here) {
      vals <- read_h5_dataset_safe(obj[["values"]])
      return(as.vector(vals))
    }
  }

  NULL
}

############################################################
## 3. Read previously parsed GEO metadata
############################################################

meta <- fread(
  PREV_METADATA_FILE,
  encoding = "UTF-8",
  na.strings = c("", "NA", "NaN")
)

required_columns <- c(
  "dataset_id",
  "sample_accession",
  "title"
)
missing_columns <- setdiff(required_columns, names(meta))
if (length(missing_columns) > 0L) {
  stop(
    "Previous metadata table is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(meta[, .(dataset_id, sample_accession)]) > 0L) {
  dup <- meta[
    duplicated(meta[, .(dataset_id, sample_accession)]) |
      duplicated(
        meta[, .(dataset_id, sample_accession)],
        fromLast = TRUE
      )
  ]
  fwrite(
    dup,
    file.path(DIRS$tables, "DUPLICATE_SAMPLE_RECORDS_FATAL.csv")
  )
  stop(
    "Duplicate dataset_id + sample_accession records were detected. ",
    "See DUPLICATE_SAMPLE_RECORDS_FATAL.csv."
  )
}

char_cols <- grep(
  "^characteristics_ch1",
  names(meta),
  value = TRUE,
  ignore.case = TRUE
)
supp_cols <- grep(
  "^supplementary_file",
  names(meta),
  value = TRUE,
  ignore.case = TRUE
)

if (length(char_cols) == 0L) {
  stop("No characteristics_ch1 columns were found in the parsed metadata.")
}

meta[, raw_characteristics := apply(
  .SD,
  1L,
  function(z) paste(normalize_text(z), collapse = " | ")
), .SDcols = char_cols]

if (length(supp_cols) > 0L) {
  meta[, raw_supplementary_files := apply(
    .SD,
    1L,
    function(z) paste(normalize_text(z), collapse = " | ")
  ), .SDcols = supp_cols]
} else {
  meta[, raw_supplementary_files := ""]
}

meta[, title_clean := normalize_text(title)]
meta[, title_lower := lower_text(title_clean)]
meta[, characteristics_lower := lower_text(raw_characteristics)]
meta[, supplementary_lower := lower_text(raw_supplementary_files)]

log_msg(
  "Loaded ",
  nrow(meta),
  " GEO sample records across ",
  uniqueN(meta$dataset_id),
  " datasets."
)

############################################################
## 4. Base locked-manifest fields
############################################################

base_cols <- c(
  "dataset_id",
  "sample_accession",
  "title_clean",
  "source_name_ch1",
  "organism_ch1",
  "library_strategy",
  "type",
  "raw_characteristics",
  "raw_supplementary_files"
)
for (nm in setdiff(base_cols, names(meta))) {
  meta[, (nm) := NA_character_]
}

locked <- meta[, .(
  dataset_id,
  sample_accession,
  original_title = title_clean,
  source_name = normalize_text(source_name_ch1),
  organism = normalize_text(organism_ch1),
  library_strategy = normalize_text(library_strategy),
  record_type = normalize_text(type),
  raw_characteristics,
  raw_supplementary_files,
  condition = NA_character_,
  genotype_or_background = NA_character_,
  diet_or_model = NA_character_,
  drug = NA_character_,
  macrophage_subset = NA_character_,
  surgery_or_stressor = NA_character_,
  assay_modality = NA_character_,
  technology = NA_character_,
  sex = NA_character_,
  reported_replicate = NA_character_,
  source_file_replicate = NA_character_,
  sample_subject_id = NA_character_,
  group_id = NA_character_,
  planned_role = NA_character_,
  primary_analysis_unit = NA_character_,
  include_in_expression_analysis = TRUE,
  claim_boundary = NA_character_,
  metadata_rule = NA_character_,
  lock_status = "UNLOCKED"
)]

locked[, title_lower := lower_text(original_title)]
locked[, char_lower := lower_text(raw_characteristics)]
locked[, supp_lower := lower_text(raw_supplementary_files)]

############################################################
## 5. Dataset-specific deterministic metadata rules
############################################################

## 5.1 GSE237156: sorted Ccr2 macrophage pharmacotranscriptomics
idx <- locked$dataset_id == "GSE237156"
if (any(idx)) {
  locked[idx, macrophage_subset := fcase(
    grepl("^ccr2\\+", title_lower), "Ccr2_positive",
    grepl("^ccr2-", title_lower), "Ccr2_negative",
    default = NA_character_
  )]
  locked[idx, diet_or_model := fcase(
    grepl("\\bhfd\\b", title_lower), "HFD",
    grepl("\\bcd\\b", title_lower), "CD",
    default = NA_character_
  )]
  locked[idx, drug := fcase(
    grepl("\\bdapa\\b|dapagliflozin", title_lower), "Dapagliflozin",
    grepl("\\bvehicle\\b|\\bveh\\b", title_lower), "Vehicle",
    default = NA_character_
  )]
  locked[idx, genotype_or_background := "ApoE_KO"]
  locked[idx, condition := fifelse(
    diet_or_model == "HFD",
    "Metabolic_diet_stress",
    "Diet_control"
  )]
  locked[idx, reported_replicate := extract_first_regex(
    original_title,
    "([0-9]+)\\s*$"
  )]
  locked[idx, source_file_replicate := extract_first_regex(
    raw_supplementary_files,
    "rep([0-9]+)"
  )]
  locked[idx, assay_modality := "Bulk_RNA_seq"]
  locked[idx, technology := "Sorted_macrophage_RNA_seq"]
  locked[idx, group_id := paste(
    macrophage_subset,
    diet_or_model,
    drug,
    sep = "__"
  )]
  locked[idx, planned_role := "Primary_drug_opposed_program_discovery"]
  locked[idx, primary_analysis_unit := "Biological_sample"]
  locked[idx, claim_boundary :=
    "Drug-opposed transcriptomic discovery; not direct molecular causality"]
  locked[idx, metadata_rule :=
    "Title-locked Ccr2 status, diet, and drug; GSM accession is authoritative"]
}

## 5.2 GSE208425: cardiac CD45+ immune single-cell context
idx <- locked$dataset_id == "GSE208425"
if (any(idx)) {
  locked[idx, genotype_or_background := fcase(
    grepl("apoe[_ ]ko", title_lower), "ApoE_KO",
    grepl("^wt|wild type", title_lower), "WT",
    default = NA_character_
  )]
  locked[idx, diet_or_model := fcase(
    grepl("_hfd_|\\bhfd\\b", title_lower), "HFD",
    grepl("_cd_|\\bcd\\b", title_lower), "CD",
    default = NA_character_
  )]
  locked[idx, condition := paste(
    genotype_or_background,
    diet_or_model,
    sep = "__"
  )]
  locked[idx, reported_replicate := extract_first_regex(
    original_title,
    "replicate\\s*([0-9]+)"
  )]
  locked[idx, assay_modality := "scRNA_seq"]
  locked[idx, technology := "Cardiac_CD45_positive_cells"]
  locked[idx, group_id := condition]
  locked[idx, planned_role := "Internal_immune_context"]
  locked[idx, primary_analysis_unit := "Biological_sample"]
  locked[idx, claim_boundary :=
    "Same project family as drug anchor; not independent validation"]
  locked[idx, metadata_rule :=
    "Title-locked genotype and diet"]
}

## 5.3 GSE209548: cardiomyocyte bulk transcriptomic support
idx <- locked$dataset_id == "GSE209548"
if (any(idx)) {
  locked[idx, genotype_or_background := fcase(
    grepl("apoe\\s*ko", title_lower), "ApoE_KO",
    grepl("^wt|^wt\\b", title_lower), "WT",
    default = NA_character_
  )]
  locked[idx, diet_or_model := fcase(
    grepl("\\bhfd\\b", title_lower), "HFD",
    grepl("\\bcd\\b", title_lower), "CD",
    default = NA_character_
  )]
  locked[idx, condition := paste(
    genotype_or_background,
    diet_or_model,
    sep = "__"
  )]
  locked[idx, reported_replicate := extract_first_regex(
    original_title,
    "replicate\\s*([0-9]+)"
  )]
  locked[idx, assay_modality := "Bulk_RNA_seq"]
  locked[idx, technology := "Cardiomyocyte_RNA_seq"]
  locked[idx, group_id := condition]
  locked[idx, planned_role := "Internal_cardiomyocyte_support"]
  locked[idx, primary_analysis_unit := "Biological_sample"]
  locked[idx, claim_boundary :=
    "Auxiliary myocardial response; not macrophage causality or independent validation"]
  locked[idx, metadata_rule :=
    "Title-locked genotype and diet"]
}

## 5.4 GSE236584: matched cardiac bulk cohort
idx <- locked$dataset_id == "GSE236584"
if (any(idx)) {
  locked[idx, condition := fifelse(
    grepl("^hfpef", title_lower),
    "HFpEF",
    fifelse(grepl("^con", title_lower), "Control", NA_character_)
  )]
  locked[idx, diet_or_model := fifelse(
    condition == "HFpEF",
    "HFD_LNAME",
    fifelse(condition == "Control", "Normal_diet_water", NA_character_)
  )]
  locked[idx, genotype_or_background := "C57BL6N_WT"]
  locked[idx, reported_replicate := extract_first_regex(
    original_title,
    "replicate[_ ]([0-9]+)"
  )]
  locked[idx, assay_modality := "Bulk_RNA_seq"]
  locked[idx, technology := "Ventricular_tissue_bulk"]
  locked[idx, group_id := condition]
  locked[idx, planned_role := "Matched_bulk_orthogonal_support"]
  locked[idx, primary_analysis_unit := "Biological_sample"]
  locked[idx, claim_boundary :=
    "Same-study orthogonal bulk support; not independent validation"]
  locked[idx, metadata_rule :=
    "Title-locked HFpEF versus CON"]
}

## 5.5 GSE236585: primary cardiac single-cell discovery cohort
idx <- locked$dataset_id == "GSE236585"
if (any(idx)) {
  locked[idx, condition := fifelse(
    grepl("^hfpef", title_lower),
    "HFpEF",
    fifelse(grepl("^con", title_lower), "Control", NA_character_)
  )]
  locked[idx, diet_or_model := fifelse(
    condition == "HFpEF",
    "HFD_LNAME",
    fifelse(condition == "Control", "Normal_diet_water", NA_character_)
  )]
  locked[idx, genotype_or_background := "C57BL6N_WT"]
  locked[idx, reported_replicate := extract_first_regex(
    original_title,
    "replicate[_ ]([0-9]+)"
  )]
  locked[idx, assay_modality := "scRNA_seq"]
  locked[idx, technology := "Ventricular_tissue_single_cell"]
  locked[idx, group_id := condition]
  locked[idx, planned_role := "Primary_cardiac_single_cell_discovery"]
  locked[idx, primary_analysis_unit := "Biological_sample_pseudobulk"]
  locked[idx, claim_boundary :=
    "Discovery cohort; cells are not independent biological replicates"]
  locked[idx, metadata_rule :=
    "Title-locked HFpEF versus CON"]
}

## 5.6 GSE223527: public group-level human PBMC GEX/ADT
idx <- locked$dataset_id == "GSE223527"
if (any(idx)) {
  locked[idx, condition := fcase(
    grepl("hfpef", title_lower), "HFpEF",
    grepl("control", title_lower), "Control",
    default = NA_character_
  )]
  locked[idx, assay_modality := fcase(
    grepl("\\bgex\\b", title_lower), "GEX",
    grepl("\\badt\\b", title_lower), "ADT",
    default = NA_character_
  )]
  locked[idx, technology := "PBMC_group_level_matrix"]
  locked[idx, group_id := paste(condition, assay_modality, sep = "__")]
  locked[idx, planned_role := fifelse(
    assay_modality == "GEX",
    "Auxiliary_human_circulating_immune_support",
    "ADT_metadata_only"
  )]
  locked[idx, primary_analysis_unit := "Public_group_level_matrix"]
  locked[idx, include_in_expression_analysis := assay_modality == "GEX"]
  locked[idx, claim_boundary :=
    "Donor-level replication is unavailable unless recoverable from additional metadata; cells are not patients"]
  locked[idx, metadata_rule :=
    "Title-locked condition and modality"]
}

## 5.7 GSE245034: external bulk SGLT2-inhibitor response validation
idx <- locked$dataset_id == "GSE245034"
if (any(idx)) {
  locked[idx, condition := fifelse(
    grepl("^non_failing", title_lower),
    "Control",
    fifelse(grepl("^hfpef", title_lower), "HFpEF", NA_character_)
  )]
  locked[idx, drug := fcase(
    grepl("tya[-_ ]?018", title_lower), "TYA_018",
    grepl("sglt2i|empagliflozin", title_lower), "Empagliflozin",
    grepl("non_failing|\\+veh", title_lower), "Vehicle",
    default = NA_character_
  )]
  locked[idx, sample_subject_id := toupper(
    extract_first_regex(original_title, "(V[0-9]+)")
  )]
  locked[idx, genotype_or_background := "C57BL6_WT"]
  locked[idx, assay_modality := "Bulk_RNA_seq"]
  locked[idx, technology := "Left_ventricle_bulk"]
  locked[idx, group_id := paste(condition, drug, sep = "__")]
  locked[idx, planned_role := "External_sample_level_SGLT2i_validation"]
  locked[idx, primary_analysis_unit := "Biological_sample"]
  locked[idx, claim_boundary :=
    "External sample-level drug-response validation within one study system"]
  locked[idx, metadata_rule :=
    "Title-locked disease group and treatment arm"]
}

## 5.8 GSE249412: cell-type-resolved external drug validation
idx <- locked$dataset_id == "GSE249412"
if (any(idx)) {
  locked[idx, condition := fifelse(
    grepl("^ctrl", title_lower),
    "Control",
    fifelse(grepl("^hfpef", title_lower), "HFpEF", NA_character_)
  )]
  locked[idx, drug := fcase(
    grepl("tya018", title_lower), "TYA_018",
    grepl("empagliflozin", title_lower), "Empagliflozin",
    grepl("vehicle", title_lower), "Vehicle",
    default = NA_character_
  )]
  locked[idx, sample_subject_id := toupper(
    extract_first_regex(original_title, "(V[0-9]+)")
  )]
  locked[idx, genotype_or_background := "C57BL6_WT"]
  locked[idx, assay_modality := "snRNA_seq"]
  locked[idx, technology := "Left_ventricle_single_nucleus"]
  locked[idx, group_id := paste(condition, drug, sep = "__")]
  locked[idx, planned_role := "Cell_type_resolved_SGLT2i_validation"]
  locked[idx, primary_analysis_unit := "Biological_sample_pseudobulk"]
  locked[idx, claim_boundary :=
    "Paired study with GSE245034; orthogonal cellular support, not a second independent cohort"]
  locked[idx, metadata_rule :=
    "Title-locked disease group and treatment arm"]
}

## 5.9 GSE270896: model-background x surgery specificity analysis
idx <- locked$dataset_id == "GSE270896"
if (any(idx)) {
  locked[idx, genotype_or_background := fcase(
    grepl("c57bl/6j|\\bwt_", title_lower), "C57BL6J_WT",
    grepl("b6\\.cg-lep|\\bob_", title_lower), "B6_Cg_Lepob_J",
    default = NA_character_
  )]
  locked[idx, surgery_or_stressor := fcase(
    grepl("\\btac\\b", title_lower), "TAC",
    grepl("\\bsham\\b", title_lower), "Sham",
    default = NA_character_
  )]
  locked[idx, condition := "Model_group_not_HF_classified"]
  locked[idx, diet_or_model := paste(
    genotype_or_background,
    surgery_or_stressor,
    sep = "__"
  )]
  locked[idx, reported_replicate := extract_first_regex(
    original_title,
    "replicate\\s*([0-9]+)"
  )]
  locked[idx, assay_modality := "snRNA_seq"]
  locked[idx, technology := "Whole_heart_single_nucleus"]
  locked[idx, group_id := diet_or_model]
  locked[idx, planned_role := "Model_specificity_comparison"]
  locked[idx, primary_analysis_unit := "Biological_sample_pseudobulk"]
  locked[idx, claim_boundary :=
    "No HFpEF or HFrEF label is assigned from metadata alone; interpret genotype-by-surgery contrasts only"]
  locked[idx, metadata_rule :=
    "Title-locked genetic background and Sham/TAC surgery"]
}

## 5.10 GSE275031: independent mouse HFpEF validation
idx <- locked$dataset_id == "GSE275031"
if (any(idx)) {
  locked[idx, condition := fifelse(
    grepl("^control", title_lower),
    "Control",
    fifelse(grepl("^heartfailure", title_lower), "HFpEF", NA_character_)
  )]
  locked[idx, diet_or_model := fifelse(
    condition == "HFpEF",
    "HFD_LNAME_7wk",
    fifelse(condition == "Control", "Untreated_control", NA_character_)
  )]
  locked[idx, genotype_or_background := "C57BL6N_WT"]
  locked[idx, reported_replicate := extract_first_regex(
    original_title,
    "_([0-9]+)\\s*$"
  )]
  locked[idx, assay_modality := "scRNA_seq"]
  locked[idx, technology := "Ventricular_tissue_single_cell"]
  locked[idx, group_id := condition]
  locked[idx, planned_role := "Independent_external_mouse_HFpEF_validation"]
  locked[idx, primary_analysis_unit := "Biological_sample_pseudobulk"]
  locked[idx, claim_boundary :=
    "Small independent cohort; emphasize effect size and direction consistency"]
  locked[idx, metadata_rule :=
    "Title-locked control versus heartfailure; treatment confirms HFD/L-NAME"]
}

## 5.11 GSE183852: human DCM/non-diseased disease comparator
idx <- locked$dataset_id == "GSE183852"
if (any(idx)) {
  locked[idx, condition := fcase(
    grepl("dilated cardiomyopathy", char_lower), "DCM",
    grepl("non-diseased", char_lower), "Non_diseased",
    default = NA_character_
  )]
  locked[idx, technology := fcase(
    grepl("technology: single cell", char_lower), "Single_Cell",
    grepl("technology: single nuclei", char_lower), "Single_Nuclei",
    default = NA_character_
  )]
  locked[idx, assay_modality := fifelse(
    technology == "Single_Cell",
    "scRNA_seq",
    fifelse(technology == "Single_Nuclei", "snRNA_seq", NA_character_)
  )]
  locked[idx, sample_subject_id := original_title]
  locked[idx, group_id := paste(condition, technology, sep = "__")]
  locked[idx, planned_role := "Human_non_HFpEF_disease_comparator"]
  locked[idx, primary_analysis_unit := "Donor_or_sample_as_reported"]
  locked[idx, claim_boundary :=
    "DCM specificity comparator only; cannot be described as human HFpEF validation"]
  locked[idx, metadata_rule :=
    "Characteristics-locked disease state and technology"]
}

############################################################
## 6. Common fields and deterministic replicate rank
############################################################

## Extract sex where it is explicitly available.
locked[grepl("sex:\\s*male", char_lower), sex := "Male"]
locked[grepl("sex:\\s*female", char_lower), sex := "Female"]

## Create a deterministic within-group rank based on GSM accession.
setorder(locked, dataset_id, group_id, sample_accession)
locked[, replicate_within_group := seq_len(.N), by = .(dataset_id, group_id)]

## Record title-versus-source-file replicate discrepancies without using
## source filename numbering as a biological grouping variable.
locked[, replicate_label_discrepancy := (
  !is.na(reported_replicate) &
  nzchar(reported_replicate) &
  !is.na(source_file_replicate) &
  nzchar(source_file_replicate) &
  reported_replicate != source_file_replicate
)]

locked[, replicate_note := fifelse(
  replicate_label_discrepancy,
  "Title and supplementary filename replicate labels differ; GSM accession and deterministic within-group rank are authoritative.",
  NA_character_
)]

## Required fields for locking.
locked[, required_lock_fields_complete := (
  !is.na(group_id) &
  nzchar(group_id) &
  !is.na(planned_role) &
  nzchar(planned_role) &
  !is.na(primary_analysis_unit) &
  nzchar(primary_analysis_unit)
)]

locked[, lock_status := fcase(
  !required_lock_fields_complete, "FAILED_REQUIRED_FIELD",
  dataset_id == "GSE270896", "LOCKED_MODEL_GROUP_NO_HF_PHENOTYPE",
  dataset_id == "GSE223527", "LOCKED_GROUP_LEVEL_PUBLIC_DATA",
  dataset_id == "GSE183852", "LOCKED_DISEASE_COMPARATOR",
  default = "LOCKED"
)]

############################################################
## 7. Expected group counts and validation
############################################################

expected_groups <- rbindlist(list(
  data.table(
    dataset_id = "GSE237156",
    group_id = as.vector(outer(
      c("Ccr2_positive", "Ccr2_negative"),
      as.vector(outer(
        c("CD", "HFD"),
        c("Vehicle", "Dapagliflozin"),
        paste,
        sep = "__"
      )),
      paste,
      sep = "__"
    )),
    expected_n = 2L
  ),
  data.table(
    dataset_id = "GSE208425",
    group_id = c(
      "WT__CD", "WT__HFD",
      "ApoE_KO__CD", "ApoE_KO__HFD"
    ),
    expected_n = 2L
  ),
  data.table(
    dataset_id = "GSE209548",
    group_id = c(
      "WT__CD", "WT__HFD",
      "ApoE_KO__CD", "ApoE_KO__HFD"
    ),
    expected_n = 3L
  ),
  data.table(
    dataset_id = "GSE236584",
    group_id = c("HFpEF", "Control"),
    expected_n = c(6L, 6L)
  ),
  data.table(
    dataset_id = "GSE236585",
    group_id = c("HFpEF", "Control"),
    expected_n = c(3L, 3L)
  ),
  data.table(
    dataset_id = "GSE223527",
    group_id = c(
      "Control__GEX", "Control__ADT",
      "HFpEF__GEX", "HFpEF__ADT"
    ),
    expected_n = 1L
  ),
  data.table(
    dataset_id = "GSE245034",
    group_id = c(
      "Control__Vehicle",
      "HFpEF__Vehicle",
      "HFpEF__TYA_018",
      "HFpEF__Empagliflozin"
    ),
    expected_n = c(8L, 8L, 11L, 11L)
  ),
  data.table(
    dataset_id = "GSE249412",
    group_id = c(
      "Control__Vehicle",
      "HFpEF__Vehicle",
      "HFpEF__TYA_018",
      "HFpEF__Empagliflozin"
    ),
    expected_n = 2L
  ),
  data.table(
    dataset_id = "GSE270896",
    group_id = c(
      "C57BL6J_WT__Sham",
      "C57BL6J_WT__TAC",
      "B6_Cg_Lepob_J__Sham",
      "B6_Cg_Lepob_J__TAC"
    ),
    expected_n = 2L
  ),
  data.table(
    dataset_id = "GSE275031",
    group_id = c("Control", "HFpEF"),
    expected_n = 2L
  ),
  data.table(
    dataset_id = "GSE183852",
    group_id = c(
      "DCM__Single_Cell",
      "DCM__Single_Nuclei",
      "Non_diseased__Single_Cell",
      "Non_diseased__Single_Nuclei"
    ),
    expected_n = c(5L, 13L, 2L, 25L)
  )
), fill = TRUE)

observed_groups <- locked[, .(
  observed_n = .N,
  sample_accessions = paste(sample_accession, collapse = "; ")
), by = .(dataset_id, group_id)]

group_validation <- merge(
  expected_groups,
  observed_groups,
  by = c("dataset_id", "group_id"),
  all = TRUE
)

group_validation[, validation_status := fcase(
  is.na(expected_n), "UNEXPECTED_GROUP",
  is.na(observed_n), "MISSING_GROUP",
  expected_n == observed_n, "PASS",
  default = "COUNT_MISMATCH"
)]

setorder(group_validation, dataset_id, group_id)

unexpected_datasets <- setdiff(
  unique(locked$dataset_id),
  unique(expected_groups$dataset_id)
)
if (length(unexpected_datasets) > 0L) {
  add_warning(
    "DATASET",
    paste(unexpected_datasets, collapse = ";"),
    "Dataset(s) were present in parsed metadata but had no deterministic locking rule."
  )
}

############################################################
## 8. Dataset-level lock summary
############################################################

dataset_lock_summary <- locked[, .(
  sample_count = .N,
  group_count = uniqueN(group_id),
  unlocked_sample_count = sum(
    !grepl("^LOCKED", lock_status),
    na.rm = TRUE
  ),
  missing_condition = sum(is.na(condition) | !nzchar(condition)),
  missing_group_id = sum(is.na(group_id) | !nzchar(group_id)),
  expression_included_samples = sum(
    include_in_expression_analysis,
    na.rm = TRUE
  ),
  replicate_label_discrepancies = sum(
    replicate_label_discrepancy,
    na.rm = TRUE
  ),
  lock_statuses = paste(sort(unique(lock_status)), collapse = "; "),
  planned_roles = paste(sort(unique(planned_role)), collapse = "; ")
), by = dataset_id]

group_pass_by_dataset <- group_validation[, .(
  group_validation_pass = all(validation_status == "PASS"),
  nonpass_group_count = sum(validation_status != "PASS")
), by = dataset_id]

dataset_lock_summary <- merge(
  dataset_lock_summary,
  group_pass_by_dataset,
  by = "dataset_id",
  all = TRUE
)
dataset_lock_summary[
  is.na(group_validation_pass),
  `:=`(
    group_validation_pass = FALSE,
    nonpass_group_count = NA_integer_
  )
]

dataset_lock_summary[, dataset_lock_status := fcase(
  unlocked_sample_count > 0L, "FAILED_SAMPLE_LOCK",
  !group_validation_pass, "FAILED_GROUP_VALIDATION",
  default = "LOCKED_AND_VALIDATED"
)]

setorder(dataset_lock_summary, dataset_id)

############################################################
## 9. SCP3342 donor-level metadata lock
############################################################

scp_donor_manifest <- data.table()
scp_donor_celltype <- data.table()
scp_validation <- data.table(
  item = c(
    "Total donors",
    "HFpEF donors",
    "Control donors",
    "Total nuclei",
    "Unknown disease-code donors"
  ),
  expected = c(43L, 19L, 24L, 48866L, 0L),
  observed = NA_integer_,
  validation_status = "NOT_RUN"
)

h5ad_files <- list.files(
  DATA_DIR,
  pattern = "\\.h5ad$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)
h5ad_file <- h5ad_files[
  grepl(H5AD_PATTERN, basename(h5ad_files), ignore.case = TRUE)
][1L]

if (
  INSPECT_SCP3342 &&
  length(h5ad_file) == 1L &&
  !is.na(h5ad_file) &&
  file.exists(h5ad_file)
) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    add_warning(
      "SCP3342",
      "hdf5r",
      "hdf5r is unavailable; SCP3342 donor metadata was not locked."
    )
    scp_validation[, validation_status := "HDF5R_UNAVAILABLE"]
  } else {
    log_msg("Locking SCP3342 donor metadata: ", basename(h5ad_file))

    tryCatch({
      h5 <- hdf5r::H5File$new(h5ad_file, mode = "r")
      obs <- h5[["obs"]]
      obs_fields <- obs$ls(recursive = FALSE)$name

      required_obs <- c("donor_id", "disease", "cell_type")
      missing_obs <- setdiff(required_obs, obs_fields)
      if (length(missing_obs) > 0L) {
        stop(
          "SCP3342 obs is missing required field(s): ",
          paste(missing_obs, collapse = ", ")
        )
      }

      donor <- read_anndata_obs_column(obs, "donor_id")
      disease_code <- read_anndata_obs_column(obs, "disease")
      cell_type <- read_anndata_obs_column(obs, "cell_type")
      sex_vec <- if ("sex" %in% obs_fields) {
        read_anndata_obs_column(obs, "sex")
      } else {
        rep(NA_character_, length(donor))
      }

      n_obs <- length(donor)
      if (
        length(disease_code) != n_obs ||
        length(cell_type) != n_obs
      ) {
        stop("SCP3342 donor, disease, and cell-type columns have incompatible lengths.")
      }

      obs_dt <- data.table(
        donor_id = as.character(donor),
        disease_code = as.character(disease_code),
        cell_type = as.character(cell_type),
        sex = as.character(sex_vec)
      )

      obs_dt[, condition := fcase(
        disease_code == "MONDO_0005252", "HFpEF",
        disease_code == "PATO_0000461", "Control",
        default = NA_character_
      )]

      donor_condition_conflicts <- obs_dt[, .(
        n_disease_codes = uniqueN(disease_code),
        n_conditions = uniqueN(condition)
      ), by = donor_id][
        n_disease_codes != 1L |
          n_conditions != 1L
      ]

      if (nrow(donor_condition_conflicts) > 0L) {
        fwrite(
          donor_condition_conflicts,
          file.path(
            DIRS$tables,
            "SCP3342_DONOR_CONDITION_CONFLICTS_FATAL.csv"
          )
        )
        stop(
          "SCP3342 donor-to-condition conflicts were detected. ",
          "See SCP3342_DONOR_CONDITION_CONFLICTS_FATAL.csv."
        )
      }

      scp_donor_manifest <- obs_dt[, .(
        disease_code = unique(disease_code),
        condition = unique(condition),
        sex = paste(sort(unique(sex[!is.na(sex) & nzchar(sex)])), collapse = "; "),
        n_nuclei = .N,
        n_cell_types = uniqueN(cell_type)
      ), by = donor_id]

      scp_donor_manifest[
        sex == "",
        sex := NA_character_
      ]
      scp_donor_manifest[, primary_analysis_unit := "Donor_pseudobulk"]
      scp_donor_manifest[, planned_role :=
        "Primary_independent_human_HFpEF_myocardial_validation"]
      scp_donor_manifest[, claim_boundary :=
        "Donor-level human HFpEF validation; nuclei are not independent replicates"]
      scp_donor_manifest[, lock_status := fifelse(
        is.na(condition),
        "FAILED_UNKNOWN_DISEASE_CODE",
        "LOCKED"
      )]

      setorder(scp_donor_manifest, condition, donor_id)

      scp_donor_celltype <- obs_dt[, .(
        n_nuclei = .N
      ), by = .(
        donor_id,
        disease_code,
        condition,
        cell_type
      )]
      setorder(
        scp_donor_celltype,
        condition,
        donor_id,
        -n_nuclei
      )

      scp_validation[, observed := c(
        uniqueN(scp_donor_manifest$donor_id),
        uniqueN(
          scp_donor_manifest[
            condition == "HFpEF",
            donor_id
          ]
        ),
        uniqueN(
          scp_donor_manifest[
            condition == "Control",
            donor_id
          ]
        ),
        nrow(obs_dt),
        uniqueN(
          scp_donor_manifest[
            is.na(condition),
            donor_id
          ]
        )
      )]
      scp_validation[, validation_status := fifelse(
        expected == observed,
        "PASS",
        "FAIL"
      )]

      try(h5$close_all(), silent = TRUE)
    }, error = function(e) {
      if (exists("h5", inherits = FALSE)) {
        try(h5$close_all(), silent = TRUE)
      }
      add_warning(
        "SCP3342_FATAL",
        basename(h5ad_file),
        conditionMessage(e)
      )
      scp_validation[, validation_status := "FAILED_TO_READ"]
    })
  }
} else {
  add_warning(
    "SCP3342",
    "H5AD",
    "The expected SCP3342 H5AD file was not found."
  )
  scp_validation[, validation_status := "FILE_NOT_FOUND"]
}

############################################################
## 10. Analysis-role and contrast plan
############################################################

contrast_plan <- rbindlist(list(
  data.table(
    dataset_id = "GSE237156",
    analytical_role = "Primary drug-opposed macrophage-program discovery",
    analysis_unit = "Biological sample",
    primary_contrast = "Within Ccr2_positive: HFD_Vehicle vs CD_Vehicle",
    secondary_contrasts = paste(
      "Within Ccr2_positive: HFD_Dapagliflozin vs HFD_Vehicle;",
      "within Ccr2_negative: corresponding diet and drug contrasts;",
      "factorial sensitivity model: macrophage_subset x diet x drug"
    ),
    validation_or_claim_boundary =
      "Do not claim direct dapagliflozin inhibition of a specific TF without experimental evidence"
  ),
  data.table(
    dataset_id = "GSE208425",
    analytical_role = "Internal immune-cell contextualization",
    analysis_unit = "Biological sample pseudobulk after cell annotation",
    primary_contrast = "ApoE_KO_HFD vs ApoE_KO_CD",
    secondary_contrasts =
      "WT_HFD vs WT_CD; genotype-by-diet interaction where estimable",
    validation_or_claim_boundary =
      "Same project family as the drug anchor; not independent replication"
  ),
  data.table(
    dataset_id = "GSE209548",
    analytical_role = "Auxiliary cardiomyocyte response",
    analysis_unit = "Biological sample",
    primary_contrast = "ApoE_KO_HFD vs ApoE_KO_CD",
    secondary_contrasts =
      "WT_HFD vs WT_CD; genotype-by-diet interaction",
    validation_or_claim_boundary =
      "Does not establish macrophage-to-cardiomyocyte causality"
  ),
  data.table(
    dataset_id = "GSE236585",
    analytical_role = "Primary cardiac single-cell discovery",
    analysis_unit = "Biological sample x cell type pseudobulk",
    primary_contrast = "HFpEF vs Control within each major cell type",
    secondary_contrasts =
      "Macrophage-subcluster pseudobulk and sample-level module scoring",
    validation_or_claim_boundary =
      "Cells must not be treated as independent biological replicates"
  ),
  data.table(
    dataset_id = "GSE236584",
    analytical_role = "Matched orthogonal cardiac bulk support",
    analysis_unit = "Biological sample",
    primary_contrast = "HFpEF vs Control",
    secondary_contrasts =
      "Program-level enrichment and direction concordance with GSE236585",
    validation_or_claim_boundary =
      "Same study as GSE236585; not independent validation"
  ),
  data.table(
    dataset_id = "GSE223527",
    analytical_role = "Auxiliary human circulating-immune support",
    analysis_unit = "Public group-level GEX matrix",
    primary_contrast = "HFpEF_GEX vs Control_GEX",
    secondary_contrasts =
      "ADT is excluded from transcriptomic differential expression",
    validation_or_claim_boundary =
      "No patient-level inference unless donor identity can be recovered"
  ),
  data.table(
    dataset_id = "GSE245034",
    analytical_role = "External sample-level SGLT2-inhibitor validation",
    analysis_unit = "Biological sample",
    primary_contrast = "HFpEF_Vehicle vs Control_Vehicle",
    secondary_contrasts = paste(
      "HFpEF_Empagliflozin vs HFpEF_Vehicle;",
      "HFpEF_TYA_018 vs HFpEF_Vehicle"
    ),
    validation_or_claim_boundary =
      "Tests external drug reversal; TYA-018 and empagliflozin are distinct interventions"
  ),
  data.table(
    dataset_id = "GSE249412",
    analytical_role = "Cell-type-resolved external drug validation",
    analysis_unit = "Biological sample x cell type pseudobulk",
    primary_contrast = "HFpEF_Vehicle vs Control_Vehicle",
    secondary_contrasts = paste(
      "HFpEF_Empagliflozin vs HFpEF_Vehicle;",
      "HFpEF_TYA_018 vs HFpEF_Vehicle"
    ),
    validation_or_claim_boundary =
      "Paired study with GSE245034; not a second independent cohort"
  ),
  data.table(
    dataset_id = "GSE270896",
    analytical_role = "Model specificity comparison",
    analysis_unit = "Biological sample x cell type pseudobulk",
    primary_contrast = "Genetic_background x Sham/TAC factorial analysis",
    secondary_contrasts =
      "Within-background TAC vs Sham and interaction testing",
    validation_or_claim_boundary =
      "Do not label groups HFpEF or HFrEF solely from GEO metadata"
  ),
  data.table(
    dataset_id = "GSE275031",
    analytical_role = "Independent external mouse HFpEF validation",
    analysis_unit = "Biological sample x cell type pseudobulk",
    primary_contrast = "HFpEF vs Control",
    secondary_contrasts =
      "Direction concordance, effect size, module score, and communication-pathway replication",
    validation_or_claim_boundary =
      "Small n; prioritize effect size and direction over isolated P values"
  ),
  data.table(
    dataset_id = "SCP3342",
    analytical_role = "Primary independent human myocardial validation",
    analysis_unit = "Donor x cell type pseudobulk",
    primary_contrast = "HFpEF vs Control within macrophage and vascular/stromal compartments",
    secondary_contrasts =
      "Donor-level module scores, TF activity, and candidate-pathway concordance",
    validation_or_claim_boundary =
      "Use 43 donors; never treat 48,866 nuclei as independent patients"
  ),
  data.table(
    dataset_id = "GSE183852",
    analytical_role = "Human non-HFpEF disease comparator",
    analysis_unit = "Donor/sample within technology",
    primary_contrast =
      "DCM vs Non_diseased, stratified by Single_Cell/Single_Nuclei technology",
    secondary_contrasts =
      "Candidate cell-type localization and disease-specificity assessment",
    validation_or_claim_boundary =
      "Cannot be described as human HFpEF validation"
  )
), fill = TRUE)

############################################################
## 11. Final validation status
############################################################

replicate_discrepancies <- locked[
  replicate_label_discrepancy == TRUE,
  .(
    dataset_id,
    sample_accession,
    original_title,
    reported_replicate,
    source_file_replicate,
    group_id,
    replicate_note
  )
]

warnings_dt <- if (length(warning_records) > 0L) {
  rbindlist(warning_records, fill = TRUE)
} else {
  data.table(
    timestamp = character(),
    category = character(),
    item = character(),
    message = character()
  )
}

nonpass_groups <- group_validation[
  validation_status != "PASS"
]
failed_sample_locks <- locked[
  !grepl("^LOCKED", lock_status)
]
failed_scp <- scp_validation[
  validation_status != "PASS"
]

script_copy_methods <- file.path(
  DIRS$methods,
  "HFpEF_Stage1_Metadata_Lock_FIXED_v3.R"
)
script_copy_check <- file.path(
  DIRS$check,
  "HFpEF_Stage1_Metadata_Lock_FIXED_v3.R"
)

script_copy_status <- "UNAVAILABLE"
if (SCRIPT_ARCHIVE_AVAILABLE) {
  copy_methods_ok <- file.copy(
    SCRIPT_FILE,
    script_copy_methods,
    overwrite = TRUE
  )
  copy_check_ok <- file.copy(
    SCRIPT_FILE,
    script_copy_check,
    overwrite = TRUE
  )

  if (isTRUE(copy_methods_ok) && isTRUE(copy_check_ok)) {
    script_copy_status <- "COPIED"
  } else {
    script_copy_status <- "FAILED"
    add_warning(
      "SCRIPT_ARCHIVE",
      basename(SCRIPT_FILE),
      "The exact executed script could not be copied into both output locations."
    )
  }
} else {
  add_warning(
    "SCRIPT_ARCHIVE",
    "running_script",
    paste0(
      "The script path was not available. Run this .R file with source() ",
      "instead of pasting the complete script into the R console."
    )
  )
}

critical_issue_count <- nrow(nonpass_groups) +
  nrow(failed_sample_locks) +
  nrow(failed_scp) +
  as.integer(script_copy_status != "COPIED")

overall_status <- if (critical_issue_count == 0L) {
  "COMPLETED_METADATA_LOCKED"
} else {
  "COMPLETED_REVIEW_REQUIRED"
}

END_TIME <- Sys.time()

run_status <- data.table(
  stage = STAGE_NAME,
  start_time = format(START_TIME, "%Y-%m-%d %H:%M:%S"),
  end_time = format(END_TIME, "%Y-%m-%d %H:%M:%S"),
  elapsed_minutes = round(
    as.numeric(difftime(END_TIME, START_TIME, units = "mins")),
    2
  ),
  geo_sample_count = nrow(locked),
  geo_dataset_count = uniqueN(locked$dataset_id),
  group_validation_nonpass = nrow(nonpass_groups),
  failed_sample_lock_count = nrow(failed_sample_locks),
  scp3342_donor_count = if (nrow(scp_donor_manifest) > 0L) {
    uniqueN(scp_donor_manifest$donor_id)
  } else {
    NA_integer_
  },
  scp3342_validation_nonpass = nrow(failed_scp),
  replicate_label_discrepancy_count = nrow(replicate_discrepancies),
  warning_count = nrow(warnings_dt),
  script_copy_status = script_copy_status,
  critical_issue_count = critical_issue_count,
  overall_status = overall_status
)

############################################################
## 12. Write tables and workbook
############################################################

## Remove internal lowercase helper columns from final manifest.
locked_final <- copy(locked)
locked_final[, c(
  "title_lower",
  "char_lower",
  "supp_lower"
) := NULL]

setcolorder(
  locked_final,
  c(
    "dataset_id",
    "sample_accession",
    "original_title",
    "group_id",
    "condition",
    "genotype_or_background",
    "diet_or_model",
    "drug",
    "macrophage_subset",
    "surgery_or_stressor",
    "assay_modality",
    "technology",
    "sex",
    "reported_replicate",
    "source_file_replicate",
    "replicate_within_group",
    "replicate_label_discrepancy",
    "replicate_note",
    "sample_subject_id",
    "source_name",
    "organism",
    "library_strategy",
    "record_type",
    "planned_role",
    "primary_analysis_unit",
    "include_in_expression_analysis",
    "claim_boundary",
    "metadata_rule",
    "required_lock_fields_complete",
    "lock_status",
    "raw_characteristics",
    "raw_supplementary_files"
  )
)

write_csv_safe(
  locked_final,
  file.path(DIRS$tables, "01_locked_sample_manifest.csv")
)
write_csv_safe(
  group_validation,
  file.path(DIRS$tables, "03_group_count_validation.csv")
)
write_csv_safe(
  scp_donor_manifest,
  file.path(DIRS$tables, "04_SCP3342_locked_donor_manifest.csv")
)
write_csv_safe(
  contrast_plan,
  file.path(DIRS$tables, "05_analysis_role_and_contrast_plan.csv")
)
write_csv_safe(
  scp_donor_celltype,
  file.path(DIRS$tables, "06_SCP3342_donor_celltype_counts.csv")
)
write_csv_safe(
  dataset_lock_summary,
  file.path(DIRS$tables, "07_dataset_lock_summary.csv")
)
write_csv_safe(
  scp_validation,
  file.path(DIRS$tables, "08_SCP3342_validation.csv")
)
write_csv_safe(
  replicate_discrepancies,
  file.path(DIRS$tables, "09_replicate_label_discrepancies.csv")
)
write_csv_safe(
  warnings_dt,
  file.path(DIRS$tables, "10_warnings_and_nonfatal_issues.csv")
)
write_csv_safe(
  run_status,
  file.path(DIRS$tables, "11_run_status.csv")
)

saveRDS(
  locked_final,
  file.path(DIRS$objects, "locked_GEO_sample_manifest.rds")
)
saveRDS(
  scp_donor_manifest,
  file.path(DIRS$objects, "locked_SCP3342_donor_manifest.rds")
)

workbook_path <- file.path(
  DIRS$tables,
  "02_locked_sample_manifest.xlsx"
)
wb <- createWorkbook()
write_sheet_safe(wb, "Locked_GEO_samples", locked_final)
write_sheet_safe(wb, "Group_validation", group_validation)
write_sheet_safe(wb, "Dataset_summary", dataset_lock_summary)
write_sheet_safe(wb, "Contrast_plan", contrast_plan)
write_sheet_safe(wb, "SCP3342_donors", scp_donor_manifest)
write_sheet_safe(wb, "SCP3342_validation", scp_validation)
write_sheet_safe(wb, "Replicate_notes", replicate_discrepancies)
write_sheet_safe(wb, "Warnings", warnings_dt)
write_sheet_safe(wb, "Run_status", run_status)
saveWorkbook(wb, workbook_path, overwrite = TRUE)

############################################################
## 13. Methods and README
############################################################

methods_text <- c(
  "HFpEF Stage 1 Metadata Lock FIXED v3",
  "",
  "Scope:",
  "- Metadata locking only; no expression analysis.",
  "- Dataset-specific rules replaced broad full-text keyword matching.",
  "",
  "Locked GEO datasets:",
  paste0("- ", sort(unique(locked_final$dataset_id))),
  "",
  "Key safeguards:",
  "- Each GEO sample is represented once by dataset ID and GSM accession.",
  "- Expected group sizes are explicitly encoded and validated.",
  "- GSE270896 is represented as genetic background x Sham/TAC without assigning HFpEF/HFrEF labels.",
  "- GSE223527 GEX and ADT are separated; ADT is excluded from transcriptomic expression analysis.",
  "- SCP3342 maps MONDO_0005252 to the HFpEF study cohort and PATO_0000461 to control.",
  "- SCP3342 inference is donor-level; nuclei are not independent patients.",
  "- GSE237156 title labels are used for biological grouping. Any supplementary-filename replicate mismatch is recorded but does not change sample identity.",
  "",
  "SCP3342 ontology mapping:",
  "- MONDO_0005252 -> HFpEF cohort",
  "- PATO_0000461 -> Control",
  "",
  "Pass criterion:",
  "- Every expected GEO group count must pass.",
  "- Every GEO sample must have a locked status.",
  "- SCP3342 must contain 43 donors: 19 HFpEF and 24 controls, totaling 48,866 nuclei.",
  "- The exact executed R script must be archived.",
  "",
  "Interpretation boundary:",
  "- Metadata lock does not establish biological validity or causality."
)
writeLines(
  methods_text,
  file.path(
    DIRS$methods,
    "metadata_lock_methods_and_boundaries.txt"
  ),
  useBytes = TRUE
)

readme_text <- c(
  "HFpEF Reanalysis Project - Stage 1 Metadata Lock FIXED v3",
  "",
  paste0("Overall status: ", overall_status),
  paste0("GEO samples locked: ", nrow(locked_final)),
  paste0("GEO datasets locked: ", uniqueN(locked_final$dataset_id)),
  paste0("Non-passing GEO groups: ", nrow(nonpass_groups)),
  paste0("Failed GEO sample locks: ", nrow(failed_sample_locks)),
  paste0(
    "SCP3342 donors locked: ",
    ifelse(
      nrow(scp_donor_manifest) > 0L,
      uniqueN(scp_donor_manifest$donor_id),
      NA_integer_
    )
  ),
  paste0("SCP3342 non-passing validation items: ", nrow(failed_scp)),
  paste0("Replicate-label discrepancies recorded: ", nrow(replicate_discrepancies)),
  paste0("Exact script archive status: ", script_copy_status),
  "",
  "Main deliverables:",
  "01_locked_sample_manifest.csv",
  "02_locked_sample_manifest.xlsx",
  "03_group_count_validation.csv",
  "04_SCP3342_locked_donor_manifest.csv",
  "05_analysis_role_and_contrast_plan.csv",
  "06_SCP3342_donor_celltype_counts.csv",
  "07_dataset_lock_summary.csv",
  "08_SCP3342_validation.csv",
  "09_replicate_label_discrepancies.csv",
  "10_warnings_and_nonfatal_issues.csv",
  "11_run_status.csv",
  "",
  "Proceed to biological analysis only when overall_status is:",
  "COMPLETED_METADATA_LOCKED"
)
writeLines(
  readme_text,
  file.path(OUT_DIR, "README_stage1_metadata_lock.txt"),
  useBytes = TRUE
)

############################################################
## 14. Build compact review package
############################################################

review_files <- c(
  file.path(DIRS$tables, "01_locked_sample_manifest.csv"),
  workbook_path,
  file.path(DIRS$tables, "03_group_count_validation.csv"),
  file.path(DIRS$tables, "04_SCP3342_locked_donor_manifest.csv"),
  file.path(DIRS$tables, "05_analysis_role_and_contrast_plan.csv"),
  file.path(DIRS$tables, "07_dataset_lock_summary.csv"),
  file.path(DIRS$tables, "08_SCP3342_validation.csv"),
  file.path(DIRS$tables, "09_replicate_label_discrepancies.csv"),
  file.path(DIRS$tables, "10_warnings_and_nonfatal_issues.csv"),
  file.path(DIRS$tables, "11_run_status.csv"),
  file.path(DIRS$methods, "metadata_lock_methods_and_boundaries.txt"),
  file.path(OUT_DIR, "README_stage1_metadata_lock.txt"),
  LOG_FILE,
  script_copy_check
)
review_files <- review_files[file.exists(review_files)]

for (f in review_files) {
  destination <- file.path(DIRS$check, basename(f))
  if (normalizePath(f, winslash = "/", mustWork = FALSE) !=
      normalizePath(destination, winslash = "/", mustWork = FALSE)) {
    file.copy(f, destination, overwrite = TRUE)
  }
}

check_manifest <- data.table(
  filename = basename(list.files(DIRS$check, full.names = TRUE)),
  size_bytes = file.info(
    list.files(DIRS$check, full.names = TRUE)
  )$size
)
fwrite(
  check_manifest,
  file.path(DIRS$check, "CHECK_package_file_manifest.csv")
)

zip::zipr(
  zipfile = CHECK_ZIP,
  files = list.files(DIRS$check, full.names = TRUE),
  root = DIRS$check
)

############################################################
## 15. Final console summary
############################################################

log_msg("Stage 1 metadata lock finished.")
log_msg("Overall status: ", overall_status)
log_msg("GEO samples locked: ", nrow(locked_final))
log_msg("GEO datasets locked: ", uniqueN(locked_final$dataset_id))
log_msg("Non-passing GEO groups: ", nrow(nonpass_groups))
log_msg("Failed GEO sample locks: ", nrow(failed_sample_locks))
log_msg("SCP3342 validation non-pass items: ", nrow(failed_scp))
log_msg("Replicate-label discrepancies: ", nrow(replicate_discrepancies))
log_msg("Script copy status: ", script_copy_status)
log_msg("Output directory: ", OUT_DIR)
log_msg("CHECK package: ", CHECK_ZIP)

cat("\n============================================================\n")
cat("HFpEF Stage 1 Metadata Lock FIXED v3 completed\n")
cat("Status: ", overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat("Proceed only if status is COMPLETED_METADATA_LOCKED.\n")
cat("============================================================\n")
