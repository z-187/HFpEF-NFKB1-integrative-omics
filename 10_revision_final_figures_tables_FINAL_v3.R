############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 8 MULTICOHORT FINAL v6
## Unified, time-efficient multi-cohort validation module
##
## Project root:
##   <HFPEF_PROJECT_DIR>
##
## Core datasets completed in this single stage:
##   1) GSE236584  matched cardiac bulk support
##   2) GSE208425  internal cardiac immune-cell context
##   3) GSE245034  external bulk SGLT2i response validation
##   4) GSE249412  cell-type-resolved SGLT2i validation
##   5) SCP3342    independent human HFpEF myocardial validation
##
## Explicit exclusions from the default Stage 8 run:
##   - GSE275031: deleted/abandoned prior Stage 8; not read here.
##   - GSE223527: group-level public matrices without recoverable
##     patient-level replication; unsuitable as a core validation cohort.
##   - GSE270896: model-specificity resource without a locked HFpEF/HFrEF
##     phenotype assignment; reserved for an optional supplementary analysis.
##
## Design principle:
##   FINAL_v3 is scientifically invalidated and must not be reused because its
##   gene-key function removed lowercase letters before uppercasing. FINAL_v6
##   also corrects FINAL_v5, which removed biologically meaningful punctuation
##   such as hyphens and therefore falsely merged distinct symbols (for example,
##   Mir194-1 with Mir1941 and mt-Tp with Mttp). FINAL_v6 rebuilds every
##   Stage 8 module from raw/frozen inputs under a new schema.
##   Stage 8 is targeted validation of hypotheses frozen in Stages 2-6.
##   It does NOT repeat full discovery clustering, FindAllMarkers, UMAP,
##   scDblFinder, candidate selection, or communication inference.
##
## Time-saving strategy:
##   - bulk cohorts are analyzed from provided count matrices;
##   - GSE208425 and GSE249412 are processed one sample at a time;
##   - fixed marker panels are used for deterministic broad-cell annotation;
##   - only frozen program genes, TF targets, ligands, receptors, and marker
##     genes are retained after sample-level aggregation;
##   - SCP3342 uses its deposited donor and cell-type annotations and is
##     aggregated in backed/chunked Python mode without loading the full
##     3.5-GB AnnData object into R memory;
##   - completed modules are checkpointed and skipped on rerun.
##
## Interpretation boundary:
##   - biological sample/donor is the inferential unit;
##   - cells/nuclei are never treated as independent replicates;
##   - GSE208425 is internal context, not independent replication;
##   - GSE236584 is same-study orthogonal support, not external validation;
##   - GSE245034 and GSE249412 are paired modalities from one study;
##   - SCP3342 is the primary independent human validation cohort;
##   - communication validation is expression support for frozen axes,
##     not proof of physical signaling or causality.
##
## Save this exact file as:
##   <HFPEF_PROJECT_DIR>/
##   HFpEF_Stage8_Multicohort_Validation_FINAL_v6.R
##
## Run from a fresh R session; do not paste line-by-line:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage8_Multicohort_Validation_FINAL_v6.R",
##     encoding = "UTF-8"
##   )
## The script automatically uses the ASCII junction when it exists and otherwise
## falls back to the real project directory. The junction is no longer mandatory.
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(warn = 1)
options(encoding = "UTF-8")
options(timeout = 7200)
options(future.globals.maxSize = 12 * 1024^3)
set.seed(20260714)

############################################################
## 0. Locked paths and run settings
############################################################

DIRECT_PROJECT_DIR <- Sys.getenv("HFPEF_PROJECT_DIR", unset = "")
ASCII_PROJECT_LINK <- Sys.getenv(
  "HFPEF_ASCII_PROJECT_LINK",
  unset = file.path(tempdir(), "HFPEF_STAGE8_ASCII_LINK")
)

## Detect the executing script before any project-path decision. This supports
## source(), Rscript --file, and a saved script pasted into an interactive session.
detect_invoked_script_early <- function() {
  candidates <- character()
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- tryCatch(frames[[i]]$ofile, error = function(e) NULL)
    if (!is.null(ofile) && length(ofile) == 1L && nzchar(ofile)) {
      candidates <- c(candidates, ofile)
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_args <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  candidates <- unique(c(candidates, file_args))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) return(NA_character_)
  gsub("\\\\", "/", path.expand(candidates[1L]))
}

EARLY_SCRIPT_FILE <- detect_invoked_script_early()
EARLY_SCRIPT_DIR <- if (
  length(EARLY_SCRIPT_FILE) == 1L &&
    !is.na(EARLY_SCRIPT_FILE) &&
    nzchar(EARLY_SCRIPT_FILE)
) {
  dirname(EARLY_SCRIPT_FILE)
} else {
  NA_character_
}

## Best-effort creation of an ASCII alias. A directory symlink is attempted
## first; on Windows, a junction is attempted next. Failure is harmless because
## FINAL_v6 can operate directly on the real Unicode project path, while all
## library-sensitive extraction and ZIP work remains under ASCII_TEMP_ROOT.
try_create_ascii_project_alias <- function(link_path, target_path) {
  if (dir.exists(link_path)) return(TRUE)

  suppressWarnings(
    try(file.symlink(target_path, link_path), silent = TRUE)
  )
  if (dir.exists(link_path)) return(TRUE)

  if (.Platform$OS.type == "windows" && nzchar(Sys.which("powershell.exe"))) {
    ps_escape <- function(x) gsub("'", "''", as.character(x), fixed = TRUE)
    ps_command <- paste0(
      "$ErrorActionPreference='Stop'; ",
      "$target='", ps_escape(target_path), "'; ",
      "$link='", ps_escape(link_path), "'; ",
      "if ((Test-Path -LiteralPath $target) -and -not (Test-Path -LiteralPath $link)) {",
      "New-Item -ItemType Junction -Path $link -Target $target | Out-Null",
      "}"
    )
    suppressWarnings(
      try(
        system2(
          "powershell.exe",
          c(
            "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass",
            "-Command", shQuote(ps_command, type = "cmd")
          ),
          stdout = FALSE,
          stderr = FALSE
        ),
        silent = TRUE
      )
    )
  }
  dir.exists(link_path)
}

ASCII_ALIAS_CREATED <- FALSE
if (nzchar(DIRECT_PROJECT_DIR) && dir.exists(DIRECT_PROJECT_DIR)) {
  ASCII_ALIAS_CREATED <- try_create_ascii_project_alias(
    ASCII_PROJECT_LINK,
    DIRECT_PROJECT_DIR
  )
}

project_dir_is_valid <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return(FALSE)
  dir.exists(path) &&
    dir.exists(file.path(path, "0.GEO")) &&
    dir.exists(file.path(path, "01_stage1_metadata_lock_FIXED_v3")) &&
    dir.exists(file.path(path, "02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2"))
}

PROJECT_DIR <- local({
  env_project <- Sys.getenv("HFPEF_PROJECT_DIR", unset = "")
  candidates <- unique(c(
    env_project,
    ASCII_PROJECT_LINK,
    DIRECT_PROJECT_DIR,
    EARLY_SCRIPT_DIR,
    getwd()
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  valid <- vapply(candidates, project_dir_is_valid, logical(1))
  if (!any(valid)) {
    stop(
      "HFpEF project root could not be located. Checked:\n",
      paste(paste0("- ", candidates), collapse = "\n"),
      "\nExpected a directory containing 0.GEO and the completed Stage 1-2 folders."
    )
  }
  gsub("\\\\", "/", path.expand(candidates[which(valid)[1L]]))
})

DATA_DIR <- file.path(PROJECT_DIR, "0.GEO")

STAGE1_DIR <- file.path(PROJECT_DIR, "01_stage1_metadata_lock_FIXED_v3")
STAGE1_MANIFEST_FILE <- file.path(
  STAGE1_DIR, "01_tables", "01_locked_sample_manifest.csv"
)
STAGE1_SCP_DONOR_FILE <- file.path(
  STAGE1_DIR, "01_tables", "04_SCP3342_locked_donor_manifest.csv"
)

STAGE2_DIR <- file.path(
  PROJECT_DIR, "02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2"
)
STAGE2_STATUS_FILE <- file.path(STAGE2_DIR, "01_tables", "24_stage2_run_status.csv")
STAGE2_CHECKS_FILE <- file.path(STAGE2_DIR, "01_tables", "23_scientific_completion_checks.csv")
STAGE2_POS_FILE <- file.path(STAGE2_DIR, "01_tables", "13_opposition_rank_Ccr2_positive.csv.gz")
STAGE2_NEG_FILE <- file.path(STAGE2_DIR, "01_tables", "14_opposition_rank_Ccr2_negative.csv.gz")
STAGE2_CROSS_FILE <- file.path(STAGE2_DIR, "01_tables", "16_cross_subset_consensus_ranking.csv.gz")

STAGE3_DIR <- file.path(
  PROJECT_DIR, "03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH"
)
STAGE3_STATUS_FILE <- file.path(STAGE3_DIR, "01_tables", "39_stage3_run_status.csv")
STAGE3_CHECKS_FILE <- file.path(STAGE3_DIR, "01_tables", "38_scientific_completion_checks.csv")
STAGE3_PROGRAM_STATS_FILE <- file.path(STAGE3_DIR, "01_tables", "27_pseudobulk_program_statistics.csv")

STAGE4_DIR <- file.path(
  PROJECT_DIR, "04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1"
)
STAGE4_STATUS_FILE <- file.path(STAGE4_DIR, "01_tables", "22_stage4_run_status.csv")
STAGE4_CHECKS_FILE <- file.path(STAGE4_DIR, "01_tables", "20_stage4_scientific_completion_checks.csv")
STAGE4_NETWORK_FILE <- file.path(STAGE4_DIR, "01_tables", "05_stage4_full_TF_target_links.csv")
STAGE4_NETWORK_FILE_GZ <- paste0(STAGE4_NETWORK_FILE, ".gz")
STAGE4_ACTIVITY_STATS_FILE <- file.path(
  STAGE4_DIR, "01_tables", "09_stage4_weighted_regulon_activity_HFpEF_vs_Control.csv"
)

STAGE5B_DIR <- file.path(
  PROJECT_DIR, "05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1"
)
STAGE5B_STATUS_FILE <- file.path(STAGE5B_DIR, "01_tables", "17_stage5B_run_status.csv")
STAGE5B_CHECKS_FILE <- file.path(STAGE5B_DIR, "01_tables", "16_stage5B_scientific_completion_checks.csv")
STAGE5B_RANK_FILE <- file.path(STAGE5B_DIR, "01_tables", "13_stage5B_final_candidate_robustness_rank.csv")

STAGE6_DIR <- file.path(
  PROJECT_DIR, "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3"
)
STAGE6_STATUS_FILE <- file.path(STAGE6_DIR, "01_tables", "23_stage6_run_status.csv")
STAGE6_CHECKS_FILE <- file.path(STAGE6_DIR, "01_tables", "22_stage6_scientific_completion_checks.csv")
STAGE6_STABILITY_FILE <- file.path(
  STAGE6_DIR, "01_tables", "18_stage6_axis_ranking_stability_summary.csv"
)

## Core validation inputs.
GSE236584_EXPR_FILE <- file.path(DATA_DIR, "GSE236584_HFpEF_CON_bulk_RNA-seq.txt.gz")
GSE208425_RAW_TAR <- file.path(DATA_DIR, "GSE208425_RAW.tar")
GSE245034_COUNTS_FILE <- file.path(DATA_DIR, "GSE245034_TENAYA0046.counts.txt.gz")
GSE249412_RAW_TAR <- file.path(DATA_DIR, "GSE249412_RAW.tar")
SCP3342_H5AD_FILE <- file.path(DATA_DIR, "HFpEF_snRNAseq_single_cell_portal_10.14.2025.h5ad")

STAGE_NAME <- "08_stage8_multicohort_validation_FINAL_v6"
OUT_DIR <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP <- file.path(PROJECT_DIR, paste0(STAGE_NAME, "_CHECK.zip"))
EXPECTED_SCRIPT_FILE <- file.path(
  PROJECT_DIR,
  "R",
  "08_stage8_multicohort_validation_FINAL_v6.R"
)

## Resume is enabled within FINAL_v6 only. FINAL_v3 checkpoints are never reused.
## The new FINAL_v6 output directory ensures a clean rebuild after the gene-key fix.
ANALYSIS_SCHEMA_VERSION <- "stage8_v6_punctuation_preserving_gene_key_20260715"
RESUME_COMPLETED_MODULES <- TRUE
FORCE_REBUILD_ALL <- FALSE

## Use a stable ASCII-only temporary directory because the project folder contains
## a non-ASCII punctuation character that some Windows/R libraries normalize poorly.
ASCII_TEMP_ROOT <- Sys.getenv(
  "HFPEF_TEMP_DIR",
  unset = file.path(tempdir(), "HFPEF_STAGE8_TEMP")
)

## Frozen hypothesis scope.
SIGNATURE_SIZES <- c(50L, 100L, 150L, 200L)
PRIMARY_SIGNATURE_SIZE <- 150L
CANDIDATE_TF_LIMIT <- 3L
MAX_STABLE_AXES <- 20L

## Targeted validation thresholds.
MIN_FEATURES_PER_CELL <- 200L
MAX_FEATURES_PER_CELL <- 9000L
MAX_COUNTS_PER_CELL <- 120000L
MAX_PERCENT_MT <- 25
MIN_CELLS_PER_SAMPLE_CELLTYPE <- 20L
MIN_RECEIVER_CELLS_PER_SAMPLE <- 30L
MIN_GENE_COVERAGE_FRACTION <- 0.30

## Python-backed SCP3342 settings.
PYTHON_AUTO_INSTALL <- TRUE
PYTHON_CHUNK_CELLS <- 500L

## Reporting.
FORMAL_FDR <- 0.05
EXPLORATORY_FDR <- 0.10
FIGURE_DPI <- 300L

############################################################
## 1. Packages, output directories, and logging
############################################################

ensure_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  missing_after <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_after) > 0L) {
    stop("Required CRAN package(s) unavailable: ", paste(missing_after, collapse = ", "))
  }
}

ensure_bioc <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
  missing_after <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_after) > 0L) {
    stop("Required Bioconductor package(s) unavailable: ", paste(missing_after, collapse = ", "))
  }
}

ensure_cran(c(
  "data.table", "Matrix", "Seurat", "SeuratObject", "ggplot2",
  "openxlsx", "zip", "digest", "scales"
))
ensure_bioc(c("edgeR", "limma"))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

DIRS <- list(
  logs = file.path(OUT_DIR, "00_logs"),
  tables = file.path(OUT_DIR, "01_tables"),
  objects = file.path(OUT_DIR, "02_objects"),
  figures = file.path(OUT_DIR, "03_figures"),
  source = file.path(OUT_DIR, "04_source_data"),
  methods = file.path(OUT_DIR, "05_methods"),
  check = file.path(OUT_DIR, "06_review_check"),
  temp = file.path(ASCII_TEMP_ROOT, STAGE_NAME, "07_temp_ascii")
)

if (!project_dir_is_valid(PROJECT_DIR)) {
  stop("Resolved PROJECT_DIR is not a valid HFpEF project root: ", PROJECT_DIR)
}
if (!dir.create(ASCII_TEMP_ROOT, recursive = TRUE, showWarnings = FALSE) &&
    !dir.exists(ASCII_TEMP_ROOT)) {
  stop("Could not create the required ASCII temporary root: ", ASCII_TEMP_ROOT)
}

if (FORCE_REBUILD_ALL && dir.exists(OUT_DIR)) {
  unlink(OUT_DIR, recursive = TRUE, force = TRUE)
}
if (FORCE_REBUILD_ALL && file.exists(CHECK_ZIP)) {
  unlink(CHECK_ZIP, force = TRUE)
}
for (d in c(OUT_DIR, unlist(DIRS, use.names = FALSE))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(DIRS$logs, "stage8_multicohort_validation.log")
WARN_FILE <- file.path(DIRS$logs, "stage8_warnings.log")
START_TIME <- Sys.time()

log_msg <- function(..., level = "INFO") {
  txt <- paste0(..., collapse = "")
  line <- sprintf(
    "[%s] [%s] %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, txt
  )
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
  invisible(line)
}

warning_records <- list()
add_warning <- function(category, item, message) {
  rec <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    category = as.character(category),
    item = as.character(item),
    message = as.character(message)
  )
  warning_records[[length(warning_records) + 1L]] <<- rec
  cat(
    sprintf("[%s] [%s] %s: %s\n", rec$timestamp, category, item, message),
    file = WARN_FILE, append = TRUE
  )
  log_msg(category, " | ", item, " | ", message, level = "WARN")
  invisible(rec)
}

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(x), path, na = "", compress = "auto")
  invisible(path)
}

safe_fread <- function(path) {
  data.table::fread(path, encoding = "UTF-8", showProgress = FALSE)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) return(NA_character_)
  hit[1L]
}

SCRIPT_FILE <- local({
  candidates <- unique(c(EARLY_SCRIPT_FILE, EXPECTED_SCRIPT_FILE))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) {
    NA_character_
  } else {
    gsub("\\\\", "/", path.expand(candidates[1L]))
  }
})

log_msg("Stage 8 multicohort validation started.")
log_msg("PROJECT_DIR resolved to: ", PROJECT_DIR)
log_msg("ASCII project link available: ", dir.exists(ASCII_PROJECT_LINK))
log_msg("OUT_DIR: ", OUT_DIR)
log_msg("Analysis schema: ", ANALYSIS_SCHEMA_VERSION)
log_msg("Resume completed modules: ", RESUME_COMPLETED_MODULES)
log_msg("Script: ", ifelse(is.na(SCRIPT_FILE), "NOT_DETECTED", SCRIPT_FILE))

############################################################
## 2. General statistical and expression utilities
############################################################

normalize_gene_punctuation <- function(x) {
  y <- as.character(x)
  ## Convert common Unicode dash/minus characters to the ASCII hyphen while
  ## preserving the hyphen itself. Hyphens are biologically meaningful in
  ## symbols such as Mir194-1, mt-Tp, Krtap1-5, and immunoglobulin genes.
  for (cp in c(0x2010, 0x2011, 0x2012, 0x2013, 0x2014, 0x2212)) {
    y <- gsub(intToUtf8(cp), "-", y, fixed = TRUE)
  }
  y
}

strip_gene_version <- function(x) {
  y <- trimws(as.character(x))
  y <- normalize_gene_punctuation(y)

  ## Remove a terminal numeric version only for Ensembl-style stable IDs.
  ## Do not strip ".1" or "-1" from ordinary gene symbols.
  upper_y <- toupper(y)
  is_ensembl_versioned <- grepl(
    "^ENS[A-Z0-9]*[GTP][0-9]+\\.[0-9]+$",
    upper_y
  )
  y[is_ensembl_versioned] <- sub(
    "\\.[0-9]+$",
    "",
    y[is_ensembl_versioned]
  )
  y[is.na(x)] <- NA_character_
  y
}

gene_key <- function(x) {
  y <- strip_gene_version(x)
  ## Matching is case-insensitive, but punctuation is retained. Only whitespace
  ## is removed because it is never a valid part of a gene symbol.
  y <- gsub("[[:space:]]+", "", y)
  y <- toupper(y)
  y[!is.na(y) & !nzchar(y)] <- NA_character_
  y[is.na(x)] <- NA_character_
  y
}

normalized_symbol_label <- function(x) {
  ## This label retains the same biologically meaningful punctuation as gene_key
  ## and is used to distinguish benign case/spacing differences from true
  ## symbol collisions.
  y <- strip_gene_version(x)
  y <- gsub("[[:space:]]+", "", y)
  y <- toupper(y)
  y[!is.na(y) & !nzchar(y)] <- NA_character_
  y[is.na(x)] <- NA_character_
  y
}

safe_filename_token <- function(x) {
  z <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  z <- gsub("^_+|_+$", "", z)
  ifelse(nzchar(z), z, "unnamed")
}

gene_key_audit_file <- file.path(DIRS$tables, "08_gene_key_audit_summary.csv")
existing_gene_audit <- if (file.exists(gene_key_audit_file)) {
  safe_fread(gene_key_audit_file)
} else {
  data.table()
}
gene_key_audit_records <- if (nrow(existing_gene_audit) > 0L) {
  split(existing_gene_audit, as.character(existing_gene_audit$source))
} else {
  list()
}

audit_gene_symbols <- function(symbols, source_label, fatal = TRUE) {
  dt <- data.table(
    source = as.character(source_label),
    original_symbol = as.character(symbols)
  )
  dt <- dt[!is.na(original_symbol) & nzchar(trimws(original_symbol))]
  dt[, normalized_symbol := normalized_symbol_label(original_symbol)]
  dt[, gene_key := gene_key(original_symbol)]
  dt <- dt[nzchar(gene_key)]

  collisions <- dt[, .(
    distinct_normalized_symbols = uniqueN(normalized_symbol),
    normalized_symbols = paste(sort(unique(normalized_symbol)), collapse = ";"),
    example_original_symbols = paste(head(sort(unique(original_symbol)), 20L), collapse = ";")
  ), by = .(source, gene_key)][distinct_normalized_symbols > 1L]

  summary_row <- data.table(
    source = as.character(source_label),
    input_rows = length(symbols),
    nonempty_rows = nrow(dt),
    distinct_normalized_symbols = uniqueN(dt$normalized_symbol),
    distinct_gene_keys = uniqueN(dt$gene_key),
    collision_keys = nrow(collisions),
    status = ifelse(nrow(collisions) == 0L, "PASS", "FAIL")
  )
  gene_key_audit_records[[as.character(source_label)]] <<- summary_row

  dir.create(DIRS$tables, recursive = TRUE, showWarnings = FALSE)
  audit_all <- rbindlist(gene_key_audit_records, use.names = TRUE, fill = TRUE)
  audit_all <- audit_all[, .SD[.N], by = source]
  write_csv_safe(audit_all, gene_key_audit_file)
  if (nrow(collisions) > 0L) {
    collision_file <- file.path(
      DIRS$tables,
      paste0("08_gene_key_collisions_", safe_filename_token(source_label), ".csv")
    )
    write_csv_safe(collisions, collision_file)
    if (fatal) {
      stop(
        "Gene-key collision detected in ", source_label,
        ". Distinct symbols would collapse to the same key. Inspect: ",
        collision_file
      )
    }
  }
  invisible(summary_row)
}

make_axis_key <- function(tf_symbol, ligand, receptor, receiver) {
  paste(
    gene_key(tf_symbol),
    gene_key(ligand),
    gene_key(receptor),
    toupper(gsub("[^A-Za-z0-9]", "", as.character(receiver))),
    sep = "__"
  )
}


## Regression guard for the exact punctuation collisions that invalidated
## FINAL_v5. These pairs must remain distinct under the canonical key.
gene_key_regression_pairs <- data.table(
  symbol_a = c(
    "Mir194-1", "mt-Tp", "Mir692-2", "Mir9-3",
    "Krtap1-5", "Ighv1-41"
  ),
  symbol_b = c(
    "Mir1941", "Mttp", "Mir6922", "Mir93",
    "Krtap15", "Ighv14-1"
  )
)
gene_key_regression_pairs[, key_a := gene_key(symbol_a)]
gene_key_regression_pairs[, key_b := gene_key(symbol_b)]
gene_key_regression_pairs[, distinct := key_a != key_b]
write_csv_safe(
  gene_key_regression_pairs,
  file.path(DIRS$tables, "08A_gene_key_punctuation_regression_guard.csv")
)
if (any(gene_key_regression_pairs$distinct != TRUE)) {
  stop(
    "Gene-key punctuation regression guard failed. FINAL_v5-style false ",
    "symbol merging would recur."
  )
}

hedges_g <- function(case, reference) {
  x <- as.numeric(case); y <- as.numeric(reference)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2L || n2 < 2L) return(NA_real_)
  pooled <- ((n1 - 1) * stats::var(x) + (n2 - 1) * stats::var(y)) /
    (n1 + n2 - 2)
  if (!is.finite(pooled) || pooled <= 0) return(NA_real_)
  d <- (mean(x) - mean(y)) / sqrt(pooled)
  j <- 1 - 3 / (4 * (n1 + n2) - 9)
  j * d
}

cliffs_delta <- function(case, reference) {
  x <- as.numeric(case); y <- as.numeric(reference)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) == 0L || length(y) == 0L) return(NA_real_)
  cmp <- outer(x, y, "-")
  (sum(cmp > 0) - sum(cmp < 0)) / length(cmp)
}

safe_wilcox_p <- function(case, reference) {
  x <- as.numeric(case); y <- as.numeric(reference)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  tryCatch(
    stats::wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

summarize_contrast <- function(
  score_dt,
  item_cols,
  group_col,
  case_group,
  reference_group,
  contrast_name
) {
  dt <- copy(as.data.table(score_dt))
  if (!all(c(item_cols, group_col, "score") %in% names(dt))) {
    stop("summarize_contrast received invalid columns.")
  }
  dt <- dt[get(group_col) %in% c(case_group, reference_group)]
  if (nrow(dt) == 0L) return(data.table())

  out <- dt[, {
    x <- score[get(group_col) == case_group]
    y <- score[get(group_col) == reference_group]
    .(
      case_group = case_group,
      reference_group = reference_group,
      case_n = sum(is.finite(x)),
      reference_n = sum(is.finite(y)),
      case_mean = mean(x, na.rm = TRUE),
      reference_mean = mean(y, na.rm = TRUE),
      effect_case_minus_reference = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      hedges_g = hedges_g(x, y),
      cliffs_delta = cliffs_delta(x, y),
      wilcoxon_p = safe_wilcox_p(x, y)
    )
  }, by = item_cols]
  out[, contrast := contrast_name]
  out[, fdr := p.adjust(wilcoxon_p, method = "BH")]

  if ("gene_coverage" %in% names(dt)) {
    cov <- dt[, .(
      median_gene_coverage = median(gene_coverage, na.rm = TRUE),
      minimum_gene_coverage = min(gene_coverage, na.rm = TRUE)
    ), by = item_cols]
    out <- merge(out, cov, by = item_cols, all.x = TRUE, sort = FALSE)
  }
  if ("target_coverage" %in% names(dt)) {
    cov <- dt[, .(
      median_target_coverage = median(target_coverage, na.rm = TRUE),
      minimum_target_coverage = min(target_coverage, na.rm = TRUE)
    ), by = item_cols]
    out <- merge(out, cov, by = item_cols, all.x = TRUE, sort = FALSE)
  }
  out
}

apply_program_coverage_guard <- function(x) {
  dt <- copy(as.data.table(x))
  dt[, coverage_pass := is.finite(median_gene_coverage) &
       median_gene_coverage >= MIN_GENE_COVERAGE_FRACTION]
  dt[coverage_pass != TRUE, direction_supported := FALSE]
  dt
}

apply_TF_coverage_guard <- function(x) {
  dt <- copy(as.data.table(x))
  dt[, coverage_pass := is.finite(median_target_coverage) &
       median_target_coverage >= MIN_GENE_COVERAGE_FRACTION]
  dt[coverage_pass != TRUE, direction_supported := FALSE]
  dt
}

zscore_within <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  z <- (x - mean(x, na.rm = TRUE)) / s
  z[!is.finite(z)] <- 0
  z
}

aggregate_rows_by_symbol <- function(raw_dt, symbol_col, sample_cols, source_label) {
  dt <- copy(as.data.table(raw_dt))
  if (!symbol_col %in% names(dt)) {
    stop("Missing symbol column in ", source_label, ": ", symbol_col)
  }
  audit_gene_symbols(dt[[symbol_col]], source_label, fatal = TRUE)
  setnames(dt, symbol_col, "gene_symbol")
  dt <- dt[!is.na(gene_symbol) & nzchar(trimws(gene_symbol))]
  dt[, gene_key := gene_key(gene_symbol)]
  dt <- dt[nzchar(gene_key)]
  dt[, lapply(.SD, function(v) sum(as.numeric(v), na.rm = TRUE)),
     by = .(gene_key), .SDcols = sample_cols]
}

counts_to_logcpm_long <- function(count_dt, sample_meta, sample_id_col = "sample_id") {
  sample_cols <- setdiff(names(count_dt), "gene_key")
  sample_meta <- as.data.table(sample_meta)
  if (!sample_id_col %in% names(sample_meta)) {
    stop("Sample metadata lacks identifier column: ", sample_id_col)
  }
  sample_meta <- sample_meta[match(sample_cols, get(sample_id_col))]
  if (any(is.na(sample_meta[[sample_id_col]])) ||
      any(as.character(sample_meta[[sample_id_col]]) != sample_cols)) {
    stop("Sample metadata could not be aligned exactly to count-matrix columns.")
  }
  m <- as.matrix(count_dt[, ..sample_cols])
  storage.mode(m) <- "numeric"
  rownames(m) <- count_dt$gene_key
  y <- edgeR::DGEList(counts = m)
  keep <- edgeR::filterByExpr(y, group = sample_meta$group_id)
  if (sum(keep) >= 50L) y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  dt <- as.data.table(logcpm, keep.rownames = "gene_key")
  dt <- melt(dt, id.vars = "gene_key", variable.name = sample_id_col,
             value.name = "logcpm")
  dt <- merge(dt, sample_meta, by = sample_id_col, all.x = TRUE, sort = FALSE)
  dt
}

candidate_counts_to_logcpm <- function(candidate_counts_long) {
  dt <- copy(as.data.table(candidate_counts_long))
  required <- c("sample_id", "group_id", "cell_type", "gene_key", "count", "total_library")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop("Candidate count table missing: ", paste(missing, collapse = ", "))
  }
  dt[, logcpm := log2((as.numeric(count) / pmax(as.numeric(total_library), 1)) * 1e6 + 1)]
  dt
}

compute_program_scores <- function(expr_long, program_manifest) {
  expr <- copy(as.data.table(expr_long))
  prog <- copy(as.data.table(program_manifest))
  required_expr <- c("sample_id", "group_id", "cell_type", "gene_key", "logcpm")
  if (!all(required_expr %in% names(expr))) {
    stop("Expression table is invalid for program scoring.")
  }
  expr[, gene_z := zscore_within(logcpm), by = .(cell_type, gene_key)]
  joined <- merge(
    prog[, .(program_name, program_size, direction, gene_key)],
    expr,
    by = "gene_key",
    allow.cartesian = TRUE,
    sort = FALSE
  )
  if (nrow(joined) == 0L) return(data.table())
  joined[, signed_z := fifelse(direction == "Disease_up_Drug_down", gene_z, -gene_z)]
  scores <- joined[, .(
    score = mean(signed_z, na.rm = TRUE),
    genes_detected = uniqueN(gene_key),
    requested_genes = uniqueN(prog[program_name == .BY$program_name, gene_key]),
    gene_coverage = uniqueN(gene_key) /
      max(1, uniqueN(prog[program_name == .BY$program_name, gene_key]))
  ), by = .(sample_id, group_id, cell_type, program_name, program_size)]
  scores
}

compute_tf_activity <- function(expr_long, tf_network, candidate_tfs) {
  expr <- copy(as.data.table(expr_long))
  net <- copy(as.data.table(tf_network))
  net <- net[source_symbol %in% candidate_tfs]
  if (nrow(net) == 0L) return(data.table())
  expr[, gene_z := zscore_within(logcpm), by = .(cell_type, gene_key)]
  net[, gene_key := gene_key(target_symbol)]
  net[, signed_weight := as.numeric(mor) * as.numeric(weight)]
  joined <- merge(net, expr, by = "gene_key", allow.cartesian = TRUE, sort = FALSE)
  joined <- joined[is.finite(gene_z) & is.finite(signed_weight)]
  if (nrow(joined) == 0L) return(data.table())
  joined[, contribution := gene_z * signed_weight]
  joined[, .(
    score = sum(contribution, na.rm = TRUE) / sum(abs(signed_weight), na.rm = TRUE),
    targets_detected = uniqueN(gene_key),
    regulon_size_requested = uniqueN(net[source_symbol == .BY$tf_symbol, gene_key]),
    target_coverage = uniqueN(gene_key) /
      max(1, uniqueN(net[source_symbol == .BY$tf_symbol, gene_key]))
  ), by = .(
    sample_id, group_id, cell_type,
    tf_symbol = source_symbol
  )]
}

compute_axis_scores <- function(expr_long, axis_manifest) {
  expr <- copy(as.data.table(expr_long))
  axes <- copy(as.data.table(axis_manifest))
  if (nrow(axes) == 0L) return(data.table())

  expr[, gene_z := zscore_within(logcpm), by = .(cell_type, gene_key)]
  sender <- expr[cell_type == "Macrophage_Monocyte",
                 .(sample_id, group_id, ligand_key = gene_key, ligand_z = gene_z)]
  receptor <- expr[cell_type %in% unique(axes$receiver),
                   .(sample_id, group_id, receiver = cell_type,
                     receptor_key = gene_key, receptor_z = gene_z)]

  axes[, ligand_key := gene_key(nichenet_ligand)]
  axes[, receptor_key := gene_key(receptor)]
  if (!"axis_key" %in% names(axes)) {
    axes[, axis_key := make_axis_key(tf_symbol, nichenet_ligand, receptor, receiver)]
  }
  tmp <- merge(axes, sender, by = "ligand_key", allow.cartesian = TRUE, sort = FALSE)
  tmp <- merge(
    tmp, receptor,
    by = c("sample_id", "group_id", "receiver", "receptor_key"),
    allow.cartesian = TRUE, sort = FALSE
  )
  if (nrow(tmp) == 0L) return(data.table())
  tmp[, .(
    score = mean(c(ligand_z, receptor_z), na.rm = TRUE),
    ligand_z = ligand_z[1L],
    receptor_z = receptor_z[1L]
  ), by = .(
    sample_id, group_id, axis_key, axis_id, tf_symbol,
    nichenet_ligand, receptor, receiver,
    median_scenario_rank, top10_scenario_frequency
  )]
}

############################################################
## 3. Upstream validation and frozen manifests
############################################################

required_inputs <- c(
  PROJECT_DIR, DATA_DIR,
  STAGE1_MANIFEST_FILE, STAGE1_SCP_DONOR_FILE,
  STAGE2_STATUS_FILE, STAGE2_CHECKS_FILE,
  STAGE2_POS_FILE, STAGE2_NEG_FILE, STAGE2_CROSS_FILE,
  STAGE3_STATUS_FILE, STAGE3_CHECKS_FILE, STAGE3_PROGRAM_STATS_FILE,
  STAGE4_STATUS_FILE, STAGE4_CHECKS_FILE, STAGE4_ACTIVITY_STATS_FILE,
  STAGE5B_STATUS_FILE, STAGE5B_CHECKS_FILE, STAGE5B_RANK_FILE,
  STAGE6_STATUS_FILE, STAGE6_CHECKS_FILE, STAGE6_STABILITY_FILE,
  GSE236584_EXPR_FILE, GSE208425_RAW_TAR,
  GSE245034_COUNTS_FILE, GSE249412_RAW_TAR, SCP3342_H5AD_FILE
)
network_input <- first_existing(c(STAGE4_NETWORK_FILE_GZ, STAGE4_NETWORK_FILE))
if (is.na(network_input)) required_inputs <- c(required_inputs, STAGE4_NETWORK_FILE)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required Stage 8 input(s):\n", paste(missing_inputs, collapse = "\n"))
}

validate_upstream <- function(status_file, checks_file, allowed, label) {
  s <- safe_fread(status_file)
  cks <- safe_fread(checks_file)
  if (!"overall_status" %in% names(s)) stop(label, " lacks overall_status.")
  if (!s$overall_status[1L] %in% allowed) {
    stop(label, " is not completed: ", s$overall_status[1L])
  }
  if (!all(c("check", "status") %in% names(cks)) || any(cks$status != "PASS")) {
    stop(label, " contains non-PASS completion checks.")
  }
  data.table(stage = label, overall_status = s$overall_status[1L], checks = nrow(cks))
}

upstream_audit <- rbindlist(list(
  validate_upstream(STAGE2_STATUS_FILE, STAGE2_CHECKS_FILE,
                    "COMPLETED_STAGE2_READY_FOR_REVIEW", "Stage2"),
  validate_upstream(STAGE3_STATUS_FILE, STAGE3_CHECKS_FILE,
                    "COMPLETED_STAGE3_READY_FOR_REVIEW", "Stage3"),
  validate_upstream(STAGE4_STATUS_FILE, STAGE4_CHECKS_FILE,
                    c("COMPLETED_STAGE4_READY_FOR_REVIEW",
                      "COMPLETED_STAGE4_READY_WITH_METHOD_CAUTION"), "Stage4"),
  validate_upstream(STAGE5B_STATUS_FILE, STAGE5B_CHECKS_FILE,
                    c("COMPLETED_STAGE5B_OFFLINE_READY_FOR_REVIEW",
                      "COMPLETED_STAGE5B_READY_FOR_REVIEW"), "Stage5B"),
  validate_upstream(STAGE6_STATUS_FILE, STAGE6_CHECKS_FILE,
                    "COMPLETED_STAGE6_READY_FOR_REVIEW", "Stage6")
))
write_csv_safe(upstream_audit, file.path(DIRS$tables, "01_upstream_status_audit.csv"))

locked_manifest <- safe_fread(STAGE1_MANIFEST_FILE)
core_dataset_ids <- c("GSE236584", "GSE208425", "GSE245034", "GSE249412")
core_meta <- locked_manifest[
  dataset_id %in% core_dataset_ids & include_in_expression_analysis == TRUE
]
if (any(core_meta$lock_status != "LOCKED")) {
  stop("At least one core Stage 8 sample is not metadata-locked.")
}
write_csv_safe(core_meta, file.path(DIRS$tables, "02_core_dataset_locked_metadata.csv"))

## Freeze Stage 2 signed programs.
stage2_pos <- safe_fread(STAGE2_POS_FILE)
stage2_neg <- safe_fread(STAGE2_NEG_FILE)
stage2_cross <- safe_fread(STAGE2_CROSS_FILE)

audit_gene_symbols(stage2_pos$symbol, "Stage2_Ccr2_positive", fatal = TRUE)
audit_gene_symbols(stage2_neg$symbol, "Stage2_Ccr2_negative", fatal = TRUE)
audit_gene_symbols(stage2_cross$symbol, "Stage2_cross_subset", fatal = TRUE)

build_subset_base <- function(dt, source_label) {
  x <- copy(dt)
  required <- c(
    "symbol", "disease_lfc", "drug_lfc", "deseq_opposed", "edger_opposed",
    "four_effect_signs_consistent", "opposition_tier", "within_subset_rank"
  )
  if (!all(required %in% names(x))) {
    stop("Stage 2 subset ranking lacks required columns.")
  }
  x <- x[
    deseq_opposed == TRUE & edger_opposed == TRUE &
      four_effect_signs_consistent == TRUE &
      is.finite(disease_lfc) & is.finite(drug_lfc) &
      !is.na(symbol) & nzchar(symbol)
  ]
  x[, gene_key := gene_key(symbol)]
  setorder(x, within_subset_rank)
  x <- x[, .SD[1L], by = gene_key]
  x[, direction := fcase(
    disease_lfc > 0 & drug_lfc < 0, "Disease_up_Drug_down",
    disease_lfc < 0 & drug_lfc > 0, "Disease_down_Drug_up",
    default = "Other"
  )]
  x <- x[direction != "Other"]
  x[, source := source_label]
  x
}

pos_base <- build_subset_base(stage2_pos, "Ccr2pos")
neg_base <- build_subset_base(stage2_neg, "Ccr2neg")

cross_base <- copy(stage2_cross)
cross_required <- c(
  "symbol", "pos_disease_lfc", "pos_drug_lfc", "neg_disease_lfc",
  "neg_drug_lfc", "consensus_category", "overall_consensus_rank"
)
if (!all(cross_required %in% names(cross_base))) {
  stop("Stage 2 cross-subset ranking lacks required columns.")
}
cross_base <- cross_base[
  consensus_category == "Cross_subset_full_directional_consensus" &
    !is.na(symbol) & nzchar(symbol)
]
cross_base[, gene_key := gene_key(symbol)]
setorder(cross_base, overall_consensus_rank)
cross_base <- cross_base[, .SD[1L], by = gene_key]
cross_base[, disease_lfc := rowMeans(cbind(pos_disease_lfc, neg_disease_lfc), na.rm = TRUE)]
cross_base[, drug_lfc := rowMeans(cbind(pos_drug_lfc, neg_drug_lfc), na.rm = TRUE)]
cross_base[, within_subset_rank := overall_consensus_rank]
cross_base[, direction := fcase(
  pos_disease_lfc > 0 & pos_drug_lfc < 0 &
    neg_disease_lfc > 0 & neg_drug_lfc < 0,
  "Disease_up_Drug_down",
  pos_disease_lfc < 0 & pos_drug_lfc > 0 &
    neg_disease_lfc < 0 & neg_drug_lfc > 0,
  "Disease_down_Drug_up",
  default = "Other"
)]
cross_base <- cross_base[direction != "Other"]
cross_base[, source := "CrossSubset"]

make_program_manifest <- function(base_dt, prefix) {
  records <- list()
  for (n in SIGNATURE_SIZES) {
    up <- base_dt[direction == "Disease_up_Drug_down"][order(within_subset_rank)]
    dn <- base_dt[direction == "Disease_down_Drug_up"][order(within_subset_rank)]
    up <- head(up, n); dn <- head(dn, n)
    program_name <- paste0(prefix, "_Top", n)
    if (nrow(up) > 0L) {
      records[[length(records) + 1L]] <- up[, .(
        program_name, program_size = n, direction, gene_key,
        original_symbol = symbol, stage2_disease_lfc = disease_lfc,
        stage2_drug_lfc = drug_lfc, stage2_rank = within_subset_rank
      )]
    }
    if (nrow(dn) > 0L) {
      records[[length(records) + 1L]] <- dn[, .(
        program_name, program_size = n, direction, gene_key,
        original_symbol = symbol, stage2_disease_lfc = disease_lfc,
        stage2_drug_lfc = drug_lfc, stage2_rank = within_subset_rank
      )]
    }
  }
  rbindlist(records, use.names = TRUE, fill = TRUE)
}

program_manifest <- rbindlist(list(
  make_program_manifest(pos_base, "Ccr2pos"),
  make_program_manifest(neg_base, "Ccr2neg"),
  make_program_manifest(cross_base, "CrossSubset")
), use.names = TRUE, fill = TRUE)
write_csv_safe(program_manifest, file.path(DIRS$tables, "03_frozen_program_manifest.csv"))
program_manifest_size_audit <- program_manifest[, .(
  actual_genes = uniqueN(gene_key),
  duplicate_gene_rows = .N - uniqueN(gene_key)
), by = .(program_name, program_size, direction)]
program_manifest_size_audit[, requested_genes := program_size]
program_manifest_size_audit[, status := fifelse(
  actual_genes > 0L & duplicate_gene_rows == 0L & actual_genes <= requested_genes,
  "PASS", "FAIL"
)]
write_csv_safe(
  program_manifest_size_audit,
  file.path(DIRS$tables, "09_program_manifest_size_audit.csv")
)
if (any(program_manifest_size_audit$status == "FAIL")) {
  stop("Frozen program manifest failed its size/duplicate audit.")
}

## Freeze candidate TFs and TF-target links.
stage5b_rank <- safe_fread(STAGE5B_RANK_FILE)
if (!all(c("tf_symbol", "final_robustness_rank") %in% names(stage5b_rank))) {
  stop("Stage 5B candidate rank file is invalid.")
}
setorder(stage5b_rank, final_robustness_rank)
candidate_tfs <- head(stage5b_rank$tf_symbol, CANDIDATE_TF_LIMIT)

stage4_network <- safe_fread(network_input)
required_net <- c("source_symbol", "target_symbol", "mor", "weight")
if (!all(required_net %in% names(stage4_network))) {
  stop("Stage 4 TF-target network lacks required columns.")
}
audit_gene_symbols(stage4_network$target_symbol, "Stage4_TF_targets", fatal = TRUE)
audit_gene_symbols(stage4_network$source_symbol, "Stage4_TF_sources", fatal = TRUE)
stage4_network <- stage4_network[
  source_symbol %in% candidate_tfs & is.finite(mor) & is.finite(weight)
]
stage4_network[, gene_key := gene_key(target_symbol)]
tf_regulon_size_audit <- stage4_network[, .(
  target_rows = .N,
  unique_target_keys = uniqueN(gene_key),
  duplicate_target_rows = .N - uniqueN(gene_key)
), by = source_symbol]
tf_regulon_size_audit[, status := fifelse(unique_target_keys > 0L, "PASS", "FAIL")]
write_csv_safe(
  tf_regulon_size_audit,
  file.path(DIRS$tables, "09B_TF_regulon_size_audit.csv")
)
if (any(tf_regulon_size_audit$status == "FAIL")) {
  stop("At least one frozen candidate TF has no valid target key.")
}

stage4_activity <- safe_fread(STAGE4_ACTIVITY_STATS_FILE)
activity_effect_col <- intersect(
  c("hfpef_minus_control", "weighted_effect", "effect"), names(stage4_activity)
)[1L]
if (is.na(activity_effect_col)) activity_effect_col <- NA_character_

tf_manifest <- unique(stage5b_rank[tf_symbol %in% candidate_tfs])
if (!is.na(activity_effect_col)) {
  tf_effects <- stage4_activity[, .(
    tf_symbol,
    discovery_HFpEF_minus_Control = get(activity_effect_col)
  )]
  tf_manifest <- merge(tf_manifest, tf_effects, by = "tf_symbol", all.x = TRUE)
}
write_csv_safe(tf_manifest, file.path(DIRS$tables, "04_frozen_TF_manifest.csv"))
write_csv_safe(stage4_network, file.path(DIRS$tables, "05_frozen_TF_target_links.csv.gz"))

## Freeze the most stable Stage 6 axes.
axis_manifest <- safe_fread(STAGE6_STABILITY_FILE)
required_axis <- c(
  "axis_id", "tf_symbol", "nichenet_ligand", "receptor", "receiver",
  "median_scenario_rank", "top10_scenario_frequency"
)
if (!all(required_axis %in% names(axis_manifest))) {
  stop("Stage 6 stability table lacks required axis columns.")
}
axis_manifest <- axis_manifest[
  tf_symbol %in% candidate_tfs &
    receiver %in% c("Endothelial", "Fibroblast", "Pericyte", "Smooth_muscle")
]
audit_gene_symbols(axis_manifest$tf_symbol, "Stage6_axis_TFs", fatal = TRUE)
audit_gene_symbols(axis_manifest$nichenet_ligand, "Stage6_axis_ligands", fatal = TRUE)
audit_gene_symbols(axis_manifest$receptor, "Stage6_axis_receptors", fatal = TRUE)
axis_manifest[, axis_key := make_axis_key(tf_symbol, nichenet_ligand, receptor, receiver)]
setorder(axis_manifest, median_scenario_rank, -top10_scenario_frequency)
axis_manifest <- head(axis_manifest, MAX_STABLE_AXES)
if (anyDuplicated(axis_manifest$axis_key)) {
  dup_axes <- axis_manifest[duplicated(axis_key) | duplicated(axis_key, fromLast = TRUE)]
  write_csv_safe(dup_axes, file.path(DIRS$tables, "08_duplicate_frozen_axis_keys.csv"))
  stop("Duplicate canonical axis_key values were found in the frozen Stage 6 panel.")
}
write_csv_safe(axis_manifest, file.path(DIRS$tables, "06_frozen_axis_manifest.csv"))
axis_key_audit <- axis_manifest[, .(
  axis_id, axis_key, tf_symbol, nichenet_ligand, receptor, receiver,
  key_nonempty = nzchar(axis_key),
  key_matches_definition = axis_key ==
    make_axis_key(tf_symbol, nichenet_ligand, receptor, receiver)
)]
axis_key_audit[, status := fifelse(
  key_nonempty & key_matches_definition, "PASS", "FAIL"
)]
write_csv_safe(axis_key_audit, file.path(DIRS$tables, "09C_axis_key_audit.csv"))
if (any(axis_key_audit$status == "FAIL")) {
  stop("Frozen axis panel failed canonical axis-key validation.")
}

## Candidate gene universe used by all targeted modules.
marker_sets_cardiac <- list(
  Cardiomyocyte = c("Tnnt2", "Tnni3", "Myh6", "Myh7", "Actc1", "Ryr2"),
  Fibroblast = c("Dcn", "Lum", "Col1a1", "Col1a2", "Col3a1", "Pdgfra", "Col14a1"),
  Endothelial = c("Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1", "Emcn", "Ptprb", "Cldn5"),
  Pericyte = c("Pdgfrb", "Rgs5", "Cspg4", "Notch3", "Kcnj8", "Abcc9"),
  Smooth_muscle = c("Acta2", "Tagln", "Myh11", "Cnn1", "Lmod1", "Smtn"),
  Macrophage_Monocyte = c("Lyz2", "Adgre1", "Csf1r", "Cd68", "Fcgr1", "Tyrobp", "C1qa", "C1qb", "C1qc", "Ms4a7"),
  Dendritic_cell = c("Flt3", "Itgax", "Clec10a", "Clec4a1", "H2-Ab1", "Cd74"),
  Neutrophil = c("S100a8", "S100a9", "Mpo", "Elane", "Retnlg", "Ly6g", "Csf3r"),
  T_NK = c("Cd3d", "Cd3e", "Trac", "Nkg7", "Prf1", "Gzmb"),
  B_cell = c("Cd79a", "Cd79b", "Ms4a1", "Cd37", "Cd74", "H2-Ab1"),
  Mast_cell = c("Kit", "Ms4a2", "Cpa3", "Tpsb2", "Mcpt4")
)
marker_sets_immune <- marker_sets_cardiac[c(
  "Macrophage_Monocyte", "Dendritic_cell", "Neutrophil", "T_NK", "B_cell", "Mast_cell"
)]

audit_gene_symbols(
  unlist(marker_sets_cardiac, use.names = FALSE),
  "Fixed_cardiac_marker_panel", fatal = TRUE
)

candidate_mouse_keys <- unique(c(
  program_manifest$gene_key,
  stage4_network$gene_key,
  gene_key(axis_manifest$nichenet_ligand),
  gene_key(axis_manifest$receptor),
  gene_key(unlist(marker_sets_cardiac, use.names = FALSE))
))
candidate_mouse_keys <- candidate_mouse_keys[nzchar(candidate_mouse_keys)]
write_csv_safe(
  data.table(gene_key = candidate_mouse_keys),
  file.path(DIRS$tables, "07_targeted_gene_universe.csv")
)

############################################################
## 4. Fixed marker-score annotation and targeted aggregation
############################################################

get_assay_layer <- function(object, layer = "counts") {
  tryCatch(
    SeuratObject::LayerData(object, assay = "RNA", layer = layer),
    error = function(e) Seurat::GetAssayData(object, assay = "RNA", layer = layer)
  )
}

restore_make_unique_feature_names <- function(features) {
  ## Seurat/ReadMtx may append .1, .2, ... to duplicate row names. Restore the
  ## base symbol only when that base symbol is also present in the same feature
  ## vector. This avoids stripping a legitimate terminal ".1" from a gene name.
  f <- as.character(features)
  base <- sub("\\.[0-9]+$", "", f)
  has_suffix <- grepl("\\.[0-9]+$", f)
  base_is_present <- base %in% f
  restore <- has_suffix & base_is_present
  f[restore] <- base[restore]
  f
}

build_feature_key_map <- function(features, source_label) {
  feature_map <- data.table(
    feature = as.character(features),
    canonical_feature = restore_make_unique_feature_names(features),
    feature_index = seq_along(features)
  )
  audit_gene_symbols(
    feature_map$canonical_feature,
    source_label,
    fatal = TRUE
  )
  feature_map[, gene_key := gene_key(canonical_feature)]
  feature_map <- feature_map[!is.na(gene_key) & nzchar(gene_key)]
  feature_map
}

annotate_by_fixed_markers <- function(
  object, marker_sets, min_raw_score = 0.03, source_label = "single_cell_features"
) {
  data_mat <- get_assay_layer(object, "data")
  feature_map <- build_feature_key_map(rownames(data_mat), source_label)

  score_list <- lapply(names(marker_sets), function(ct) {
    keys <- unique(gene_key(marker_sets[[ct]]))
    keys <- keys[nzchar(keys)]
    key_scores <- lapply(keys, function(k) {
      idx <- feature_map[gene_key == k, feature_index]
      if (length(idx) == 0L) return(NULL)
      Matrix::colSums(data_mat[idx, , drop = FALSE])
    })
    key_scores <- key_scores[!vapply(key_scores, is.null, logical(1))]
    if (length(key_scores) == 0L) {
      rep(0, ncol(data_mat))
    } else {
      Reduce(`+`, key_scores) / length(key_scores)
    }
  })
  score_mat <- do.call(cbind, score_list)
  colnames(score_mat) <- names(marker_sets)
  rownames(score_mat) <- colnames(data_mat)

  score_z <- apply(score_mat, 2L, zscore_within)
  if (is.vector(score_z)) score_z <- matrix(score_z, ncol = 1L)
  colnames(score_z) <- colnames(score_mat)
  rownames(score_z) <- rownames(score_mat)

  best_index <- max.col(score_z, ties.method = "first")
  best_type <- colnames(score_z)[best_index]
  best_raw <- score_mat[cbind(seq_len(nrow(score_mat)), best_index)]
  best_type[!is.finite(best_raw) | best_raw < min_raw_score] <- "Unresolved"

  list(cell_type = best_type, raw_scores = score_mat)
}

aggregate_targeted_from_seurat <- function(
  object, sample_id, group_id, candidate_keys, min_cells
) {
  counts <- get_assay_layer(object, "counts")
  feature_map <- build_feature_key_map(
    rownames(counts), paste0("single_cell_counts_", sample_id)
  )
  selected <- feature_map[gene_key %in% candidate_keys]

  meta <- as.data.table(object@meta.data, keep.rownames = "cell")
  meta <- meta[match(colnames(counts), cell)]
  cell_types <- unique(meta$cell_type)
  records <- list()
  cell_summary <- list()

  for (ct in cell_types) {
    cells <- meta[cell_type == ct, cell]
    cell_summary[[length(cell_summary) + 1L]] <- data.table(
      sample_id = sample_id, group_id = group_id,
      cell_type = ct, n_cells = length(cells)
    )
    threshold <- if (ct %in% c("Endothelial", "Fibroblast", "Pericyte", "Smooth_muscle")) {
      MIN_RECEIVER_CELLS_PER_SAMPLE
    } else {
      min_cells
    }
    if (length(cells) < threshold || nrow(selected) == 0L) next
    idx <- match(cells, colnames(counts))
    total_library <- sum(counts[, idx, drop = FALSE])
    sub <- counts[selected$feature_index, idx, drop = FALSE]
    raw_sums <- Matrix::rowSums(sub)
    summed <- data.table(
      gene_key = selected$gene_key,
      count = as.numeric(raw_sums)
    )[, .(count = sum(count, na.rm = TRUE)), by = gene_key]
    records[[length(records) + 1L]] <- summed[, .(
      sample_id = sample_id,
      group_id = group_id,
      cell_type = ct,
      gene_key,
      count,
      total_library = as.numeric(total_library),
      n_cells = length(cells)
    )]
  }
  list(
    counts = rbindlist(records, use.names = TRUE, fill = TRUE),
    cells = rbindlist(cell_summary, use.names = TRUE, fill = TRUE)
  )
}

process_single_10x_sample <- function(
  matrix_input, sample_id, group_id, marker_sets, candidate_keys,
  input_type = c("mtx", "h5"), features_file = NULL, barcodes_file = NULL
) {
  input_type <- match.arg(input_type)
  if (input_type == "h5") {
    mat <- Seurat::Read10X_h5(matrix_input, use.names = TRUE, unique.features = TRUE)
    if (is.list(mat)) {
      preferred <- intersect(c("Gene Expression", "RNA"), names(mat))
      mat <- if (length(preferred) > 0L) mat[[preferred[1L]]] else mat[[1L]]
    }
  } else {
    feature_table <- safe_fread(features_file)
    feature_column <- if (ncol(feature_table) >= 2L) 2L else 1L
    mat <- Seurat::ReadMtx(
      mtx = matrix_input, features = features_file, cells = barcodes_file,
      feature.column = feature_column
    )
    rownames(mat) <- make.unique(rownames(mat))
  }
  mat <- methods::as(mat, "dgCMatrix")
  obj <- Seurat::CreateSeuratObject(mat, min.cells = 3L, min.features = 100L)
  obj$sample_id <- sample_id
  obj$group_id <- group_id
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-|^Mt-|^MT-")

  keep <- rownames(obj@meta.data)[
    obj$nFeature_RNA >= MIN_FEATURES_PER_CELL &
      obj$nFeature_RNA <= MAX_FEATURES_PER_CELL &
      obj$nCount_RNA <= MAX_COUNTS_PER_CELL &
      obj$percent.mt <= MAX_PERCENT_MT
  ]
  before <- ncol(obj)
  obj <- subset(obj, cells = keep)
  after <- ncol(obj)
  if (after < 50L) stop("Too few cells after QC for sample ", sample_id)

  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  ann <- annotate_by_fixed_markers(
    obj, marker_sets, source_label = paste0("single_cell_normalized_", sample_id)
  )
  obj$cell_type <- ann$cell_type
  agg <- aggregate_targeted_from_seurat(
    obj, sample_id, group_id, candidate_keys, MIN_CELLS_PER_SAMPLE_CELLTYPE
  )
  qc <- data.table(
    sample_id = sample_id, group_id = group_id,
    cells_before_qc = before, cells_after_qc = after,
    retained_fraction = after / before
  )
  rm(obj, mat)
  gc()
  list(counts = agg$counts, cells = agg$cells, qc = qc)
}

module_is_complete <- function(status_path, required_outputs) {
  if (!RESUME_COMPLETED_MODULES || !file.exists(status_path)) return(FALSE)
  if (!all(file.exists(required_outputs))) return(FALSE)
  if (any(file.info(required_outputs)$size <= 0L)) return(FALSE)
  s <- tryCatch(safe_fread(status_path), error = function(e) data.table())
  if (!all(c("status", "analysis_schema") %in% names(s)) || nrow(s) == 0L) return(FALSE)
  identical(s$status[1L], "COMPLETED") &&
    identical(s$analysis_schema[1L], ANALYSIS_SCHEMA_VERSION)
}

write_module_status <- function(path, module, status = "COMPLETED") {
  write_csv_safe(data.table(
    module = module, status = status,
    analysis_schema = ANALYSIS_SCHEMA_VERSION,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ), path)
}

############################################################
## 5. Module A: GSE236584 matched cardiac bulk support
############################################################

module_a_status <- file.path(DIRS$tables, "10A_GSE236584_module_status.csv")
module_a_program_file <- file.path(DIRS$tables, "10B_GSE236584_program_validation.csv")
module_a_tf_file <- file.path(DIRS$tables, "10C_GSE236584_TF_validation.csv")
module_a_expr_file <- file.path(DIRS$source, "10D_GSE236584_targeted_logCPM.csv.gz")

if (!module_is_complete(module_a_status, c(module_a_program_file, module_a_tf_file, module_a_expr_file))) {
  log_msg("Module A: GSE236584 matched bulk support.")
  raw <- safe_fread(GSE236584_EXPR_FILE)
  count_cols <- grep("^Count_", names(raw), value = TRUE)
  if (length(count_cols) != 12L || !"symbol" %in% names(raw)) {
    stop("GSE236584 expression table does not contain the expected 12 count columns.")
  }
  counts_dt <- aggregate_rows_by_symbol(raw, "symbol", count_cols, "GSE236584_bulk_symbols")
  setnames(counts_dt, count_cols, gsub("^Count_\\s*", "", count_cols))
  sample_ids <- setdiff(names(counts_dt), "gene_key")
  meta_a <- data.table(sample_id = sample_ids)
  meta_a[, group_id := fifelse(grepl("^con-", sample_id, ignore.case = TRUE),
                               "Control", "HFpEF")]
  expr_a <- counts_to_logcpm_long(counts_dt, meta_a)
  expr_a[, cell_type := "Whole_heart"]
  write_csv_safe(expr_a[gene_key %in% candidate_mouse_keys], module_a_expr_file)

  prog_scores <- compute_program_scores(
    expr_a[gene_key %in% unique(program_manifest$gene_key)], program_manifest
  )
  tf_scores <- compute_tf_activity(
    expr_a[gene_key %in% unique(stage4_network$gene_key)], stage4_network, candidate_tfs
  )
  prog_res <- summarize_contrast(
    prog_scores, c("cell_type", "program_name", "program_size"),
    "group_id", "HFpEF", "Control", "HFpEF_vs_Control"
  )
  prog_res[, dataset_id := "GSE236584"]
  prog_res[, evidence_role := "Matched_bulk_orthogonal_support"]
  prog_res[, expected_direction := "positive"]
  prog_res[, direction_supported := effect_case_minus_reference > 0]
  prog_res <- apply_program_coverage_guard(prog_res)

  tf_res <- summarize_contrast(
    tf_scores, c("cell_type", "tf_symbol"),
    "group_id", "HFpEF", "Control", "HFpEF_vs_Control"
  )
  tf_res[, dataset_id := "GSE236584"]
  tf_res <- merge(
    tf_res,
    tf_manifest[, .(tf_symbol, discovery_HFpEF_minus_Control)],
    by = "tf_symbol", all.x = TRUE
  )
  tf_res[, direction_supported :=
           sign(effect_case_minus_reference) == sign(discovery_HFpEF_minus_Control)]
  tf_res <- apply_TF_coverage_guard(tf_res)

  write_csv_safe(prog_res, module_a_program_file)
  write_csv_safe(tf_res, module_a_tf_file)
  write_module_status(module_a_status, "GSE236584")
} else {
  log_msg("Module A skipped: completed outputs detected.")
}

############################################################
## 6. Module B: GSE245034 external bulk drug validation
############################################################

module_b_status <- file.path(DIRS$tables, "20A_GSE245034_module_status.csv")
module_b_program_file <- file.path(DIRS$tables, "20B_GSE245034_program_validation.csv")
module_b_tf_file <- file.path(DIRS$tables, "20C_GSE245034_TF_validation.csv")
module_b_expr_file <- file.path(DIRS$source, "20D_GSE245034_targeted_logCPM.csv.gz")

if (!module_is_complete(module_b_status, c(module_b_program_file, module_b_tf_file, module_b_expr_file))) {
  log_msg("Module B: GSE245034 external bulk SGLT2i validation.")
  raw <- safe_fread(GSE245034_COUNTS_FILE)
  if (!"Gene" %in% names(raw)) stop("GSE245034 count table lacks Gene column.")
  count_cols <- setdiff(names(raw), "Gene")
  counts_dt <- aggregate_rows_by_symbol(raw, "Gene", count_cols, "GSE245034_bulk_symbols")

  meta_b <- core_meta[dataset_id == "GSE245034"]
  meta_b[, sample_code := sub(".*_(V[0-9]+)$", "\\1", original_title)]
  meta_b <- meta_b[sample_code %in% count_cols]
  if (nrow(meta_b) != length(count_cols)) {
    missing_codes <- setdiff(count_cols, meta_b$sample_code)
    stop("GSE245034 metadata mapping failed for: ", paste(missing_codes, collapse = ", "))
  }
  sample_meta_b <- meta_b[, .(sample_id = sample_code, group_id)]
  expr_b <- counts_to_logcpm_long(counts_dt, sample_meta_b)
  expr_b[, cell_type := "Whole_heart"]
  write_csv_safe(expr_b[gene_key %in% candidate_mouse_keys], module_b_expr_file)

  prog_scores <- compute_program_scores(
    expr_b[gene_key %in% unique(program_manifest$gene_key)], program_manifest
  )
  tf_scores <- compute_tf_activity(
    expr_b[gene_key %in% unique(stage4_network$gene_key)], stage4_network, candidate_tfs
  )

  contrasts_b <- list(
    disease = c("HFpEF__Vehicle", "Control__Vehicle", "HFpEF_Vehicle_vs_Control_Vehicle"),
    empa = c("HFpEF__Empagliflozin", "HFpEF__Vehicle", "Empagliflozin_vs_HFpEF_Vehicle"),
    tya = c("HFpEF__TYA_018", "HFpEF__Vehicle", "TYA_018_vs_HFpEF_Vehicle")
  )
  prog_list <- lapply(contrasts_b, function(z) {
    summarize_contrast(
      prog_scores, c("cell_type", "program_name", "program_size"),
      "group_id", z[1L], z[2L], z[3L]
    )
  })
  prog_res <- rbindlist(prog_list, fill = TRUE)
  prog_res[, dataset_id := "GSE245034"]
  prog_res[, evidence_role := "External_bulk_SGLT2i_validation"]
  prog_res[, expected_direction := fifelse(
    contrast == "HFpEF_Vehicle_vs_Control_Vehicle", "positive", "negative"
  )]
  prog_res[, direction_supported := fifelse(
    expected_direction == "positive",
    effect_case_minus_reference > 0,
    effect_case_minus_reference < 0
  )]
  prog_res <- apply_program_coverage_guard(prog_res)

  tf_list <- lapply(contrasts_b, function(z) {
    summarize_contrast(
      tf_scores, c("cell_type", "tf_symbol"),
      "group_id", z[1L], z[2L], z[3L]
    )
  })
  tf_res <- rbindlist(tf_list, fill = TRUE)
  tf_res[, dataset_id := "GSE245034"]
  tf_res <- merge(
    tf_res,
    tf_manifest[, .(tf_symbol, discovery_HFpEF_minus_Control)],
    by = "tf_symbol", all.x = TRUE
  )
  tf_res[, direction_supported := fifelse(
    contrast == "HFpEF_Vehicle_vs_Control_Vehicle",
    sign(effect_case_minus_reference) == sign(discovery_HFpEF_minus_Control),
    sign(effect_case_minus_reference) == -sign(discovery_HFpEF_minus_Control)
  )]
  tf_res <- apply_TF_coverage_guard(tf_res)

  write_csv_safe(prog_res, module_b_program_file)
  write_csv_safe(tf_res, module_b_tf_file)
  write_module_status(module_b_status, "GSE245034")
} else {
  log_msg("Module B skipped: completed outputs detected.")
}

############################################################
## 7. Module C: GSE208425 internal immune context
############################################################

module_c_status <- file.path(DIRS$tables, "30A_GSE208425_module_status.csv")
module_c_program_file <- file.path(DIRS$tables, "30B_GSE208425_program_validation.csv")
module_c_tf_file <- file.path(DIRS$tables, "30C_GSE208425_TF_validation.csv")
module_c_counts_file <- file.path(DIRS$source, "30D_GSE208425_targeted_pseudobulk_counts.csv.gz")
module_c_cells_file <- file.path(DIRS$tables, "30E_GSE208425_cell_counts.csv")
module_c_qc_file <- file.path(DIRS$tables, "30F_GSE208425_sample_QC.csv")

if (!module_is_complete(
  module_c_status,
  c(module_c_program_file, module_c_tf_file, module_c_counts_file, module_c_cells_file)
)) {
  log_msg("Module C: GSE208425 internal cardiac immune context.")
  extract_c <- file.path(DIRS$temp, "GSE208425")
  dir.create(extract_c, recursive = TRUE, showWarnings = FALSE)
  members <- utils::untar(GSE208425_RAW_TAR, list = TRUE)
  selected_members <- members[grepl("_(matrix\\.mtx|genes\\.tsv|barcodes\\.tsv)\\.gz$", members)]
  utils::untar(GSE208425_RAW_TAR, files = selected_members, exdir = extract_c)

  meta_c <- core_meta[dataset_id == "GSE208425"]
  count_records <- list(); cell_records <- list(); qc_records <- list()
  for (i in seq_len(nrow(meta_c))) {
    gsm <- meta_c$sample_accession[i]
    prefix_name <- list.files(extract_c, pattern = paste0("^", gsm, ".*_matrix\\.mtx\\.gz$"), full.names = FALSE)
    if (length(prefix_name) != 1L) stop("GSE208425 matrix mapping failed for ", gsm)
    matrix_file <- file.path(extract_c, prefix_name[1L])
    features_file <- sub("_matrix\\.mtx\\.gz$", "_genes.tsv.gz", matrix_file)
    barcodes_file <- sub("_matrix\\.mtx\\.gz$", "_barcodes.tsv.gz", matrix_file)
    log_msg("  GSE208425 sample ", i, "/", nrow(meta_c), ": ", gsm)
    res <- process_single_10x_sample(
      matrix_file, gsm, meta_c$group_id[i], marker_sets_immune,
      candidate_mouse_keys, input_type = "mtx",
      features_file = features_file, barcodes_file = barcodes_file
    )
    count_records[[gsm]] <- res$counts
    cell_records[[gsm]] <- res$cells
    qc_records[[gsm]] <- res$qc
  }
  counts_c <- rbindlist(count_records, fill = TRUE)
  cells_c <- rbindlist(cell_records, fill = TRUE)
  qc_c <- rbindlist(qc_records, fill = TRUE)
  write_csv_safe(counts_c, module_c_counts_file)
  write_csv_safe(cells_c, module_c_cells_file)
  write_csv_safe(qc_c, module_c_qc_file)

  expr_c <- candidate_counts_to_logcpm(counts_c)
  prog_scores <- compute_program_scores(expr_c, program_manifest)
  tf_scores <- compute_tf_activity(expr_c, stage4_network, candidate_tfs)

  contrasts_c <- list(
    ApoE = c("ApoE_KO__HFD", "ApoE_KO__CD", "ApoE_KO_HFD_vs_CD"),
    WT = c("WT__HFD", "WT__CD", "WT_HFD_vs_CD")
  )
  prog_res <- rbindlist(lapply(contrasts_c, function(z) {
    summarize_contrast(
      prog_scores, c("cell_type", "program_name", "program_size"),
      "group_id", z[1L], z[2L], z[3L]
    )
  }), fill = TRUE)
  prog_res[, dataset_id := "GSE208425"]
  prog_res[, evidence_role := "Internal_immune_context"]
  prog_res[, expected_direction := "positive"]
  prog_res[, direction_supported := effect_case_minus_reference > 0]
  prog_res <- apply_program_coverage_guard(prog_res)

  tf_res <- rbindlist(lapply(contrasts_c, function(z) {
    summarize_contrast(
      tf_scores, c("cell_type", "tf_symbol"),
      "group_id", z[1L], z[2L], z[3L]
    )
  }), fill = TRUE)
  tf_res[, dataset_id := "GSE208425"]
  tf_res <- merge(
    tf_res,
    tf_manifest[, .(tf_symbol, discovery_HFpEF_minus_Control)],
    by = "tf_symbol", all.x = TRUE
  )
  tf_res[, direction_supported :=
           sign(effect_case_minus_reference) == sign(discovery_HFpEF_minus_Control)]
  tf_res <- apply_TF_coverage_guard(tf_res)

  write_csv_safe(prog_res, module_c_program_file)
  write_csv_safe(tf_res, module_c_tf_file)
  write_module_status(module_c_status, "GSE208425")
} else {
  log_msg("Module C skipped: completed outputs detected.")
}

############################################################
## 8. Module D: GSE249412 cell-type-resolved drug validation
############################################################

module_d_status <- file.path(DIRS$tables, "40A_GSE249412_module_status.csv")
module_d_program_file <- file.path(DIRS$tables, "40B_GSE249412_program_validation.csv")
module_d_tf_file <- file.path(DIRS$tables, "40C_GSE249412_TF_validation.csv")
module_d_axis_file <- file.path(DIRS$tables, "40D_GSE249412_axis_validation.csv")
module_d_counts_file <- file.path(DIRS$source, "40E_GSE249412_targeted_pseudobulk_counts.csv.gz")
module_d_cells_file <- file.path(DIRS$tables, "40F_GSE249412_cell_counts.csv")
module_d_qc_file <- file.path(DIRS$tables, "40G_GSE249412_sample_QC.csv")

if (!module_is_complete(
  module_d_status,
  c(module_d_program_file, module_d_tf_file, module_d_axis_file,
    module_d_counts_file, module_d_cells_file)
)) {
  log_msg("Module D: GSE249412 cell-type-resolved SGLT2i validation.")
  extract_d <- file.path(DIRS$temp, "GSE249412")
  dir.create(extract_d, recursive = TRUE, showWarnings = FALSE)
  members <- utils::untar(GSE249412_RAW_TAR, list = TRUE)
  h5_members <- members[grepl("filtered_feature_bc_matrix\\.h5$", members)]
  if (length(h5_members) != 8L) stop("Expected 8 GSE249412 filtered H5 matrices.")
  utils::untar(GSE249412_RAW_TAR, files = h5_members, exdir = extract_d)

  meta_d <- core_meta[dataset_id == "GSE249412"]
  count_records <- list(); cell_records <- list(); qc_records <- list()
  for (i in seq_len(nrow(meta_d))) {
    gsm <- meta_d$sample_accession[i]
    h5_name <- list.files(extract_d, pattern = paste0("^", gsm, ".*filtered_feature_bc_matrix\\.h5$"), full.names = FALSE)
    if (length(h5_name) != 1L) stop("GSE249412 H5 mapping failed for ", gsm)
    h5_file <- file.path(extract_d, h5_name[1L])
    log_msg("  GSE249412 sample ", i, "/", nrow(meta_d), ": ", gsm)
    res <- process_single_10x_sample(
      h5_file, gsm, meta_d$group_id[i], marker_sets_cardiac,
      candidate_mouse_keys, input_type = "h5"
    )
    count_records[[gsm]] <- res$counts
    cell_records[[gsm]] <- res$cells
    qc_records[[gsm]] <- res$qc
  }
  counts_d <- rbindlist(count_records, fill = TRUE)
  cells_d <- rbindlist(cell_records, fill = TRUE)
  qc_d <- rbindlist(qc_records, fill = TRUE)
  write_csv_safe(counts_d, module_d_counts_file)
  write_csv_safe(cells_d, module_d_cells_file)
  write_csv_safe(qc_d, module_d_qc_file)

  expr_d <- candidate_counts_to_logcpm(counts_d)
  prog_scores <- compute_program_scores(expr_d, program_manifest)
  tf_scores <- compute_tf_activity(expr_d, stage4_network, candidate_tfs)
  axis_scores <- compute_axis_scores(expr_d, axis_manifest)

  contrasts_d <- list(
    disease = c("HFpEF__Vehicle", "Control__Vehicle", "HFpEF_Vehicle_vs_Control_Vehicle"),
    empa = c("HFpEF__Empagliflozin", "HFpEF__Vehicle", "Empagliflozin_vs_HFpEF_Vehicle"),
    tya = c("HFpEF__TYA_018", "HFpEF__Vehicle", "TYA_018_vs_HFpEF_Vehicle")
  )
  prog_res <- rbindlist(lapply(contrasts_d, function(z) {
    summarize_contrast(
      prog_scores, c("cell_type", "program_name", "program_size"),
      "group_id", z[1L], z[2L], z[3L]
    )
  }), fill = TRUE)
  prog_res[, dataset_id := "GSE249412"]
  prog_res[, evidence_role := "Cell_type_resolved_SGLT2i_validation"]
  prog_res[, expected_direction := fifelse(
    contrast == "HFpEF_Vehicle_vs_Control_Vehicle", "positive", "negative"
  )]
  prog_res[, direction_supported := fifelse(
    expected_direction == "positive", effect_case_minus_reference > 0,
    effect_case_minus_reference < 0
  )]
  prog_res <- apply_program_coverage_guard(prog_res)

  tf_res <- rbindlist(lapply(contrasts_d, function(z) {
    summarize_contrast(
      tf_scores, c("cell_type", "tf_symbol"),
      "group_id", z[1L], z[2L], z[3L]
    )
  }), fill = TRUE)
  tf_res[, dataset_id := "GSE249412"]
  tf_res <- merge(
    tf_res,
    tf_manifest[, .(tf_symbol, discovery_HFpEF_minus_Control)],
    by = "tf_symbol", all.x = TRUE
  )
  tf_res[, direction_supported := fifelse(
    contrast == "HFpEF_Vehicle_vs_Control_Vehicle",
    sign(effect_case_minus_reference) == sign(discovery_HFpEF_minus_Control),
    sign(effect_case_minus_reference) == -sign(discovery_HFpEF_minus_Control)
  )]
  tf_res <- apply_TF_coverage_guard(tf_res)

  axis_res <- rbindlist(lapply(contrasts_d, function(z) {
    summarize_contrast(
      axis_scores,
      c("axis_key", "axis_id", "tf_symbol", "nichenet_ligand", "receptor", "receiver",
        "median_scenario_rank", "top10_scenario_frequency"),
      "group_id", z[1L], z[2L], z[3L]
    )
  }), fill = TRUE)
  axis_res[, dataset_id := "GSE249412"]
  axis_res[, expected_direction := fifelse(
    contrast == "HFpEF_Vehicle_vs_Control_Vehicle", "positive", "negative"
  )]
  axis_res[, direction_supported := fifelse(
    expected_direction == "positive", effect_case_minus_reference > 0,
    effect_case_minus_reference < 0
  )]

  write_csv_safe(prog_res, module_d_program_file)
  write_csv_safe(tf_res, module_d_tf_file)
  write_csv_safe(axis_res, module_d_axis_file)
  write_module_status(module_d_status, "GSE249412")
} else {
  log_msg("Module D skipped: completed outputs detected.")
}

############################################################
## 9. Module E: SCP3342 independent human donor validation
############################################################

module_e_status <- file.path(DIRS$tables, "50A_SCP3342_module_status.csv")
module_e_program_file <- file.path(DIRS$tables, "50B_SCP3342_program_validation.csv")
module_e_tf_file <- file.path(DIRS$tables, "50C_SCP3342_TF_validation.csv")
module_e_axis_file <- file.path(DIRS$tables, "50D_SCP3342_axis_validation.csv")
module_e_counts_file <- file.path(DIRS$source, "50E_SCP3342_targeted_donor_pseudobulk_counts.csv.gz")
module_e_cellmap_file <- file.path(DIRS$tables, "50F_SCP3342_cell_type_mapping.csv")
module_e_genemap_file <- file.path(DIRS$tables, "50G_SCP3342_gene_mapping_audit.csv")

find_python <- function() {
  py <- Sys.which("python")
  if (nzchar(py)) return(list(exe = py, prefix = character()))
  py3 <- Sys.which("python3")
  if (nzchar(py3)) return(list(exe = py3, prefix = character()))
  py_launcher <- Sys.which("py")
  if (nzchar(py_launcher)) return(list(exe = py_launcher, prefix = c("-3.12")))
  stop("Python was not found. Install Python 3.12 or make it available in PATH.")
}

run_python <- function(py, args, stdout = "", stderr = "") {
  system2(py$exe, c(py$prefix, args), stdout = stdout, stderr = stderr)
}

if (!module_is_complete(
  module_e_status,
  c(module_e_program_file, module_e_tf_file, module_e_axis_file,
    module_e_counts_file, module_e_cellmap_file, module_e_genemap_file)
)) {
  log_msg("Module E: SCP3342 independent human donor validation.")

  ## Conserved-symbol projection. This is deliberately transparent: mouse
  ## symbols are converted to uppercase human-style symbols. The retained
  ## mapping and coverage are reported; no unmatched gene is imputed.
  human_key <- gene_key
  candidate_human_symbols <- unique(c(
    human_key(program_manifest$original_symbol),
    human_key(stage4_network$target_symbol),
    human_key(axis_manifest$nichenet_ligand),
    human_key(axis_manifest$receptor)
  ))
  candidate_human_symbols <- candidate_human_symbols[nzchar(candidate_human_symbols) & !is.na(candidate_human_symbols)]
  human_gene_manifest <- data.table(
    conserved_mouse_symbol_key = candidate_human_symbols,
    human_gene_symbol = candidate_human_symbols
  )
  human_gene_file <- file.path(DIRS$temp, "SCP3342_candidate_human_genes.csv")
  write_csv_safe(human_gene_manifest, human_gene_file)

  py_script <- file.path(DIRS$temp, "aggregate_scp3342_targeted.py")
  py_lines <- c(
    "import argparse, os, sys",
    "import numpy as np",
    "import pandas as pd",
    "import scipy.sparse as sp",
    "import anndata as ad",
    "",
    "p=argparse.ArgumentParser()",
    "p.add_argument('--h5ad', required=True)",
    "p.add_argument('--genes', required=True)",
    "p.add_argument('--out_counts', required=True)",
    "p.add_argument('--out_map', required=True)",
    "p.add_argument('--out_gene_map', required=True)",
    "p.add_argument('--chunk', type=int, default=500)",
    "a=p.parse_args()",
    "adata=ad.read_h5ad(a.h5ad, backed='r')",
    "obs=adata.obs.copy()",
    "if 'donor_id' not in obs.columns: raise RuntimeError('Missing obs column: donor_id')",
    "disease_col='original_disease' if 'original_disease' in obs.columns else ('disease__ontology_label' if 'disease__ontology_label' in obs.columns else 'disease')",
    "cell_col='original_cell_type' if 'original_cell_type' in obs.columns else ('cell_type__ontology_label' if 'cell_type__ontology_label' in obs.columns else 'cell_type')",
    "for col in [disease_col, cell_col]:",
    "    if col not in obs.columns: raise RuntimeError(f'Missing obs column: {col}')",
    "obs['donor_id']=obs['donor_id'].astype(str)",
    "obs['disease_raw']=obs[disease_col].astype(str)",
    "obs['cell_type_raw']=obs[cell_col].astype(str)",
    "obs['sex']=obs['sex'].astype(str) if 'sex' in obs.columns else 'Unknown'",
    "def broad_ct(x):",
    "    s=x.lower()",
    "    if 'macroph' in s or 'monocyte' in s or s.startswith('mac'): return 'Macrophage_Monocyte'",
    "    if 'endocardial' in s or ('endothelial' in s and 'lymph' not in s) or s in {'ec','endo','endothelial cell','endothelial1','endothelial2'}: return 'Endothelial'",
    "    if 'fibroblast' in s or s in {'fb','fib'}: return 'Fibroblast'",
    "    if 'pericyte' in s or s in {'pc','peri'}: return 'Pericyte'",
    "    if 'smooth muscle' in s or 'vascular smooth' in s or s in {'smc','vsmc'}: return 'Smooth_muscle'",
    "    return 'Other'",
    "obs['cell_type']=obs['cell_type_raw'].map(broad_ct)",
    "disease_lower=obs['disease_raw'].str.lower()",
    "obs['condition']=np.where(disease_lower.str.contains('hfpef|heart failure|preserved'), 'HFpEF', np.where(disease_lower.str.contains('control|normal'), 'Control', 'Unknown'))",
    "mapping=obs[['cell_type_raw','cell_type']].drop_duplicates().sort_values(['cell_type','cell_type_raw'])",
    "mapping.to_csv(a.out_map, index=False)",
    "keep_ct={'Macrophage_Monocyte','Endothelial','Fibroblast','Pericyte','Smooth_muscle'}",
    "candidate=pd.read_csv(a.genes)",
    "genes=list(dict.fromkeys(candidate['human_gene_symbol'].astype(str).tolist()))",
    "import re",
    "def canonical_gene_key(x):",
    "    s=str(x).strip()",
    "    for ch in ['\\u2010','\\u2011','\\u2012','\\u2013','\\u2014','\\u2212']:",
    "        s=s.replace(ch,'-')",
    "    u=s.upper()",
    "    if re.match(r'^ENS[A-Z0-9]*[GTP][0-9]+\\.[0-9]+$', u):",
    "        s=re.sub(r'\\.[0-9]+$', '', s)",
    "    s=re.sub(r'\\s+', '', s)",
    "    return s.upper()",
    "var_names=[str(x) for x in adata.var_names]",
    "var_keys=[canonical_gene_key(x) for x in var_names]",
    "key_to_names={}",
    "key_to_index={}",
    "for i,(k,nm) in enumerate(zip(var_keys,var_names)):",
    "    if not k: continue",
    "    key_to_names.setdefault(k,set()).add(canonical_gene_key(nm))",
    "    key_to_index.setdefault(k,i)",
    "collision_keys={k:v for k,v in key_to_names.items() if len(v)>1 and k in set(genes)}",
    "if collision_keys: raise RuntimeError('Candidate gene-key collisions in SCP3342 var_names: '+str(collision_keys))",
    "gene_idx=[]; genes_found=[]; map_rows=[]",
    "for g in genes:",
    "    k=canonical_gene_key(g)",
    "    if k in key_to_index:",
    "        gene_idx.append(key_to_index[k]); genes_found.append(k)",
    "        map_rows.append((g,k,var_names[key_to_index[k]],'MATCHED'))",
    "    else:",
    "        map_rows.append((g,k,'','UNMATCHED'))",
    "pd.DataFrame(map_rows, columns=['requested_symbol','gene_key','matched_var_name','status']).to_csv(a.out_gene_map,index=False)",
    "gene_idx=np.asarray(gene_idx,dtype=int)",
    "if len(gene_idx)==0: raise RuntimeError('No candidate genes matched SCP3342 var_names')",
    "layer='raw_count_cellbender' if 'raw_count_cellbender' in adata.layers.keys() else None",
    "mat=adata.layers[layer] if layer is not None else adata.X",
    "keys=(obs['donor_id']+'||'+obs['condition']+'||'+obs['sex']+'||'+obs['cell_type']).to_numpy()",
    "relevant=obs['cell_type'].isin(keep_ct).to_numpy()",
    "acc={}; libs={}; ncells={}",
    "n=adata.n_obs",
    "for start in range(0,n,a.chunk):",
    "    end=min(n,start+a.chunk)",
    "    mask=relevant[start:end]",
    "    if not mask.any(): continue",
    "    block=mat[start:end,:]",
    "    if not sp.issparse(block): block=sp.csr_matrix(block)",
    "    lib=np.asarray(block.sum(axis=1)).ravel()",
    "    sub=block[:,gene_idx]",
    "    local_keys=keys[start:end]",
    "    for k in np.unique(local_keys[mask]):",
    "        rows=np.where((local_keys==k)&mask)[0]",
    "        v=np.asarray(sub[rows,:].sum(axis=0)).ravel()",
    "        acc[k]=acc.get(k, np.zeros(len(genes_found), dtype=np.float64))+v",
    "        libs[k]=libs.get(k,0.0)+float(lib[rows].sum())",
    "        ncells[k]=ncells.get(k,0)+int(len(rows))",
    "rows_out=[]",
    "for k,v in acc.items():",
    "    donor,condition,sex,ct=k.split('||',3)",
    "    for g,c in zip(genes_found,v):",
    "        rows_out.append((donor,condition,sex,ct,g,float(c),float(libs[k]),int(ncells[k])))",
    "out=pd.DataFrame(rows_out, columns=['sample_id','group_id','sex','cell_type','gene_key','count','total_library','n_cells'])",
    "out.to_csv(a.out_counts, index=False, compression='gzip')",
    "adata.file.close()"
  )
  writeLines(py_lines, py_script, useBytes = TRUE)

  py <- find_python()
  check_code <- run_python(
    py,
    c("-c", shQuote("import anndata, scipy, pandas, h5py, numpy")),
    stdout = FALSE, stderr = FALSE
  )
  if (!identical(check_code, 0L) && PYTHON_AUTO_INSTALL) {
    log_msg("Installing required Python packages for backed AnnData aggregation.")
    install_code <- run_python(
      py,
      c("-m", "pip", "install", "--user", "anndata", "scipy", "pandas", "h5py", "numpy")
    )
    if (!identical(install_code, 0L)) {
      stop("Python dependency installation failed for SCP3342.")
    }
  }

  py_log <- file.path(DIRS$logs, "SCP3342_python_aggregation.log")
  exit_code <- run_python(
    py,
    c(
      shQuote(py_script),
      "--h5ad", shQuote(SCP3342_H5AD_FILE),
      "--genes", shQuote(human_gene_file),
      "--out_counts", shQuote(module_e_counts_file),
      "--out_map", shQuote(module_e_cellmap_file),
      "--out_gene_map", shQuote(module_e_genemap_file),
      "--chunk", as.character(PYTHON_CHUNK_CELLS)
    ),
    stdout = py_log, stderr = py_log
  )
  if (!identical(exit_code, 0L) || !file.exists(module_e_counts_file)) {
    stop("SCP3342 Python aggregation failed. Inspect: ", py_log)
  }

  counts_e <- safe_fread(module_e_counts_file)
  counts_e[, gene_key := human_key(gene_key)]
  donor_check <- unique(counts_e[, .(sample_id, group_id)])
  if (uniqueN(donor_check$sample_id) != 43L ||
      sum(donor_check$group_id == "HFpEF") != 19L ||
      sum(donor_check$group_id == "Control") != 24L) {
    stop(
      "SCP3342 donor aggregation must resolve 43 donors (19 HFpEF, 24 Control); observed ",
      uniqueN(donor_check$sample_id), " donors, ",
      sum(donor_check$group_id == "HFpEF"), " HFpEF and ",
      sum(donor_check$group_id == "Control"), " Control."
    )
  }
  expr_e <- candidate_counts_to_logcpm(counts_e)

  ## Human conserved-symbol program and regulon manifests.
  human_program <- copy(program_manifest)
  human_program[, gene_key := human_key(original_symbol)]
  human_network <- copy(stage4_network)
  human_network[, target_symbol := human_key(target_symbol)]
  human_network[, gene_key := target_symbol]
  human_axes <- copy(axis_manifest)
  human_axes[, nichenet_ligand := human_key(nichenet_ligand)]
  human_axes[, receptor := human_key(receptor)]
  human_axes[, axis_key := make_axis_key(tf_symbol, nichenet_ligand, receptor, receiver)]

  prog_scores <- compute_program_scores(expr_e, human_program)
  tf_scores <- compute_tf_activity(expr_e, human_network, candidate_tfs)
  axis_scores <- compute_axis_scores(expr_e, human_axes)

  prog_res <- summarize_contrast(
    prog_scores,
    c("cell_type", "program_name", "program_size"),
    "group_id", "HFpEF", "Control", "HFpEF_vs_Control"
  )
  prog_res[, dataset_id := "SCP3342"]
  prog_res[, evidence_role := "Primary_independent_human_myocardial_validation"]
  prog_res[, expected_direction := "positive"]
  prog_res[, direction_supported := effect_case_minus_reference > 0]
  prog_res <- apply_program_coverage_guard(prog_res)

  tf_res <- summarize_contrast(
    tf_scores, c("cell_type", "tf_symbol"),
    "group_id", "HFpEF", "Control", "HFpEF_vs_Control"
  )
  tf_res[, dataset_id := "SCP3342"]
  tf_res <- merge(
    tf_res,
    tf_manifest[, .(tf_symbol, discovery_HFpEF_minus_Control)],
    by = "tf_symbol", all.x = TRUE
  )
  tf_res[, direction_supported :=
           sign(effect_case_minus_reference) == sign(discovery_HFpEF_minus_Control)]
  tf_res <- apply_TF_coverage_guard(tf_res)

  axis_res <- summarize_contrast(
    axis_scores,
    c("axis_key", "axis_id", "tf_symbol", "nichenet_ligand", "receptor", "receiver",
      "median_scenario_rank", "top10_scenario_frequency"),
    "group_id", "HFpEF", "Control", "HFpEF_vs_Control"
  )
  axis_res[, dataset_id := "SCP3342"]
  axis_res[, expected_direction := "positive"]
  axis_res[, direction_supported := effect_case_minus_reference > 0]

  write_csv_safe(prog_res, module_e_program_file)
  write_csv_safe(tf_res, module_e_tf_file)
  write_csv_safe(axis_res, module_e_axis_file)
  write_module_status(module_e_status, "SCP3342")
} else {
  log_msg("Module E skipped: completed outputs detected.")
}

############################################################
## 10. Unified multicohort evidence integration
############################################################

log_msg("Integrating multicohort validation evidence.")
program_files <- c(
  module_a_program_file, module_b_program_file,
  module_c_program_file, module_d_program_file, module_e_program_file
)
tf_files <- c(
  module_a_tf_file, module_b_tf_file,
  module_c_tf_file, module_d_tf_file, module_e_tf_file
)
axis_files <- c(module_d_axis_file, module_e_axis_file)

program_evidence <- rbindlist(lapply(program_files, safe_fread), fill = TRUE)
tf_evidence <- rbindlist(lapply(tf_files, safe_fread), fill = TRUE)
axis_evidence <- rbindlist(lapply(axis_files, safe_fread), fill = TRUE)

## Prioritize the predeclared Top150 programs for integrated reporting;
## Top50/100/200 remain as sensitivity results in the full tables.
program_primary <- program_evidence[program_size == PRIMARY_SIGNATURE_SIZE]

program_summary <- program_primary[, .(
  evidence_rows = .N,
  datasets = uniqueN(dataset_id),
  compartments = uniqueN(cell_type),
  supported_rows = sum(direction_supported %in% TRUE, na.rm = TRUE),
  support_fraction = mean(direction_supported %in% TRUE, na.rm = TRUE),
  median_abs_hedges_g = median(abs(hedges_g), na.rm = TRUE),
  formal_fdr_rows = sum(fdr < FORMAL_FDR, na.rm = TRUE),
  exploratory_fdr_rows = sum(fdr < EXPLORATORY_FDR, na.rm = TRUE)
), by = program_name]
setorder(program_summary, -support_fraction, -median_abs_hedges_g)
program_summary[, integrated_rank := seq_len(.N)]

tf_summary <- tf_evidence[, .(
  evidence_rows = .N,
  datasets = uniqueN(dataset_id),
  compartments = uniqueN(cell_type),
  supported_rows = sum(direction_supported %in% TRUE, na.rm = TRUE),
  support_fraction = mean(direction_supported %in% TRUE, na.rm = TRUE),
  median_abs_hedges_g = median(abs(hedges_g), na.rm = TRUE),
  formal_fdr_rows = sum(fdr < FORMAL_FDR, na.rm = TRUE)
), by = tf_symbol]
setorder(tf_summary, -support_fraction, -median_abs_hedges_g)
tf_summary[, integrated_rank := seq_len(.N)]

axis_evidence[, axis_key := make_axis_key(tf_symbol, nichenet_ligand, receptor, receiver)]
axis_summary <- axis_evidence[, .(
  axis_id = axis_key[1L],
  tf_symbol = gene_key(tf_symbol[1L]),
  nichenet_ligand = gene_key(nichenet_ligand[1L]),
  receptor = gene_key(receptor[1L]),
  receiver = toupper(receiver[1L]),
  evidence_rows = .N,
  datasets = uniqueN(dataset_id),
  dataset_list = paste(sort(unique(dataset_id)), collapse = ";"),
  supported_rows = sum(direction_supported %in% TRUE, na.rm = TRUE),
  support_fraction = mean(direction_supported %in% TRUE, na.rm = TRUE),
  median_abs_hedges_g = median(abs(hedges_g), na.rm = TRUE),
  formal_fdr_rows = sum(fdr < FORMAL_FDR, na.rm = TRUE)
), by = axis_key]
setorder(axis_summary, -support_fraction, -datasets, -median_abs_hedges_g)
axis_summary[, integrated_rank := seq_len(.N)]

write_csv_safe(program_evidence, file.path(DIRS$tables, "60_multicohort_program_evidence.csv.gz"))
write_csv_safe(tf_evidence, file.path(DIRS$tables, "61_multicohort_TF_evidence.csv.gz"))
write_csv_safe(axis_evidence, file.path(DIRS$tables, "62_multicohort_axis_evidence.csv.gz"))
write_csv_safe(program_summary, file.path(DIRS$tables, "63_program_integrated_summary.csv"))
write_csv_safe(tf_summary, file.path(DIRS$tables, "64_TF_integrated_summary.csv"))
write_csv_safe(axis_summary, file.path(DIRS$tables, "65_axis_integrated_summary.csv"))

claim_boundary <- data.table(
  dataset_id = c("GSE236584", "GSE208425", "GSE245034", "GSE249412", "SCP3342"),
  role = c(
    "Matched orthogonal cardiac bulk support",
    "Internal immune-cell context",
    "External sample-level SGLT2i response validation",
    "Cell-type-resolved SGLT2i validation paired with GSE245034",
    "Primary independent human HFpEF myocardial validation"
  ),
  independent_external_validation = c(FALSE, FALSE, TRUE, FALSE, TRUE),
  primary_analysis_unit = c(
    "biological sample", "biological sample x immune cell type",
    "biological sample", "biological sample x cardiac cell type",
    "donor x deposited cell type"
  ),
  claim_boundary = c(
    "Same study as GSE236585; not external replication",
    "Same project family as GSE237156; not independent replication",
    "Cross-SGLT2 inhibitor validation; empagliflozin is not dapagliflozin",
    "Paired modality with GSE245034; not a second independent cohort",
    "Human donor-level validation; conserved-symbol projection is not proof of direct orthology for every gene"
  )
)
write_csv_safe(claim_boundary, file.path(DIRS$tables, "66_dataset_roles_and_claim_boundaries.csv"))

############################################################
## 11. Figures and key-results workbook
############################################################

plot_program <- program_primary[
  cell_type %in% c("Whole_heart", "Macrophage_Monocyte") &
    contrast %in% c(
      "HFpEF_vs_Control", "HFpEF_Vehicle_vs_Control_Vehicle",
      "Empagliflozin_vs_HFpEF_Vehicle", "ApoE_KO_HFD_vs_CD"
    )
]
plot_program[, column_label := paste(dataset_id, cell_type, contrast, sep = " | ")]
plot_program[, display_effect := pmax(pmin(hedges_g, 4), -4)]

if (nrow(plot_program) > 0L) {
  p1 <- ggplot(plot_program, aes(x = column_label, y = program_name, fill = display_effect)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_point(
      data = plot_program[direction_supported == TRUE],
      aes(x = column_label, y = program_name),
      inherit.aes = FALSE, shape = 21, size = 1.8, fill = "white"
    ) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B33A3A",
                         midpoint = 0, limits = c(-4, 4), oob = scales::squish) +
    labs(
      title = "Frozen drug-opposed programs across validation cohorts",
      x = NULL, y = NULL, fill = "Hedges' g\n(clipped)",
      caption = "White circles denote prespecified directional support."
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1))
  ggsave(file.path(DIRS$figures, "Fig8A_multicohort_program_validation.png"),
         p1, width = 13, height = 6.5, dpi = FIGURE_DPI)
  ggsave(file.path(DIRS$figures, "Fig8A_multicohort_program_validation.pdf"),
         p1, width = 13, height = 6.5)
}

plot_tf <- tf_evidence[
  cell_type %in% c("Whole_heart", "Macrophage_Monocyte") &
    contrast %in% c(
      "HFpEF_vs_Control", "HFpEF_Vehicle_vs_Control_Vehicle",
      "Empagliflozin_vs_HFpEF_Vehicle", "ApoE_KO_HFD_vs_CD"
    )
]
plot_tf[, column_label := paste(dataset_id, cell_type, contrast, sep = " | ")]
plot_tf[, display_effect := pmax(pmin(hedges_g, 4), -4)]
if (nrow(plot_tf) > 0L) {
  p2 <- ggplot(plot_tf, aes(x = column_label, y = tf_symbol, fill = display_effect)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_point(
      data = plot_tf[direction_supported == TRUE],
      aes(x = column_label, y = tf_symbol),
      inherit.aes = FALSE, shape = 21, size = 2, fill = "white"
    ) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B33A3A",
                         midpoint = 0, limits = c(-4, 4), oob = scales::squish) +
    labs(title = "Frozen TF activities across validation cohorts",
         x = NULL, y = NULL, fill = "Hedges' g\n(clipped)") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1))
  ggsave(file.path(DIRS$figures, "Fig8B_multicohort_TF_validation.png"),
         p2, width = 13, height = 4.8, dpi = FIGURE_DPI)
  ggsave(file.path(DIRS$figures, "Fig8B_multicohort_TF_validation.pdf"),
         p2, width = 13, height = 4.8)
}

plot_axes <- axis_evidence[
  contrast %in% c("HFpEF_vs_Control", "HFpEF_Vehicle_vs_Control_Vehicle",
                  "Empagliflozin_vs_HFpEF_Vehicle")
]
plot_axes[, column_label := paste(dataset_id, contrast, sep = " | ")]
plot_axes[, axis_short := paste(
  gene_key(tf_symbol), gene_key(nichenet_ligand),
  gene_key(receptor), toupper(receiver), sep = "-"
)]
plot_axes[, display_effect := pmax(pmin(hedges_g, 4), -4)]
plot_axes <- plot_axes[axis_key %in% head(axis_summary$axis_key, 15L)]
if (nrow(plot_axes) > 0L) {
  p3 <- ggplot(plot_axes, aes(x = column_label, y = axis_short, fill = display_effect)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_point(
      data = plot_axes[direction_supported == TRUE],
      aes(x = column_label, y = axis_short),
      inherit.aes = FALSE, shape = 21, size = 1.8, fill = "white"
    ) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B33A3A",
                         midpoint = 0, limits = c(-4, 4), oob = scales::squish) +
    labs(title = "Frozen macrophage-to-vascular/stromal axes",
         x = NULL, y = NULL, fill = "Hedges' g\n(clipped)") +
    theme_bw(base_size = 8.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(DIRS$figures, "Fig8C_multicohort_axis_validation.png"),
         p3, width = 11, height = 7.5, dpi = FIGURE_DPI)
  ggsave(file.path(DIRS$figures, "Fig8C_multicohort_axis_validation.pdf"),
         p3, width = 11, height = 7.5)
}

workbook_file <- file.path(DIRS$tables, "70_stage8_multicohort_key_results.xlsx")
wb <- openxlsx::createWorkbook()
add_sheet <- function(name, x) {
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, as.data.frame(x), withFilter = TRUE)
  openxlsx::freezePane(wb, name, firstRow = TRUE)
  openxlsx::setColWidths(wb, name, cols = 1:min(ncol(x), 25), widths = "auto")
}
add_sheet("dataset_roles", claim_boundary)
add_sheet("program_summary", program_summary)
add_sheet("TF_summary", tf_summary)
add_sheet("axis_summary", axis_summary)
add_sheet("program_evidence", program_evidence)
add_sheet("TF_evidence", tf_evidence)
add_sheet("axis_evidence", axis_evidence)
add_sheet("frozen_programs", program_manifest)
add_sheet("frozen_TFs", tf_manifest)
add_sheet("frozen_axes", axis_manifest)
add_sheet("gene_key_audit", safe_fread(gene_key_audit_file))
add_sheet("program_size_audit", program_manifest_size_audit)
add_sheet("TF_regulon_audit", tf_regulon_size_audit)
add_sheet("axis_key_audit", axis_key_audit)
add_sheet("SCP_gene_mapping", safe_fread(module_e_genemap_file))
openxlsx::saveWorkbook(wb, workbook_file, overwrite = TRUE)

############################################################
## 12. Methods, checks, status, and CHECK package
############################################################

methods_text <- c(
  "Stage 8 was a targeted multicohort validation module. All programs, TF candidates,",
  "TF-target links, and ligand-receptor axes were frozen from Stages 2-6 before any",
  "Stage 8 outcome was inspected. GSE236584 provided same-study cardiac bulk support;",
  "GSE208425 provided same-project immune-cell contextualization; GSE245034 provided",
  "external sample-level empagliflozin and TYA-018 response validation; GSE249412",
  "provided cell-type-resolved validation from the paired study; and SCP3342 provided",
  "independent human donor-level myocardial validation.",
  "",
  "For time-efficient validation, GSE208425 and GSE249412 were processed independently",
  "by biological sample using filtered 10x matrices, fixed QC thresholds, log normalization,",
  "and deterministic broad-cell marker scores. No discovery clustering, UMAP, marker",
  "selection, or outcome-driven annotation was performed. Candidate-gene counts were",
  "aggregated by sample and cell type, with total library size retained for targeted logCPM.",
  "SCP3342 was read in backed, chunked Python mode using deposited donor, disease, sex,",
  "and cell-type annotations. Raw CellBender counts were preferred when available.",
  "",
  "Signed program scores were calculated from within-compartment gene-standardized",
  "expression, with disease-up/drug-down genes contributing positively and disease-down/",
  "drug-up genes contributing negatively. TF activity used the locked Stage 4 signed,",
  "weighted regulon formula. Communication support used the mean standardized expression",
  "of the frozen macrophage ligand and receiver-cell receptor within each biological sample.",
  "All gene identifiers were canonicalized only after preserving both uppercase and lowercase",
  "letters; version suffixes were removed, symbols were converted to uppercase, and fatal",
  "collision audits prevented distinct symbols from being merged. Ligand-receptor axes used",
  "a canonical cross-species axis_key for mouse-human integration.",
  "All inferential summaries used biological samples or human donors; cells were never",
  "treated as independent replicates. Direction, effect size, and cross-cohort consistency",
  "were primary; P values and false-discovery rates were secondary."
)
writeLines(methods_text, file.path(DIRS$methods, "stage8_methods_and_claim_boundaries.txt"))

parameters <- data.table(
  parameter = c(
    "SIGNATURE_SIZES", "PRIMARY_SIGNATURE_SIZE", "CANDIDATE_TF_LIMIT",
    "MAX_STABLE_AXES", "MIN_FEATURES_PER_CELL", "MAX_FEATURES_PER_CELL",
    "MAX_COUNTS_PER_CELL", "MAX_PERCENT_MT", "MIN_CELLS_PER_SAMPLE_CELLTYPE",
    "MIN_RECEIVER_CELLS_PER_SAMPLE", "PYTHON_CHUNK_CELLS",
    "ANALYSIS_SCHEMA_VERSION", "PROJECT_DIR", "DIRECT_PROJECT_DIR",
    "ASCII_PROJECT_LINK", "ASCII_ALIAS_CREATED", "ASCII_TEMP_ROOT",
    "RESUME_COMPLETED_MODULES", "FORCE_REBUILD_ALL"
  ),
  value = c(
    paste(SIGNATURE_SIZES, collapse = ";"), PRIMARY_SIGNATURE_SIZE,
    CANDIDATE_TF_LIMIT, MAX_STABLE_AXES, MIN_FEATURES_PER_CELL,
    MAX_FEATURES_PER_CELL, MAX_COUNTS_PER_CELL, MAX_PERCENT_MT,
    MIN_CELLS_PER_SAMPLE_CELLTYPE, MIN_RECEIVER_CELLS_PER_SAMPLE,
    PYTHON_CHUNK_CELLS, ANALYSIS_SCHEMA_VERSION, PROJECT_DIR,
    DIRECT_PROJECT_DIR, ASCII_PROJECT_LINK, ASCII_ALIAS_CREATED, ASCII_TEMP_ROOT,
    RESUME_COMPLETED_MODULES, FORCE_REBUILD_ALL
  )
)
write_csv_safe(parameters, file.path(DIRS$methods, "stage8_parameters.csv"))

if (!is.na(SCRIPT_FILE) && file.exists(SCRIPT_FILE)) {
  file.copy(
    SCRIPT_FILE,
    file.path(DIRS$methods, basename(SCRIPT_FILE)),
    overwrite = TRUE
  )
} else {
  add_warning(
    "REPRODUCIBILITY", "script_archive",
    "The running script path was not detected. Save the exact script under the expected project path and rerun with source()."
  )
}

warnings_dt <- if (length(warning_records) > 0L) {
  rbindlist(warning_records, fill = TRUE)
} else {
  data.table(
    timestamp = character(), category = character(),
    item = character(), message = character()
  )
}
write_csv_safe(warnings_dt, file.path(DIRS$tables, "71_warnings_and_nonfatal_issues.csv"))

## Re-read persisted module outputs so completion checks also work after checkpoint resume.
module_status_files <- c(
  module_a_status, module_b_status, module_c_status, module_d_status, module_e_status
)
module_status_dt <- rbindlist(lapply(module_status_files, safe_fread), fill = TRUE)
counts_c_check <- safe_fread(module_c_counts_file)
counts_d_check <- safe_fread(module_d_counts_file)
counts_e_check <- safe_fread(module_e_counts_file)
cells_c_check <- safe_fread(module_c_cells_file)
cells_d_check <- safe_fread(module_d_cells_file)
scp_gene_map_check <- safe_fread(module_e_genemap_file)

gene_audit_check <- if (file.exists(gene_key_audit_file)) {
  safe_fread(gene_key_audit_file)
} else {
  rbindlist(gene_key_audit_records, use.names = TRUE, fill = TRUE)
}
gene_audit_check <- gene_audit_check[, .SD[.N], by = source]

expected_meta_counts <- data.table(
  dataset_id = c(
    rep("GSE236584", 2L), rep("GSE208425", 4L),
    rep("GSE245034", 4L), rep("GSE249412", 4L)
  ),
  group_id = c(
    "Control", "HFpEF",
    "ApoE_KO__CD", "ApoE_KO__HFD", "WT__CD", "WT__HFD",
    "Control__Vehicle", "HFpEF__Vehicle", "HFpEF__Empagliflozin", "HFpEF__TYA_018",
    "Control__Vehicle", "HFpEF__Vehicle", "HFpEF__Empagliflozin", "HFpEF__TYA_018"
  ),
  expected_n = c(6L, 6L, 2L, 2L, 2L, 2L, 8L, 8L, 11L, 11L, 2L, 2L, 2L, 2L)
)
observed_meta_counts <- core_meta[, .(observed_n = uniqueN(sample_accession)), by = .(dataset_id, group_id)]
meta_count_audit <- merge(
  expected_meta_counts, observed_meta_counts,
  by = c("dataset_id", "group_id"), all.x = TRUE
)
meta_count_audit[is.na(observed_n), observed_n := 0L]
meta_count_audit[, status := fifelse(observed_n == expected_n, "PASS", "FAIL")]
write_csv_safe(meta_count_audit, file.path(DIRS$tables, "74_biological_sample_count_audit.csv"))

scp_donor_audit <- unique(counts_e_check[, .(sample_id, group_id)])[, .(
  donors = uniqueN(sample_id)
), by = group_id]
write_csv_safe(scp_donor_audit, file.path(DIRS$tables, "75_SCP3342_donor_count_audit.csv"))

pseudobulk_duplicate_rows <- sum(duplicated(counts_c_check[, .(sample_id, group_id, cell_type, gene_key)])) +
  sum(duplicated(counts_d_check[, .(sample_id, group_id, cell_type, gene_key)])) +
  sum(duplicated(counts_e_check[, .(sample_id, group_id, cell_type, gene_key)]))

axis_key_consistent <- all(
  axis_evidence$axis_key ==
    make_axis_key(
      axis_evidence$tf_symbol,
      axis_evidence$nichenet_ligand,
      axis_evidence$receptor,
      axis_evidence$receiver
    )
)

checks <- data.table(
  check = c(
    "all_five_core_modules_completed_under_v6_schema",
    "GSE236584_results_present",
    "GSE208425_results_present",
    "GSE245034_results_present",
    "GSE249412_results_present",
    "SCP3342_results_present",
    "all_locked_sample_counts_match_design",
    "SCP3342_has_19_HFpEF_and_24_Control_donors",
    "gene_key_collision_audits_all_pass",
    "gene_key_punctuation_regression_guard_pass",
    "program_manifest_size_audit_pass",
    "program_manifest_has_no_duplicate_keys",
    "TF_regulons_have_targets",
    "frozen_axis_keys_unique_and_nonempty",
    "canonical_axis_keys_consistent_in_evidence",
    "pseudobulk_tables_have_no_duplicate_sample_celltype_gene_keys",
    "single_cell_modules_use_expected_biological_samples",
    "SCP3342_candidate_gene_mapping_nonempty",
    "program_support_requires_minimum_gene_coverage",
    "TF_support_requires_minimum_target_coverage",
    "program_integration_nonempty",
    "TF_integration_nonempty",
    "axis_integration_nonempty",
    "GSE275031_not_used",
    "running_script_archived",
    "workbook_present"
  ),
  pass = c(
    nrow(module_status_dt) == 5L &&
      all(module_status_dt$status == "COMPLETED") &&
      all(module_status_dt$analysis_schema == ANALYSIS_SCHEMA_VERSION),
    file.exists(module_a_program_file) && file.info(module_a_program_file)$size > 0L,
    file.exists(module_c_program_file) && file.info(module_c_program_file)$size > 0L,
    file.exists(module_b_program_file) && file.info(module_b_program_file)$size > 0L,
    file.exists(module_d_program_file) && file.info(module_d_program_file)$size > 0L,
    file.exists(module_e_program_file) && file.info(module_e_program_file)$size > 0L,
    all(meta_count_audit$status == "PASS"),
    scp_donor_audit[group_id == "HFpEF", donors] == 19L &&
      scp_donor_audit[group_id == "Control", donors] == 24L,
    nrow(gene_audit_check) > 0L && all(gene_audit_check$status == "PASS") &&
      sum(gene_audit_check$collision_keys) == 0L,
    nrow(gene_key_regression_pairs) > 0L &&
      all(gene_key_regression_pairs$distinct == TRUE),
    all(program_manifest_size_audit$status == "PASS"),
    !anyDuplicated(program_manifest[, .(program_name, direction, gene_key)]),
    all(tf_regulon_size_audit$status == "PASS"),
    nrow(axis_manifest) > 0L && !anyDuplicated(axis_manifest$axis_key) &&
      all(nzchar(axis_manifest$axis_key)),
    isTRUE(axis_key_consistent),
    pseudobulk_duplicate_rows == 0L,
    uniqueN(cells_c_check$sample_id) == 8L && uniqueN(cells_d_check$sample_id) == 8L,
    nrow(scp_gene_map_check[status == "MATCHED"]) > 0L,
    all(program_evidence[direction_supported == TRUE,
      median_gene_coverage >= MIN_GENE_COVERAGE_FRACTION]),
    all(tf_evidence[direction_supported == TRUE,
      median_target_coverage >= MIN_GENE_COVERAGE_FRACTION]),
    nrow(program_evidence) > 0L,
    nrow(tf_evidence) > 0L,
    nrow(axis_evidence) > 0L,
    !any(grepl("GSE275031", c(required_inputs, program_evidence$dataset_id,
                              tf_evidence$dataset_id, axis_evidence$dataset_id))),
    !is.na(SCRIPT_FILE) && file.exists(file.path(DIRS$methods, basename(SCRIPT_FILE))),
    file.exists(workbook_file) && file.info(workbook_file)$size > 0L
  )
)
checks[, status := fifelse(pass, "PASS", "FAIL")]
write_csv_safe(checks, file.path(DIRS$tables, "72_stage8_scientific_completion_checks.csv"))

if (any(checks$status == "FAIL")) {
  stop(
    "Stage 8 finished calculations but failed completion checks: ",
    paste(checks[status == "FAIL", check], collapse = ", ")
  )
}

run_status <- data.table(
  overall_status = "COMPLETED_STAGE8_MULTICOHORT_READY_FOR_REVIEW",
  analysis_schema = ANALYSIS_SCHEMA_VERSION,
  started_at = format(START_TIME, "%Y-%m-%d %H:%M:%S"),
  completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  elapsed_minutes = round(as.numeric(difftime(Sys.time(), START_TIME, units = "mins")), 2),
  core_datasets = "GSE236584;GSE208425;GSE245034;GSE249412;SCP3342",
  abandoned_dataset = "GSE275031_NOT_USED",
  primary_human_validation = "SCP3342",
  program_rows = nrow(program_evidence),
  TF_rows = nrow(tf_evidence),
  axis_rows = nrow(axis_evidence),
  cross_dataset_axis_keys = sum(axis_summary$datasets >= 2L),
  gene_key_collision_keys = sum(gene_audit_check$collision_keys),
  warnings = nrow(warnings_dt)
)
write_csv_safe(run_status, file.path(DIRS$tables, "73_stage8_run_status.csv"))

readme <- c(
  "HFpEF Stage 8 multicohort validation FINAL v6",
  "",
  "Core datasets analyzed:",
  "- GSE236584: matched cardiac bulk support",
  "- GSE208425: internal cardiac immune context",
  "- GSE245034: external bulk SGLT2i response",
  "- GSE249412: cell-type-resolved SGLT2i response",
  "- SCP3342: independent human donor-level myocardial validation",
  "",
  "GSE275031 is not used in this stage.",
  "GSE223527 and GSE270896 are intentionally excluded from the core module because",
  "their metadata boundaries do not support the same level of patient-level HFpEF validation.",
  "",
  paste0("Overall status: ", run_status$overall_status),
  paste0("Elapsed minutes in this run: ", run_status$elapsed_minutes)
)
writeLines(readme, file.path(OUT_DIR, "README_stage8.txt"))

## Compact review package.
if (dir.exists(DIRS$check)) unlink(DIRS$check, recursive = TRUE, force = TRUE)
dir.create(DIRS$check, recursive = TRUE, showWarnings = FALSE)

check_files <- c(
  file.path(DIRS$tables, c(
    "01_upstream_status_audit.csv",
    "02_core_dataset_locked_metadata.csv",
    "03_frozen_program_manifest.csv",
    "04_frozen_TF_manifest.csv",
    "06_frozen_axis_manifest.csv",
    "08_gene_key_audit_summary.csv",
    "08A_gene_key_punctuation_regression_guard.csv",
    "09_program_manifest_size_audit.csv",
    "09B_TF_regulon_size_audit.csv",
    "09C_axis_key_audit.csv",
    "10B_GSE236584_program_validation.csv",
    "10C_GSE236584_TF_validation.csv",
    "20B_GSE245034_program_validation.csv",
    "20C_GSE245034_TF_validation.csv",
    "30B_GSE208425_program_validation.csv",
    "30C_GSE208425_TF_validation.csv",
    "30E_GSE208425_cell_counts.csv",
    "30F_GSE208425_sample_QC.csv",
    "40B_GSE249412_program_validation.csv",
    "40C_GSE249412_TF_validation.csv",
    "40D_GSE249412_axis_validation.csv",
    "40F_GSE249412_cell_counts.csv",
    "40G_GSE249412_sample_QC.csv",
    "50B_SCP3342_program_validation.csv",
    "50C_SCP3342_TF_validation.csv",
    "50D_SCP3342_axis_validation.csv",
    "50F_SCP3342_cell_type_mapping.csv",
    "50G_SCP3342_gene_mapping_audit.csv",
    "63_program_integrated_summary.csv",
    "64_TF_integrated_summary.csv",
    "65_axis_integrated_summary.csv",
    "66_dataset_roles_and_claim_boundaries.csv",
    "70_stage8_multicohort_key_results.xlsx",
    "71_warnings_and_nonfatal_issues.csv",
    "72_stage8_scientific_completion_checks.csv",
    "73_stage8_run_status.csv",
    "74_biological_sample_count_audit.csv",
    "75_SCP3342_donor_count_audit.csv"
  )),
  list.files(DIRS$figures, full.names = TRUE),
  file.path(DIRS$methods, c(
    "stage8_methods_and_claim_boundaries.txt",
    "stage8_parameters.csv"
  )),
  file.path(OUT_DIR, "README_stage8.txt")
)
if (!is.na(SCRIPT_FILE)) {
  check_files <- c(check_files, file.path(DIRS$methods, basename(SCRIPT_FILE)))
}
check_files <- unique(check_files[file.exists(check_files)])
file.copy(check_files, DIRS$check, overwrite = TRUE)

manifest_check <- data.table(
  filename = basename(check_files),
  source_path = gsub("\\\\", "/", path.expand(check_files)),
  size_bytes = file.info(check_files)$size,
  md5 = unname(tools::md5sum(check_files))
)
write_csv_safe(manifest_check, file.path(DIRS$check, "CHECK_package_file_manifest.csv"))

## Build and validate the archive in an ASCII-only staging path, then copy it to
## the requested project location. This avoids malformed ZIP entry names on Windows.
check_staging_ascii <- file.path(DIRS$temp, "CHECK_staging_ascii")
if (dir.exists(check_staging_ascii)) unlink(check_staging_ascii, recursive = TRUE, force = TRUE)
dir.create(check_staging_ascii, recursive = TRUE, showWarnings = FALSE)
check_stage_files <- list.files(DIRS$check, full.names = TRUE)
file.copy(check_stage_files, check_staging_ascii, overwrite = TRUE)

check_zip_temp <- file.path(DIRS$temp, paste0(STAGE_NAME, "_CHECK_ASCII.zip"))
if (file.exists(check_zip_temp)) unlink(check_zip_temp, force = TRUE)
zip::zipr(
  zipfile = check_zip_temp,
  files = list.files(check_staging_ascii, full.names = TRUE),
  root = check_staging_ascii,
  include_directories = FALSE
)

zip_listing <- utils::unzip(check_zip_temp, list = TRUE)
expected_zip_names <- sort(basename(list.files(check_staging_ascii, full.names = TRUE)))
observed_zip_names <- sort(basename(zip_listing$Name))
if (!identical(expected_zip_names, observed_zip_names) ||
    any(zip_listing$Length <= 0L & !grepl("/$", zip_listing$Name))) {
  stop("CHECK ZIP validation failed: archive entries do not match staged files.")
}
if (file.exists(CHECK_ZIP)) unlink(CHECK_ZIP, force = TRUE)
if (!file.copy(check_zip_temp, CHECK_ZIP, overwrite = TRUE) ||
    !file.exists(CHECK_ZIP) || file.info(CHECK_ZIP)$size <= 0L) {
  stop("Failed to copy the validated CHECK ZIP to the project root.")
}

log_msg("Stage 8 completed successfully.")
log_msg("Status: ", run_status$overall_status)
log_msg("CHECK package: ", CHECK_ZIP)

cat("\n============================================================\n")
cat("Stage 8 multicohort validation completed.\n")
cat("Status: ", run_status$overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat("============================================================\n")


