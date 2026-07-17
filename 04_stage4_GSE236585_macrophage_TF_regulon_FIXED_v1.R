############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 3 FIXED v2 REPLACE
## GSE236585 cardiac scRNA-seq discovery and projection of
## Stage 2 drug-opposed macrophage programs
##
## FIXED v2 changes:
##   - Deletes incomplete Stage 3 FIXED_v1 and prior v2 outputs before rerun.
##   - Fixes data.table::setorder() failure caused by -abs(expression).
##   - Uses explicit absolute-effect helper columns for signature deduplication.
##   - Fixes pseudobulk metadata construction in aggregate_sparse_counts().
##   - Saves post-QC and post-doublet checkpoint objects.
##   - Adds a fixed-path fallback for script self-archiving.
##
## Project root:
##   <HFPEF_PROJECT_DIR>
##
## Required inputs:
##   1) 0.GEO/GSE236585_RAW.tar
##   2) Stage 1 locked sample manifest
##   3) Stage 2 FIXED v2 drug-opposition ranking tables
##
## Output:
##   <HFPEF_PROJECT_DIR>/
##   03_stage3_GSE236585_scRNA_projection_FIXED_v2
##
## CHECK package:
##   <HFPEF_PROJECT_DIR>/
##   03_stage3_GSE236585_scRNA_projection_FIXED_v2_CHECK.zip
##
## Primary objectives:
##   1) Reconstruct six independent 10x scRNA-seq samples.
##   2) Perform sample-aware QC and optional doublet removal.
##   3) Cluster and annotate major cardiac cell populations.
##   4) Aggregate raw counts by biological sample and cell type.
##   5) Test HFpEF versus control by pseudobulk edgeR and
##      limma-voom sensitivity analysis.
##   6) Project Stage 2 Ccr2-positive, Ccr2-negative, and
##      cross-subset drug-opposed programs without forcing Nfkb1.
##   7) Quantify sample-level program localization and sensitivity
##      to Top50/Top100/Top150/Top200 signature size.
##   8) Recluster macrophage/monocyte cells and identify
##      reproducible macrophage-state candidates.
##
## Statistical boundary:
##   - Biological sample, not cell, is the inferential unit.
##   - Cell-level plots are descriptive.
##   - Pseudobulk and sample-level scores are used for comparisons.
##   - n = 3 HFpEF and n = 3 control samples.
##   - Results are discovery evidence, not causal validation.
##
## Recommended run:
##   source(
##     "<HFPEF_PROJECT_DIR>/
##      HFpEF_Stage3_GSE236585_scRNA_Projection_FIXED_v2_REPLACE.R",
##     encoding = "UTF-8"
##   )
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(warn = 1)
options(encoding = "UTF-8")
options(future.globals.maxSize = 12 * 1024^3)
set.seed(20260714)

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
DATA_DIR <- file.path(PROJECT_DIR, "0.GEO")

STAGE1_DIR <- file.path(
  PROJECT_DIR,
  "01_stage1_metadata_lock_FIXED_v3"
)
STAGE1_MANIFEST <- file.path(
  STAGE1_DIR,
  "01_tables",
  "01_locked_sample_manifest.csv"
)

STAGE2_DIR <- file.path(
  PROJECT_DIR,
  "02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2"
)
STAGE2_STATUS_FILE <- file.path(
  STAGE2_DIR,
  "01_tables",
  "24_stage2_run_status.csv"
)
STAGE2_POS_FILE <- file.path(
  STAGE2_DIR,
  "01_tables",
  "13_opposition_rank_Ccr2_positive.csv.gz"
)
STAGE2_NEG_FILE <- file.path(
  STAGE2_DIR,
  "01_tables",
  "14_opposition_rank_Ccr2_negative.csv.gz"
)
STAGE2_CROSS_FILE <- file.path(
  STAGE2_DIR,
  "01_tables",
  "16_cross_subset_consensus_ranking.csv.gz"
)

RAW_TAR_FILE <- file.path(
  DATA_DIR,
  "GSE236585_RAW.tar"
)

OLD_STAGE_NAME <- "03_stage3_GSE236585_scRNA_projection_FIXED_v1"
OLD_OUT_DIR <- file.path(PROJECT_DIR, OLD_STAGE_NAME)
OLD_CHECK_ZIP <- file.path(
  PROJECT_DIR,
  paste0(OLD_STAGE_NAME, "_CHECK.zip")
)

STAGE_NAME <- "03_stage3_GSE236585_scRNA_projection_FIXED_v2"
OUT_DIR <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP <- file.path(
  PROJECT_DIR,
  paste0(STAGE_NAME, "_CHECK.zip")
)

EXPECTED_SCRIPT_FILE <- file.path(
  PROJECT_DIR,
  "R",
  "03a_stage3_GSE236585_scRNA_projection_FIXED_v2.R"
)

REPLACE_EXISTING_STAGE3 <- TRUE
RUN_SCDOUBLETFINDER <- TRUE

## QC parameters.
MIN_FEATURES_HARD <- 200L
MIN_CELLS_PER_GENE <- 3L
MAX_FEATURES_HARD <- 9000L
MAX_COUNTS_HARD <- 120000L
MAX_PERCENT_MT_HARD <- 25
QC_MAD_MULTIPLIER <- 4
QC_UPPER_QUANTILE <- 0.995

## Clustering parameters.
N_VARIABLE_FEATURES <- 3000L
N_PCS <- 40L
DIMS_USE <- 1:30
MAJOR_CLUSTER_RESOLUTION <- 0.55
MACROPHAGE_CLUSTER_RESOLUTION <- 0.65
MAX_CELLS_PER_CLUSTER_MARKER_TEST <- 3000L

## Pseudobulk eligibility.
MIN_CELLS_PER_SAMPLE_CELLTYPE <- 20L
MIN_SAMPLES_PER_CONDITION <- 3L
PSEUDOBULK_MIN_COUNT <- 10L

## Stage 2 signature construction.
SIGNATURE_SIZES <- c(50L, 100L, 150L, 200L)
PRIMARY_SIGNATURE_SIZE <- 150L
SIGNATURE_TIERS <- c(
  "Tier_A_both_DESeq2_FDR_and_edgeR_direction",
  "Tier_B_one_DESeq2_FDR_effect_supported",
  "Tier_C_effect_and_method_supported"
)

## Reporting.
FORMAL_FDR <- 0.05
EXPLORATORY_FDR <- 0.10
TOP_DE_GENES_PER_CELLTYPE_CHECK <- 100L
TOP_MARKERS_PER_CLUSTER_CHECK <- 20L

############################################################
## 1. Script detection, validation, and output setup
############################################################

detect_script_file <- function() {
  candidates <- character()

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- tryCatch(
      frames[[i]]$ofile,
      error = function(e) NULL
    )
    if (
      !is.null(ofile) &&
      length(ofile) == 1L &&
      nzchar(ofile)
    ) {
      candidates <- c(candidates, ofile)
    }
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    candidates <- c(
      candidates,
      sub("^--file=", "", file_arg[1L])
    )
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
if (
  (
    length(SCRIPT_FILE) != 1L ||
    is.na(SCRIPT_FILE) ||
    !file.exists(SCRIPT_FILE)
  ) &&
  file.exists(EXPECTED_SCRIPT_FILE)
) {
  SCRIPT_FILE <- normalizePath(
    EXPECTED_SCRIPT_FILE,
    winslash = "/",
    mustWork = TRUE
  )
}

required_inputs <- c(
  PROJECT_DIR,
  DATA_DIR,
  STAGE1_MANIFEST,
  STAGE2_STATUS_FILE,
  STAGE2_POS_FILE,
  STAGE2_NEG_FILE,
  STAGE2_CROSS_FILE,
  RAW_TAR_FILE
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Required input path(s) are missing:\n",
    paste(missing_inputs, collapse = "\n")
  )
}

stage2_status <- data.table::fread(
  STAGE2_STATUS_FILE,
  encoding = "UTF-8"
)
if (
  !"overall_status" %in% names(stage2_status) ||
  stage2_status$overall_status[1L] !=
    "COMPLETED_STAGE2_READY_FOR_REVIEW"
) {
  stop(
    "Stage 2 is not locked as COMPLETED_STAGE2_READY_FOR_REVIEW. ",
    "Review 24_stage2_run_status.csv before Stage 3."
  )
}

replacement_audit <- data.frame(
  path = c(
    OLD_OUT_DIR,
    OLD_CHECK_ZIP,
    OUT_DIR,
    CHECK_ZIP
  ),
  path_type = c(
    "prior_incomplete_stage3_v1_output_directory",
    "prior_incomplete_stage3_v1_check_zip",
    "current_stage3_v2_output_directory",
    "current_stage3_v2_check_zip"
  ),
  existed_before = FALSE,
  deletion_attempted = FALSE,
  deletion_succeeded = FALSE,
  stringsAsFactors = FALSE
)

if (REPLACE_EXISTING_STAGE3) {
  for (i in seq_len(nrow(replacement_audit))) {
    target <- replacement_audit$path[i]
    existed <- dir.exists(target) || file.exists(target)
    replacement_audit$existed_before[i] <- existed

    if (existed) {
      replacement_audit$deletion_attempted[i] <- TRUE
      unlink(
        target,
        recursive = dir.exists(target),
        force = TRUE
      )
      replacement_audit$deletion_succeeded[i] <- !(
        dir.exists(target) || file.exists(target)
      )
      if (!replacement_audit$deletion_succeeded[i]) {
        stop(
          "Failed to remove previous Stage 3 path:\n",
          target
        )
      }
    } else {
      replacement_audit$deletion_succeeded[i] <- TRUE
    }
  }
} else if (
  dir.exists(OLD_OUT_DIR) ||
  file.exists(OLD_CHECK_ZIP) ||
  dir.exists(OUT_DIR) ||
  file.exists(CHECK_ZIP)
) {
  stop(
    "Existing Stage 3 output was detected while replacement is disabled."
  )
}

DIRS <- list(
  logs = file.path(OUT_DIR, "00_logs"),
  tables = file.path(OUT_DIR, "01_tables"),
  objects = file.path(OUT_DIR, "02_objects"),
  figures = file.path(OUT_DIR, "03_figures"),
  source = file.path(OUT_DIR, "04_source_data"),
  methods = file.path(OUT_DIR, "05_methods"),
  check = file.path(OUT_DIR, "06_review_check"),
  extracted = file.path(OUT_DIR, "07_extracted_10x")
)
for (d in c(OUT_DIR, unlist(DIRS, use.names = FALSE))) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

write.csv(
  replacement_audit,
  file.path(
    DIRS$logs,
    "stage3_replacement_audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

LOG_FILE <- file.path(
  DIRS$logs,
  "stage3_GSE236585.log"
)
WARN_FILE <- file.path(
  DIRS$logs,
  "stage3_warnings.log"
)
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
  cat(
    line,
    "\n",
    file = LOG_FILE,
    append = TRUE
  )
  invisible(line)
}

warning_records <- list()
add_warning <- function(category, item, message) {
  rec <- data.frame(
    timestamp = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
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

  log_msg(
    category,
    " | ",
    item,
    " | ",
    message,
    level = "WARN"
  )
  invisible(rec)
}

log_msg("Stage 3 GSE236585 analysis started.")
log_msg(
  "Replacement of incomplete Stage 3 outputs enabled: ",
  REPLACE_EXISTING_STAGE3
)
for (i in seq_len(nrow(replacement_audit))) {
  log_msg(
    "Replacement audit | ",
    replacement_audit$path_type[i],
    " | existed_before=",
    replacement_audit$existed_before[i],
    " | deletion_succeeded=",
    replacement_audit$deletion_succeeded[i],
    " | path=",
    replacement_audit$path[i]
  )
}
log_msg("PROJECT_DIR: ", PROJECT_DIR)
log_msg("RAW_TAR_FILE: ", RAW_TAR_FILE)
log_msg("STAGE2_DIR: ", STAGE2_DIR)
log_msg("OUT_DIR: ", OUT_DIR)
log_msg(
  "Detected SCRIPT_FILE: ",
  ifelse(
    length(SCRIPT_FILE) == 1L &&
      !is.na(SCRIPT_FILE) &&
      file.exists(SCRIPT_FILE),
    SCRIPT_FILE,
    "NOT_DETECTED"
  )
)

############################################################
## 2. Package setup
############################################################

ensure_cran <- function(pkgs, required = TRUE) {
  missing <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
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
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(still_missing) > 0L) {
    msg <- paste(
      "CRAN package(s) unavailable:",
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

ensure_bioc <- function(pkgs, required = TRUE) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }

  missing <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing) > 0L) {
    log_msg(
      "Installing missing Bioconductor package(s): ",
      paste(missing, collapse = ", ")
    )
    try(
      BiocManager::install(
        missing,
        ask = FALSE,
        update = FALSE
      ),
      silent = TRUE
    )
  }

  still_missing <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(still_missing) > 0L) {
    msg <- paste(
      "Bioconductor package(s) unavailable:",
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

ensure_cran(
  c(
    "Seurat",
    "SeuratObject",
    "data.table",
    "Matrix",
    "ggplot2",
    "ggrepel",
    "patchwork",
    "pheatmap",
    "openxlsx",
    "scales",
    "zip",
    "digest",
    "matrixStats"
  ),
  required = TRUE
)

ensure_bioc(
  c(
    "edgeR",
    "limma"
  ),
  required = TRUE
)

if (RUN_SCDOUBLETFINDER) {
  ensure_bioc(
    c(
      "SingleCellExperiment",
      "scDblFinder"
    ),
    required = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(openxlsx)
  library(edgeR)
  library(limma)
})

############################################################
## 3. General utilities
############################################################

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\r|\\n", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

gene_key <- function(x) {
  x <- normalize_text(x)
  x <- sub(
    "([._-][0-9]+)$",
    "",
    x
  )
  toupper(x)
}

extract_gsm <- function(x) {
  x <- as.character(x)
  m <- regexpr(
    "GSM[0-9]+",
    x,
    ignore.case = TRUE
  )
  out <- rep(NA_character_, length(x))
  valid <- m > 0L
  out[valid] <- toupper(
    regmatches(x, m)[valid]
  )
  out
}

write_csv_safe <- function(
  x,
  path,
  compress = FALSE
) {
  if (
    is.null(x) ||
    ncol(x) == 0L
  ) {
    fwrite(
      data.table(
        note = "No records generated."
      ),
      path
    )
  } else {
    fwrite(
      x,
      path,
      compress = if (compress) "gzip" else "none"
    )
  }
}

sanitize_sheet_name <- function(x) {
  x <- gsub(
    "[\\[\\]:*?/\\\\]",
    "_",
    x
  )
  substr(x, 1L, 31L)
}

write_sheet_safe <- function(wb, sheet, x) {
  sheet <- sanitize_sheet_name(sheet)
  addWorksheet(wb, sheet)

  if (
    is.null(x) ||
    nrow(x) == 0L ||
    ncol(x) == 0L
  ) {
    writeData(
      wb,
      sheet,
      data.frame(
        note = "No records generated."
      )
    )
    return(invisible(NULL))
  }

  y <- as.data.frame(
    x,
    stringsAsFactors = FALSE
  )

  char_cols <- vapply(
    y,
    is.character,
    logical(1)
  )
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
    cols = seq_len(min(ncol(y), 35L)),
    widths = "auto"
  )
  invisible(NULL)
}

save_plot_bundle <- function(
  plot_object,
  stem,
  width,
  height
) {
  png_path <- file.path(
    DIRS$figures,
    paste0(stem, ".png")
  )
  pdf_path <- file.path(
    DIRS$figures,
    paste0(stem, ".pdf")
  )
  tiff_path <- file.path(
    DIRS$figures,
    paste0(stem, ".tiff")
  )

  ggsave(
    png_path,
    plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    pdf_path,
    plot_object,
    width = width,
    height = height,
    bg = "white"
  )
  ggsave(
    tiff_path,
    plot_object,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )

  invisible(
    c(
      png = png_path,
      pdf = pdf_path,
      tiff = tiff_path
    )
  )
}

join_layers_safe <- function(
  object,
  assay = "RNA"
) {
  DefaultAssay(object) <- assay

  out <- tryCatch(
    {
      if (
        "JoinLayers" %in%
          getNamespaceExports("SeuratObject")
      ) {
        SeuratObject::JoinLayers(
          object,
          assay = assay
        )
      } else if (
        exists(
          "JoinLayers",
          mode = "function"
        )
      ) {
        JoinLayers(
          object,
          assay = assay
        )
      } else {
        object
      }
    },
    error = function(e) {
      add_warning(
        "SEURAT_LAYER",
        assay,
        paste0(
          "JoinLayers failed; original object retained: ",
          conditionMessage(e)
        )
      )
      object
    }
  )

  DefaultAssay(out) <- assay
  out
}

get_assay_matrix <- function(
  object,
  layer = "counts",
  assay = "RNA"
) {
  DefaultAssay(object) <- assay

  mat <- tryCatch(
    SeuratObject::LayerData(
      object,
      assay = assay,
      layer = layer
    ),
    error = function(e1) {
      tryCatch(
        GetAssayData(
          object,
          assay = assay,
          layer = layer
        ),
        error = function(e2) {
          GetAssayData(
            object,
            assay = assay,
            slot = layer
          )
        }
      )
    }
  )

  mat
}

count_lines_file <- function(path) {
  con <- if (
    grepl("\\.gz$", path, ignore.case = TRUE)
  ) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(
    try(close(con), silent = TRUE),
    add = TRUE
  )

  n <- 0L
  repeat {
    chunk <- readLines(
      con,
      n = 100000L,
      warn = FALSE
    )
    if (length(chunk) == 0L) break
    n <- n + length(chunk)
  }
  n
}

read_mtx_dimensions <- function(path) {
  con <- if (
    grepl("\\.gz$", path, ignore.case = TRUE)
  ) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(
    try(close(con), silent = TRUE),
    add = TRUE
  )

  repeat {
    line <- readLines(
      con,
      n = 1L,
      warn = FALSE
    )
    if (length(line) == 0L) {
      return(
        c(
          rows = NA_integer_,
          cols = NA_integer_,
          nnz = NA_real_
        )
      )
    }
    if (!grepl("^%", line)) {
      fields <- suppressWarnings(
        as.numeric(
          strsplit(
            trimws(line),
            "[[:space:]]+"
          )[[1L]]
        )
      )
      if (length(fields) >= 3L) {
        return(
          c(
            rows = fields[1L],
            cols = fields[2L],
            nnz = fields[3L]
          )
        )
      }
    }
  }
}

adaptive_upper_threshold <- function(
  x,
  minimum_upper,
  maximum_upper,
  mad_multiplier = QC_MAD_MULTIPLIER,
  quantile_probability = QC_UPPER_QUANTILE
) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(maximum_upper)

  med <- median(x)
  mad_value <- mad(
    x,
    center = med,
    constant = 1.4826
  )
  q_value <- as.numeric(
    quantile(
      x,
      probs = quantile_probability,
      na.rm = TRUE,
      names = FALSE
    )
  )

  candidate <- max(
    q_value,
    med + mad_multiplier * mad_value,
    minimum_upper,
    na.rm = TRUE
  )

  min(maximum_upper, candidate)
}

hedges_g <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  n1 <- length(x)
  n2 <- length(y)

  if (n1 < 2L || n2 < 2L) {
    return(NA_real_)
  }

  pooled_var <- (
    (n1 - 1L) * stats::var(x) +
      (n2 - 1L) * stats::var(y)
  ) / (n1 + n2 - 2L)

  if (
    !is.finite(pooled_var) ||
    pooled_var <= 0
  ) {
    return(NA_real_)
  }

  d <- (
    mean(x) - mean(y)
  ) / sqrt(pooled_var)

  correction <- 1 -
    3 / (
      4 * (n1 + n2) - 9
    )

  correction * d
}

safe_wilcox_p <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (
    length(x) < 1L ||
    length(y) < 1L
  ) {
    return(NA_real_)
  }

  tryCatch(
    wilcox.test(
      x,
      y,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)

  suppressWarnings(
    cor(
      x[keep],
      y[keep],
      method = "spearman"
    )
  )
}

safe_pearson <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)

  suppressWarnings(
    cor(
      x[keep],
      y[keep],
      method = "pearson"
    )
  )
}

edgeR_norm_lib_sizes <- function(y) {
  if (
    "normLibSizes" %in%
      getNamespaceExports("edgeR")
  ) {
    edgeR::normLibSizes(y)
  } else {
    edgeR::calcNormFactors(y)
  }
}

############################################################
## 4. Read Stage 1 metadata and Stage 2 ranked programs
############################################################

locked_manifest <- fread(
  STAGE1_MANIFEST,
  encoding = "UTF-8",
  na.strings = c("", "NA", "NaN")
)

required_manifest_columns <- c(
  "dataset_id",
  "sample_accession",
  "original_title",
  "condition",
  "group_id",
  "lock_status"
)
missing_manifest_columns <- setdiff(
  required_manifest_columns,
  names(locked_manifest)
)
if (length(missing_manifest_columns) > 0L) {
  stop(
    "Stage 1 manifest is missing column(s): ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}

sample_meta <- locked_manifest[
  dataset_id == "GSE236585"
]
if (
  nrow(sample_meta) != 6L ||
  uniqueN(sample_meta$sample_accession) != 6L
) {
  stop(
    "Expected six unique GSE236585 samples in the Stage 1 manifest."
  )
}
if (any(!grepl("^LOCKED", sample_meta$lock_status))) {
  stop(
    "At least one GSE236585 sample is not metadata-locked."
  )
}

sample_meta[, condition := factor(
  condition,
  levels = c("Control", "HFpEF")
)]
if (
  any(is.na(sample_meta$condition)) ||
  any(table(sample_meta$condition) != 3L)
) {
  stop(
    "Expected three Control and three HFpEF samples."
  )
}

setorder(
  sample_meta,
  condition,
  sample_accession
)
fwrite(
  sample_meta,
  file.path(
    DIRS$tables,
    "01_locked_GSE236585_sample_metadata.csv"
  )
)

stage2_pos <- fread(
  STAGE2_POS_FILE,
  encoding = "UTF-8"
)
stage2_neg <- fread(
  STAGE2_NEG_FILE,
  encoding = "UTF-8"
)
stage2_cross <- fread(
  STAGE2_CROSS_FILE,
  encoding = "UTF-8"
)

required_stage2_columns <- c(
  "symbol",
  "disease_lfc",
  "drug_lfc",
  "deseq_opposed",
  "edger_opposed",
  "four_effect_signs_consistent",
  "opposition_tier",
  "within_subset_rank"
)
for (nm in required_stage2_columns) {
  if (
    !nm %in% names(stage2_pos) ||
    !nm %in% names(stage2_neg)
  ) {
    stop(
      "Stage 2 opposition table is missing required column: ",
      nm
    )
  }
}

derive_subset_signature_table <- function(
  x,
  subset_label
) {
  y <- copy(x)

  y <- y[
    !is.na(symbol) &
      nzchar(symbol) &
      opposition_tier %in% SIGNATURE_TIERS &
      deseq_opposed == TRUE &
      edger_opposed == TRUE &
      four_effect_signs_consistent == TRUE &
      is.finite(disease_lfc) &
      is.finite(drug_lfc)
  ]

  y[, symbol_key := gene_key(symbol)]
  y <- y[
    nzchar(symbol_key)
  ]
  ## data.table::setorder() accepts column names, not expressions.
  ## Explicit helper columns prevent the previous -abs(expression) failure.
  y[, abs_disease_lfc_for_order := abs(disease_lfc)]
  y[, abs_drug_lfc_for_order := abs(drug_lfc)]
  setorder(
    y,
    within_subset_rank,
    -abs_disease_lfc_for_order,
    -abs_drug_lfc_for_order
  )
  y <- y[, .SD[1L], by = symbol_key]
  y[
    ,
    c(
      "abs_disease_lfc_for_order",
      "abs_drug_lfc_for_order"
    ) := NULL
  ]

  y[, direction := fcase(
    disease_lfc > 0 &
      drug_lfc < 0,
    "Disease_up_Drug_down",

    disease_lfc < 0 &
      drug_lfc > 0,
    "Disease_down_Drug_up",

    default = "Other"
  )]
  y <- y[direction != "Other"]
  y[, subset_source := subset_label]
  y
}

pos_signature_base <- derive_subset_signature_table(
  stage2_pos,
  "Ccr2_positive"
)
neg_signature_base <- derive_subset_signature_table(
  stage2_neg,
  "Ccr2_negative"
)

required_cross_columns <- c(
  "symbol",
  "pos_disease_lfc",
  "pos_drug_lfc",
  "neg_disease_lfc",
  "neg_drug_lfc",
  "consensus_category",
  "overall_consensus_rank"
)
missing_cross_columns <- setdiff(
  required_cross_columns,
  names(stage2_cross)
)
if (length(missing_cross_columns) > 0L) {
  stop(
    "Cross-subset Stage 2 table is missing column(s): ",
    paste(missing_cross_columns, collapse = ", ")
  )
}

cross_signature_base <- copy(stage2_cross)
cross_signature_base <- cross_signature_base[
  !is.na(symbol) &
    nzchar(symbol) &
    consensus_category ==
      "Cross_subset_full_directional_consensus" &
    is.finite(pos_disease_lfc) &
    is.finite(pos_drug_lfc) &
    is.finite(neg_disease_lfc) &
    is.finite(neg_drug_lfc)
]
cross_signature_base[, symbol_key := gene_key(symbol)]
setorder(
  cross_signature_base,
  overall_consensus_rank
)
cross_signature_base <- cross_signature_base[
  ,
  .SD[1L],
  by = symbol_key
]
cross_signature_base[, direction := fcase(
  pos_disease_lfc > 0 &
    pos_drug_lfc < 0 &
    neg_disease_lfc > 0 &
    neg_drug_lfc < 0,
  "Disease_up_Drug_down",

  pos_disease_lfc < 0 &
    pos_drug_lfc > 0 &
    neg_disease_lfc < 0 &
    neg_drug_lfc > 0,
  "Disease_down_Drug_up",

  default = "Other"
)]
cross_signature_base <- cross_signature_base[
  direction != "Other"
]
cross_signature_base[, subset_source := "Cross_subset"]
cross_signature_base[
  ,
  disease_lfc := rowMeans(
    cbind(
      pos_disease_lfc,
      neg_disease_lfc
    ),
    na.rm = TRUE
  )
]
cross_signature_base[
  ,
  drug_lfc := rowMeans(
    cbind(
      pos_drug_lfc,
      neg_drug_lfc
    ),
    na.rm = TRUE
  )
]
cross_signature_base[
  ,
  within_subset_rank :=
    overall_consensus_rank
]

build_signature_sets <- function(
  base_table,
  prefix
) {
  out <- list()
  manifest_records <- list()

  for (n_target in SIGNATURE_SIZES) {
    up_dt <- base_table[
      direction == "Disease_up_Drug_down"
    ][
      order(within_subset_rank)
    ][
      seq_len(min(n_target, .N))
    ]

    down_dt <- base_table[
      direction == "Disease_down_Drug_up"
    ][
      order(within_subset_rank)
    ][
      seq_len(min(n_target, .N))
    ]

    set_name <- paste0(
      prefix,
      "_Top",
      n_target
    )

    out[[set_name]] <- list(
      up = unique(up_dt$symbol),
      down = unique(down_dt$symbol),
      requested_size_per_direction = n_target,
      source = prefix
    )

    if (nrow(up_dt) > 0L) {
      manifest_records[[length(manifest_records) + 1L]] <- up_dt[
        ,
        .(
          signature_name = set_name,
          signature_source = prefix,
          requested_size_per_direction = n_target,
          direction,
          symbol,
          symbol_key,
          stage2_disease_lfc = disease_lfc,
          stage2_drug_lfc = drug_lfc,
          stage2_rank = within_subset_rank
        )
      ]
    }

    if (nrow(down_dt) > 0L) {
      manifest_records[[length(manifest_records) + 1L]] <- down_dt[
        ,
        .(
          signature_name = set_name,
          signature_source = prefix,
          requested_size_per_direction = n_target,
          direction,
          symbol,
          symbol_key,
          stage2_disease_lfc = disease_lfc,
          stage2_drug_lfc = drug_lfc,
          stage2_rank = within_subset_rank
        )
      ]
    }
  }

  list(
    sets = out,
    manifest = if (
      length(manifest_records) > 0L
    ) {
      rbindlist(
        manifest_records,
        use.names = TRUE,
        fill = TRUE
      )
    } else {
      data.table()
    }
  )
}

pos_sets <- build_signature_sets(
  pos_signature_base,
  "Ccr2pos"
)
neg_sets <- build_signature_sets(
  neg_signature_base,
  "Ccr2neg"
)
cross_sets <- build_signature_sets(
  cross_signature_base,
  "CrossSubset"
)

signature_sets <- c(
  pos_sets$sets,
  neg_sets$sets,
  cross_sets$sets
)
signature_manifest <- rbindlist(
  list(
    pos_sets$manifest,
    neg_sets$manifest,
    cross_sets$manifest
  ),
  use.names = TRUE,
  fill = TRUE
)

if (
  length(signature_sets) !=
    length(SIGNATURE_SIZES) * 3L
) {
  stop(
    "The expected signature-set collection was not generated."
  )
}

fwrite(
  signature_manifest,
  file.path(
    DIRS$tables,
    "02_stage2_signature_gene_manifest.csv"
  )
)

signature_size_summary <- signature_manifest[
  ,
  .(
    selected_genes = uniqueN(symbol_key)
  ),
  by = .(
    signature_name,
    signature_source,
    requested_size_per_direction,
    direction
  )
]
fwrite(
  signature_size_summary,
  file.path(
    DIRS$tables,
    "03_stage2_signature_size_summary.csv"
  )
)

############################################################
## 5. Map, extract, and validate 10x files
############################################################

tar_members <- tryCatch(
  utils::untar(
    RAW_TAR_FILE,
    list = TRUE
  ),
  error = function(e) {
    stop(
      "Unable to list GSE236585_RAW.tar: ",
      conditionMessage(e)
    )
  }
)

if (length(tar_members) != 18L) {
  add_warning(
    "ARCHIVE",
    "GSE236585_RAW.tar",
    paste0(
      "Expected 18 archive members but found ",
      length(tar_members),
      ". Mapping will continue using GSM accessions."
    )
  )
}

mapping_records <- list()
for (gsm in sample_meta$sample_accession) {
  matrix_hit <- tar_members[
    grepl(gsm, tar_members, fixed = TRUE) &
      grepl(
        "_matrix\\.mtx\\.gz$",
        tar_members,
        ignore.case = TRUE
      )
  ]
  feature_hit <- tar_members[
    grepl(gsm, tar_members, fixed = TRUE) &
      grepl(
        "_features\\.tsv\\.gz$",
        tar_members,
        ignore.case = TRUE
      )
  ]
  barcode_hit <- tar_members[
    grepl(gsm, tar_members, fixed = TRUE) &
      grepl(
        "_barcodes\\.tsv\\.gz$",
        tar_members,
        ignore.case = TRUE
      )
  ]

  mapping_records[[length(mapping_records) + 1L]] <- data.table(
    sample_accession = gsm,
    matrix_member = paste(
      matrix_hit,
      collapse = "; "
    ),
    feature_member = paste(
      feature_hit,
      collapse = "; "
    ),
    barcode_member = paste(
      barcode_hit,
      collapse = "; "
    ),
    matrix_match_count = length(matrix_hit),
    feature_match_count = length(feature_hit),
    barcode_match_count = length(barcode_hit),
    mapping_status = if (
      length(matrix_hit) == 1L &&
      length(feature_hit) == 1L &&
      length(barcode_hit) == 1L
    ) {
      "PASS"
    } else {
      "FAIL"
    }
  )
}

file_map <- rbindlist(
  mapping_records,
  use.names = TRUE,
  fill = TRUE
)
file_map <- merge(
  sample_meta[
    ,
    .(
      sample_accession,
      original_title,
      condition,
      group_id
    )
  ],
  file_map,
  by = "sample_accession",
  all.x = TRUE
)

if (any(file_map$mapping_status != "PASS")) {
  fwrite(
    file_map,
    file.path(
      DIRS$tables,
      "FATAL_10x_file_mapping.csv"
    )
  )
  stop(
    "At least one GSE236585 sample did not map uniquely to a 10x triple."
  )
}

members_to_extract <- unique(
  c(
    file_map$matrix_member,
    file_map$feature_member,
    file_map$barcode_member
  )
)

log_msg(
  "Extracting ",
  length(members_to_extract),
  " GSE236585 10x files."
)

utils::untar(
  RAW_TAR_FILE,
  files = members_to_extract,
  exdir = DIRS$extracted
)

file_map[
  ,
  matrix_file := file.path(
    DIRS$extracted,
    matrix_member
  )
]
file_map[
  ,
  feature_file := file.path(
    DIRS$extracted,
    feature_member
  )
]
file_map[
  ,
  barcode_file := file.path(
    DIRS$extracted,
    barcode_member
  )
]

file_map[
  ,
  files_exist := (
    file.exists(matrix_file) &
      file.exists(feature_file) &
      file.exists(barcode_file)
  )
]

if (any(!file_map$files_exist)) {
  stop(
    "At least one extracted 10x file does not exist."
  )
}

dimension_records <- list()
for (i in seq_len(nrow(file_map))) {
  dims <- read_mtx_dimensions(
    file_map$matrix_file[i]
  )
  n_features <- count_lines_file(
    file_map$feature_file[i]
  )
  n_barcodes <- count_lines_file(
    file_map$barcode_file[i]
  )

  dimension_records[[i]] <- data.table(
    sample_accession =
      file_map$sample_accession[i],
    matrix_rows = as.integer(
      dims["rows"]
    ),
    matrix_columns = as.integer(
      dims["cols"]
    ),
    matrix_nonzero_entries = as.numeric(
      dims["nnz"]
    ),
    feature_lines = n_features,
    barcode_lines = n_barcodes,
    rows_match_features = (
      as.integer(dims["rows"]) ==
        n_features
    ),
    columns_match_barcodes = (
      as.integer(dims["cols"]) ==
        n_barcodes
    )
  )
}

dimension_validation <- rbindlist(
  dimension_records
)
dimension_validation[
  ,
  dimensions_valid := (
    rows_match_features &
      columns_match_barcodes
  )
]

file_map <- merge(
  file_map,
  dimension_validation,
  by = "sample_accession",
  all.x = TRUE
)

fwrite(
  file_map,
  file.path(
    DIRS$tables,
    "04_GSE236585_10x_file_mapping_and_dimensions.csv"
  )
)

if (any(!file_map$dimensions_valid)) {
  stop(
    "At least one GSE236585 10x triple failed dimension validation."
  )
}

############################################################
## 6. Build sample Seurat objects and perform adaptive QC
############################################################

sample_objects <- list()
qc_threshold_records <- list()
qc_summary_records <- list()

for (i in seq_len(nrow(file_map))) {
  gsm <- file_map$sample_accession[i]
  condition_i <- as.character(
    file_map$condition[i]
  )

  log_msg(
    "Reading sample ",
    i,
    "/",
    nrow(file_map),
    ": ",
    gsm,
    " | ",
    condition_i
  )

  counts <- Seurat::ReadMtx(
    mtx = file_map$matrix_file[i],
    cells = file_map$barcode_file[i],
    features = file_map$feature_file[i],
    feature.column = 2L,
    cell.column = 1L,
    unique.features = TRUE,
    strip.suffix = FALSE
  )

  if (is.list(counts)) {
    counts <- counts[[1L]]
  }

  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE236585",
    min.cells = MIN_CELLS_PER_GENE,
    min.features = 100L
  )

  obj <- RenameCells(
    obj,
    add.cell.id = gsm
  )
  obj$dataset_id <- "GSE236585"
  obj$sample_accession <- gsm
  obj$condition <- condition_i
  obj$original_title <-
    file_map$original_title[i]

  obj[["percent.mt"]] <-
    PercentageFeatureSet(
      obj,
      pattern = "^mt-|^Mt-|^MT-"
    )

  feature_upper <- adaptive_upper_threshold(
    obj$nFeature_RNA,
    minimum_upper = 2500,
    maximum_upper = MAX_FEATURES_HARD
  )
  count_upper <- adaptive_upper_threshold(
    obj$nCount_RNA,
    minimum_upper = 30000,
    maximum_upper = MAX_COUNTS_HARD
  )
  mt_upper <- adaptive_upper_threshold(
    obj$percent.mt,
    minimum_upper = 12,
    maximum_upper = MAX_PERCENT_MT_HARD
  )

  qc_threshold_records[[gsm]] <- data.table(
    sample_accession = gsm,
    condition = condition_i,
    minimum_features = MIN_FEATURES_HARD,
    maximum_features = round(feature_upper, 2),
    maximum_counts = round(count_upper, 2),
    maximum_percent_mt = round(mt_upper, 2),
    threshold_method = paste0(
      "max(sample 99.5th percentile, median + ",
      QC_MAD_MULTIPLIER,
      " MAD) constrained by hard ceilings"
    )
  )

  before_n <- ncol(obj)

  keep_cells <- rownames(
    obj@meta.data
  )[
    obj$nFeature_RNA >=
      MIN_FEATURES_HARD &
      obj$nFeature_RNA <=
        feature_upper &
      obj$nCount_RNA <=
        count_upper &
      obj$percent.mt <=
        mt_upper
  ]

  obj <- subset(
    obj,
    cells = keep_cells
  )

  after_n <- ncol(obj)

  qc_summary_records[[gsm]] <- data.table(
    sample_accession = gsm,
    condition = condition_i,
    cells_before_qc = before_n,
    cells_after_qc = after_n,
    retained_fraction = after_n / before_n,
    median_features_after = median(
      obj$nFeature_RNA
    ),
    median_counts_after = median(
      obj$nCount_RNA
    ),
    median_percent_mt_after = median(
      obj$percent.mt
    )
  )

  if (after_n < 100L) {
    stop(
      "Sample ",
      gsm,
      " retained fewer than 100 cells after QC."
    )
  }

  sample_objects[[gsm]] <- obj
}

qc_thresholds <- rbindlist(
  qc_threshold_records
)
qc_summary <- rbindlist(
  qc_summary_records
)

fwrite(
  qc_thresholds,
  file.path(
    DIRS$tables,
    "05_sample_specific_QC_thresholds.csv"
  )
)
fwrite(
  qc_summary,
  file.path(
    DIRS$tables,
    "06_sample_QC_retention_summary.csv"
  )
)

if (length(sample_objects) == 1L) {
  cardiac <- sample_objects[[1L]]
} else {
  cardiac <- merge(
    sample_objects[[1L]],
    y = sample_objects[-1L],
    project = "GSE236585"
  )
}
cardiac <- join_layers_safe(
  cardiac,
  assay = "RNA"
)

saveRDS(
  cardiac,
  file.path(
    DIRS$objects,
    "CHECKPOINT_GSE236585_post_QC_pre_doublet.rds"
  )
)

############################################################
## 7. Optional scDblFinder doublet removal
############################################################

doublet_status <- "NOT_REQUESTED"
doublet_summary <- data.table()

if (RUN_SCDOUBLETFINDER) {
  if (
    requireNamespace(
      "SingleCellExperiment",
      quietly = TRUE
    ) &&
    requireNamespace(
      "scDblFinder",
      quietly = TRUE
    )
  ) {
    log_msg("Running optional scDblFinder.")

    doublet_result <- tryCatch(
      {
        sce <- Seurat::as.SingleCellExperiment(
          cardiac,
          assay = "RNA"
        )
        SingleCellExperiment::colData(
          sce
        )$sample_accession <-
          cardiac$sample_accession

        sce <- scDblFinder::scDblFinder(
          sce,
          samples = "sample_accession"
        )

        class_values <- as.character(
          SingleCellExperiment::colData(
            sce
          )$scDblFinder.class
        )
        score_values <- as.numeric(
          SingleCellExperiment::colData(
            sce
          )$scDblFinder.score
        )

        names(class_values) <- colnames(sce)
        names(score_values) <- colnames(sce)

        list(
          class = class_values,
          score = score_values
        )
      },
      error = function(e) {
        add_warning(
          "DOUBLETS",
          "scDblFinder",
          paste0(
            "scDblFinder failed and was skipped: ",
            conditionMessage(e)
          )
        )
        NULL
      }
    )

    if (!is.null(doublet_result)) {
      cardiac$doublet_class <-
        doublet_result$class[
          colnames(cardiac)
        ]
      cardiac$doublet_score <-
        doublet_result$score[
          colnames(cardiac)
        ]

      doublet_summary <- as.data.table(
        cardiac@meta.data,
        keep.rownames = "cell"
      )[
        ,
        .(
          cells = .N
        ),
        by = .(
          sample_accession,
          condition,
          doublet_class
        )
      ]
      doublet_summary[
        ,
        fraction := cells / sum(cells),
        by = sample_accession
      ]

      singlet_cells <- rownames(
        cardiac@meta.data
      )[
        cardiac$doublet_class ==
          "singlet"
      ]

      cardiac <- subset(
        cardiac,
        cells = singlet_cells
      )
      cardiac <- join_layers_safe(
        cardiac,
        assay = "RNA"
      )

      doublet_status <- "COMPLETED"
    } else {
      cardiac$doublet_class <- "not_run"
      cardiac$doublet_score <- NA_real_
      doublet_status <- "FAILED_NONFATAL"
    }
  } else {
    cardiac$doublet_class <- "package_unavailable"
    cardiac$doublet_score <- NA_real_
    doublet_status <- "PACKAGE_UNAVAILABLE"
    add_warning(
      "DOUBLETS",
      "packages",
      "SingleCellExperiment and/or scDblFinder unavailable; doublet removal was skipped."
    )
  }
} else {
  cardiac$doublet_class <- "not_requested"
  cardiac$doublet_score <- NA_real_
}

write_csv_safe(
  doublet_summary,
  file.path(
    DIRS$tables,
    "07_scDblFinder_summary.csv"
  )
)

post_doublet_counts <- as.data.table(
  cardiac@meta.data,
  keep.rownames = "cell"
)[
  ,
  .(
    retained_cells = .N
  ),
  by = .(
    sample_accession,
    condition
  )
]
fwrite(
  post_doublet_counts,
  file.path(
    DIRS$tables,
    "08_post_doublet_cell_counts.csv"
  )
)

saveRDS(
  cardiac,
  file.path(
    DIRS$objects,
    "CHECKPOINT_GSE236585_post_doublet_pre_clustering.rds"
  )
)

############################################################
## 8. Normalize, cluster, and annotate major cardiac cells
############################################################

log_msg(
  "Cardiac cells entering clustering: ",
  ncol(cardiac),
  " | genes: ",
  nrow(cardiac)
)

cardiac <- NormalizeData(
  cardiac,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
cardiac <- join_layers_safe(
  cardiac,
  assay = "RNA"
)
cardiac <- FindVariableFeatures(
  cardiac,
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES,
  verbose = FALSE
)
cardiac <- ScaleData(
  cardiac,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)
cardiac <- RunPCA(
  cardiac,
  npcs = N_PCS,
  verbose = FALSE
)

available_pcs <- ncol(
  Embeddings(
    cardiac,
    reduction = "pca"
  )
)
dims_use <- DIMS_USE[
  DIMS_USE <= available_pcs
]
if (length(dims_use) < 10L) {
  stop(
    "Fewer than 10 principal components were available."
  )
}

cardiac <- FindNeighbors(
  cardiac,
  dims = dims_use,
  verbose = FALSE
)
cardiac <- FindClusters(
  cardiac,
  resolution = MAJOR_CLUSTER_RESOLUTION,
  algorithm = 1L,
  random.seed = 20260714,
  verbose = FALSE
)
cardiac <- RunUMAP(
  cardiac,
  dims = dims_use,
  seed.use = 20260714,
  verbose = FALSE
)

pca_stdev <- Stdev(
  cardiac,
  reduction = "pca"
)
pca_variance <- (
  pca_stdev^2 /
    sum(pca_stdev^2)
)
pca_variance_table <- data.table(
  PC = seq_along(pca_stdev),
  variance_fraction = pca_variance,
  cumulative_variance =
    cumsum(pca_variance),
  used_for_graph = seq_along(
    pca_stdev
  ) %in% dims_use
)
fwrite(
  pca_variance_table,
  file.path(
    DIRS$tables,
    "09_PCA_variance_and_dimensions.csv"
  )
)

major_marker_sets <- list(
  Cardiomyocyte = c(
    "Tnnt2", "Tnni3", "Myh6",
    "Myh7", "Actc1", "Ryr2"
  ),
  Fibroblast = c(
    "Dcn", "Lum", "Col1a1",
    "Col1a2", "Col3a1", "Pdgfra"
  ),
  Endothelial = c(
    "Pecam1", "Cdh5", "Vwf",
    "Kdr", "Flt1", "Emcn"
  ),
  Lymphatic_endothelial = c(
    "Prox1", "Lyve1", "Pdpn",
    "Flt4", "Ccl21a"
  ),
  Pericyte = c(
    "Pdgfrb", "Rgs5", "Cspg4",
    "Notch3", "Kcnj8", "Abcc9"
  ),
  Smooth_muscle = c(
    "Acta2", "Tagln", "Myh11",
    "Cnn1", "Myl9"
  ),
  Macrophage_Monocyte = c(
    "Lyz2", "Adgre1", "Csf1r",
    "Cd68", "Fcgr1", "Apoe"
  ),
  Dendritic_cell = c(
    "Flt3", "Itgax", "Clec10a",
    "H2-Ab1", "Cd74"
  ),
  Neutrophil = c(
    "S100a8", "S100a9", "Mpo",
    "Elane", "Retnlg", "Ly6g"
  ),
  T_NK = c(
    "Cd3d", "Cd3e", "Trac",
    "Nkg7", "Klrb1c", "Prf1"
  ),
  B_cell = c(
    "Ms4a1", "Cd79a", "Cd79b",
    "Cd74", "Cd37"
  ),
  Mast_cell = c(
    "Kit", "Cpa3", "Tpsb2",
    "Fcer1a"
  ),
  Epicardial_Mesothelial = c(
    "Upk3b", "Msln", "Krt19",
    "Wt1", "C3"
  )
)

calculate_cluster_marker_scores <- function(
  object,
  marker_sets,
  cluster_column = "seurat_clusters"
) {
  data_matrix <- get_assay_matrix(
    object,
    layer = "data",
    assay = "RNA"
  )
  cell_clusters <- as.character(
    object@meta.data[
      colnames(data_matrix),
      cluster_column
    ]
  )
  cluster_factor <- factor(
    cell_clusters,
    levels = sort(
      unique(cell_clusters)
    )
  )

  design <- Matrix::sparse.model.matrix(
    ~ 0 + cluster_factor
  )
  colnames(design) <- levels(
    cluster_factor
  )

  cluster_sums <- data_matrix %*% design
  cluster_counts <- as.numeric(
    table(cluster_factor)[
      levels(cluster_factor)
    ]
  )
  cluster_means <- sweep(
    cluster_sums,
    2L,
    cluster_counts,
    "/"
  )

  feature_keys <- gene_key(
    rownames(cluster_means)
  )
  score_records <- list()
  availability_records <- list()

  for (set_name in names(marker_sets)) {
    requested <- marker_sets[[set_name]]
    requested_keys <- gene_key(requested)
    idx <- which(
      feature_keys %in% requested_keys
    )

    availability_records[[set_name]] <- data.table(
      marker_set = set_name,
      requested_genes =
        length(unique(requested_keys)),
      detected_genes = length(idx),
      detected_gene_symbols = paste(
        rownames(cluster_means)[idx],
        collapse = ";"
      )
    )

    if (length(idx) == 0L) {
      score_value <- rep(
        NA_real_,
        ncol(cluster_means)
      )
    } else {
      score_value <- Matrix::colMeans(
        cluster_means[
          idx,
          ,
          drop = FALSE
        ]
      )
    }

    score_records[[set_name]] <- data.table(
      seurat_cluster =
        colnames(cluster_means),
      marker_set = set_name,
      raw_marker_score =
        as.numeric(score_value)
    )
  }

  scores <- rbindlist(
    score_records
  )
  availability <- rbindlist(
    availability_records
  )

  scores[
    ,
    standardized_marker_score := {
      if (
        all(!is.finite(raw_marker_score)) ||
        stats::sd(
          raw_marker_score,
          na.rm = TRUE
        ) == 0
      ) {
        rep(0, .N)
      } else {
        as.numeric(
          scale(raw_marker_score)
        )
      }
    },
    by = marker_set
  ]

  setorder(
    scores,
    seurat_cluster,
    -standardized_marker_score,
    -raw_marker_score
  )

  annotation <- scores[
    ,
    {
      top <- .SD[1L]
      second <- if (.N >= 2L) {
        .SD[2L]
      } else {
        top
      }

      list(
        major_cell_type =
          top$marker_set,
        top_standardized_score =
          top$standardized_marker_score,
        top_raw_score =
          top$raw_marker_score,
        second_marker_set =
          second$marker_set,
        second_standardized_score =
          second$standardized_marker_score,
        score_margin =
          top$standardized_marker_score -
            second$standardized_marker_score
      )
    },
    by = seurat_cluster
  ]

  annotation[
    ,
    annotation_confidence := fcase(
      !is.finite(top_raw_score) |
        top_raw_score <= 0,
      "Low",

      score_margin >= 0.75,
      "High",

      score_margin >= 0.25,
      "Moderate",

      default = "Low"
    )
  ]

  list(
    scores = scores,
    annotation = annotation,
    availability = availability
  )
}

major_score_result <-
  calculate_cluster_marker_scores(
    cardiac,
    major_marker_sets
  )

cluster_score_long <-
  major_score_result$scores
cluster_annotation <-
  major_score_result$annotation
major_marker_availability <-
  major_score_result$availability

fwrite(
  cluster_score_long,
  file.path(
    DIRS$tables,
    "10_major_celltype_cluster_marker_scores.csv"
  )
)
fwrite(
  cluster_annotation,
  file.path(
    DIRS$tables,
    "11_major_celltype_cluster_annotation.csv"
  )
)
fwrite(
  major_marker_availability,
  file.path(
    DIRS$tables,
    "12_major_marker_gene_availability.csv"
  )
)

cluster_to_label <- setNames(
  cluster_annotation$major_cell_type,
  cluster_annotation$seurat_cluster
)
cluster_to_confidence <- setNames(
  cluster_annotation$annotation_confidence,
  cluster_annotation$seurat_cluster
)

cardiac$major_cell_type <-
  unname(
    cluster_to_label[
      as.character(
        cardiac$seurat_clusters
      )
    ]
  )
cardiac$annotation_confidence <-
  unname(
    cluster_to_confidence[
      as.character(
        cardiac$seurat_clusters
      )
    ]
  )

if (
  any(is.na(cardiac$major_cell_type)) ||
  any(!nzchar(cardiac$major_cell_type))
) {
  stop(
    "At least one cardiac cell lacked a major-cell-type annotation."
  )
}

log_msg(
  "Major cardiac cell types annotated: ",
  paste(
    sort(
      unique(cardiac$major_cell_type)
    ),
    collapse = ", "
  )
)

## Cluster markers are generated for review, not used as independent
## biological replicates.
cluster_markers <- tryCatch(
  FindAllMarkers(
    cardiac,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.20,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    max.cells.per.ident =
      MAX_CELLS_PER_CLUSTER_MARKER_TEST,
    random.seed = 20260714,
    verbose = FALSE
  ),
  error = function(e) {
    add_warning(
      "MARKERS",
      "FindAllMarkers",
      conditionMessage(e)
    )
    data.frame()
  }
)
setDT(cluster_markers)

write_csv_safe(
  cluster_markers,
  file.path(
    DIRS$tables,
    "13_cluster_positive_markers.csv.gz"
  ),
  compress = TRUE
)

if (nrow(cluster_markers) > 0L) {
  fc_column <- intersect(
    c(
      "avg_log2FC",
      "avg_logFC"
    ),
    names(cluster_markers)
  )[1L]

  if (
    length(fc_column) == 1L &&
    !is.na(fc_column)
  ) {
    top_cluster_markers <- cluster_markers[
      order(
        cluster,
        p_val_adj,
        -abs(get(fc_column))
      )
    ][
      ,
      head(.SD, TOP_MARKERS_PER_CLUSTER_CHECK),
      by = cluster
    ]
  } else {
    top_cluster_markers <- cluster_markers[
      order(
        cluster,
        p_val_adj
      )
    ][
      ,
      head(.SD, TOP_MARKERS_PER_CLUSTER_CHECK),
      by = cluster
    ]
  }
} else {
  top_cluster_markers <- data.table()
}

write_csv_safe(
  top_cluster_markers,
  file.path(
    DIRS$tables,
    "14_top_cluster_markers_for_review.csv"
  )
)

############################################################
## 9. Cell-type composition and exploratory sample-level tests
############################################################

cell_metadata <- as.data.table(
  cardiac@meta.data,
  keep.rownames = "cell"
)

celltype_counts <- cell_metadata[
  ,
  .(
    cell_count = .N
  ),
  by = .(
    sample_accession,
    condition,
    major_cell_type
  )
]
celltype_counts[
  ,
  total_cells_in_sample := sum(cell_count),
  by = sample_accession
]
celltype_counts[
  ,
  cell_fraction :=
    cell_count / total_cells_in_sample
]

all_samples <- unique(
  cell_metadata[
    ,
    .(
      sample_accession,
      condition
    )
  ]
)
all_celltypes <- sort(
  unique(
    cell_metadata$major_cell_type
  )
)

composition_complete <- CJ(
  sample_accession =
    all_samples$sample_accession,
  major_cell_type =
    all_celltypes,
  unique = TRUE
)
composition_complete <- merge(
  composition_complete,
  all_samples,
  by = "sample_accession",
  all.x = TRUE
)
composition_complete <- merge(
  composition_complete,
  celltype_counts[
    ,
    .(
      sample_accession,
      major_cell_type,
      cell_count,
      total_cells_in_sample,
      cell_fraction
    )
  ],
  by = c(
    "sample_accession",
    "major_cell_type"
  ),
  all.x = TRUE
)
composition_complete[
  is.na(cell_count),
  cell_count := 0L
]
sample_totals <- cell_metadata[
  ,
  .(
    total_cells_in_sample = .N
  ),
  by = sample_accession
]
composition_complete[
  ,
  total_cells_in_sample := NULL
]
composition_complete <- merge(
  composition_complete,
  sample_totals,
  by = "sample_accession",
  all.x = TRUE
)
composition_complete[
  is.na(cell_fraction),
  cell_fraction :=
    cell_count / total_cells_in_sample
]

fwrite(
  composition_complete,
  file.path(
    DIRS$tables,
    "15_celltype_composition_by_sample.csv"
  )
)

composition_stats <- composition_complete[
  ,
  {
    control_values <- cell_fraction[
      condition == "Control"
    ]
    hfpef_values <- cell_fraction[
      condition == "HFpEF"
    ]

    list(
      control_mean_fraction =
        mean(control_values),
      hfpef_mean_fraction =
        mean(hfpef_values),
      hfpef_minus_control =
        mean(hfpef_values) -
          mean(control_values),
      hedges_g_HFpEF_vs_Control =
        hedges_g(
          hfpef_values,
          control_values
        ),
      wilcoxon_p =
        safe_wilcox_p(
          hfpef_values,
          control_values
        )
    )
  },
  by = major_cell_type
]
composition_stats[
  ,
  wilcoxon_fdr := p.adjust(
    wilcoxon_p,
    method = "BH"
  )
]
fwrite(
  composition_stats,
  file.path(
    DIRS$tables,
    "16_celltype_composition_exploratory_statistics.csv"
  )
)

############################################################
## 10. Score Stage 2 signed programs in individual cells
############################################################

match_signature_genes <- function(
  requested_genes,
  available_features
) {
  requested_keys <- gene_key(
    requested_genes
  )
  available_keys <- gene_key(
    available_features
  )

  idx <- match(
    requested_keys,
    available_keys
  )
  idx <- idx[!is.na(idx)]
  unique(
    available_features[idx]
  )
}

add_signed_program_scores <- function(
  object,
  signature_sets
) {
  data_matrix <- get_assay_matrix(
    object,
    layer = "data",
    assay = "RNA"
  )
  available_features <- rownames(
    data_matrix
  )

  coverage_records <- list()

  for (set_name in names(signature_sets)) {
    requested_up <-
      signature_sets[[set_name]]$up
    requested_down <-
      signature_sets[[set_name]]$down

    present_up <- match_signature_genes(
      requested_up,
      available_features
    )
    present_down <- match_signature_genes(
      requested_down,
      available_features
    )

    if (length(present_up) > 0L) {
      up_score <- Matrix::colMeans(
        data_matrix[
          present_up,
          ,
          drop = FALSE
        ]
      )
    } else {
      up_score <- rep(
        NA_real_,
        ncol(data_matrix)
      )
    }

    if (length(present_down) > 0L) {
      down_score <- Matrix::colMeans(
        data_matrix[
          present_down,
          ,
          drop = FALSE
        ]
      )
    } else {
      down_score <- rep(
        NA_real_,
        ncol(data_matrix)
      )
    }

    net_score <- up_score - down_score

    up_column <- paste0(
      "score_",
      set_name,
      "_up"
    )
    down_column <- paste0(
      "score_",
      set_name,
      "_down"
    )
    net_column <- paste0(
      "score_",
      set_name,
      "_net"
    )

    object[[up_column]] <-
      as.numeric(up_score)
    object[[down_column]] <-
      as.numeric(down_score)
    object[[net_column]] <-
      as.numeric(net_score)

    coverage_records[[set_name]] <- data.table(
      signature_name = set_name,
      requested_up_genes =
        length(unique(gene_key(requested_up))),
      detected_up_genes =
        length(present_up),
      requested_down_genes =
        length(unique(gene_key(requested_down))),
      detected_down_genes =
        length(present_down),
      detected_up_symbols =
        paste(present_up, collapse = ";"),
      detected_down_symbols =
        paste(present_down, collapse = ";"),
      net_score_column = net_column
    )
  }

  list(
    object = object,
    coverage = rbindlist(
      coverage_records,
      use.names = TRUE,
      fill = TRUE
    )
  )
}

scored_result <- add_signed_program_scores(
  cardiac,
  signature_sets
)
cardiac <- scored_result$object
signature_coverage <- scored_result$coverage

fwrite(
  signature_coverage,
  file.path(
    DIRS$tables,
    "17_signature_gene_coverage_in_GSE236585.csv"
  )
)

primary_signature_names <- c(
  paste0(
    "Ccr2pos_Top",
    PRIMARY_SIGNATURE_SIZE
  ),
  paste0(
    "Ccr2neg_Top",
    PRIMARY_SIGNATURE_SIZE
  ),
  paste0(
    "CrossSubset_Top",
    PRIMARY_SIGNATURE_SIZE
  )
)
primary_score_columns <- paste0(
  "score_",
  primary_signature_names,
  "_net"
)

if (
  any(!primary_score_columns %in%
        names(cardiac@meta.data))
) {
  stop(
    "One or more primary signature score columns were not generated."
  )
}

score_columns_all <- grep(
  "^score_.*_net$",
  names(cardiac@meta.data),
  value = TRUE
)

cell_score_summary <- as.data.table(
  cardiac@meta.data,
  keep.rownames = "cell"
)[
  ,
  c(
    list(
      cell_count = .N
    ),
    lapply(
      .SD,
      mean,
      na.rm = TRUE
    )
  ),
  by = .(
    sample_accession,
    condition,
    major_cell_type
  ),
  .SDcols = score_columns_all
]
fwrite(
  cell_score_summary,
  file.path(
    DIRS$tables,
    "18_sample_celltype_mean_cell_level_scores.csv"
  )
)

cell_score_long <- melt(
  cell_score_summary,
  id.vars = c(
    "sample_accession",
    "condition",
    "major_cell_type",
    "cell_count"
  ),
  measure.vars = score_columns_all,
  variable.name = "score_column",
  value.name = "mean_cell_score"
)
cell_score_long[
  ,
  signature_name := sub(
    "^score_",
    "",
    sub(
      "_net$",
      "",
      score_column
    )
  )
]

cell_score_stats <- cell_score_long[
  ,
  {
    control_values <-
      mean_cell_score[
        condition == "Control"
      ]
    hfpef_values <-
      mean_cell_score[
        condition == "HFpEF"
      ]

    list(
      control_mean =
        mean(control_values),
      hfpef_mean =
        mean(hfpef_values),
      hfpef_minus_control =
        mean(hfpef_values) -
          mean(control_values),
      hedges_g_HFpEF_vs_Control =
        hedges_g(
          hfpef_values,
          control_values
        ),
      wilcoxon_p =
        safe_wilcox_p(
          hfpef_values,
          control_values
        ),
      control_samples =
        sum(is.finite(control_values)),
      hfpef_samples =
        sum(is.finite(hfpef_values))
    )
  },
  by = .(
    major_cell_type,
    signature_name
  )
]
cell_score_stats[
  ,
  wilcoxon_fdr := p.adjust(
    wilcoxon_p,
    method = "BH"
  ),
  by = signature_name
]
fwrite(
  cell_score_stats,
  file.path(
    DIRS$tables,
    "19_sample_level_cell_score_statistics.csv"
  )
)

############################################################
## 11. Build sample-by-cell-type pseudobulk count matrices
############################################################

aggregate_sparse_counts <- function(
  object,
  group_columns
) {
  object <- join_layers_safe(
    object,
    assay = "RNA"
  )
  counts <- get_assay_matrix(
    object,
    layer = "counts",
    assay = "RNA"
  )

  meta <- as.data.table(
    object@meta.data,
    keep.rownames = "cell"
  )
  meta <- meta[
    match(
      colnames(counts),
      cell
    )
  ]

  if (any(is.na(meta$cell))) {
    stop(
      "Cell metadata could not be aligned to the count matrix."
    )
  }

  group_text <- do.call(
    paste,
    c(
      meta[, ..group_columns],
      sep = "||"
    )
  )
  group_factor <- factor(
    group_text,
    levels = unique(group_text)
  )

  design <- Matrix::sparse.model.matrix(
    ~ 0 + group_factor
  )
  colnames(design) <- levels(
    group_factor
  )

  aggregated <- counts %*% design
  colnames(aggregated) <- levels(
    group_factor
  )

  ## Build grouping metadata explicitly. Combining a character vector of
  ## column names with list(pseudobulk_id = group_text) under with=FALSE is
  ## invalid data.table syntax.
  group_meta <- copy(
    meta[
      ,
      ..group_columns
    ]
  )
  group_meta[, pseudobulk_id := group_text]
  group_meta <- unique(group_meta)

  cell_counts <- data.table(
    pseudobulk_id = levels(
      group_factor
    ),
    cells = as.integer(
      table(group_factor)[
        levels(group_factor)
      ]
    )
  )

  group_meta <- merge(
    group_meta,
    cell_counts,
    by = "pseudobulk_id",
    all.x = TRUE
  )

  group_meta <- group_meta[
    match(
      colnames(aggregated),
      pseudobulk_id
    )
  ]

  list(
    counts = aggregated,
    metadata = group_meta
  )
}

pseudobulk_major <- aggregate_sparse_counts(
  cardiac,
  c(
    "sample_accession",
    "condition",
    "major_cell_type"
  )
)

saveRDS(
  pseudobulk_major$counts,
  file.path(
    DIRS$objects,
    "GSE236585_major_celltype_pseudobulk_counts.rds"
  )
)
fwrite(
  pseudobulk_major$metadata,
  file.path(
    DIRS$tables,
    "20_major_celltype_pseudobulk_metadata.csv"
  )
)

pseudobulk_eligibility <- pseudobulk_major$metadata[
  ,
  .(
    pseudobulk_samples = .N,
    minimum_cells_per_pseudobulk =
      min(cells),
    control_samples = sum(
      condition == "Control" &
        cells >=
          MIN_CELLS_PER_SAMPLE_CELLTYPE
    ),
    hfpef_samples = sum(
      condition == "HFpEF" &
        cells >=
          MIN_CELLS_PER_SAMPLE_CELLTYPE
    )
  ),
  by = major_cell_type
]
pseudobulk_eligibility[
  ,
  eligible := (
    control_samples >=
      MIN_SAMPLES_PER_CONDITION &
      hfpef_samples >=
        MIN_SAMPLES_PER_CONDITION
  )
]
fwrite(
  pseudobulk_eligibility,
  file.path(
    DIRS$tables,
    "21_major_celltype_pseudobulk_eligibility.csv"
  )
)

eligible_celltypes <- pseudobulk_eligibility[
  eligible == TRUE,
  major_cell_type
]

if (length(eligible_celltypes) < 3L) {
  stop(
    "Fewer than three major cell types were eligible for six-sample pseudobulk analysis."
  )
}

############################################################
## 12. Pseudobulk edgeR and limma-voom comparisons
############################################################

run_celltype_pseudobulk <- function(
  cell_type,
  pb_counts,
  pb_meta
) {
  meta_ct <- pb_meta[
    major_cell_type == cell_type &
      cells >=
        MIN_CELLS_PER_SAMPLE_CELLTYPE
  ]
  meta_ct[, condition := factor(
    condition,
    levels = c("Control", "HFpEF")
  )]
  setorder(
    meta_ct,
    condition,
    sample_accession
  )

  if (
    sum(meta_ct$condition == "Control") <
      MIN_SAMPLES_PER_CONDITION ||
    sum(meta_ct$condition == "HFpEF") <
      MIN_SAMPLES_PER_CONDITION
  ) {
    return(NULL)
  }

  counts_ct <- pb_counts[
    ,
    meta_ct$pseudobulk_id,
    drop = FALSE
  ]
  colnames(counts_ct) <-
    meta_ct$sample_accession

  y <- edgeR::DGEList(
    counts = counts_ct,
    samples = as.data.frame(meta_ct)
  )
  keep <- edgeR::filterByExpr(
    y,
    group = meta_ct$condition,
    min.count = PSEUDOBULK_MIN_COUNT
  )
  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  y <- edgeR_norm_lib_sizes(y)

  design <- model.matrix(
    ~ condition,
    data = meta_ct
  )

  y <- edgeR::estimateDisp(
    y,
    design,
    robust = TRUE
  )
  fit_edge <- edgeR::glmQLFit(
    y,
    design,
    robust = TRUE
  )
  test_edge <- edgeR::glmQLFTest(
    fit_edge,
    coef = "conditionHFpEF"
  )
  edge_tab <- edgeR::topTags(
    test_edge,
    n = Inf,
    sort.by = "none"
  )$table
  setDT(
    edge_tab,
    keep.rownames = "feature"
  )
  edge_out <- edge_tab[
    ,
    .(
      feature,
      edgeR_logFC = logFC,
      edgeR_logCPM = logCPM,
      edgeR_F = F,
      edgeR_pvalue = PValue,
      edgeR_padj = FDR
    )
  ]

  voom_object <- limma::voom(
    y,
    design,
    plot = FALSE
  )
  fit_limma <- limma::lmFit(
    voom_object,
    design
  )
  fit_limma <- limma::eBayes(
    fit_limma,
    robust = TRUE,
    trend = FALSE
  )
  limma_tab <- limma::topTable(
    fit_limma,
    coef = "conditionHFpEF",
    number = Inf,
    sort.by = "none"
  )
  setDT(
    limma_tab,
    keep.rownames = "feature"
  )
  limma_out <- limma_tab[
    ,
    .(
      feature,
      limma_logFC = logFC,
      limma_AveExpr = AveExpr,
      limma_t = t,
      limma_pvalue = P.Value,
      limma_padj = adj.P.Val
    )
  ]

  merged <- merge(
    edge_out,
    limma_out,
    by = "feature",
    all = TRUE
  )
  merged[, feature_key := gene_key(feature)]
  merged[, major_cell_type := cell_type]
  merged[, edgeR_limma_sign_agreement := (
    sign(edgeR_logFC) ==
      sign(limma_logFC)
  )]
  merged[, edgeR_formal_fdr_005 := (
    !is.na(edgeR_padj) &
      edgeR_padj <= FORMAL_FDR
  )]
  merged[, edgeR_exploratory_fdr_010 := (
    !is.na(edgeR_padj) &
      edgeR_padj <= EXPLORATORY_FDR
  )]
  merged[, limma_formal_fdr_005 := (
    !is.na(limma_padj) &
      limma_padj <= FORMAL_FDR
  )]
  merged[, limma_exploratory_fdr_010 := (
    !is.na(limma_padj) &
      limma_padj <= EXPLORATORY_FDR
  )]

  summary <- data.table(
    major_cell_type = cell_type,
    pseudobulk_samples = nrow(meta_ct),
    genes_tested = nrow(merged),
    edgeR_fdr_005_up = sum(
      merged$edgeR_padj <=
        FORMAL_FDR &
        merged$edgeR_logFC > 0,
      na.rm = TRUE
    ),
    edgeR_fdr_005_down = sum(
      merged$edgeR_padj <=
        FORMAL_FDR &
        merged$edgeR_logFC < 0,
      na.rm = TRUE
    ),
    edgeR_fdr_010_up = sum(
      merged$edgeR_padj <=
        EXPLORATORY_FDR &
        merged$edgeR_logFC > 0,
      na.rm = TRUE
    ),
    edgeR_fdr_010_down = sum(
      merged$edgeR_padj <=
        EXPLORATORY_FDR &
        merged$edgeR_logFC < 0,
      na.rm = TRUE
    ),
    limma_fdr_005_up = sum(
      merged$limma_padj <=
        FORMAL_FDR &
        merged$limma_logFC > 0,
      na.rm = TRUE
    ),
    limma_fdr_005_down = sum(
      merged$limma_padj <=
        FORMAL_FDR &
        merged$limma_logFC < 0,
      na.rm = TRUE
    ),
    lfc_pearson = safe_pearson(
      merged$edgeR_logFC,
      merged$limma_logFC
    ),
    lfc_spearman = safe_spearman(
      merged$edgeR_logFC,
      merged$limma_logFC
    ),
    sign_agreement = mean(
      merged$edgeR_limma_sign_agreement,
      na.rm = TRUE
    )
  )

  list(
    results = merged,
    summary = summary,
    y = y,
    voom = voom_object,
    metadata = meta_ct
  )
}

pseudobulk_results_list <- list()
pseudobulk_summary_list <- list()
pseudobulk_objects <- list()

for (ct in eligible_celltypes) {
  log_msg(
    "Pseudobulk analysis: ",
    ct
  )

  result_ct <- run_celltype_pseudobulk(
    cell_type = ct,
    pb_counts =
      pseudobulk_major$counts,
    pb_meta =
      pseudobulk_major$metadata
  )

  if (!is.null(result_ct)) {
    pseudobulk_results_list[[ct]] <-
      result_ct$results
    pseudobulk_summary_list[[ct]] <-
      result_ct$summary
    pseudobulk_objects[[ct]] <- list(
      y = result_ct$y,
      voom = result_ct$voom,
      metadata = result_ct$metadata
    )
  }
}

pseudobulk_results <- rbindlist(
  pseudobulk_results_list,
  use.names = TRUE,
  fill = TRUE
)
pseudobulk_summary <- rbindlist(
  pseudobulk_summary_list,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(pseudobulk_results) == 0L) {
  stop(
    "No major-cell-type pseudobulk results were generated."
  )
}

fwrite(
  pseudobulk_results,
  file.path(
    DIRS$tables,
    "22_major_celltype_pseudobulk_DE_edgeR_limma.csv.gz"
  ),
  compress = "gzip"
)
fwrite(
  pseudobulk_summary,
  file.path(
    DIRS$tables,
    "23_major_celltype_pseudobulk_DE_summary.csv"
  )
)
saveRDS(
  pseudobulk_objects,
  file.path(
    DIRS$objects,
    "GSE236585_major_celltype_pseudobulk_model_objects.rds"
  )
)

top_pseudobulk_de <- pseudobulk_results[
  order(
    major_cell_type,
    edgeR_padj,
    -abs(edgeR_logFC)
  )
][
  ,
  head(
    .SD,
    TOP_DE_GENES_PER_CELLTYPE_CHECK
  ),
  by = major_cell_type
]
fwrite(
  top_pseudobulk_de,
  file.path(
    DIRS$tables,
    "24_top_pseudobulk_DE_genes_per_celltype.csv"
  )
)

############################################################
## 13. Sample-level pseudobulk program scores
############################################################

score_pseudobulk_programs <- function(
  pb_counts,
  pb_meta,
  signature_sets
) {
  score_records <- list()
  coverage_records <- list()

  for (ct in unique(pb_meta$major_cell_type)) {
    meta_ct <- pb_meta[
      major_cell_type == ct
    ]
    counts_ct <- pb_counts[
      ,
      meta_ct$pseudobulk_id,
      drop = FALSE
    ]
    colnames(counts_ct) <-
      meta_ct$sample_accession

    y <- edgeR::DGEList(
      counts = counts_ct
    )
    y <- edgeR_norm_lib_sizes(y)
    log_cpm <- edgeR::cpm(
      y,
      log = TRUE,
      prior.count = 2
    )

    feature_keys <- gene_key(
      rownames(log_cpm)
    )

    ## Gene-wise standardization across biological samples within
    ## each cell type.
    gene_z <- t(
      scale(
        t(log_cpm)
      )
    )
    gene_z[
      !is.finite(gene_z)
    ] <- 0

    for (set_name in names(signature_sets)) {
      up_keys <- gene_key(
        signature_sets[[set_name]]$up
      )
      down_keys <- gene_key(
        signature_sets[[set_name]]$down
      )

      up_idx <- which(
        feature_keys %in% up_keys
      )
      down_idx <- which(
        feature_keys %in% down_keys
      )

      raw_up <- if (length(up_idx) > 0L) {
        Matrix::colMeans(
          log_cpm[
            up_idx,
            ,
            drop = FALSE
          ]
        )
      } else {
        rep(NA_real_, ncol(log_cpm))
      }

      raw_down <- if (length(down_idx) > 0L) {
        Matrix::colMeans(
          log_cpm[
            down_idx,
            ,
            drop = FALSE
          ]
        )
      } else {
        rep(NA_real_, ncol(log_cpm))
      }

      z_up <- if (length(up_idx) > 0L) {
        colMeans(
          gene_z[
            up_idx,
            ,
            drop = FALSE
          ]
        )
      } else {
        rep(NA_real_, ncol(gene_z))
      }

      z_down <- if (length(down_idx) > 0L) {
        colMeans(
          gene_z[
            down_idx,
            ,
            drop = FALSE
          ]
        )
      } else {
        rep(NA_real_, ncol(gene_z))
      }

      score_records[[length(score_records) + 1L]] <- data.table(
        sample_accession =
          colnames(log_cpm),
        major_cell_type = ct,
        signature_name = set_name,
        raw_up_mean_logCPM =
          as.numeric(raw_up),
        raw_down_mean_logCPM =
          as.numeric(raw_down),
        raw_net_score =
          as.numeric(raw_up - raw_down),
        z_up_score =
          as.numeric(z_up),
        z_down_score =
          as.numeric(z_down),
        z_net_score =
          as.numeric(z_up - z_down)
      )

      coverage_records[[length(coverage_records) + 1L]] <- data.table(
        major_cell_type = ct,
        signature_name = set_name,
        detected_up_genes =
          length(up_idx),
        detected_down_genes =
          length(down_idx)
      )
    }
  }

  scores <- rbindlist(
    score_records,
    use.names = TRUE,
    fill = TRUE
  )
  coverage <- unique(
    rbindlist(
      coverage_records,
      use.names = TRUE,
      fill = TRUE
    )
  )

  scores <- merge(
    scores,
    unique(
      pb_meta[
        ,
        .(
          sample_accession,
          condition
        )
      ]
    ),
    by = "sample_accession",
    all.x = TRUE
  )

  list(
    scores = scores,
    coverage = coverage
  )
}

pb_score_result <- score_pseudobulk_programs(
  pseudobulk_major$counts,
  pseudobulk_major$metadata,
  signature_sets
)
pseudobulk_program_scores <-
  pb_score_result$scores
pseudobulk_score_coverage <-
  pb_score_result$coverage

fwrite(
  pseudobulk_program_scores,
  file.path(
    DIRS$tables,
    "25_pseudobulk_sample_level_program_scores.csv"
  )
)
fwrite(
  pseudobulk_score_coverage,
  file.path(
    DIRS$tables,
    "26_pseudobulk_program_gene_coverage.csv"
  )
)

pseudobulk_program_stats <-
  pseudobulk_program_scores[
    ,
    {
      control_values <-
        z_net_score[
          condition == "Control"
        ]
      hfpef_values <-
        z_net_score[
          condition == "HFpEF"
        ]

      list(
        control_mean_z_net =
          mean(control_values),
        hfpef_mean_z_net =
          mean(hfpef_values),
        hfpef_minus_control_z_net =
          mean(hfpef_values) -
            mean(control_values),
        hedges_g_HFpEF_vs_Control =
          hedges_g(
            hfpef_values,
            control_values
          ),
        wilcoxon_p =
          safe_wilcox_p(
            hfpef_values,
            control_values
          ),
        control_samples =
          sum(is.finite(control_values)),
        hfpef_samples =
          sum(is.finite(hfpef_values))
      )
    },
    by = .(
      major_cell_type,
      signature_name
    )
  ]
pseudobulk_program_stats[
  ,
  wilcoxon_fdr := p.adjust(
    wilcoxon_p,
    method = "BH"
  ),
  by = signature_name
]
fwrite(
  pseudobulk_program_stats,
  file.path(
    DIRS$tables,
    "27_pseudobulk_program_statistics.csv"
  )
)

############################################################
## 14. Stage 2 to Stage 3 gene-level concordance
############################################################

prepare_stage2_concordance_reference <- function(
  base_table,
  subset_label
) {
  out <- base_table[
    !is.na(symbol) &
      nzchar(symbol) &
      is.finite(disease_lfc)
  ]
  out[, feature_key := gene_key(symbol)]
  setorder(
    out,
    within_subset_rank
  )
  out <- out[, .SD[1L], by = feature_key]
  out[
    ,
    .(
      feature_key,
      subset_source = subset_label,
      stage2_symbol = symbol,
      stage2_disease_lfc = disease_lfc,
      stage2_drug_lfc = drug_lfc,
      stage2_tier = opposition_tier,
      stage2_rank =
        within_subset_rank
    )
  ]
}

stage2_concordance_reference <- rbindlist(
  list(
    prepare_stage2_concordance_reference(
      stage2_pos,
      "Ccr2_positive"
    ),
    prepare_stage2_concordance_reference(
      stage2_neg,
      "Ccr2_negative"
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

concordance_records <- list()
for (ct in unique(
  pseudobulk_results$major_cell_type
)) {
  de_ct <- pseudobulk_results[
    major_cell_type == ct,
    .(
      feature_key,
      stage3_edgeR_logFC =
        edgeR_logFC,
      stage3_edgeR_padj =
        edgeR_padj,
      stage3_limma_logFC =
        limma_logFC,
      stage3_limma_padj =
        limma_padj
    )
  ]

  for (subset_label in c(
    "Ccr2_positive",
    "Ccr2_negative"
  )) {
    ref <- stage2_concordance_reference[
      subset_source == subset_label
    ]
    merged <- merge(
      ref,
      de_ct,
      by = "feature_key",
      all = FALSE
    )

    tier_abc <- merged[
      stage2_tier %in% SIGNATURE_TIERS
    ]

    concordance_records[[
      paste(ct, subset_label, sep = "__")
    ]] <- data.table(
      major_cell_type = ct,
      stage2_subset = subset_label,
      all_overlap_genes = nrow(merged),
      tier_ABC_overlap_genes =
        nrow(tier_abc),
      all_gene_spearman =
        safe_spearman(
          merged$stage2_disease_lfc,
          merged$stage3_edgeR_logFC
        ),
      all_gene_sign_agreement =
        mean(
          sign(
            merged$stage2_disease_lfc
          ) ==
            sign(
              merged$stage3_edgeR_logFC
            ),
          na.rm = TRUE
        ),
      tier_ABC_spearman =
        safe_spearman(
          tier_abc$stage2_disease_lfc,
          tier_abc$stage3_edgeR_logFC
        ),
      tier_ABC_sign_agreement =
        mean(
          sign(
            tier_abc$stage2_disease_lfc
          ) ==
            sign(
              tier_abc$stage3_edgeR_logFC
            ),
          na.rm = TRUE
        ),
      median_stage3_logFC_stage2_disease_up =
        median(
          tier_abc[
            stage2_disease_lfc > 0,
            stage3_edgeR_logFC
          ],
          na.rm = TRUE
        ),
      median_stage3_logFC_stage2_disease_down =
        median(
          tier_abc[
            stage2_disease_lfc < 0,
            stage3_edgeR_logFC
          ],
          na.rm = TRUE
        ),
      orientation_difference =
        median(
          tier_abc[
            stage2_disease_lfc > 0,
            stage3_edgeR_logFC
          ],
          na.rm = TRUE
        ) -
          median(
            tier_abc[
              stage2_disease_lfc < 0,
              stage3_edgeR_logFC
            ],
            na.rm = TRUE
          )
    )
  }
}

stage2_stage3_concordance <- rbindlist(
  concordance_records,
  use.names = TRUE,
  fill = TRUE
)
fwrite(
  stage2_stage3_concordance,
  file.path(
    DIRS$tables,
    "28_stage2_stage3_gene_level_concordance.csv"
  )
)

############################################################
## 15. Macrophage/monocyte subclustering
############################################################

macrophage_cells <- rownames(
  cardiac@meta.data
)[
  cardiac$major_cell_type ==
    "Macrophage_Monocyte"
]

macrophage_status <- "NOT_RUN"
macrophage_cluster_annotation <-
  data.table()
macrophage_state_composition <-
  data.table()
macrophage_state_scores <-
  data.table()
macrophage_state_score_stats <-
  data.table()
macrophage_markers <- data.table()

if (length(macrophage_cells) >= 200L) {
  log_msg(
    "Reclustering ",
    length(macrophage_cells),
    " macrophage/monocyte cells."
  )

  macrophage <- subset(
    cardiac,
    cells = macrophage_cells
  )
  macrophage <- join_layers_safe(
    macrophage,
    assay = "RNA"
  )
  macrophage <- NormalizeData(
    macrophage,
    verbose = FALSE
  )
  macrophage <- FindVariableFeatures(
    macrophage,
    nfeatures = min(
      2500L,
      nrow(macrophage)
    ),
    verbose = FALSE
  )
  macrophage <- ScaleData(
    macrophage,
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
  macrophage <- RunPCA(
    macrophage,
    npcs = min(
      30L,
      N_PCS
    ),
    verbose = FALSE
  )

  macrophage_dims <- seq_len(
    min(
      20L,
      ncol(
        Embeddings(
          macrophage,
          reduction = "pca"
        )
      )
    )
  )

  macrophage <- FindNeighbors(
    macrophage,
    dims = macrophage_dims,
    verbose = FALSE
  )
  macrophage <- FindClusters(
    macrophage,
    resolution =
      MACROPHAGE_CLUSTER_RESOLUTION,
    algorithm = 1L,
    random.seed = 20260714,
    verbose = FALSE
  )
  macrophage <- RunUMAP(
    macrophage,
    dims = macrophage_dims,
    seed.use = 20260714,
    verbose = FALSE
  )

  macrophage_marker_sets <- list(
    Resident_Timd4_Lyve1 = c(
      "Timd4", "Lyve1", "Folr2",
      "Mrc1", "Cd163", "Vsig4",
      "C1qa", "C1qb", "C1qc"
    ),
    Ccr2_Monocyte_like = c(
      "Ccr2", "Ly6c2", "Plac8",
      "Chil3", "Ctss", "Lgals3"
    ),
    Inflammatory_Il1b = c(
      "Il1b", "Tnf", "Nfkbia",
      "Ccl2", "Ccl3", "Cxcl2",
      "S100a8", "S100a9"
    ),
    Spp1_Trem2_Remodeling = c(
      "Spp1", "Trem2", "Gpnmb",
      "Fabp5", "Lpl", "Apoe",
      "Ctsb", "Ctsd"
    ),
    Interferon_responsive = c(
      "Isg15", "Ifit1", "Ifit2",
      "Ifit3", "Irf7", "Rsad2",
      "Oas1a", "Stat1"
    ),
    Antigen_presentation = c(
      "H2-Ab1", "H2-Aa", "Cd74",
      "Ciita", "H2-Eb1"
    ),
    Cycling = c(
      "Mki67", "Top2a", "Stmn1",
      "Tubb5", "Hmgb2"
    )
  )

  macrophage_score_result <-
    calculate_cluster_marker_scores(
      macrophage,
      macrophage_marker_sets
    )

  macrophage_cluster_annotation <-
    copy(
      macrophage_score_result$annotation
    )
  setnames(
    macrophage_cluster_annotation,
    "major_cell_type",
    "macrophage_state"
  )

  macrophage_cluster_to_state <- setNames(
    macrophage_cluster_annotation$macrophage_state,
    macrophage_cluster_annotation$seurat_cluster
  )
  macrophage$macrophage_state <-
    unname(
      macrophage_cluster_to_state[
        as.character(
          macrophage$seurat_clusters
        )
      ]
    )

  fwrite(
    macrophage_score_result$scores,
    file.path(
      DIRS$tables,
      "29_macrophage_cluster_marker_scores.csv"
    )
  )
  fwrite(
    macrophage_cluster_annotation,
    file.path(
      DIRS$tables,
      "30_macrophage_cluster_annotation.csv"
    )
  )
  fwrite(
    macrophage_score_result$availability,
    file.path(
      DIRS$tables,
      "31_macrophage_marker_gene_availability.csv"
    )
  )

  macrophage_markers <- tryCatch(
    FindAllMarkers(
      macrophage,
      assay = "RNA",
      only.pos = TRUE,
      min.pct = 0.20,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      max.cells.per.ident =
        MAX_CELLS_PER_CLUSTER_MARKER_TEST,
      random.seed = 20260714,
      verbose = FALSE
    ),
    error = function(e) {
      add_warning(
        "MACROPHAGE_MARKERS",
        "FindAllMarkers",
        conditionMessage(e)
      )
      data.frame()
    }
  )
  setDT(macrophage_markers)
  write_csv_safe(
    macrophage_markers,
    file.path(
      DIRS$tables,
      "32_macrophage_cluster_positive_markers.csv.gz"
    ),
    compress = TRUE
  )

  macrophage_meta <- as.data.table(
    macrophage@meta.data,
    keep.rownames = "cell"
  )

  macrophage_state_composition <- macrophage_meta[
    ,
    .(
      cell_count = .N
    ),
    by = .(
      sample_accession,
      condition,
      macrophage_state
    )
  ]
  macrophage_state_composition[
    ,
    macrophage_total_cells :=
      sum(cell_count),
    by = sample_accession
  ]
  macrophage_state_composition[
    ,
    state_fraction :=
      cell_count /
        macrophage_total_cells
  ]

  fwrite(
    macrophage_state_composition,
    file.path(
      DIRS$tables,
      "33_macrophage_state_composition_by_sample.csv"
    )
  )

  macrophage_score_columns <- intersect(
    score_columns_all,
    names(macrophage_meta)
  )

  macrophage_state_scores <- macrophage_meta[
    ,
    c(
      list(
        cell_count = .N
      ),
      lapply(
        .SD,
        mean,
        na.rm = TRUE
      )
    ),
    by = .(
      sample_accession,
      condition,
      macrophage_state
    ),
    .SDcols = macrophage_score_columns
  ]

  fwrite(
    macrophage_state_scores,
    file.path(
      DIRS$tables,
      "34_macrophage_state_sample_level_scores.csv"
    )
  )

  macrophage_score_long <- melt(
    macrophage_state_scores,
    id.vars = c(
      "sample_accession",
      "condition",
      "macrophage_state",
      "cell_count"
    ),
    measure.vars =
      macrophage_score_columns,
    variable.name = "score_column",
    value.name = "mean_cell_score"
  )
  macrophage_score_long[
    ,
    signature_name := sub(
      "^score_",
      "",
      sub(
        "_net$",
        "",
        score_column
      )
    )
  ]

  macrophage_state_score_stats <-
    macrophage_score_long[
      ,
      {
        control_values <-
          mean_cell_score[
            condition == "Control"
          ]
        hfpef_values <-
          mean_cell_score[
            condition == "HFpEF"
          ]

        list(
          control_mean =
            mean(control_values),
          hfpef_mean =
            mean(hfpef_values),
          hfpef_minus_control =
            mean(hfpef_values) -
              mean(control_values),
          hedges_g_HFpEF_vs_Control =
            hedges_g(
              hfpef_values,
              control_values
            ),
          wilcoxon_p =
            safe_wilcox_p(
              hfpef_values,
              control_values
            ),
          control_samples =
            sum(is.finite(control_values)),
          hfpef_samples =
            sum(is.finite(hfpef_values))
        )
      },
      by = .(
        macrophage_state,
        signature_name
      )
    ]
  macrophage_state_score_stats[
    ,
    wilcoxon_fdr := p.adjust(
      wilcoxon_p,
      method = "BH"
    ),
    by = signature_name
  ]

  fwrite(
    macrophage_state_score_stats,
    file.path(
      DIRS$tables,
      "35_macrophage_state_score_statistics.csv"
    )
  )

  saveRDS(
    macrophage,
    file.path(
      DIRS$objects,
      "GSE236585_macrophage_subclustered_seurat.rds"
    )
  )

  macrophage_status <- "COMPLETED"
} else {
  add_warning(
    "MACROPHAGE",
    "cell_count",
    paste0(
      "Only ",
      length(macrophage_cells),
      " macrophage/monocyte cells were annotated; subclustering was skipped."
    )
  )
  macrophage_status <- "INSUFFICIENT_CELLS"
}

############################################################
## 16. Save main Seurat object
############################################################

saveRDS(
  cardiac,
  file.path(
    DIRS$objects,
    "GSE236585_stage3_annotated_projected_seurat.rds"
  )
)

############################################################
## 17. Figures and source data
############################################################

major_palette <- c(
  "Cardiomyocyte" = "#8C564B",
  "Fibroblast" = "#E377C2",
  "Endothelial" = "#17BECF",
  "Lymphatic_endothelial" = "#9EDAE5",
  "Pericyte" = "#9467BD",
  "Smooth_muscle" = "#7F7F7F",
  "Macrophage_Monocyte" = "#D62728",
  "Dendritic_cell" = "#FF9896",
  "Neutrophil" = "#BCBD22",
  "T_NK" = "#1F77B4",
  "B_cell" = "#AEC7E8",
  "Mast_cell" = "#FF7F0E",
  "Epicardial_Mesothelial" = "#2CA02C"
)
condition_palette <- c(
  "Control" = "#4C78A8",
  "HFpEF" = "#E45756"
)

umap_coordinates <- as.data.table(
  Embeddings(
    cardiac,
    reduction = "umap"
  ),
  keep.rownames = "cell"
)
umap_coordinates <- merge(
  umap_coordinates,
  as.data.table(
    cardiac@meta.data,
    keep.rownames = "cell"
  )[
    ,
    c(
      "cell",
      "sample_accession",
      "condition",
      "seurat_clusters",
      "major_cell_type",
      primary_score_columns
    ),
    with = FALSE
  ],
  by = "cell",
  all.x = TRUE
)
fwrite(
  umap_coordinates,
  file.path(
    DIRS$source,
    "Fig3A_3B_3D_UMAP_source.csv.gz"
  ),
  compress = "gzip"
)

p_umap_celltype <- DimPlot(
  cardiac,
  reduction = "umap",
  group.by = "major_cell_type",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  scale_color_manual(
    values = major_palette,
    na.value = "grey70"
  ) +
  labs(
    title = "GSE236585 cardiac scRNA-seq",
    subtitle = "Major cell types assigned by cluster-level canonical-marker scoring",
    color = "Major cell type"
  ) +
  theme_bw(base_size = 10)

save_plot_bundle(
  p_umap_celltype,
  "Fig3A_GSE236585_UMAP_major_cell_types",
  9,
  7
)

p_umap_condition <- DimPlot(
  cardiac,
  reduction = "umap",
  group.by = "condition",
  split.by = "condition",
  raster = TRUE
) +
  scale_color_manual(
    values = condition_palette
  ) +
  labs(
    title = "GSE236585 cells by biological condition",
    color = "Condition"
  ) +
  theme_bw(base_size = 10)

save_plot_bundle(
  p_umap_condition,
  "Fig3B_GSE236585_UMAP_by_condition",
  11,
  5.5
)

fwrite(
  composition_complete,
  file.path(
    DIRS$source,
    "Fig3C_celltype_composition_source.csv"
  )
)
p_composition <- ggplot(
  composition_complete,
  aes(
    x = sample_accession,
    y = cell_fraction,
    fill = major_cell_type
  )
) +
  geom_col(
    width = 0.78
  ) +
  facet_grid(
    ~ condition,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = major_palette,
    na.value = "grey70"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    x = NULL,
    y = "Fraction of retained cells",
    fill = "Major cell type",
    title = "Sample-level cardiac cell-type composition",
    subtitle = "Composition is descriptive; inferential comparisons use biological-sample fractions"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot_bundle(
  p_composition,
  "Fig3C_GSE236585_celltype_composition",
  11,
  6
)

primary_pos_col <- paste0(
  "score_Ccr2pos_Top",
  PRIMARY_SIGNATURE_SIZE,
  "_net"
)
primary_neg_col <- paste0(
  "score_Ccr2neg_Top",
  PRIMARY_SIGNATURE_SIZE,
  "_net"
)

p_score_umap_pos <- FeaturePlot(
  cardiac,
  features = primary_pos_col,
  reduction = "umap",
  raster = TRUE
) +
  labs(
    title = paste0(
      "Projected Ccr2-positive disease-like program (Top",
      PRIMARY_SIGNATURE_SIZE,
      ")"
    ),
    subtitle = "Higher score indicates expression aligned with the Stage 2 disease-associated direction"
  ) +
  theme_bw(base_size = 10)

save_plot_bundle(
  p_score_umap_pos,
  "Fig3D_GSE236585_Ccr2pos_program_UMAP",
  8,
  6.5
)

p_score_umap_neg <- FeaturePlot(
  cardiac,
  features = primary_neg_col,
  reduction = "umap",
  raster = TRUE
) +
  labs(
    title = paste0(
      "Projected Ccr2-negative disease-like program (Top",
      PRIMARY_SIGNATURE_SIZE,
      ")"
    ),
    subtitle = "Higher score indicates expression aligned with the Stage 2 disease-associated direction"
  ) +
  theme_bw(base_size = 10)

save_plot_bundle(
  p_score_umap_neg,
  "FigS3A_GSE236585_Ccr2neg_program_UMAP",
  8,
  6.5
)

primary_program_stats <- pseudobulk_program_stats[
  signature_name %in%
    primary_signature_names
]
fwrite(
  primary_program_stats,
  file.path(
    DIRS$source,
    "Fig3E_program_localization_source.csv"
  )
)

p_program_localization <- ggplot(
  primary_program_stats,
  aes(
    x = hfpef_minus_control_z_net,
    y = reorder(
      major_cell_type,
      hfpef_minus_control_z_net
    ),
    shape = signature_name
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ signature_name,
    scales = "free_y"
  ) +
  labs(
    x = "HFpEF − Control sample-level pseudobulk program score",
    y = NULL,
    shape = "Program",
    title = "Cell-type localization of Stage 2 disease-like programs",
    subtitle = "Positive values indicate higher disease-like program activity in HFpEF biological samples"
  ) +
  theme_bw(base_size = 10)

save_plot_bundle(
  p_program_localization,
  "Fig3E_stage2_program_localization_by_celltype",
  12,
  7
)

fwrite(
  stage2_stage3_concordance,
  file.path(
    DIRS$source,
    "Fig3F_stage2_stage3_concordance_source.csv"
  )
)
p_concordance <- ggplot(
  stage2_stage3_concordance,
  aes(
    x = tier_ABC_spearman,
    y = reorder(
      major_cell_type,
      tier_ABC_spearman
    ),
    shape = stage2_subset
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ stage2_subset
  ) +
  labs(
    x = "Spearman correlation:\nStage 2 disease log2FC vs Stage 3 HFpEF pseudobulk log2FC",
    y = NULL,
    shape = "Stage 2 subset",
    title = "Gene-level cross-dataset concordance by cardiac cell type"
  ) +
  theme_bw(base_size = 10)

save_plot_bundle(
  p_concordance,
  "Fig3F_stage2_stage3_gene_level_concordance",
  11,
  7
)

if (exists("macrophage") &&
    macrophage_status == "COMPLETED") {
  macrophage_umap <- as.data.table(
    Embeddings(
      macrophage,
      reduction = "umap"
    ),
    keep.rownames = "cell"
  )
  macrophage_umap <- merge(
    macrophage_umap,
    as.data.table(
      macrophage@meta.data,
      keep.rownames = "cell"
    )[
      ,
      .(
        cell,
        sample_accession,
        condition,
        seurat_clusters,
        macrophage_state
      )
    ],
    by = "cell",
    all.x = TRUE
  )
  fwrite(
    macrophage_umap,
    file.path(
      DIRS$source,
      "Fig3G_macrophage_UMAP_source.csv"
    )
  )

  p_macrophage_umap <- DimPlot(
    macrophage,
    reduction = "umap",
    group.by = "macrophage_state",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) +
    labs(
      title = "GSE236585 macrophage/monocyte subclustering",
      subtitle = "State labels are marker-defined discovery candidates",
      color = "Macrophage state"
    ) +
    theme_bw(base_size = 10)

  save_plot_bundle(
    p_macrophage_umap,
    "Fig3G_GSE236585_macrophage_state_UMAP",
    9,
    7
  )

  primary_macrophage_scores <-
    macrophage_state_scores[
      ,
      c(
        "sample_accession",
        "condition",
        "macrophage_state",
        "cell_count",
        primary_pos_col,
        primary_neg_col
      ),
      with = FALSE
    ]
  primary_macrophage_long <- melt(
    primary_macrophage_scores,
    id.vars = c(
      "sample_accession",
      "condition",
      "macrophage_state",
      "cell_count"
    ),
    measure.vars = c(
      primary_pos_col,
      primary_neg_col
    ),
    variable.name = "program",
    value.name = "sample_mean_score"
  )
  fwrite(
    primary_macrophage_long,
    file.path(
      DIRS$source,
      "Fig3H_macrophage_state_program_source.csv"
    )
  )

  p_macrophage_scores <- ggplot(
    primary_macrophage_long,
    aes(
      x = condition,
      y = sample_mean_score,
      shape = condition
    )
  ) +
    geom_point(
      size = 2.8,
      position = position_jitter(
        width = 0.08,
        height = 0
      )
    ) +
    facet_grid(
      macrophage_state ~ program,
      scales = "free_y"
    ) +
    scale_shape_manual(
      values = c(
        "Control" = 16,
        "HFpEF" = 17
      )
    ) +
    labs(
      x = NULL,
      y = "Sample-level mean disease-like score",
      shape = "Condition",
      title = "Stage 2 programs across macrophage-state candidates",
      subtitle = "Each point is one biological sample"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      )
    )

  save_plot_bundle(
    p_macrophage_scores,
    "Fig3H_macrophage_state_sample_level_program_scores",
    12,
    11
  )
}

## QC figures.
qc_plot_data <- as.data.table(
  cardiac@meta.data,
  keep.rownames = "cell"
)[
  ,
  .(
    cell,
    sample_accession,
    condition,
    nFeature_RNA,
    nCount_RNA,
    percent.mt
  )
]
fwrite(
  qc_plot_data,
  file.path(
    DIRS$source,
    "FigS3B_QC_violin_source.csv.gz"
  ),
  compress = "gzip"
)

p_qc <- VlnPlot(
  cardiac,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample_accession",
  pt.size = 0,
  ncol = 3
) &
  theme_bw(base_size = 8) &
  theme(
    axis.text.x = element_text(
      angle = 50,
      hjust = 1
    )
  )

save_plot_bundle(
  p_qc,
  "FigS3B_GSE236585_QC_violin",
  14,
  5
)

## Major-cell-type marker score heatmap.
marker_heat <- dcast(
  cluster_score_long,
  marker_set ~ seurat_cluster,
  value.var = "standardized_marker_score"
)
marker_heat_matrix <- as.matrix(
  marker_heat[
    ,
    -1,
    with = FALSE
  ]
)
rownames(marker_heat_matrix) <-
  marker_heat$marker_set
write_csv_safe(
  marker_heat,
  file.path(
    DIRS$source,
    "FigS3C_major_marker_score_heatmap_source.csv"
  )
)

png(
  file.path(
    DIRS$figures,
    "FigS3C_major_marker_score_heatmap.png"
  ),
  width = 3000,
  height = 2400,
  res = 300
)
pheatmap(
  marker_heat_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "Cluster-level canonical-marker scores"
)
dev.off()

pdf(
  file.path(
    DIRS$figures,
    "FigS3C_major_marker_score_heatmap.pdf"
  ),
  width = 10,
  height = 8
)
pheatmap(
  marker_heat_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "Cluster-level canonical-marker scores"
)
dev.off()

############################################################
## 18. Workbook, methods, and parameter documentation
############################################################

workbook_path <- file.path(
  DIRS$tables,
  "36_GSE236585_stage3_key_results.xlsx"
)
wb <- createWorkbook()

write_sheet_safe(
  wb,
  "Sample_metadata",
  sample_meta
)
write_sheet_safe(
  wb,
  "10x_mapping",
  file_map
)
write_sheet_safe(
  wb,
  "QC_thresholds",
  qc_thresholds
)
write_sheet_safe(
  wb,
  "QC_summary",
  qc_summary
)
write_sheet_safe(
  wb,
  "Cluster_annotation",
  cluster_annotation
)
write_sheet_safe(
  wb,
  "Celltype_composition",
  composition_complete
)
write_sheet_safe(
  wb,
  "Composition_stats",
  composition_stats
)
write_sheet_safe(
  wb,
  "Signature_sizes",
  signature_size_summary
)
write_sheet_safe(
  wb,
  "Signature_coverage",
  signature_coverage
)
write_sheet_safe(
  wb,
  "Program_stats",
  pseudobulk_program_stats
)
write_sheet_safe(
  wb,
  "PB_eligibility",
  pseudobulk_eligibility
)
write_sheet_safe(
  wb,
  "PB_DE_summary",
  pseudobulk_summary
)
write_sheet_safe(
  wb,
  "Stage2_Stage3_concordance",
  stage2_stage3_concordance
)
write_sheet_safe(
  wb,
  "Macrophage_annotation",
  macrophage_cluster_annotation
)
write_sheet_safe(
  wb,
  "Macrophage_score_stats",
  macrophage_state_score_stats
)

saveWorkbook(
  wb,
  workbook_path,
  overwrite = TRUE
)

parameter_table <- data.table(
  parameter = c(
    "Random seed",
    "Minimum features per cell",
    "Maximum features hard ceiling",
    "Maximum counts hard ceiling",
    "Maximum mitochondrial percentage hard ceiling",
    "Adaptive QC multiplier",
    "Adaptive QC upper quantile",
    "Doublet method",
    "Variable features",
    "PCA dimensions",
    "Major clustering resolution",
    "Macrophage clustering resolution",
    "Minimum cells per sample-cell type",
    "Minimum samples per condition",
    "Pseudobulk primary method",
    "Pseudobulk sensitivity method",
    "Signature tiers",
    "Signature sizes",
    "Primary signature size",
    "Inferential unit"
  ),
  value = c(
    "20260714",
    as.character(MIN_FEATURES_HARD),
    as.character(MAX_FEATURES_HARD),
    as.character(MAX_COUNTS_HARD),
    as.character(MAX_PERCENT_MT_HARD),
    as.character(QC_MAD_MULTIPLIER),
    as.character(QC_UPPER_QUANTILE),
    doublet_status,
    as.character(N_VARIABLE_FEATURES),
    paste(dims_use, collapse = ","),
    as.character(MAJOR_CLUSTER_RESOLUTION),
    as.character(MACROPHAGE_CLUSTER_RESOLUTION),
    as.character(MIN_CELLS_PER_SAMPLE_CELLTYPE),
    as.character(MIN_SAMPLES_PER_CONDITION),
    "edgeR quasi-likelihood pseudobulk",
    "limma-voom pseudobulk",
    paste(SIGNATURE_TIERS, collapse = "; "),
    paste(SIGNATURE_SIZES, collapse = ", "),
    as.character(PRIMARY_SIGNATURE_SIZE),
    "Biological sample"
  ),
  rationale = c(
    "Reproducibility",
    "Remove low-complexity droplets",
    "Prevent extreme high-feature cells while retaining cardiac cell diversity",
    "Prevent extreme high-count cells",
    "Remove strongly mitochondrial cells",
    "Sample-specific outlier detection",
    "Sample-specific outlier detection",
    "Optional per-sample doublet classification; nonfatal if unavailable",
    "Major cell-state representation",
    "Graph construction and UMAP",
    "Major cardiac cell-type discovery",
    "Macrophage-state discovery",
    "Ensure adequate pseudobulk support",
    "Require all three biological samples per condition",
    "Primary count-based sample-level inference",
    "Independent mean-variance modeling sensitivity analysis",
    "Exclude weak directional-only Tier D genes from projected signatures",
    "Signature-size sensitivity analysis",
    "Primary visualization and localization score",
    "Avoid cell-level pseudoreplication"
  )
)
fwrite(
  parameter_table,
  file.path(
    DIRS$methods,
    "stage3_parameters_and_rationale.csv"
  )
)

methods_text <- c(
  "HFpEF Stage 3: GSE236585 cardiac scRNA-seq discovery and projection",
  "",
  "Input and biological design:",
  "- Six independent ventricular scRNA-seq samples: three HFpEF and three control mice.",
  "- Sample identities and conditions were taken from the locked Stage 1 manifest.",
  "- Eighteen 10x files were mapped by GSM accession and validated by matrix-feature-barcode dimensions.",
  "",
  "Quality control:",
  "- Cells with fewer than 200 detected genes were removed.",
  "- Sample-specific upper thresholds used the larger of the 99.5th percentile and median + 4 MAD, constrained by hard ceilings.",
  "- Hard ceilings were 9,000 genes, 120,000 counts, and 25% mitochondrial reads.",
  "- scDblFinder was applied by biological sample when available; failure or package absence was recorded and did not alter other analyses.",
  "",
  "Clustering and annotation:",
  "- LogNormalize, 3,000 variable genes, percent-mitochondrial regression, PCA, graph clustering, and UMAP were used.",
  "- Major cell types were assigned at the cluster level from prespecified canonical cardiac marker sets.",
  "- Annotation confidence and competing marker scores were retained for manual review.",
  "- Cluster marker tests are descriptive and do not constitute biological replication.",
  "",
  "Stage 2 program projection:",
  "- Stage 2 Tier A-C genes supported by opposite disease/drug directions and DESeq2-edgeR directional agreement were used.",
  "- Separate disease-up/drug-down and disease-down/drug-up components were constructed.",
  "- Top50, Top100, Top150, and Top200 versions were evaluated for Ccr2-positive, Ccr2-negative, and cross-subset programs.",
  "- The signed disease-like score was mean expression of disease-up genes minus mean expression of disease-down genes.",
  "- No gene, including Nfkb1, was forced into a signature.",
  "",
  "Pseudobulk inference:",
  "- Raw counts were summed within biological sample and major cell type.",
  "- Cell types required at least 20 cells in each of all three samples per condition.",
  "- edgeR quasi-likelihood was the primary differential-expression method.",
  "- limma-voom was used as an independent sensitivity method.",
  "- HFpEF versus control was the tested contrast.",
  "",
  "Sample-level program testing:",
  "- Pseudobulk logCPM was standardized gene-wise across the six biological samples within each cell type.",
  "- Program scores and HFpEF-control effect sizes were calculated at the sample level.",
  "- Wilcoxon P values and Hedges g were reported as exploratory measures because n = 3 per condition.",
  "",
  "Macrophage analysis:",
  "- Cells annotated as macrophage/monocyte were reclustered.",
  "- Resident, Ccr2/monocyte-like, inflammatory, Spp1/Trem2 remodeling, interferon-responsive, antigen-presentation, and cycling marker programs were compared.",
  "- State labels are discovery candidates requiring external and experimental validation.",
  "",
  "Claim boundary:",
  "- Stage 3 localizes and tests concordance of Stage 2 programs in an independent cardiac dataset.",
  "- It does not prove dapagliflozin exposure in GSE236585, direct drug targeting, transcription-factor causality, or cell-cell communication."
)
writeLines(
  methods_text,
  file.path(
    DIRS$methods,
    "stage3_methods_and_claim_boundaries.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    DIRS$methods,
    "sessionInfo.txt"
  )
)

############################################################
## 19. Completion checks and run status
############################################################

warnings_dt <- if (
  length(warning_records) > 0L
) {
  rbindlist(
    warning_records,
    use.names = TRUE,
    fill = TRUE
  )
} else {
  data.table(
    timestamp = character(),
    category = character(),
    item = character(),
    message = character()
  )
}
fwrite(
  warnings_dt,
  file.path(
    DIRS$tables,
    "37_warnings_and_nonfatal_issues.csv"
  )
)

primary_coverage_check <- signature_coverage[
  signature_name %in%
    primary_signature_names
]

scientific_checks <- data.table(
  check = c(
    "Locked biological samples",
    "Control samples",
    "HFpEF samples",
    "Valid 10x triples",
    "Samples retaining >=100 cells",
    "Major cell types annotated",
    "Eligible pseudobulk cell types",
    "Pseudobulk DE results",
    "Primary signature sets",
    "Primary signatures with >=10 detected genes in each direction",
    "Stage2-Stage3 concordance rows",
    "Macrophage cells",
    "Macrophage subclustering"
  ),
  observed = c(
    nrow(sample_meta),
    sum(sample_meta$condition == "Control"),
    sum(sample_meta$condition == "HFpEF"),
    sum(file_map$dimensions_valid),
    sum(qc_summary$cells_after_qc >= 100L),
    uniqueN(cardiac$major_cell_type),
    length(eligible_celltypes),
    nrow(pseudobulk_results),
    length(primary_signature_names),
    sum(
      primary_coverage_check$detected_up_genes >= 10L &
        primary_coverage_check$detected_down_genes >= 10L
    ),
    nrow(stage2_stage3_concordance),
    length(macrophage_cells),
    as.integer(
      macrophage_status == "COMPLETED"
    )
  ),
  expected = c(
    6L,
    3L,
    3L,
    6L,
    6L,
    5L,
    3L,
    1L,
    3L,
    3L,
    1L,
    200L,
    1L
  ),
  comparison = c(
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "at_least",
    "at_least",
    "at_least",
    "equal",
    "equal",
    "at_least",
    "at_least",
    "equal"
  )
)
scientific_checks[
  ,
  status := fcase(
    comparison == "equal" &
      observed == expected,
    "PASS",

    comparison == "at_least" &
      observed >= expected,
    "PASS",

    default = "FAIL"
  )
]

fwrite(
  scientific_checks,
  file.path(
    DIRS$tables,
    "38_scientific_completion_checks.csv"
  )
)

script_copy_status <- "NOT_DETECTED"
if (
  length(SCRIPT_FILE) == 1L &&
  !is.na(SCRIPT_FILE) &&
  file.exists(SCRIPT_FILE)
) {
  methods_script <- file.path(
    DIRS$methods,
    "HFpEF_Stage3_GSE236585_scRNA_Projection_FIXED_v2_REPLACE.R"
  )
  check_script <- file.path(
    DIRS$check,
    "HFpEF_Stage3_GSE236585_scRNA_Projection_FIXED_v2_REPLACE.R"
  )

  copy_methods <- file.copy(
    SCRIPT_FILE,
    methods_script,
    overwrite = TRUE
  )
  copy_check <- file.copy(
    SCRIPT_FILE,
    check_script,
    overwrite = TRUE
  )

  script_copy_status <- if (
    isTRUE(copy_methods) &&
    isTRUE(copy_check)
  ) {
    "COPIED"
  } else {
    "COPY_FAILED"
  }
}

END_TIME <- Sys.time()
overall_status <- if (
  all(scientific_checks$status == "PASS")
) {
  "COMPLETED_STAGE3_READY_FOR_REVIEW"
} else {
  "COMPLETED_STAGE3_REVIEW_REQUIRED"
}

run_status <- data.table(
  stage = STAGE_NAME,
  start_time = format(
    START_TIME,
    "%Y-%m-%d %H:%M:%S"
  ),
  end_time = format(
    END_TIME,
    "%Y-%m-%d %H:%M:%S"
  ),
  elapsed_minutes = round(
    as.numeric(
      difftime(
        END_TIME,
        START_TIME,
        units = "mins"
      )
    ),
    2
  ),
  biological_samples = nrow(sample_meta),
  retained_cells = ncol(cardiac),
  retained_genes = nrow(cardiac),
  major_clusters =
    uniqueN(cardiac$seurat_clusters),
  major_cell_types =
    uniqueN(cardiac$major_cell_type),
  eligible_pseudobulk_cell_types =
    length(eligible_celltypes),
  pseudobulk_tested_genes =
    nrow(pseudobulk_results),
  macrophage_cells =
    length(macrophage_cells),
  macrophage_status =
    macrophage_status,
  doublet_status =
    doublet_status,
  warnings = nrow(warnings_dt),
  script_copy_status =
    script_copy_status,
  scientific_checks_failed =
    sum(scientific_checks$status != "PASS"),
  overall_status =
    overall_status
)

fwrite(
  run_status,
  file.path(
    DIRS$tables,
    "39_stage3_run_status.csv"
  )
)

readme <- c(
  "HFpEF Reanalysis Project - Stage 3",
  "GSE236585 cardiac scRNA-seq discovery and Stage 2 program projection",
  "",
  paste0(
    "Overall status: ",
    overall_status
  ),
  paste0(
    "Biological samples: ",
    nrow(sample_meta)
  ),
  paste0(
    "Retained cells: ",
    ncol(cardiac)
  ),
  paste0(
    "Major cell types: ",
    uniqueN(cardiac$major_cell_type)
  ),
  paste0(
    "Eligible pseudobulk cell types: ",
    length(eligible_celltypes)
  ),
  paste0(
    "Macrophage cells: ",
    length(macrophage_cells)
  ),
  paste0(
    "Doublet analysis: ",
    doublet_status
  ),
  paste0(
    "Script snapshot: ",
    script_copy_status
  ),
  "",
  "Primary interpretation boundary:",
  "- Cells are descriptive observations; biological samples are the inferential units.",
  "- Stage 2 signatures are projected without forcing Nfkb1 or a fixed original candidate panel.",
  "- GSE236585 contains no dapagliflozin exposure and therefore tests disease-program localization and concordance, not drug response.",
  "",
  "Upload the CHECK package for review before Stage 4."
)
writeLines(
  readme,
  file.path(
    OUT_DIR,
    "README_stage3.txt"
  ),
  useBytes = TRUE
)

############################################################
## 20. Compact CHECK package
############################################################

## Compact pseudobulk program tables.
write_csv_safe(
  pseudobulk_program_stats[
    signature_name %in%
      primary_signature_names
  ],
  file.path(
    DIRS$check,
    "PRIMARY_program_localization_statistics.csv"
  )
)
write_csv_safe(
  top_pseudobulk_de,
  file.path(
    DIRS$check,
    "TOP_pseudobulk_DE_by_celltype.csv"
  )
)
write_csv_safe(
  stage2_stage3_concordance,
  file.path(
    DIRS$check,
    "Stage2_Stage3_concordance.csv"
  )
)
write_csv_safe(
  macrophage_cluster_annotation,
  file.path(
    DIRS$check,
    "Macrophage_cluster_annotation.csv"
  )
)
write_csv_safe(
  macrophage_state_score_stats[
    signature_name %in%
      primary_signature_names
  ],
  file.path(
    DIRS$check,
    "Macrophage_primary_program_statistics.csv"
  )
)

review_files <- c(
  file.path(
    DIRS$tables,
    "01_locked_GSE236585_sample_metadata.csv"
  ),
  file.path(
    DIRS$tables,
    "03_stage2_signature_size_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "04_GSE236585_10x_file_mapping_and_dimensions.csv"
  ),
  file.path(
    DIRS$tables,
    "05_sample_specific_QC_thresholds.csv"
  ),
  file.path(
    DIRS$tables,
    "06_sample_QC_retention_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "07_scDblFinder_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "11_major_celltype_cluster_annotation.csv"
  ),
  file.path(
    DIRS$tables,
    "14_top_cluster_markers_for_review.csv"
  ),
  file.path(
    DIRS$tables,
    "15_celltype_composition_by_sample.csv"
  ),
  file.path(
    DIRS$tables,
    "16_celltype_composition_exploratory_statistics.csv"
  ),
  file.path(
    DIRS$tables,
    "17_signature_gene_coverage_in_GSE236585.csv"
  ),
  file.path(
    DIRS$tables,
    "21_major_celltype_pseudobulk_eligibility.csv"
  ),
  file.path(
    DIRS$tables,
    "23_major_celltype_pseudobulk_DE_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "27_pseudobulk_program_statistics.csv"
  ),
  file.path(
    DIRS$tables,
    "28_stage2_stage3_gene_level_concordance.csv"
  ),
  file.path(
    DIRS$tables,
    "30_macrophage_cluster_annotation.csv"
  ),
  file.path(
    DIRS$tables,
    "33_macrophage_state_composition_by_sample.csv"
  ),
  file.path(
    DIRS$tables,
    "35_macrophage_state_score_statistics.csv"
  ),
  workbook_path,
  file.path(
    DIRS$tables,
    "37_warnings_and_nonfatal_issues.csv"
  ),
  file.path(
    DIRS$tables,
    "38_scientific_completion_checks.csv"
  ),
  file.path(
    DIRS$tables,
    "39_stage3_run_status.csv"
  ),
  file.path(
    DIRS$methods,
    "stage3_parameters_and_rationale.csv"
  ),
  file.path(
    DIRS$methods,
    "stage3_methods_and_claim_boundaries.txt"
  ),
  file.path(
    DIRS$methods,
    "sessionInfo.txt"
  ),
  file.path(
    OUT_DIR,
    "README_stage3.txt"
  ),
  LOG_FILE,
  list.files(
    DIRS$figures,
    pattern = "\\.png$",
    full.names = TRUE
  )
)
review_files <- unique(
  review_files[
    file.exists(review_files)
  ]
)

for (f in review_files) {
  target <- file.path(
    DIRS$check,
    basename(f)
  )

  if (
    normalizePath(
      f,
      winslash = "/",
      mustWork = FALSE
    ) !=
      normalizePath(
        target,
        winslash = "/",
        mustWork = FALSE
      )
  ) {
    file.copy(
      f,
      target,
      overwrite = TRUE
    )
  }
}

check_files <- list.files(
  DIRS$check,
  full.names = TRUE
)
check_manifest <- data.table(
  filename = basename(check_files),
  size_bytes = as.numeric(
    file.info(check_files)$size
  )
)
check_manifest[
  ,
  sha256 := vapply(
    check_files,
    function(f) {
      digest::digest(
        file = f,
        algo = "sha256",
        serialize = FALSE
      )
    },
    character(1)
  )
]
fwrite(
  check_manifest,
  file.path(
    DIRS$check,
    "CHECK_package_file_manifest.csv"
  )
)

if (file.exists(CHECK_ZIP)) {
  unlink(
    CHECK_ZIP,
    force = TRUE
  )
}
zip::zipr(
  zipfile = CHECK_ZIP,
  files = list.files(
    DIRS$check,
    full.names = TRUE
  ),
  root = DIRS$check
)

log_msg("Stage 3 analysis finished.")
log_msg("Overall status: ", overall_status)
log_msg(
  "Retained cells: ",
  ncol(cardiac)
)
log_msg(
  "Major cell types: ",
  uniqueN(cardiac$major_cell_type)
)
log_msg(
  "Eligible pseudobulk cell types: ",
  length(eligible_celltypes)
)
log_msg(
  "Macrophage cells: ",
  length(macrophage_cells)
)
log_msg(
  "CHECK package: ",
  CHECK_ZIP
)

cat("\n============================================================\n")
cat("HFpEF Stage 3 GSE236585 analysis completed\n")
cat("Status: ", overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat("Upload the CHECK package for review before Stage 4.\n")
cat("============================================================\n")
