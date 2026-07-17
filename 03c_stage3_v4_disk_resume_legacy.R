############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 3 FIXED v4
## DISK-RESUME CONTINUATION AFTER THE FIRST SUPPLEMENT
##
## This is the correct continuation for the situation where R was
## CLOSED after running:
##   1) 3、第三阶段代码.R
##   2) 3.1、第三阶段代码、第一次补充、后面出错.R
##
## The first supplement completed and saved:
##   - the annotated/projected cardiac Seurat object;
##   - the macrophage/monocyte subclustered Seurat object;
##   - Tables 01-35, including pseudobulk and concordance results.
##
## It then stopped at the beginning of Section 17 because this line:
##   )condition_palette <- c(
## was syntactically invalid.
##
## This script starts in a NEW R session and reloads all required
## objects and tables from disk. It does NOT depend on objects remaining
## in memory.
##
## It DOES NOT repeat:
##   - raw 10x extraction;
##   - original QC;
##   - scDblFinder;
##   - normalization, PCA, neighbors, clustering, or UMAP;
##   - major-cluster FindAllMarkers;
##   - major-cell-type annotation;
##   - cell-level program scoring;
##   - pseudobulk edgeR or limma-voom;
##   - Stage 2-Stage 3 concordance;
##   - macrophage/monocyte subclustering or macrophage markers.
##
## It ONLY performs:
##   - Section 17 figures and figure source data;
##   - Section 18 workbook, methods, and parameters;
##   - Section 19 completion checks and run status;
##   - Section 20 compact CHECK package.
##
## Recommended run:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage3_GSE236585_v4_DISK_RESUME_AFTER_FIRST_SUPPLEMENT.R",
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
## 0. Paths and downstream-only replacement boundary
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
STAGE_NAME <- "03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH"
OUT_DIR <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP <- file.path(
  PROJECT_DIR,
  paste0(STAGE_NAME, "_CHECK.zip")
)

EXPECTED_SCRIPT_FILE <- file.path(PROJECT_DIR, "03c_stage3_v4_disk_resume_legacy.R")

SOURCE_STAGE3_V2_DIR <- file.path(
  PROJECT_DIR,
  "03_stage3_GSE236585_scRNA_projection_FIXED_v2"
)
SOURCE_CHECKPOINT <- file.path(
  SOURCE_STAGE3_V2_DIR,
  "02_objects",
  "CHECKPOINT_GSE236585_post_QC_pre_doublet.rds"
)

DIRS <- list(
  logs = file.path(OUT_DIR, "00_logs"),
  tables = file.path(OUT_DIR, "01_tables"),
  objects = file.path(OUT_DIR, "02_objects"),
  figures = file.path(OUT_DIR, "03_figures"),
  source = file.path(OUT_DIR, "04_source_data"),
  methods = file.path(OUT_DIR, "05_methods"),
  check = file.path(OUT_DIR, "06_review_check")
)

for (d in c(OUT_DIR, unlist(DIRS, use.names = FALSE))) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

CARDIAC_RDS <- file.path(
  DIRS$objects,
  "GSE236585_stage3_annotated_projected_seurat.rds"
)
MACROPHAGE_RDS <- file.path(
  DIRS$objects,
  "GSE236585_macrophage_subclustered_seurat.rds"
)

## Only downstream products are rebuilt. Tables 01-35 and the two RDS
## objects are never deleted or overwritten by this cleanup.
REBUILD_DOWNSTREAM_ONLY <- TRUE

downstream_cleanup_paths <- c(
  DIRS$figures,
  DIRS$source,
  DIRS$methods,
  DIRS$check,
  file.path(
    DIRS$tables,
    "36_GSE236585_stage3_key_results.xlsx"
  ),
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
    DIRS$tables,
    "40_disk_resume_provenance.csv"
  ),
  file.path(
    OUT_DIR,
    "README_stage3.txt"
  ),
  CHECK_ZIP
)

cleanup_audit <- data.table::data.table(
  path = downstream_cleanup_paths,
  existed_before = vapply(
    downstream_cleanup_paths,
    function(x) dir.exists(x) || file.exists(x),
    logical(1)
  ),
  deleted = FALSE
)

if (REBUILD_DOWNSTREAM_ONLY) {
  for (i in seq_len(nrow(cleanup_audit))) {
    target <- cleanup_audit$path[i]

    if (dir.exists(target)) {
      existing_children <- list.files(
        target,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
      )
      if (length(existing_children) > 0L) {
        unlink(
          existing_children,
          recursive = TRUE,
          force = TRUE
        )
      }
      cleanup_audit$deleted[i] <- TRUE
    } else if (file.exists(target)) {
      unlink(
        target,
        recursive = FALSE,
        force = TRUE
      )
      cleanup_audit$deleted[i] <- !file.exists(target)
    } else {
      cleanup_audit$deleted[i] <- TRUE
    }
  }
}

for (d in c(
  DIRS$figures,
  DIRS$source,
  DIRS$methods,
  DIRS$check
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

START_TIME <- Sys.time()
LOG_FILE <- file.path(
  DIRS$logs,
  "stage3_disk_resume_downstream.log"
)
UPSTREAM_LOG_FILE <- file.path(
  DIRS$logs,
  "stage3_GSE236585.log"
)
UPSTREAM_WARNING_FILE <- file.path(
  DIRS$logs,
  "stage3_warnings.log"
)

if (file.exists(LOG_FILE)) {
  unlink(LOG_FILE, force = TRUE)
}

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

ORIGINAL_ANALYSIS_START_TIME <- NA_character_
if (file.exists(UPSTREAM_LOG_FILE)) {
  upstream_log_lines <- readLines(
    UPSTREAM_LOG_FILE,
    warn = FALSE,
    encoding = "UTF-8"
  )
  timestamp_match <- regmatches(
    upstream_log_lines,
    regexpr(
      "\\[20[0-9]{2}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}\\]",
      upstream_log_lines
    )
  )
  timestamp_match <- timestamp_match[
    nzchar(timestamp_match)
  ]
  if (length(timestamp_match) > 0L) {
    ORIGINAL_ANALYSIS_START_TIME <- gsub(
      "^\\[|\\]$",
      "",
      timestamp_match[1L]
    )
  }
}

if (is.na(ORIGINAL_ANALYSIS_START_TIME)) {
  ORIGINAL_ANALYSIS_START_TIME <- "NOT_RECOVERED_FROM_LOG"
}

log_msg(
  "Stage 3 disk-resume continuation started after R restart."
)
log_msg(
  "Upstream Seurat objects and Tables 01-35 will be reloaded from disk."
)
log_msg(
  "No upstream analysis will be repeated."
)

############################################################
## 1. Package validation and utility functions
############################################################

required_packages <- c(
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
  "digest"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Required package(s) are missing: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running the disk-resume script."
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
})

write_csv_safe <- function(
  x,
  path,
  compress = FALSE
) {
  if (
    is.null(x) ||
    ncol(x) == 0L
  ) {
    data.table::fwrite(
      data.table::data.table(
        note = "No records generated."
      ),
      path
    )
  } else {
    data.table::fwrite(
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
  openxlsx::addWorksheet(wb, sheet)

  if (
    is.null(x) ||
    nrow(x) == 0L ||
    ncol(x) == 0L
  ) {
    openxlsx::writeData(
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

  openxlsx::writeData(wb, sheet, y)
  openxlsx::freezePane(
    wb,
    sheet,
    firstRow = TRUE
  )
  openxlsx::setColWidths(
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

  ggplot2::ggsave(
    png_path,
    plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    pdf_path,
    plot_object,
    width = width,
    height = height,
    bg = "white"
  )
  ggplot2::ggsave(
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

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes")] <- TRUE
  out[y %in% c("false", "f", "0", "no")] <- FALSE
  out
}

read_table_required <- function(filename) {
  path <- file.path(DIRS$tables, filename)
  if (!file.exists(path)) {
    stop("Required upstream table is missing: ", path)
  }
  data.table::fread(
    path,
    encoding = "UTF-8",
    na.strings = c("", "NA", "NaN")
  )
}

############################################################
## 2. Reload saved upstream objects and Tables 01-35
############################################################

required_disk_inputs <- c(
  CARDIAC_RDS,
  MACROPHAGE_RDS,
  SOURCE_CHECKPOINT,
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
    "04A_checkpoint_resume_audit.csv"
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
    "06A_checkpoint_cell_count_validation.csv"
  ),
  file.path(
    DIRS$tables,
    "07_scDblFinder_cell_calls.csv.gz"
  ),
  file.path(
    DIRS$tables,
    "07B_scDblFinder_rate_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "08_post_doublet_cell_counts.csv"
  ),
  file.path(
    DIRS$tables,
    "09_PCA_variance_and_dimensions.csv"
  ),
  file.path(
    DIRS$tables,
    "10_major_celltype_cluster_marker_scores.csv"
  ),
  file.path(
    DIRS$tables,
    "11_major_celltype_cluster_annotation.csv"
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
    "22_major_celltype_pseudobulk_DE_edgeR_limma.csv.gz"
  ),
  file.path(
    DIRS$tables,
    "23_major_celltype_pseudobulk_DE_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "24_top_pseudobulk_DE_genes_per_celltype.csv"
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
    "34_macrophage_state_sample_level_scores.csv"
  ),
  file.path(
    DIRS$tables,
    "35_macrophage_state_score_statistics.csv"
  )
)

missing_disk_inputs <- required_disk_inputs[
  !file.exists(required_disk_inputs)
]
if (length(missing_disk_inputs) > 0L) {
  stop(
    "The first supplement did not leave all required disk outputs. ",
    "Missing path(s):\n",
    paste(missing_disk_inputs, collapse = "\n")
  )
}

log_msg("Loading annotated cardiac Seurat object.")
cardiac <- readRDS(CARDIAC_RDS)
log_msg("Loading subclustered macrophage Seurat object.")
macrophage <- readRDS(MACROPHAGE_RDS)

if (!inherits(cardiac, "Seurat")) {
  stop("CARDIAC_RDS is not a Seurat object.")
}
if (!inherits(macrophage, "Seurat")) {
  stop("MACROPHAGE_RDS is not a Seurat object.")
}

sample_meta <- read_table_required(
  "01_locked_GSE236585_sample_metadata.csv"
)
signature_size_summary <- read_table_required(
  "03_stage2_signature_size_summary.csv"
)
checkpoint_audit <- read_table_required(
  "04A_checkpoint_resume_audit.csv"
)
file_map <- read_table_required(
  "04_GSE236585_10x_file_mapping_and_dimensions.csv"
)
qc_thresholds <- read_table_required(
  "05_sample_specific_QC_thresholds.csv"
)
qc_summary <- read_table_required(
  "06_sample_QC_retention_summary.csv"
)
qc_validation <- read_table_required(
  "06A_checkpoint_cell_count_validation.csv"
)
doublet_calls <- read_table_required(
  "07_scDblFinder_cell_calls.csv.gz"
)
doublet_rates <- read_table_required(
  "07B_scDblFinder_rate_summary.csv"
)
post_doublet_counts <- read_table_required(
  "08_post_doublet_cell_counts.csv"
)
pca_variance_table <- read_table_required(
  "09_PCA_variance_and_dimensions.csv"
)
cluster_score_long <- read_table_required(
  "10_major_celltype_cluster_marker_scores.csv"
)
cluster_annotation <- read_table_required(
  "11_major_celltype_cluster_annotation.csv"
)
composition_complete <- read_table_required(
  "15_celltype_composition_by_sample.csv"
)
composition_stats <- read_table_required(
  "16_celltype_composition_exploratory_statistics.csv"
)
signature_coverage <- read_table_required(
  "17_signature_gene_coverage_in_GSE236585.csv"
)
pseudobulk_eligibility <- read_table_required(
  "21_major_celltype_pseudobulk_eligibility.csv"
)
pseudobulk_results <- read_table_required(
  "22_major_celltype_pseudobulk_DE_edgeR_limma.csv.gz"
)
pseudobulk_summary <- read_table_required(
  "23_major_celltype_pseudobulk_DE_summary.csv"
)
top_pseudobulk_de <- read_table_required(
  "24_top_pseudobulk_DE_genes_per_celltype.csv"
)
pseudobulk_program_stats <- read_table_required(
  "27_pseudobulk_program_statistics.csv"
)
stage2_stage3_concordance <- read_table_required(
  "28_stage2_stage3_gene_level_concordance.csv"
)
macrophage_cluster_annotation <- read_table_required(
  "30_macrophage_cluster_annotation.csv"
)
macrophage_state_composition <- read_table_required(
  "33_macrophage_state_composition_by_sample.csv"
)
macrophage_state_scores <- read_table_required(
  "34_macrophage_state_sample_level_scores.csv"
)
macrophage_state_score_stats <- read_table_required(
  "35_macrophage_state_score_statistics.csv"
)

sample_meta[
  ,
  condition := factor(
    condition,
    levels = c("Control", "HFpEF")
  )
]
composition_complete[
  ,
  condition := factor(
    condition,
    levels = c("Control", "HFpEF")
  )
]
macrophage_state_scores[
  ,
  condition := factor(
    condition,
    levels = c("Control", "HFpEF")
  )
]

if (
  "checkpoint_matches_v2_QC" %in%
    names(qc_validation)
) {
  qc_validation[
    ,
    checkpoint_matches_v2_QC :=
      as_logical_safe(
        checkpoint_matches_v2_QC
      )
  ]
}
if ("eligible" %in% names(pseudobulk_eligibility)) {
  pseudobulk_eligibility[
    ,
    eligible :=
      as_logical_safe(eligible)
  ]
}

############################################################
## 3. Reconstruct compact variables and validate disk state
############################################################

MIN_FEATURES_HARD <- 200L
MAX_FEATURES_HARD <- 9000L
MAX_COUNTS_HARD <- 120000L
MAX_PERCENT_MT_HARD <- 25
QC_MAD_MULTIPLIER <- 4
QC_UPPER_QUANTILE <- 0.995

SCDOUBLETFINDER_DBR_PER_1K <- 0.008
doublet_status <- "COMPLETED_REQUIRED"

N_VARIABLE_FEATURES <- 3000L
MAJOR_CLUSTER_RESOLUTION <- 0.55
MACROPHAGE_CLUSTER_RESOLUTION <- 0.65

MIN_CELLS_PER_SAMPLE_CELLTYPE <- 20L
MIN_SAMPLES_PER_CONDITION <- 3L
EXCLUDED_FROM_PSEUDOBULK <- c(
  "Low_quality_mitochondrial",
  "Cycling_unresolved",
  "Unresolved"
)

SIGNATURE_TIERS <- c(
  "Tier_A_both_DESeq2_FDR_and_edgeR_direction",
  "Tier_B_one_DESeq2_FDR_effect_supported",
  "Tier_C_effect_and_method_supported"
)
SIGNATURE_SIZES <- c(50L, 100L, 150L, 200L)
PRIMARY_SIGNATURE_SIZE <- 150L

dims_use <- pca_variance_table[
  used_for_graph == TRUE,
  PC
]
if (length(dims_use) == 0L) {
  dims_use <- 1:30
}

primary_signature_names <- unique(
  signature_size_summary[
    requested_size_per_direction ==
      PRIMARY_SIGNATURE_SIZE,
    signature_name
  ]
)
expected_primary_signatures <- c(
  "Ccr2pos_Top150",
  "Ccr2neg_Top150",
  "CrossSubset_Top150"
)
if (
  length(primary_signature_names) != 3L ||
  !setequal(
    primary_signature_names,
    expected_primary_signatures
  )
) {
  stop(
    "Primary signature reconstruction failed. Observed: ",
    paste(primary_signature_names, collapse = ", ")
  )
}
primary_signature_names <-
  expected_primary_signatures

primary_score_columns <- paste0(
  "score_",
  primary_signature_names,
  "_net"
)

eligible_celltypes <- pseudobulk_eligibility[
  eligible == TRUE,
  major_cell_type
]
eligible_celltypes <- unique(
  as.character(eligible_celltypes)
)

macrophage_cells <- colnames(macrophage)
macrophage_status <- "COMPLETED"

expected_singlets <- sum(
  post_doublet_counts$retained_singlet_cells,
  na.rm = TRUE
)
expected_macrophage_cells <- sum(
  cardiac$major_cell_type ==
    "Macrophage_Monocyte",
  na.rm = TRUE
)

required_cardiac_metadata <- c(
  "sample_accession",
  "condition",
  "seurat_clusters",
  "major_cell_type",
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  primary_score_columns
)
missing_cardiac_metadata <- setdiff(
  required_cardiac_metadata,
  names(cardiac@meta.data)
)
if (length(missing_cardiac_metadata) > 0L) {
  stop(
    "The saved cardiac object is missing metadata column(s): ",
    paste(missing_cardiac_metadata, collapse = ", ")
  )
}

required_macrophage_metadata <- c(
  "sample_accession",
  "condition",
  "seurat_clusters",
  "macrophage_state"
)
missing_macrophage_metadata <- setdiff(
  required_macrophage_metadata,
  names(macrophage@meta.data)
)
if (length(missing_macrophage_metadata) > 0L) {
  stop(
    "The saved macrophage object is missing metadata column(s): ",
    paste(missing_macrophage_metadata, collapse = ", ")
  )
}

if (!"umap" %in% names(cardiac@reductions)) {
  stop("The saved cardiac object does not contain UMAP.")
}
if (!"umap" %in% names(macrophage@reductions)) {
  stop("The saved macrophage object does not contain UMAP.")
}
if (ncol(cardiac) != expected_singlets) {
  stop(
    "Cardiac object/table mismatch: object cells=",
    ncol(cardiac),
    "; expected singlets=",
    expected_singlets
  )
}
if (ncol(macrophage) != expected_macrophage_cells) {
  stop(
    "Macrophage object/annotation mismatch: object cells=",
    ncol(macrophage),
    "; cardiac macrophage labels=",
    expected_macrophage_cells
  )
}
if (
  data.table::uniqueN(
    sample_meta$sample_accession
  ) != 6L ||
  sum(sample_meta$condition == "Control") != 3L ||
  sum(sample_meta$condition == "HFpEF") != 3L
) {
  stop(
    "Sample metadata is not the expected 3 Control + 3 HFpEF design."
  )
}
if (
  data.table::uniqueN(
    doublet_calls$sample_accession
  ) != 6L
) {
  stop("Doublet-call table does not contain all six samples.")
}
if (
  any(
    !primary_score_columns %in%
      names(macrophage_state_scores)
  )
) {
  stop(
    "Macrophage score table is missing primary score columns."
  )
}

## Reconstruct a compact sample-level checkpoint metadata object solely
## for completion checks. No cell-level checkpoint is reloaded.
checkpoint_meta <- data.table::copy(sample_meta)

## Preserve upstream warnings, if any.
warning_records <- list()
if (
  file.exists(UPSTREAM_WARNING_FILE) &&
  file.info(UPSTREAM_WARNING_FILE)$size > 0
) {
  upstream_warning_lines <- readLines(
    UPSTREAM_WARNING_FILE,
    warn = FALSE,
    encoding = "UTF-8"
  )
  upstream_warning_lines <- upstream_warning_lines[
    nzchar(trimws(upstream_warning_lines))
  ]

  if (length(upstream_warning_lines) > 0L) {
    warning_records[[1L]] <- data.frame(
      timestamp = NA_character_,
      category = "UPSTREAM_WARNING_LOG",
      item = basename(UPSTREAM_WARNING_FILE),
      message = upstream_warning_lines,
      stringsAsFactors = FALSE
    )
  }
}

disk_resume_provenance <- data.table::data.table(
  continuation_type =
    "Independent disk resume after R restart",
  based_on_script_1 =
    "3、第三阶段代码.R",
  based_on_script_2 =
    "3.1、第三阶段代码、第一次补充、后面出错.R",
  resume_section =
    "Section 17: Figures and source data",
  cardiac_object =
    normalizePath(
      CARDIAC_RDS,
      winslash = "/",
      mustWork = TRUE
    ),
  macrophage_object =
    normalizePath(
      MACROPHAGE_RDS,
      winslash = "/",
      mustWork = TRUE
    ),
  cardiac_singlets = ncol(cardiac),
  macrophage_monocyte_cells =
    ncol(macrophage),
  upstream_analysis_repeated = FALSE,
  R_same_session_required = FALSE,
  original_analysis_start_time =
    ORIGINAL_ANALYSIS_START_TIME,
  downstream_resume_start_time =
    format(START_TIME, "%Y-%m-%d %H:%M:%S")
)

data.table::fwrite(
  disk_resume_provenance,
  file.path(
    DIRS$tables,
    "40_disk_resume_provenance.csv"
  )
)
data.table::fwrite(
  cleanup_audit,
  file.path(
    DIRS$logs,
    "stage3_disk_resume_cleanup_audit.csv"
  )
)

log_msg(
  "Disk state validated: ",
  ncol(cardiac),
  " cardiac singlets; ",
  ncol(macrophage),
  " macrophage/monocyte cells; ",
  length(eligible_celltypes),
  " eligible pseudobulk cell types."
)
log_msg(
  "Continuing directly at Section 17."
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
  "Epicardial_Mesothelial" = "#2CA02C",
  "Schwann_Glial" = "#8C6D31",
  "Platelet_Megakaryocyte" = "#C49C94",
  "Erythroid" = "#AD494A",
  "Cycling_unresolved" = "#637939",
  "Low_quality_mitochondrial" = "#BDBDBD",
  "Unresolved" = "#525252"
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
  "Disk_resume_provenance",
  disk_resume_provenance
)
write_sheet_safe(
  wb,
  "Sample_metadata",
  sample_meta
)
write_sheet_safe(
  wb,
  "Checkpoint_audit",
  checkpoint_audit
)
write_sheet_safe(
  wb,
  "Checkpoint_validation",
  qc_validation
)
write_sheet_safe(
  wb,
  "10x_mapping_inherited",
  file_map
)
write_sheet_safe(
  wb,
  "QC_thresholds",
  qc_thresholds
)
write_sheet_safe(
  wb,
  "QC_summary_inherited",
  qc_summary
)
write_sheet_safe(
  wb,
  "Doublet_rates",
  doublet_rates
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
    "Repair source checkpoint",
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
    paste0(
      doublet_status,
      "; per-sample scDblFinder; dbr.per1k=",
      SCDOUBLETFINDER_DBR_PER_1K
    ),
    normalizePath(
      SOURCE_CHECKPOINT,
      winslash = "/",
      mustWork = TRUE
    ),
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
    "Mandatory per-sample doublet classification; the repair stops if it fails",
    "Reuses completed v2 post-QC cells and does not repeat raw 10x extraction or original QC",
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
  "HFpEF Stage 3 FIXED v4: checkpoint-based doublet and annotation repair",
  "",
  "Input and biological design:",
  "- This downstream continuation was executed after R had been closed and restarted.",
  "- It loaded the annotated cardiac and macrophage Seurat objects saved by the first supplement.",
  "- It did not repeat raw 10x extraction, QC, scDblFinder, clustering, annotation, pseudobulk inference, or macrophage subclustering.",
  "- Six independent ventricular scRNA-seq samples: three HFpEF and three control mice.",
  "- Sample identities and conditions were taken from the locked Stage 1 manifest.",
  "- The repair loaded the completed v2 post-QC/pre-doublet Seurat checkpoint.",
  "- Raw 10x matrices were not re-read and original per-sample QC was not repeated.",
  "- Checkpoint sample-level cell counts were required to match the completed v2 QC summary.",
  "",
  "Quality control and doublets:",
  "- The completed v2 sample-specific QC results were inherited from the validated checkpoint.",
  "- scDblFinder was then run separately for each biological sample with a fixed random seed.",
  "- dbr.per1k was set to 0.008 (0.8% per 1,000 captured cells) when supported by the installed scDblFinder version.",
  "- Predicted doublets were removed before normalization, clustering, annotation, pseudobulk analysis, and program projection.",
  "- scDblFinder failure is fatal in this repair version.",
  "",
  "Clustering and annotation:",
  "- LogNormalize, 3,000 variable genes, percent-mitochondrial regression, PCA, graph clustering, and UMAP were repeated using singlets only.",
  "- Cluster annotation combined prespecified marker-program scores with deterministic top-marker rules.",
  "- The marker vocabulary explicitly included Schwann/glial, platelet/megakaryocyte, erythroid, cycling, and low-quality mitochondrial classes.",
  "- Low-quality mitochondrial, cycling-unresolved, and unresolved clusters were excluded from inferential pseudobulk analyses.",
  "- Annotation evidence, priority-marker hits, cluster QC, sample dominance, confidence, and manual-review flags were retained.",
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
    "Checkpoint biological samples",
    "Checkpoint matches completed v2 QC",
    "Mandatory scDblFinder completed",
    "Samples with valid doublet calls",
    "Samples retaining >=100 singlets",
    "Major cell types annotated",
    "No missing major-cell-type labels",
    "Excluded labels absent from eligible pseudobulk set",
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
    nrow(qc_validation),
    sum(
      qc_validation$
        checkpoint_matches_v2_QC
    ),
    as.integer(
      doublet_status ==
        "COMPLETED_REQUIRED"
    ),
    uniqueN(
      doublet_calls$sample_accession
    ),
    sum(
      post_doublet_counts$
        retained_singlet_cells >=
        100L
    ),
    uniqueN(cardiac$major_cell_type),
    sum(
      is.na(cardiac$major_cell_type) |
        !nzchar(cardiac$major_cell_type)
    ),
    sum(
      eligible_celltypes %in%
        EXCLUDED_FROM_PSEUDOBULK
    ),
    length(eligible_celltypes),
    nrow(pseudobulk_results),
    length(primary_signature_names),
    sum(
      primary_coverage_check$
        detected_up_genes >= 10L &
        primary_coverage_check$
          detected_down_genes >= 10L
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
    1L,
    6L,
    6L,
    5L,
    0L,
    0L,
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
    "equal",
    "equal",
    "equal",
    "at_least",
    "equal",
    "equal",
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
    "HFpEF_Stage3_GSE236585_v4_DISK_RESUME_AFTER_FIRST_SUPPLEMENT.R"
  )
  check_script <- file.path(
    DIRS$check,
    "HFpEF_Stage3_GSE236585_v4_DISK_RESUME_AFTER_FIRST_SUPPLEMENT.R"
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
  original_analysis_start_time =
    ORIGINAL_ANALYSIS_START_TIME,
  downstream_resume_start_time = format(
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
  "HFpEF Reanalysis Project - Stage 3 FIXED v4 PATCH",
  "GSE236585 checkpoint-based doublet and annotation repair",
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
  "Repair provenance:",
  "- R was closed after the first supplement stopped at the first line of Section 17.",
  "- This continuation independently reloaded the saved cardiac and macrophage Seurat objects from disk.",
  "- No completed upstream analysis was repeated.",
  "- Loaded the completed Stage 3 v2 post-QC/pre-doublet checkpoint.",
  "- Did not re-read raw 10x matrices or repeat original QC.",
  "- Recomputed all downstream results after mandatory doublet removal and revised annotation.",
  "",
  "Primary interpretation boundary:",
  "- Cells are descriptive observations; biological samples are the inferential units.",
  "- Stage 2 signatures are projected without forcing Nfkb1 or a fixed original candidate panel.",
  "- GSE236585 contains no dapagliflozin exposure and therefore tests disease-program localization and concordance, not drug response.",
  "",
  "Upload the v4 PATCH CHECK package for review before Stage 4."
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
    "40_disk_resume_provenance.csv"
  ),
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
    "04A_checkpoint_resume_audit.csv"
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
    "06A_checkpoint_cell_count_validation.csv"
  ),
  file.path(
    DIRS$tables,
    "07A_scDblFinder_class_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "07B_scDblFinder_rate_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "11_major_celltype_cluster_annotation.csv"
  ),
  file.path(
    DIRS$tables,
    "14A_cluster_priority_marker_evidence.csv"
  ),
  file.path(
    DIRS$tables,
    "14B_cluster_QC_and_sample_balance.csv"
  ),
  file.path(
    DIRS$tables,
    "14C_known_misannotation_class_audit.csv"
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

log_msg("Stage 3 disk-resume downstream analysis finished.")
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
cat("HFpEF Stage 3 FIXED v4 disk-resume continuation completed\n")
cat("Status: ", overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat("Upload the v4 PATCH CHECK package for review before Stage 4.\n")
cat("============================================================\n")
