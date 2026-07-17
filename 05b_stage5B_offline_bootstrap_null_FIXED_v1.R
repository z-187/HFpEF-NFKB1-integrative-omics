############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 4 FIXED v1
## GSE236585 macrophage TF-regulon consensus prioritization
##
## This stage is rebuilt for the CURRENT project:
##   <HFPEF_PROJECT_DIR>
##
## Required completed inputs:
##   Stage 2:
##     02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2
##   Stage 3:
##     03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH
##
## Scientific objectives:
##   1) Load the locked Stage 3 macrophage/monocyte object without
##      repeating QC, doublet detection, clustering, or annotation.
##   2) obtain a broad mouse TF-target prior without forcing Nfkb1;
##   3) estimate signed regulon activity using an always-available
##      weighted target score and, when available, AUCell;
##   4) test HFpEF versus Control at the BIOLOGICAL-SAMPLE level;
##   5) compare TF activity with TF-expression-only ranking;
##   6) quantify leave-one-Control/leave-one-HFpEF-pair robustness;
##   7) test regulon overlap with Stage 2 drug-opposed programs and
##      their Stage 3 macrophage-supported subsets;
##   8) produce an auditable TF priority table without preselecting
##      Nfkb1 or any fixed TF candidate panel.
##
## Interpretation boundary:
##   - GSE236585 has no dapagliflozin treatment.
##   - Stage 4 performs regulatory PRIORITIZATION, not causal proof.
##   - Biological samples, not cells, are the inferential units.
##   - Cell-level activity maps are descriptive only.
##   - A prior-based regulon network is preferred. A clearly labelled
##     data-driven correlation fallback is used only if no prior network
##     can be loaded.
##
## Output:
##   <HFPEF_PROJECT_DIR>/
##   04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1
##
## CHECK package:
##   <HFPEF_PROJECT_DIR>/
##   04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1_CHECK.zip
##
## Recommended run:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage4_GSE236585_Macrophage_TF_Regulon_FIXED_v1.R",
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
## 0. User settings and exact input paths
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

STAGE3_DIR <- file.path(
  PROJECT_DIR,
  "03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH"
)
STAGE3_STATUS_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "39_stage3_run_status.csv"
)
STAGE3_CHECKS_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "38_scientific_completion_checks.csv"
)
STAGE3_MACROPHAGE_RDS <- file.path(
  STAGE3_DIR,
  "02_objects",
  "GSE236585_macrophage_subclustered_seurat.rds"
)
STAGE3_CARDIAC_RDS <- file.path(
  STAGE3_DIR,
  "02_objects",
  "GSE236585_stage3_annotated_projected_seurat.rds"
)
STAGE3_SAMPLE_META_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "01_locked_GSE236585_sample_metadata.csv"
)
STAGE3_MAC_STATE_ANNOTATION_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "30_macrophage_cluster_annotation.csv"
)
STAGE3_MAC_STATE_COMPOSITION_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "33_macrophage_state_composition_by_sample.csv"
)
STAGE3_MAC_STATE_SCORE_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "35_macrophage_state_score_statistics.csv"
)
STAGE3_MAJOR_PB_DE_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "22_major_celltype_pseudobulk_DE_edgeR_limma.csv.gz"
)
STAGE3_CONCORDANCE_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "28_stage2_stage3_gene_level_concordance.csv"
)

STAGE_NAME <-
  "04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1"
OUT_DIR <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP <- file.path(
  PROJECT_DIR,
  paste0(STAGE_NAME, "_CHECK.zip")
)
EXPECTED_SCRIPT_FILE <- file.path(
  PROJECT_DIR,
  "R",
  "04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1.R"
)

REPLACE_EXISTING_STAGE4 <- TRUE
INSTALL_OPTIONAL_PACKAGES <- TRUE
RUN_AUCELL_SENSITIVITY <- TRUE

## Program construction.
SIGNATURE_TIERS <- c(
  "Tier_A_both_DESeq2_FDR_and_edgeR_direction",
  "Tier_B_one_DESeq2_FDR_effect_supported",
  "Tier_C_effect_and_method_supported"
)
SIGNATURE_SIZES <- c(50L, 100L, 150L, 200L)
PRIMARY_SIGNATURE_SIZE <- 150L

## Regulon filtering.
DOROTHEA_CONFIDENCE_LEVELS <- c("A", "B", "C")
MIN_TARGETS_PER_REGULON <- 10L
MAX_TARGETS_PER_REGULON <- 250L
MIN_TARGET_DETECTION_FRACTION <- 0.01
MIN_TF_DETECTION_FRACTION <- 0.01
MIN_STATE_SAMPLE_CELLS <- 20L

## Data-driven fallback parameters.
FALLBACK_MIN_ABS_COR <- 0.25
FALLBACK_TOP_TARGETS_PER_TF <- 120L
FALLBACK_MIN_PROFILES <- 6L

## Reporting.
TOP_TFS_REPORT <- 25L
TOP_TFS_HEATMAP <- 20L
FORMAL_FDR <- 0.05
EXPLORATORY_FDR <- 0.10

############################################################
## 1. Preflight, output replacement, and logging
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
  STAGE2_STATUS_FILE,
  STAGE2_POS_FILE,
  STAGE2_NEG_FILE,
  STAGE2_CROSS_FILE,
  STAGE3_STATUS_FILE,
  STAGE3_CHECKS_FILE,
  STAGE3_MACROPHAGE_RDS,
  STAGE3_CARDIAC_RDS,
  STAGE3_SAMPLE_META_FILE,
  STAGE3_MAC_STATE_ANNOTATION_FILE,
  STAGE3_MAC_STATE_COMPOSITION_FILE,
  STAGE3_MAC_STATE_SCORE_FILE,
  STAGE3_MAJOR_PB_DE_FILE,
  STAGE3_CONCORDANCE_FILE
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Required Stage 2/3 input path(s) are missing:\n",
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
    "Stage 2 is not locked as COMPLETED_STAGE2_READY_FOR_REVIEW."
  )
}

stage3_status <- data.table::fread(
  STAGE3_STATUS_FILE,
  encoding = "UTF-8"
)
if (
  !"overall_status" %in% names(stage3_status) ||
  stage3_status$overall_status[1L] !=
    "COMPLETED_STAGE3_READY_FOR_REVIEW"
) {
  stop(
    "Stage 3 is not locked as COMPLETED_STAGE3_READY_FOR_REVIEW."
  )
}

stage3_checks <- data.table::fread(
  STAGE3_CHECKS_FILE,
  encoding = "UTF-8"
)
if (
  !all(c("check", "status") %in% names(stage3_checks)) ||
  any(stage3_checks$status != "PASS")
) {
  stop(
    "At least one Stage 3 scientific completion check is not PASS."
  )
}

replacement_audit <- data.table::data.table(
  path = c(OUT_DIR, CHECK_ZIP),
  path_type = c(
    "stage4_output_directory",
    "stage4_check_zip"
  ),
  existed_before = FALSE,
  deletion_attempted = FALSE,
  deletion_succeeded = FALSE
)

if (REPLACE_EXISTING_STAGE4) {
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
        stop("Failed to remove previous Stage 4 path: ", target)
      }
    } else {
      replacement_audit$deletion_succeeded[i] <- TRUE
    }
  }
} else if (dir.exists(OUT_DIR) || file.exists(CHECK_ZIP)) {
  stop(
    "Existing Stage 4 output detected while replacement is disabled."
  )
}

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
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

START_TIME <- Sys.time()
LOG_FILE <- file.path(DIRS$logs, "stage4_GSE236585_TF_regulon.log")
WARN_FILE <- file.path(DIRS$logs, "stage4_warnings.log")

data.table::fwrite(
  replacement_audit,
  file.path(DIRS$logs, "stage4_replacement_audit.csv")
)

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

log_msg("Stage 4 TF-regulon consensus analysis started.")
log_msg("PROJECT_DIR: ", PROJECT_DIR)
log_msg("STAGE2_DIR: ", STAGE2_DIR)
log_msg("STAGE3_DIR: ", STAGE3_DIR)
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

  if (
    length(missing) > 0L &&
    INSTALL_OPTIONAL_PACKAGES
  ) {
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
    "digest"
  ),
  required = TRUE
)

ensure_bioc(
  c("edgeR", "limma"),
  required = TRUE
)

## These are optional because a deterministic fallback is implemented.
ensure_bioc(
  c(
    "dorothea",
    "decoupleR",
    "AUCell",
    "AnnotationDbi",
    "org.Mm.eg.db"
  ),
  required = FALSE
)

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

package_status <- data.table::data.table(
  package = c(
    "dorothea",
    "decoupleR",
    "AUCell",
    "AnnotationDbi",
    "org.Mm.eg.db"
  ),
  available = vapply(
    c(
      "dorothea",
      "decoupleR",
      "AUCell",
      "AnnotationDbi",
      "org.Mm.eg.db"
    ),
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)

data.table::fwrite(
  package_status,
  file.path(DIRS$tables, "00_stage4_optional_package_status.csv")
)

############################################################
## 3. Utility functions
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
  x <- sub("([._-][0-9]+)$", "", x)
  toupper(x)
}

rescale01 <- function(x, neutral_if_constant = 0) {
  x <- as.numeric(x)
  out <- rep(neutral_if_constant, length(x))
  finite <- is.finite(x)
  if (!any(finite)) return(out)
  rng <- range(x[finite], na.rm = TRUE)
  if (!is.finite(rng[1L]) || !is.finite(rng[2L])) {
    return(out)
  }
  if (diff(rng) == 0) {
    out[finite] <- neutral_if_constant
    return(out)
  }
  out[finite] <- (x[finite] - rng[1L]) / diff(rng)
  out
}

safe_neglog10 <- function(p, cap = 10) {
  p <- as.numeric(p)
  p[!is.finite(p) | is.na(p)] <- 1
  p[p <= 0] <- .Machine$double.xmin
  pmin(-log10(p), cap)
}

hedges_g <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  n1 <- length(x)
  n2 <- length(y)
  if (n1 < 2L || n2 < 2L) return(NA_real_)

  pooled_var <- (
    (n1 - 1L) * stats::var(x) +
      (n2 - 1L) * stats::var(y)
  ) / (n1 + n2 - 2L)

  if (!is.finite(pooled_var) || pooled_var <= 0) {
    return(NA_real_)
  }

  d <- (mean(x) - mean(y)) / sqrt(pooled_var)
  correction <- 1 - 3 / (4 * (n1 + n2) - 9)
  correction * d
}

safe_wilcox_p <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 1L || length(y) < 1L) {
    return(NA_real_)
  }
  tryCatch(
    stats::wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  suppressWarnings(
    stats::cor(
      x[keep],
      y[keep],
      method = "spearman"
    )
  )
}

write_csv_safe <- function(x, path, compress = FALSE) {
  if (is.null(x) || ncol(x) == 0L) {
    data.table::fwrite(
      data.table::data.table(note = "No records generated."),
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
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
  substr(x, 1L, 31L)
}

write_sheet_safe <- function(wb, sheet, x) {
  sheet <- sanitize_sheet_name(sheet)
  openxlsx::addWorksheet(wb, sheet)

  if (is.null(x) || nrow(x) == 0L || ncol(x) == 0L) {
    openxlsx::writeData(
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

  openxlsx::writeData(wb, sheet, y)
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  openxlsx::setColWidths(
    wb,
    sheet,
    cols = seq_len(min(ncol(y), 35L)),
    widths = "auto"
  )
  invisible(NULL)
}

save_plot_bundle <- function(plot_object, stem, width, height) {
  png_path <- file.path(DIRS$figures, paste0(stem, ".png"))
  pdf_path <- file.path(DIRS$figures, paste0(stem, ".pdf"))
  tiff_path <- file.path(DIRS$figures, paste0(stem, ".tiff"))

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

  invisible(c(png = png_path, pdf = pdf_path, tiff = tiff_path))
}

join_layers_safe <- function(object, assay = "RNA") {
  DefaultAssay(object) <- assay

  out <- tryCatch(
    {
      if (
        "JoinLayers" %in%
          getNamespaceExports("SeuratObject")
      ) {
        SeuratObject::JoinLayers(object, assay = assay)
      } else if (exists("JoinLayers", mode = "function")) {
        JoinLayers(object, assay = assay)
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

  tryCatch(
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
}

make_feature_map <- function(features) {
  dt <- data.table::data.table(
    feature = as.character(features),
    feature_key = gene_key(features)
  )
  data.table::setorder(dt, feature_key, feature)
  dt <- dt[, .SD[1L], by = feature_key]
  dt
}

map_symbols_to_features <- function(symbols, feature_map) {
  query <- data.table::data.table(
    symbol = as.character(symbols),
    feature_key = gene_key(symbols)
  )
  out <- merge(
    query,
    feature_map,
    by = "feature_key",
    all.x = TRUE,
    sort = FALSE
  )
  out
}

scale_rows <- function(m) {
  m <- as.matrix(m)
  row_means <- rowMeans(m, na.rm = TRUE)
  row_sds <- apply(m, 1L, stats::sd, na.rm = TRUE)
  row_sds[!is.finite(row_sds) | row_sds == 0] <- 1
  z <- sweep(m, 1L, row_means, "-")
  z <- sweep(z, 1L, row_sds, "/")
  z[!is.finite(z)] <- 0
  z
}

aggregate_sparse_counts <- function(counts, groups) {
  groups <- as.character(groups)
  if (length(groups) != ncol(counts)) {
    stop("Grouping vector length does not match count-matrix columns.")
  }

  group_levels <- sort(unique(groups))
  group_factor <- factor(groups, levels = group_levels)
  design <- Matrix::sparse.model.matrix(
    ~ 0 + group_factor
  )
  colnames(design) <- group_levels
  aggregated <- counts %*% design
  colnames(aggregated) <- group_levels
  aggregated
}

logcpm_from_counts <- function(counts) {
  y <- edgeR::DGEList(counts = as.matrix(counts))
  y <- edgeR::calcNormFactors(y)
  edgeR::cpm(y, log = TRUE, prior.count = 2)
}

weighted_regulon_activity <- function(expr_matrix, network_dt) {
  required_cols <- c(
    "source_symbol",
    "target_feature",
    "mor",
    "weight"
  )
  missing_cols <- setdiff(required_cols, names(network_dt))
  if (length(missing_cols) > 0L) {
    stop(
      "Network is missing column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  target_features <- intersect(
    unique(network_dt$target_feature),
    rownames(expr_matrix)
  )
  if (length(target_features) < 10L) {
    stop("Too few network targets are present in the expression matrix.")
  }

  z <- scale_rows(expr_matrix[target_features, , drop = FALSE])
  tfs <- sort(unique(network_dt$source_symbol))
  out <- matrix(
    NA_real_,
    nrow = length(tfs),
    ncol = ncol(z),
    dimnames = list(tfs, colnames(z))
  )

  for (tf in tfs) {
    net_tf <- network_dt[
      source_symbol == tf &
        target_feature %in% rownames(z)
    ]
    if (nrow(net_tf) < MIN_TARGETS_PER_REGULON) next

    net_tf <- net_tf[
      match(unique(target_feature), target_feature)
    ]
    weights <- as.numeric(net_tf$mor) * as.numeric(net_tf$weight)
    weights[!is.finite(weights)] <- 0
    denom <- sum(abs(weights))
    if (!is.finite(denom) || denom <= 0) next

    target_z <- z[
      net_tf$target_feature,
      ,
      drop = FALSE
    ]
    out[tf, ] <- colSums(
      sweep(target_z, 1L, weights, "*")
    ) / denom
  }

  out <- out[rowSums(is.finite(out)) > 0L, , drop = FALSE]
  out
}

sample_level_activity_test <- function(
  activity_matrix,
  sample_metadata,
  method_label
) {
  sample_metadata <- data.table::copy(sample_metadata)
  sample_metadata <- sample_metadata[
    match(colnames(activity_matrix), sample_accession)
  ]

  if (
    any(is.na(sample_metadata$sample_accession)) ||
    any(sample_metadata$sample_accession != colnames(activity_matrix))
  ) {
    stop("Sample metadata could not be aligned to activity matrix.")
  }

  sample_metadata[, condition := factor(
    condition,
    levels = c("Control", "HFpEF")
  )]

  design <- model.matrix(
    ~ condition,
    data = sample_metadata
  )

  fit <- limma::lmFit(activity_matrix, design)
  fit <- limma::eBayes(fit, robust = TRUE)
  tab <- limma::topTable(
    fit,
    coef = "conditionHFpEF",
    number = Inf,
    sort.by = "none"
  )
  data.table::setDT(tab, keep.rownames = "tf_symbol")

  records <- lapply(
    rownames(activity_matrix),
    function(tf) {
      values <- as.numeric(activity_matrix[tf, ])
      control <- values[sample_metadata$condition == "Control"]
      hfpef <- values[sample_metadata$condition == "HFpEF"]

      data.table::data.table(
        tf_symbol = tf,
        method = method_label,
        control_mean = mean(control, na.rm = TRUE),
        hfpef_mean = mean(hfpef, na.rm = TRUE),
        hfpef_minus_control =
          mean(hfpef, na.rm = TRUE) -
          mean(control, na.rm = TRUE),
        hedges_g_HFpEF_vs_Control = hedges_g(hfpef, control),
        wilcoxon_p = safe_wilcox_p(hfpef, control),
        control_samples = sum(is.finite(control)),
        hfpef_samples = sum(is.finite(hfpef))
      )
    }
  )

  stats_dt <- data.table::rbindlist(records)
  stats_dt[, wilcoxon_fdr := p.adjust(wilcoxon_p, method = "BH")]

  limma_dt <- tab[
    ,
    .(
      tf_symbol,
      limma_logFC = logFC,
      limma_AveExpr = AveExpr,
      limma_t = t,
      limma_pvalue = P.Value,
      limma_padj = adj.P.Val
    )
  ]

  merge(
    stats_dt,
    limma_dt,
    by = "tf_symbol",
    all.x = TRUE
  )
}

save_heatmap_bundle <- function(
  mat,
  stem,
  annotation_col = NULL,
  width = 10,
  height = 8
) {
  png_path <- file.path(DIRS$figures, paste0(stem, ".png"))
  pdf_path <- file.path(DIRS$figures, paste0(stem, ".pdf"))
  tiff_path <- file.path(DIRS$figures, paste0(stem, ".tiff"))

  grDevices::png(
    png_path,
    width = width * 300,
    height = height * 300,
    res = 300
  )
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    border_color = NA,
    main = "Sample-level macrophage TF activity"
  )
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = width, height = height)
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    border_color = NA,
    main = "Sample-level macrophage TF activity"
  )
  grDevices::dev.off()

  grDevices::tiff(
    tiff_path,
    width = width,
    height = height,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    border_color = NA,
    main = "Sample-level macrophage TF activity"
  )
  grDevices::dev.off()

  invisible(c(png = png_path, pdf = pdf_path, tiff = tiff_path))
}

############################################################
## 4. Load and validate Stage 3 objects
############################################################

sample_meta <- data.table::fread(
  STAGE3_SAMPLE_META_FILE,
  encoding = "UTF-8"
)
sample_meta[, condition := factor(
  condition,
  levels = c("Control", "HFpEF")
)]
data.table::setorder(sample_meta, condition, sample_accession)

if (
  data.table::uniqueN(sample_meta$sample_accession) != 6L ||
  sum(sample_meta$condition == "Control") != 3L ||
  sum(sample_meta$condition == "HFpEF") != 3L
) {
  stop("Stage 3 sample metadata is not the expected 3 + 3 design.")
}

log_msg("Loading Stage 3 macrophage/monocyte Seurat object.")
macrophage <- readRDS(STAGE3_MACROPHAGE_RDS)
if (!inherits(macrophage, "Seurat")) {
  stop("Stage 3 macrophage RDS is not a Seurat object.")
}
macrophage <- join_layers_safe(macrophage, assay = "RNA")

log_msg("Loading Stage 3 annotated cardiac Seurat object.")
cardiac <- readRDS(STAGE3_CARDIAC_RDS)
if (!inherits(cardiac, "Seurat")) {
  stop("Stage 3 cardiac RDS is not a Seurat object.")
}
cardiac <- join_layers_safe(cardiac, assay = "RNA")

required_mac_meta <- c(
  "sample_accession",
  "condition",
  "seurat_clusters",
  "macrophage_state"
)
missing_mac_meta <- setdiff(
  required_mac_meta,
  names(macrophage@meta.data)
)
if (length(missing_mac_meta) > 0L) {
  stop(
    "Macrophage object is missing metadata column(s): ",
    paste(missing_mac_meta, collapse = ", ")
  )
}

if (ncol(macrophage) != 1822L) {
  stop(
    "Unexpected Stage 3 macrophage cell count. Expected 1822; observed ",
    ncol(macrophage),
    "."
  )
}

mac_meta <- data.table::as.data.table(
  macrophage@meta.data,
  keep.rownames = "cell"
)
if (
  data.table::uniqueN(mac_meta$sample_accession) != 6L ||
  any(!mac_meta$sample_accession %in% sample_meta$sample_accession)
) {
  stop("Macrophage object does not contain the six locked samples.")
}

mac_meta[, condition := factor(
  condition,
  levels = c("Control", "HFpEF")
)]

mac_state_annotation <- data.table::fread(
  STAGE3_MAC_STATE_ANNOTATION_FILE,
  encoding = "UTF-8"
)
mac_state_composition <- data.table::fread(
  STAGE3_MAC_STATE_COMPOSITION_FILE,
  encoding = "UTF-8"
)
mac_state_score_stats <- data.table::fread(
  STAGE3_MAC_STATE_SCORE_FILE,
  encoding = "UTF-8"
)
stage3_concordance <- data.table::fread(
  STAGE3_CONCORDANCE_FILE,
  encoding = "UTF-8"
)

input_audit <- data.table::data.table(
  item = c(
    "Stage2_status",
    "Stage3_status",
    "Stage3_failed_checks",
    "Biological_samples",
    "Control_samples",
    "HFpEF_samples",
    "Macrophage_cells",
    "Macrophage_states",
    "Cardiac_cells"
  ),
  value = c(
    stage2_status$overall_status[1L],
    stage3_status$overall_status[1L],
    sum(stage3_checks$status != "PASS"),
    data.table::uniqueN(sample_meta$sample_accession),
    sum(sample_meta$condition == "Control"),
    sum(sample_meta$condition == "HFpEF"),
    ncol(macrophage),
    data.table::uniqueN(macrophage$macrophage_state),
    ncol(cardiac)
  )
)

data.table::fwrite(
  input_audit,
  file.path(DIRS$tables, "01_stage4_input_validation.csv")
)

############################################################
## 5. Build Stage 2 programs and Stage 3-supported subsets
############################################################

stage2_pos <- data.table::fread(
  STAGE2_POS_FILE,
  encoding = "UTF-8"
)
stage2_neg <- data.table::fread(
  STAGE2_NEG_FILE,
  encoding = "UTF-8"
)
stage2_cross <- data.table::fread(
  STAGE2_CROSS_FILE,
  encoding = "UTF-8"
)

required_subset_columns <- c(
  "symbol",
  "disease_lfc",
  "drug_lfc",
  "deseq_opposed",
  "edger_opposed",
  "four_effect_signs_consistent",
  "opposition_tier",
  "within_subset_rank"
)

for (nm in required_subset_columns) {
  if (
    !nm %in% names(stage2_pos) ||
    !nm %in% names(stage2_neg)
  ) {
    stop("Stage 2 subset table is missing column: ", nm)
  }
}

prepare_subset_program <- function(x, subset_name) {
  y <- data.table::copy(x)
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
  y[, abs_disease_lfc_order := abs(disease_lfc)]
  y[, abs_drug_lfc_order := abs(drug_lfc)]
  data.table::setorder(
    y,
    within_subset_rank,
    -abs_disease_lfc_order,
    -abs_drug_lfc_order
  )
  y <- y[, .SD[1L], by = symbol_key]
  y[, direction := data.table::fcase(
    disease_lfc > 0 & drug_lfc < 0,
    "Disease_up_Drug_down",
    disease_lfc < 0 & drug_lfc > 0,
    "Disease_down_Drug_up",
    default = "Other"
  )]
  y <- y[direction != "Other"]
  y[, subset_name := subset_name]
  y
}

pos_base <- prepare_subset_program(stage2_pos, "Ccr2pos")
neg_base <- prepare_subset_program(stage2_neg, "Ccr2neg")

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
    "Stage 2 cross-subset table is missing column(s): ",
    paste(missing_cross_columns, collapse = ", ")
  )
}

cross_base <- data.table::copy(stage2_cross)
cross_base <- cross_base[
  !is.na(symbol) &
    nzchar(symbol) &
    consensus_category ==
      "Cross_subset_full_directional_consensus" &
    is.finite(pos_disease_lfc) &
    is.finite(pos_drug_lfc) &
    is.finite(neg_disease_lfc) &
    is.finite(neg_drug_lfc)
]
cross_base[, symbol_key := gene_key(symbol)]
data.table::setorder(cross_base, overall_consensus_rank)
cross_base <- cross_base[, .SD[1L], by = symbol_key]
cross_base[, disease_lfc := rowMeans(
  cbind(pos_disease_lfc, neg_disease_lfc),
  na.rm = TRUE
)]
cross_base[, drug_lfc := rowMeans(
  cbind(pos_drug_lfc, neg_drug_lfc),
  na.rm = TRUE
)]
cross_base[, within_subset_rank := overall_consensus_rank]
cross_base[, direction := data.table::fcase(
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
cross_base <- cross_base[direction != "Other"]
cross_base[, subset_name := "CrossSubset"]

build_program_manifest <- function(base_dt) {
  records <- list()
  for (n_target in SIGNATURE_SIZES) {
    for (dir_value in c(
      "Disease_up_Drug_down",
      "Disease_down_Drug_up"
    )) {
      selected <- base_dt[
        direction == dir_value
      ][
        order(within_subset_rank)
      ][
        seq_len(min(n_target, .N))
      ]
      if (nrow(selected) == 0L) next
      selected[, signature_size := n_target]
      selected[, program_name := paste(
        subset_name,
        dir_value,
        paste0("Top", n_target),
        sep = "_"
      )]
      records[[length(records) + 1L]] <- selected[
        ,
        .(
          program_name,
          subset_name,
          direction,
          signature_size,
          symbol,
          symbol_key,
          stage2_disease_lfc = disease_lfc,
          stage2_drug_lfc = drug_lfc,
          stage2_rank = within_subset_rank
        )
      ]
    }
  }
  data.table::rbindlist(records, use.names = TRUE, fill = TRUE)
}

program_manifest <- data.table::rbindlist(
  list(
    build_program_manifest(pos_base),
    build_program_manifest(neg_base),
    build_program_manifest(cross_base)
  ),
  use.names = TRUE,
  fill = TRUE
)

stage3_pb_de <- data.table::fread(
  STAGE3_MAJOR_PB_DE_FILE,
  encoding = "UTF-8"
)
required_pb_cols <- c(
  "feature",
  "feature_key",
  "major_cell_type",
  "edgeR_logFC",
  "limma_logFC",
  "edgeR_limma_sign_agreement"
)
missing_pb_cols <- setdiff(required_pb_cols, names(stage3_pb_de))
if (length(missing_pb_cols) > 0L) {
  stop(
    "Stage 3 pseudobulk DE table is missing column(s): ",
    paste(missing_pb_cols, collapse = ", ")
  )
}

mac_pb_de <- stage3_pb_de[
  major_cell_type == "Macrophage_Monocyte"
]
mac_pb_de <- mac_pb_de[
  match(unique(feature_key), feature_key)
]

primary_manifest <- program_manifest[
  signature_size == PRIMARY_SIGNATURE_SIZE
]
primary_manifest <- merge(
  primary_manifest,
  mac_pb_de[
    ,
    .(
      symbol_key = feature_key,
      stage3_feature = feature,
      stage3_edgeR_logFC = edgeR_logFC,
      stage3_limma_logFC = limma_logFC,
      stage3_edgeR_limma_sign_agreement =
        edgeR_limma_sign_agreement
    )
  ],
  by = "symbol_key",
  all.x = TRUE
)

primary_manifest[, stage3_direction_supported := (
  !is.na(stage3_edgeR_logFC) &
    !is.na(stage3_limma_logFC) &
    stage3_edgeR_limma_sign_agreement == TRUE &
    sign(stage3_edgeR_logFC) == sign(stage2_disease_lfc) &
    sign(stage3_limma_logFC) == sign(stage2_disease_lfc)
)]

supported_programs <- primary_manifest[
  stage3_direction_supported == TRUE
]
supported_programs[, program_name := paste0(
  program_name,
  "_Stage3Supported"
)]

program_manifest_all <- data.table::rbindlist(
  list(program_manifest, supported_programs),
  use.names = TRUE,
  fill = TRUE
)

data.table::fwrite(
  program_manifest_all,
  file.path(DIRS$tables, "02_stage4_program_gene_manifest.csv")
)

program_size_summary <- program_manifest_all[
  ,
  .(
    genes = data.table::uniqueN(symbol_key),
    stage3_supported_genes = sum(
      stage3_direction_supported == TRUE,
      na.rm = TRUE
    )
  ),
  by = .(
    program_name,
    subset_name,
    direction,
    signature_size
  )
]

data.table::fwrite(
  program_size_summary,
  file.path(DIRS$tables, "03_stage4_program_size_summary.csv")
)

############################################################
## 6. Construct macrophage pseudobulk expression matrices
############################################################

counts_mac <- get_assay_matrix(
  macrophage,
  layer = "counts",
  assay = "RNA"
)

mac_meta <- mac_meta[
  match(colnames(counts_mac), cell)
]
if (
  any(is.na(mac_meta$cell)) ||
  any(mac_meta$cell != colnames(counts_mac))
) {
  stop("Macrophage metadata could not be aligned to the count matrix.")
}

sample_counts <- aggregate_sparse_counts(
  counts_mac,
  mac_meta$sample_accession
)
sample_logcpm <- logcpm_from_counts(sample_counts)

sample_order <- sample_meta$sample_accession
if (!all(sample_order %in% colnames(sample_logcpm))) {
  stop("Not all locked samples are present in macrophage pseudobulk matrix.")
}
sample_logcpm <- sample_logcpm[
  ,
  sample_order,
  drop = FALSE
]

state_group <- paste(
  mac_meta$sample_accession,
  mac_meta$macrophage_state,
  sep = "__"
)
state_group_counts <- data.table::data.table(
  state_group = state_group
)[, .N, by = state_group]
eligible_state_groups <- state_group_counts[
  N >= MIN_STATE_SAMPLE_CELLS,
  state_group
]

state_keep <- state_group %in% eligible_state_groups
state_counts <- aggregate_sparse_counts(
  counts_mac[, state_keep, drop = FALSE],
  state_group[state_keep]
)
state_logcpm <- logcpm_from_counts(state_counts)

pseudobulk_audit <- data.table::data.table(
  item = c(
    "macrophage_genes",
    "macrophage_cells",
    "sample_pseudobulks",
    "eligible_sample_state_pseudobulks",
    "minimum_cells_per_sample_state"
  ),
  value = c(
    nrow(counts_mac),
    ncol(counts_mac),
    ncol(sample_counts),
    ncol(state_counts),
    MIN_STATE_SAMPLE_CELLS
  )
)

data.table::fwrite(
  pseudobulk_audit,
  file.path(DIRS$tables, "04_stage4_macrophage_pseudobulk_audit.csv")
)

saveRDS(
  list(
    sample_counts = sample_counts,
    sample_logcpm = sample_logcpm,
    state_counts = state_counts,
    state_logcpm = state_logcpm
  ),
  file.path(
    DIRS$objects,
    "GSE236585_stage4_macrophage_pseudobulk_matrices.rds"
  )
)

############################################################
## 7. Load a broad mouse TF-target network
############################################################

load_dorothea_network <- function() {
  if (!requireNamespace("dorothea", quietly = TRUE)) {
    return(NULL)
  }

  env <- new.env(parent = emptyenv())
  loaded <- tryCatch(
    {
      utils::data(
        "dorothea_mm",
        package = "dorothea",
        envir = env
      )
      TRUE
    },
    error = function(e) FALSE
  )

  if (!loaded || !exists("dorothea_mm", envir = env)) {
    return(NULL)
  }

  net <- data.table::as.data.table(
    get("dorothea_mm", envir = env)
  )
  required <- c("tf", "target", "mor")
  if (!all(required %in% names(net))) return(NULL)

  if (!"confidence" %in% names(net)) {
    net[, confidence := "U"]
  }
  net <- net[confidence %in% DOROTHEA_CONFIDENCE_LEVELS]
  if (nrow(net) == 0L) return(NULL)

  confidence_weight <- c(
    A = 1.00,
    B = 0.85,
    C = 0.70,
    D = 0.50,
    E = 0.30,
    U = 0.50
  )

  net[, source_symbol := as.character(tf)]
  net[, target_symbol := as.character(target)]
  net[, mor := as.numeric(mor)]
  net[, weight := unname(confidence_weight[as.character(confidence)])]
  net[!is.finite(weight), weight := 0.50]
  net[, network_source := "DoRothEA_mouse_ABC"]

  net[
    ,
    .(
      source_symbol,
      target_symbol,
      mor,
      weight,
      confidence,
      network_source
    )
  ]
}

load_collectri_network <- function() {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    return(NULL)
  }

  fun <- tryCatch(
    get("get_collectri", envir = asNamespace("decoupleR")),
    error = function(e) NULL
  )
  if (is.null(fun)) return(NULL)

  formals_names <- names(formals(fun))
  args <- list(organism = "mouse")
  if ("split_complexes" %in% formals_names) {
    args$split_complexes <- FALSE
  }

  net <- tryCatch(
    data.table::as.data.table(do.call(fun, args)),
    error = function(e) {
      add_warning(
        "REGULON_NETWORK",
        "CollecTRI",
        conditionMessage(e)
      )
      NULL
    }
  )
  if (is.null(net) || nrow(net) == 0L) return(NULL)

  required <- c("source", "target", "mor")
  if (!all(required %in% names(net))) return(NULL)

  net[, source_symbol := as.character(source)]
  net[, target_symbol := as.character(target)]
  net[, mor := as.numeric(mor)]
  if ("likelihood" %in% names(net)) {
    net[, weight := abs(as.numeric(likelihood))]
  } else {
    net[, weight := 1]
  }
  net[!is.finite(weight) | weight <= 0, weight := 1]
  net[, confidence := "CollecTRI"]
  net[, network_source := "CollecTRI_mouse"]

  net[
    ,
    .(
      source_symbol,
      target_symbol,
      mor,
      weight,
      confidence,
      network_source
    )
  ]
}

get_mouse_tf_symbols <- function(features) {
  tf_symbols <- character()

  if (
    requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("org.Mm.eg.db", quietly = TRUE)
  ) {
    go_terms <- c(
      "GO:0003700",
      "GO:0000981",
      "GO:0001227",
      "GO:0001228"
    )

    tf_annot <- tryCatch(
      suppressMessages(
        AnnotationDbi::select(
          org.Mm.eg.db::org.Mm.eg.db,
          keys = go_terms,
          columns = "SYMBOL",
          keytype = "GOALL"
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(tf_annot) && "SYMBOL" %in% names(tf_annot)) {
      tf_symbols <- unique(
        as.character(tf_annot$SYMBOL)
      )
      tf_symbols <- tf_symbols[
        !is.na(tf_symbols) & nzchar(tf_symbols)
      ]
    }
  }

  ## Last-resort embedded immune/cardiovascular TF vocabulary. This is
  ## used only if annotation resources are unavailable.
  if (length(tf_symbols) < 20L) {
    tf_symbols <- unique(c(
      "Nfkb1", "Nfkb2", "Rela", "Relb", "Rel", "Irf1", "Irf2",
      "Irf3", "Irf4", "Irf5", "Irf7", "Irf8", "Irf9", "Stat1",
      "Stat2", "Stat3", "Stat4", "Stat5a", "Stat5b", "Stat6",
      "Spi1", "Cebpa", "Cebpb", "Cebpd", "Jun", "Junb", "Jund",
      "Fos", "Fosb", "Fosl1", "Fosl2", "Atf2", "Atf3", "Atf4",
      "Atf6", "Klf2", "Klf4", "Klf5", "Klf6", "Klf9", "Klf10",
      "Runx1", "Runx2", "Runx3", "Maf", "Mafb", "Nfe2l2",
      "Hif1a", "Hif3a", "Egr1", "Egr2", "Egr3", "Myc", "Max",
      "Bcl3", "Bcl6", "Prdm1", "Prdm2", "Prdm8", "Batf", "Batf2",
      "Batf3", "Tfeb", "Tfe3", "Mitf", "Foxo1", "Foxo3", "Foxp1",
      "Foxp3", "Srebf1", "Srebf2", "Pparg", "Ppara", "Ppard",
      "Rxra", "Rxrb", "Nr1h3", "Nr1h2", "Nr4a1", "Nr4a2", "Nr4a3",
      "Smad1", "Smad2", "Smad3", "Smad4", "Smad5", "Smad7",
      "Sox4", "Sox5", "Sox9", "Yy1", "Gata2", "Gata3", "Gata6",
      "Ets1", "Ets2", "Etv3", "Etv5", "Fli1", "Elk1", "Elk3",
      "Creb1", "Creb3", "Xbp1", "Chop", "Ddit3"
    ))
  }

  fmap <- make_feature_map(features)
  mapped <- map_symbols_to_features(tf_symbols, fmap)
  unique(mapped[!is.na(feature), feature])
}

build_correlation_fallback <- function(
  expression_profiles,
  tf_features
) {
  if (ncol(expression_profiles) < FALLBACK_MIN_PROFILES) {
    stop(
      "Too few pseudobulk profiles for correlation fallback: ",
      ncol(expression_profiles)
    )
  }

  tf_features <- intersect(tf_features, rownames(expression_profiles))
  if (length(tf_features) < 10L) {
    stop("Too few detected TFs for correlation fallback.")
  }

  records <- list()
  for (tf in tf_features) {
    tf_values <- as.numeric(expression_profiles[tf, ])
    correlations <- apply(
      expression_profiles,
      1L,
      function(gene_values) {
        suppressWarnings(
          stats::cor(
            tf_values,
            as.numeric(gene_values),
            method = "spearman",
            use = "pairwise.complete.obs"
          )
        )
      }
    )
    correlations[!is.finite(correlations)] <- 0
    correlations[tf] <- 0

    order_idx <- order(abs(correlations), decreasing = TRUE)
    selected <- names(correlations)[
      order_idx[
        seq_len(
          min(FALLBACK_TOP_TARGETS_PER_TF, length(order_idx))
        )
      ]
    ]
    selected <- selected[
      abs(correlations[selected]) >= FALLBACK_MIN_ABS_COR
    ]

    if (length(selected) < MIN_TARGETS_PER_REGULON) {
      selected <- names(correlations)[
        order_idx[
          seq_len(
            min(
              max(MIN_TARGETS_PER_REGULON, 30L),
              length(order_idx)
            )
          )
        ]
      ]
    }

    records[[tf]] <- data.table::data.table(
      source_feature = tf,
      target_feature = selected,
      mor = sign(correlations[selected]),
      weight = abs(correlations[selected]),
      confidence = "data_driven",
      network_source =
        "sample_state_pseudobulk_spearman_fallback"
    )
  }

  data.table::rbindlist(records, use.names = TRUE, fill = TRUE)
}

feature_map <- make_feature_map(rownames(macrophage))
network_raw <- load_dorothea_network()

if (is.null(network_raw) || nrow(network_raw) == 0L) {
  network_raw <- load_collectri_network()
}

network_mode <- "PRIOR_NETWORK"
if (!is.null(network_raw) && nrow(network_raw) > 0L) {
  source_map <- map_symbols_to_features(
    network_raw$source_symbol,
    feature_map
  )
  target_map <- map_symbols_to_features(
    network_raw$target_symbol,
    feature_map
  )

  source_map <- unique(
    source_map[
      ,
      .(source_symbol = symbol, source_feature = feature)
    ],
    by = "source_symbol"
  )
  target_map <- unique(
    target_map[
      ,
      .(target_symbol = symbol, target_feature = feature)
    ],
    by = "target_symbol"
  )

  network_dt <- merge(
    network_raw,
    source_map,
    by = "source_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  network_dt <- merge(
    network_dt,
    target_map,
    by = "target_symbol",
    all.x = TRUE,
    sort = FALSE
  )
} else {
  network_mode <- "DATA_DRIVEN_FALLBACK"
  tf_features <- get_mouse_tf_symbols(rownames(macrophage))
  fallback_profiles <- if (
    ncol(state_logcpm) >= FALLBACK_MIN_PROFILES
  ) {
    state_logcpm
  } else {
    sample_logcpm
  }

  network_dt <- build_correlation_fallback(
    fallback_profiles,
    tf_features
  )

  network_dt[, source_symbol := source_feature]
  network_dt[, target_symbol := target_feature]
}

network_dt <- data.table::as.data.table(network_dt)
required_network_cols <- c(
  "source_symbol",
  "target_symbol",
  "source_feature",
  "target_feature",
  "mor",
  "weight",
  "network_source"
)
missing_network_cols <- setdiff(required_network_cols, names(network_dt))
if (length(missing_network_cols) > 0L) {
  stop(
    "Constructed network is missing column(s): ",
    paste(missing_network_cols, collapse = ", ")
  )
}

network_dt <- network_dt[
  !is.na(source_feature) &
    !is.na(target_feature) &
    source_feature %in% rownames(macrophage) &
    target_feature %in% rownames(macrophage) &
    source_feature != target_feature &
    is.finite(mor) &
    mor != 0 &
    is.finite(weight) &
    weight > 0
]

normalized_mac <- get_assay_matrix(
  macrophage,
  layer = "data",
  assay = "RNA"
)

target_detection <- Matrix::rowMeans(normalized_mac > 0)
source_detection <- target_detection

network_dt <- network_dt[
  target_feature %in%
    names(target_detection)[
      target_detection >= MIN_TARGET_DETECTION_FRACTION
    ] &
    source_feature %in%
    names(source_detection)[
      source_detection >= MIN_TF_DETECTION_FRACTION
    ]
]

## Deduplicate identical TF-target edges and retain strongest evidence.
data.table::setorder(
  network_dt,
  source_symbol,
  target_feature,
  -weight
)
network_dt <- network_dt[
  ,
  .SD[1L],
  by = .(source_symbol, target_feature)
]

network_dt[, rank_within_regulon := data.table::frank(
  -weight,
  ties.method = "first"
), by = source_symbol]
network_dt <- network_dt[
  rank_within_regulon <= MAX_TARGETS_PER_REGULON
]

regulon_sizes <- network_dt[
  ,
  .(
    regulon_size = data.table::uniqueN(target_feature),
    activator_targets = data.table::uniqueN(
      target_feature[mor > 0]
    ),
    repressor_targets = data.table::uniqueN(
      target_feature[mor < 0]
    ),
    source_feature = source_feature[1L],
    network_source = network_source[1L]
  ),
  by = source_symbol
]
valid_tfs <- regulon_sizes[
  regulon_size >= MIN_TARGETS_PER_REGULON,
  source_symbol
]
network_dt <- network_dt[source_symbol %in% valid_tfs]
regulon_sizes <- regulon_sizes[source_symbol %in% valid_tfs]

if (length(valid_tfs) < 10L) {
  stop(
    "Too few valid TF regulons after filtering: ",
    length(valid_tfs)
  )
}

## Compatibility names used by later legacy stages.
network_dt[, regulatoryGene := source_feature]
network_dt[, targetGene := target_feature]
network_dt[, tf_symbol := source_symbol]

write_csv_safe(
  network_dt,
  file.path(DIRS$tables, "05_stage4_full_TF_target_links.csv"),
  compress = TRUE
)
write_csv_safe(
  network_dt,
  file.path(DIRS$tables, "06_stage4_top_regulon_targets.csv")
)
write_csv_safe(
  regulon_sizes,
  file.path(DIRS$tables, "07_stage4_regulon_size_summary.csv")
)

regulon_list <- split(
  network_dt$target_feature,
  network_dt$source_symbol
)
regulon_list <- lapply(regulon_list, unique)
saveRDS(
  regulon_list,
  file.path(DIRS$objects, "stage4_regulon_target_list.rds")
)

network_audit <- data.table::data.table(
  network_mode = network_mode,
  network_source = paste(
    unique(network_dt$network_source),
    collapse = ";"
  ),
  raw_or_filtered_edges = nrow(network_dt),
  valid_regulons = length(valid_tfs),
  median_regulon_size = stats::median(
    regulon_sizes$regulon_size
  ),
  min_regulon_size = min(regulon_sizes$regulon_size),
  max_regulon_size = max(regulon_sizes$regulon_size),
  Nfkb1_forced = FALSE,
  Nfkb1_present = "NFKB1" %in% gene_key(valid_tfs)
)

data.table::fwrite(
  network_audit,
  file.path(DIRS$tables, "08_stage4_regulon_network_audit.csv")
)

log_msg(
  "Regulon network ready: mode=",
  network_mode,
  "; TFs=",
  length(valid_tfs),
  "; edges=",
  nrow(network_dt)
)

############################################################
## 8. Signed weighted regulon activity
############################################################

weighted_sample_activity <- weighted_regulon_activity(
  sample_logcpm,
  network_dt
)

sample_activity_stats_weighted <- sample_level_activity_test(
  weighted_sample_activity,
  sample_meta,
  "signed_weighted_target_zscore"
)

## Cell-level activity is descriptive and is used for localization only.
cell_target_features <- intersect(
  unique(network_dt$target_feature),
  rownames(normalized_mac)
)
cell_activity <- weighted_regulon_activity(
  normalized_mac[cell_target_features, , drop = FALSE],
  network_dt
)

saveRDS(
  weighted_sample_activity,
  file.path(
    DIRS$objects,
    "stage4_weighted_regulon_activity_sample_matrix.rds"
  )
)
saveRDS(
  cell_activity,
  file.path(
    DIRS$objects,
    "stage4_weighted_regulon_activity_cell_matrix.rds"
  )
)

write_csv_safe(
  sample_activity_stats_weighted,
  file.path(
    DIRS$tables,
    "09_stage4_weighted_regulon_activity_HFpEF_vs_Control.csv"
  )
)

############################################################
## 9. AUCell sensitivity analysis
############################################################

aucell_cell_activity <- NULL
aucell_sample_activity <- NULL
sample_activity_stats_aucell <- data.table::data.table()
aucell_status <- "NOT_RUN"

run_aucell_signed <- function(count_matrix, network) {
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    return(NULL)
  }

  positive_sets <- split(
    network[mor > 0, target_feature],
    network[mor > 0, source_symbol]
  )
  negative_sets <- split(
    network[mor < 0, target_feature],
    network[mor < 0, source_symbol]
  )

  positive_sets <- lapply(positive_sets, unique)
  negative_sets <- lapply(negative_sets, unique)
  positive_sets <- positive_sets[
    lengths(positive_sets) >= 5L
  ]
  negative_sets <- negative_sets[
    lengths(negative_sets) >= 5L
  ]

  gene_sets <- c(
    stats::setNames(
      positive_sets,
      paste0(names(positive_sets), "__POS")
    ),
    stats::setNames(
      negative_sets,
      paste0(names(negative_sets), "__NEG")
    )
  )

  if (length(gene_sets) < 5L) return(NULL)

  build_fun <- get(
    "AUCell_buildRankings",
    envir = asNamespace("AUCell")
  )
  build_args <- list(count_matrix)
  build_formals <- names(formals(build_fun))
  if ("plotStats" %in% build_formals) {
    build_args$plotStats <- FALSE
  }
  if ("verbose" %in% build_formals) {
    build_args$verbose <- FALSE
  }
  if ("nCores" %in% build_formals) {
    build_args$nCores <- 1L
  }

  rankings <- do.call(build_fun, build_args)

  auc_fun <- get(
    "AUCell_calcAUC",
    envir = asNamespace("AUCell")
  )
  auc_args <- list(
    geneSets = gene_sets,
    rankings = rankings
  )
  auc_formals <- names(formals(auc_fun))
  if ("aucMaxRank" %in% auc_formals) {
    auc_args$aucMaxRank <- max(
      50L,
      ceiling(0.05 * nrow(count_matrix))
    )
  }
  if ("verbose" %in% auc_formals) {
    auc_args$verbose <- FALSE
  }
  if ("nCores" %in% auc_formals) {
    auc_args$nCores <- 1L
  }

  auc_result <- do.call(auc_fun, auc_args)
  auc_matrix <- as.matrix(AUCell::getAUC(auc_result))

  tfs <- sort(unique(network$source_symbol))
  net_matrix <- matrix(
    NA_real_,
    nrow = length(tfs),
    ncol = ncol(auc_matrix),
    dimnames = list(tfs, colnames(auc_matrix))
  )

  for (tf in tfs) {
    pos_name <- paste0(tf, "__POS")
    neg_name <- paste0(tf, "__NEG")
    pos <- if (pos_name %in% rownames(auc_matrix)) {
      auc_matrix[pos_name, ]
    } else {
      rep(0, ncol(auc_matrix))
    }
    neg <- if (neg_name %in% rownames(auc_matrix)) {
      auc_matrix[neg_name, ]
    } else {
      rep(0, ncol(auc_matrix))
    }
    net_matrix[tf, ] <- pos - neg
  }

  net_matrix <- net_matrix[
    rowSums(is.finite(net_matrix)) > 0,
    ,
    drop = FALSE
  ]
  net_matrix
}

if (RUN_AUCELL_SENSITIVITY) {
  aucell_cell_activity <- tryCatch(
    {
      log_msg("Running signed AUCell sensitivity analysis.")
      run_aucell_signed(counts_mac, network_dt)
    },
    error = function(e) {
      add_warning(
        "AUCELL",
        "signed_activity",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(aucell_cell_activity) && nrow(aucell_cell_activity) > 0L) {
    aucell_status <- "COMPLETED"
    aucell_sample_activity <- matrix(
      NA_real_,
      nrow = nrow(aucell_cell_activity),
      ncol = nrow(sample_meta),
      dimnames = list(
        rownames(aucell_cell_activity),
        sample_meta$sample_accession
      )
    )

    for (sample_id in sample_meta$sample_accession) {
      sample_cells <- mac_meta[
        sample_accession == sample_id,
        cell
      ]
      sample_cells <- intersect(
        sample_cells,
        colnames(aucell_cell_activity)
      )
      if (length(sample_cells) > 0L) {
        aucell_sample_activity[, sample_id] <- rowMeans(
          aucell_cell_activity[, sample_cells, drop = FALSE],
          na.rm = TRUE
        )
      }
    }

    sample_activity_stats_aucell <- sample_level_activity_test(
      aucell_sample_activity,
      sample_meta,
      "signed_AUCell"
    )

    saveRDS(
      aucell_cell_activity,
      file.path(
        DIRS$objects,
        "stage4_AUCell_regulon_activity_cell_matrix.rds"
      )
    )
    saveRDS(
      aucell_sample_activity,
      file.path(
        DIRS$objects,
        "stage4_AUCell_regulon_activity_sample_matrix.rds"
      )
    )
    write_csv_safe(
      sample_activity_stats_aucell,
      file.path(
        DIRS$tables,
        "10_stage4_AUCell_regulon_activity_HFpEF_vs_Control.csv"
      )
    )
  } else {
    aucell_status <- "UNAVAILABLE_OR_FAILED_NONFATAL"
  }
}

############################################################
## 10. TF-expression-only comparator
############################################################

source_map_unique <- network_dt[
  ,
  .(source_feature = source_feature[1L]),
  by = source_symbol
]
source_map_unique <- source_map_unique[
  source_feature %in% rownames(sample_logcpm)
]

tf_expression_matrix <- sample_logcpm[
  source_map_unique$source_feature,
  ,
  drop = FALSE
]
rownames(tf_expression_matrix) <- source_map_unique$source_symbol

tf_expression_stats <- sample_level_activity_test(
  tf_expression_matrix,
  sample_meta,
  "TF_expression_logCPM"
)

write_csv_safe(
  tf_expression_stats,
  file.path(
    DIRS$tables,
    "11_stage4_TF_expression_HFpEF_vs_Control.csv"
  )
)

############################################################
## 11. Sample-by-state activity localization
############################################################

cell_activity_dt <- data.table::as.data.table(
  t(cell_activity),
  keep.rownames = "cell"
)
cell_activity_long <- data.table::melt(
  cell_activity_dt,
  id.vars = "cell",
  variable.name = "tf_symbol",
  value.name = "weighted_activity"
)
cell_activity_long <- merge(
  cell_activity_long,
  mac_meta[
    ,
    .(
      cell,
      sample_accession,
      condition,
      macrophage_state
    )
  ],
  by = "cell",
  all.x = TRUE
)

state_activity_summary <- cell_activity_long[
  ,
  .(
    cell_count = .N,
    sample_state_mean_activity = mean(
      weighted_activity,
      na.rm = TRUE
    )
  ),
  by = .(
    tf_symbol,
    sample_accession,
    condition,
    macrophage_state
  )
]
state_activity_summary <- state_activity_summary[
  cell_count >= MIN_STATE_SAMPLE_CELLS
]

state_activity_tests <- state_activity_summary[
  ,
  {
    control <- sample_state_mean_activity[
      condition == "Control"
    ]
    hfpef <- sample_state_mean_activity[
      condition == "HFpEF"
    ]
    .(
      control_mean = if (length(control) > 0L) {
        mean(control, na.rm = TRUE)
      } else {
        NA_real_
      },
      hfpef_mean = if (length(hfpef) > 0L) {
        mean(hfpef, na.rm = TRUE)
      } else {
        NA_real_
      },
      hfpef_minus_control = if (
        length(control) > 0L &&
        length(hfpef) > 0L
      ) {
        mean(hfpef, na.rm = TRUE) -
          mean(control, na.rm = TRUE)
      } else {
        NA_real_
      },
      hedges_g_HFpEF_vs_Control = hedges_g(hfpef, control),
      wilcoxon_p = safe_wilcox_p(hfpef, control),
      control_samples = length(control),
      hfpef_samples = length(hfpef)
    )
  },
  by = .(tf_symbol, macrophage_state)
]
state_activity_tests[, wilcoxon_fdr := p.adjust(
  wilcoxon_p,
  method = "BH"
)]

write_csv_safe(
  state_activity_summary,
  file.path(
    DIRS$tables,
    "12_stage4_TF_activity_by_sample_and_macrophage_state.csv"
  ),
  compress = TRUE
)
write_csv_safe(
  state_activity_tests,
  file.path(
    DIRS$tables,
    "13_stage4_TF_activity_state_HFpEF_vs_Control.csv"
  )
)

############################################################
## 12. Leave-one-Control/leave-one-HFpEF-pair robustness
############################################################

full_weighted_stats <- data.table::copy(
  sample_activity_stats_weighted
)
full_effect_map <- full_weighted_stats[
  ,
  .(tf_symbol, full_effect = hfpef_minus_control)
]

control_samples <- sample_meta[
  condition == "Control",
  sample_accession
]
hfpef_samples <- sample_meta[
  condition == "HFpEF",
  sample_accession
]

loo_records <- list()
for (control_removed in control_samples) {
  for (hfpef_removed in hfpef_samples) {
    keep_samples <- setdiff(
      sample_meta$sample_accession,
      c(control_removed, hfpef_removed)
    )
    keep_meta <- sample_meta[
      sample_accession %in% keep_samples
    ]
    keep_meta <- keep_meta[
      match(keep_samples, sample_accession)
    ]

    effects <- vapply(
      rownames(weighted_sample_activity),
      function(tf) {
        values <- weighted_sample_activity[
          tf,
          keep_samples
        ]
        mean(
          values[keep_meta$condition == "HFpEF"],
          na.rm = TRUE
        ) - mean(
          values[keep_meta$condition == "Control"],
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    ranks <- rank(-abs(effects), ties.method = "average")
    loo_records[[length(loo_records) + 1L]] <-
      data.table::data.table(
        tf_symbol = names(effects),
        removed_control = control_removed,
        removed_hfpef = hfpef_removed,
        loo_effect = as.numeric(effects),
        loo_abs_effect_rank = as.numeric(ranks)
      )
  }
}

loo_results <- data.table::rbindlist(loo_records)
loo_results <- merge(
  loo_results,
  full_effect_map,
  by = "tf_symbol",
  all.x = TRUE
)
loo_results[, sign_consistent_with_full := (
  sign(loo_effect) == sign(full_effect)
)]

loo_summary <- loo_results[
  ,
  .(
    leave_pair_runs = .N,
    sign_stability = mean(
      sign_consistent_with_full,
      na.rm = TRUE
    ),
    median_loo_effect = stats::median(
      loo_effect,
      na.rm = TRUE
    ),
    min_loo_effect = min(loo_effect, na.rm = TRUE),
    max_loo_effect = max(loo_effect, na.rm = TRUE),
    median_abs_effect_rank = stats::median(
      loo_abs_effect_rank,
      na.rm = TRUE
    ),
    top10_frequency = mean(
      loo_abs_effect_rank <= 10,
      na.rm = TRUE
    ),
    top20_frequency = mean(
      loo_abs_effect_rank <= 20,
      na.rm = TRUE
    )
  ),
  by = tf_symbol
]

write_csv_safe(
  loo_results,
  file.path(
    DIRS$tables,
    "14_stage4_leave_one_pair_out_TF_activity_results.csv"
  ),
  compress = TRUE
)
write_csv_safe(
  loo_summary,
  file.path(
    DIRS$tables,
    "15_stage4_leave_one_pair_out_TF_robustness_summary.csv"
  )
)

############################################################
## 13. Regulon overlap with Stage 2 and Stage 3-supported programs
############################################################

program_sets <- split(
  program_manifest_all$symbol_key,
  program_manifest_all$program_name
)
program_sets <- lapply(program_sets, unique)

network_dt[, target_key := gene_key(target_symbol)]
network_universe <- unique(network_dt$target_key)
network_universe <- network_universe[nzchar(network_universe)]

regulon_overlap_records <- list()
for (tf in sort(unique(network_dt$source_symbol))) {
  tf_targets <- unique(
    network_dt[
      source_symbol == tf,
      target_key
    ]
  )
  tf_targets <- intersect(tf_targets, network_universe)

  for (program_name in names(program_sets)) {
    program_genes <- intersect(
      program_sets[[program_name]],
      network_universe
    )
    overlap_genes <- intersect(tf_targets, program_genes)

    hypergeom_p <- if (
      length(program_genes) > 0L &&
      length(tf_targets) > 0L
    ) {
      stats::phyper(
        q = length(overlap_genes) - 1L,
        m = length(program_genes),
        n = length(network_universe) - length(program_genes),
        k = length(tf_targets),
        lower.tail = FALSE
      )
    } else {
      NA_real_
    }

    regulon_overlap_records[[
      length(regulon_overlap_records) + 1L
    ]] <- data.table::data.table(
      tf_symbol = tf,
      program_name = program_name,
      regulon_size = length(tf_targets),
      program_genes_in_universe = length(program_genes),
      overlap_count = length(overlap_genes),
      overlap_fraction_of_regulon = if (
        length(tf_targets) > 0L
      ) {
        length(overlap_genes) / length(tf_targets)
      } else {
        NA_real_
      },
      overlap_fraction_of_program = if (
        length(program_genes) > 0L
      ) {
        length(overlap_genes) / length(program_genes)
      } else {
        NA_real_
      },
      overlap_genes = paste(overlap_genes, collapse = ";"),
      hypergeom_p = hypergeom_p
    )
  }
}

regulon_overlap <- data.table::rbindlist(
  regulon_overlap_records,
  use.names = TRUE,
  fill = TRUE
)
regulon_overlap[, hypergeom_fdr := p.adjust(
  hypergeom_p,
  method = "BH"
), by = program_name]
regulon_overlap[, neglog10_fdr := safe_neglog10(hypergeom_fdr)]

write_csv_safe(
  regulon_overlap,
  file.path(
    DIRS$tables,
    "16_stage4_regulon_overlap_with_stage2_stage3_programs.csv"
  ),
  compress = TRUE
)

## Select a primary overlap set without forcing any TF.
preferred_primary_programs <- c(
  paste0(
    "Ccr2pos_Disease_up_Drug_down_Top",
    PRIMARY_SIGNATURE_SIZE,
    "_Stage3Supported"
  ),
  paste0(
    "CrossSubset_Disease_up_Drug_down_Top",
    PRIMARY_SIGNATURE_SIZE,
    "_Stage3Supported"
  ),
  paste0(
    "Ccr2pos_Disease_up_Drug_down_Top",
    PRIMARY_SIGNATURE_SIZE
  ),
  paste0(
    "CrossSubset_Disease_up_Drug_down_Top",
    PRIMARY_SIGNATURE_SIZE
  )
)

available_programs <- unique(regulon_overlap$program_name)
primary_overlap_program <- preferred_primary_programs[
  preferred_primary_programs %in% available_programs
][1L]
if (
  length(primary_overlap_program) == 0L ||
  is.na(primary_overlap_program)
) {
  primary_overlap_program <- available_programs[1L]
}

primary_overlap <- regulon_overlap[
  program_name == primary_overlap_program
]
primary_overlap <- merge(
  primary_overlap,
  unique(
    network_dt[
      ,
      .(
        tf_symbol = source_symbol,
        regulatoryGene = source_feature
      )
    ],
    by = "tf_symbol"
  ),
  by = "tf_symbol",
  all.x = TRUE
)
primary_overlap[, target_program := primary_overlap_program]
primary_overlap[, compatibility_alias_only := TRUE]

## Compatibility filename retained for rewritten Stage 5 integration.
write_csv_safe(
  primary_overlap,
  file.path(
    DIRS$tables,
    "11_stage4_regulon_overlap_with_stage3b_clean_genes.csv"
  )
)

############################################################
## 14. Method comparison and composite TF priority
############################################################

weighted_priority <- sample_activity_stats_weighted[
  ,
  .(
    tf_symbol,
    weighted_effect = hfpef_minus_control,
    weighted_hedges_g = hedges_g_HFpEF_vs_Control,
    weighted_limma_padj = limma_padj,
    weighted_wilcoxon_fdr = wilcoxon_fdr
  )
]

priority_dt <- merge(
  weighted_priority,
  tf_expression_stats[
    ,
    .(
      tf_symbol,
      expression_effect = hfpef_minus_control,
      expression_hedges_g = hedges_g_HFpEF_vs_Control,
      expression_limma_padj = limma_padj
    )
  ],
  by = "tf_symbol",
  all.x = TRUE
)

priority_dt <- merge(
  priority_dt,
  loo_summary,
  by = "tf_symbol",
  all.x = TRUE
)

max_overlap <- regulon_overlap[
  ,
  .(
    max_overlap_count = max(overlap_count, na.rm = TRUE),
    max_overlap_fraction_program = max(
      overlap_fraction_of_program,
      na.rm = TRUE
    ),
    max_overlap_neglog10_fdr = max(
      neglog10_fdr,
      na.rm = TRUE
    ),
    best_overlap_program = program_name[
      which.max(neglog10_fdr)
    ]
  ),
  by = tf_symbol
]
priority_dt <- merge(
  priority_dt,
  max_overlap,
  by = "tf_symbol",
  all.x = TRUE
)

supported_overlap <- regulon_overlap[
  grepl("Stage3Supported$", program_name)
][
  ,
  .(
    supported_overlap_count = max(overlap_count, na.rm = TRUE),
    supported_overlap_neglog10_fdr = max(
      neglog10_fdr,
      na.rm = TRUE
    )
  ),
  by = tf_symbol
]
priority_dt <- merge(
  priority_dt,
  supported_overlap,
  by = "tf_symbol",
  all.x = TRUE
)

if (
  aucell_status == "COMPLETED" &&
  nrow(sample_activity_stats_aucell) > 0L
) {
  priority_dt <- merge(
    priority_dt,
    sample_activity_stats_aucell[
      ,
      .(
        tf_symbol,
        aucell_effect = hfpef_minus_control,
        aucell_hedges_g = hedges_g_HFpEF_vs_Control,
        aucell_limma_padj = limma_padj
      )
    ],
    by = "tf_symbol",
    all.x = TRUE
  )
  priority_dt[, activity_method_sign_agreement := (
    sign(weighted_effect) == sign(aucell_effect)
  )]
} else {
  priority_dt[, aucell_effect := NA_real_]
  priority_dt[, aucell_hedges_g := NA_real_]
  priority_dt[, aucell_limma_padj := NA_real_]
  priority_dt[, activity_method_sign_agreement := NA]
}

priority_dt[!is.finite(max_overlap_count), max_overlap_count := 0]
priority_dt[
  !is.finite(max_overlap_fraction_program),
  max_overlap_fraction_program := 0
]
priority_dt[
  !is.finite(max_overlap_neglog10_fdr),
  max_overlap_neglog10_fdr := 0
]
priority_dt[
  !is.finite(supported_overlap_count),
  supported_overlap_count := 0
]
priority_dt[
  !is.finite(supported_overlap_neglog10_fdr),
  supported_overlap_neglog10_fdr := 0
]
priority_dt[!is.finite(sign_stability), sign_stability := 0]
priority_dt[!is.finite(top20_frequency), top20_frequency := 0]

priority_dt[, score_activity_effect := rescale01(
  abs(weighted_hedges_g)
)]
priority_dt[, score_activity_evidence := rescale01(
  safe_neglog10(weighted_limma_padj)
)]
priority_dt[, score_method_agreement := data.table::fcase(
  activity_method_sign_agreement == TRUE,
  1,
  activity_method_sign_agreement == FALSE,
  0,
  default = 0.5
)]
priority_dt[, score_overlap := rescale01(
  max_overlap_neglog10_fdr +
    log1p(max_overlap_count)
)]
priority_dt[, score_supported_overlap := rescale01(
  supported_overlap_neglog10_fdr +
    log1p(supported_overlap_count)
)]
priority_dt[, score_robustness := pmin(
  1,
  0.65 * sign_stability +
    0.35 * top20_frequency
)]
priority_dt[, score_expression := rescale01(
  abs(expression_hedges_g)
)]

priority_dt[, priority_score :=
  0.25 * score_activity_effect +
  0.15 * score_activity_evidence +
  0.10 * score_method_agreement +
  0.20 * score_overlap +
  0.10 * score_supported_overlap +
  0.15 * score_robustness +
  0.05 * score_expression
]

data.table::setorder(
  priority_dt,
  -priority_score,
  -score_robustness,
  -score_overlap
)
priority_dt[, priority_rank := seq_len(.N)]
priority_dt[, Nfkb1_forced := FALSE]

write_csv_safe(
  priority_dt,
  file.path(
    DIRS$tables,
    "12_stage4_candidate_TF_priority_score.csv"
  )
)

method_comparison <- merge(
  weighted_priority,
  tf_expression_stats[
    ,
    .(
      tf_symbol,
      expression_effect = hfpef_minus_control
    )
  ],
  by = "tf_symbol",
  all.x = TRUE
)
if (aucell_status == "COMPLETED") {
  method_comparison <- merge(
    method_comparison,
    sample_activity_stats_aucell[
      ,
      .(
        tf_symbol,
        aucell_effect = hfpef_minus_control
      )
    ],
    by = "tf_symbol",
    all.x = TRUE
  )
} else {
  method_comparison[, aucell_effect := NA_real_]
}

method_summary <- data.table::data.table(
  comparison = c(
    "Weighted activity vs TF expression",
    "Weighted activity vs AUCell",
    "Priority rank vs TF-expression absolute-effect rank"
  ),
  spearman = c(
    safe_spearman(
      method_comparison$weighted_effect,
      method_comparison$expression_effect
    ),
    safe_spearman(
      method_comparison$weighted_effect,
      method_comparison$aucell_effect
    ),
    safe_spearman(
      priority_dt$priority_rank,
      rank(
        -abs(priority_dt$expression_effect),
        ties.method = "average"
      )
    )
  ),
  n_TFs = c(
    sum(
      is.finite(method_comparison$weighted_effect) &
        is.finite(method_comparison$expression_effect)
    ),
    sum(
      is.finite(method_comparison$weighted_effect) &
        is.finite(method_comparison$aucell_effect)
    ),
    nrow(priority_dt)
  )
)

write_csv_safe(
  method_comparison,
  file.path(
    DIRS$tables,
    "17_stage4_TF_method_comparison_per_TF.csv"
  )
)
write_csv_safe(
  method_summary,
  file.path(
    DIRS$tables,
    "18_stage4_TF_method_comparison_summary.csv"
  )
)

############################################################
## 15. Add top TF activities to the macrophage object
############################################################

top_tfs <- head(
  priority_dt$tf_symbol,
  min(TOP_TFS_REPORT, nrow(priority_dt))
)

for (tf in top_tfs) {
  if (tf %in% rownames(cell_activity)) {
    column_name <- paste0(
      "stage4_TF_activity_",
      make.names(tf)
    )
    macrophage@meta.data[[column_name]] <- as.numeric(
      cell_activity[tf, colnames(macrophage)]
    )
  }
}

macrophage$stage4_condition <- factor(
  macrophage$condition,
  levels = c("Control", "HFpEF")
)

saveRDS(
  macrophage,
  file.path(
    DIRS$objects,
    "GSE236585_stage4_macrophage_regulon_scored.rds"
  )
)

## Compatibility object name for later rewritten stages.
saveRDS(
  macrophage,
  file.path(
    DIRS$objects,
    "GSE236585_stage4_macrophage_subset_reprocessed.rds"
  )
)

## Restart checkpoint: all scientific calculations are complete before
## figure/report generation begins.
saveRDS(
  list(
    macrophage = macrophage,
    sample_meta = sample_meta,
    network_dt = network_dt,
    regulon_sizes = regulon_sizes,
    program_manifest_all = program_manifest_all,
    weighted_sample_activity = weighted_sample_activity,
    sample_activity_stats_weighted = sample_activity_stats_weighted,
    aucell_status = aucell_status,
    aucell_sample_activity = aucell_sample_activity,
    sample_activity_stats_aucell = sample_activity_stats_aucell,
    tf_expression_stats = tf_expression_stats,
    state_activity_tests = state_activity_tests,
    loo_summary = loo_summary,
    regulon_overlap = regulon_overlap,
    priority_dt = priority_dt,
    method_comparison = method_comparison,
    method_summary = method_summary,
    network_mode = network_mode
  ),
  file.path(
    DIRS$objects,
    "CHECKPOINT_stage4_scientific_results_pre_figures.rds"
  ),
  compress = FALSE
)

log_msg(
  "Stage 4 scientific calculations completed and checkpointed before figures."
)

############################################################
## 16. Figures and source data
############################################################

condition_palette <- c(
  Control = "#4C78A8",
  HFpEF = "#E45756"
)

## Figure 4A: sample-level weighted activity heatmap.
heat_tfs <- head(
  priority_dt$tf_symbol,
  min(TOP_TFS_HEATMAP, nrow(priority_dt))
)
heat_tfs <- intersect(heat_tfs, rownames(weighted_sample_activity))
heat_matrix <- weighted_sample_activity[
  heat_tfs,
  sample_meta$sample_accession,
  drop = FALSE
]
heat_matrix <- scale_rows(heat_matrix)
annotation_col <- data.frame(
  Condition = as.character(sample_meta$condition),
  row.names = sample_meta$sample_accession,
  stringsAsFactors = FALSE
)

write_csv_safe(
  data.table::as.data.table(
    heat_matrix,
    keep.rownames = "tf_symbol"
  ),
  file.path(
    DIRS$source,
    "Fig4A_sample_level_TF_activity_heatmap_source.csv"
  )
)

save_heatmap_bundle(
  heat_matrix,
  "Fig4A_sample_level_macrophage_TF_activity_heatmap",
  annotation_col = annotation_col,
  width = 9,
  height = 8
)

## Figure 4B: composite priority ranking.
plot_priority <- priority_dt[
  seq_len(min(TOP_TFS_REPORT, .N))
]
plot_priority <- plot_priority[order(priority_score)]
plot_priority[, tf_symbol := factor(
  tf_symbol,
  levels = tf_symbol
)]

write_csv_safe(
  plot_priority,
  file.path(
    DIRS$source,
    "Fig4B_TF_priority_source.csv"
  )
)

p_priority <- ggplot2::ggplot(
  plot_priority,
  ggplot2::aes(
    x = priority_score,
    y = tf_symbol
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = priority_score,
      y = tf_symbol,
      yend = tf_symbol
    ),
    linewidth = 0.5,
    color = "grey65"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = max_overlap_count,
      fill = weighted_effect
    ),
    shape = 21,
    color = "black",
    stroke = 0.3
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0
  ) +
  ggplot2::scale_size_continuous(range = c(2.5, 8)) +
  ggplot2::labs(
    title = "Macrophage TF-regulon consensus prioritization",
    subtitle = paste0(
      "No TF was forced; network mode: ",
      network_mode
    ),
    x = "Composite priority score",
    y = NULL,
    fill = "HFpEF - Control\nTF activity",
    size = "Maximum\nprogram overlap"
  ) +
  ggplot2::theme_bw(base_size = 10)

save_plot_bundle(
  p_priority,
  "Fig4B_macrophage_TF_priority_ranking",
  9,
  8
)

## Figure 4C: overlap landscape.
## Add program metadata before filtering the overlap landscape.
overlap_program_meta <- unique(
  program_manifest_all[
    ,
    .(
      program_name,
      subset_name,
      direction,
      signature_size
    )
  ]
)
overlap_plot_dt <- merge(
  regulon_overlap,
  overlap_program_meta,
  by = "program_name",
  all.x = TRUE
)
overlap_plot_dt <- overlap_plot_dt[
  tf_symbol %in% head(priority_dt$tf_symbol, 15L) &
    (
      signature_size == PRIMARY_SIGNATURE_SIZE |
      grepl("Stage3Supported$", program_name)
    )
]
overlap_plot_dt[, tf_symbol := factor(
  tf_symbol,
  levels = rev(head(priority_dt$tf_symbol, 15L))
)]

write_csv_safe(
  overlap_plot_dt,
  file.path(
    DIRS$source,
    "Fig4C_regulon_program_overlap_source.csv"
  )
)

p_overlap <- ggplot2::ggplot(
  overlap_plot_dt,
  ggplot2::aes(
    x = program_name,
    y = tf_symbol
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = overlap_count,
      fill = neglog10_fdr
    ),
    shape = 21,
    color = "black",
    stroke = 0.25
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#B2182B"
  ) +
  ggplot2::scale_size_continuous(range = c(1.5, 7)) +
  ggplot2::labs(
    title = "Regulon overlap with drug-opposed programs",
    subtitle = "Stage 3-supported subsets are explicitly labelled",
    x = NULL,
    y = NULL,
    fill = "-log10(FDR)",
    size = "Overlap genes"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 55,
      hjust = 1
    )
  )

save_plot_bundle(
  p_overlap,
  "Fig4C_regulon_overlap_with_stage2_stage3_programs",
  13,
  7
)

## Figure 4D: leave-one-pair robustness.
robust_plot_dt <- priority_dt[
  seq_len(min(TOP_TFS_REPORT, .N))
]

write_csv_safe(
  robust_plot_dt,
  file.path(
    DIRS$source,
    "Fig4D_leave_one_pair_robustness_source.csv"
  )
)

p_robust <- ggplot2::ggplot(
  robust_plot_dt,
  ggplot2::aes(
    x = sign_stability,
    y = top20_frequency,
    label = tf_symbol
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = abs(weighted_hedges_g),
      fill = priority_score
    ),
    shape = 21,
    color = "black",
    stroke = 0.25
  ) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = Inf
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::scale_fill_gradient(
    low = "#FEE8C8",
    high = "#B30000"
  ) +
  ggplot2::labs(
    title = "Leave-one-Control/leave-one-HFpEF-pair robustness",
    x = "Direction stability across 9 leave-pair-out runs",
    y = "Frequency ranked among top 20 TF effects",
    size = "|Hedges g|",
    fill = "Priority score"
  ) +
  ggplot2::theme_bw(base_size = 10)

save_plot_bundle(
  p_robust,
  "Fig4D_leave_one_pair_out_TF_robustness",
  9,
  7
)

## Figure 4E: top-ranked TF activity UMAP.
top_tf <- priority_dt$tf_symbol[1L]
top_tf_col <- paste0(
  "stage4_TF_activity_",
  make.names(top_tf)
)
if (top_tf_col %in% names(macrophage@meta.data)) {
  p_top_umap <- Seurat::FeaturePlot(
    macrophage,
    features = top_tf_col,
    reduction = "umap",
    raster = TRUE
  ) +
    ggplot2::labs(
      title = paste0(
        "Top-ranked TF activity: ",
        top_tf
      ),
      subtitle = "Cell-level map is descriptive; inference uses six biological samples"
    ) +
    ggplot2::theme_bw(base_size = 10)

  save_plot_bundle(
    p_top_umap,
    "Fig4E_top_ranked_TF_activity_UMAP",
    8,
    6.5
  )
}

## Supplementary method-comparison figure.
if (aucell_status == "COMPLETED") {
  method_plot_dt <- method_comparison[
    is.finite(weighted_effect) &
      is.finite(aucell_effect)
  ]

  write_csv_safe(
    method_plot_dt,
    file.path(
      DIRS$source,
      "FigS4A_weighted_vs_AUCell_source.csv"
    )
  )

  p_method <- ggplot2::ggplot(
    method_plot_dt,
    ggplot2::aes(
      x = weighted_effect,
      y = aucell_effect,
      label = tf_symbol
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point(size = 2.5) +
    ggrepel::geom_text_repel(
      data = method_plot_dt[
        order(-abs(weighted_effect))[1:min(15L, .N)]
      ],
      size = 3,
      max.overlaps = Inf
    ) +
    ggplot2::labs(
      title = "TF-activity method comparison",
      subtitle = paste0(
        "Spearman rho = ",
        round(
          safe_spearman(
            method_plot_dt$weighted_effect,
            method_plot_dt$aucell_effect
          ),
          3
        )
      ),
      x = "Weighted regulon score: HFpEF - Control",
      y = "Signed AUCell: HFpEF - Control"
    ) +
    ggplot2::theme_bw(base_size = 10)

  save_plot_bundle(
    p_method,
    "FigS4A_TF_activity_method_comparison",
    8,
    7
  )
}

## Supplementary state-localization heatmap.
state_top <- state_activity_tests[
  tf_symbol %in% head(priority_dt$tf_symbol, 15L)
]
state_heat <- data.table::dcast(
  state_top,
  tf_symbol ~ macrophage_state,
  value.var = "hfpef_minus_control"
)
if (nrow(state_heat) > 1L && ncol(state_heat) > 2L) {
  state_matrix <- as.matrix(
    state_heat[, -1L, with = FALSE]
  )
  rownames(state_matrix) <- state_heat$tf_symbol
  state_matrix[!is.finite(state_matrix)] <- 0

  write_csv_safe(
    state_heat,
    file.path(
      DIRS$source,
      "FigS4B_TF_activity_by_macrophage_state_source.csv"
    )
  )

  save_heatmap_bundle(
    state_matrix,
    "FigS4B_TF_activity_by_macrophage_state",
    annotation_col = NULL,
    width = 10,
    height = 7
  )
}

############################################################
## 17. Workbook, methods, and parameter documentation
############################################################

workbook_path <- file.path(
  DIRS$tables,
  "19_GSE236585_stage4_TF_regulon_key_results.xlsx"
)
wb <- openxlsx::createWorkbook()

write_sheet_safe(wb, "Input_validation", input_audit)
write_sheet_safe(wb, "Program_sizes", program_size_summary)
write_sheet_safe(wb, "Network_audit", network_audit)
write_sheet_safe(wb, "Regulon_sizes", regulon_sizes)
write_sheet_safe(wb, "Weighted_activity", sample_activity_stats_weighted)
write_sheet_safe(wb, "AUCell_activity", sample_activity_stats_aucell)
write_sheet_safe(wb, "TF_expression", tf_expression_stats)
write_sheet_safe(wb, "LOO_robustness", loo_summary)
write_sheet_safe(wb, "Program_overlap", regulon_overlap)
write_sheet_safe(wb, "TF_priority", priority_dt)
write_sheet_safe(wb, "Method_comparison", method_summary)
write_sheet_safe(wb, "State_activity", state_activity_tests)

openxlsx::saveWorkbook(
  wb,
  workbook_path,
  overwrite = TRUE
)

parameter_table <- data.table::data.table(
  parameter = c(
    "Random seed",
    "Stage 3 macrophage cells",
    "Inferential unit",
    "Primary signature size",
    "Allowed Stage 2 tiers",
    "Minimum targets per regulon",
    "Maximum targets per regulon",
    "Minimum target detection fraction",
    "Minimum TF detection fraction",
    "Minimum cells per sample-state",
    "Primary activity method",
    "Sensitivity activity method",
    "Network mode",
    "AUCell status",
    "Leave-pair-out runs",
    "Nfkb1 forced"
  ),
  value = c(
    "20260714",
    as.character(ncol(macrophage)),
    "Biological sample",
    as.character(PRIMARY_SIGNATURE_SIZE),
    paste(SIGNATURE_TIERS, collapse = "; "),
    as.character(MIN_TARGETS_PER_REGULON),
    as.character(MAX_TARGETS_PER_REGULON),
    as.character(MIN_TARGET_DETECTION_FRACTION),
    as.character(MIN_TF_DETECTION_FRACTION),
    as.character(MIN_STATE_SAMPLE_CELLS),
    "Signed weighted target z-score",
    "Signed AUCell when available",
    network_mode,
    aucell_status,
    as.character(length(control_samples) * length(hfpef_samples)),
    "FALSE"
  ),
  rationale = c(
    "Reproducibility",
    "Use the locked Stage 3 macrophage subset",
    "Avoid cell-level pseudoreplication",
    "Primary cross-stage program while retaining size sensitivity",
    "Exclude weak directional-only Tier D genes",
    "Avoid unstable very small regulons",
    "Limit dominance by very large prior regulons",
    "Remove nearly undetected targets",
    "Require TF detectability in macrophages",
    "Require minimally supported state-level summaries",
    "Always available and direction-aware",
    "Independent rank-based activity sensitivity analysis",
    "Prior network preferred; explicit fallback otherwise",
    "Nonfatal sensitivity method",
    "Every Control-HFpEF removal combination",
    "No candidate TF was manually inserted or promoted"
  )
)

data.table::fwrite(
  parameter_table,
  file.path(DIRS$methods, "stage4_parameters_and_rationale.csv")
)

priority_formula <- c(
  "Stage 4 TF priority score",
  "",
  "priority_score =",
  "  0.25 * activity_effect +",
  "  0.15 * activity_statistical_evidence +",
  "  0.10 * activity_method_sign_agreement +",
  "  0.20 * regulon_program_overlap +",
  "  0.10 * Stage3_supported_program_overlap +",
  "  0.15 * leave_pair_out_robustness +",
  "  0.05 * TF_expression_effect",
  "",
  "All continuous components are rescaled to [0,1].",
  "AUCell unavailability is assigned a neutral method-agreement score of 0.5.",
  "No TF, including Nfkb1, is forced into the network or ranking."
)
writeLines(
  priority_formula,
  file.path(DIRS$methods, "stage4_priority_score_definition.txt"),
  useBytes = TRUE
)

methods_text <- c(
  "HFpEF Stage 4 FIXED v1: macrophage TF-regulon consensus prioritization",
  "",
  "Input boundary:",
  "- The completed Stage 3 GSE236585 macrophage/monocyte Seurat object was loaded from disk.",
  "- Stage 4 did not repeat raw-data import, QC, scDblFinder, clustering, annotation, or Stage 3 pseudobulk analysis.",
  "- The locked input contained 1,822 macrophage/monocyte cells from three HFpEF and three control biological samples.",
  "",
  "Program definition:",
  "- Ccr2-positive, Ccr2-negative, and cross-subset Stage 2 drug-opposed programs were reconstructed at Top50, Top100, Top150, and Top200 sizes.",
  "- The primary size was Top150, but all sizes were retained in auditable tables.",
  "- Stage 3-supported subsets required concordant disease direction in macrophage pseudobulk edgeR and limma-voom results.",
  "",
  "Regulon network:",
  paste0("- Actual network mode: ", network_mode, "."),
  paste0("- Actual network source: ", paste(unique(network_dt$network_source), collapse = "; "), "."),
  "- A broad mouse TF-target prior was preferred; no fixed candidate-TF list was used.",
  "- If no prior resource was available, a clearly labelled sample-state pseudobulk Spearman fallback was used.",
  "- Regulons required at least 10 detected targets and were capped at 250 targets.",
  "",
  "Activity inference and statistics:",
  "- The primary signed activity score was a weighted mean of gene-wise standardized target expression, incorporating activation or repression direction.",
  paste0("- AUCell sensitivity status: ", aucell_status, "."),
  "- HFpEF versus control comparisons used six biological-sample pseudobulks.",
  "- limma moderated statistics, Wilcoxon P values, Hedges g, and FDR values were retained.",
  "- Cell-level UMAP activity maps were descriptive and were not used as independent replicates.",
  "",
  "Robustness and comparison:",
  "- Nine leave-one-Control/leave-one-HFpEF-pair analyses quantified sign stability and top-rank frequency.",
  "- Regulon activity was compared with TF-expression-only ranking.",
  "- When available, weighted regulon effects were compared with signed AUCell effects.",
  "",
  "Claim boundary:",
  "- Stage 4 prioritizes regulatory hypotheses; it does not demonstrate TF binding or causal perturbation.",
  "- GSE236585 contains no dapagliflozin exposure and therefore does not directly test drug action.",
  "- Nfkb1 was not forced into the analysis; its final position, if present, is data-derived."
)
writeLines(
  methods_text,
  file.path(DIRS$methods, "stage4_methods_and_claim_boundaries.txt"),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(DIRS$methods, "sessionInfo.txt")
)

############################################################
## 18. Completion checks, status, and CHECK package
############################################################

nfkb1_row <- priority_dt[gene_key(tf_symbol) == "NFKB1"]
nfkb1_rank <- if (nrow(nfkb1_row) > 0L) {
  nfkb1_row$priority_rank[1L]
} else {
  NA_integer_
}

scientific_checks <- data.table::data.table(
  check = c(
    "Stage 2 ready",
    "Stage 3 ready",
    "Stage 3 checks failed",
    "Biological samples",
    "Control samples",
    "HFpEF samples",
    "Macrophage cells",
    "Valid TF regulons",
    "Weighted sample activity TFs",
    "Primary program sets",
    "Stage 3-supported program genes",
    "Leave-pair-out runs",
    "Priority table TFs",
    "Nfkb1 forced",
    "Output workbook"
  ),
  observed = c(
    as.integer(
      stage2_status$overall_status[1L] ==
        "COMPLETED_STAGE2_READY_FOR_REVIEW"
    ),
    as.integer(
      stage3_status$overall_status[1L] ==
        "COMPLETED_STAGE3_READY_FOR_REVIEW"
    ),
    sum(stage3_checks$status != "PASS"),
    data.table::uniqueN(sample_meta$sample_accession),
    sum(sample_meta$condition == "Control"),
    sum(sample_meta$condition == "HFpEF"),
    ncol(macrophage),
    length(valid_tfs),
    nrow(weighted_sample_activity),
    data.table::uniqueN(
      primary_manifest$program_name
    ),
    nrow(supported_programs),
    data.table::uniqueN(
      paste(
        loo_results$removed_control,
        loo_results$removed_hfpef,
        sep = "__"
      )
    ),
    nrow(priority_dt),
    0L,
    as.integer(file.exists(workbook_path))
  ),
  expected = c(
    1L,
    1L,
    0L,
    6L,
    3L,
    3L,
    1822L,
    10L,
    10L,
    6L,
    10L,
    9L,
    10L,
    0L,
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
    "at_least",
    "at_least",
    "equal",
    "at_least",
    "equal",
    "at_least",
    "equal",
    "equal"
  )
)

scientific_checks[, status := data.table::fcase(
  comparison == "equal" & observed == expected,
  "PASS",
  comparison == "at_least" & observed >= expected,
  "PASS",
  default = "FAIL"
)]

data.table::fwrite(
  scientific_checks,
  file.path(DIRS$tables, "20_stage4_scientific_completion_checks.csv")
)

warnings_dt <- if (length(warning_records) > 0L) {
  data.table::rbindlist(
    warning_records,
    use.names = TRUE,
    fill = TRUE
  )
} else {
  data.table::data.table(
    timestamp = character(),
    category = character(),
    item = character(),
    message = character()
  )
}

data.table::fwrite(
  warnings_dt,
  file.path(DIRS$tables, "21_stage4_warnings_and_nonfatal_issues.csv")
)

script_copy_status <- "NOT_DETECTED"
if (
  length(SCRIPT_FILE) == 1L &&
  !is.na(SCRIPT_FILE) &&
  file.exists(SCRIPT_FILE)
) {
  methods_script <- file.path(
    DIRS$methods,
    basename(EXPECTED_SCRIPT_FILE)
  )
  check_script <- file.path(
    DIRS$check,
    basename(EXPECTED_SCRIPT_FILE)
  )
  copied_methods <- file.copy(
    SCRIPT_FILE,
    methods_script,
    overwrite = TRUE
  )
  copied_check <- file.copy(
    SCRIPT_FILE,
    check_script,
    overwrite = TRUE
  )
  script_copy_status <- if (
    isTRUE(copied_methods) &&
    isTRUE(copied_check)
  ) {
    "COPIED"
  } else {
    "COPY_FAILED"
  }
}

END_TIME <- Sys.time()
base_ready <- all(scientific_checks$status == "PASS")
overall_status <- if (!base_ready) {
  "COMPLETED_STAGE4_REVIEW_REQUIRED"
} else if (network_mode == "DATA_DRIVEN_FALLBACK") {
  "COMPLETED_STAGE4_READY_WITH_METHOD_CAUTION"
} else {
  "COMPLETED_STAGE4_READY_FOR_REVIEW"
}

run_status <- data.table::data.table(
  stage = STAGE_NAME,
  start_time = format(START_TIME, "%Y-%m-%d %H:%M:%S"),
  end_time = format(END_TIME, "%Y-%m-%d %H:%M:%S"),
  elapsed_minutes = round(
    as.numeric(
      difftime(END_TIME, START_TIME, units = "mins")
    ),
    2
  ),
  biological_samples = nrow(sample_meta),
  macrophage_cells = ncol(macrophage),
  valid_regulons = length(valid_tfs),
  network_mode = network_mode,
  network_source = paste(
    unique(network_dt$network_source),
    collapse = ";"
  ),
  aucell_status = aucell_status,
  leave_pair_out_runs = 9L,
  prioritized_TFs = nrow(priority_dt),
  top_TF = priority_dt$tf_symbol[1L],
  top_TF_priority_score = priority_dt$priority_score[1L],
  Nfkb1_rank = nfkb1_rank,
  Nfkb1_forced = FALSE,
  warnings = nrow(warnings_dt),
  script_copy_status = script_copy_status,
  scientific_checks_failed = sum(
    scientific_checks$status != "PASS"
  ),
  overall_status = overall_status
)

data.table::fwrite(
  run_status,
  file.path(DIRS$tables, "22_stage4_run_status.csv")
)

readme <- c(
  "HFpEF Reanalysis Project - Stage 4 FIXED v1",
  "GSE236585 macrophage TF-regulon consensus prioritization",
  "",
  paste0("Overall status: ", overall_status),
  paste0("Biological samples: ", nrow(sample_meta)),
  paste0("Macrophage cells: ", ncol(macrophage)),
  paste0("Valid regulons: ", length(valid_tfs)),
  paste0("Network mode: ", network_mode),
  paste0("Network source: ", paste(unique(network_dt$network_source), collapse = ";")),
  paste0("AUCell status: ", aucell_status),
  paste0("Top TF: ", priority_dt$tf_symbol[1L]),
  paste0("Nfkb1 rank: ", ifelse(is.na(nfkb1_rank), "not ranked", nfkb1_rank)),
  paste0("Script snapshot: ", script_copy_status),
  "",
  "Primary boundaries:",
  "- No TF was forced into the ranking.",
  "- Biological samples are the inferential units.",
  "- Cell-level maps are descriptive.",
  "- Regulatory activity and target links are computational prioritization evidence.",
  "- GSE236585 does not contain dapagliflozin treatment.",
  "",
  "Upload the Stage 4 CHECK package before starting Stage 5."
)
writeLines(
  readme,
  file.path(OUT_DIR, "README_stage4.txt"),
  useBytes = TRUE
)

## Compact review files.
review_files <- c(
  LOG_FILE,
  file.path(DIRS$tables, "00_stage4_optional_package_status.csv"),
  file.path(DIRS$tables, "01_stage4_input_validation.csv"),
  file.path(DIRS$tables, "03_stage4_program_size_summary.csv"),
  file.path(DIRS$tables, "04_stage4_macrophage_pseudobulk_audit.csv"),
  file.path(DIRS$tables, "07_stage4_regulon_size_summary.csv"),
  file.path(DIRS$tables, "08_stage4_regulon_network_audit.csv"),
  file.path(DIRS$tables, "09_stage4_weighted_regulon_activity_HFpEF_vs_Control.csv"),
  file.path(DIRS$tables, "10_stage4_AUCell_regulon_activity_HFpEF_vs_Control.csv"),
  file.path(DIRS$tables, "11_stage4_TF_expression_HFpEF_vs_Control.csv"),
  file.path(DIRS$tables, "12_stage4_candidate_TF_priority_score.csv"),
  file.path(DIRS$tables, "15_stage4_leave_one_pair_out_TF_robustness_summary.csv"),
  file.path(DIRS$tables, "18_stage4_TF_method_comparison_summary.csv"),
  file.path(DIRS$tables, "19_GSE236585_stage4_TF_regulon_key_results.xlsx"),
  file.path(DIRS$tables, "20_stage4_scientific_completion_checks.csv"),
  file.path(DIRS$tables, "21_stage4_warnings_and_nonfatal_issues.csv"),
  file.path(DIRS$tables, "22_stage4_run_status.csv"),
  file.path(DIRS$methods, "stage4_parameters_and_rationale.csv"),
  file.path(DIRS$methods, "stage4_priority_score_definition.txt"),
  file.path(DIRS$methods, "stage4_methods_and_claim_boundaries.txt"),
  file.path(DIRS$methods, "sessionInfo.txt"),
  file.path(OUT_DIR, "README_stage4.txt"),
  list.files(DIRS$figures, pattern = "\\.png$", full.names = TRUE)
)

review_files <- unique(review_files[file.exists(review_files)])
for (f in review_files) {
  target <- file.path(DIRS$check, basename(f))
  if (
    normalizePath(f, winslash = "/", mustWork = FALSE) !=
      normalizePath(target, winslash = "/", mustWork = FALSE)
  ) {
    file.copy(f, target, overwrite = TRUE)
  }
}

check_files <- list.files(DIRS$check, full.names = TRUE)
check_manifest <- data.table::data.table(
  filename = basename(check_files),
  size_bytes = as.numeric(file.info(check_files)$size)
)
check_manifest[, sha256 := vapply(
  check_files,
  function(f) {
    digest::digest(
      file = f,
      algo = "sha256",
      serialize = FALSE
    )
  },
  character(1)
)]

data.table::fwrite(
  check_manifest,
  file.path(DIRS$check, "CHECK_package_file_manifest.csv")
)

if (file.exists(CHECK_ZIP)) {
  unlink(CHECK_ZIP, force = TRUE)
}
zip::zipr(
  zipfile = CHECK_ZIP,
  files = list.files(DIRS$check, full.names = TRUE),
  root = DIRS$check
)

log_msg("Stage 4 analysis finished.")
log_msg("Overall status: ", overall_status)
log_msg("Network mode: ", network_mode)
log_msg("Valid regulons: ", length(valid_tfs))
log_msg("Top TF: ", priority_dt$tf_symbol[1L])
log_msg("Nfkb1 rank: ", ifelse(is.na(nfkb1_rank), "not ranked", nfkb1_rank))
log_msg("CHECK package: ", CHECK_ZIP)

cat("\n============================================================\n")
cat("HFpEF Stage 4 TF-regulon analysis completed\n")
cat("Status: ", overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat("Top TF: ", priority_dt$tf_symbol[1L], "\n", sep = "")
cat("Nfkb1 was not forced.\n")
cat("Upload the CHECK package before Stage 5.\n")
cat("============================================================\n")
