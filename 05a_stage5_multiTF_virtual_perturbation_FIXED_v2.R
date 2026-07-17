############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 5 FIXED v2
## Comparative multi-candidate in-silico TF perturbation
##
## Current project:
##   <HFPEF_PROJECT_DIR>
##
## Required completed input:
##   04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1
##
## Critical v2 correction:
##   Stage 5 v1 contained a data.table scoping collision because the
##   Stage 4 network had both source_symbol and tf_symbol columns.
##   Inside the perturbation function, source_symbol == tf_symbol
##   compared two columns and selected the whole network for every TF.
##   v2 removes the alias column, uses a non-colliding tf_requested
##   variable, and verifies every TF-specific target count against the
##   locked Stage 4 regulon-size table before any perturbation.
##
## Prespecified biological candidates:
##   Bhlhe40, Runx1, Spi1, Rel, Nfkb1, Rela
##
## Primary scientific question:
##   Which candidate TF activity normalization most consistently moves
##   HFpEF macrophage transcription toward the Control state across:
##     1) Stage 2 drug-opposed programs;
##     2) Stage 3-supported program subsets;
##     3) inflammatory and macrophage-state programs;
##     4) biological samples;
##     5) macrophage states;
##     6) two transparent perturbation formulations?
##
## Important design choice:
##   The primary intervention is DISEASE-NORMALIZING activity adjustment,
##   not uniform knockout. This is necessary because some candidate TF
##   activities are higher in HFpEF, whereas others are lower. Uniform
##   inhibition is retained as a secondary directionality test.
##
## Primary perturbation formulations:
##   A) weighted minimum-norm target adjustment;
##   B) equal-signed target adjustment as an alternative formulation.
##
## Sensitivity dimensions:
##   - 25%, 50%, 75%, and 100% activity adjustment;
##   - disease normalization versus activity attenuation;
##   - matched low-priority TF controls;
##   - sample-level and macrophage-state-level analyses;
##   - rank aggregation rather than a single opaque weighted score.
##
## Interpretation boundary:
##   - This is computational perturbation prioritization.
##   - It is not experimental TF knockdown or knockout.
##   - It does not prove direct dapagliflozin action on any TF.
##   - Biological samples, not cells, are the inferential units.
##   - Predicted expression shifts are network-constrained estimates.
##
## Output:
##   <HFPEF_PROJECT_DIR>/
##   05_stage5_multiTF_virtual_perturbation_FIXED_v2
##
## CHECK:
##   <HFPEF_PROJECT_DIR>/
##   05_stage5_multiTF_virtual_perturbation_FIXED_v2_CHECK.zip
##
## Recommended run:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage5_MultiTF_Virtual_Perturbation_FIXED_v2.R",
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
## 0. Paths and analysis settings
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

STAGE4_DIR <- file.path(
  PROJECT_DIR,
  "04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1"
)

STAGE4_STATUS_FILE <- file.path(
  STAGE4_DIR,
  "01_tables",
  "22_stage4_run_status.csv"
)
STAGE4_CHECKS_FILE <- file.path(
  STAGE4_DIR,
  "01_tables",
  "20_stage4_scientific_completion_checks.csv"
)
STAGE4_PRIORITY_FILE <- file.path(
  STAGE4_DIR,
  "01_tables",
  "12_stage4_candidate_TF_priority_score.csv"
)
STAGE4_NETWORK_FILE <- file.path(
  STAGE4_DIR,
  "01_tables",
  "05_stage4_full_TF_target_links.csv"
)
STAGE4_REGULON_SIZE_FILE <- file.path(
  STAGE4_DIR,
  "01_tables",
  "07_stage4_regulon_size_summary.csv"
)
STAGE4_PROGRAM_FILE <- file.path(
  STAGE4_DIR,
  "01_tables",
  "02_stage4_program_gene_manifest.csv"
)
STAGE4_PSEUDOBULK_RDS <- file.path(
  STAGE4_DIR,
  "02_objects",
  "GSE236585_stage4_macrophage_pseudobulk_matrices.rds"
)
STAGE4_ACTIVITY_RDS <- file.path(
  STAGE4_DIR,
  "02_objects",
  "stage4_weighted_regulon_activity_sample_matrix.rds"
)
STAGE4_MACROPHAGE_RDS <- file.path(
  STAGE4_DIR,
  "02_objects",
  "GSE236585_stage4_macrophage_regulon_scored.rds"
)

STAGE3_DIR <- file.path(
  PROJECT_DIR,
  "03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH"
)
STAGE3_SAMPLE_META_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "01_locked_GSE236585_sample_metadata.csv"
)

STAGE_NAME <- "05_stage5_multiTF_virtual_perturbation_FIXED_v2"
OUT_DIR <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP <- file.path(
  PROJECT_DIR,
  paste0(STAGE_NAME, "_CHECK.zip")
)

EXPECTED_SCRIPT_FILE <- file.path(PROJECT_DIR, "05a_stage5_multiTF_virtual_perturbation_FIXED_v2.R")

REPLACE_EXISTING_STAGE5 <- TRUE

CANDIDATE_TFS_REQUESTED <- c(
  "Bhlhe40",
  "Runx1",
  "Spi1",
  "Rel",
  "Nfkb1",
  "Rela"
)

N_MATCHED_CONTROL_TFS <- 3L
CONTROL_PRIORITY_QUANTILE <- 0.60
CONTROL_MAX_PROGRAM_OVERLAP <- 1L
CONTROL_MAX_SUPPORTED_OVERLAP <- 0L

PERTURBATION_METHODS <- c(
  "weighted_minimum_norm",
  "equal_signed_targets"
)
PERTURBATION_MODES <- c(
  "disease_normalization",
  "activity_attenuation"
)
PERTURBATION_STRENGTHS <- c(0.25, 0.50, 0.75, 1.00)

PRIMARY_METHOD <- "weighted_minimum_norm"
PRIMARY_MODE <- "disease_normalization"
PRIMARY_STRENGTH <- 1.00
PRIMARY_SIGNATURE_SIZE <- 150L

MAX_ABS_GENE_SHIFT_SD <- 2.50
MIN_ABS_OBSERVED_PROGRAM_GAP <- 0.10
MIN_TARGETS_PER_TF <- 10L
MIN_STATE_PROFILES_PER_CONDITION <- 2L

TOP_GENES_PER_TF_REPORT <- 100L
TOP_LIGANDS_PER_TF_REPORT <- 30L
TOP_TFS_FOR_FIGURES <- 9L

############################################################
## 1. Preflight, output setup, and logging
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
  STAGE4_STATUS_FILE,
  STAGE4_CHECKS_FILE,
  STAGE4_PRIORITY_FILE,
  STAGE4_NETWORK_FILE,
  STAGE4_REGULON_SIZE_FILE,
  STAGE4_PROGRAM_FILE,
  STAGE4_PSEUDOBULK_RDS,
  STAGE4_ACTIVITY_RDS,
  STAGE4_MACROPHAGE_RDS,
  STAGE3_SAMPLE_META_FILE
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Required Stage 3/4 input path(s) are missing:\n",
    paste(missing_inputs, collapse = "\n")
  )
}

stage4_status <- data.table::fread(
  STAGE4_STATUS_FILE,
  encoding = "UTF-8"
)

allowed_stage4_status <- c(
  "COMPLETED_STAGE4_READY_FOR_REVIEW",
  "COMPLETED_STAGE4_READY_WITH_METHOD_CAUTION"
)

if (
  !"overall_status" %in% names(stage4_status) ||
  !stage4_status$overall_status[1L] %in%
    allowed_stage4_status
) {
  stop(
    "Stage 4 is not in an allowed completed state. Observed: ",
    ifelse(
      "overall_status" %in% names(stage4_status),
      stage4_status$overall_status[1L],
      "missing overall_status"
    )
  )
}

stage4_checks <- data.table::fread(
  STAGE4_CHECKS_FILE,
  encoding = "UTF-8"
)

if (
  !all(c("check", "status") %in% names(stage4_checks)) ||
  any(stage4_checks$status != "PASS")
) {
  stop(
    "At least one Stage 4 scientific completion check is not PASS."
  )
}

replacement_audit <- data.table::data.table(
  path = c(OUT_DIR, CHECK_ZIP),
  path_type = c(
    "stage5_output_directory",
    "stage5_check_zip"
  ),
  existed_before = FALSE,
  deletion_attempted = FALSE,
  deletion_succeeded = FALSE
)

if (REPLACE_EXISTING_STAGE5) {
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
        stop("Failed to remove previous Stage 5 path: ", target)
      }
    } else {
      replacement_audit$deletion_succeeded[i] <- TRUE
    }
  }
} else if (dir.exists(OUT_DIR) || file.exists(CHECK_ZIP)) {
  stop(
    "Existing Stage 5 output detected while replacement is disabled."
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
LOG_FILE <- file.path(
  DIRS$logs,
  "stage5_multiTF_virtual_perturbation.log"
)
WARN_FILE <- file.path(
  DIRS$logs,
  "stage5_warnings.log"
)

data.table::fwrite(
  replacement_audit,
  file.path(DIRS$logs, "stage5_replacement_audit.csv")
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

log_msg("Stage 5 multi-TF virtual perturbation started.")
log_msg("Stage 4 status: ", stage4_status$overall_status[1L])
log_msg("Output: ", OUT_DIR)
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
      "Required package(s) unavailable:",
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
    "writexl",
    "scales",
    "zip",
    "digest"
  ),
  required = TRUE
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
})

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

read_table_auto <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2L)
  close(con)
  on.exit(NULL, add = FALSE)

  is_gzip <- length(magic) == 2L &&
    identical(
      as.integer(magic),
      c(31L, 139L)
    )

  if (is_gzip) {
    data.table::as.data.table(
      utils::read.csv(
        gzfile(path, open = "rt"),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )
  } else {
    data.table::fread(path, encoding = "UTF-8")
  }
}

rescale01 <- function(x, neutral_if_constant = 0.5) {
  x <- as.numeric(x)
  out <- rep(neutral_if_constant, length(x))
  finite <- is.finite(x)

  if (!any(finite)) return(out)

  rng <- range(x[finite], na.rm = TRUE)

  if (
    !is.finite(rng[1L]) ||
    !is.finite(rng[2L]) ||
    diff(rng) == 0
  ) {
    out[finite] <- neutral_if_constant
    return(out)
  }

  out[finite] <- (x[finite] - rng[1L]) / diff(rng)
  out
}

safe_median <- function(x, default = NA_real_) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(default)
  stats::median(x)
}

safe_mean <- function(x, default = NA_real_) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(default)
  mean(x)
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

safe_wilcox_p <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 1L || length(y) < 1L) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(
      x,
      y,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

write_csv_safe <- function(x, path, compress = FALSE) {
  if (is.null(x) || ncol(x) == 0L) {
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
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
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

save_plot_bundle <- function(
  plot_object,
  stem,
  width,
  height
) {
  paths <- c(
    png = file.path(DIRS$figures, paste0(stem, ".png")),
    pdf = file.path(DIRS$figures, paste0(stem, ".pdf")),
    tiff = file.path(DIRS$figures, paste0(stem, ".tiff"))
  )

  ggplot2::ggsave(
    paths["png"],
    plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    paths["pdf"],
    plot_object,
    width = width,
    height = height,
    bg = "white"
  )

  ggplot2::ggsave(
    paths["tiff"],
    plot_object,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )

  invisible(paths)
}

save_heatmap_bundle <- function(
  mat,
  stem,
  annotation_col = NULL,
  width = 10,
  height = 8,
  main = NULL
) {
  paths <- c(
    png = file.path(DIRS$figures, paste0(stem, ".png")),
    pdf = file.path(DIRS$figures, paste0(stem, ".pdf")),
    tiff = file.path(DIRS$figures, paste0(stem, ".tiff"))
  )

  grDevices::png(
    paths["png"],
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
    main = main
  )
  grDevices::dev.off()

  grDevices::pdf(
    paths["pdf"],
    width = width,
    height = height
  )
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    border_color = NA,
    main = main
  )
  grDevices::dev.off()

  grDevices::tiff(
    paths["tiff"],
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
    main = main
  )
  grDevices::dev.off()

  invisible(paths)
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

make_feature_map <- function(features) {
  dt <- data.table::data.table(
    feature = as.character(features),
    feature_key = gene_key(features)
  )

  data.table::setorder(dt, feature_key, feature)
  dt[, .SD[1L], by = feature_key]
}

map_keys_to_features <- function(keys, feature_map) {
  keys <- unique(gene_key(keys))
  keys <- keys[nzchar(keys)]

  out <- merge(
    data.table::data.table(feature_key = keys),
    feature_map,
    by = "feature_key",
    all.x = TRUE,
    sort = FALSE
  )

  unique(out[!is.na(feature), feature])
}

parse_state_profile_names <- function(x) {
  parts <- strsplit(
    as.character(x),
    "__",
    fixed = TRUE
  )

  sample_accession <- vapply(
    parts,
    function(z) z[1L],
    character(1)
  )

  macrophage_state <- vapply(
    parts,
    function(z) {
      if (length(z) <= 1L) {
        return("UNKNOWN_STATE")
      }
      paste(z[-1L], collapse = "__")
    },
    character(1)
  )

  data.table::data.table(
    profile = as.character(x),
    sample_accession = sample_accession,
    macrophage_state = macrophage_state
  )
}

weighted_regulon_activity <- function(
  expr_matrix,
  network_dt,
  tf_subset = NULL
) {
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

  net <- data.table::copy(network_dt)

  if (!is.null(tf_subset)) {
    net <- net[source_symbol %in% tf_subset]
  }

  target_features <- intersect(
    unique(net$target_feature),
    rownames(expr_matrix)
  )

  if (length(target_features) < 10L) {
    stop(
      "Too few network targets are present in the expression matrix."
    )
  }

  z <- scale_rows(
    expr_matrix[target_features, , drop = FALSE]
  )

  tfs <- sort(unique(net$source_symbol))

  out <- matrix(
    NA_real_,
    nrow = length(tfs),
    ncol = ncol(z),
    dimnames = list(tfs, colnames(z))
  )

  for (tf in tfs) {
    net_tf <- net[
      source_symbol == tf &
        target_feature %in% rownames(z)
    ]

    net_tf <- net_tf[
      match(unique(target_feature), target_feature)
    ]

    if (nrow(net_tf) < MIN_TARGETS_PER_TF) next

    signed_weights <- as.numeric(net_tf$mor) *
      as.numeric(net_tf$weight)
    signed_weights[!is.finite(signed_weights)] <- 0

    denom <- sum(abs(signed_weights))
    if (!is.finite(denom) || denom <= 0) next

    target_z <- z[
      net_tf$target_feature,
      ,
      drop = FALSE
    ]

    out[tf, ] <- colSums(
      sweep(
        target_z,
        1L,
        signed_weights,
        "*"
      )
    ) / denom
  }

  out[
    rowSums(is.finite(out)) > 0L,
    ,
    drop = FALSE
  ]
}

make_program_definitions <- function(
  stage4_manifest,
  feature_map
) {
  manifest <- data.table::copy(stage4_manifest)

  required <- c(
    "program_name",
    "subset_name",
    "direction",
    "signature_size",
    "symbol_key"
  )

  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0L) {
    stop(
      "Stage 4 program manifest is missing column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  manifest[, support_class := ifelse(
    grepl("Stage3Supported$", program_name),
    "Stage3_supported",
    "Full_Stage2"
  )]

  definitions <- list()
  summary_records <- list()

  grouping <- unique(
    manifest[
      ,
      .(
        subset_name,
        signature_size,
        support_class
      )
    ]
  )

  for (i in seq_len(nrow(grouping))) {
    subset_i <- grouping$subset_name[i]
    size_i <- grouping$signature_size[i]
    support_i <- grouping$support_class[i]

    part <- manifest[
      subset_name == subset_i &
        signature_size == size_i &
        support_class == support_i
    ]

    up_keys <- unique(
      part[
        direction == "Disease_up_Drug_down",
        symbol_key
      ]
    )

    down_keys <- unique(
      part[
        direction == "Disease_down_Drug_up",
        symbol_key
      ]
    )

    up_features <- map_keys_to_features(
      up_keys,
      feature_map
    )
    down_features <- map_keys_to_features(
      down_keys,
      feature_map
    )

    if (
      length(up_features) + length(down_features) <
        5L
    ) {
      next
    }

    program_id <- paste(
      "DrugOpposedNet",
      subset_i,
      support_i,
      paste0("Top", size_i),
      sep = "__"
    )

    definitions[[program_id]] <- list(
      program_id = program_id,
      program_category = "Stage2_drug_opposed_net",
      subset_name = subset_i,
      signature_size = as.integer(size_i),
      support_class = support_i,
      up_features = up_features,
      down_features = down_features,
      primary_program = (
        as.integer(size_i) == PRIMARY_SIGNATURE_SIZE
      )
    )

    summary_records[[length(summary_records) + 1L]] <-
      data.table::data.table(
        program_id = program_id,
        program_category = "Stage2_drug_opposed_net",
        subset_name = subset_i,
        signature_size = as.integer(size_i),
        support_class = support_i,
        detected_up_genes = length(up_features),
        detected_down_genes = length(down_features),
        primary_program = (
          as.integer(size_i) == PRIMARY_SIGNATURE_SIZE
        )
      )
  }

  functional_sets <- list(
    Inflammatory_Il1b = c(
      "Il1b", "Tnf", "Nfkbia", "Ccl2",
      "Ccl3", "Ccl4", "Cxcl2", "S100a8",
      "S100a9", "Ptgs2"
    ),
    NFkB_TNF_response = c(
      "Tnf", "Nfkbia", "Nfkbiz", "Rel",
      "Rela", "Nfkb1", "Tnfaip3", "Icam1",
      "Ccl2", "Ccl3", "Il1b"
    ),
    Inflammasome_pyroptosis = c(
      "Nlrp3", "Pycard", "Casp1", "Gsdmd",
      "Il1b", "Il18", "Txnip", "P2rx7"
    ),
    Interferon_response = c(
      "Isg15", "Ifit1", "Ifit2", "Ifit3",
      "Irf7", "Rsad2", "Oas1a", "Stat1",
      "Cxcl9", "Cxcl10"
    ),
    Antigen_presentation = c(
      "H2-Ab1", "H2-Aa", "Cd74",
      "Ciita", "H2-Eb1", "Tap1", "B2m"
    ),
    Spp1_Trem2_remodeling = c(
      "Spp1", "Trem2", "Gpnmb", "Fabp5",
      "Lpl", "Apoe", "Ctsb", "Ctsd",
      "Lgals3"
    ),
    Resident_Timd4_Lyve1 = c(
      "Timd4", "Lyve1", "Folr2", "Mrc1",
      "Cd163", "Vsig4", "C1qa", "C1qb",
      "C1qc"
    ),
    Ccr2_monocyte_like = c(
      "Ccr2", "Ly6c2", "Plac8",
      "Chil3", "Ctss", "Lgals3"
    ),
    Myeloid_identity = c(
      "Spi1", "Csf1r", "Lyz2", "Aif1",
      "Tyrobp", "Fcerg", "Ctss", "Laptm5",
      "C1qa", "C1qb", "C1qc"
    ),
    Lipid_cholesterol = c(
      "Apoe", "Lpl", "Abca1", "Abcg1",
      "Pparg", "Nr1h3", "Soat1", "Lipa",
      "Fabp5"
    ),
    Oxidative_stress = c(
      "Nfe2l2", "Hmox1", "Nqo1", "Gclc",
      "Gclm", "Sod2", "Txnrd1", "Prdx1"
    ),
    Cycling = c(
      "Mki67", "Top2a", "Stmn1",
      "Tubb5", "Hmgb2"
    )
  )

  for (set_name in names(functional_sets)) {
    features <- map_keys_to_features(
      functional_sets[[set_name]],
      feature_map
    )

    if (length(features) < 3L) next

    program_id <- paste0(
      "Functional__",
      set_name
    )

    definitions[[program_id]] <- list(
      program_id = program_id,
      program_category = "Functional_state",
      subset_name = set_name,
      signature_size = NA_integer_,
      support_class = "Curated",
      up_features = features,
      down_features = character(),
      primary_program = set_name %in% c(
        "Inflammatory_Il1b",
        "NFkB_TNF_response",
        "Inflammasome_pyroptosis",
        "Interferon_response",
        "Spp1_Trem2_remodeling",
        "Myeloid_identity"
      )
    )

    summary_records[[length(summary_records) + 1L]] <-
      data.table::data.table(
        program_id = program_id,
        program_category = "Functional_state",
        subset_name = set_name,
        signature_size = NA_integer_,
        support_class = "Curated",
        detected_up_genes = length(features),
        detected_down_genes = 0L,
        primary_program = set_name %in% c(
          "Inflammatory_Il1b",
          "NFkB_TNF_response",
          "Inflammasome_pyroptosis",
          "Interferon_response",
          "Spp1_Trem2_remodeling",
          "Myeloid_identity"
        )
      )
  }

  list(
    definitions = definitions,
    summary = data.table::rbindlist(
      summary_records,
      use.names = TRUE,
      fill = TRUE
    )
  )
}

score_programs <- function(z_matrix, definitions) {
  out <- matrix(
    NA_real_,
    nrow = length(definitions),
    ncol = ncol(z_matrix),
    dimnames = list(
      names(definitions),
      colnames(z_matrix)
    )
  )

  for (program_id in names(definitions)) {
    def <- definitions[[program_id]]

    up <- intersect(
      def$up_features,
      rownames(z_matrix)
    )
    down <- intersect(
      def$down_features,
      rownames(z_matrix)
    )

    up_score <- if (length(up) > 0L) {
      colMeans(
        z_matrix[up, , drop = FALSE],
        na.rm = TRUE
      )
    } else {
      rep(0, ncol(z_matrix))
    }

    down_score <- if (length(down) > 0L) {
      colMeans(
        z_matrix[down, , drop = FALSE],
        na.rm = TRUE
      )
    } else {
      rep(0, ncol(z_matrix))
    }

    out[program_id, ] <- up_score - down_score
  }

  out
}

simulate_tf_activity_adjustment <- function(
  observed_z,
  activity_matrix,
  tf_requested,
  network_dt,
  sample_meta,
  perturbation_mode,
  perturbation_strength,
  perturbation_method
) {
  if (!tf_requested %in% rownames(activity_matrix)) {
    stop("TF is missing from the activity matrix: ", tf_requested)
  }

  sample_ids <- colnames(observed_z)

  meta <- data.table::copy(sample_meta)
  meta <- meta[match(sample_ids, sample_accession)]

  if (
    any(is.na(meta$sample_accession)) ||
    any(meta$sample_accession != sample_ids)
  ) {
    stop("Sample metadata could not be aligned in perturbation.")
  }

  observed_activity <- as.numeric(
    activity_matrix[tf_requested, sample_ids]
  )
  names(observed_activity) <- sample_ids

  control_reference <- mean(
    observed_activity[meta$condition == "Control"],
    na.rm = TRUE
  )

  desired_activity <- observed_activity
  hfpef_idx <- which(meta$condition == "HFpEF")

  if (perturbation_mode == "disease_normalization") {
    desired_activity[hfpef_idx] <-
      observed_activity[hfpef_idx] +
      perturbation_strength *
        (
          control_reference -
            observed_activity[hfpef_idx]
        )
  } else if (
    perturbation_mode == "activity_attenuation"
  ) {
    desired_activity[hfpef_idx] <-
      observed_activity[hfpef_idx] *
      (1 - perturbation_strength)
  } else {
    stop(
      "Unknown perturbation_mode: ",
      perturbation_mode
    )
  }

  delta_activity <- desired_activity - observed_activity

  net_tf <- data.table::copy(
    network_dt[
      source_symbol == tf_requested &
        target_feature %in% rownames(observed_z)
    ]
  )

  data.table::setorder(
    net_tf,
    target_feature,
    -weight
  )
  net_tf <- net_tf[, .SD[1L], by = target_feature]

  if (nrow(net_tf) < MIN_TARGETS_PER_TF) {
    stop(
      "Too few detected targets for ",
      tf_requested,
      ": ",
      nrow(net_tf)
    )
  }

  signed_weights <- as.numeric(net_tf$mor) *
    as.numeric(net_tf$weight)
  signed_weights[!is.finite(signed_weights)] <- 0

  denom_abs <- sum(abs(signed_weights))
  denom_sq <- sum(signed_weights^2)

  if (
    !is.finite(denom_abs) ||
    denom_abs <= 0 ||
    !is.finite(denom_sq) ||
    denom_sq <= 0
  ) {
    stop(
      "Invalid target weights for ",
      tf_requested
    )
  }

  if (
    perturbation_method ==
      "weighted_minimum_norm"
  ) {
    gene_coefficients <-
      signed_weights *
      denom_abs /
      denom_sq
  } else if (
    perturbation_method ==
      "equal_signed_targets"
  ) {
    gene_coefficients <- sign(signed_weights)
  } else {
    stop(
      "Unknown perturbation_method: ",
      perturbation_method
    )
  }

  target_delta <- outer(
    gene_coefficients,
    delta_activity,
    "*"
  )

  target_delta <- pmax(
    pmin(
      target_delta,
      MAX_ABS_GENE_SHIFT_SD
    ),
    -MAX_ABS_GENE_SHIFT_SD
  )

  rownames(target_delta) <- net_tf$target_feature
  colnames(target_delta) <- sample_ids

  perturbed_z <- observed_z
  perturbed_z[
    net_tf$target_feature,
    sample_ids
  ] <- perturbed_z[
    net_tf$target_feature,
    sample_ids,
    drop = FALSE
  ] + target_delta

  realized_delta_activity <- colSums(
    sweep(
      target_delta,
      1L,
      signed_weights,
      "*"
    )
  ) / denom_abs

  gene_effect <- data.table::data.table(
    tf_symbol = tf_requested,
    perturbation_mode = perturbation_mode,
    perturbation_strength = perturbation_strength,
    perturbation_method = perturbation_method,
    target_feature = net_tf$target_feature,
    target_key = gene_key(net_tf$target_feature),
    mor = as.numeric(net_tf$mor),
    edge_weight = as.numeric(net_tf$weight),
    gene_coefficient = gene_coefficients,
    mean_delta_z_HFpEF = rowMeans(
      target_delta[
        ,
        meta$condition == "HFpEF",
        drop = FALSE
      ],
      na.rm = TRUE
    ),
    max_abs_delta_z = apply(
      abs(target_delta),
      1L,
      max,
      na.rm = TRUE
    )
  )

  list(
    perturbed_z = perturbed_z,
    observed_activity = observed_activity,
    desired_activity = desired_activity,
    delta_activity = delta_activity,
    realized_delta_activity = realized_delta_activity,
    control_reference = control_reference,
    target_count = nrow(net_tf),
    global_rms_shift_HFpEF = sqrt(
      mean(
        target_delta[
          ,
          meta$condition == "HFpEF",
          drop = FALSE
        ]^2,
        na.rm = TRUE
      )
    ),
    gene_effect = gene_effect
  )
}

summarize_program_effects <- function(
  observed_scores,
  perturbed_scores,
  sample_meta,
  program_summary,
  tf_symbol,
  perturbation_mode,
  perturbation_strength,
  perturbation_method,
  target_count,
  global_rms_shift,
  delta_activity,
  realized_delta_activity
) {
  sample_ids <- colnames(observed_scores)

  meta <- data.table::copy(sample_meta)
  meta <- meta[match(sample_ids, sample_accession)]

  control_idx <- which(meta$condition == "Control")
  hfpef_idx <- which(meta$condition == "HFpEF")

  records <- lapply(
    rownames(observed_scores),
    function(program_id) {
      observed_control <- as.numeric(
        observed_scores[program_id, control_idx]
      )
      observed_hfpef <- as.numeric(
        observed_scores[program_id, hfpef_idx]
      )
      perturbed_hfpef <- as.numeric(
        perturbed_scores[program_id, hfpef_idx]
      )

      control_mean <- mean(
        observed_control,
        na.rm = TRUE
      )
      observed_hfpef_mean <- mean(
        observed_hfpef,
        na.rm = TRUE
      )
      perturbed_hfpef_mean <- mean(
        perturbed_hfpef,
        na.rm = TRUE
      )

      observed_gap <- observed_hfpef_mean -
        control_mean
      perturbed_gap <- perturbed_hfpef_mean -
        control_mean

      absolute_gap_reduction <-
        abs(observed_gap) -
        abs(perturbed_gap)

      gap_eligible <- is.finite(observed_gap) &&
        abs(observed_gap) >=
          MIN_ABS_OBSERVED_PROGRAM_GAP

      recovery_fraction <- if (
        gap_eligible
      ) {
        absolute_gap_reduction /
          abs(observed_gap)
      } else {
        NA_real_
      }

      sample_improved <- abs(
        perturbed_hfpef - control_mean
      ) < abs(
        observed_hfpef - control_mean
      )

      data.table::data.table(
        tf_symbol = tf_symbol,
        perturbation_mode = perturbation_mode,
        perturbation_strength =
          perturbation_strength,
        perturbation_method =
          perturbation_method,
        program_id = program_id,
        observed_control_mean = control_mean,
        observed_hfpef_mean =
          observed_hfpef_mean,
        perturbed_hfpef_mean =
          perturbed_hfpef_mean,
        observed_hfpef_minus_control =
          observed_gap,
        perturbed_hfpef_minus_control =
          perturbed_gap,
        hfpef_score_change =
          perturbed_hfpef_mean -
          observed_hfpef_mean,
        absolute_gap_reduction =
          absolute_gap_reduction,
        recovery_fraction =
          recovery_fraction,
        observed_gap_eligible =
          gap_eligible,
        sample_improvement_fraction =
          mean(sample_improved, na.rm = TRUE),
        target_count = target_count,
        global_rms_shift_HFpEF =
          global_rms_shift,
        mean_requested_delta_activity_HFpEF =
          mean(
            delta_activity[hfpef_idx],
            na.rm = TRUE
          ),
        mean_realized_delta_activity_HFpEF =
          mean(
            realized_delta_activity[hfpef_idx],
            na.rm = TRUE
          )
      )
    }
  )

  out <- data.table::rbindlist(records)

  merge(
    out,
    program_summary,
    by = "program_id",
    all.x = TRUE
  )
}

select_matched_control_tfs <- function(
  priority_dt,
  network_dt,
  candidate_tfs
) {
  regulon_sizes <- network_dt[
    ,
    .(
      regulon_size =
        data.table::uniqueN(target_feature)
    ),
    by = source_symbol
  ]

  pool <- merge(
    priority_dt,
    regulon_sizes,
    by.x = "tf_symbol",
    by.y = "source_symbol",
    all.x = TRUE
  )

  candidate_rows <- pool[
    tf_symbol %in% candidate_tfs
  ]

  cutoff_rank <- stats::quantile(
    pool$priority_rank,
    probs = CONTROL_PRIORITY_QUANTILE,
    na.rm = TRUE,
    names = FALSE
  )

  eligible <- pool[
    !tf_symbol %in% candidate_tfs &
      priority_rank >= cutoff_rank &
      max_overlap_count <=
        CONTROL_MAX_PROGRAM_OVERLAP &
      supported_overlap_count <=
        CONTROL_MAX_SUPPORTED_OVERLAP &
      is.finite(regulon_size) &
      regulon_size >= MIN_TARGETS_PER_TF
  ]

  excluded_family_keys <- c(
    "NFKB1",
    "NFKB2",
    "RELA",
    "RELB",
    "REL"
  )

  eligible <- eligible[
    !gene_key(tf_symbol) %in%
      excluded_family_keys
  ]

  if (nrow(eligible) < N_MATCHED_CONTROL_TFS) {
    eligible <- pool[
      !tf_symbol %in% candidate_tfs &
        is.finite(regulon_size) &
        regulon_size >= MIN_TARGETS_PER_TF
    ]
    eligible <- eligible[
      !gene_key(tf_symbol) %in%
        excluded_family_keys
    ]
  }

  if (nrow(eligible) < N_MATCHED_CONTROL_TFS) {
    stop(
      "Fewer than ",
      N_MATCHED_CONTROL_TFS,
      " eligible matched control TFs."
    )
  }

  candidate_medians <- c(
    log_regulon = safe_median(
      log1p(candidate_rows$regulon_size),
      0
    ),
    abs_activity = safe_median(
      abs(candidate_rows$weighted_hedges_g),
      0
    ),
    abs_expression = safe_median(
      abs(candidate_rows$expression_hedges_g),
      0
    )
  )

  scale_activity <- stats::mad(
    abs(pool$weighted_hedges_g),
    na.rm = TRUE
  )
  scale_expression <- stats::mad(
    abs(pool$expression_hedges_g),
    na.rm = TRUE
  )

  if (!is.finite(scale_activity) || scale_activity == 0) {
    scale_activity <- 1
  }
  if (!is.finite(scale_expression) || scale_expression == 0) {
    scale_expression <- 1
  }

  eligible[, matching_distance :=
    abs(
      log1p(regulon_size) -
        candidate_medians["log_regulon"]
    ) +
    abs(
      abs(weighted_hedges_g) -
        candidate_medians["abs_activity"]
    ) /
      scale_activity +
    abs(
      abs(expression_hedges_g) -
        candidate_medians["abs_expression"]
    ) /
      scale_expression
  ]

  data.table::setorder(
    eligible,
    matching_distance,
    -priority_rank
  )

  eligible[
    seq_len(
      min(N_MATCHED_CONTROL_TFS, .N)
    )
  ]
}

rank_metric <- function(
  x,
  higher_is_better = TRUE
) {
  x <- as.numeric(x)
  finite <- is.finite(x)

  ranks <- rep(
    max(1L, sum(finite)) + 1,
    length(x)
  )

  if (any(finite)) {
    ranks[finite] <- rank(
      if (higher_is_better) {
        -x[finite]
      } else {
        x[finite]
      },
      ties.method = "average"
    )
  }

  ranks
}

############################################################
## 4. Load and validate Stage 4 results
############################################################

priority_dt <- read_table_auto(
  STAGE4_PRIORITY_FILE
)

network_dt <- read_table_auto(
  STAGE4_NETWORK_FILE
)

stage4_regulon_sizes <- read_table_auto(
  STAGE4_REGULON_SIZE_FILE
)

program_manifest <- read_table_auto(
  STAGE4_PROGRAM_FILE
)

sample_meta <- data.table::fread(
  STAGE3_SAMPLE_META_FILE,
  encoding = "UTF-8"
)

sample_meta[, condition := factor(
  condition,
  levels = c("Control", "HFpEF")
)]

data.table::setorder(
  sample_meta,
  condition,
  sample_accession
)

if (
  data.table::uniqueN(
    sample_meta$sample_accession
  ) != 6L ||
  sum(sample_meta$condition == "Control") !=
    3L ||
  sum(sample_meta$condition == "HFpEF") !=
    3L
) {
  stop(
    "Stage 3 sample metadata is not the expected 3 + 3 design."
  )
}

required_priority_columns <- c(
  "tf_symbol",
  "priority_rank",
  "priority_score",
  "weighted_effect",
  "weighted_hedges_g",
  "expression_hedges_g",
  "max_overlap_count",
  "supported_overlap_count",
  "sign_stability"
)

missing_priority_columns <- setdiff(
  required_priority_columns,
  names(priority_dt)
)

if (length(missing_priority_columns) > 0L) {
  stop(
    "Stage 4 priority table is missing column(s): ",
    paste(missing_priority_columns, collapse = ", ")
  )
}

required_network_columns <- c(
  "source_symbol",
  "target_feature",
  "mor",
  "weight"
)

missing_network_columns <- setdiff(
  required_network_columns,
  names(network_dt)
)

if (length(missing_network_columns) > 0L) {
  stop(
    "Stage 4 network table is missing column(s): ",
    paste(missing_network_columns, collapse = ", ")
  )
}

network_dt[, mor := as.numeric(mor)]
network_dt[, weight := as.numeric(weight)]

## Remove the Stage 4 compatibility alias that caused the v1
## source_symbol == tf_symbol column-to-column comparison.
if ("tf_symbol" %in% names(network_dt)) {
  network_dt[, tf_symbol := NULL]
}

network_dt <- network_dt[
  !is.na(source_symbol) &
    nzchar(source_symbol) &
    !is.na(target_feature) &
    nzchar(target_feature) &
    is.finite(mor) &
    mor != 0 &
    is.finite(weight) &
    weight > 0
]

data.table::setorder(
  network_dt,
  source_symbol,
  target_feature,
  -weight
)

network_dt <- network_dt[
  ,
  .SD[1L],
  by = .(
    source_symbol,
    target_feature
  )
]

required_regulon_size_columns <- c(
  "source_symbol",
  "regulon_size"
)

missing_regulon_size_columns <- setdiff(
  required_regulon_size_columns,
  names(stage4_regulon_sizes)
)

if (length(missing_regulon_size_columns) > 0L) {
  stop(
    "Stage 4 regulon-size table is missing column(s): ",
    paste(missing_regulon_size_columns, collapse = ", ")
  )
}

stage5_regulon_sizes <- network_dt[
  ,
  .(
    stage5_regulon_size =
      data.table::uniqueN(target_feature)
  ),
  by = source_symbol
]

regulon_integrity_audit <- merge(
  stage4_regulon_sizes[
    ,
    .(
      source_symbol,
      stage4_regulon_size =
        as.integer(regulon_size)
    )
  ],
  stage5_regulon_sizes,
  by = "source_symbol",
  all = TRUE
)

regulon_integrity_audit[
  ,
  target_count_match := (
    !is.na(stage4_regulon_size) &
      !is.na(stage5_regulon_size) &
      stage4_regulon_size ==
        stage5_regulon_size
  )
]

write_csv_safe(
  regulon_integrity_audit,
  file.path(
    DIRS$tables,
    "00_stage5_regulon_integrity_audit.csv"
  )
)

if (any(regulon_integrity_audit$target_count_match != TRUE)) {
  failed_tfs <- regulon_integrity_audit[
    target_count_match != TRUE,
    source_symbol
  ]
  stop(
    "Stage 4 network target counts do not match the locked ",
    "regulon-size table for: ",
    paste(failed_tfs, collapse = ", ")
  )
}

pseudobulk_objects <- readRDS(
  STAGE4_PSEUDOBULK_RDS
)

required_pb_objects <- c(
  "sample_logcpm",
  "state_logcpm"
)

missing_pb_objects <- setdiff(
  required_pb_objects,
  names(pseudobulk_objects)
)

if (length(missing_pb_objects) > 0L) {
  stop(
    "Stage 4 pseudobulk RDS is missing object(s): ",
    paste(missing_pb_objects, collapse = ", ")
  )
}

sample_logcpm <- as.matrix(
  pseudobulk_objects$sample_logcpm
)
state_logcpm <- as.matrix(
  pseudobulk_objects$state_logcpm
)

weighted_sample_activity <- as.matrix(
  readRDS(STAGE4_ACTIVITY_RDS)
)

log_msg("Loading Stage 4 macrophage Seurat object.")
macrophage <- readRDS(
  STAGE4_MACROPHAGE_RDS
)

if (!inherits(macrophage, "Seurat")) {
  stop(
    "Stage 4 macrophage RDS is not a Seurat object."
  )
}

if (ncol(macrophage) != 1822L) {
  stop(
    "Unexpected Stage 4 macrophage cell count. Expected 1822; observed ",
    ncol(macrophage)
  )
}

sample_order <- sample_meta$sample_accession

if (
  !all(sample_order %in% colnames(sample_logcpm)) ||
  !all(sample_order %in% colnames(weighted_sample_activity))
) {
  stop(
    "Not all six locked samples are present in Stage 4 matrices."
  )
}

sample_logcpm <- sample_logcpm[
  ,
  sample_order,
  drop = FALSE
]

weighted_sample_activity <-
  weighted_sample_activity[
    ,
    sample_order,
    drop = FALSE
  ]

feature_map <- make_feature_map(
  rownames(sample_logcpm)
)

############################################################
## 5. Resolve candidates and select matched TF controls
############################################################

network_tf_map <- unique(
  network_dt[
    ,
    .(
      tf_key = gene_key(source_symbol),
      tf_symbol = source_symbol
    )
  ],
  by = "tf_key"
)

candidate_manifest <- merge(
  data.table::data.table(
    requested_tf = CANDIDATE_TFS_REQUESTED,
    requested_order = seq_along(
      CANDIDATE_TFS_REQUESTED
    ),
    tf_key = gene_key(
      CANDIDATE_TFS_REQUESTED
    )
  ),
  network_tf_map,
  by = "tf_key",
  all.x = TRUE,
  sort = FALSE
)

data.table::setorder(
  candidate_manifest,
  requested_order
)

missing_candidates <- candidate_manifest[
  is.na(tf_symbol),
  requested_tf
]

if (length(missing_candidates) > 0L) {
  stop(
    "Candidate TF(s) are absent from the Stage 4 network: ",
    paste(missing_candidates, collapse = ", ")
  )
}

candidate_tfs <- candidate_manifest$tf_symbol

if (
  any(!candidate_tfs %in%
    rownames(weighted_sample_activity))
) {
  stop(
    "At least one candidate TF is absent from the Stage 4 activity matrix: ",
    paste(
      setdiff(
        candidate_tfs,
        rownames(weighted_sample_activity)
      ),
      collapse = ", "
    )
  )
}

matched_controls <- select_matched_control_tfs(
  priority_dt,
  network_dt,
  candidate_tfs
)

control_tfs <- matched_controls$tf_symbol
analysis_tfs <- c(candidate_tfs, control_tfs)

analysis_integrity <- regulon_integrity_audit[
  source_symbol %in% analysis_tfs
]

if (
  nrow(analysis_integrity) != length(analysis_tfs) ||
  any(analysis_integrity$target_count_match != TRUE)
) {
  stop(
    "Candidate/control TF-specific network integrity failed."
  )
}

analysis_tf_manifest <- merge(
  data.table::data.table(
    tf_symbol = analysis_tfs,
    analysis_role = c(
      rep("Biological_candidate", length(candidate_tfs)),
      rep("Matched_low_priority_control", length(control_tfs))
    ),
    requested_order = seq_along(analysis_tfs)
  ),
  priority_dt,
  by = "tf_symbol",
  all.x = TRUE,
  sort = FALSE
)

analysis_tf_manifest <- merge(
  analysis_tf_manifest,
  network_dt[
    ,
    .(
      regulon_size =
        data.table::uniqueN(target_feature)
    ),
    by = source_symbol
  ],
  by.x = "tf_symbol",
  by.y = "source_symbol",
  all.x = TRUE,
  sort = FALSE
)

data.table::setorder(
  analysis_tf_manifest,
  requested_order
)

write_csv_safe(
  candidate_manifest,
  file.path(
    DIRS$tables,
    "01_stage5_candidate_TF_resolution.csv"
  )
)

write_csv_safe(
  analysis_tf_manifest,
  file.path(
    DIRS$tables,
    "02_stage5_candidate_and_matched_control_TFs.csv"
  )
)

log_msg(
  "Biological candidates: ",
  paste(candidate_tfs, collapse = ", ")
)
log_msg(
  "Matched controls: ",
  paste(control_tfs, collapse = ", ")
)

############################################################
## 6. Program definitions and observed scores
############################################################

program_objects <- make_program_definitions(
  program_manifest,
  feature_map
)

program_definitions <- program_objects$definitions
program_summary <- program_objects$summary

if (length(program_definitions) < 10L) {
  stop(
    "Too few evaluable transcriptional programs: ",
    length(program_definitions)
  )
}

sample_z <- scale_rows(sample_logcpm)

observed_program_scores <- score_programs(
  sample_z,
  program_definitions
)

program_observed_summary <- data.table::rbindlist(
  lapply(
    rownames(observed_program_scores),
    function(program_id) {
      control_values <- observed_program_scores[
        program_id,
        sample_meta$condition == "Control"
      ]
      hfpef_values <- observed_program_scores[
        program_id,
        sample_meta$condition == "HFpEF"
      ]

      data.table::data.table(
        program_id = program_id,
        control_mean = mean(
          control_values,
          na.rm = TRUE
        ),
        hfpef_mean = mean(
          hfpef_values,
          na.rm = TRUE
        ),
        hfpef_minus_control =
          mean(hfpef_values, na.rm = TRUE) -
          mean(control_values, na.rm = TRUE),
        wilcoxon_p = safe_wilcox_p(
          hfpef_values,
          control_values
        ),
        observed_gap_eligible = abs(
          mean(hfpef_values, na.rm = TRUE) -
            mean(control_values, na.rm = TRUE)
        ) >= MIN_ABS_OBSERVED_PROGRAM_GAP
      )
    }
  )
)

program_observed_summary[, wilcoxon_fdr := p.adjust(
  wilcoxon_p,
  method = "BH"
)]

program_observed_summary <- merge(
  program_observed_summary,
  program_summary,
  by = "program_id",
  all.x = TRUE
)

eligible_primary_stage2_programs <-
  program_observed_summary[
    program_category ==
      "Stage2_drug_opposed_net" &
      primary_program == TRUE &
      observed_gap_eligible == TRUE,
    data.table::uniqueN(program_id)
  ]

if (eligible_primary_stage2_programs < 2L) {
  stop(
    "Fewer than two primary Stage 2 drug-opposed programs have ",
    "a stable observed HFpEF-Control gap for perturbation testing. ",
    "Observed count: ",
    eligible_primary_stage2_programs
  )
}

write_csv_safe(
  program_summary,
  file.path(
    DIRS$tables,
    "03_stage5_program_definition_summary.csv"
  )
)

write_csv_safe(
  program_observed_summary,
  file.path(
    DIRS$tables,
    "04_stage5_observed_program_scores.csv"
  )
)

write_csv_safe(
  data.table::as.data.table(
    observed_program_scores,
    keep.rownames = "program_id"
  ),
  file.path(
    DIRS$source,
    "observed_sample_level_program_scores.csv"
  )
)

############################################################
## 7. Sample-level comparative perturbation
############################################################

program_effect_records <- list()
gene_effect_records <- list()
activity_audit_records <- list()

total_runs <- length(analysis_tfs) *
  length(PERTURBATION_MODES) *
  length(PERTURBATION_STRENGTHS) *
  length(PERTURBATION_METHODS)

run_index <- 0L

for (tf_symbol in analysis_tfs) {
  for (perturbation_mode in PERTURBATION_MODES) {
    for (
      perturbation_strength in
        PERTURBATION_STRENGTHS
    ) {
      for (
        perturbation_method in
          PERTURBATION_METHODS
      ) {
        run_index <- run_index + 1L

        log_msg(
          "Perturbation ",
          run_index,
          "/",
          total_runs,
          " | TF=",
          tf_symbol,
          " | mode=",
          perturbation_mode,
          " | strength=",
          perturbation_strength,
          " | method=",
          perturbation_method
        )

        sim <- simulate_tf_activity_adjustment(
          observed_z = sample_z,
          activity_matrix =
            weighted_sample_activity,
          tf_requested = tf_symbol,
          network_dt = network_dt,
          sample_meta = sample_meta,
          perturbation_mode =
            perturbation_mode,
          perturbation_strength =
            perturbation_strength,
          perturbation_method =
            perturbation_method
        )

        perturbed_program_scores <- score_programs(
          sim$perturbed_z,
          program_definitions
        )

        effect_dt <- summarize_program_effects(
          observed_scores =
            observed_program_scores,
          perturbed_scores =
            perturbed_program_scores,
          sample_meta = sample_meta,
          program_summary =
            program_summary,
          tf_symbol = tf_symbol,
          perturbation_mode =
            perturbation_mode,
          perturbation_strength =
            perturbation_strength,
          perturbation_method =
            perturbation_method,
          target_count =
            sim$target_count,
          global_rms_shift =
            sim$global_rms_shift_HFpEF,
          delta_activity =
            sim$delta_activity,
          realized_delta_activity =
            sim$realized_delta_activity
        )

        program_effect_records[[
          length(program_effect_records) + 1L
        ]] <- effect_dt

        gene_effect_records[[
          length(gene_effect_records) + 1L
        ]] <- sim$gene_effect

        activity_audit_records[[
          length(activity_audit_records) + 1L
        ]] <- data.table::data.table(
          tf_symbol = tf_symbol,
          perturbation_mode =
            perturbation_mode,
          perturbation_strength =
            perturbation_strength,
          perturbation_method =
            perturbation_method,
          sample_accession =
            colnames(sample_z),
          condition = as.character(
            sample_meta$condition
          ),
          observed_activity =
            sim$observed_activity,
          desired_activity =
            sim$desired_activity,
          requested_delta_activity =
            sim$delta_activity,
          realized_delta_activity =
            sim$realized_delta_activity,
          control_reference =
            sim$control_reference,
          target_count =
            sim$target_count
        )
      }
    }
  }
}

program_effects <- data.table::rbindlist(
  program_effect_records,
  use.names = TRUE,
  fill = TRUE
)

gene_effects <- data.table::rbindlist(
  gene_effect_records,
  use.names = TRUE,
  fill = TRUE
)

activity_audit <- data.table::rbindlist(
  activity_audit_records,
  use.names = TRUE,
  fill = TRUE
)

program_effects <- merge(
  program_effects,
  analysis_tf_manifest[
    ,
    .(
      tf_symbol,
      analysis_role,
      stage4_priority_rank =
        priority_rank,
      stage4_priority_score =
        priority_score
    )
  ],
  by = "tf_symbol",
  all.x = TRUE
)

gene_effects <- merge(
  gene_effects,
  analysis_tf_manifest[
    ,
    .(
      tf_symbol,
      analysis_role,
      stage4_priority_rank =
        priority_rank
    )
  ],
  by = "tf_symbol",
  all.x = TRUE
)

write_csv_safe(
  program_effects,
  file.path(
    DIRS$tables,
    "05_stage5_sample_level_program_perturbation_results.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  activity_audit,
  file.path(
    DIRS$tables,
    "06_stage5_requested_and_realized_TF_activity_changes.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  gene_effects,
  file.path(
    DIRS$tables,
    "07_stage5_all_predicted_target_gene_changes.csv"
  ),
  compress = TRUE
)

############################################################
## 8. Gene-level reporting and ligand preparation
############################################################

ligand_panel <- unique(c(
  "Tnf", "Il1b", "Il18", "Ccl2", "Ccl3",
  "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Ccl12",
  "Cxcl1", "Cxcl2", "Cxcl9", "Cxcl10",
  "Spp1", "Tgfb1", "Vegfa", "Vegfb",
  "Pdgfa", "Pdgfb", "Apoe", "Lgals3",
  "Osm", "Il6", "Il10", "Il15", "Il27",
  "Kitl", "Hbegf", "Areg", "Ereg",
  "Gas6", "Pros1", "Mif", "Nampt",
  "Fn1", "Thbs1", "Csf1", "Csf2",
  "Csf3", "Bmp2", "Wnt5a"
))

ligand_keys <- gene_key(ligand_panel)

gene_effects[, is_candidate_ligand := (
  target_key %in% ligand_keys
)]

primary_gene_effects <- gene_effects[
  perturbation_mode == PRIMARY_MODE &
    perturbation_strength == PRIMARY_STRENGTH
]

primary_gene_effects[
  ,
  absolute_mean_delta_z :=
    abs(mean_delta_z_HFpEF)
]

data.table::setorder(
  primary_gene_effects,
  tf_symbol,
  perturbation_method,
  -absolute_mean_delta_z
)

top_gene_effects <- primary_gene_effects[
  ,
  head(.SD, TOP_GENES_PER_TF_REPORT),
  by = .(
    tf_symbol,
    analysis_role,
    perturbation_method
  )
]

ligand_effects <- primary_gene_effects[
  is_candidate_ligand == TRUE
]

data.table::setorder(
  ligand_effects,
  tf_symbol,
  perturbation_method,
  -absolute_mean_delta_z
)

ligand_effects <- ligand_effects[
  ,
  head(.SD, TOP_LIGANDS_PER_TF_REPORT),
  by = .(
    tf_symbol,
    analysis_role,
    perturbation_method
  )
]

write_csv_safe(
  top_gene_effects,
  file.path(
    DIRS$tables,
    "08_stage5_primary_top_predicted_genes_per_TF.csv"
  )
)

write_csv_safe(
  ligand_effects,
  file.path(
    DIRS$tables,
    "09_stage5_candidate_ligand_changes_for_stage6.csv"
  )
)

############################################################
## 9. Macrophage-state-level perturbation
############################################################

state_profile_meta <- parse_state_profile_names(
  colnames(state_logcpm)
)

state_profile_meta <- merge(
  state_profile_meta,
  sample_meta[
    ,
    .(
      sample_accession,
      condition
    )
  ],
  by = "sample_accession",
  all.x = TRUE
)

state_profile_meta <- state_profile_meta[
  match(
    colnames(state_logcpm),
    profile
  )
]

if (any(is.na(state_profile_meta$condition))) {
  stop(
    "State-level profiles could not be mapped to sample conditions."
  )
}

state_z <- scale_rows(state_logcpm)

state_activity <- weighted_regulon_activity(
  state_logcpm,
  network_dt,
  tf_subset = analysis_tfs
)

missing_state_candidates <- setdiff(
  candidate_tfs,
  rownames(state_activity)
)

if (length(missing_state_candidates) > 0L) {
  stop(
    "Candidate TF activity could not be reconstructed in the ",
    "sample-state pseudobulk matrix: ",
    paste(missing_state_candidates, collapse = ", ")
  )
}

state_program_scores <- score_programs(
  state_z,
  program_definitions
)

eligible_states <- state_profile_meta[
  ,
  .(
    control_profiles = sum(
      condition == "Control"
    ),
    hfpef_profiles = sum(
      condition == "HFpEF"
    )
  ),
  by = macrophage_state
][
  control_profiles >=
    MIN_STATE_PROFILES_PER_CONDITION &
    hfpef_profiles >=
      MIN_STATE_PROFILES_PER_CONDITION,
  macrophage_state
]

state_effect_records <- list()

for (state_name in eligible_states) {
  state_profiles <- state_profile_meta[
    macrophage_state == state_name,
    profile
  ]

  state_meta_i <- state_profile_meta[
    macrophage_state == state_name
  ]
  state_meta_i <- state_meta_i[
    match(state_profiles, profile)
  ]

  state_z_i <- state_z[
    ,
    state_profiles,
    drop = FALSE
  ]

  state_activity_i <- state_activity[
    ,
    state_profiles,
    drop = FALSE
  ]

  state_scores_i <- state_program_scores[
    ,
    state_profiles,
    drop = FALSE
  ]

  state_sample_meta_i <- data.table::data.table(
    sample_accession =
      state_meta_i$profile,
    condition = factor(
      state_meta_i$condition,
      levels = c("Control", "HFpEF")
    )
  )

  for (tf_symbol in candidate_tfs) {
    for (
      perturbation_strength in
        c(0.50, 1.00)
    ) {
      for (
        perturbation_method in
          PERTURBATION_METHODS
      ) {
        if (
          !tf_symbol %in%
            rownames(state_activity_i)
        ) {
          next
        }

        sim_state <- simulate_tf_activity_adjustment(
          observed_z = state_z_i,
          activity_matrix =
            state_activity_i,
          tf_requested = tf_symbol,
          network_dt = network_dt,
          sample_meta =
            state_sample_meta_i,
          perturbation_mode =
            "disease_normalization",
          perturbation_strength =
            perturbation_strength,
          perturbation_method =
            perturbation_method
        )

        perturbed_state_scores <- score_programs(
          sim_state$perturbed_z,
          program_definitions
        )

        state_effect <- summarize_program_effects(
          observed_scores =
            state_scores_i,
          perturbed_scores =
            perturbed_state_scores,
          sample_meta =
            state_sample_meta_i,
          program_summary =
            program_summary,
          tf_symbol = tf_symbol,
          perturbation_mode =
            "disease_normalization",
          perturbation_strength =
            perturbation_strength,
          perturbation_method =
            perturbation_method,
          target_count =
            sim_state$target_count,
          global_rms_shift =
            sim_state$global_rms_shift_HFpEF,
          delta_activity =
            sim_state$delta_activity,
          realized_delta_activity =
            sim_state$realized_delta_activity
        )

        state_effect[, macrophage_state := state_name]

        state_effect_records[[
          length(state_effect_records) + 1L
        ]] <- state_effect
      }
    }
  }
}

state_program_effects <- if (
  length(state_effect_records) > 0L
) {
  data.table::rbindlist(
    state_effect_records,
    use.names = TRUE,
    fill = TRUE
  )
} else {
  data.table::data.table(
    tf_symbol = character(),
    perturbation_mode = character(),
    perturbation_strength = numeric(),
    perturbation_method = character(),
    program_id = character(),
    program_category = character(),
    subset_name = character(),
    primary_program = logical(),
    observed_gap_eligible = logical(),
    absolute_gap_reduction = numeric(),
    recovery_fraction = numeric(),
    sample_improvement_fraction = numeric(),
    macrophage_state = character()
  )
}

write_csv_safe(
  state_program_effects,
  file.path(
    DIRS$tables,
    "10_stage5_macrophage_state_specific_perturbation_results.csv"
  ),
  compress = TRUE
)

############################################################
## 10. Method concordance and candidate summaries
############################################################

primary_results <- program_effects[
  perturbation_mode == PRIMARY_MODE &
    perturbation_strength ==
      PRIMARY_STRENGTH
]

method_wide <- data.table::dcast(
  primary_results[
    analysis_role == "Biological_candidate"
  ],
  tf_symbol + program_id ~
    perturbation_method,
  value.var = "absolute_gap_reduction"
)

method_concordance <- method_wide[
  ,
  .(
    programs_compared = sum(
      is.finite(weighted_minimum_norm) &
        is.finite(equal_signed_targets)
    ),
    spearman_recovery =
      safe_spearman(
        weighted_minimum_norm,
        equal_signed_targets
      ),
    direction_agreement = mean(
      sign(weighted_minimum_norm) ==
        sign(equal_signed_targets),
      na.rm = TRUE
    )
  ),
  by = tf_symbol
]

candidate_summary_records <- lapply(
  candidate_tfs,
  function(tf_symbol) {
    tf_i <- tf_symbol

    tf_primary <- primary_results[
      tf_symbol == tf_i
    ]

    tf_stage2_primary <- tf_primary[
      program_category ==
        "Stage2_drug_opposed_net" &
        primary_program == TRUE &
        observed_gap_eligible == TRUE
    ]

    tf_stage2_all <- tf_primary[
      program_category ==
        "Stage2_drug_opposed_net" &
        observed_gap_eligible == TRUE
    ]

    tf_functional <- tf_primary[
      program_category ==
        "Functional_state" &
        observed_gap_eligible == TRUE
    ]

    inflammatory_names <- c(
      "Inflammatory_Il1b",
      "NFkB_TNF_response",
      "Inflammasome_pyroptosis",
      "Interferon_response"
    )

    tf_inflammation <- tf_functional[
      subset_name %in%
        inflammatory_names
    ]

    tf_myeloid <- tf_functional[
      subset_name == "Myeloid_identity"
    ]

    tf_state <- state_program_effects[
      tf_symbol == tf_i &
        perturbation_strength ==
          PRIMARY_STRENGTH &
        program_category ==
          "Stage2_drug_opposed_net" &
        primary_program == TRUE &
        observed_gap_eligible == TRUE
    ]

    method_row <- method_concordance[
      tf_symbol == tf_i
    ]

    data.table::data.table(
      tf_symbol = tf_symbol,
      stage2_primary_median_gap_reduction =
        safe_median(
          tf_stage2_primary$absolute_gap_reduction,
          0
        ),
      stage2_primary_median_recovery_fraction =
        safe_median(
          tf_stage2_primary$recovery_fraction,
          0
        ),
      stage2_primary_positive_fraction =
        safe_mean(
          tf_stage2_primary$absolute_gap_reduction > 0,
          0
        ),
      stage2_allsize_positive_fraction =
        safe_mean(
          tf_stage2_all$absolute_gap_reduction > 0,
          0
        ),
      biological_sample_improvement_fraction =
        safe_median(
          tf_stage2_primary$
            sample_improvement_fraction,
          0
        ),
      inflammation_median_gap_reduction =
        safe_median(
          tf_inflammation$
            absolute_gap_reduction,
          0
        ),
      state_primary_positive_fraction =
        safe_mean(
          tf_state$absolute_gap_reduction > 0,
          0
        ),
      state_sample_improvement_fraction =
        safe_median(
          tf_state$
            sample_improvement_fraction,
          0
        ),
      method_spearman =
        if (nrow(method_row) > 0L) {
          method_row$spearman_recovery[1L]
        } else {
          NA_real_
        },
      method_direction_agreement =
        if (nrow(method_row) > 0L) {
          method_row$direction_agreement[1L]
        } else {
          NA_real_
        },
      median_global_rms_shift =
        safe_median(
          tf_primary$global_rms_shift_HFpEF,
          NA_real_
        ),
      myeloid_identity_absolute_change =
        safe_median(
          abs(tf_myeloid$hfpef_score_change),
          0
        )
    )
  }
)

candidate_summary <- data.table::rbindlist(
  candidate_summary_records,
  use.names = TRUE,
  fill = TRUE
)

control_primary_summary <- primary_results[
  analysis_role ==
    "Matched_low_priority_control" &
    program_category ==
      "Stage2_drug_opposed_net" &
    primary_program == TRUE &
    observed_gap_eligible == TRUE,
  .(
    control_median_gap_reduction =
      safe_median(
        absolute_gap_reduction,
        0
      ),
    control_positive_fraction =
      safe_mean(
        absolute_gap_reduction > 0,
        0
      )
  ),
  by = tf_symbol
]

control_reference_recovery <- safe_median(
  control_primary_summary$
    control_median_gap_reduction,
  0
)

candidate_summary[
  ,
  matched_control_recovery_margin :=
    stage2_primary_median_gap_reduction -
    control_reference_recovery
]

candidate_summary <- merge(
  candidate_summary,
  priority_dt[
    tf_symbol %in% candidate_tfs,
    .(
      tf_symbol,
      stage4_priority_rank =
        priority_rank,
      stage4_priority_score =
        priority_score,
      stage4_weighted_effect =
        weighted_effect,
      stage4_weighted_hedges_g =
        weighted_hedges_g,
      stage4_aucell_effect =
        aucell_effect,
      stage4_supported_overlap_count =
        supported_overlap_count
    )
  ],
  by = "tf_symbol",
  all.x = TRUE
)

candidate_summary[
  ,
  specificity_score :=
    1 /
    (
      1 +
        pmax(
          median_global_rms_shift,
          0
        ) +
        pmax(
          myeloid_identity_absolute_change,
          0
        )
    )
]

candidate_summary[
  ,
  rank_stage2_recovery := rank_metric(
    stage2_primary_median_gap_reduction,
    higher_is_better = TRUE
  )
]
candidate_summary[
  ,
  rank_cross_size_robustness := rank_metric(
    stage2_allsize_positive_fraction,
    higher_is_better = TRUE
  )
]
candidate_summary[
  ,
  rank_sample_stability := rank_metric(
    biological_sample_improvement_fraction,
    higher_is_better = TRUE
  )
]
candidate_summary[
  ,
  rank_state_robustness := rank_metric(
    state_primary_positive_fraction,
    higher_is_better = TRUE
  )
]
candidate_summary[
  ,
  rank_method_concordance := rank_metric(
    method_spearman,
    higher_is_better = TRUE
  )
]
candidate_summary[
  ,
  rank_inflammation_recovery := rank_metric(
    inflammation_median_gap_reduction,
    higher_is_better = TRUE
  )
]
candidate_summary[
  ,
  rank_specificity := rank_metric(
    specificity_score,
    higher_is_better = TRUE
  )
]
candidate_summary[
  ,
  rank_control_margin := rank_metric(
    matched_control_recovery_margin,
    higher_is_better = TRUE
  )
]
candidate_summary[
  ,
  rank_stage4_support := rank_metric(
    -stage4_priority_rank,
    higher_is_better = TRUE
  )
]

rank_columns <- c(
  "rank_stage2_recovery",
  "rank_cross_size_robustness",
  "rank_sample_stability",
  "rank_state_robustness",
  "rank_method_concordance",
  "rank_inflammation_recovery",
  "rank_specificity",
  "rank_control_margin",
  "rank_stage4_support"
)

candidate_summary[
  ,
  mean_evidence_rank := rowMeans(
    .SD,
    na.rm = TRUE
  ),
  .SDcols = rank_columns
]

candidate_summary[
  ,
  median_evidence_rank := apply(
    .SD,
    1L,
    stats::median,
    na.rm = TRUE
  ),
  .SDcols = rank_columns
]

candidate_summary[
  ,
  final_rank_aggregation_score :=
    0.5 * mean_evidence_rank +
    0.5 * median_evidence_rank
]

data.table::setorder(
  candidate_summary,
  final_rank_aggregation_score,
  rank_stage2_recovery,
  rank_specificity
)

candidate_summary[, final_candidate_rank := seq_len(.N)]
candidate_summary[, Nfkb1_forced := FALSE]

write_csv_safe(
  method_concordance,
  file.path(
    DIRS$tables,
    "11_stage5_perturbation_method_concordance.csv"
  )
)

write_csv_safe(
  control_primary_summary,
  file.path(
    DIRS$tables,
    "12_stage5_matched_control_TF_reference_results.csv"
  )
)

write_csv_safe(
  candidate_summary,
  file.path(
    DIRS$tables,
    "13_stage5_candidate_TF_rank_aggregation.csv"
  )
)

############################################################
## 11. Ranking sensitivity scenarios
############################################################

ranking_scenarios <- list(
  all_evidence = rank_columns,
  without_stage4_prior = setdiff(
    rank_columns,
    "rank_stage4_support"
  ),
  program_and_robustness = c(
    "rank_stage2_recovery",
    "rank_cross_size_robustness",
    "rank_sample_stability",
    "rank_state_robustness",
    "rank_method_concordance"
  ),
  specificity_emphasis = c(
    "rank_stage2_recovery",
    "rank_specificity",
    "rank_control_margin",
    "rank_method_concordance"
  ),
  inflammatory_emphasis = c(
    "rank_stage2_recovery",
    "rank_inflammation_recovery",
    "rank_sample_stability",
    "rank_method_concordance"
  )
)

scenario_records <- list()

for (scenario_name in names(ranking_scenarios)) {
  cols <- ranking_scenarios[[scenario_name]]

  scenario_dt <- candidate_summary[
    ,
    .(
      tf_symbol,
      scenario_score = rowMeans(
        .SD,
        na.rm = TRUE
      )
    ),
    .SDcols = cols
  ]

  data.table::setorder(
    scenario_dt,
    scenario_score
  )

  scenario_dt[, scenario_rank := seq_len(.N)]
  scenario_dt[, scenario := scenario_name]
  scenario_dt[, metrics_used := paste(cols, collapse = ";")]

  scenario_records[[
    length(scenario_records) + 1L
  ]] <- scenario_dt
}

ranking_sensitivity <- data.table::rbindlist(
  scenario_records,
  use.names = TRUE,
  fill = TRUE
)

ranking_stability <- ranking_sensitivity[
  ,
  .(
    scenarios = .N,
    median_scenario_rank =
      stats::median(scenario_rank),
    best_scenario_rank = min(scenario_rank),
    worst_scenario_rank = max(scenario_rank),
    top3_frequency = mean(scenario_rank <= 3L)
  ),
  by = tf_symbol
]

write_csv_safe(
  ranking_sensitivity,
  file.path(
    DIRS$tables,
    "14_stage5_candidate_ranking_sensitivity_scenarios.csv"
  )
)

write_csv_safe(
  ranking_stability,
  file.path(
    DIRS$tables,
    "15_stage5_candidate_ranking_stability_summary.csv"
  )
)

############################################################
## 12. Normalization-versus-attenuation directionality
############################################################

mode_comparison <- program_effects[
  tf_symbol %in% candidate_tfs &
    perturbation_strength ==
      PRIMARY_STRENGTH &
    perturbation_method ==
      PRIMARY_METHOD &
    program_category ==
      "Stage2_drug_opposed_net" &
    primary_program == TRUE &
    observed_gap_eligible == TRUE,
  .(
    median_gap_reduction =
      safe_median(
        absolute_gap_reduction,
        0
      ),
    positive_fraction =
      safe_mean(
        absolute_gap_reduction > 0,
        0
      ),
    sample_improvement_fraction =
      safe_median(
        sample_improvement_fraction,
        0
      )
  ),
  by = .(
    tf_symbol,
    perturbation_mode
  )
]

mode_wide <- data.table::dcast(
  mode_comparison,
  tf_symbol ~ perturbation_mode,
  value.var = "median_gap_reduction"
)

if (
  all(
    c(
      "disease_normalization",
      "activity_attenuation"
    ) %in% names(mode_wide)
  )
) {
  mode_wide[
    ,
    normalization_advantage :=
      disease_normalization -
      activity_attenuation
  ]
}

write_csv_safe(
  mode_comparison,
  file.path(
    DIRS$tables,
    "16_stage5_normalization_vs_attenuation_results.csv"
  )
)

write_csv_safe(
  mode_wide,
  file.path(
    DIRS$tables,
    "17_stage5_normalization_vs_attenuation_summary.csv"
  )
)

############################################################
## 13. Save scientific checkpoint before figures
############################################################

scientific_checkpoint <- list(
  analysis_tf_manifest =
    analysis_tf_manifest,
  program_summary =
    program_summary,
  observed_program_scores =
    observed_program_scores,
  program_effects =
    program_effects,
  gene_effects =
    gene_effects,
  activity_audit =
    activity_audit,
  state_program_effects =
    state_program_effects,
  method_concordance =
    method_concordance,
  candidate_summary =
    candidate_summary,
  ranking_sensitivity =
    ranking_sensitivity,
  ranking_stability =
    ranking_stability,
  mode_comparison =
    mode_comparison,
  ligand_effects =
    ligand_effects
)

saveRDS(
  scientific_checkpoint,
  file.path(
    DIRS$objects,
    "CHECKPOINT_stage5_scientific_results_pre_figures.rds"
  ),
  compress = FALSE
)

log_msg(
  "Stage 5 scientific calculations completed and checkpointed before figures."
)

############################################################
## 14. Figures
############################################################

candidate_order <- candidate_summary[
  order(final_candidate_rank),
  tf_symbol
]

## Figure 5A: primary program recovery heatmap.
heat_dt <- primary_results[
  tf_symbol %in% candidate_tfs &
    program_category ==
      "Stage2_drug_opposed_net" &
    primary_program == TRUE &
    observed_gap_eligible == TRUE,
  .(
    median_gap_reduction =
      safe_median(
        absolute_gap_reduction,
        0
      )
  ),
  by = .(
    tf_symbol,
    subset_name,
    support_class
  )
]

heat_dt[, program_label := paste(
  subset_name,
  support_class,
  sep = "__"
)]

heat_wide <- data.table::dcast(
  heat_dt,
  tf_symbol ~ program_label,
  value.var = "median_gap_reduction",
  fill = 0
)

heat_wide <- heat_wide[
  match(candidate_order, tf_symbol)
]

heat_matrix <- as.matrix(
  heat_wide[
    ,
    -1L,
    with = FALSE
  ]
)
rownames(heat_matrix) <- heat_wide$tf_symbol
heat_matrix[!is.finite(heat_matrix)] <- 0

write_csv_safe(
  heat_wide,
  file.path(
    DIRS$source,
    "Fig5A_primary_program_recovery_heatmap_source.csv"
  )
)

save_heatmap_bundle(
  heat_matrix,
  "Fig5A_multiTF_primary_program_recovery_heatmap",
  width = 10,
  height = 7,
  main =
    "Disease-normalizing perturbation: primary program recovery"
)

## Figure 5B: final rank aggregation.
rank_plot_dt <- data.table::copy(
  candidate_summary
)

rank_plot_dt[, tf_symbol := factor(
  tf_symbol,
  levels = rev(candidate_order)
)]

write_csv_safe(
  rank_plot_dt,
  file.path(
    DIRS$source,
    "Fig5B_candidate_rank_aggregation_source.csv"
  )
)

p_rank <- ggplot2::ggplot(
  rank_plot_dt,
  ggplot2::aes(
    x = -final_rank_aggregation_score,
    y = tf_symbol
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = -max(final_rank_aggregation_score),
      xend = -final_rank_aggregation_score,
      y = tf_symbol,
      yend = tf_symbol
    ),
    color = "grey65",
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = stage2_primary_positive_fraction,
      fill =
        stage2_primary_median_gap_reduction
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
  ggplot2::scale_size_continuous(
    limits = c(0, 1),
    range = c(3, 8)
  ) +
  ggplot2::labs(
    title =
      "Comparative multi-TF perturbation prioritization",
    subtitle =
      "Rank aggregation across program recovery, robustness, method agreement, state effects, specificity, and Stage 4 support",
    x =
      "Higher rank-aggregation performance",
    y = NULL,
    fill =
      "Median primary\nprogram recovery",
    size =
      "Positive primary\nprogram fraction"
  ) +
  ggplot2::theme_bw(base_size = 10)

save_plot_bundle(
  p_rank,
  "Fig5B_multiTF_rank_aggregation",
  9,
  6.5
)

## Figure 5C: disease normalization versus activity attenuation.
if (
  all(
    c(
      "disease_normalization",
      "activity_attenuation"
    ) %in% names(mode_wide)
  )
) {
  mode_plot <- merge(
    mode_wide,
    candidate_summary[
      ,
      .(
        tf_symbol,
        final_candidate_rank
      )
    ],
    by = "tf_symbol",
    all.x = TRUE
  )

  write_csv_safe(
    mode_plot,
    file.path(
      DIRS$source,
      "Fig5C_normalization_vs_attenuation_source.csv"
    )
  )

  p_mode <- ggplot2::ggplot(
    mode_plot,
    ggplot2::aes(
      x = activity_attenuation,
      y = disease_normalization,
      label = tf_symbol
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = final_candidate_rank
      ),
      shape = 21,
      size = 4,
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      size = 3.5,
      max.overlaps = Inf
    ) +
    ggplot2::scale_fill_gradient(
      low = "#D73027",
      high = "#4575B4",
      trans = "reverse"
    ) +
    ggplot2::labs(
      title =
        "Disease normalization is not equivalent to uniform TF inhibition",
      subtitle =
        "Positive values indicate reduced HFpEF–Control program separation",
      x =
        "Activity attenuation: median primary-program recovery",
      y =
        "Disease normalization: median primary-program recovery",
      fill = "Final rank"
    ) +
    ggplot2::theme_bw(base_size = 10)

  save_plot_bundle(
    p_mode,
    "Fig5C_normalization_vs_activity_attenuation",
    8,
    7
  )
}

## Figure 5D: macrophage-state-specific recovery.
if (nrow(state_program_effects) > 0L) {
  state_plot_dt <- state_program_effects[
    tf_symbol %in% candidate_tfs &
      perturbation_strength ==
        PRIMARY_STRENGTH &
      program_category ==
        "Stage2_drug_opposed_net" &
      primary_program == TRUE &
      observed_gap_eligible == TRUE,
    .(
      median_gap_reduction =
        safe_median(
          absolute_gap_reduction,
          0
        )
    ),
    by = .(
      tf_symbol,
      macrophage_state
    )
  ]

  state_wide <- data.table::dcast(
    state_plot_dt,
    tf_symbol ~ macrophage_state,
    value.var = "median_gap_reduction",
    fill = 0
  )

  state_wide <- state_wide[
    match(candidate_order, tf_symbol)
  ]

  state_matrix <- as.matrix(
    state_wide[
      ,
      -1L,
      with = FALSE
    ]
  )
  rownames(state_matrix) <- state_wide$tf_symbol
  state_matrix[!is.finite(state_matrix)] <- 0

  write_csv_safe(
    state_wide,
    file.path(
      DIRS$source,
      "Fig5D_state_specific_recovery_source.csv"
    )
  )

  save_heatmap_bundle(
    state_matrix,
    "Fig5D_macrophage_state_specific_program_recovery",
    width = 10,
    height = 7,
    main =
      "State-specific disease-program recovery"
  )
}

## Figure 5E: predicted candidate-ligand changes.
ligand_plot_dt <- ligand_effects[
  tf_symbol %in% candidate_tfs &
    perturbation_method ==
      PRIMARY_METHOD
]

if (nrow(ligand_plot_dt) > 0L) {
  ligand_plot_dt <- ligand_plot_dt[
    ,
    head(
      .SD[
        order(-absolute_mean_delta_z)
      ],
      12L
    ),
    by = tf_symbol
  ]

  ligand_plot_dt[, tf_symbol := factor(
    tf_symbol,
    levels = candidate_order
  )]

  write_csv_safe(
    ligand_plot_dt,
    file.path(
      DIRS$source,
      "Fig5E_candidate_ligand_change_source.csv"
    )
  )

  p_ligand <- ggplot2::ggplot(
    ligand_plot_dt,
    ggplot2::aes(
      x = target_feature,
      y = tf_symbol
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = absolute_mean_delta_z,
        fill = mean_delta_z_HFpEF
      ),
      shape = 21,
      color = "black",
      stroke = 0.25
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0
    ) +
    ggplot2::scale_size_continuous(
      range = c(2, 8)
    ) +
    ggplot2::labs(
      title =
        "Predicted ligand consequences of candidate TF normalization",
      subtitle =
        "Inputs for the subsequent macrophage-to-vascular communication stage",
      x = "Candidate macrophage ligand",
      y = NULL,
      fill =
        "Predicted mean\ndelta-z in HFpEF",
      size = "|Predicted delta-z|"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 55,
        hjust = 1
      )
    )

  save_plot_bundle(
    p_ligand,
    "Fig5E_candidate_TF_ligand_effects",
    12,
    6.5
  )
}

## Supplementary: perturbation strength curves.
strength_plot_dt <- program_effects[
  tf_symbol %in% candidate_tfs &
    perturbation_mode ==
      PRIMARY_MODE &
    perturbation_method ==
      PRIMARY_METHOD &
    program_category ==
      "Stage2_drug_opposed_net" &
    primary_program == TRUE &
    observed_gap_eligible == TRUE,
  .(
    median_gap_reduction =
      safe_median(
        absolute_gap_reduction,
        0
      )
  ),
  by = .(
    tf_symbol,
    perturbation_strength
  )
]

write_csv_safe(
  strength_plot_dt,
  file.path(
    DIRS$source,
    "FigS5A_perturbation_strength_source.csv"
  )
)

p_strength <- ggplot2::ggplot(
  strength_plot_dt,
  ggplot2::aes(
    x = perturbation_strength,
    y = median_gap_reduction,
    group = tf_symbol,
    label = tf_symbol
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  ggplot2::geom_line(
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    size = 2.5
  ) +
  ggrepel::geom_text_repel(
    data = strength_plot_dt[
      perturbation_strength ==
        max(PERTURBATION_STRENGTHS)
    ],
    size = 3,
    max.overlaps = Inf
  ) +
  ggplot2::scale_x_continuous(
    breaks = PERTURBATION_STRENGTHS
  ) +
  ggplot2::labs(
    title =
      "Perturbation-strength sensitivity",
    x = "Fraction of HFpEF TF activity normalized toward Control",
    y =
      "Median primary-program gap reduction"
  ) +
  ggplot2::theme_bw(base_size = 10)

save_plot_bundle(
  p_strength,
  "FigS5A_perturbation_strength_sensitivity",
  9,
  6.5
)

## Supplementary: candidate versus matched controls.
control_plot_dt <- primary_results[
  program_category ==
    "Stage2_drug_opposed_net" &
    primary_program == TRUE &
    observed_gap_eligible == TRUE,
  .(
    median_gap_reduction =
      safe_median(
        absolute_gap_reduction,
        0
      )
  ),
  by = .(
    tf_symbol,
    analysis_role
  )
]

control_plot_dt[, tf_symbol := factor(
  tf_symbol,
  levels = c(candidate_order, control_tfs)
)]

write_csv_safe(
  control_plot_dt,
  file.path(
    DIRS$source,
    "FigS5B_candidate_vs_control_source.csv"
  )
)

p_control <- ggplot2::ggplot(
  control_plot_dt,
  ggplot2::aes(
    x = tf_symbol,
    y = median_gap_reduction,
    fill = analysis_role
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  ggplot2::geom_col(
    width = 0.72
  ) +
  ggplot2::labs(
    title =
      "Biological candidates versus matched low-priority TF controls",
    x = NULL,
    y =
      "Median primary-program gap reduction",
    fill = "Analysis role"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot_bundle(
  p_control,
  "FigS5B_candidate_vs_matched_control_TFs",
  10,
  6
)

############################################################
## 15. Workbook, methods, and parameters
############################################################

workbook_path <- file.path(
  DIRS$tables,
  "18_stage5_multiTF_virtual_perturbation_key_results.xlsx"
)

xlsx_safe_df <- function(x) {
  y <- as.data.frame(
    x,
    stringsAsFactors = FALSE
  )

  factor_cols <- vapply(
    y,
    is.factor,
    logical(1)
  )
  if (any(factor_cols)) {
    y[factor_cols] <- lapply(
      y[factor_cols],
      as.character
    )
  }

  list_cols <- vapply(
    y,
    is.list,
    logical(1)
  )
  if (any(list_cols)) {
    y[list_cols] <- lapply(
      y[list_cols],
      function(z) {
        vapply(
          z,
          function(v) paste(v, collapse = ";"),
          character(1)
        )
      }
    )
  }

  y
}

workbook_sheets <- list(
  TF_manifest =
    xlsx_safe_df(analysis_tf_manifest),
  Program_definitions =
    xlsx_safe_df(program_summary),
  Observed_programs =
    xlsx_safe_df(program_observed_summary),
  Candidate_summary =
    xlsx_safe_df(candidate_summary),
  Ranking_sensitivity =
    xlsx_safe_df(ranking_sensitivity),
  Ranking_stability =
    xlsx_safe_df(ranking_stability),
  Method_concordance =
    xlsx_safe_df(method_concordance),
  Mode_comparison =
    xlsx_safe_df(mode_comparison),
  Matched_controls =
    xlsx_safe_df(control_primary_summary),
  Ligand_effects =
    xlsx_safe_df(ligand_effects),
  Top_gene_effects =
    xlsx_safe_df(top_gene_effects),
  State_effects =
    xlsx_safe_df(
      state_program_effects[
        perturbation_strength ==
          PRIMARY_STRENGTH &
          primary_program == TRUE
      ]
    )
)

writexl::write_xlsx(
  workbook_sheets,
  path = workbook_path
)

workbook_zip_list <- tryCatch(
  utils::unzip(
    workbook_path,
    list = TRUE
  ),
  error = function(e) NULL
)

workbook_integrity_ok <- (
  file.exists(workbook_path) &&
    file.info(workbook_path)$size > 10000 &&
    !is.null(workbook_zip_list) &&
    all(
      c(
        "[Content_Types].xml",
        "xl/workbook.xml",
        "xl/styles.xml"
      ) %in% workbook_zip_list$Name
    )
)

if (!workbook_integrity_ok) {
  stop(
    "The Stage 5 results workbook failed XLSX package-integrity checks."
  )
}

parameter_table <- data.table::data.table(
  parameter = c(
    "Random seed",
    "Biological candidate TFs",
    "Matched control TF count",
    "Primary perturbation mode",
    "Secondary perturbation mode",
    "Perturbation strengths",
    "Primary perturbation method",
    "Sensitivity perturbation method",
    "Maximum absolute gene shift",
    "Minimum observed program gap",
    "Primary signature size",
    "Minimum TF targets",
    "Minimum state profiles per condition",
    "Inferential unit",
    "Ranking strategy",
    "Nfkb1 forced"
  ),
  value = c(
    "20260714",
    paste(candidate_tfs, collapse = "; "),
    as.character(length(control_tfs)),
    PRIMARY_MODE,
    "activity_attenuation",
    paste(PERTURBATION_STRENGTHS, collapse = ", "),
    PRIMARY_METHOD,
    "equal_signed_targets",
    as.character(MAX_ABS_GENE_SHIFT_SD),
    as.character(MIN_ABS_OBSERVED_PROGRAM_GAP),
    as.character(PRIMARY_SIGNATURE_SIZE),
    as.character(MIN_TARGETS_PER_TF),
    as.character(MIN_STATE_PROFILES_PER_CONDITION),
    "Biological sample",
    "Mean and median rank aggregation across transparent evidence dimensions",
    "FALSE"
  ),
  rationale = c(
    "Reproducibility",
    "Prespecified from Stage 4 multi-method evidence",
    "Estimate specificity relative to low-priority TFs",
    "Moves each HFpEF TF activity toward the Control reference, regardless of whether it is disease-up or disease-down",
    "Tests whether uniform attenuation would improve or worsen the disease-associated state",
    "Assess monotonicity and strength sensitivity",
    "Minimum-norm expression adjustment consistent with the Stage 4 signed regulon score",
    "Alternative equal-signed target formulation",
    "Prevent implausibly large predicted target shifts",
    "Avoid unstable recovery fractions for near-zero baseline gaps",
    "Primary Stage 2/3 cross-stage program size",
    "Avoid unstable small regulons",
    "Require minimally represented state-specific comparisons",
    "Avoid cell-level pseudoreplication",
    "Avoid dependence on one opaque weighted composite score",
    "No candidate was manually promoted"
  )
)

write_csv_safe(
  parameter_table,
  file.path(
    DIRS$methods,
    "stage5_parameters_and_rationale.csv"
  )
)

methods_text <- c(
  "HFpEF Stage 5 FIXED v2: comparative multi-candidate in-silico TF perturbation",
  "",
  "Input boundary:",
  "- Stage 5 loaded the completed Stage 4 DoRothEA-based TF-target network, sample-level TF activity matrix, macrophage pseudobulk expression matrices, program manifest, and scored macrophage object.",
  "- Stage 5 did not repeat raw-data processing, quality control, doublet detection, clustering, cell annotation, pseudobulk differential expression, or Stage 4 TF activity inference.",
  "",
  "Candidate design:",
  paste0(
    "- Prespecified biological candidates: ",
    paste(candidate_tfs, collapse = ", "),
    "."
  ),
  paste0(
    "- Matched low-priority TF controls: ",
    paste(control_tfs, collapse = ", "),
    "."
  ),
  "- Controls were selected from the lower Stage 4 priority distribution while matching regulon size and TF activity/expression effect magnitudes.",
  "- Nfkb1 was not forced into any output rank.",
  "",
  "Perturbation formulation:",
  "- The primary perturbation was disease normalization: each HFpEF sample's candidate TF activity was shifted toward the mean Control activity.",
  "- This formulation can attenuate a disease-increased TF or restore a disease-decreased TF.",
  "- Uniform activity attenuation toward zero was analysed separately and was not assumed to be therapeutic.",
  "- Four strengths were tested: 25%, 50%, 75%, and 100%.",
  "- The weighted minimum-norm formulation distributed the required TF-activity change over signed targets in proportion to prior edge weight while minimizing squared target shifts.",
  "- The equal-signed formulation assigned the same signed shift to each target and served as an alternative perturbation sensitivity formulation.",
  "- Target shifts were capped at 2.5 standardized-expression units.",
  "",
  "Outcome programs:",
  "- Stage 2 Ccr2-positive, Ccr2-negative, and cross-subset drug-opposed programs were evaluated at Top50, Top100, Top150, and Top200 sizes.",
  "- Stage 3-supported subsets were evaluated separately.",
  "- Functional programs included inflammatory, NF-kB/TNF, inflammasome, interferon, antigen-presentation, remodeling, resident, monocyte-like, myeloid-identity, lipid/cholesterol, oxidative-stress, and cycling programs.",
  "",
  "Recovery metrics:",
  "- Program recovery was defined as reduction in the absolute HFpEF-Control program-score gap.",
  "- Program-specific recovery fractions were reported only when the observed absolute gap was at least 0.10 standardized units.",
  "- Biological-sample stability was the fraction of HFpEF samples moved closer to the Control mean.",
  "- Macrophage-state analyses required at least two Control and two HFpEF sample-state profiles.",
  "",
  "Priority strategy:",
  "- Candidate ranking used mean and median rank aggregation across primary program recovery, signature-size robustness, biological-sample stability, macrophage-state robustness, perturbation-method concordance, inflammatory-program recovery, specificity, matched-control margin, and Stage 4 support.",
  "- Five alternative ranking scenarios quantified dependence on evidence selection.",
  "- The ranking is a prioritization tool, not a statistical significance test.",
  "",
  "Claim boundary:",
  "- Results represent network-constrained predicted consequences of TF activity adjustment.",
  "- They are not experimental knockout, knockdown, or pharmacological inhibition results.",
  "- Stage 5 does not prove that dapagliflozin directly acts on any candidate TF.",
  "- Predicted ligand changes are inputs for a separate communication analysis and do not establish ligand-receptor signaling."
)

writeLines(
  methods_text,
  file.path(
    DIRS$methods,
    "stage5_methods_and_claim_boundaries.txt"
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
## 16. Completion checks and run status
############################################################

warnings_dt <- if (
  length(warning_records) > 0L
) {
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

write_csv_safe(
  warnings_dt,
  file.path(
    DIRS$tables,
    "19_stage5_warnings_and_nonfatal_issues.csv"
  )
)

expected_sample_runs <- length(analysis_tfs) *
  length(PERTURBATION_MODES) *
  length(PERTURBATION_STRENGTHS) *
  length(PERTURBATION_METHODS)

observed_sample_runs <- data.table::uniqueN(
  paste(
    program_effects$tf_symbol,
    program_effects$perturbation_mode,
    program_effects$perturbation_strength,
    program_effects$perturbation_method,
    sep = "__"
  )
)

scientific_checks <- data.table::data.table(
  check = c(
    "Stage 4 completed status",
    "Stage 4 failed checks",
    "Biological samples",
    "Control samples",
    "HFpEF samples",
    "Macrophage cells",
    "Prespecified candidate TFs",
    "Matched control TFs",
    "Candidate TFs in network",
    "Candidate TFs in activity matrix",
    "Evaluable programs",
    "Eligible primary Stage 2 programs",
    "Candidate TFs in state activity matrix",
    "Candidate/control regulon target counts match Stage 4",
    "Expected sample-level perturbation runs",
    "Observed sample-level perturbation runs",
    "Two perturbation formulations",
    "Two perturbation modes",
    "Four perturbation strengths",
    "Candidate summary rows",
    "Ranking sensitivity scenarios",
    "Scientific checkpoint",
    "Workbook"
  ),
  observed = c(
    as.integer(
      stage4_status$overall_status[1L] %in%
        allowed_stage4_status
    ),
    sum(stage4_checks$status != "PASS"),
    data.table::uniqueN(
      sample_meta$sample_accession
    ),
    sum(sample_meta$condition == "Control"),
    sum(sample_meta$condition == "HFpEF"),
    ncol(macrophage),
    length(candidate_tfs),
    length(control_tfs),
    sum(candidate_tfs %in%
      unique(network_dt$source_symbol)),
    sum(candidate_tfs %in%
      rownames(weighted_sample_activity)),
    length(program_definitions),
    eligible_primary_stage2_programs,
    sum(candidate_tfs %in%
      rownames(state_activity)),
    sum(
      analysis_integrity$target_count_match == TRUE
    ),
    expected_sample_runs,
    observed_sample_runs,
    data.table::uniqueN(
      program_effects$perturbation_method
    ),
    data.table::uniqueN(
      program_effects$perturbation_mode
    ),
    data.table::uniqueN(
      program_effects$perturbation_strength
    ),
    nrow(candidate_summary),
    data.table::uniqueN(
      ranking_sensitivity$scenario
    ),
    as.integer(
      file.exists(
        file.path(
          DIRS$objects,
          "CHECKPOINT_stage5_scientific_results_pre_figures.rds"
        )
      )
    ),
    as.integer(file.exists(workbook_path))
  ),
  expected = c(
    1L,
    0L,
    6L,
    3L,
    3L,
    1822L,
    6L,
    N_MATCHED_CONTROL_TFS,
    6L,
    6L,
    10L,
    2L,
    6L,
    length(analysis_tfs),
    expected_sample_runs,
    expected_sample_runs,
    2L,
    2L,
    4L,
    6L,
    5L,
    1L,
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
    "equal",
    "equal",
    "at_least",
    "at_least",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal"
  )
)

scientific_checks[, status := data.table::fcase(
  comparison == "equal" &
    observed == expected,
  "PASS",
  comparison == "at_least" &
    observed >= expected,
  "PASS",
  default = "FAIL"
)]

write_csv_safe(
  scientific_checks,
  file.path(
    DIRS$tables,
    "20_stage5_scientific_completion_checks.csv"
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

overall_status <- if (
  all(scientific_checks$status == "PASS")
) {
  "COMPLETED_STAGE5_READY_FOR_REVIEW"
} else {
  "COMPLETED_STAGE5_REVIEW_REQUIRED"
}

nfkb1_rank <- candidate_summary[
  gene_key(tf_symbol) == "NFKB1",
  final_candidate_rank
]

run_status <- data.table::data.table(
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
  biological_samples =
    nrow(sample_meta),
  macrophage_cells =
    ncol(macrophage),
  candidate_TFs =
    paste(candidate_tfs, collapse = ";"),
  matched_control_TFs =
    paste(control_tfs, collapse = ";"),
  perturbation_methods =
    paste(PERTURBATION_METHODS, collapse = ";"),
  perturbation_modes =
    paste(PERTURBATION_MODES, collapse = ";"),
  perturbation_strengths =
    paste(PERTURBATION_STRENGTHS, collapse = ";"),
  sample_level_runs =
    observed_sample_runs,
  regulon_integrity_failures =
    sum(
      regulon_integrity_audit$
        target_count_match != TRUE
    ),
  state_level_rows =
    nrow(state_program_effects),
  top_candidate_TF =
    candidate_summary$tf_symbol[1L],
  Nfkb1_rank =
    ifelse(
      length(nfkb1_rank) == 1L,
      nfkb1_rank,
      NA_integer_
    ),
  Nfkb1_forced = FALSE,
  warnings = nrow(warnings_dt),
  script_copy_status =
    script_copy_status,
  scientific_checks_failed =
    sum(scientific_checks$status != "PASS"),
  overall_status =
    overall_status
)

write_csv_safe(
  run_status,
  file.path(
    DIRS$tables,
    "21_stage5_run_status.csv"
  )
)

readme <- c(
  "HFpEF Reanalysis Project - Stage 5 FIXED v2",
  "Comparative multi-candidate in-silico TF perturbation",
  "",
  paste0("Overall status: ", overall_status),
  paste0(
    "Biological candidates: ",
    paste(candidate_tfs, collapse = ", ")
  ),
  paste0(
    "Matched controls: ",
    paste(control_tfs, collapse = ", ")
  ),
  paste0(
    "Top candidate by rank aggregation: ",
    candidate_summary$tf_symbol[1L]
  ),
  paste0(
    "Nfkb1 rank: ",
    ifelse(
      length(nfkb1_rank) == 1L,
      nfkb1_rank,
      "not ranked"
    )
  ),
  "",
  "Primary interpretation:",
  "- v2 verified that every perturbation used only the requested TF-specific regulon.",
  "- Disease-normalizing activity adjustment is the principal analysis.",
  "- Uniform activity attenuation is a separate directionality sensitivity analysis.",
  "- Two perturbation formulations and four strengths were compared.",
  "- Candidate TFs were compared with matched low-priority controls.",
  "- Ranking used transparent rank aggregation and five sensitivity scenarios.",
  "- Predicted ligand changes are prepared for the next communication stage.",
  "",
  "Causal boundary:",
  "- These are network-predicted perturbation consequences.",
  "- They are not experimental knockdown or knockout results.",
  "- Nfkb1 was not forced into the ranking.",
  "",
  "Upload the Stage 5 CHECK package before starting Stage 6."
)

writeLines(
  readme,
  file.path(OUT_DIR, "README_stage5.txt"),
  useBytes = TRUE
)

############################################################
## 17. CHECK package
############################################################

review_files <- c(
  LOG_FILE,
  file.path(
    DIRS$tables,
    "00_stage5_regulon_integrity_audit.csv"
  ),
  file.path(
    DIRS$tables,
    "01_stage5_candidate_TF_resolution.csv"
  ),
  file.path(
    DIRS$tables,
    "02_stage5_candidate_and_matched_control_TFs.csv"
  ),
  file.path(
    DIRS$tables,
    "03_stage5_program_definition_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "04_stage5_observed_program_scores.csv"
  ),
  file.path(
    DIRS$tables,
    "08_stage5_primary_top_predicted_genes_per_TF.csv"
  ),
  file.path(
    DIRS$tables,
    "09_stage5_candidate_ligand_changes_for_stage6.csv"
  ),
  file.path(
    DIRS$tables,
    "11_stage5_perturbation_method_concordance.csv"
  ),
  file.path(
    DIRS$tables,
    "12_stage5_matched_control_TF_reference_results.csv"
  ),
  file.path(
    DIRS$tables,
    "13_stage5_candidate_TF_rank_aggregation.csv"
  ),
  file.path(
    DIRS$tables,
    "14_stage5_candidate_ranking_sensitivity_scenarios.csv"
  ),
  file.path(
    DIRS$tables,
    "15_stage5_candidate_ranking_stability_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "16_stage5_normalization_vs_attenuation_results.csv"
  ),
  file.path(
    DIRS$tables,
    "17_stage5_normalization_vs_attenuation_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "18_stage5_multiTF_virtual_perturbation_key_results.xlsx"
  ),
  file.path(
    DIRS$tables,
    "19_stage5_warnings_and_nonfatal_issues.csv"
  ),
  file.path(
    DIRS$tables,
    "20_stage5_scientific_completion_checks.csv"
  ),
  file.path(
    DIRS$tables,
    "21_stage5_run_status.csv"
  ),
  file.path(
    DIRS$methods,
    "stage5_parameters_and_rationale.csv"
  ),
  file.path(
    DIRS$methods,
    "stage5_methods_and_claim_boundaries.txt"
  ),
  file.path(
    DIRS$methods,
    "sessionInfo.txt"
  ),
  file.path(
    OUT_DIR,
    "README_stage5.txt"
  ),
  list.files(
    DIRS$figures,
    pattern = "\\.png$",
    full.names = TRUE
  )
)

review_files <- unique(
  review_files[file.exists(review_files)]
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

check_manifest <- data.table::data.table(
  filename = basename(check_files),
  size_bytes =
    as.numeric(file.info(check_files)$size)
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

write_csv_safe(
  check_manifest,
  file.path(
    DIRS$check,
    "CHECK_package_file_manifest.csv"
  )
)

if (file.exists(CHECK_ZIP)) {
  unlink(CHECK_ZIP, force = TRUE)
}

zip::zipr(
  zipfile = CHECK_ZIP,
  files = list.files(
    DIRS$check,
    full.names = TRUE
  ),
  root = DIRS$check
)

log_msg("Stage 5 analysis finished.")
log_msg("Overall status: ", overall_status)
log_msg(
  "Top candidate TF: ",
  candidate_summary$tf_symbol[1L]
)
log_msg(
  "Nfkb1 rank: ",
  ifelse(
    length(nfkb1_rank) == 1L,
    nfkb1_rank,
    "not ranked"
  )
)
log_msg("CHECK package: ", CHECK_ZIP)

cat("\n============================================================\n")
cat("HFpEF Stage 5 multi-TF perturbation completed\n")
cat("Status: ", overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat(
  "Top candidate: ",
  candidate_summary$tf_symbol[1L],
  "\n",
  sep = ""
)
cat("Nfkb1 was not forced.\n")
cat("Upload the CHECK package before Stage 6.\n")
cat("============================================================\n")
