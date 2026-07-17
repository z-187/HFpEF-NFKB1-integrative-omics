############################################################
## HFpEF Reanalysis Project
## Stage 6 FINAL v3
## Candidate-TF-dependent macrophage-to-vascular/stromal
## communication analysis using the official NicheNet-v2
## mouse ligand-receptor network and ligand-target matrix
##
## Project:
##   <HFPEF_PROJECT_DIR>
##
## Required completed inputs:
##   Stage 3:
##     03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH
##   Stage 4:
##     04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1
##   Stage 5:
##     05_stage5_multiTF_virtual_perturbation_FIXED_v2
##   Stage 5B:
##     05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1
##
## FINAL v3 corrections and safeguards:
##   1) Retains the v2 correction that replaces calculated expressions
##      inside data.table::setorder() with explicit columns and setorderv().
##   2) Fixes the v2 figure-export failure caused by a dynamically
##      calculated ggplot height exceeding the 50-inch ggsave limit.
##   3) Restricts manuscript figures to prespecified top-ranked axes while
##      preserving complete source tables for all communication axes.
##   4) Caps every ggplot and heatmap export dimension to a safe range and
##      verifies that every exported PNG/PDF/TIFF file exists and is nonempty.
##   5) Removes the data.table pronoun ambiguity in ranking-sensitivity
##      calculations by explicit base data-frame column selection.
##   6) Adds figure-export auditing and completion checks.
##
## Candidate roles:
##   Nfkb1   primary inflammation/communication candidate
##   Rela    NF-kB family sensitivity candidate
##   Bhlhe40 program-recovery contrast candidate
##
## Sender:
##   Macrophage_Monocyte
##
## Receivers:
##   Endothelial
##   Fibroblast
##   Pericyte
##   Smooth_muscle
##
## Cross-stage gate:
##   Stage 5 TF-sensitive ligand
##   + macrophage expression
##   + Stage 2 drug-opposed direction where available
##   + Stage 3 macrophage pseudobulk direction support
##   + receiver receptor expression
##   + NicheNet-v2 receiver target activity
##
## Resource policy:
##   - Uses only the fixed official NicheNet-v2 Zenodo record
##     DOI: 10.5281/zenodo.7074291.
##   - Searches for verified local copies first.
##   - Downloads only when the files are absent.
##   - Verifies the published MD5 hashes before analysis.
##   - Does not use decoupleR, OmnipathR, CellChat, or an
##     unversioned live API.
##
## Interpretation boundary:
##   - This is ligand-receptor-target prioritization, not proof
##     of physical signaling or TF causality.
##   - Stage 5 ligand changes are predicted consequences of
##     TF-activity normalization, not experimental perturbations.
##   - Receiver differential expression is sample-level
##     pseudobulk evidence; cell-level expression is descriptive.
##   - Sample ligand-receptor correlations use six paired samples
##     and are descriptive only.
##
## Output:
##   <HFPEF_PROJECT_DIR>/
##   06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3
##
## CHECK:
##   <HFPEF_PROJECT_DIR>/
##   06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3_CHECK.zip
##
## Recommended run:
##   source(
##     "<HFPEF_PROJECT_DIR>/06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3.R",
##     encoding = "UTF-8"
##   )
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(warn = 1)
options(encoding = "UTF-8")
options(timeout = 7200)
options(future.globals.maxSize = 16 * 1024^3)

set.seed(20260714)

############################################################
## 0. Locked paths and settings
############################################################

PROJECT_DIR <- Sys.getenv("HFPEF_PROJECT_DIR", unset = "")
if (!nzchar(PROJECT_DIR)) {
  stop(
    "HFPEF_PROJECT_DIR is not set. Define it as the local project root ",
    "containing 0.GEO and the completed stage-output folders."
  )
}
PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)

STAGE3_DIR <- file.path(
  PROJECT_DIR,
  "03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH"
)

STAGE4_DIR <- file.path(
  PROJECT_DIR,
  "04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1"
)

STAGE5_DIR <- file.path(
  PROJECT_DIR,
  "05_stage5_multiTF_virtual_perturbation_FIXED_v2"
)

STAGE5B_DIR <- file.path(
  PROJECT_DIR,
  "05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1"
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

STAGE3_SAMPLE_META_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "01_locked_GSE236585_sample_metadata.csv"
)

STAGE3_SEURAT_FILE <- file.path(
  STAGE3_DIR,
  "02_objects",
  "GSE236585_stage3_annotated_projected_seurat.rds"
)

STAGE3_PSEUDOBULK_DE_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "22_major_celltype_pseudobulk_DE_edgeR_limma.csv.gz"
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

STAGE4_PROGRAM_FILE <- file.path(
  STAGE4_DIR,
  "01_tables",
  "02_stage4_program_gene_manifest.csv"
)

STAGE5_STATUS_FILE <- file.path(
  STAGE5_DIR,
  "01_tables",
  "21_stage5_run_status.csv"
)

STAGE5_CHECKS_FILE <- file.path(
  STAGE5_DIR,
  "01_tables",
  "20_stage5_scientific_completion_checks.csv"
)

STAGE5_ALL_TARGET_FILE <- file.path(
  STAGE5_DIR,
  "01_tables",
  "07_stage5_all_predicted_target_gene_changes.csv"
)

STAGE5_LIGAND_PANEL_FILE <- file.path(
  STAGE5_DIR,
  "01_tables",
  "09_stage5_candidate_ligand_changes_for_stage6.csv"
)

STAGE5B_STATUS_FILE <- file.path(
  STAGE5B_DIR,
  "01_tables",
  "17_stage5B_run_status.csv"
)

STAGE5B_CHECKS_FILE <- file.path(
  STAGE5B_DIR,
  "01_tables",
  "16_stage5B_scientific_completion_checks.csv"
)

STAGE5B_RANK_FILE <- file.path(
  STAGE5B_DIR,
  "01_tables",
  "13_stage5B_final_candidate_robustness_rank.csv"
)

STAGE_NAME <-
  "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3"

OUT_DIR <- file.path(
  PROJECT_DIR,
  STAGE_NAME
)

CHECK_ZIP <- file.path(
  PROJECT_DIR,
  paste0(
    STAGE_NAME,
    "_CHECK.zip"
  )
)

EXPECTED_SCRIPT_FILE <- file.path(
  PROJECT_DIR,
  "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3.R"
)

REPLACE_EXISTING_STAGE6 <- TRUE

RESOURCE_DIR <- file.path(
  PROJECT_DIR,
  "00_external_resources",
  "NicheNet_v2_7074291"
)

NICHENET_MATRIX_FILE <- file.path(
  RESOURCE_DIR,
  "ligand_target_matrix_nsga2r_final_mouse.rds"
)

NICHENET_LR_FILE <- file.path(
  RESOURCE_DIR,
  "lr_network_mouse_21122021.rds"
)

NICHENET_MATRIX_MD5 <-
  "ac80d846fe0bfc4879a5b52ca85ffeb9"

NICHENET_LR_MD5 <-
  "cf33ee8b6bf84bdf2d11cab9c8f94b9e"

NICHENET_RECORD_ID <- "7074291"

NICHENET_FILES <- list(
  ligand_target_matrix = list(
    filename =
      "ligand_target_matrix_nsga2r_final_mouse.rds",
    destination =
      NICHENET_MATRIX_FILE,
    md5 =
      NICHENET_MATRIX_MD5,
    minimum_size =
      150 * 1024^2
  ),
  ligand_receptor_network = list(
    filename =
      "lr_network_mouse_21122021.rds",
    destination =
      NICHENET_LR_FILE,
    md5 =
      NICHENET_LR_MD5,
    minimum_size =
      20000
  )
)

NICHENET_SEARCH_DIRS <- unique(c(
  RESOURCE_DIR,
  PROJECT_DIR,
  file.path(
    PROJECT_DIR,
    "00_resources"
  ),
  file.path(
    PROJECT_DIR,
    "NicheNet_resources"
  ),
  file.path(
    PROJECT_DIR,
    "nichenet_resources"
  ),
  file.path(PROJECT_DIR, "00_external_resources"),
  file.path(PROJECT_DIR, "00_resources")
))

CANDIDATE_TFS <- c(
  "Nfkb1",
  "Rela",
  "Bhlhe40"
)

CANDIDATE_ROLE_MAP <- c(
  Nfkb1 =
    "Primary_inflammation_communication_candidate",
  Rela =
    "NFkB_family_sensitivity_candidate",
  Bhlhe40 =
    "Program_recovery_contrast_candidate"
)

SENDER_CELL_TYPE <-
  "Macrophage_Monocyte"

RECEIVER_CELL_TYPES <- c(
  "Endothelial",
  "Fibroblast",
  "Pericyte",
  "Smooth_muscle"
)

PRIMARY_PERTURBATION_METHOD <-
  "weighted_minimum_norm"

PRIMARY_PERTURBATION_MODE <-
  "disease_normalization"

PRIMARY_PERTURBATION_STRENGTH <-
  1

PRIMARY_STAGE2_SIGNATURE_SIZE <-
  150L

MIN_SENDER_CELLS_PER_SAMPLE <-
  20L

MIN_RECEIVER_CELLS_PER_SAMPLE <-
  50L

MIN_EXPRESSED_CELL_FRACTION <-
  0.05

MIN_HFPEF_SAMPLES_EXPRESSED <-
  2L

RECEIVER_PRIMARY_FDR <-
  0.10

RECEIVER_MIN_ABS_LOGFC <-
  0.10

MIN_RECEIVER_TARGET_GENES <-
  20L

FALLBACK_RECEIVER_TARGET_GENES <-
  150L

MAX_RECEIVER_TARGET_GENES <-
  300L

TOP_TARGET_LINKS_PER_LIGAND <-
  50L

TOP_AXES_FOR_REPORT <-
  100L

TOP_AXES_FOR_FIGURE <-
  30L

TOP_TARGETS_FOR_HEATMAP <-
  30L

MAX_RECEPTOR_AXES_FOR_FIGURE <-
  40L

MAX_ACTIVITY_ROWS_FOR_HEATMAP <-
  50L

MAX_TARGET_ROWS_FOR_HEATMAP <-
  60L

MAX_GGPLOT_WIDTH_IN <-
  18

MAX_GGPLOT_HEIGHT_IN <-
  18

MAX_HEATMAP_WIDTH_IN <-
  18

MAX_HEATMAP_HEIGHT_IN <-
  18

RESOURCE_DOWNLOAD_ATTEMPTS <-
  4L

if (
  any(
    c(
      MAX_GGPLOT_WIDTH_IN,
      MAX_GGPLOT_HEIGHT_IN,
      MAX_HEATMAP_WIDTH_IN,
      MAX_HEATMAP_HEIGHT_IN
    ) >=
      50
  )
) {
  stop(
    "Configured figure dimensions must remain below the 50-inch ggsave limit."
  )
}

############################################################
## 1. Preflight, replacement, directories, and logging
############################################################

detect_script_file <- function() {
  candidates <- character()

  frames <- sys.frames()

  for (frame_index in rev(seq_along(frames))) {
    ofile <- tryCatch(
      frames[[frame_index]]$ofile,
      error = function(e) NULL
    )

    if (
      !is.null(ofile) &&
      length(ofile) == 1L &&
      nzchar(ofile)
    ) {
      candidates <- c(
        candidates,
        ofile
      )
    }
  }

  arguments <- commandArgs(
    trailingOnly = FALSE
  )

  file_argument <- grep(
    "^--file=",
    arguments,
    value = TRUE
  )

  if (length(file_argument) > 0L) {
    candidates <- c(
      candidates,
      sub(
        "^--file=",
        "",
        file_argument[1L]
      )
    )
  }

  candidates <- unique(candidates)
  candidates <- candidates[
    file.exists(candidates)
  ]

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
  STAGE3_STATUS_FILE,
  STAGE3_CHECKS_FILE,
  STAGE3_SAMPLE_META_FILE,
  STAGE3_SEURAT_FILE,
  STAGE3_PSEUDOBULK_DE_FILE,
  STAGE4_STATUS_FILE,
  STAGE4_CHECKS_FILE,
  STAGE4_PROGRAM_FILE,
  STAGE5_STATUS_FILE,
  STAGE5_CHECKS_FILE,
  STAGE5_ALL_TARGET_FILE,
  STAGE5_LIGAND_PANEL_FILE,
  STAGE5B_STATUS_FILE,
  STAGE5B_CHECKS_FILE,
  STAGE5B_RANK_FILE
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs) > 0L) {
  stop(
    "Required Stage 3/4/5/5B input path(s) are missing:\n",
    paste(
      missing_inputs,
      collapse = "\n"
    )
  )
}

stage3_status <- data.table::fread(
  STAGE3_STATUS_FILE,
  encoding = "UTF-8"
)

stage4_status <- data.table::fread(
  STAGE4_STATUS_FILE,
  encoding = "UTF-8"
)

stage5_status <- data.table::fread(
  STAGE5_STATUS_FILE,
  encoding = "UTF-8"
)

stage5b_status <- data.table::fread(
  STAGE5B_STATUS_FILE,
  encoding = "UTF-8"
)

stage3_checks <- data.table::fread(
  STAGE3_CHECKS_FILE,
  encoding = "UTF-8"
)

stage4_checks <- data.table::fread(
  STAGE4_CHECKS_FILE,
  encoding = "UTF-8"
)

stage5_checks <- data.table::fread(
  STAGE5_CHECKS_FILE,
  encoding = "UTF-8"
)

stage5b_checks <- data.table::fread(
  STAGE5B_CHECKS_FILE,
  encoding = "UTF-8"
)

expected_statuses <- c(
  Stage3 =
    "COMPLETED_STAGE3_READY_FOR_REVIEW",
  Stage4 =
    "COMPLETED_STAGE4_READY_FOR_REVIEW",
  Stage5 =
    "COMPLETED_STAGE5_READY_FOR_REVIEW",
  Stage5B =
    "COMPLETED_STAGE5B_OFFLINE_READY_FOR_REVIEW"
)

observed_statuses <- c(
  Stage3 =
    stage3_status$overall_status[1L],
  Stage4 =
    stage4_status$overall_status[1L],
  Stage5 =
    stage5_status$overall_status[1L],
  Stage5B =
    stage5b_status$overall_status[1L]
)

if (
  !identical(
    unname(observed_statuses),
    unname(expected_statuses)
  )
) {
  stop(
    "One or more required stages are not in the locked completed state:\n",
    paste(
      names(observed_statuses),
      observed_statuses,
      sep = "=",
      collapse = "\n"
    )
  )
}

for (
  check_table in list(
    stage3_checks,
    stage4_checks,
    stage5_checks,
    stage5b_checks
  )
) {
  if (
    !all(
      c(
        "check",
        "status"
      ) %in%
        names(check_table)
    ) ||
    any(
      check_table$status !=
        "PASS"
    )
  ) {
    stop(
      "At least one required upstream scientific completion check is not PASS."
    )
  }
}

replacement_audit <- data.table::data.table(
  path = c(
    OUT_DIR,
    CHECK_ZIP
  ),
  path_type = c(
    "stage6_output_directory",
    "stage6_check_zip"
  ),
  existed_before = FALSE,
  deletion_attempted = FALSE,
  deletion_succeeded = FALSE
)

if (REPLACE_EXISTING_STAGE6) {
  for (
    audit_index in seq_len(
      nrow(replacement_audit)
    )
  ) {
    target_path <-
      replacement_audit$path[
        audit_index
      ]

    existed <- (
      dir.exists(target_path) ||
        file.exists(target_path)
    )

    replacement_audit$
      existed_before[
        audit_index
      ] <- existed

    if (existed) {
      replacement_audit$
        deletion_attempted[
          audit_index
        ] <- TRUE

      unlink(
        target_path,
        recursive =
          dir.exists(target_path),
        force = TRUE
      )

      replacement_audit$
        deletion_succeeded[
          audit_index
        ] <- !(
          dir.exists(target_path) ||
            file.exists(target_path)
        )

      if (
        !replacement_audit$
          deletion_succeeded[
            audit_index
          ]
      ) {
        stop(
          "Failed to remove previous Stage 6 path: ",
          target_path
        )
      }
    } else {
      replacement_audit$
        deletion_succeeded[
          audit_index
        ] <- TRUE
    }
  }
} else if (
  dir.exists(OUT_DIR) ||
  file.exists(CHECK_ZIP)
) {
  stop(
    "Existing Stage 6 output detected while replacement is disabled."
  )
}

DIRS <- list(
  logs = file.path(
    OUT_DIR,
    "00_logs"
  ),
  tables = file.path(
    OUT_DIR,
    "01_tables"
  ),
  objects = file.path(
    OUT_DIR,
    "02_objects"
  ),
  figures = file.path(
    OUT_DIR,
    "03_figures"
  ),
  source = file.path(
    OUT_DIR,
    "04_source_data"
  ),
  methods = file.path(
    OUT_DIR,
    "05_methods"
  ),
  check = file.path(
    OUT_DIR,
    "06_review_check"
  )
)

dir.create(
  RESOURCE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

for (
  directory_i in c(
    OUT_DIR,
    unlist(
      DIRS,
      use.names = FALSE
    )
  )
) {
  dir.create(
    directory_i,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

START_TIME <- Sys.time()

LOG_FILE <- file.path(
  DIRS$logs,
  "stage6_TF_dependent_communication.log"
)

WARN_FILE <- file.path(
  DIRS$logs,
  "stage6_warnings.log"
)

data.table::fwrite(
  replacement_audit,
  file.path(
    DIRS$logs,
    "stage6_replacement_audit.csv"
  )
)

log_msg <- function(
  ...,
  level = "INFO"
) {
  text <- paste0(
    ...,
    collapse = ""
  )

  line <- sprintf(
    "[%s] [%s] %s",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    level,
    text
  )

  cat(
    line,
    "\n"
  )

  cat(
    line,
    "\n",
    file = LOG_FILE,
    append = TRUE
  )

  invisible(line)
}

warning_records <- list()

add_warning <- function(
  category,
  item,
  message
) {
  record <- data.table::data.table(
    timestamp = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    category = as.character(category),
    item = as.character(item),
    message = as.character(message)
  )

  warning_records[[length(warning_records) + 1L]] <<- record

  cat(
    sprintf(
      "[%s] [%s] %s: %s\n",
      record$timestamp,
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

  invisible(record)
}

log_msg(
  "Stage 6 analysis started."
)

log_msg(
  "Upstream statuses: ",
  paste(
    names(observed_statuses),
    observed_statuses,
    sep = "=",
    collapse = "; "
  )
)

log_msg(
  "Output: ",
  OUT_DIR
)

############################################################
## 2. Packages
############################################################

ensure_cran <- function(
  packages,
  required = TRUE
) {
  missing <- packages[
    !vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing) > 0L) {
    log_msg(
      "Installing missing CRAN package(s): ",
      paste(
        missing,
        collapse = ", "
      )
    )

    try(
      install.packages(
        missing,
        repos =
          "https://cloud.r-project.org",
        dependencies = TRUE
      ),
      silent = TRUE
    )
  }

  still_missing <- packages[
    !vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(still_missing) > 0L) {
    package_message <- paste(
      "Required CRAN package(s) unavailable:",
      paste(
        still_missing,
        collapse = ", "
      )
    )

    if (required) {
      stop(package_message)
    } else {
      add_warning(
        "PACKAGE",
        paste(
          still_missing,
          collapse = ";"
        ),
        package_message
      )
    }
  }

  invisible(
    setdiff(
      packages,
      still_missing
    )
  )
}

ensure_cran(
  c(
    "Seurat",
    "SeuratObject",
    "data.table",
    "Matrix",
    "ggplot2",
    "pheatmap",
    "writexl",
    "zip",
    "digest",
    "scales"
  ),
  required = TRUE
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(pheatmap)
})

ordering_self_test <- data.table::data.table(
  combined_p = c(
    0.01,
    0.01,
    0.02
  ),
  abs_mean_logFC = c(
    0.2,
    0.8,
    1.0
  )
)

data.table::setorderv(
  ordering_self_test,
  cols = c(
    "combined_p",
    "abs_mean_logFC"
  ),
  order = c(
    1L,
    -1L
  ),
  na.last = TRUE
)

if (
  !identical(
    ordering_self_test$
      abs_mean_logFC,
    c(
      0.8,
      0.2,
      1.0
    )
  )
) {
  stop(
    "Internal data.table ordering self-test failed."
  )
}

rm(ordering_self_test)

scenario_column_self_test <- data.frame(
  rank_a = c(
    1,
    2
  ),
  rank_b = c(
    3,
    4
  )
)

scenario_columns_self_test <- c(
  "rank_a",
  "rank_b",
  "rank_b"
)

scenario_matrix_self_test <- as.matrix(
  scenario_column_self_test[
    ,
    scenario_columns_self_test,
    drop = FALSE
  ]
)

if (
  !identical(
    dim(
      scenario_matrix_self_test
    ),
    c(
      2L,
      3L
    )
  )
) {
  stop(
    "Internal scenario-column duplication self-test failed."
  )
}

rm(
  scenario_column_self_test,
  scenario_columns_self_test,
  scenario_matrix_self_test
)

############################################################
## 3. General utilities
############################################################

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""

  x <- gsub(
    "\\r|\\n",
    " ",
    x
  )

  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )

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

read_table_auto <- function(path) {
  binary_connection <- file(
    path,
    open = "rb"
  )

  magic <- readBin(
    binary_connection,
    what = "raw",
    n = 2L
  )

  close(binary_connection)

  is_gzip <- (
    length(magic) == 2L &&
      identical(
        as.integer(magic),
        c(
          31L,
          139L
        )
      )
  )

  if (is_gzip) {
    gzip_connection <- gzfile(
      path,
      open = "rt"
    )

    on.exit(
      close(gzip_connection),
      add = TRUE
    )

    output <- data.table::as.data.table(
      utils::read.csv(
        gzip_connection,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )

    close(gzip_connection)

    on.exit(
      NULL,
      add = FALSE
    )

    output
  } else {
    data.table::fread(
      path,
      encoding = "UTF-8"
    )
  }
}

write_csv_safe <- function(
  table_object,
  path,
  compress = FALSE
) {
  if (
    is.null(table_object) ||
    ncol(table_object) == 0L
  ) {
    data.table::fwrite(
      data.table::data.table(
        note =
          "No records generated."
      ),
      path
    )
  } else {
    data.table::fwrite(
      table_object,
      path,
      compress = if (compress) {
        "gzip"
      } else {
        "none"
      }
    )
  }
}

safe_mean <- function(
  x,
  default = NA_real_
) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(default)
  }

  mean(x)
}

safe_median <- function(
  x,
  default = NA_real_
) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(default)
  }

  stats::median(x)
}

safe_max <- function(
  x,
  default = NA_real_
) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(default)
  }

  max(x)
}

safe_quantile <- function(
  x,
  probability,
  default = NA_real_
) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(default)
  }

  as.numeric(
    stats::quantile(
      x,
      probs = probability,
      names = FALSE,
      na.rm = TRUE,
      type = 7
    )
  )
}

safe_spearman <- function(
  x,
  y
) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  keep <- (
    is.finite(x) &
      is.finite(y)
  )

  if (sum(keep) < 4L) {
    return(NA_real_)
  }

  suppressWarnings(
    stats::cor(
      x[keep],
      y[keep],
      method = "spearman"
    )
  )
}

safe_pearson <- function(
  x,
  y
) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  keep <- (
    is.finite(x) &
      is.finite(y)
  )

  if (sum(keep) < 4L) {
    return(NA_real_)
  }

  suppressWarnings(
    stats::cor(
      x[keep],
      y[keep],
      method = "pearson"
    )
  )
}

rank_metric <- function(
  x,
  higher_is_better = TRUE
) {
  x <- as.numeric(x)
  finite <- is.finite(x)

  output <- rep(
    sum(finite) + 1,
    length(x)
  )

  if (any(finite)) {
    output[finite] <- rank(
      if (higher_is_better) {
        -x[finite]
      } else {
        x[finite]
      },
      ties.method = "average"
    )
  }

  output
}

make_feature_map <- function(features) {
  feature_table <- data.table::data.table(
    feature = as.character(features),
    feature_key =
      gene_key(features)
  )

  data.table::setorder(
    feature_table,
    feature_key,
    feature
  )

  feature_table[
    ,
    .SD[1L],
    by = feature_key
  ]
}

map_symbols_to_features <- function(
  symbols,
  feature_map
) {
  query <- data.table::data.table(
    symbol = as.character(symbols),
    feature_key =
      gene_key(symbols)
  )

  merge(
    query,
    feature_map,
    by = "feature_key",
    all.x = TRUE,
    sort = FALSE
  )
}

get_assay_matrix <- function(
  object,
  assay = "RNA",
  layer = "counts",
  features = NULL,
  cells = NULL
) {
  Seurat::DefaultAssay(object) <-
    assay

  matrix_object <- tryCatch(
    SeuratObject::LayerData(
      object,
      assay = assay,
      layer = layer,
      features = features,
      cells = cells
    ),
    error = function(error_layer) {
      tryCatch(
        Seurat::GetAssayData(
          object,
          assay = assay,
          layer = layer
        ),
        error = function(error_get_layer) {
          Seurat::GetAssayData(
            object,
            assay = assay,
            slot = layer
          )
        }
      )
    }
  )

  if (!is.null(features)) {
    retained_features <- intersect(
      features,
      rownames(matrix_object)
    )

    matrix_object <- matrix_object[
      retained_features,
      ,
      drop = FALSE
    ]
  }

  if (!is.null(cells)) {
    retained_cells <- intersect(
      cells,
      colnames(matrix_object)
    )

    matrix_object <- matrix_object[
      ,
      retained_cells,
      drop = FALSE
    ]
  }

  matrix_object
}

join_layers_compat <- function(
  object,
  assay = "RNA"
) {
  Seurat::DefaultAssay(object) <-
    assay

  output <- tryCatch(
    {
      if (
        "JoinLayers" %in%
          getNamespaceExports(
            "SeuratObject"
          )
      ) {
        SeuratObject::JoinLayers(
          object,
          assay = assay
        )
      } else if (
        "JoinLayers" %in%
          getNamespaceExports(
            "Seurat"
          )
      ) {
        Seurat::JoinLayers(
          object,
          assay = assay
        )
      } else {
        object
      }
    },
    error = function(error_join) {
      add_warning(
        "SEURAT_LAYER",
        assay,
        paste0(
          "JoinLayers was not applied: ",
          conditionMessage(
            error_join
          )
        )
      )

      object
    }
  )

  Seurat::DefaultAssay(output) <-
    assay

  output
}

simple_auc <- function(
  labels,
  scores
) {
  labels <- as.integer(labels)
  scores <- as.numeric(scores)

  keep <- (
    !is.na(labels) &
      is.finite(scores)
  )

  labels <- labels[keep]
  scores <- scores[keep]

  positive_n <- sum(
    labels == 1L
  )

  negative_n <- sum(
    labels == 0L
  )

  if (
    positive_n == 0L ||
    negative_n == 0L
  ) {
    return(NA_real_)
  }

  score_ranks <- rank(
    scores,
    ties.method = "average"
  )

  (
    sum(
      score_ranks[
        labels == 1L
      ]
    ) -
      positive_n *
        (
          positive_n + 1
        ) /
        2
  ) /
    (
      positive_n *
        negative_n
    )
}

simple_aupr <- function(
  labels,
  scores
) {
  labels <- as.integer(labels)
  scores <- as.numeric(scores)

  keep <- (
    !is.na(labels) &
      is.finite(scores)
  )

  labels <- labels[keep]
  scores <- scores[keep]

  positive_n <- sum(
    labels == 1L
  )

  if (positive_n == 0L) {
    return(NA_real_)
  }

  ordering <- order(
    scores,
    decreasing = TRUE
  )

  labels <- labels[ordering]

  true_positive <- cumsum(
    labels == 1L
  )

  false_positive <- cumsum(
    labels == 0L
  )

  recall <- c(
    0,
    true_positive /
      positive_n
  )

  precision <- c(
    1,
    true_positive /
      pmax(
        true_positive +
          false_positive,
        1
      )
  )

  sum(
    diff(recall) *
      precision[-1L],
    na.rm = TRUE
  )
}

split_lr_entity <- function(entity) {
  entity <- normalize_text(entity)

  entity <- gsub(
    "[()]",
    "",
    entity
  )

  components <- unlist(
    strsplit(
      entity,
      "[_&+;:|]",
      perl = TRUE
    ),
    use.names = FALSE
  )

  components <- trimws(components)

  unique(
    components[
      nzchar(components)
    ]
  )
}

figure_export_records <- list()

sanitize_figure_dimension <- function(
  value,
  minimum,
  maximum,
  dimension_name,
  stem
) {
  value <- as.numeric(value)

  if (
    length(value) != 1L ||
    !is.finite(value) ||
    value <= 0
  ) {
    stop(
      "Invalid ",
      dimension_name,
      " for figure ",
      stem,
      ": ",
      paste(value, collapse = ";")
    )
  }

  effective_value <- min(
    max(
      value,
      minimum
    ),
    maximum
  )

  if (
    abs(
      effective_value -
        value
    ) >
      1e-12
  ) {
    add_warning(
      "FIGURE_DIMENSION",
      stem,
      paste0(
        dimension_name,
        " was capped from ",
        round(value, 3),
        " to ",
        round(
          effective_value,
          3
        ),
        " inches."
      )
    )
  }

  effective_value
}

record_figure_export <- function(
  stem,
  plot_type,
  requested_width,
  requested_height,
  effective_width,
  effective_height,
  paths
) {
  paths_exist <- file.exists(paths)
  sizes <- rep(
    NA_real_,
    length(paths)
  )

  sizes[paths_exist] <- as.numeric(
    file.info(
      paths[paths_exist]
    )$size
  )

  files_valid <- (
    length(paths) == 3L &&
      all(paths_exist) &&
      all(
        is.finite(sizes) &
          sizes > 0
      )
  )

  figure_export_records[[length(figure_export_records) + 1L]] <-
    data.table::data.table(
    stem = stem,
    plot_type = plot_type,
    requested_width_in =
      as.numeric(
        requested_width
      ),
    requested_height_in =
      as.numeric(
        requested_height
      ),
    effective_width_in =
      as.numeric(
        effective_width
      ),
    effective_height_in =
      as.numeric(
        effective_height
      ),
    png_exists =
      paths_exist["png"],
    pdf_exists =
      paths_exist["pdf"],
    tiff_exists =
      paths_exist["tiff"],
    png_size_bytes =
      sizes["png"],
    pdf_size_bytes =
      sizes["pdf"],
    tiff_size_bytes =
      sizes["tiff"],
    files_valid =
      files_valid
  )

  if (!files_valid) {
    stop(
      "Figure export validation failed for ",
      stem,
      "."
    )
  }

  invisible(TRUE)
}

save_plot_bundle <- function(
  plot_object,
  stem,
  width,
  height
) {
  if (
    !inherits(
      plot_object,
      "ggplot"
    )
  ) {
    stop(
      "save_plot_bundle received a non-ggplot object for ",
      stem,
      "."
    )
  }

  requested_width <- width
  requested_height <- height

  width <- sanitize_figure_dimension(
    value = width,
    minimum = 3,
    maximum =
      MAX_GGPLOT_WIDTH_IN,
    dimension_name = "width",
    stem = stem
  )

  height <- sanitize_figure_dimension(
    value = height,
    minimum = 3,
    maximum =
      MAX_GGPLOT_HEIGHT_IN,
    dimension_name = "height",
    stem = stem
  )

  paths <- c(
    png = file.path(
      DIRS$figures,
      paste0(
        stem,
        ".png"
      )
    ),
    pdf = file.path(
      DIRS$figures,
      paste0(
        stem,
        ".pdf"
      )
    ),
    tiff = file.path(
      DIRS$figures,
      paste0(
        stem,
        ".tiff"
      )
    )
  )

  ggplot2::ggsave(
    filename = paths["png"],
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = paths["pdf"],
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    limitsize = TRUE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = paths["tiff"],
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    limitsize = TRUE,
    bg = "white"
  )

  record_figure_export(
    stem = stem,
    plot_type = "ggplot",
    requested_width =
      requested_width,
    requested_height =
      requested_height,
    effective_width = width,
    effective_height = height,
    paths = paths
  )

  invisible(paths)
}

save_heatmap_bundle <- function(
  matrix_object,
  stem,
  width,
  height,
  title = NULL
) {
  matrix_object <- as.matrix(
    matrix_object
  )

  if (
    nrow(matrix_object) < 1L ||
    ncol(matrix_object) < 1L
  ) {
    stop(
      "save_heatmap_bundle received an empty matrix for ",
      stem,
      "."
    )
  }

  if (!is.numeric(matrix_object)) {
    stop(
      "save_heatmap_bundle requires a numeric matrix for ",
      stem,
      "."
    )
  }

  if (
    any(
      !is.finite(
        matrix_object
      )
    )
  ) {
    add_warning(
      "HEATMAP_NONFINITE",
      stem,
      "Non-finite values were replaced with zero before plotting."
    )

    matrix_object[
      !is.finite(
        matrix_object
      )
    ] <- 0
  }

  requested_width <- width
  requested_height <- height

  width <- sanitize_figure_dimension(
    value = width,
    minimum = 3,
    maximum =
      MAX_HEATMAP_WIDTH_IN,
    dimension_name = "width",
    stem = stem
  )

  height <- sanitize_figure_dimension(
    value = height,
    minimum = 3,
    maximum =
      MAX_HEATMAP_HEIGHT_IN,
    dimension_name = "height",
    stem = stem
  )

  paths <- c(
    png = file.path(
      DIRS$figures,
      paste0(
        stem,
        ".png"
      )
    ),
    pdf = file.path(
      DIRS$figures,
      paste0(
        stem,
        ".pdf"
      )
    ),
    tiff = file.path(
      DIRS$figures,
      paste0(
        stem,
        ".tiff"
      )
    )
  )

  draw_heatmap <- function() {
    pheatmap::pheatmap(
      matrix_object,
      cluster_rows =
        nrow(matrix_object) > 1L,
      cluster_cols =
        ncol(matrix_object) > 1L,
      border_color = NA,
      main = title
    )
  }

  grDevices::png(
    filename = paths["png"],
    width = round(
      width * 300
    ),
    height = round(
      height * 300
    ),
    res = 300
  )

  tryCatch(
    draw_heatmap(),
    finally = {
      grDevices::dev.off()
    }
  )

  grDevices::pdf(
    file = paths["pdf"],
    width = width,
    height = height
  )

  tryCatch(
    draw_heatmap(),
    finally = {
      grDevices::dev.off()
    }
  )

  grDevices::tiff(
    filename = paths["tiff"],
    width = width,
    height = height,
    units = "in",
    res = 600,
    compression = "lzw"
  )

  tryCatch(
    draw_heatmap(),
    finally = {
      grDevices::dev.off()
    }
  )

  record_figure_export(
    stem = stem,
    plot_type = "heatmap",
    requested_width =
      requested_width,
    requested_height =
      requested_height,
    effective_width = width,
    effective_height = height,
    paths = paths
  )

  invisible(paths)
}


############################################################
## 4. Verified NicheNet-v2 resources
############################################################

file_md5 <- function(path) {
  unname(
    tools::md5sum(path)
  )
}

resource_is_valid <- function(
  path,
  expected_md5,
  minimum_size
) {
  if (!file.exists(path)) {
    return(FALSE)
  }

  file_size <- as.numeric(
    file.info(path)$size
  )

  if (
    !is.finite(file_size) ||
    file_size <
      minimum_size
  ) {
    return(FALSE)
  }

  observed_md5 <- tryCatch(
    file_md5(path),
    error = function(e) NA_character_
  )

  isTRUE(
    identical(
      tolower(observed_md5),
      tolower(expected_md5)
    )
  )
}

find_verified_local_resource <- function(
  filename,
  expected_md5,
  minimum_size,
  search_directories
) {
  search_directories <- unique(
    search_directories[
      !is.na(search_directories) &
        dir.exists(
          search_directories
        )
    ]
  )

  for (
    directory_i in search_directories
  ) {
    candidates <- list.files(
      directory_i,
      pattern = paste0(
        "^",
        gsub(
          "\\.",
          "\\\\.",
          filename
        ),
        "$"
      ),
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(candidates) == 0L) {
      next
    }

    for (
      candidate_path in candidates
    ) {
      if (
        resource_is_valid(
          candidate_path,
          expected_md5,
          minimum_size
        )
      ) {
        return(
          normalizePath(
            candidate_path,
            winslash = "/",
            mustWork = TRUE
          )
        )
      }
    }
  }

  NA_character_
}

download_verified_resource <- function(
  filename,
  destination,
  expected_md5,
  minimum_size
) {
  url_candidates <- c(
    paste0(
      "https://zenodo.org/api/records/",
      NICHENET_RECORD_ID,
      "/files/",
      filename,
      "/content"
    ),
    paste0(
      "https://zenodo.org/records/",
      NICHENET_RECORD_ID,
      "/files/",
      filename,
      "?download=1"
    ),
    paste0(
      "https://zenodo.org/record/",
      NICHENET_RECORD_ID,
      "/files/",
      filename,
      "?download=1"
    )
  )

  temporary_file <- paste0(
    destination,
    ".part"
  )

  if (file.exists(temporary_file)) {
    unlink(
      temporary_file,
      force = TRUE
    )
  }

  for (
    url_i in url_candidates
  ) {
    for (
      attempt_i in seq_len(
        RESOURCE_DOWNLOAD_ATTEMPTS
      )
    ) {
      log_msg(
        "Downloading NicheNet resource: ",
        filename,
        " | attempt ",
        attempt_i,
        "/",
        RESOURCE_DOWNLOAD_ATTEMPTS,
        " | ",
        url_i
      )

      success <- tryCatch(
        {
          utils::download.file(
            url_i,
            temporary_file,
            mode = "wb",
            method = "libcurl",
            quiet = FALSE
          )

          resource_is_valid(
            temporary_file,
            expected_md5,
            minimum_size
          )
        },
        error = function(download_error) {
          add_warning(
            "RESOURCE_DOWNLOAD",
            filename,
            conditionMessage(
              download_error
            )
          )

          FALSE
        }
      )

      if (success) {
        if (file.exists(destination)) {
          unlink(
            destination,
            force = TRUE
          )
        }

        moved <- file.rename(
          temporary_file,
          destination
        )

        if (!moved) {
          copied <- file.copy(
            temporary_file,
            destination,
            overwrite = TRUE
          )

          unlink(
            temporary_file,
            force = TRUE
          )

          moved <- copied
        }

        if (
          isTRUE(moved) &&
          resource_is_valid(
            destination,
            expected_md5,
            minimum_size
          )
        ) {
          return(TRUE)
        }
      }

      if (file.exists(temporary_file)) {
        unlink(
          temporary_file,
          force = TRUE
        )
      }

      Sys.sleep(
        min(
          5 * attempt_i,
          20
        )
      )
    }
  }

  FALSE
}

ensure_nichenet_resource <- function(
  resource_definition,
  search_directories
) {
  filename <-
    resource_definition$filename

  destination <-
    resource_definition$destination

  expected_md5 <-
    resource_definition$md5

  minimum_size <-
    resource_definition$minimum_size

  if (
    resource_is_valid(
      destination,
      expected_md5,
      minimum_size
    )
  ) {
    return(
      data.table::data.table(
        filename = filename,
        final_path =
          normalizePath(
            destination,
            winslash = "/",
            mustWork = TRUE
          ),
        acquisition_mode =
          "EXISTING_VERIFIED",
        size_bytes =
          as.numeric(
            file.info(destination)$size
          ),
        expected_md5 =
          expected_md5,
        observed_md5 =
          file_md5(destination),
        valid = TRUE
      )
    )
  }

  local_match <-
    find_verified_local_resource(
      filename,
      expected_md5,
      minimum_size,
      search_directories
    )

  if (
    length(local_match) == 1L &&
    !is.na(local_match)
  ) {
    if (
      normalizePath(
        local_match,
        winslash = "/",
        mustWork = TRUE
      ) !=
        normalizePath(
          destination,
          winslash = "/",
          mustWork = FALSE
        )
    ) {
      copied <- file.copy(
        local_match,
        destination,
        overwrite = TRUE
      )

      if (!copied) {
        stop(
          "Could not copy verified local NicheNet resource: ",
          local_match
        )
      }
    }

    return(
      data.table::data.table(
        filename = filename,
        final_path =
          normalizePath(
            destination,
            winslash = "/",
            mustWork = TRUE
          ),
        acquisition_mode =
          "COPIED_LOCAL_VERIFIED",
        size_bytes =
          as.numeric(
            file.info(destination)$size
          ),
        expected_md5 =
          expected_md5,
        observed_md5 =
          file_md5(destination),
        valid = TRUE
      )
    )
  }

  downloaded <-
    download_verified_resource(
      filename,
      destination,
      expected_md5,
      minimum_size
    )

  if (!downloaded) {
    stop(
      "Unable to obtain the verified NicheNet-v2 resource ",
      filename,
      ". Place the exact Zenodo record 7074291 file at:\n",
      destination,
      "\nExpected MD5: ",
      expected_md5
    )
  }

  data.table::data.table(
    filename = filename,
    final_path =
      normalizePath(
        destination,
        winslash = "/",
        mustWork = TRUE
      ),
    acquisition_mode =
      "DOWNLOADED_VERIFIED",
    size_bytes =
      as.numeric(
        file.info(destination)$size
      ),
    expected_md5 =
      expected_md5,
    observed_md5 =
      file_md5(destination),
    valid = TRUE
  )
}

resource_audit <- data.table::rbindlist(
  lapply(
    NICHENET_FILES,
    ensure_nichenet_resource,
    search_directories =
      NICHENET_SEARCH_DIRS
  ),
  use.names = TRUE,
  fill = TRUE
)

if (
  nrow(resource_audit) !=
    length(NICHENET_FILES) ||
  any(resource_audit$valid != TRUE)
) {
  stop(
    "NicheNet-v2 resource validation failed."
  )
}

write_csv_safe(
  resource_audit,
  file.path(
    DIRS$tables,
    "00_stage6_nichenet_resource_audit.csv"
  )
)

log_msg(
  "Verified NicheNet-v2 resources: ",
  paste(
    resource_audit$filename,
    resource_audit$acquisition_mode,
    sep = "=",
    collapse = "; "
  )
)

############################################################
## 5. Load locked cross-stage inputs
############################################################

sample_meta <- data.table::fread(
  STAGE3_SAMPLE_META_FILE,
  encoding = "UTF-8"
)

sample_meta[
  ,
  condition := factor(
    condition,
    levels = c(
      "Control",
      "HFpEF"
    )
  )
]

data.table::setorder(
  sample_meta,
  condition,
  sample_accession
)

if (
  data.table::uniqueN(
    sample_meta$sample_accession
  ) != 6L ||
  sum(
    sample_meta$condition ==
      "Control"
  ) != 3L ||
  sum(
    sample_meta$condition ==
      "HFpEF"
  ) != 3L
) {
  stop(
    "The locked sample metadata is not the expected 3 + 3 design."
  )
}

stage3_de <- read_table_auto(
  STAGE3_PSEUDOBULK_DE_FILE
)

stage4_programs <- read_table_auto(
  STAGE4_PROGRAM_FILE
)

stage5_all_targets <- read_table_auto(
  STAGE5_ALL_TARGET_FILE
)

stage5_ligand_panel <- data.table::fread(
  STAGE5_LIGAND_PANEL_FILE,
  encoding = "UTF-8"
)

stage5b_rank <- data.table::fread(
  STAGE5B_RANK_FILE,
  encoding = "UTF-8"
)

required_stage3_de_columns <- c(
  "feature",
  "feature_key",
  "major_cell_type",
  "edgeR_logFC",
  "edgeR_pvalue",
  "edgeR_padj",
  "limma_logFC",
  "limma_pvalue",
  "limma_padj",
  "edgeR_limma_sign_agreement"
)

missing_stage3_de_columns <- setdiff(
  required_stage3_de_columns,
  names(stage3_de)
)

if (
  length(
    missing_stage3_de_columns
  ) > 0L
) {
  stop(
    "Stage 3 pseudobulk DE table is missing column(s): ",
    paste(
      missing_stage3_de_columns,
      collapse = ", "
    )
  )
}

required_stage4_program_columns <- c(
  "program_name",
  "subset_name",
  "direction",
  "signature_size",
  "symbol",
  "symbol_key",
  "stage2_disease_lfc",
  "stage2_drug_lfc"
)

missing_stage4_program_columns <- setdiff(
  required_stage4_program_columns,
  names(stage4_programs)
)

if (
  length(
    missing_stage4_program_columns
  ) > 0L
) {
  stop(
    "Stage 4 program manifest is missing column(s): ",
    paste(
      missing_stage4_program_columns,
      collapse = ", "
    )
  )
}

required_stage5_columns <- c(
  "tf_symbol",
  "analysis_role",
  "perturbation_method",
  "perturbation_mode",
  "perturbation_strength",
  "target_feature",
  "target_key",
  "mor",
  "edge_weight",
  "mean_delta_z_HFpEF"
)

missing_stage5_columns <- setdiff(
  required_stage5_columns,
  names(stage5_all_targets)
)

if (
  length(
    missing_stage5_columns
  ) > 0L
) {
  stop(
    "Stage 5 target table is missing column(s): ",
    paste(
      missing_stage5_columns,
      collapse = ", "
    )
  )
}

required_stage5b_columns <- c(
  "tf_symbol",
  "final_robustness_rank",
  "final_robustness_score",
  "positive_recovery_probability",
  "candidate_percentile",
  "inflammation_median_gap_reduction"
)

missing_stage5b_columns <- setdiff(
  required_stage5b_columns,
  names(stage5b_rank)
)

if (
  length(
    missing_stage5b_columns
  ) > 0L
) {
  stop(
    "Stage 5B rank table is missing column(s): ",
    paste(
      missing_stage5b_columns,
      collapse = ", "
    )
  )
}

candidate_rank <- stage5b_rank[
  tf_symbol %in%
    CANDIDATE_TFS
]

if (
  nrow(candidate_rank) !=
    length(CANDIDATE_TFS)
) {
  stop(
    "The three prespecified Stage 6 candidates are not all present in Stage 5B."
  )
}

candidate_rank[
  ,
  candidate_role :=
    unname(
      CANDIDATE_ROLE_MAP[
        tf_symbol
      ]
    )
]

data.table::setorder(
  candidate_rank,
  final_robustness_rank
)

write_csv_safe(
  candidate_rank,
  file.path(
    DIRS$tables,
    "01_stage6_candidate_TF_manifest.csv"
  )
)

############################################################
## 6. Candidate TF-sensitive ligand universe
############################################################

lr_network_raw <- data.table::as.data.table(
  readRDS(
    NICHENET_LR_FILE
  )
)

required_lr_columns <- c(
  "from",
  "to"
)

missing_lr_columns <- setdiff(
  required_lr_columns,
  names(lr_network_raw)
)

if (
  length(
    missing_lr_columns
  ) > 0L
) {
  stop(
    "NicheNet ligand-receptor network is missing column(s): ",
    paste(
      missing_lr_columns,
      collapse = ", "
    )
  )
}

lr_network <- unique(
  lr_network_raw[
    ,
    .(
      ligand =
        as.character(from),
      receptor =
        as.character(to)
    )
  ]
)

lr_network[
  ,
  ligand_key :=
    gene_key(ligand)
]

lr_network[
  ,
  receptor_key :=
    gene_key(receptor)
]

lr_network <- lr_network[
  nzchar(ligand) &
    nzchar(receptor)
]

stage5_primary_targets <- data.table::copy(
  stage5_all_targets[
    tf_symbol %in%
      CANDIDATE_TFS &
      analysis_role ==
        "Biological_candidate" &
      perturbation_method ==
        PRIMARY_PERTURBATION_METHOD &
      perturbation_mode ==
        PRIMARY_PERTURBATION_MODE &
      abs(
        perturbation_strength -
          PRIMARY_PERTURBATION_STRENGTH
      ) <
        1e-12
  ]
)

stage5_primary_targets[
  ,
  absolute_mean_delta_z :=
    abs(mean_delta_z_HFpEF)
]

data.table::setorder(
  stage5_primary_targets,
  tf_symbol,
  target_key,
  -absolute_mean_delta_z
)

stage5_primary_targets <- stage5_primary_targets[
  ,
  .SD[1L],
  by = .(
    tf_symbol,
    target_key
  )
]

stage5_panel_keys <- unique(
  stage5_ligand_panel[
    tf_symbol %in%
      CANDIDATE_TFS &
      perturbation_method ==
        PRIMARY_PERTURBATION_METHOD &
      perturbation_mode ==
        PRIMARY_PERTURBATION_MODE &
      abs(
        perturbation_strength -
          PRIMARY_PERTURBATION_STRENGTH
      ) <
        1e-12,
    .(
      tf_symbol,
      target_key
    )
  ]
)

lr_ligand_map <- unique(
  lr_network[
    ,
    .(
      target_key =
        ligand_key,
      nichenet_ligand =
        ligand
    )
  ],
  by = "target_key"
)

candidate_ligands <- merge(
  stage5_primary_targets,
  lr_ligand_map,
  by = "target_key",
  all = FALSE
)

candidate_ligands <- merge(
  candidate_ligands,
  candidate_rank[
    ,
    .(
      tf_symbol,
      candidate_role,
      final_robustness_rank,
      final_robustness_score,
      positive_recovery_probability,
      candidate_percentile,
      inflammation_median_gap_reduction
    )
  ],
  by = "tf_symbol",
  all.x = TRUE
)

candidate_ligands <- merge(
  candidate_ligands,
  stage5_panel_keys[
    ,
    `:=`(
      in_stage5_curated_ligand_panel =
        TRUE
    )
  ],
  by = c(
    "tf_symbol",
    "target_key"
  ),
  all.x = TRUE
)

candidate_ligands[
  is.na(
    in_stage5_curated_ligand_panel
  ),
  in_stage5_curated_ligand_panel :=
    FALSE
]

candidate_ligands[
  ,
  predicted_receiver_direction :=
    ifelse(
      mean_delta_z_HFpEF < 0,
      "HFpEF_up",
      ifelse(
        mean_delta_z_HFpEF > 0,
        "HFpEF_down",
        "No_direction"
      )
    )
]

candidate_ligands <- candidate_ligands[
  predicted_receiver_direction !=
    "No_direction"
]

data.table::setorder(
  candidate_ligands,
  final_robustness_rank,
  -absolute_mean_delta_z,
  nichenet_ligand
)

candidate_ligand_coverage <- candidate_ligands[
  ,
  .(
    direct_NicheNet_LR_ligands =
      data.table::uniqueN(
        target_key
      ),
    curated_panel_ligands =
      data.table::uniqueN(
        target_key[
          in_stage5_curated_ligand_panel
        ]
      ),
    predicted_decrease_ligands =
      data.table::uniqueN(
        target_key[
          mean_delta_z_HFpEF < 0
        ]
      ),
    predicted_increase_ligands =
      data.table::uniqueN(
        target_key[
          mean_delta_z_HFpEF > 0
        ]
      )
  ),
  by = .(
    tf_symbol,
    candidate_role,
    final_robustness_rank
  )
]

candidate_ligand_coverage <- merge(
  candidate_rank[
    ,
    .(
      tf_symbol,
      candidate_role,
      final_robustness_rank
    )
  ],
  candidate_ligand_coverage,
  by = c(
    "tf_symbol",
    "candidate_role",
    "final_robustness_rank"
  ),
  all.x = TRUE
)

for (
  column_i in c(
    "direct_NicheNet_LR_ligands",
    "curated_panel_ligands",
    "predicted_decrease_ligands",
    "predicted_increase_ligands"
  )
) {
  data.table::set(
    candidate_ligand_coverage,
    which(
      is.na(
        candidate_ligand_coverage[[column_i]]
      )
    ),
    column_i,
    0L
  )
}

if (nrow(candidate_ligands) == 0L) {
  stop(
    "None of the three Stage 6 candidates has a direct Stage 5 target in the verified NicheNet ligand-receptor network."
  )
}

zero_ligand_candidates <- candidate_ligand_coverage[
  direct_NicheNet_LR_ligands < 1L,
  tf_symbol
]

if (length(zero_ligand_candidates) > 0L) {
  add_warning(
    "CANDIDATE_LIGAND_COVERAGE",
    paste(
      zero_ligand_candidates,
      collapse = ";"
    ),
    "No direct Stage 5 target was present in the verified NicheNet ligand-receptor network; zero communication support will be reported rather than forcing an axis."
  )
}

write_csv_safe(
  candidate_ligands,
  file.path(
    DIRS$tables,
    "02_stage6_candidate_TF_sensitive_ligands_pre_expression.csv"
  )
)

write_csv_safe(
  candidate_ligand_coverage,
  file.path(
    DIRS$tables,
    "03_stage6_candidate_ligand_coverage.csv"
  )
)

############################################################
## 7. Seurat object and expression summaries
############################################################

cardiac <- readRDS(
  STAGE3_SEURAT_FILE
)

cardiac <- join_layers_compat(
  cardiac,
  assay = "RNA"
)

cell_metadata <- data.table::as.data.table(
  cardiac@meta.data,
  keep.rownames =
    "cell_id"
)

required_metadata_columns <- c(
  "cell_id",
  "sample_accession",
  "condition",
  "major_cell_type"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(cell_metadata)
)

if (
  length(
    missing_metadata_columns
  ) > 0L
) {
  stop(
    "Stage 3 Seurat metadata is missing column(s): ",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}

cell_metadata[
  ,
  condition := factor(
    condition,
    levels = c(
      "Control",
      "HFpEF"
    )
  )
]

selected_cell_types <- c(
  SENDER_CELL_TYPE,
  RECEIVER_CELL_TYPES
)

cell_count_audit <- cell_metadata[
  major_cell_type %in%
    selected_cell_types,
  .(
    cells = .N
  ),
  by = .(
    sample_accession,
    condition,
    major_cell_type
  )
]

complete_cell_grid <- data.table::CJ(
  sample_accession =
    sample_meta$sample_accession,
  major_cell_type =
    selected_cell_types,
  unique = TRUE
)

complete_cell_grid <- merge(
  complete_cell_grid,
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

cell_count_audit <- merge(
  complete_cell_grid,
  cell_count_audit,
  by = c(
    "sample_accession",
    "condition",
    "major_cell_type"
  ),
  all.x = TRUE
)

cell_count_audit[
  is.na(cells),
  cells := 0L
]

cell_count_audit[
  ,
  minimum_required := ifelse(
    major_cell_type ==
      SENDER_CELL_TYPE,
    MIN_SENDER_CELLS_PER_SAMPLE,
    MIN_RECEIVER_CELLS_PER_SAMPLE
  )
]

cell_count_audit[
  ,
  sufficient := (
    cells >=
      minimum_required
  )
]

if (
  any(
    cell_count_audit$
      sufficient != TRUE
  )
) {
  failed_groups <- cell_count_audit[
    sufficient != TRUE,
    paste0(
      sample_accession,
      ":",
      major_cell_type,
      "=",
      cells
    )
  ]

  stop(
    "Insufficient cells for one or more sender/receiver sample groups: ",
    paste(
      failed_groups,
      collapse = "; "
    )
  )
}

write_csv_safe(
  cell_count_audit,
  file.path(
    DIRS$tables,
    "04_stage6_sender_receiver_cell_count_audit.csv"
  )
)

counts_matrix <- get_assay_matrix(
  cardiac,
  assay = "RNA",
  layer = "counts"
)

if (
  ncol(counts_matrix) !=
    nrow(cell_metadata) ||
  !all(
    colnames(counts_matrix) %in%
      cell_metadata$cell_id
  )
) {
  stop(
    "Seurat count matrix and metadata could not be aligned."
  )
}

cell_metadata <- cell_metadata[
  match(
    colnames(counts_matrix),
    cell_id
  )
]

if (
  any(
    is.na(
      cell_metadata$cell_id
    )
  ) ||
  any(
    cell_metadata$cell_id !=
      colnames(counts_matrix)
  )
) {
  stop(
    "Seurat count matrix cell order does not match metadata."
  )
}

if (
  !"nCount_RNA" %in%
    names(cell_metadata)
) {
  cell_metadata[
    ,
    nCount_RNA :=
      as.numeric(
        Matrix::colSums(
          counts_matrix
        )
      )
  ]
}

feature_map <- make_feature_map(
  rownames(counts_matrix)
)

candidate_ligand_feature_map <-
  map_symbols_to_features(
    candidate_ligands$
      nichenet_ligand,
    feature_map
  )

candidate_ligand_feature_map <- unique(
  candidate_ligand_feature_map[
    !is.na(feature),
    .(
      target_key =
        feature_key,
      sender_ligand_feature =
        feature
    )
  ],
  by = "target_key"
)

candidate_ligands <- merge(
  candidate_ligands,
  candidate_ligand_feature_map,
  by = "target_key",
  all.x = TRUE
)

candidate_lr_network <- merge(
  candidate_ligands[
    ,
    .(
      tf_symbol,
      candidate_role,
      final_robustness_rank,
      final_robustness_score,
      positive_recovery_probability,
      candidate_percentile,
      inflammation_median_gap_reduction,
      target_feature,
      target_key,
      nichenet_ligand,
      sender_ligand_feature,
      mor,
      edge_weight,
      mean_delta_z_HFpEF,
      absolute_mean_delta_z,
      in_stage5_curated_ligand_panel,
      predicted_receiver_direction
    )
  ],
  lr_network,
  by.x = c(
    "nichenet_ligand",
    "target_key"
  ),
  by.y = c(
    "ligand",
    "ligand_key"
  ),
  all.x = FALSE,
  allow.cartesian = TRUE
)

candidate_lr_network[
  ,
  lr_pair_id :=
    paste(
      nichenet_ligand,
      receptor,
      sep = "__"
    )
]

receptor_component_records <- list()

for (
  lr_index in seq_len(
    nrow(candidate_lr_network)
  )
) {
  components <- split_lr_entity(
    candidate_lr_network$
      receptor[
        lr_index
      ]
  )

  if (length(components) == 0L) {
    next
  }

  receptor_component_records[[length(receptor_component_records) + 1L]] <-
    data.table::data.table(
      lr_pair_id =
        candidate_lr_network$
          lr_pair_id[
            lr_index
          ],
      receptor =
        candidate_lr_network$
          receptor[
            lr_index
          ],
      receptor_component =
        components,
      receptor_component_key =
        gene_key(components)
    )
}

receptor_components <- data.table::rbindlist(
  receptor_component_records,
  use.names = TRUE,
  fill = TRUE
)

receptor_component_feature_map <-
  map_symbols_to_features(
    receptor_components$
      receptor_component,
    feature_map
  )

receptor_component_feature_map <- unique(
  receptor_component_feature_map[
    ,
    .(
      receptor_component_key =
        feature_key,
      receptor_component_feature =
        feature
    )
  ],
  by =
    "receptor_component_key"
)

receptor_components <- merge(
  receptor_components,
  receptor_component_feature_map,
  by = "receptor_component_key",
  all.x = TRUE
)

receptor_entity_mapping <- receptor_components[
  ,
  .(
    receptor_components_required =
      data.table::uniqueN(
        receptor_component_key
      ),
    receptor_components_detected =
      data.table::uniqueN(
        receptor_component_key[
          !is.na(
            receptor_component_feature
          )
        ]
      ),
    receptor_component_features =
      paste(
        sort(
          unique(
            receptor_component_feature[
              !is.na(
                receptor_component_feature
              )
            ]
          )
        ),
        collapse = ";"
      ),
    receptor_component_keys =
      paste(
        sort(
          unique(
            receptor_component_key
          )
        ),
        collapse = ";"
      )
  ),
  by = .(
    lr_pair_id,
    receptor
  )
]

receptor_entity_mapping[
  ,
  receptor_mapping_complete := (
    receptor_components_required ==
      receptor_components_detected &
      receptor_components_required >
        0L
  )
]

candidate_lr_network <- merge(
  candidate_lr_network,
  receptor_entity_mapping,
  by = c(
    "lr_pair_id",
    "receptor"
  ),
  all.x = TRUE
)

candidate_lr_network <- candidate_lr_network[
  !is.na(sender_ligand_feature) &
    receptor_mapping_complete ==
      TRUE
]

if (nrow(candidate_lr_network) == 0L) {
  stop(
    "No expressed-feature-mappable NicheNet ligand-receptor pair remains for any Stage 6 candidate."
  )
}

candidate_ligand_features <- sort(
  unique(
    candidate_lr_network$
      sender_ligand_feature
  )
)

receptor_component_features <- sort(
  unique(
    receptor_components[
      !is.na(
        receptor_component_feature
      ),
      receptor_component_feature
    ]
  )
)

summarize_selected_gene_expression <- function(
  count_matrix,
  metadata,
  selected_features,
  selected_cell_types_arg
) {
  records <- list()

  for (
    cell_type_i in
      selected_cell_types_arg
  ) {
    for (
      sample_i in
        sample_meta$sample_accession
    ) {
      cells_i <- metadata[
        major_cell_type ==
          cell_type_i &
          sample_accession ==
            sample_i,
        cell_id
      ]

      cells_i <- intersect(
        cells_i,
        colnames(count_matrix)
      )

      if (length(cells_i) == 0L) {
        next
      }

      matrix_i <- count_matrix[
        selected_features,
        cells_i,
        drop = FALSE
      ]

      total_umi <- sum(
        metadata[
          cell_id %in% cells_i,
          nCount_RNA
        ],
        na.rm = TRUE
      )

      gene_sum <- as.numeric(
        Matrix::rowSums(
          matrix_i
        )
      )

      pct_expressed <- as.numeric(
        Matrix::rowMeans(
          matrix_i > 0
        )
      )

      condition_i <- as.character(
        metadata[
          cell_id ==
            cells_i[1L],
          condition
        ][1L]
      )

      records[[length(records) + 1L]] <-
        data.table::data.table(
          feature =
            rownames(matrix_i),
          feature_key =
            gene_key(
              rownames(matrix_i)
            ),
          major_cell_type =
            cell_type_i,
          sample_accession =
            sample_i,
          condition =
            condition_i,
          cells =
            length(cells_i),
          pct_expressed =
            pct_expressed,
          summed_counts =
            gene_sum,
          total_umi =
            total_umi,
          log2_cpm =
            log2(
              gene_sum /
                pmax(
                  total_umi,
                  1
                ) *
                1e6 +
                1
            )
        )
    }
  }

  data.table::rbindlist(
    records,
    use.names = TRUE,
    fill = TRUE
  )
}

selected_expression_features <- sort(
  unique(
    c(
      candidate_ligand_features,
      receptor_component_features
    )
  )
)

selected_expression_by_sample <-
  summarize_selected_gene_expression(
    counts_matrix,
    cell_metadata,
    selected_expression_features,
    selected_cell_types
  )

write_csv_safe(
  selected_expression_by_sample,
  file.path(
    DIRS$tables,
    "05_stage6_selected_ligand_receptor_expression_by_sample.csv"
  ),
  compress = TRUE
)

receiver_background_records <- list()

for (
  receiver_i in
    RECEIVER_CELL_TYPES
) {
  receiver_cells <- cell_metadata[
    major_cell_type ==
      receiver_i,
    cell_id
  ]

  receiver_cells <- intersect(
    receiver_cells,
    colnames(counts_matrix)
  )

  receiver_matrix <- counts_matrix[
    ,
    receiver_cells,
    drop = FALSE
  ]

  receiver_background_records[[length(receiver_background_records) + 1L]] <-
    data.table::data.table(
      receiver =
        receiver_i,
      feature =
        rownames(
          receiver_matrix
        ),
      feature_key =
        gene_key(
          rownames(
            receiver_matrix
          )
        ),
      pct_expressed =
        as.numeric(
          Matrix::rowMeans(
            receiver_matrix > 0
          )
        )
    )
}

receiver_background_expression <-
  data.table::rbindlist(
    receiver_background_records,
    use.names = TRUE,
    fill = TRUE
  )

write_csv_safe(
  receiver_background_expression[
    pct_expressed >=
      MIN_EXPRESSED_CELL_FRACTION
  ],
  file.path(
    DIRS$tables,
    "06_stage6_receiver_expressed_background_genes.csv"
  ),
  compress = TRUE
)

rm(
  cardiac,
  counts_matrix
)

gc()

############################################################
## 8. Sender ligand and receiver receptor expression gates
############################################################

sender_sample_expression <-
  selected_expression_by_sample[
    major_cell_type ==
      SENDER_CELL_TYPE &
      feature %in%
        candidate_ligand_features
  ]

sender_condition_expression <-
  sender_sample_expression[
    ,
    .(
      median_pct_expressed =
        safe_median(
          pct_expressed,
          0
        ),
      median_log2_cpm =
        safe_median(
          log2_cpm,
          0
        ),
      samples_expressed =
        sum(
          pct_expressed >=
            MIN_EXPRESSED_CELL_FRACTION
        )
    ),
    by = .(
      feature,
      feature_key,
      condition
    )
  ]

sender_condition_wide <- data.table::dcast(
  sender_condition_expression,
  feature +
    feature_key ~
    condition,
  value.var = c(
    "median_pct_expressed",
    "median_log2_cpm",
    "samples_expressed"
  ),
  fill = 0
)

required_sender_wide_columns <- c(
  "median_pct_expressed_Control",
  "median_pct_expressed_HFpEF",
  "median_log2_cpm_Control",
  "median_log2_cpm_HFpEF",
  "samples_expressed_Control",
  "samples_expressed_HFpEF"
)

for (
  missing_column_i in setdiff(
    required_sender_wide_columns,
    names(sender_condition_wide)
  )
) {
  sender_condition_wide[
    ,
    (missing_column_i) := 0
  ]
}

sender_condition_wide[
  ,
  sender_expression_gate := (
    median_pct_expressed_HFpEF >=
      MIN_EXPRESSED_CELL_FRACTION &
      samples_expressed_HFpEF >=
        MIN_HFPEF_SAMPLES_EXPRESSED
  )
]

candidate_ligands <- merge(
  candidate_ligands,
  sender_condition_wide[
    ,
    .(
      sender_ligand_feature =
        feature,
      sender_pct_Control =
        median_pct_expressed_Control,
      sender_pct_HFpEF =
        median_pct_expressed_HFpEF,
      sender_log2CPM_Control =
        median_log2_cpm_Control,
      sender_log2CPM_HFpEF =
        median_log2_cpm_HFpEF,
      sender_samples_expressed_Control =
        samples_expressed_Control,
      sender_samples_expressed_HFpEF =
        samples_expressed_HFpEF,
      sender_expression_gate
    )
  ],
  by = "sender_ligand_feature",
  all.x = TRUE
)

candidate_ligands[
  is.na(sender_expression_gate),
  sender_expression_gate :=
    FALSE
]

receptor_component_sample_expression <-
  selected_expression_by_sample[
    feature %in%
      receptor_component_features &
      major_cell_type %in%
        RECEIVER_CELL_TYPES
  ]

receptor_component_sample_expression <- merge(
  receptor_component_sample_expression,
  receptor_components[
    !is.na(
      receptor_component_feature
    ),
    .(
      lr_pair_id,
      receptor,
      receptor_component_feature
    )
  ],
  by.x = "feature",
  by.y =
    "receptor_component_feature",
  allow.cartesian = TRUE
)

receptor_entity_sample_records <- list()

receptor_entity_groups <- unique(
  candidate_lr_network[
    ,
    .(
      lr_pair_id,
      receptor,
      receptor_components_required
    )
  ]
)

for (
  entity_index in seq_len(
    nrow(receptor_entity_groups)
  )
) {
  entity_i <-
    receptor_entity_groups[
      entity_index
    ]

  expression_i <-
    receptor_component_sample_expression[
      lr_pair_id ==
        entity_i$lr_pair_id
    ]

  if (nrow(expression_i) == 0L) {
    next
  }

  summarized_i <- expression_i[
    ,
    .(
      detected_components =
        data.table::uniqueN(
          feature
        ),
      receptor_pct_expressed =
        min(
          pct_expressed,
          na.rm = TRUE
        ),
      receptor_log2_cpm =
        min(
          log2_cpm,
          na.rm = TRUE
        )
    ),
    by = .(
      lr_pair_id,
      receptor,
      major_cell_type,
      sample_accession,
      condition
    )
  ]

  summarized_i <- summarized_i[
    detected_components ==
      entity_i$
        receptor_components_required
  ]

  receptor_entity_sample_records[[length(receptor_entity_sample_records) + 1L]] <-
    summarized_i
}

receptor_entity_sample_expression <-
  data.table::rbindlist(
    receptor_entity_sample_records,
    use.names = TRUE,
    fill = TRUE
  )

receptor_condition_expression <-
  receptor_entity_sample_expression[
    ,
    .(
      median_receptor_pct =
        safe_median(
          receptor_pct_expressed,
          0
        ),
      median_receptor_log2_cpm =
        safe_median(
          receptor_log2_cpm,
          0
        ),
      receptor_samples_expressed =
        sum(
          receptor_pct_expressed >=
            MIN_EXPRESSED_CELL_FRACTION
        )
    ),
    by = .(
      lr_pair_id,
      receptor,
      receiver =
        major_cell_type,
      condition
    )
  ]

receptor_condition_wide <-
  data.table::dcast(
    receptor_condition_expression,
    lr_pair_id +
      receptor +
      receiver ~
      condition,
    value.var = c(
      "median_receptor_pct",
      "median_receptor_log2_cpm",
      "receptor_samples_expressed"
    ),
    fill = 0
  )

required_receptor_wide_columns <- c(
  "median_receptor_pct_Control",
  "median_receptor_pct_HFpEF",
  "median_receptor_log2_cpm_Control",
  "median_receptor_log2_cpm_HFpEF",
  "receptor_samples_expressed_Control",
  "receptor_samples_expressed_HFpEF"
)

for (
  missing_column_i in setdiff(
    required_receptor_wide_columns,
    names(receptor_condition_wide)
  )
) {
  receptor_condition_wide[
    ,
    (missing_column_i) := 0
  ]
}

receptor_condition_wide[
  ,
  receptor_expression_gate := (
    median_receptor_pct_HFpEF >=
      MIN_EXPRESSED_CELL_FRACTION &
      receptor_samples_expressed_HFpEF >=
        MIN_HFPEF_SAMPLES_EXPRESSED
  )
]

write_csv_safe(
  sender_condition_wide,
  file.path(
    DIRS$tables,
    "07_stage6_macrophage_ligand_expression_summary.csv"
  )
)

write_csv_safe(
  receptor_condition_wide,
  file.path(
    DIRS$tables,
    "08_stage6_receiver_receptor_expression_summary.csv"
  )
)

############################################################
## 9. Stage 2 and Stage 3 ligand-direction support
############################################################

primary_stage2_programs <- stage4_programs[
  signature_size ==
    PRIMARY_STAGE2_SIGNATURE_SIZE &
    !grepl(
      "_Stage3Supported$",
      program_name
    )
]

stage3_supported_ligand_keys <- unique(
  stage4_programs[
    signature_size ==
      PRIMARY_STAGE2_SIGNATURE_SIZE &
      grepl(
        "_Stage3Supported$",
        program_name
      ),
    symbol_key
  ]
)

stage2_ligand_support <- primary_stage2_programs[
  ,
  .(
    stage2_program_count =
      data.table::uniqueN(
        program_name
      ),
    stage2_direction_count =
      data.table::uniqueN(
        direction
      ),
    stage2_direction = if (
      data.table::uniqueN(
        direction
      ) == 1L
    ) {
      unique(direction)[1L]
    } else {
      "Ambiguous"
    },
    stage2_disease_lfc =
      safe_median(
        stage2_disease_lfc
      ),
    stage2_drug_lfc =
      safe_median(
        stage2_drug_lfc
      ),
    stage2_subsets =
      paste(
        sort(
          unique(subset_name)
        ),
        collapse = ";"
      )
  ),
  by = .(
    target_key =
      symbol_key
  )
]

stage2_ligand_support[
  ,
  stage3_supported_any := (
    target_key %in%
      stage3_supported_ligand_keys
  )
]

macrophage_de <- stage3_de[
  major_cell_type ==
    SENDER_CELL_TYPE
]

data.table::setorder(
  macrophage_de,
  feature_key,
  edgeR_padj,
  limma_padj
)

macrophage_de <- macrophage_de[
  ,
  .SD[1L],
  by = feature_key
]

candidate_ligands <- merge(
  candidate_ligands,
  stage2_ligand_support,
  by = "target_key",
  all.x = TRUE
)

candidate_ligands <- merge(
  candidate_ligands,
  macrophage_de[
    ,
    .(
      target_key =
        feature_key,
      stage3_sender_feature =
        feature,
      stage3_sender_edgeR_logFC =
        edgeR_logFC,
      stage3_sender_edgeR_padj =
        edgeR_padj,
      stage3_sender_limma_logFC =
        limma_logFC,
      stage3_sender_limma_padj =
        limma_padj,
      stage3_sender_sign_agreement =
        edgeR_limma_sign_agreement
    )
  ],
  by = "target_key",
  all.x = TRUE
)

candidate_ligands[
  is.na(stage2_program_count),
  stage2_program_count := 0L
]

candidate_ligands[
  is.na(stage3_supported_any),
  stage3_supported_any := FALSE
]

candidate_ligands[
  ,
  stage5_expected_disease_sign :=
    -sign(
      mean_delta_z_HFpEF
    )
]

candidate_ligands[
  ,
  stage2_expected_disease_sign :=
    sign(
      stage2_disease_lfc
    )
]

candidate_ligands[
  ,
  stage3_sender_disease_sign := ifelse(
    stage3_sender_sign_agreement ==
      TRUE &
      is.finite(
        stage3_sender_edgeR_logFC
      ) &
      is.finite(
        stage3_sender_limma_logFC
      ),
    sign(
      (
        stage3_sender_edgeR_logFC +
          stage3_sender_limma_logFC
      ) /
        2
    ),
    NA_real_
  )
]

candidate_ligands[
  ,
  stage5_stage2_direction_match := (
    is.finite(
      stage2_expected_disease_sign
    ) &
      stage2_expected_disease_sign !=
        0 &
      stage5_expected_disease_sign ==
        stage2_expected_disease_sign
  )
]

candidate_ligands[
  ,
  stage3_stage5_direction_match := (
    is.finite(
      stage3_sender_disease_sign
    ) &
      stage3_sender_disease_sign !=
        0 &
      stage3_sender_disease_sign ==
        stage5_expected_disease_sign
  )
]

candidate_ligands[
  is.na(stage5_stage2_direction_match),
  stage5_stage2_direction_match := FALSE
]

candidate_ligands[
  is.na(stage3_stage5_direction_match),
  stage3_stage5_direction_match := FALSE
]

candidate_ligands[
  ,
  strict_cross_stage_support := (
    sender_expression_gate ==
      TRUE &
      stage5_stage2_direction_match ==
        TRUE &
      stage3_stage5_direction_match ==
        TRUE &
      stage3_supported_any ==
        TRUE
  )
]

candidate_ligands[
  ,
  support_tier := data.table::fcase(
    strict_cross_stage_support ==
      TRUE,
    "Tier_A_strict_cross_stage",
    sender_expression_gate ==
      TRUE &
      stage3_stage5_direction_match ==
        TRUE,
    "Tier_B_Stage3_direction_supported",
    sender_expression_gate ==
      TRUE,
    "Tier_C_sender_expression_only",
    default =
      "Tier_D_not_sender_expressed"
  )
]

write_csv_safe(
  candidate_ligands,
  file.path(
    DIRS$tables,
    "09_stage6_candidate_ligand_cross_stage_support.csv"
  )
)

############################################################
## 10. Receiver disease-response gene sets
############################################################

receiver_de <- stage3_de[
  major_cell_type %in%
    RECEIVER_CELL_TYPES
]

receiver_de[
  ,
  mean_logFC :=
    rowMeans(
      cbind(
        edgeR_logFC,
        limma_logFC
      ),
      na.rm = TRUE
    )
]

receiver_de[
  ,
  abs_mean_logFC :=
    abs(mean_logFC)
]

receiver_de[
  ,
  combined_p := pmin(
    edgeR_pvalue,
    limma_pvalue,
    na.rm = TRUE
  )
]

receiver_de[
  !is.finite(combined_p),
  combined_p := 1
]

select_receiver_geneset <- function(
  receiver_arg,
  direction_arg,
  receiver_de_table,
  background_table
) {
  sign_required <- if (
    direction_arg ==
      "HFpEF_up"
  ) {
    1
  } else {
    -1
  }

  background_keys <- unique(
    background_table[
      receiver ==
        receiver_arg &
        pct_expressed >=
          MIN_EXPRESSED_CELL_FRACTION,
      feature_key
    ]
  )

  de_i <- data.table::copy(
    receiver_de_table[
      major_cell_type ==
        receiver_arg &
        feature_key %in%
          background_keys &
        edgeR_limma_sign_agreement ==
          TRUE &
        sign(mean_logFC) ==
          sign_required
    ]
  )

  primary_i <- de_i[
    abs(mean_logFC) >=
      RECEIVER_MIN_ABS_LOGFC &
      (
        (
          is.finite(edgeR_padj) &
            edgeR_padj <=
              RECEIVER_PRIMARY_FDR
        ) |
          (
            is.finite(limma_padj) &
              limma_padj <=
                RECEIVER_PRIMARY_FDR
          )
      )
  ]

  data.table::setorderv(
    primary_i,
    cols = c(
      "combined_p",
      "abs_mean_logFC"
    ),
    order = c(
      1L,
      -1L
    ),
    na.last = TRUE
  )

  selection_tier <-
    "FDR_0.10_and_logFC"

  selected_i <- primary_i

  if (
    nrow(selected_i) <
      MIN_RECEIVER_TARGET_GENES
  ) {
    fallback_i <- de_i[
      abs(mean_logFC) >=
        RECEIVER_MIN_ABS_LOGFC
    ]

    data.table::setorderv(
      fallback_i,
      cols = c(
        "combined_p",
        "abs_mean_logFC"
      ),
      order = c(
        1L,
        -1L
      ),
      na.last = TRUE
    )

    selected_i <- head(
      fallback_i,
      min(
        FALLBACK_RECEIVER_TARGET_GENES,
        nrow(fallback_i)
      )
    )

    selection_tier <-
      "Exploratory_sign_agree_logFC_ranked"
  }

  if (
    nrow(selected_i) <
      MIN_RECEIVER_TARGET_GENES
  ) {
    data.table::setorderv(
      de_i,
      cols = c(
        "combined_p",
        "abs_mean_logFC"
      ),
      order = c(
        1L,
        -1L
      ),
      na.last = TRUE
    )

    selected_i <- head(
      de_i,
      min(
        FALLBACK_RECEIVER_TARGET_GENES,
        nrow(de_i)
      )
    )

    selection_tier <-
      "Exploratory_sign_agree_ranked"
  }

  selected_i <- head(
    selected_i,
    min(
      MAX_RECEIVER_TARGET_GENES,
      nrow(selected_i)
    )
  )

  selected_i[
    ,
    `:=`(
      receiver =
        receiver_arg,
      receiver_direction =
        direction_arg,
      selection_tier =
        selection_tier
    )
  ]

  selected_i
}

receiver_gene_set_records <- list()

for (
  receiver_i in
    RECEIVER_CELL_TYPES
) {
  for (
    direction_i in c(
      "HFpEF_up",
      "HFpEF_down"
    )
  ) {
    selected_i <- select_receiver_geneset(
      receiver_i,
      direction_i,
      receiver_de,
      receiver_background_expression
    )

    receiver_gene_set_records[[length(receiver_gene_set_records) + 1L]] <-
      selected_i
  }
}

receiver_gene_sets <- data.table::rbindlist(
  receiver_gene_set_records,
  use.names = TRUE,
  fill = TRUE
)

receiver_gene_set_summary <- receiver_gene_sets[
  ,
  .(
    selected_genes =
      data.table::uniqueN(
        feature_key
      ),
    median_abs_logFC =
      safe_median(
        abs(mean_logFC),
        0
      ),
    minimum_edgeR_padj =
      suppressWarnings(
        min(
          edgeR_padj,
          na.rm = TRUE
        )
      ),
    minimum_limma_padj =
      suppressWarnings(
        min(
          limma_padj,
          na.rm = TRUE
        )
      ),
    selection_tier =
      unique(selection_tier)[1L]
  ),
  by = .(
    receiver,
    receiver_direction
  )
]

receiver_gene_set_grid <- data.table::CJ(
  receiver =
    RECEIVER_CELL_TYPES,
  receiver_direction = c(
    "HFpEF_up",
    "HFpEF_down"
  ),
  unique = TRUE
)

receiver_gene_set_summary <- merge(
  receiver_gene_set_grid,
  receiver_gene_set_summary,
  by = c(
    "receiver",
    "receiver_direction"
  ),
  all.x = TRUE
)

receiver_gene_set_summary[
  is.na(selected_genes),
  selected_genes := 0L
]

receiver_gene_set_summary[
  is.na(selection_tier),
  selection_tier :=
    "No_eligible_receiver_genes"
]

receiver_gene_set_summary[
  !is.finite(
    minimum_edgeR_padj
  ),
  minimum_edgeR_padj :=
    NA_real_
]

receiver_gene_set_summary[
  !is.finite(
    minimum_limma_padj
  ),
  minimum_limma_padj :=
    NA_real_
]

write_csv_safe(
  receiver_gene_sets,
  file.path(
    DIRS$tables,
    "10_stage6_receiver_disease_response_gene_sets.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  receiver_gene_set_summary,
  file.path(
    DIRS$tables,
    "11_stage6_receiver_gene_set_summary.csv"
  )
)

############################################################
## 11. Load and subset the NicheNet ligand-target matrix
############################################################

log_msg(
  "Loading verified NicheNet-v2 mouse ligand-target matrix."
)

ligand_target_matrix <- readRDS(
  NICHENET_MATRIX_FILE
)

if (
  is.null(
    rownames(
      ligand_target_matrix
    )
  ) ||
  is.null(
    colnames(
      ligand_target_matrix
    )
  ) ||
  nrow(ligand_target_matrix) <
    1000L ||
  ncol(ligand_target_matrix) <
    100L
) {
  stop(
    "The NicheNet ligand-target matrix has an invalid structure."
  )
}

ligand_matrix_map <- data.table::data.table(
  nichenet_ligand =
    colnames(
      ligand_target_matrix
    ),
  ligand_matrix_key =
    gene_key(
      colnames(
        ligand_target_matrix
      )
    )
)

data.table::setorder(
  ligand_matrix_map,
  ligand_matrix_key,
  nichenet_ligand
)

ligand_matrix_map <- ligand_matrix_map[
  ,
  .SD[1L],
  by = ligand_matrix_key
]

candidate_ligands[
  ,
  ligand_matrix_key :=
    gene_key(
      nichenet_ligand
    )
]

candidate_ligands <- merge(
  candidate_ligands,
  ligand_matrix_map,
  by = "ligand_matrix_key",
  all.x = TRUE,
  suffixes = c(
    "",
    "_matrix"
  )
)

if (
  "nichenet_ligand_matrix" %in%
    names(candidate_ligands)
) {
  candidate_ligands[
    ,
    ligand_matrix_symbol :=
      nichenet_ligand_matrix
  ]
} else {
  candidate_ligands[
    ,
    ligand_matrix_symbol :=
      nichenet_ligand
  ]
}

candidate_ligands <- candidate_ligands[
  !is.na(
    ligand_matrix_symbol
  )
]

matrix_candidate_ligands <- sort(
  unique(
    candidate_ligands$
      ligand_matrix_symbol
  )
)

if (
  length(
    matrix_candidate_ligands
  ) < 1L
) {
  stop(
    "No Stage 5 candidate ligand is present in the verified NicheNet ligand-target matrix."
  )
}

ligand_target_subset <- ligand_target_matrix[
  ,
  matrix_candidate_ligands,
  drop = FALSE
]

matrix_target_map <- data.table::data.table(
  matrix_target =
    rownames(
      ligand_target_subset
    ),
  target_key =
    gene_key(
      rownames(
        ligand_target_subset
      )
    )
)

data.table::setorder(
  matrix_target_map,
  target_key,
  matrix_target
)

matrix_target_map <- matrix_target_map[
  ,
  .SD[1L],
  by = target_key
]

rm(
  ligand_target_matrix
)

gc()

############################################################
## 12. NicheNet-v2 ligand activity
############################################################

receiver_background_matrix_map <- merge(
  receiver_background_expression[
    pct_expressed >=
      MIN_EXPRESSED_CELL_FRACTION,
    .(
      receiver,
      target_key =
        feature_key,
      receiver_pct_expressed =
        pct_expressed
    )
  ],
  matrix_target_map,
  by = "target_key",
  all = FALSE
)

receiver_gene_sets_matrix <- merge(
  receiver_gene_sets,
  matrix_target_map,
  by.x = "feature_key",
  by.y = "target_key",
  all = FALSE
)

candidate_ligand_activity_input <- unique(
  candidate_ligands[
    sender_expression_gate ==
      TRUE,
    .(
      ligand_matrix_symbol,
      predicted_receiver_direction
    )
  ]
)

ligand_activity_records <- list()

for (
  receiver_i in
    RECEIVER_CELL_TYPES
) {
  for (
    direction_i in c(
      "HFpEF_up",
      "HFpEF_down"
    )
  ) {
    background_i <- unique(
      receiver_background_matrix_map[
        receiver ==
          receiver_i,
        matrix_target
      ]
    )

    gene_set_i <- unique(
      receiver_gene_sets_matrix[
        receiver ==
          receiver_i &
          receiver_direction ==
            direction_i,
        matrix_target
      ]
    )

    gene_set_i <- intersect(
      gene_set_i,
      background_i
    )

    ligand_input_i <-
      candidate_ligand_activity_input[
        predicted_receiver_direction ==
          direction_i,
        ligand_matrix_symbol
      ]

    ligand_input_i <- intersect(
      ligand_input_i,
      colnames(
        ligand_target_subset
      )
    )

    if (
      length(background_i) <
        20L ||
      length(gene_set_i) <
        2L ||
      length(ligand_input_i) <
        1L
    ) {
      add_warning(
        "NICHENET_ACTIVITY",
        paste(
          receiver_i,
          direction_i,
          sep = ":"
        ),
        paste0(
          "Insufficient background, gene-set, or candidate-ligand coverage: background=",
          length(background_i),
          "; gene_set=",
          length(gene_set_i),
          "; ligands=",
          length(ligand_input_i)
        )
      )

      next
    }

    labels <- as.integer(
      background_i %in%
        gene_set_i
    )

    baseline_aupr <- mean(
      labels
    )

    for (
      ligand_i in
        ligand_input_i
    ) {
      scores_i <- as.numeric(
        ligand_target_subset[
          background_i,
          ligand_i
        ]
      )

      auroc_i <- simple_auc(
        labels,
        scores_i
      )

      aupr_i <- simple_aupr(
        labels,
        scores_i
      )

      pearson_i <- safe_pearson(
        scores_i,
        labels
      )

      spearman_i <- safe_spearman(
        scores_i,
        labels
      )

      ligand_activity_records[[length(ligand_activity_records) + 1L]] <-
        data.table::data.table(
          receiver =
            receiver_i,
          receiver_direction =
            direction_i,
          ligand_matrix_symbol =
            ligand_i,
          background_genes =
            length(
              background_i
            ),
          receiver_gene_set_genes =
            length(
              gene_set_i
            ),
          positive_fraction =
            baseline_aupr,
          auroc =
            auroc_i,
          aupr =
            aupr_i,
          aupr_corrected =
            aupr_i -
            baseline_aupr,
          pearson =
            pearson_i,
          spearman =
            spearman_i
        )
    }
  }
}

ligand_activity <- data.table::rbindlist(
  ligand_activity_records,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(ligand_activity) == 0L) {
  stop(
    "No NicheNet ligand activity result was generated."
  )
}

ligand_activity[
  ,
  rank_aupr_corrected :=
    rank_metric(
      aupr_corrected,
      higher_is_better = TRUE
    ),
  by = .(
    receiver,
    receiver_direction
  )
]

ligand_activity[
  ,
  rank_pearson :=
    rank_metric(
      pearson,
      higher_is_better = TRUE
    ),
  by = .(
    receiver,
    receiver_direction
  )
]

ligand_activity[
  ,
  rank_auroc :=
    rank_metric(
      auroc,
      higher_is_better = TRUE
    ),
  by = .(
    receiver,
    receiver_direction
  )
]

ligand_activity[
  ,
  nichenet_mean_rank :=
    rowMeans(
      .SD,
      na.rm = TRUE
    ),
  .SDcols = c(
    "rank_aupr_corrected",
    "rank_pearson",
    "rank_auroc"
  )
]

ligand_activity[
  ,
  nichenet_activity_rank :=
    rank(
      nichenet_mean_rank,
      ties.method = "average"
    ),
  by = .(
    receiver,
    receiver_direction
  )
]

data.table::setorder(
  ligand_activity,
  receiver,
  receiver_direction,
  nichenet_activity_rank
)

write_csv_safe(
  ligand_activity,
  file.path(
    DIRS$tables,
    "12_stage6_nichenet_ligand_activity.csv"
  )
)

############################################################
## 13. NicheNet ligand-target links
############################################################

ligand_target_link_records <- list()

for (
  activity_index in seq_len(
    nrow(ligand_activity)
  )
) {
  activity_i <-
    ligand_activity[
      activity_index
    ]

  target_table_i <- receiver_gene_sets_matrix[
    receiver ==
      activity_i$receiver &
      receiver_direction ==
        activity_i$
          receiver_direction
  ]

  target_table_i <- target_table_i[
    matrix_target %in%
      rownames(
        ligand_target_subset
      )
  ]

  if (nrow(target_table_i) == 0L) {
    next
  }

  target_table_i[
    ,
    regulatory_potential :=
      as.numeric(
        ligand_target_subset[
          matrix_target,
          activity_i$
            ligand_matrix_symbol
        ]
      )
  ]

  target_table_i <- target_table_i[
    is.finite(
      regulatory_potential
    ) &
      regulatory_potential > 0
  ]

  data.table::setorderv(
    target_table_i,
    cols = c(
      "regulatory_potential",
      "combined_p",
      "abs_mean_logFC"
    ),
    order = c(
      -1L,
      1L,
      -1L
    ),
    na.last = TRUE
  )

  target_table_i <- head(
    target_table_i,
    min(
      TOP_TARGET_LINKS_PER_LIGAND,
      nrow(target_table_i)
    )
  )

  target_table_i[
    ,
    `:=`(
      ligand_matrix_symbol =
        activity_i$
          ligand_matrix_symbol,
      ligand_aupr_corrected =
        activity_i$
          aupr_corrected,
      ligand_pearson =
        activity_i$
          pearson,
      ligand_activity_rank =
        activity_i$
          nichenet_activity_rank
    )
  ]

  ligand_target_link_records[[length(ligand_target_link_records) + 1L]] <-
    target_table_i
}

ligand_target_links <- data.table::rbindlist(
  ligand_target_link_records,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(ligand_target_links) == 0L) {
  stop(
    "No positive NicheNet ligand-target link was generated."
  )
}

ligand_target_summary <- ligand_target_links[
  ,
  .(
    target_links =
      data.table::uniqueN(
        feature_key
      ),
    maximum_regulatory_potential =
      max(
        regulatory_potential,
        na.rm = TRUE
      ),
    median_regulatory_potential =
      safe_median(
        regulatory_potential,
        0
      ),
    median_receiver_abs_logFC =
      safe_median(
        abs(mean_logFC),
        0
      )
  ),
  by = .(
    ligand_matrix_symbol,
    receiver,
    receiver_direction
  )
]

write_csv_safe(
  ligand_target_links,
  file.path(
    DIRS$tables,
    "13_stage6_nichenet_ligand_target_links.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  ligand_target_summary,
  file.path(
    DIRS$tables,
    "14_stage6_nichenet_ligand_target_summary.csv"
  )
)

############################################################
## 14. Ligand-receptor axes and sample-level co-variation
############################################################

candidate_ligands_for_axes <- candidate_ligands[
  sender_expression_gate ==
    TRUE
]

candidate_axes <- merge(
  candidate_lr_network,
  candidate_ligands_for_axes[
    ,
    .(
      tf_symbol,
      target_key,
      nichenet_ligand,
      ligand_matrix_symbol,
      sender_pct_Control,
      sender_pct_HFpEF,
      sender_log2CPM_Control,
      sender_log2CPM_HFpEF,
      sender_samples_expressed_Control,
      sender_samples_expressed_HFpEF,
      sender_expression_gate,
      stage2_program_count,
      stage2_direction,
      stage2_disease_lfc,
      stage2_drug_lfc,
      stage3_supported_any,
      stage3_sender_edgeR_logFC,
      stage3_sender_edgeR_padj,
      stage3_sender_limma_logFC,
      stage3_sender_limma_padj,
      stage3_sender_sign_agreement,
      stage5_expected_disease_sign,
      stage2_expected_disease_sign,
      stage3_sender_disease_sign,
      stage5_stage2_direction_match,
      stage3_stage5_direction_match,
      strict_cross_stage_support,
      support_tier
    )
  ],
  by = c(
    "tf_symbol",
    "target_key",
    "nichenet_ligand"
  ),
  all.x = FALSE
)

candidate_axes <- merge(
  candidate_axes,
  receptor_condition_wide,
  by = c(
    "lr_pair_id",
    "receptor"
  ),
  allow.cartesian = TRUE
)

candidate_axes <- candidate_axes[
  receptor_expression_gate ==
    TRUE
]

candidate_axes <- merge(
  candidate_axes,
  ligand_activity,
  by.x = c(
    "ligand_matrix_symbol",
    "predicted_receiver_direction",
    "receiver"
  ),
  by.y = c(
    "ligand_matrix_symbol",
    "receiver_direction",
    "receiver"
  ),
  allow.cartesian = TRUE
)

if (nrow(candidate_axes) == 0L) {
  stop(
    "No expression- and NicheNet-supported candidate communication axis was generated."
  )
}

missing_axis_candidates <- setdiff(
  CANDIDATE_TFS,
  unique(
    candidate_axes$tf_symbol
  )
)

if (length(missing_axis_candidates) > 0L) {
  add_warning(
    "CANDIDATE_AXIS_COVERAGE",
    paste(
      missing_axis_candidates,
      collapse = ";"
    ),
    "No expression- and NicheNet-supported communication axis was identified; the candidate will remain in the summary with zero axes."
  )
}

candidate_axes <- merge(
  candidate_axes,
  ligand_target_summary,
  by.x = c(
    "ligand_matrix_symbol",
    "receiver",
    "predicted_receiver_direction"
  ),
  by.y = c(
    "ligand_matrix_symbol",
    "receiver",
    "receiver_direction"
  ),
  all.x = TRUE
)

sender_sample_for_correlation <-
  sender_sample_expression[
    ,
    .(
      sender_ligand_feature =
        feature,
      sample_accession,
      condition,
      sender_sample_pct =
        pct_expressed,
      sender_sample_log2_cpm =
        log2_cpm
    )
  ]

axis_correlation_keys <- unique(
  candidate_axes[
    ,
    .(
      sender_ligand_feature,
      lr_pair_id,
      receptor,
      receiver
    )
  ]
)

axis_correlation_records <- list()

for (
  correlation_index in seq_len(
    nrow(axis_correlation_keys)
  )
) {
  key_i <-
    axis_correlation_keys[
      correlation_index
    ]

  sender_i <-
    sender_sample_for_correlation[
      sender_ligand_feature ==
        key_i$
          sender_ligand_feature
    ]

  receptor_i <-
    receptor_entity_sample_expression[
      lr_pair_id ==
        key_i$lr_pair_id &
      receptor ==
        key_i$receptor &
      major_cell_type ==
        key_i$receiver,
      .(
        sample_accession,
        condition,
        receptor_sample_pct =
          receptor_pct_expressed,
        receptor_sample_log2_cpm =
          receptor_log2_cpm
      )
    ]

  paired_i <- merge(
    sender_i,
    receptor_i,
    by = c(
      "sample_accession",
      "condition"
    ),
    all = FALSE
  )

  axis_correlation_records[[length(axis_correlation_records) + 1L]] <-
    data.table::data.table(
      sender_ligand_feature =
        key_i$
          sender_ligand_feature,
      lr_pair_id =
        key_i$lr_pair_id,
      receptor =
        key_i$receptor,
      receiver =
        key_i$receiver,
      paired_samples =
        nrow(paired_i),
      sample_log2CPM_spearman =
        safe_spearman(
          paired_i$
            sender_sample_log2_cpm,
          paired_i$
            receptor_sample_log2_cpm
        ),
      sample_pct_spearman =
        safe_spearman(
          paired_i$
            sender_sample_pct,
          paired_i$
            receptor_sample_pct
        ),
      sender_HFpEF_minus_Control_log2CPM =
        safe_mean(
          paired_i[
            condition ==
              "HFpEF",
            sender_sample_log2_cpm
          ]
        ) -
        safe_mean(
          paired_i[
            condition ==
              "Control",
            sender_sample_log2_cpm
          ]
        ),
      receptor_HFpEF_minus_Control_log2CPM =
        safe_mean(
          paired_i[
            condition ==
              "HFpEF",
            receptor_sample_log2_cpm
          ]
        ) -
        safe_mean(
          paired_i[
            condition ==
              "Control",
            receptor_sample_log2_cpm
          ]
        )
    )
}

axis_sample_correlation <- data.table::rbindlist(
  axis_correlation_records,
  use.names = TRUE,
  fill = TRUE
)

candidate_axes <- merge(
  candidate_axes,
  axis_sample_correlation,
  by = c(
    "sender_ligand_feature",
    "lr_pair_id",
    "receptor",
    "receiver"
  ),
  all.x = TRUE
)

candidate_axes[
  ,
  support_tier_numeric :=
    data.table::fcase(
      strict_cross_stage_support ==
        TRUE,
      1,
      stage3_stage5_direction_match ==
        TRUE,
      2,
      sender_expression_gate ==
        TRUE,
      3,
      default = 4
    )
]

candidate_axes[
  ,
  rank_TF_robustness :=
    rank_metric(
      final_robustness_rank,
      higher_is_better = FALSE
    )
]

candidate_axes[
  ,
  rank_ligand_change :=
    rank_metric(
      absolute_mean_delta_z,
      higher_is_better = TRUE
    )
]

candidate_axes[
  ,
  rank_sender_expression :=
    rank_metric(
      sender_pct_HFpEF,
      higher_is_better = TRUE
    )
]

candidate_axes[
  ,
  rank_receptor_expression :=
    rank_metric(
      median_receptor_pct_HFpEF,
      higher_is_better = TRUE
    )
]

candidate_axes[
  ,
  rank_NicheNet_AUPR :=
    rank_metric(
      aupr_corrected,
      higher_is_better = TRUE
    )
]

candidate_axes[
  ,
  rank_NicheNet_pearson :=
    rank_metric(
      pearson,
      higher_is_better = TRUE
    )
]

candidate_axes[
  ,
  rank_target_support :=
    rank_metric(
      maximum_regulatory_potential,
      higher_is_better = TRUE
    )
]

candidate_axes[
  ,
  rank_cross_stage_support :=
    rank_metric(
      support_tier_numeric,
      higher_is_better = FALSE
    )
]

primary_rank_columns <- c(
  "rank_TF_robustness",
  "rank_ligand_change",
  "rank_sender_expression",
  "rank_receptor_expression",
  "rank_NicheNet_AUPR",
  "rank_NicheNet_pearson",
  "rank_target_support",
  "rank_cross_stage_support"
)

candidate_axes[
  ,
  mean_evidence_rank :=
    rowMeans(
      .SD,
      na.rm = TRUE
    ),
  .SDcols =
    primary_rank_columns
]

candidate_axes[
  ,
  median_evidence_rank :=
    apply(
      .SD,
      1L,
      stats::median,
      na.rm = TRUE
    ),
  .SDcols =
    primary_rank_columns
]

candidate_axes[
  ,
  final_axis_score :=
    0.5 *
      mean_evidence_rank +
    0.5 *
      median_evidence_rank
]

data.table::setorder(
  candidate_axes,
  final_axis_score,
  rank_NicheNet_AUPR,
  rank_cross_stage_support,
  final_robustness_rank
)

candidate_axes[
  ,
  final_axis_rank :=
    seq_len(.N)
]

candidate_axes[
  ,
  axis_id :=
    paste(
      tf_symbol,
      nichenet_ligand,
      receptor,
      receiver,
      sep = "__"
    )
]

write_csv_safe(
  candidate_axes,
  file.path(
    DIRS$tables,
    "15_stage6_candidate_TF_ligand_receptor_axes.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  head(
    candidate_axes,
    min(
      TOP_AXES_FOR_REPORT,
      nrow(candidate_axes)
    )
  ),
  file.path(
    DIRS$tables,
    "16_stage6_top_candidate_axes.csv"
  )
)

############################################################
## 15. Axis ranking sensitivity and candidate summary
############################################################

scenario_definitions <- list(
  All_evidence =
    primary_rank_columns,
  Without_TF_prior = c(
    "rank_ligand_change",
    "rank_sender_expression",
    "rank_receptor_expression",
    "rank_NicheNet_AUPR",
    "rank_NicheNet_pearson",
    "rank_target_support",
    "rank_cross_stage_support"
  ),
  NicheNet_expression = c(
    "rank_sender_expression",
    "rank_receptor_expression",
    "rank_NicheNet_AUPR",
    "rank_NicheNet_pearson",
    "rank_target_support"
  ),
  Cross_stage_emphasis = c(
    "rank_TF_robustness",
    "rank_ligand_change",
    "rank_sender_expression",
    "rank_receptor_expression",
    "rank_NicheNet_AUPR",
    "rank_NicheNet_pearson",
    "rank_target_support",
    "rank_cross_stage_support",
    "rank_cross_stage_support"
  )
)

scenario_records <- list()

for (
  scenario_name in names(
    scenario_definitions
  )
) {
  columns_i <-
    scenario_definitions[[scenario_name]]

  scenario_table <- candidate_axes[
    ,
    .(
      axis_id,
      tf_symbol,
      nichenet_ligand,
      receptor,
      receiver
    )
  ]

  scenario_values <- as.matrix(
    as.data.frame(
      candidate_axes
    )[
      ,
      columns_i,
      drop = FALSE
    ]
  )

  if (
    nrow(scenario_values) !=
      nrow(scenario_table)
  ) {
    stop(
      "Scenario-ranking matrix row mismatch for ",
      scenario_name,
      "."
    )
  }

  scenario_table[
    ,
    scenario_score :=
      rowMeans(
        scenario_values,
        na.rm = TRUE
      )
  ]

  data.table::setorder(
    scenario_table,
    scenario_score,
    axis_id
  )

  scenario_table[
    ,
    scenario_rank :=
      seq_len(.N)
  ]

  scenario_table[
    ,
    scenario :=
      scenario_name
  ]

  scenario_records[[length(scenario_records) + 1L]] <-
    scenario_table
}

axis_ranking_sensitivity <- data.table::rbindlist(
  scenario_records,
  use.names = TRUE,
  fill = TRUE
)

axis_ranking_stability <- axis_ranking_sensitivity[
  ,
  .(
    scenarios =
      .N,
    median_scenario_rank =
      safe_median(
        scenario_rank
      ),
    best_scenario_rank =
      min(
        scenario_rank,
        na.rm = TRUE
      ),
    worst_scenario_rank =
      max(
        scenario_rank,
        na.rm = TRUE
      ),
    top10_scenario_frequency =
      safe_mean(
        scenario_rank <=
          10
      ),
    top20_scenario_frequency =
      safe_mean(
        scenario_rank <=
          20
      )
  ),
  by = .(
    axis_id,
    tf_symbol,
    nichenet_ligand,
    receptor,
    receiver
  )
]

candidate_summary <- candidate_axes[
  ,
  .(
    candidate_role =
      unique(candidate_role)[1L],
    Stage5B_rank =
      unique(
        final_robustness_rank
      )[1L],
    candidate_ligands =
      data.table::uniqueN(
        nichenet_ligand
      ),
    candidate_receptors =
      data.table::uniqueN(
        receptor
      ),
    receivers_with_axes =
      data.table::uniqueN(
        receiver
      ),
    total_axes =
      .N,
    strict_cross_stage_axes =
      sum(
        strict_cross_stage_support ==
          TRUE,
        na.rm = TRUE
      ),
    best_axis_rank =
      min(
        final_axis_rank,
        na.rm = TRUE
      ),
    median_AUPR_corrected =
      safe_median(
        aupr_corrected
      ),
    median_NicheNet_pearson =
      safe_median(
        pearson
      ),
    median_sender_pct_HFpEF =
      safe_median(
        sender_pct_HFpEF
      ),
    median_receptor_pct_HFpEF =
      safe_median(
        median_receptor_pct_HFpEF
      ),
    median_sample_log2CPM_spearman =
      safe_median(
        sample_log2CPM_spearman
      )
  ),
  by = tf_symbol
]

candidate_summary <- merge(
  candidate_rank[
    ,
    .(
      tf_symbol,
      candidate_role,
      Stage5B_rank =
        final_robustness_rank
    )
  ],
  candidate_summary,
  by = c(
    "tf_symbol",
    "candidate_role",
    "Stage5B_rank"
  ),
  all.x = TRUE
)

for (
  column_i in c(
    "candidate_ligands",
    "candidate_receptors",
    "receivers_with_axes",
    "total_axes",
    "strict_cross_stage_axes"
  )
) {
  data.table::set(
    candidate_summary,
    which(
      is.na(
        candidate_summary[[column_i]]
      )
    ),
    column_i,
    0L
  )
}

write_csv_safe(
  axis_ranking_sensitivity,
  file.path(
    DIRS$tables,
    "17_stage6_axis_ranking_sensitivity_scenarios.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  axis_ranking_stability,
  file.path(
    DIRS$tables,
    "18_stage6_axis_ranking_stability_summary.csv"
  )
)

write_csv_safe(
  candidate_summary,
  file.path(
    DIRS$tables,
    "19_stage6_candidate_TF_communication_summary.csv"
  )
)

############################################################
## 16. Scientific checkpoint before figures
############################################################

saveRDS(
  list(
    resource_audit =
      resource_audit,
    candidate_rank =
      candidate_rank,
    candidate_ligand_coverage =
      candidate_ligand_coverage,
    candidate_ligands =
      candidate_ligands,
    cell_count_audit =
      cell_count_audit,
    sender_condition_wide =
      sender_condition_wide,
    receptor_condition_wide =
      receptor_condition_wide,
    receiver_gene_sets =
      receiver_gene_sets,
    receiver_gene_set_summary =
      receiver_gene_set_summary,
    ligand_activity =
      ligand_activity,
    ligand_target_links =
      ligand_target_links,
    candidate_axes =
      candidate_axes,
    axis_ranking_sensitivity =
      axis_ranking_sensitivity,
    axis_ranking_stability =
      axis_ranking_stability,
    candidate_summary =
      candidate_summary
  ),
  file.path(
    DIRS$objects,
    "CHECKPOINT_stage6_scientific_results_pre_figures.rds"
  ),
  compress = FALSE
)

log_msg(
  "Stage 6 scientific calculations checkpointed before figures."
)

############################################################
## 17. Figures and source data
############################################################

role_labels <- c(
  Nfkb1 = "Nfkb1",
  Rela = "Rela",
  Bhlhe40 = "Bhlhe40"
)

top_axes_figure <- head(
  candidate_axes,
  min(
    TOP_AXES_FOR_FIGURE,
    nrow(candidate_axes)
  )
)

top_axes_figure[
  ,
  axis_label :=
    paste0(
      nichenet_ligand,
      " - ",
      receptor,
      " -> ",
      receiver
    )
]

top_axes_figure[
  ,
  axis_label := factor(
    axis_label,
    levels = rev(
      unique(axis_label)
    )
  )
]

write_csv_safe(
  top_axes_figure,
  file.path(
    DIRS$source,
    "Fig6A_top_candidate_axes_source.csv"
  )
)

plot_axes <- ggplot2::ggplot(
  top_axes_figure,
  ggplot2::aes(
    x =
      -final_axis_score,
    y =
      axis_label,
    fill =
      tf_symbol,
    size =
      pmax(
        aupr_corrected,
        0
      )
  )
) +
  ggplot2::geom_point(
    shape = 21,
    color = "black"
  ) +
  ggplot2::labs(
    title =
      "Candidate-TF-dependent macrophage-to-vascular/stromal axes",
    subtitle =
      "Official NicheNet-v2 prior; larger points indicate higher corrected AUPR",
    x =
      "Higher integrated axis priority",
    y = NULL,
    fill = "Candidate TF",
    size =
      "Corrected AUPR"
  ) +
  ggplot2::theme_bw(
    base_size = 9
  ) +
  ggplot2::theme(
    legend.position = "right"
  )

save_plot_bundle(
  plot_axes,
  "Fig6A_top_candidate_TF_ligand_receptor_axes",
  11,
  9
)

activity_plot_data <- merge(
  ligand_activity,
  unique(
    candidate_ligands[
      ,
      .(
        tf_symbol,
        ligand_matrix_symbol
      )
    ]
  ),
  by = "ligand_matrix_symbol",
  allow.cartesian = TRUE
)

activity_plot_data[
  ,
  activity_column :=
    paste(
      receiver,
      receiver_direction,
      sep = "__"
    )
]

activity_plot_data[
  ,
  ligand_TF :=
    paste(
      tf_symbol,
      ligand_matrix_symbol,
      sep = ":"
    )
]

activity_plot_data_all <- data.table::copy(
  activity_plot_data
)

activity_row_priority <- activity_plot_data_all[
  ,
  .(
    maximum_corrected_AUPR =
      safe_max(
        aupr_corrected,
        0
      ),
    supported_receiver_directions =
      data.table::uniqueN(
        activity_column
      )
  ),
  by = ligand_TF
]

activity_row_priority[
  !is.finite(
    maximum_corrected_AUPR
  ),
  maximum_corrected_AUPR := 0
]

data.table::setorderv(
  activity_row_priority,
  cols = c(
    "maximum_corrected_AUPR",
    "supported_receiver_directions",
    "ligand_TF"
  ),
  order = c(
    -1L,
    -1L,
    1L
  ),
  na.last = TRUE
)

activity_ligand_rows_for_figure <- head(
  activity_row_priority$ligand_TF,
  min(
    MAX_ACTIVITY_ROWS_FOR_HEATMAP,
    nrow(
      activity_row_priority
    )
  )
)

activity_plot_data <- activity_plot_data_all[
  ligand_TF %in%
    activity_ligand_rows_for_figure
]

activity_matrix <- data.table::dcast(
  activity_plot_data,
  ligand_TF ~
    activity_column,
  value.var =
    "aupr_corrected",
  fill = 0,
  fun.aggregate = max
)

activity_matrix_values <- as.matrix(
  activity_matrix[
    ,
    setdiff(
      names(activity_matrix),
      "ligand_TF"
    ),
    with = FALSE
  ]
)

rownames(activity_matrix_values) <-
  activity_matrix$ligand_TF

write_csv_safe(
  activity_plot_data_all,
  file.path(
    DIRS$source,
    "Fig6B_NicheNet_activity_heatmap_ALL_source.csv"
  )
)

write_csv_safe(
  activity_plot_data,
  file.path(
    DIRS$source,
    "Fig6B_NicheNet_activity_heatmap_DISPLAY_source.csv"
  )
)

save_heatmap_bundle(
  activity_matrix_values,
  "Fig6B_candidate_ligand_NicheNet_activity_heatmap",
  11,
  max(
    6,
    0.28 *
      nrow(activity_matrix_values) +
      2
  ),
  title =
    "Corrected AUPR by receiver and disease-response direction"
)

receptor_plot_data_all <- unique(
  candidate_axes[
    ,
    .(
      axis_id,
      final_axis_rank,
      tf_symbol,
      nichenet_ligand,
      receptor,
      receiver,
      median_receptor_pct_HFpEF,
      median_receptor_log2_cpm_HFpEF,
      strict_cross_stage_support
    )
  ]
)

top_receptor_axis_ids <- head(
  candidate_axes[
    order(
      final_axis_rank
    ),
    axis_id
  ],
  min(
    MAX_RECEPTOR_AXES_FOR_FIGURE,
    nrow(candidate_axes)
  )
)

receptor_plot_data <- receptor_plot_data_all[
  axis_id %in%
    top_receptor_axis_ids
]

receptor_plot_data[
  ,
  axis_label :=
    paste(
      tf_symbol,
      nichenet_ligand,
      receptor,
      sep = ":"
    )
]

receptor_axis_order <- receptor_plot_data[
  order(
    final_axis_rank,
    -median_receptor_pct_HFpEF,
    axis_label
  ),
  unique(
    axis_label
  )
]

receptor_plot_data[
  ,
  axis_label := factor(
    axis_label,
    levels = rev(
      receptor_axis_order
    )
  )
]

write_csv_safe(
  receptor_plot_data_all,
  file.path(
    DIRS$source,
    "Fig6C_receiver_receptor_expression_ALL_source.csv"
  )
)

write_csv_safe(
  receptor_plot_data,
  file.path(
    DIRS$source,
    "Fig6C_receiver_receptor_expression_DISPLAY_source.csv"
  )
)

plot_receptor <- ggplot2::ggplot(
  receptor_plot_data,
  ggplot2::aes(
    x = receiver,
    y = axis_label,
    size =
      median_receptor_pct_HFpEF,
    fill =
      median_receptor_log2_cpm_HFpEF,
    shape =
      strict_cross_stage_support
  )
) +
  ggplot2::geom_point(
    color = "black"
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      `FALSE` = 21,
      `TRUE` = 24
    )
  ) +
  ggplot2::labs(
    title =
      "Receiver receptor expression support",
    subtitle =
      paste0(
        "Top ",
        data.table::uniqueN(
          receptor_plot_data$axis_id
        ),
        " integrated axes; complete axes retained in source data"
      ),
    x = NULL,
    y =
      "TF:ligand:receptor",
    size =
      "HFpEF receptor\nfraction",
    fill =
      "HFpEF receptor\nlog2 CPM",
    shape =
      "Strict cross-stage\nligand support"
  ) +
  ggplot2::theme_bw(
    base_size = 9
  ) +
  ggplot2::theme(
    axis.text.y =
      ggplot2::element_text(
        size = 7
      )
  )

save_plot_bundle(
  plot_receptor,
  "Fig6C_receiver_receptor_expression_support",
  10,
  max(
    7,
    0.25 *
      data.table::uniqueN(
        receptor_plot_data$
          axis_label
      ) +
      2
  )
)


target_axis_manifest <- candidate_axes[
  is.finite(
    target_links
  ) &
    target_links > 0
]

data.table::setorder(
  target_axis_manifest,
  final_axis_rank,
  axis_id
)

top_axis_count_for_targets <- min(
  10L,
  nrow(
    target_axis_manifest
  )
)

target_axis_manifest <- head(
  target_axis_manifest,
  top_axis_count_for_targets
)

target_heatmap_data <- merge(
  ligand_target_links,
  unique(
    target_axis_manifest[
      ,
      .(
        tf_symbol,
        ligand_matrix_symbol,
        receiver,
        predicted_receiver_direction
      )
    ]
  ),
  by.x = c(
    "ligand_matrix_symbol",
    "receiver",
    "receiver_direction"
  ),
  by.y = c(
    "ligand_matrix_symbol",
    "receiver",
    "predicted_receiver_direction"
  ),
  allow.cartesian = TRUE
)

if (nrow(target_heatmap_data) > 0L) {
  data.table::setorder(
    target_heatmap_data,
    tf_symbol,
    ligand_matrix_symbol,
    receiver,
    -regulatory_potential
  )

  target_heatmap_data <- target_heatmap_data[
    ,
    head(
      .SD,
      TOP_TARGETS_FOR_HEATMAP
    ),
    by = .(
      tf_symbol,
      ligand_matrix_symbol,
      receiver
    )
  ]

  target_heatmap_data[
    ,
    axis_column :=
      paste(
        tf_symbol,
        ligand_matrix_symbol,
        receiver,
        sep = "__"
      )
  ]

  target_heatmap_data_all <- data.table::copy(
    target_heatmap_data
  )

  target_row_priority <- target_heatmap_data_all[
    ,
    .(
      maximum_regulatory_potential =
        safe_max(
          regulatory_potential,
          0
        ),
      supported_axes =
        data.table::uniqueN(
          axis_column
        )
    ),
    by = feature_key
  ]

  target_row_priority[
    !is.finite(
      maximum_regulatory_potential
    ),
    maximum_regulatory_potential := 0
  ]

  data.table::setorderv(
    target_row_priority,
    cols = c(
      "maximum_regulatory_potential",
      "supported_axes",
      "feature_key"
    ),
    order = c(
      -1L,
      -1L,
      1L
    ),
    na.last = TRUE
  )

  target_features_for_figure <- head(
    target_row_priority$feature_key,
    min(
      MAX_TARGET_ROWS_FOR_HEATMAP,
      nrow(
        target_row_priority
      )
    )
  )

  target_heatmap_data <- target_heatmap_data_all[
    feature_key %in%
      target_features_for_figure
  ]

  target_heatmap_matrix <- data.table::dcast(
    target_heatmap_data,
    feature_key ~
      axis_column,
    value.var =
      "regulatory_potential",
    fill = 0,
    fun.aggregate = max
  )

  target_heatmap_values <- as.matrix(
    target_heatmap_matrix[
      ,
      setdiff(
        names(
          target_heatmap_matrix
        ),
        "feature_key"
      ),
      with = FALSE
    ]
  )

  rownames(target_heatmap_values) <-
    target_heatmap_matrix$
      feature_key
} else {
  target_heatmap_data_all <-
    data.table::data.table()

  target_heatmap_values <- matrix(
    numeric(),
    nrow = 0L,
    ncol = 0L
  )
}

write_csv_safe(
  target_heatmap_data_all,
  file.path(
    DIRS$source,
    "Fig6D_ligand_target_heatmap_ALL_source.csv"
  )
)

write_csv_safe(
  target_heatmap_data,
  file.path(
    DIRS$source,
    "Fig6D_ligand_target_heatmap_DISPLAY_source.csv"
  )
)

if (
  nrow(target_heatmap_values) >
    0L &&
  ncol(target_heatmap_values) >
    0L
) {
  save_heatmap_bundle(
    target_heatmap_values,
    "Fig6D_top_axis_ligand_target_regulatory_potential",
    11,
    max(
      7,
      0.20 *
        nrow(
          target_heatmap_values
        ) +
        2
    ),
    title =
      "Top NicheNet ligand-target regulatory potentials"
  )
}

candidate_summary_plot <- data.table::copy(
  candidate_summary
)

candidate_summary_plot[
  ,
  tf_symbol := factor(
    tf_symbol,
    levels =
      candidate_rank[
        order(
          final_robustness_rank
        ),
        tf_symbol
      ]
  )
]

write_csv_safe(
  candidate_summary_plot,
  file.path(
    DIRS$source,
    "Fig6E_candidate_summary_source.csv"
  )
)

plot_candidate_summary <- ggplot2::ggplot(
  candidate_summary_plot,
  ggplot2::aes(
    x = tf_symbol,
    y = total_axes,
    fill = tf_symbol
  )
) +
  ggplot2::geom_col(
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "strict=",
        strict_cross_stage_axes
      )
    ),
    vjust = -0.3,
    size = 3
  ) +
  ggplot2::labs(
    title =
      "Communication-axis coverage by candidate TF",
    subtitle =
      "Bhlhe40 is retained as a program-recovery contrast and is not required to produce a communication axis",
    x = NULL,
    y =
      "Expression- and NicheNet-supported axes"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  ) +
  ggplot2::expand_limits(
    y =
      max(
        candidate_summary_plot$
          total_axes,
        na.rm = TRUE
      ) *
      1.15 +
      1
  )

save_plot_bundle(
  plot_candidate_summary,
  "Fig6E_candidate_TF_communication_axis_coverage",
  8,
  6
)

############################################################
## 18. Figure export audit
############################################################

figure_export_audit <- data.table::rbindlist(
  figure_export_records,
  use.names = TRUE,
  fill = TRUE
)

expected_main_figure_stems <- c(
  "Fig6A_top_candidate_TF_ligand_receptor_axes",
  "Fig6B_candidate_ligand_NicheNet_activity_heatmap",
  "Fig6C_receiver_receptor_expression_support",
  "Fig6E_candidate_TF_communication_axis_coverage"
)

if (
  nrow(target_heatmap_values) > 0L &&
  ncol(target_heatmap_values) > 0L
) {
  expected_main_figure_stems <- append(
    expected_main_figure_stems,
    "Fig6D_top_axis_ligand_target_regulatory_potential",
    after = 3L
  )
}

figure_export_audit[
  ,
  expected_main_figure :=
    stem %in%
      expected_main_figure_stems
]

write_csv_safe(
  figure_export_audit,
  file.path(
    DIRS$tables,
    "20A_stage6_figure_export_audit.csv"
  )
)

if (
  nrow(figure_export_audit) <
    length(
      expected_main_figure_stems
    ) ||
  any(
    figure_export_audit$
      files_valid != TRUE
  ) ||
  !all(
    expected_main_figure_stems %in%
      figure_export_audit$stem
  )
) {
  stop(
    "One or more required Stage 6 figures were not exported and validated."
  )
}

############################################################
## 19. Workbook, methods, and parameters
############################################################

workbook_path <- file.path(
  DIRS$tables,
  "20_stage6_TF_dependent_communication_key_results.xlsx"
)

workbook_sheets <- list(
  Resources =
    as.data.frame(
      resource_audit
    ),
  Candidate_TFs =
    as.data.frame(
      candidate_rank
    ),
  Ligand_coverage =
    as.data.frame(
      candidate_ligand_coverage
    ),
  Ligand_support =
    as.data.frame(
      candidate_ligands
    ),
  Cell_counts =
    as.data.frame(
      cell_count_audit
    ),
  Receiver_gene_sets =
    as.data.frame(
      receiver_gene_set_summary
    ),
  NicheNet_activity =
    as.data.frame(
      ligand_activity
    ),
  Top_axes =
    as.data.frame(
      head(
        candidate_axes,
        min(
          500L,
          nrow(candidate_axes)
        )
      )
    ),
  Candidate_summary =
    as.data.frame(
      candidate_summary
    ),
  Rank_stability =
    as.data.frame(
      axis_ranking_stability
    ),
  Figure_audit =
    as.data.frame(
      figure_export_audit
    )
)

writexl::write_xlsx(
  workbook_sheets,
  workbook_path
)

xlsx_contents <- utils::unzip(
  workbook_path,
  list = TRUE
)

xlsx_required_files <- c(
  "[Content_Types].xml",
  "xl/workbook.xml",
  "xl/worksheets/sheet1.xml"
)

xlsx_structure_ok <- all(
  xlsx_required_files %in%
    xlsx_contents$Name
)

if (!xlsx_structure_ok) {
  stop(
    "Generated Stage 6 XLSX failed internal structure validation."
  )
}

parameter_table <- data.table::data.table(
  parameter = c(
    "Random seed",
    "Sender cell type",
    "Receiver cell types",
    "Stage 6 candidate TFs",
    "Primary Stage 5 perturbation method",
    "Primary Stage 5 perturbation mode",
    "Primary Stage 5 perturbation strength",
    "NicheNet record",
    "Ligand-target matrix MD5",
    "Ligand-receptor network MD5",
    "Minimum sender cells per sample",
    "Minimum receiver cells per sample",
    "Minimum expressed-cell fraction",
    "Minimum HFpEF samples expressed",
    "Receiver primary FDR",
    "Receiver minimum absolute logFC",
    "Minimum receiver target genes",
    "Fallback receiver target genes",
    "Maximum receiver target genes",
    "Top target links per ligand",
    "Inferential unit",
    "Cell-level expression role",
    "Sample ligand-receptor correlation role",
    "Maximum receptor axes displayed",
    "Maximum activity heatmap rows",
    "Maximum target heatmap rows",
    "Maximum ggplot dimension",
    "Maximum heatmap dimension",
    "Nfkb1 forced"
  ),
  value = c(
    "20260714",
    SENDER_CELL_TYPE,
    paste(
      RECEIVER_CELL_TYPES,
      collapse = "; "
    ),
    paste(
      CANDIDATE_TFS,
      collapse = "; "
    ),
    PRIMARY_PERTURBATION_METHOD,
    PRIMARY_PERTURBATION_MODE,
    as.character(
      PRIMARY_PERTURBATION_STRENGTH
    ),
    paste0(
      "Zenodo ",
      NICHENET_RECORD_ID,
      "; DOI 10.5281/zenodo.7074291"
    ),
    NICHENET_MATRIX_MD5,
    NICHENET_LR_MD5,
    as.character(
      MIN_SENDER_CELLS_PER_SAMPLE
    ),
    as.character(
      MIN_RECEIVER_CELLS_PER_SAMPLE
    ),
    as.character(
      MIN_EXPRESSED_CELL_FRACTION
    ),
    as.character(
      MIN_HFPEF_SAMPLES_EXPRESSED
    ),
    as.character(
      RECEIVER_PRIMARY_FDR
    ),
    as.character(
      RECEIVER_MIN_ABS_LOGFC
    ),
    as.character(
      MIN_RECEIVER_TARGET_GENES
    ),
    as.character(
      FALLBACK_RECEIVER_TARGET_GENES
    ),
    as.character(
      MAX_RECEIVER_TARGET_GENES
    ),
    as.character(
      TOP_TARGET_LINKS_PER_LIGAND
    ),
    "Biological sample",
    "Descriptive expression gate only",
    "Descriptive co-variation only; n=6",
    as.character(
      MAX_RECEPTOR_AXES_FOR_FIGURE
    ),
    as.character(
      MAX_ACTIVITY_ROWS_FOR_HEATMAP
    ),
    as.character(
      MAX_TARGET_ROWS_FOR_HEATMAP
    ),
    paste0(
      MAX_GGPLOT_WIDTH_IN,
      " x ",
      MAX_GGPLOT_HEIGHT_IN,
      " inches"
    ),
    paste0(
      MAX_HEATMAP_WIDTH_IN,
      " x ",
      MAX_HEATMAP_HEIGHT_IN,
      " inches"
    ),
    "FALSE"
  ),
  rationale = c(
    "Reproducibility",
    "Locked Stage 3 macrophage/monocyte annotation",
    "Eligible vascular and stromal populations present in all six samples",
    "Stage 5B prioritization: NF-kB communication module plus Bhlhe40 contrast",
    "Primary Stage 5 v2 formulation",
    "Moves HFpEF TF activity toward the Control reference",
    "Primary Stage 5 v2 strength",
    "Fixed official mouse NicheNet-v2 prior",
    "Verify the exact 191-MB mouse ligand-target matrix",
    "Verify the exact mouse ligand-receptor network",
    "Prevent unstable sender expression estimates",
    "Prevent unstable receiver expression estimates",
    "Require detectable transcript support",
    "Require reproducibility in at least two HFpEF samples",
    "Primary receiver target selection where available",
    "Avoid trivial receiver effects",
    "Minimum informative NicheNet gene set",
    "Fallback for small n=3 per group without claiming formal significance",
    "Limit broad nonspecific response sets",
    "Restrict target-link reports to the strongest prior links",
    "Avoid cell-level pseudoreplication",
    "Cell-level percentages do not generate inferential P values",
    "Six-sample correlations are not used as formal significance tests",
    "Keep the receptor dot plot readable while complete axes remain in source data",
    "Keep the ligand-activity heatmap readable while complete activity data remain in source data",
    "Keep the ligand-target heatmap readable while complete links remain in source data",
    "Prevent ggsave from receiving an invalid or excessive dimension",
    "Prevent excessive raster/vector heatmap dimensions",
    "No candidate is manually promoted in axis ranking"
  )
)

write_csv_safe(
  parameter_table,
  file.path(
    DIRS$methods,
    "stage6_parameters_and_rationale.csv"
  )
)

methods_text <- c(
  "HFpEF Stage 6 FINAL v3",
  "Candidate-TF-dependent macrophage-to-vascular/stromal communication analysis",
  "",
  "Input boundary:",
  "- Stage 6 used the completed Stage 3 annotated Seurat object and sample-level major-cell-type pseudobulk differential-expression table.",
  "- Stage 4 supplied the Top150 drug-opposed macrophage program directions and Stage 3 directional-support flags.",
  "- Stage 5 supplied gene-level predicted consequences of disease-normalizing TF-activity adjustment.",
  "- Stage 5B supplied the final candidate robustness ranks.",
  "- No raw single-cell processing, clustering, cell-type annotation, TF inference, or virtual perturbation was repeated.",
  "",
  "Candidates and cell populations:",
  "- Nfkb1 was the primary inflammation/communication candidate.",
  "- Rela was the NF-kB family sensitivity candidate.",
  "- Bhlhe40 was retained as a program-recovery contrast candidate and was not required to generate a communication axis.",
  "- Macrophage_Monocyte was the sender; Endothelial, Fibroblast, Pericyte, and Smooth_muscle were receivers.",
  "",
  "Prior resources:",
  "- The exact mouse NicheNet-v2 ligand-target matrix and ligand-receptor network from Zenodo record 7074291 were used.",
  "- Published MD5 hashes were verified before analysis.",
  "- The workflow did not use a live interaction API, decoupleR, OmnipathR, or an unversioned network.",
  "",
  "Candidate ligand gate:",
  "- Stage 5 weighted minimum-norm, disease-normalization, strength-1 gene-level predictions were intersected with NicheNet ligands.",
  "- Ligands had to be expressed in macrophages in at least two HFpEF samples and at a median cell fraction of at least 5%.",
  "- Strict Tier-A support additionally required the Stage 5 predicted normalization direction, Top150 Stage 2 drug-opposed direction, and Stage 3 macrophage pseudobulk direction to agree.",
  "",
  "Receiver gene sets:",
  "- Receiver disease-response genes came from the locked Stage 3 sample-level edgeR and limma-voom pseudobulk results.",
  "- The primary set required edgeR/limma sign agreement, absolute mean logFC of at least 0.10, and FDR <=0.10 in either method.",
  "- When fewer than 20 genes met that criterion, a prespecified exploratory sign-consistent ranked set was used and explicitly labeled.",
  "- Background genes were expressed in at least 5% of receiver cells.",
  "",
  "NicheNet activity:",
  "- For every candidate ligand and receiver direction, AUROC, AUPR, corrected AUPR, Pearson correlation, and Spearman correlation were computed from the official ligand-target regulatory-potential matrix.",
  "- Corrected AUPR was the primary ligand-activity metric; Pearson and AUROC were retained as sensitivity metrics.",
  "- Positive ligand-target links were restricted to the receiver disease-response gene set and ranked by regulatory potential.",
  "",
  "Ligand-receptor expression support:",
  "- Receptor components were conservatively required to be detectable; for receptor complexes, the minimum component expression was used.",
  "- A receptor had to be expressed in at least 5% of receiver cells and in at least two HFpEF samples.",
  "",
  "Ranking:",
  "- Axis priority used transparent rank aggregation across Stage 5B TF robustness, Stage 5 ligand-change magnitude, sender expression, receiver receptor expression, NicheNet corrected AUPR, NicheNet Pearson correlation, ligand-target support, and cross-stage direction support.",
  "- Sensitivity scenarios omitted the TF prior, emphasized NicheNet/expression, or emphasized cross-stage support.",
  "- Nfkb1 was not forced.",
  "",
  "Figure policy:",
  "- Complete communication-axis, receptor-expression, ligand-activity, and ligand-target tables were retained as source data.",
  "- Manuscript figures displayed prespecified top-ranked subsets to remain legible.",
  "- Plot dimensions were capped and every PNG/PDF/TIFF export was checked for existence and nonzero file size.",
  "",
  "Claim boundary:",
  "- The results prioritize candidate TF-ligand-receptor-target axes.",
  "- They do not establish physical ligand-receptor binding, direct TF control in vivo, or a causal effect of dapagliflozin.",
  "- Receiver pseudobulk statistics use biological samples; cell-level percentages and six-sample correlations are descriptive."
)

writeLines(
  methods_text,
  file.path(
    DIRS$methods,
    "stage6_methods_and_claim_boundaries.txt"
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
## 20. Completion checks and run status
############################################################

warnings_table <- if (
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
  warnings_table,
  file.path(
    DIRS$tables,
    "21_stage6_warnings_and_nonfatal_issues.csv"
  )
)

scientific_check_records <- list()

add_scientific_check <- function(
  check_name,
  observed,
  expected,
  comparison = "equal"
) {
  status <- if (
    comparison == "equal"
  ) {
    ifelse(
      observed == expected,
      "PASS",
      "FAIL"
    )
  } else if (
    comparison == "at_least"
  ) {
    ifelse(
      observed >= expected,
      "PASS",
      "FAIL"
    )
  } else if (
    comparison == "at_most"
  ) {
    ifelse(
      observed <= expected,
      "PASS",
      "FAIL"
    )
  } else {
    stop(
      "Unknown scientific-check comparison: ",
      comparison
    )
  }

  scientific_check_records[[length(scientific_check_records) + 1L]] <<-
    data.table::data.table(
      check = check_name,
      observed =
        as.numeric(observed),
      expected =
        as.numeric(expected),
      comparison =
        comparison,
      status =
        status
    )

  invisible(status)
}

add_scientific_check(
  "Stage 3 completed",
  as.integer(
    observed_statuses["Stage3"] ==
      expected_statuses["Stage3"]
  ),
  1
)

add_scientific_check(
  "Stage 4 completed",
  as.integer(
    observed_statuses["Stage4"] ==
      expected_statuses["Stage4"]
  ),
  1
)

add_scientific_check(
  "Stage 5 completed",
  as.integer(
    observed_statuses["Stage5"] ==
      expected_statuses["Stage5"]
  ),
  1
)

add_scientific_check(
  "Stage 5B completed",
  as.integer(
    observed_statuses["Stage5B"] ==
      expected_statuses["Stage5B"]
  ),
  1
)

add_scientific_check(
  "Upstream failed checks",
  sum(
    stage3_checks$status != "PASS"
  ) +
    sum(
      stage4_checks$status != "PASS"
    ) +
    sum(
      stage5_checks$status != "PASS"
    ) +
    sum(
      stage5b_checks$status != "PASS"
    ),
  0
)

add_scientific_check(
  "Biological samples",
  nrow(sample_meta),
  6
)

add_scientific_check(
  "Control samples",
  sum(
    sample_meta$condition ==
      "Control"
  ),
  3
)

add_scientific_check(
  "HFpEF samples",
  sum(
    sample_meta$condition ==
      "HFpEF"
  ),
  3
)

add_scientific_check(
  "Verified NicheNet resources",
  sum(
    resource_audit$valid ==
      TRUE
  ),
  2
)

add_scientific_check(
  "Sender and receiver sample-cell groups",
  sum(
    cell_count_audit$sufficient ==
      TRUE
  ),
  nrow(cell_count_audit)
)

add_scientific_check(
  "Stage 6 candidate TFs",
  nrow(candidate_rank),
  length(CANDIDATE_TFS)
)

add_scientific_check(
  "Candidates represented in ligand coverage",
  nrow(candidate_ligand_coverage),
  length(CANDIDATE_TFS)
)

add_scientific_check(
  "Candidates with at least one NicheNet ligand",
  sum(
    candidate_ligand_coverage$
      direct_NicheNet_LR_ligands >
      0L
  ),
  1,
  "at_least"
)

add_scientific_check(
  "Receiver gene sets",
  nrow(receiver_gene_set_summary),
  length(RECEIVER_CELL_TYPES) *
    2
)

add_scientific_check(
  "Receiver-direction sets with at least two genes",
  sum(
    receiver_gene_set_summary$
      selected_genes >=
      2
  ),
  1,
  "at_least"
)

add_scientific_check(
  "NicheNet ligand activities",
  nrow(ligand_activity),
  1,
  "at_least"
)

add_scientific_check(
  "Positive ligand-target links",
  nrow(ligand_target_links),
  1,
  "at_least"
)

add_scientific_check(
  "Candidate communication axes",
  nrow(candidate_axes),
  1,
  "at_least"
)

add_scientific_check(
  "Candidates represented in communication summary",
  nrow(candidate_summary),
  length(CANDIDATE_TFS)
)

add_scientific_check(
  "Candidates with at least one communication axis",
  sum(
    candidate_summary$total_axes >
      0L
  ),
  1,
  "at_least"
)

add_scientific_check(
  "Axis ranking scenarios",
  data.table::uniqueN(
    axis_ranking_sensitivity$
      scenario
  ),
  length(
    scenario_definitions
  )
)

add_scientific_check(
  "Scientific checkpoint",
  as.integer(
    file.exists(
      file.path(
        DIRS$objects,
        "CHECKPOINT_stage6_scientific_results_pre_figures.rds"
      )
    )
  ),
  1
)

add_scientific_check(
  "Workbook",
  as.integer(
    file.exists(
      workbook_path
    )
  ),
  1
)

add_scientific_check(
  "Workbook structure",
  as.integer(
    xlsx_structure_ok
  ),
  1
)

add_scientific_check(
  "Required main figures exported",
  sum(
    expected_main_figure_stems %in%
      figure_export_audit$stem
  ),
  length(
    expected_main_figure_stems
  )
)

add_scientific_check(
  "Figure export failures",
  sum(
    figure_export_audit$
      files_valid != TRUE
  ),
  0
)

add_scientific_check(
  "Figure dimension violations",
  sum(
    figure_export_audit$
      effective_width_in >
      ifelse(
        figure_export_audit$
          plot_type == "ggplot",
        MAX_GGPLOT_WIDTH_IN,
        MAX_HEATMAP_WIDTH_IN
      ) |
      figure_export_audit$
        effective_height_in >
      ifelse(
        figure_export_audit$
          plot_type == "ggplot",
        MAX_GGPLOT_HEIGHT_IN,
        MAX_HEATMAP_HEIGHT_IN
      )
  ),
  0
)

scientific_checks <- data.table::rbindlist(
  scientific_check_records,
  use.names = TRUE,
  fill = TRUE
)

write_csv_safe(
  scientific_checks,
  file.path(
    DIRS$tables,
    "22_stage6_scientific_completion_checks.csv"
  )
)

script_copy_status <-
  "NOT_DETECTED"

if (
  length(SCRIPT_FILE) == 1L &&
  !is.na(SCRIPT_FILE) &&
  file.exists(SCRIPT_FILE)
) {
  methods_script <- file.path(
    DIRS$methods,
    basename(
      EXPECTED_SCRIPT_FILE
    )
  )

  check_script <- file.path(
    DIRS$check,
    basename(
      EXPECTED_SCRIPT_FILE
    )
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
  all(
    scientific_checks$status ==
      "PASS"
  )
) {
  "COMPLETED_STAGE6_READY_FOR_REVIEW"
} else {
  "COMPLETED_STAGE6_REVIEW_REQUIRED"
}

top_axis <- if (
  nrow(candidate_axes) > 0L
) {
  candidate_axes$axis_id[1L]
} else {
  NA_character_
}

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
  sender_cell_type =
    SENDER_CELL_TYPE,
  receiver_cell_types =
    paste(
      RECEIVER_CELL_TYPES,
      collapse = ";"
    ),
  candidate_TFs =
    paste(
      CANDIDATE_TFS,
      collapse = ";"
    ),
  verified_NicheNet_resources =
    sum(
      resource_audit$valid ==
        TRUE
    ),
  candidate_ligands =
    data.table::uniqueN(
      candidate_ligands$
        target_key
    ),
  NicheNet_activity_rows =
    nrow(ligand_activity),
  ligand_target_links =
    nrow(ligand_target_links),
  communication_axes =
    nrow(candidate_axes),
  top_axis =
    top_axis,
  Nfkb1_axes =
    candidate_axes[
      tf_symbol == "Nfkb1",
      .N
    ],
  Rela_axes =
    candidate_axes[
      tf_symbol == "Rela",
      .N
    ],
  Bhlhe40_axes =
    candidate_axes[
      tf_symbol == "Bhlhe40",
      .N
    ],
  Nfkb1_forced = FALSE,
  warnings =
    nrow(warnings_table),
  script_copy_status =
    script_copy_status,
  scientific_checks_failed =
    sum(
      scientific_checks$status !=
        "PASS"
    ),
  overall_status =
    overall_status
)

write_csv_safe(
  run_status,
  file.path(
    DIRS$tables,
    "23_stage6_run_status.csv"
  )
)

readme <- c(
  "HFpEF Reanalysis Project - Stage 6 FINAL v3",
  "Candidate-TF-dependent macrophage-to-vascular/stromal communication analysis",
  "",
  paste0(
    "Overall status: ",
    overall_status
  ),
  paste0(
    "Verified NicheNet-v2 resources: ",
    sum(
      resource_audit$valid ==
        TRUE
    ),
    "/2"
  ),
  paste0(
    "Communication axes: ",
    nrow(candidate_axes)
  ),
  paste0(
    "Top axis: ",
    top_axis
  ),
  "",
  "Candidates:",
  "- Nfkb1: primary inflammation/communication candidate.",
  "- Rela: NF-kB family sensitivity candidate.",
  "- Bhlhe40: program-recovery contrast candidate.",
  "",
  "Primary interpretation:",
  "- The output prioritizes TF-sensitive macrophage ligands, expressed receiver receptors, and NicheNet-supported receiver targets.",
  "- It does not establish causal or physical signaling.",
  "",
  "Figure export:",
  "- Complete source tables are retained; manuscript figures use readable top-ranked subsets.",
  "- Every PNG/PDF/TIFF export is audited for existence, nonzero size, and bounded dimensions.",
  "",
  "Upload the Stage 6 CHECK package before proceeding to external validation."
)

writeLines(
  readme,
  file.path(
    OUT_DIR,
    "README_stage6.txt"
  ),
  useBytes = TRUE
)

############################################################
## 21. CHECK package and hashes
############################################################

review_files <- c(
  LOG_FILE,
  file.path(
    DIRS$tables,
    "00_stage6_nichenet_resource_audit.csv"
  ),
  file.path(
    DIRS$tables,
    "01_stage6_candidate_TF_manifest.csv"
  ),
  file.path(
    DIRS$tables,
    "03_stage6_candidate_ligand_coverage.csv"
  ),
  file.path(
    DIRS$tables,
    "04_stage6_sender_receiver_cell_count_audit.csv"
  ),
  file.path(
    DIRS$tables,
    "07_stage6_macrophage_ligand_expression_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "08_stage6_receiver_receptor_expression_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "09_stage6_candidate_ligand_cross_stage_support.csv"
  ),
  file.path(
    DIRS$tables,
    "11_stage6_receiver_gene_set_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "12_stage6_nichenet_ligand_activity.csv"
  ),
  file.path(
    DIRS$tables,
    "14_stage6_nichenet_ligand_target_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "16_stage6_top_candidate_axes.csv"
  ),
  file.path(
    DIRS$tables,
    "18_stage6_axis_ranking_stability_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "19_stage6_candidate_TF_communication_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "20A_stage6_figure_export_audit.csv"
  ),
  workbook_path,
  file.path(
    DIRS$tables,
    "21_stage6_warnings_and_nonfatal_issues.csv"
  ),
  file.path(
    DIRS$tables,
    "22_stage6_scientific_completion_checks.csv"
  ),
  file.path(
    DIRS$tables,
    "23_stage6_run_status.csv"
  ),
  file.path(
    DIRS$methods,
    "stage6_parameters_and_rationale.csv"
  ),
  file.path(
    DIRS$methods,
    "stage6_methods_and_claim_boundaries.txt"
  ),
  file.path(
    DIRS$methods,
    "sessionInfo.txt"
  ),
  file.path(
    OUT_DIR,
    "README_stage6.txt"
  ),
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

for (
  file_i in review_files
) {
  target_file <- file.path(
    DIRS$check,
    basename(file_i)
  )

  if (
    normalizePath(
      file_i,
      winslash = "/",
      mustWork = FALSE
    ) !=
      normalizePath(
        target_file,
        winslash = "/",
        mustWork = FALSE
      )
  ) {
    file.copy(
      file_i,
      target_file,
      overwrite = TRUE
    )
  }
}

check_files <- list.files(
  DIRS$check,
  full.names = TRUE
)

check_manifest <- data.table::data.table(
  filename =
    basename(check_files),
  size_bytes =
    as.numeric(
      file.info(
        check_files
      )$size
    )
)

check_manifest[
  ,
  sha256 := vapply(
    check_files,
    function(file_i) {
      digest::digest(
        file = file_i,
        algo = "sha256",
        serialize = FALSE
      )
    },
    character(1)
  )
]

write_csv_safe(
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

log_msg(
  "Stage 6 analysis finished."
)

log_msg(
  "Overall status: ",
  overall_status
)

log_msg(
  "Communication axes: ",
  nrow(candidate_axes)
)

log_msg(
  "Top axis: ",
  top_axis
)

log_msg(
  "CHECK package: ",
  CHECK_ZIP
)

cat(
  "\n============================================================\n"
)

cat(
  "HFpEF Stage 6 communication analysis completed\n"
)

cat(
  "Status: ",
  overall_status,
  "\n",
  sep = ""
)

cat(
  "Output: ",
  OUT_DIR,
  "\n",
  sep = ""
)

cat(
  "CHECK: ",
  CHECK_ZIP,
  "\n",
  sep = ""
)

cat(
  "Communication axes: ",
  nrow(candidate_axes),
  "\n",
  sep = ""
)

cat(
  "Top axis: ",
  top_axis,
  "\n",
  sep = ""
)

cat(
  "Nfkb1 was not forced.\n"
)

cat(
  "Upload the CHECK package before the next stage.\n"
)

cat(
  "============================================================\n"
)
