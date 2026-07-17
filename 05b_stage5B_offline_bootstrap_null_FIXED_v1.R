############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 5B OFFLINE FIXED v1
## DoRothEA regulon bootstrap and unbiased matched-TF null
##
## Project:
##   <HFPEF_PROJECT_DIR>
##
## Required completed inputs:
##   Stage 4:
##     04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1
##   Stage 5:
##     05_stage5_multiTF_virtual_perturbation_FIXED_v2
##
## Prespecified candidate TFs:
##   Bhlhe40, Runx1, Spi1, Rel, Nfkb1, Rela
##
## Stage 5B OFFLINE objectives:
##   1) Recalculate the six full-regulon DoRothEA perturbations and
##      verify exact consistency with the completed Stage 5 v2 output.
##   2) Bootstrap 80% of each candidate regulon for 500 replicates to
##      quantify target-set dependence and rank stability.
##   3) Construct 1,000-draw matched random-TF null distributions for
##      each candidate without using Stage 2/3 program overlap to
##      select null TFs.
##   4) Produce empirical candidate percentiles and one-sided empirical
##      P values for prioritization before Stage 6 communication analysis.
##
## This version is deliberately offline:
##   - no CollecTRI;
##   - no decoupleR;
##   - no OmnipathR;
##   - no run-time internet access;
##   - no reprocessing of Stage 1-5.
##
## Interpretation boundary:
##   - This is a computational robustness analysis.
##   - It is not experimental knockdown, knockout, or inhibition.
##   - Bootstrap intervals quantify dependence on regulon targets, not
##     population-level biological confidence intervals.
##   - Empirical null P values are prioritization statistics, not causal
##     significance tests.
##   - Biological samples, not cells, are the inferential units.
##
## Output:
##   <HFPEF_PROJECT_DIR>/
##   05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1
##
## CHECK:
##   <HFPEF_PROJECT_DIR>/
##   05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1_CHECK.zip
##
## Recommended run:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage5B_OFFLINE_Bootstrap_Null_FIXED_v1.R",
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
## 0. Locked paths and analysis settings
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

STAGE3_SAMPLE_META_FILE <- file.path(
  STAGE3_DIR,
  "01_tables",
  "01_locked_GSE236585_sample_metadata.csv"
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

STAGE5_CANDIDATE_FILE <- file.path(
  STAGE5_DIR,
  "01_tables",
  "13_stage5_candidate_TF_rank_aggregation.csv"
)

STAGE_NAME <-
  "05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1"

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

EXPECTED_SCRIPT_FILE <- file.path(PROJECT_DIR, "05b_stage5B_offline_bootstrap_null_FIXED_v1.R")

REPLACE_EXISTING_STAGE5B <- TRUE

CANDIDATE_TFS_REQUESTED <- c(
  "Bhlhe40",
  "Runx1",
  "Spi1",
  "Rel",
  "Nfkb1",
  "Rela"
)

PRIMARY_SIGNATURE_SIZE <- 150L
MIN_ABS_OBSERVED_PROGRAM_GAP <- 0.10
MIN_TARGETS_PER_REGULON <- 10L
MAX_ABS_GENE_SHIFT_SD <- 2.50

BOOTSTRAP_REPLICATES <- 500L
BOOTSTRAP_TARGET_FRACTION <- 0.80

NULL_DRAWS_PER_CANDIDATE <- 1000L
NULL_NEAREST_POOL_SIZE <- 60L
NULL_MIN_POOL_SIZE <- 20L

STAGE5_AUDIT_TOLERANCE <- 1e-6

############################################################
## 1. Preflight and output setup
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
      candidates <- c(
        candidates,
        ofile
      )
    }
  }

  args <- commandArgs(
    trailingOnly = FALSE
  )

  file_arg <- grep(
    "^--file=",
    args,
    value = TRUE
  )

  if (length(file_arg) > 0L) {
    candidates <- c(
      candidates,
      sub(
        "^--file=",
        "",
        file_arg[1L]
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
  STAGE3_SAMPLE_META_FILE,
  STAGE4_STATUS_FILE,
  STAGE4_CHECKS_FILE,
  STAGE4_NETWORK_FILE,
  STAGE4_REGULON_SIZE_FILE,
  STAGE4_PROGRAM_FILE,
  STAGE4_PSEUDOBULK_RDS,
  STAGE4_ACTIVITY_RDS,
  STAGE5_STATUS_FILE,
  STAGE5_CHECKS_FILE,
  STAGE5_CANDIDATE_FILE
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs) > 0L) {
  stop(
    "Required Stage 3/4/5 input path(s) are missing:\n",
    paste(
      missing_inputs,
      collapse = "\n"
    )
  )
}

stage4_status <- data.table::fread(
  STAGE4_STATUS_FILE,
  encoding = "UTF-8"
)

stage5_status <- data.table::fread(
  STAGE5_STATUS_FILE,
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
    "Stage 4 is not in an allowed completed state."
  )
}

if (
  !"overall_status" %in% names(stage5_status) ||
  stage5_status$overall_status[1L] !=
    "COMPLETED_STAGE5_READY_FOR_REVIEW"
) {
  stop(
    "Stage 5 v2 is not locked as ",
    "COMPLETED_STAGE5_READY_FOR_REVIEW."
  )
}

if (
  !all(
    c(
      "check",
      "status"
    ) %in% names(stage4_checks)
  ) ||
  any(stage4_checks$status != "PASS")
) {
  stop(
    "At least one Stage 4 completion check is not PASS."
  )
}

if (
  !all(
    c(
      "check",
      "status"
    ) %in% names(stage5_checks)
  ) ||
  any(stage5_checks$status != "PASS")
) {
  stop(
    "At least one Stage 5 completion check is not PASS."
  )
}

replacement_audit <- data.table::data.table(
  path = c(
    OUT_DIR,
    CHECK_ZIP
  ),
  path_type = c(
    "stage5B_output_directory",
    "stage5B_check_zip"
  ),
  existed_before = FALSE,
  deletion_attempted = FALSE,
  deletion_succeeded = FALSE
)

if (REPLACE_EXISTING_STAGE5B) {
  for (
    i in seq_len(
      nrow(replacement_audit)
    )
  ) {
    target_path <- replacement_audit$path[i]

    existed <- (
      dir.exists(target_path) ||
        file.exists(target_path)
    )

    replacement_audit$existed_before[i] <- existed

    if (existed) {
      replacement_audit$deletion_attempted[i] <- TRUE

      unlink(
        target_path,
        recursive = dir.exists(target_path),
        force = TRUE
      )

      replacement_audit$deletion_succeeded[i] <- !(
        dir.exists(target_path) ||
          file.exists(target_path)
      )

      if (
        !replacement_audit$
          deletion_succeeded[i]
      ) {
        stop(
          "Failed to remove previous Stage 5B path: ",
          target_path
        )
      }
    } else {
      replacement_audit$deletion_succeeded[i] <- TRUE
    }
  }
} else if (
  dir.exists(OUT_DIR) ||
  file.exists(CHECK_ZIP)
) {
  stop(
    "Existing Stage 5B output detected while replacement is disabled."
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
  "stage5B_OFFLINE_bootstrap_null.log"
)

WARN_FILE <- file.path(
  DIRS$logs,
  "stage5B_OFFLINE_warnings.log"
)

data.table::fwrite(
  replacement_audit,
  file.path(
    DIRS$logs,
    "stage5B_replacement_audit.csv"
  )
)

log_msg <- function(
  ...,
  level = "INFO"
) {
  message_text <- paste0(
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
    message_text
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
  "Stage 5B OFFLINE robustness analysis started."
)

log_msg(
  "Stage 4 status: ",
  stage4_status$overall_status[1L]
)

log_msg(
  "Stage 5 status: ",
  stage5_status$overall_status[1L]
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
    "data.table",
    "Matrix",
    "ggplot2",
    "pheatmap",
    "writexl",
    "zip",
    "digest"
  ),
  required = TRUE
)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(pheatmap)
})

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
        note =
          "No records generated."
      ),
      path
    )
  } else {
    data.table::fwrite(
      x,
      path,
      compress = if (compress) {
        "gzip"
      } else {
        "none"
      }
    )
  }
}

save_plot_bundle <- function(
  plot_object,
  stem,
  width,
  height
) {
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
  matrix_object,
  stem,
  width,
  height,
  title = NULL
) {
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

  grDevices::png(
    paths["png"],
    width = width * 300,
    height = height * 300,
    res = 300
  )

  pheatmap::pheatmap(
    matrix_object,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    main = title
  )

  grDevices::dev.off()

  grDevices::pdf(
    paths["pdf"],
    width = width,
    height = height
  )

  pheatmap::pheatmap(
    matrix_object,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    main = title
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
    matrix_object,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    main = title
  )

  grDevices::dev.off()

  invisible(paths)
}

scale_rows <- function(matrix_object) {
  matrix_object <- as.matrix(
    matrix_object
  )

  row_means <- rowMeans(
    matrix_object,
    na.rm = TRUE
  )

  row_standard_deviations <- apply(
    matrix_object,
    1L,
    stats::sd,
    na.rm = TRUE
  )

  row_standard_deviations[
    !is.finite(
      row_standard_deviations
    ) |
      row_standard_deviations == 0
  ] <- 1

  output <- sweep(
    matrix_object,
    1L,
    row_means,
    "-"
  )

  output <- sweep(
    output,
    1L,
    row_standard_deviations,
    "/"
  )

  output[!is.finite(output)] <- 0

  output
}

make_feature_map <- function(features) {
  feature_table <- data.table::data.table(
    feature = as.character(features),
    feature_key = gene_key(features)
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

############################################################
## 4. Program and perturbation functions
############################################################

make_program_definitions <- function(
  stage4_manifest,
  feature_map
) {
  manifest <- data.table::copy(
    stage4_manifest
  )

  required_columns <- c(
    "program_name",
    "subset_name",
    "direction",
    "signature_size",
    "symbol_key"
  )

  missing_columns <- setdiff(
    required_columns,
    names(manifest)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Stage 4 program manifest is missing column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  manifest[
    ,
    support_class := ifelse(
      grepl(
        "Stage3Supported$",
        program_name
      ),
      "Stage3_supported",
      "Full_Stage2"
    )
  ]

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

  for (
    group_index in seq_len(
      nrow(grouping)
    )
  ) {
    subset_i <-
      grouping$subset_name[
        group_index
      ]

    signature_size_i <-
      grouping$signature_size[
        group_index
      ]

    support_class_i <-
      grouping$support_class[
        group_index
      ]

    part <- manifest[
      subset_name == subset_i &
        signature_size ==
          signature_size_i &
        support_class ==
          support_class_i
    ]

    up_map <- map_symbols_to_features(
      part[
        direction ==
          "Disease_up_Drug_down",
        symbol_key
      ],
      feature_map
    )

    down_map <- map_symbols_to_features(
      part[
        direction ==
          "Disease_down_Drug_up",
        symbol_key
      ],
      feature_map
    )

    up_features <- unique(
      up_map[
        !is.na(feature),
        feature
      ]
    )

    down_features <- unique(
      down_map[
        !is.na(feature),
        feature
      ]
    )

    if (
      length(up_features) +
        length(down_features) <
        5L
    ) {
      next
    }

    program_id <- paste(
      "DrugOpposedNet",
      subset_i,
      support_class_i,
      paste0(
        "Top",
        signature_size_i
      ),
      sep = "__"
    )

    definitions[[program_id]] <- list(
      program_id = program_id,
      program_category =
        "Stage2_drug_opposed_net",
      subset_name = subset_i,
      signature_size =
        as.integer(
          signature_size_i
        ),
      support_class =
        support_class_i,
      up_features = up_features,
      down_features = down_features,
      primary_program = (
        as.integer(
          signature_size_i
        ) ==
          PRIMARY_SIGNATURE_SIZE
      )
    )

    summary_records[[length(summary_records) + 1L]] <- data.table::data.table(
      program_id = program_id,
      program_category =
        "Stage2_drug_opposed_net",
      subset_name = subset_i,
      signature_size =
        as.integer(
          signature_size_i
        ),
      support_class =
        support_class_i,
      detected_up_genes =
        length(up_features),
      detected_down_genes =
        length(down_features),
      primary_program = (
        as.integer(
          signature_size_i
        ) ==
          PRIMARY_SIGNATURE_SIZE
      )
    )
  }

  functional_sets <- list(
    Inflammatory_Il1b = c(
      "Il1b",
      "Tnf",
      "Nfkbia",
      "Ccl2",
      "Ccl3",
      "Ccl4",
      "Cxcl2",
      "S100a8",
      "S100a9",
      "Ptgs2"
    ),
    NFkB_TNF_response = c(
      "Tnf",
      "Nfkbia",
      "Nfkbiz",
      "Rel",
      "Rela",
      "Nfkb1",
      "Tnfaip3",
      "Icam1",
      "Ccl2",
      "Ccl3",
      "Il1b"
    ),
    Inflammasome_pyroptosis = c(
      "Nlrp3",
      "Pycard",
      "Casp1",
      "Gsdmd",
      "Il1b",
      "Il18",
      "Txnip",
      "P2rx7"
    ),
    Interferon_response = c(
      "Isg15",
      "Ifit1",
      "Ifit2",
      "Ifit3",
      "Irf7",
      "Rsad2",
      "Oas1a",
      "Stat1",
      "Cxcl9",
      "Cxcl10"
    ),
    Antigen_presentation = c(
      "H2-Ab1",
      "H2-Aa",
      "Cd74",
      "Ciita",
      "H2-Eb1",
      "Tap1",
      "B2m"
    ),
    Spp1_Trem2_remodeling = c(
      "Spp1",
      "Trem2",
      "Gpnmb",
      "Fabp5",
      "Lpl",
      "Apoe",
      "Ctsb",
      "Ctsd",
      "Lgals3"
    ),
    Resident_Timd4_Lyve1 = c(
      "Timd4",
      "Lyve1",
      "Folr2",
      "Mrc1",
      "Cd163",
      "Vsig4",
      "C1qa",
      "C1qb",
      "C1qc"
    ),
    Ccr2_monocyte_like = c(
      "Ccr2",
      "Ly6c2",
      "Plac8",
      "Chil3",
      "Ctss",
      "Lgals3"
    ),
    Myeloid_identity = c(
      "Spi1",
      "Csf1r",
      "Lyz2",
      "Aif1",
      "Tyrobp",
      "Fcerg",
      "Ctss",
      "Laptm5",
      "C1qa",
      "C1qb",
      "C1qc"
    ),
    Lipid_cholesterol = c(
      "Apoe",
      "Lpl",
      "Abca1",
      "Abcg1",
      "Pparg",
      "Nr1h3",
      "Soat1",
      "Lipa",
      "Fabp5"
    ),
    Oxidative_stress = c(
      "Nfe2l2",
      "Hmox1",
      "Nqo1",
      "Gclc",
      "Gclm",
      "Sod2",
      "Txnrd1",
      "Prdx1"
    ),
    Cycling = c(
      "Mki67",
      "Top2a",
      "Stmn1",
      "Tubb5",
      "Hmgb2"
    )
  )

  for (
    set_name in names(
      functional_sets
    )
  ) {
    mapped <- map_symbols_to_features(
      functional_sets[[set_name]],
      feature_map
    )

    features <- unique(
      mapped[
        !is.na(feature),
        feature
      ]
    )

    if (length(features) < 3L) {
      next
    }

    program_id <- paste0(
      "Functional__",
      set_name
    )

    definitions[[program_id]] <- list(
      program_id = program_id,
      program_category =
        "Functional_state",
      subset_name = set_name,
      signature_size = NA_integer_,
      support_class = "Curated",
      up_features = features,
      down_features = character(),
      primary_program = set_name %in%
        c(
          "Inflammatory_Il1b",
          "NFkB_TNF_response",
          "Inflammasome_pyroptosis",
          "Interferon_response",
          "Spp1_Trem2_remodeling",
          "Myeloid_identity"
        )
    )

    summary_records[[length(summary_records) + 1L]] <- data.table::data.table(
      program_id = program_id,
      program_category =
        "Functional_state",
      subset_name = set_name,
      signature_size = NA_integer_,
      support_class = "Curated",
      detected_up_genes =
        length(features),
      detected_down_genes = 0L,
      primary_program = set_name %in%
        c(
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

score_programs <- function(
  z_matrix,
  definitions
) {
  output <- matrix(
    NA_real_,
    nrow = length(definitions),
    ncol = ncol(z_matrix),
    dimnames = list(
      names(definitions),
      colnames(z_matrix)
    )
  )

  for (
    program_id in names(
      definitions
    )
  ) {
    definition <- definitions[[program_id]]

    up_features <- intersect(
      definition$up_features,
      rownames(z_matrix)
    )

    down_features <- intersect(
      definition$down_features,
      rownames(z_matrix)
    )

    up_score <- if (
      length(up_features) > 0L
    ) {
      colMeans(
        z_matrix[
          up_features,
          ,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    } else {
      rep(
        0,
        ncol(z_matrix)
      )
    }

    down_score <- if (
      length(down_features) > 0L
    ) {
      colMeans(
        z_matrix[
          down_features,
          ,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    } else {
      rep(
        0,
        ncol(z_matrix)
      )
    }

    output[program_id, ] <-
      up_score - down_score
  }

  output
}

weighted_regulon_activity <- function(
  expression_matrix,
  network_subset
) {
  target_features <- intersect(
    unique(
      network_subset$target_feature
    ),
    rownames(expression_matrix)
  )

  if (
    length(target_features) <
      MIN_TARGETS_PER_REGULON
  ) {
    return(NULL)
  }

  scaled_expression <- scale_rows(
    expression_matrix[
      target_features,
      ,
      drop = FALSE
    ]
  )

  source_symbols <- sort(
    unique(
      network_subset$source_symbol
    )
  )

  output <- matrix(
    NA_real_,
    nrow = length(source_symbols),
    ncol = ncol(
      scaled_expression
    ),
    dimnames = list(
      source_symbols,
      colnames(
        scaled_expression
      )
    )
  )

  for (
    source_i in source_symbols
  ) {
    network_i <- data.table::copy(
      network_subset[
        source_symbol == source_i &
          target_feature %in%
            rownames(
              scaled_expression
            )
      ]
    )

    data.table::setorder(
      network_i,
      target_feature,
      -weight
    )

    network_i <- network_i[
      ,
      .SD[1L],
      by = target_feature
    ]

    if (
      nrow(network_i) <
        MIN_TARGETS_PER_REGULON
    ) {
      next
    }

    signed_weights <-
      as.numeric(
        network_i$mor
      ) *
      as.numeric(
        network_i$weight
      )

    signed_weights[
      !is.finite(signed_weights)
    ] <- 0

    denominator <- sum(
      abs(signed_weights)
    )

    if (
      !is.finite(denominator) ||
      denominator <= 0
    ) {
      next
    }

    target_z <- scaled_expression[
      network_i$target_feature,
      ,
      drop = FALSE
    ]

    output[source_i, ] <- colSums(
      sweep(
        target_z,
        1L,
        signed_weights,
        "*"
      )
    ) / denominator
  }

  output[
    rowSums(
      is.finite(output)
    ) > 0L,
    ,
    drop = FALSE
  ]
}

simulate_normalization <- function(
  observed_z,
  activity_matrix,
  tf_requested,
  network_table,
  sample_meta
) {
  if (
    !tf_requested %in%
      rownames(activity_matrix)
  ) {
    stop(
      "TF activity is unavailable: ",
      tf_requested
    )
  }

  sample_ids <- colnames(
    observed_z
  )

  metadata <- data.table::copy(
    sample_meta
  )

  metadata <- metadata[
    match(
      sample_ids,
      sample_accession
    )
  ]

  if (
    any(
      is.na(
        metadata$sample_accession
      )
    ) ||
    any(
      metadata$sample_accession !=
        sample_ids
    )
  ) {
    stop(
      "Sample metadata could not be aligned."
    )
  }

  network_tf <- data.table::copy(
    network_table[
      source_symbol ==
        tf_requested &
        target_feature %in%
          rownames(observed_z)
    ]
  )

  data.table::setorder(
    network_tf,
    target_feature,
    -weight
  )

  network_tf <- network_tf[
    ,
    .SD[1L],
    by = target_feature
  ]

  if (
    nrow(network_tf) <
      MIN_TARGETS_PER_REGULON
  ) {
    stop(
      "Too few targets for ",
      tf_requested,
      ": ",
      nrow(network_tf)
    )
  }

  observed_activity <- as.numeric(
    activity_matrix[
      tf_requested,
      sample_ids
    ]
  )

  names(observed_activity) <-
    sample_ids

  control_reference <- mean(
    observed_activity[
      metadata$condition ==
        "Control"
    ],
    na.rm = TRUE
  )

  desired_activity <-
    observed_activity

  hfpef_index <- which(
    metadata$condition ==
      "HFpEF"
  )

  desired_activity[hfpef_index] <-
    control_reference

  delta_activity <-
    desired_activity -
    observed_activity

  signed_weights <-
    as.numeric(
      network_tf$mor
    ) *
    as.numeric(
      network_tf$weight
    )

  signed_weights[
    !is.finite(signed_weights)
  ] <- 0

  denominator_absolute <- sum(
    abs(signed_weights)
  )

  denominator_squared <- sum(
    signed_weights^2
  )

  if (
    !is.finite(
      denominator_absolute
    ) ||
    denominator_absolute <= 0 ||
    !is.finite(
      denominator_squared
    ) ||
    denominator_squared <= 0
  ) {
    stop(
      "Invalid target weights for ",
      tf_requested
    )
  }

  gene_coefficients <-
    signed_weights *
    denominator_absolute /
    denominator_squared

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

  rownames(target_delta) <-
    network_tf$target_feature

  colnames(target_delta) <-
    sample_ids

  perturbed_z <- observed_z

  perturbed_z[
    network_tf$target_feature,
    sample_ids
  ] <- perturbed_z[
    network_tf$target_feature,
    sample_ids,
    drop = FALSE
  ] + target_delta

  list(
    perturbed_z = perturbed_z,
    target_count =
      nrow(network_tf),
    global_rms_shift_HFpEF =
      sqrt(
        mean(
          target_delta[
            ,
            metadata$condition ==
              "HFpEF",
            drop = FALSE
          ]^2,
          na.rm = TRUE
        )
      )
  )
}

summarize_perturbation <- function(
  observed_scores,
  perturbed_scores,
  sample_meta,
  program_summary
) {
  metadata <- data.table::copy(
    sample_meta
  )

  metadata <- metadata[
    match(
      colnames(
        observed_scores
      ),
      sample_accession
    )
  ]

  control_index <- which(
    metadata$condition ==
      "Control"
  )

  hfpef_index <- which(
    metadata$condition ==
      "HFpEF"
  )

  records <- lapply(
    rownames(observed_scores),
    function(program_id) {
      control_values <- as.numeric(
        observed_scores[
          program_id,
          control_index
        ]
      )

      observed_hfpef <- as.numeric(
        observed_scores[
          program_id,
          hfpef_index
        ]
      )

      perturbed_hfpef <- as.numeric(
        perturbed_scores[
          program_id,
          hfpef_index
        ]
      )

      control_mean <- mean(
        control_values,
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

      observed_gap <-
        observed_hfpef_mean -
        control_mean

      perturbed_gap <-
        perturbed_hfpef_mean -
        control_mean

      gap_reduction <-
        abs(observed_gap) -
        abs(perturbed_gap)

      data.table::data.table(
        program_id = program_id,
        observed_control_mean =
          control_mean,
        observed_hfpef_mean =
          observed_hfpef_mean,
        perturbed_hfpef_mean =
          perturbed_hfpef_mean,
        observed_gap =
          observed_gap,
        perturbed_gap =
          perturbed_gap,
        absolute_gap_reduction =
          gap_reduction,
        recovery_fraction = if (
          is.finite(observed_gap) &&
          abs(observed_gap) >=
            MIN_ABS_OBSERVED_PROGRAM_GAP
        ) {
          gap_reduction /
            abs(observed_gap)
        } else {
          NA_real_
        },
        observed_gap_eligible = (
          is.finite(observed_gap) &&
            abs(observed_gap) >=
              MIN_ABS_OBSERVED_PROGRAM_GAP
        ),
        sample_improvement_fraction =
          mean(
            abs(
              perturbed_hfpef -
                control_mean
            ) <
              abs(
                observed_hfpef -
                  control_mean
              ),
            na.rm = TRUE
          )
      )
    }
  )

  detail <- data.table::rbindlist(
    records,
    use.names = TRUE,
    fill = TRUE
  )

  merge(
    detail,
    program_summary,
    by = "program_id",
    all.x = TRUE
  )
}

result_summary <- function(
  detail_table
) {
  primary_stage2 <- detail_table[
    program_category ==
      "Stage2_drug_opposed_net" &
      primary_program == TRUE &
      observed_gap_eligible == TRUE
  ]

  all_stage2 <- detail_table[
    program_category ==
      "Stage2_drug_opposed_net" &
      observed_gap_eligible == TRUE
  ]

  supported_primary <- primary_stage2[
    support_class ==
      "Stage3_supported"
  ]

  inflammation <- detail_table[
    program_category ==
      "Functional_state" &
      subset_name %in%
        c(
          "Inflammatory_Il1b",
          "NFkB_TNF_response",
          "Inflammasome_pyroptosis",
          "Interferon_response"
        ) &
      observed_gap_eligible == TRUE
  ]

  myeloid <- detail_table[
    program_category ==
      "Functional_state" &
      subset_name ==
        "Myeloid_identity"
  ]

  data.table::data.table(
    primary_programs_evaluated =
      nrow(primary_stage2),
    primary_median_gap_reduction =
      safe_median(
        primary_stage2$
          absolute_gap_reduction,
        0
      ),
    primary_median_recovery_fraction =
      safe_median(
        primary_stage2$
          recovery_fraction,
        0
      ),
    primary_positive_fraction =
      safe_mean(
        primary_stage2$
          absolute_gap_reduction > 0,
        0
      ),
    allsize_positive_fraction =
      safe_mean(
        all_stage2$
          absolute_gap_reduction > 0,
        0
      ),
    sample_improvement_fraction =
      safe_median(
        primary_stage2$
          sample_improvement_fraction,
        0
      ),
    supported_primary_median_gap_reduction =
      safe_median(
        supported_primary$
          absolute_gap_reduction,
        0
      ),
    inflammation_median_gap_reduction =
      safe_median(
        inflammation$
          absolute_gap_reduction,
        0
      ),
    myeloid_identity_absolute_change =
      safe_median(
        abs(
          myeloid$
            perturbed_hfpef_mean -
            myeloid$
              observed_hfpef_mean
        ),
        0
      )
  )
}

evaluate_tf <- function(
  tf_requested,
  network_table,
  activity_matrix,
  sample_z,
  sample_meta,
  program_definitions,
  program_summary,
  observed_program_scores
) {
  simulation <- simulate_normalization(
    observed_z = sample_z,
    activity_matrix =
      activity_matrix,
    tf_requested =
      tf_requested,
    network_table =
      network_table,
    sample_meta =
      sample_meta
  )

  perturbed_scores <- score_programs(
    simulation$perturbed_z,
    program_definitions
  )

  detail <- summarize_perturbation(
    observed_scores =
      observed_program_scores,
    perturbed_scores =
      perturbed_scores,
    sample_meta =
      sample_meta,
    program_summary =
      program_summary
  )

  summary <- result_summary(
    detail
  )

  summary[
    ,
    `:=`(
      tf_symbol =
        tf_requested,
      target_count =
        simulation$target_count,
      global_rms_shift_HFpEF =
        simulation$
          global_rms_shift_HFpEF,
      specificity_score =
        1 /
        (
          1 +
            simulation$
              global_rms_shift_HFpEF +
            myeloid_identity_absolute_change
        )
    )
  ]

  list(
    detail = detail,
    summary = summary
  )
}

############################################################
## 5. Load and lock input objects
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

network_table <- read_table_auto(
  STAGE4_NETWORK_FILE
)

locked_regulon_sizes <- read_table_auto(
  STAGE4_REGULON_SIZE_FILE
)

program_manifest <- read_table_auto(
  STAGE4_PROGRAM_FILE
)

stage5_candidate_table <-
  data.table::fread(
    STAGE5_CANDIDATE_FILE,
    encoding = "UTF-8"
  )

pseudobulk_objects <- readRDS(
  STAGE4_PSEUDOBULK_RDS
)

if (
  !"sample_logcpm" %in%
    names(pseudobulk_objects)
) {
  stop(
    "Stage 4 pseudobulk object does not contain sample_logcpm."
  )
}

sample_logcpm <- as.matrix(
  pseudobulk_objects$
    sample_logcpm
)

stage4_activity <- as.matrix(
  readRDS(
    STAGE4_ACTIVITY_RDS
  )
)

sample_order <-
  sample_meta$sample_accession

if (
  !all(
    sample_order %in%
      colnames(sample_logcpm)
  ) ||
  !all(
    sample_order %in%
      colnames(stage4_activity)
  )
) {
  stop(
    "One or more locked samples are missing from the Stage 4 matrices."
  )
}

sample_logcpm <- sample_logcpm[
  ,
  sample_order,
  drop = FALSE
]

stage4_activity <- stage4_activity[
  ,
  sample_order,
  drop = FALSE
]

if ("tf_symbol" %in% names(network_table)) {
  network_table[
    ,
    tf_symbol := NULL
  ]
}

required_network_columns <- c(
  "source_symbol",
  "target_feature",
  "mor",
  "weight"
)

missing_network_columns <- setdiff(
  required_network_columns,
  names(network_table)
)

if (
  length(
    missing_network_columns
  ) > 0L
) {
  stop(
    "Stage 4 network is missing column(s): ",
    paste(
      missing_network_columns,
      collapse = ", "
    )
  )
}

network_table[
  ,
  mor := as.numeric(mor)
]

network_table[
  ,
  weight := as.numeric(weight)
]

network_table <- network_table[
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
  network_table,
  source_symbol,
  target_feature,
  -weight
)

network_table <- network_table[
  ,
  .SD[1L],
  by = .(
    source_symbol,
    target_feature
  )
]

observed_regulon_sizes <-
  network_table[
    ,
    .(
      observed_regulon_size =
        data.table::uniqueN(
          target_feature
        )
    ),
    by = source_symbol
  ]

network_integrity <- merge(
  locked_regulon_sizes[
    ,
    .(
      source_symbol,
      locked_regulon_size =
        as.integer(
          regulon_size
        )
    )
  ],
  observed_regulon_sizes,
  by = "source_symbol",
  all = TRUE
)

network_integrity[
  ,
  target_count_match := (
    !is.na(
      locked_regulon_size
    ) &
      !is.na(
        observed_regulon_size
      ) &
      locked_regulon_size ==
        observed_regulon_size
  )
]

if (
  any(
    network_integrity$
      target_count_match != TRUE
  )
) {
  failed_sources <- network_integrity[
    target_count_match != TRUE,
    source_symbol
  ]

  stop(
    "Stage 4 network integrity mismatch for: ",
    paste(
      failed_sources,
      collapse = ", "
    )
  )
}

candidate_key_table <-
  data.table::data.table(
    requested_tf =
      CANDIDATE_TFS_REQUESTED,
    requested_order =
      seq_along(
        CANDIDATE_TFS_REQUESTED
      ),
    tf_key =
      gene_key(
        CANDIDATE_TFS_REQUESTED
      )
  )

network_source_map <- unique(
  network_table[
    ,
    .(
      tf_key =
        gene_key(source_symbol),
      source_symbol
    )
  ],
  by = "tf_key"
)

candidate_resolution <- merge(
  candidate_key_table,
  network_source_map,
  by = "tf_key",
  all.x = TRUE,
  sort = FALSE
)

data.table::setorder(
  candidate_resolution,
  requested_order
)

candidate_resolution[
  ,
  activity_available := (
    !is.na(source_symbol) &
      source_symbol %in%
        rownames(stage4_activity)
  )
]

if (
  nrow(candidate_resolution) !=
    length(
      CANDIDATE_TFS_REQUESTED
    ) ||
  any(
    candidate_resolution$
      activity_available != TRUE
  )
) {
  stop(
    "At least one prespecified candidate TF is unavailable in the locked Stage 4 network or activity matrix."
  )
}

candidate_tfs <-
  candidate_resolution$
    source_symbol

feature_map <- make_feature_map(
  rownames(sample_logcpm)
)

program_objects <- make_program_definitions(
  program_manifest,
  feature_map
)

program_definitions <-
  program_objects$definitions

program_summary <-
  program_objects$summary

sample_z <- scale_rows(
  sample_logcpm
)

observed_program_scores <-
  score_programs(
    sample_z,
    program_definitions
  )

eligible_primary_programs <- program_summary[
  program_category ==
    "Stage2_drug_opposed_net" &
    primary_program == TRUE,
  data.table::uniqueN(
    program_id
  )
]

if (
  eligible_primary_programs <
    2L
) {
  stop(
    "Fewer than two primary Stage 2 programs are available."
  )
}

write_csv_safe(
  network_integrity,
  file.path(
    DIRS$tables,
    "00_stage5B_network_integrity_audit.csv"
  )
)

write_csv_safe(
  candidate_resolution,
  file.path(
    DIRS$tables,
    "01_stage5B_candidate_resolution.csv"
  )
)

write_csv_safe(
  program_summary,
  file.path(
    DIRS$tables,
    "02_stage5B_program_definition_summary.csv"
  )
)

############################################################
## 6. Recalculate full candidate effects and audit Stage 5
############################################################

full_detail_records <- list()
full_summary_records <- list()

for (
  candidate_i in candidate_tfs
) {
  log_msg(
    "Full-regulon candidate evaluation: ",
    candidate_i
  )

  candidate_result <- evaluate_tf(
    tf_requested =
      candidate_i,
    network_table =
      network_table,
    activity_matrix =
      stage4_activity,
    sample_z =
      sample_z,
    sample_meta =
      sample_meta,
    program_definitions =
      program_definitions,
    program_summary =
      program_summary,
    observed_program_scores =
      observed_program_scores
  )

  detail_i <- data.table::copy(
    candidate_result$detail
  )

  detail_i[
    ,
    tf_symbol := candidate_i
  ]

  full_detail_records[[length(full_detail_records) + 1L]] <- detail_i

  full_summary_records[[length(full_summary_records) + 1L]] <- candidate_result$summary
}

full_candidate_detail <-
  data.table::rbindlist(
    full_detail_records,
    use.names = TRUE,
    fill = TRUE
  )

full_candidate_summary <-
  data.table::rbindlist(
    full_summary_records,
    use.names = TRUE,
    fill = TRUE
  )

stage5_required_columns <- c(
  "tf_symbol",
  "stage2_primary_median_gap_reduction"
)

missing_stage5_columns <- setdiff(
  stage5_required_columns,
  names(stage5_candidate_table)
)

if (
  length(
    missing_stage5_columns
  ) > 0L
) {
  stop(
    "Stage 5 candidate table is missing audit column(s): ",
    paste(
      missing_stage5_columns,
      collapse = ", "
    )
  )
}

stage5_audit <- merge(
  full_candidate_summary[
    ,
    .(
      tf_symbol,
      recalculated_primary_gap_reduction =
        primary_median_gap_reduction,
      recalculated_target_count =
        target_count
    )
  ],
  stage5_candidate_table[
    tf_symbol %in%
      candidate_tfs,
    .(
      tf_symbol,
      stage5_primary_gap_reduction =
        stage2_primary_median_gap_reduction
    )
  ],
  by = "tf_symbol",
  all = TRUE
)

stage5_audit[
  ,
  absolute_difference := abs(
    recalculated_primary_gap_reduction -
      stage5_primary_gap_reduction
  )
]

stage5_audit[
  ,
  within_tolerance := (
    is.finite(
      absolute_difference
    ) &
      absolute_difference <=
        STAGE5_AUDIT_TOLERANCE
  )
]

if (
  nrow(stage5_audit) !=
    length(candidate_tfs) ||
  any(
    stage5_audit$
      within_tolerance != TRUE
  )
) {
  stop(
    "Recalculated full-regulon candidate results do not reproduce Stage 5 v2 within tolerance."
  )
}

write_csv_safe(
  full_candidate_detail,
  file.path(
    DIRS$tables,
    "03_stage5B_full_candidate_program_effects.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  full_candidate_summary,
  file.path(
    DIRS$tables,
    "04_stage5B_full_candidate_summary.csv"
  )
)

write_csv_safe(
  stage5_audit,
  file.path(
    DIRS$tables,
    "05_stage5B_stage5_v2_reproduction_audit.csv"
  )
)

############################################################
## 7. Candidate regulon bootstrap
############################################################

bootstrap_records <- list()

for (
  candidate_i in candidate_tfs
) {
  candidate_network <- data.table::copy(
    network_table[
      source_symbol ==
        candidate_i
    ]
  )

  data.table::setorder(
    candidate_network,
    target_feature,
    -weight
  )

  candidate_network <-
    candidate_network[
      ,
      .SD[1L],
      by = target_feature
    ]

  full_target_count <-
    nrow(candidate_network)

  bootstrap_target_count <- max(
    MIN_TARGETS_PER_REGULON,
    floor(
      BOOTSTRAP_TARGET_FRACTION *
        full_target_count
    )
  )

  bootstrap_target_count <- min(
    bootstrap_target_count,
    full_target_count
  )

  log_msg(
    "Bootstrap: ",
    candidate_i,
    " | full targets=",
    full_target_count,
    " | retained=",
    bootstrap_target_count,
    " | replicates=",
    BOOTSTRAP_REPLICATES
  )

  candidate_seed_offset <- match(
    candidate_i,
    candidate_tfs
  )

  for (
    bootstrap_i in seq_len(
      BOOTSTRAP_REPLICATES
    )
  ) {
    set.seed(
      20260714 +
        candidate_seed_offset *
          100000L +
        bootstrap_i
    )

    sampled_rows <- sample.int(
      full_target_count,
      size =
        bootstrap_target_count,
      replace = FALSE
    )

    bootstrap_network <-
      candidate_network[
        sampled_rows
      ]

    bootstrap_activity <-
      weighted_regulon_activity(
        expression_matrix =
          sample_logcpm,
        network_subset =
          bootstrap_network
      )

    if (
      is.null(
        bootstrap_activity
      ) ||
      !candidate_i %in%
        rownames(
          bootstrap_activity
        )
    ) {
      stop(
        "Bootstrap TF activity could not be computed for ",
        candidate_i,
        " replicate ",
        bootstrap_i,
        "."
      )
    }

    bootstrap_result <- evaluate_tf(
      tf_requested =
        candidate_i,
      network_table =
        bootstrap_network,
      activity_matrix =
        bootstrap_activity,
      sample_z =
        sample_z,
      sample_meta =
        sample_meta,
      program_definitions =
        program_definitions,
      program_summary =
        program_summary,
      observed_program_scores =
        observed_program_scores
    )

    summary_i <- bootstrap_result$summary

    bootstrap_records[[length(bootstrap_records) + 1L]] <- data.table::data.table(
      tf_symbol =
        candidate_i,
      bootstrap_replicate =
        bootstrap_i,
      full_target_count =
        full_target_count,
      bootstrap_target_count =
        bootstrap_target_count,
      primary_median_gap_reduction =
        summary_i$
          primary_median_gap_reduction,
      primary_median_recovery_fraction =
        summary_i$
          primary_median_recovery_fraction,
      primary_positive_fraction =
        summary_i$
          primary_positive_fraction,
      allsize_positive_fraction =
        summary_i$
          allsize_positive_fraction,
      sample_improvement_fraction =
        summary_i$
          sample_improvement_fraction,
      supported_primary_median_gap_reduction =
        summary_i$
          supported_primary_median_gap_reduction,
      inflammation_median_gap_reduction =
        summary_i$
          inflammation_median_gap_reduction,
      global_rms_shift_HFpEF =
        summary_i$
          global_rms_shift_HFpEF,
      specificity_score =
        summary_i$
          specificity_score
    )
  }
}

bootstrap_results <-
  data.table::rbindlist(
    bootstrap_records,
    use.names = TRUE,
    fill = TRUE
  )

bootstrap_results[
  ,
  bootstrap_candidate_rank :=
    rank(
      -primary_median_gap_reduction,
      ties.method = "average"
    ),
  by = bootstrap_replicate
]

bootstrap_summary <- bootstrap_results[
  ,
  .(
    bootstrap_replicates = .N,
    median_primary_gap_reduction =
      safe_median(
        primary_median_gap_reduction
      ),
    q025_primary_gap_reduction =
      safe_quantile(
        primary_median_gap_reduction,
        0.025
      ),
    q975_primary_gap_reduction =
      safe_quantile(
        primary_median_gap_reduction,
        0.975
      ),
    positive_recovery_probability =
      safe_mean(
        primary_median_gap_reduction >
          0
      ),
    median_primary_positive_fraction =
      safe_median(
        primary_positive_fraction
      ),
    median_allsize_positive_fraction =
      safe_median(
        allsize_positive_fraction
      ),
    median_sample_improvement_fraction =
      safe_median(
        sample_improvement_fraction
      ),
    median_supported_recovery =
      safe_median(
        supported_primary_median_gap_reduction
      ),
    median_inflammation_recovery =
      safe_median(
        inflammation_median_gap_reduction
      ),
    median_global_rms_shift =
      safe_median(
        global_rms_shift_HFpEF
      ),
    median_specificity_score =
      safe_median(
        specificity_score
      ),
    median_bootstrap_rank =
      safe_median(
        bootstrap_candidate_rank
      ),
    top1_frequency =
      safe_mean(
        bootstrap_candidate_rank <=
          1
      ),
    top3_frequency =
      safe_mean(
        bootstrap_candidate_rank <=
          3
      )
  ),
  by = tf_symbol
]

if (
  nrow(bootstrap_results) !=
    BOOTSTRAP_REPLICATES *
      length(candidate_tfs)
) {
  stop(
    "Bootstrap output row count is incomplete."
  )
}

if (
  min(
    bootstrap_results[
      ,
      .N,
      by = tf_symbol
    ]$N
  ) !=
    BOOTSTRAP_REPLICATES
) {
  stop(
    "At least one candidate does not have all bootstrap replicates."
  )
}

write_csv_safe(
  bootstrap_results,
  file.path(
    DIRS$tables,
    "06_stage5B_candidate_regulon_bootstrap_replicates.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  bootstrap_summary,
  file.path(
    DIRS$tables,
    "07_stage5B_candidate_regulon_bootstrap_summary.csv"
  )
)

############################################################
## 8. Build TF covariates for unbiased null matching
############################################################

regulon_covariates <-
  network_table[
    ,
    .(
      regulon_size =
        data.table::uniqueN(
          target_feature
        )
    ),
    by = source_symbol
  ]

activity_effect_records <- lapply(
  rownames(stage4_activity),
  function(source_i) {
    values <- as.numeric(
      stage4_activity[
        source_i,
        sample_order
      ]
    )

    data.table::data.table(
      source_symbol =
        source_i,
      activity_effect =
        mean(
          values[
            sample_meta$condition ==
              "HFpEF"
          ],
          na.rm = TRUE
        ) -
        mean(
          values[
            sample_meta$condition ==
              "Control"
          ],
          na.rm = TRUE
        )
    )
  }
)

activity_covariates <-
  data.table::rbindlist(
    activity_effect_records
  )

source_feature_map <-
  map_symbols_to_features(
    regulon_covariates$
      source_symbol,
    feature_map
  )

source_feature_map <- unique(
  source_feature_map[
    ,
    .(
      source_symbol =
        symbol,
      source_feature =
        feature
    )
  ],
  by = "source_symbol"
)

tf_covariates <- merge(
  regulon_covariates,
  activity_covariates,
  by = "source_symbol",
  all = FALSE
)

tf_covariates <- merge(
  tf_covariates,
  source_feature_map,
  by = "source_symbol",
  all.x = TRUE
)

tf_covariates[
  ,
  expression_effect := NA_real_
]

detected_sources <- tf_covariates[
  !is.na(source_feature) &
    source_feature %in%
      rownames(sample_logcpm)
]

if (
  nrow(detected_sources) > 0L
) {
  expression_matrix <- sample_logcpm[
    detected_sources$
      source_feature,
    sample_order,
    drop = FALSE
  ]

  expression_effects <- apply(
    expression_matrix,
    1L,
    function(values) {
      mean(
        values[
          sample_meta$condition ==
            "HFpEF"
        ],
        na.rm = TRUE
      ) -
        mean(
          values[
            sample_meta$condition ==
              "Control"
        ],
        na.rm = TRUE
      )
    }
  )

  tf_covariates[
    match(
      detected_sources$
        source_symbol,
      source_symbol
    ),
    expression_effect :=
      as.numeric(
        expression_effects
      )
  ]
}

tf_covariates[
  ,
  source_key :=
    gene_key(source_symbol)
]

tf_covariates[
  ,
  log_regulon_size :=
    log1p(regulon_size)
]

write_csv_safe(
  tf_covariates,
  file.path(
    DIRS$tables,
    "08_stage5B_all_TF_matching_covariates.csv"
  )
)

############################################################
## 9. Create matched null pools without program-overlap filter
############################################################

candidate_keys <- gene_key(
  candidate_tfs
)

eligible_null_tfs <- tf_covariates[
  !source_key %in%
    candidate_keys &
    regulon_size >=
      MIN_TARGETS_PER_REGULON &
    is.finite(activity_effect)
]

if (
  nrow(eligible_null_tfs) <
    NULL_MIN_POOL_SIZE
) {
  stop(
    "Insufficient eligible TFs for null matching."
  )
}

scale_log_size <- stats::mad(
  tf_covariates$
    log_regulon_size,
  na.rm = TRUE
)

scale_activity <- stats::mad(
  abs(
    tf_covariates$
      activity_effect
  ),
  na.rm = TRUE
)

scale_expression <- stats::mad(
  abs(
    tf_covariates$
      expression_effect
  ),
  na.rm = TRUE
)

if (
  !is.finite(scale_log_size) ||
  scale_log_size <= 0
) {
  scale_log_size <- 1
}

if (
  !is.finite(scale_activity) ||
  scale_activity <= 0
) {
  scale_activity <- 1
}

if (
  !is.finite(scale_expression) ||
  scale_expression <= 0
) {
  scale_expression <- 1
}

null_pool_records <- list()
required_null_tfs <- character()

for (
  candidate_i in candidate_tfs
) {
  candidate_covariate <- tf_covariates[
    source_symbol ==
      candidate_i
  ]

  if (
    nrow(candidate_covariate) != 1L
  ) {
    stop(
      "Candidate covariates are not uniquely available for ",
      candidate_i,
      "."
    )
  }

  pool_i <- data.table::copy(
    eligible_null_tfs
  )

  pool_i[
    ,
    size_distance :=
      abs(
        log_regulon_size -
          candidate_covariate$
            log_regulon_size[1L]
      ) /
      scale_log_size
  ]

  pool_i[
    ,
    activity_distance :=
      abs(
        abs(activity_effect) -
          abs(
            candidate_covariate$
              activity_effect[1L]
          )
      ) /
      scale_activity
  ]

  candidate_expression_available <-
    is.finite(
      candidate_covariate$
        expression_effect[1L]
    )

  if (
    candidate_expression_available
  ) {
    pool_i[
      ,
      expression_distance := ifelse(
        is.finite(expression_effect),
        abs(
          abs(expression_effect) -
            abs(
              candidate_covariate$
                expression_effect[1L]
            )
        ) /
          scale_expression,
        2
      )
    ]
  } else {
    pool_i[
      ,
      expression_distance := 0
    ]
  }

  pool_i[
    ,
    matching_distance :=
      size_distance +
      activity_distance +
      expression_distance
  ]

  data.table::setorder(
    pool_i,
    matching_distance,
    source_symbol
  )

  pool_i <- head(
    pool_i,
    min(
      NULL_NEAREST_POOL_SIZE,
      nrow(pool_i)
    )
  )

  if (
    nrow(pool_i) <
      NULL_MIN_POOL_SIZE
  ) {
    stop(
      "Matched null pool is too small for ",
      candidate_i,
      "."
    )
  }

  positive_distances <-
    pool_i$matching_distance[
      is.finite(
        pool_i$matching_distance
      ) &
        pool_i$matching_distance > 0
    ]

  distance_scale <- safe_median(
    positive_distances,
    1
  )

  if (
    !is.finite(distance_scale) ||
    distance_scale <= 0
  ) {
    distance_scale <- 1
  }

  pool_i[
    ,
    sampling_weight :=
      exp(
        -matching_distance /
          distance_scale
      )
  ]

  pool_i[
    ,
    sampling_probability :=
      sampling_weight /
      sum(sampling_weight)
  ]

  pool_i[
    ,
    `:=`(
      candidate_tf =
        candidate_i,
      candidate_regulon_size =
        candidate_covariate$
          regulon_size[1L],
      candidate_activity_effect =
        candidate_covariate$
          activity_effect[1L],
      candidate_expression_effect =
        candidate_covariate$
          expression_effect[1L],
      expression_dimension_used =
        candidate_expression_available,
      program_overlap_filter_used =
        FALSE
    )
  ]

  null_pool_records[[length(null_pool_records) + 1L]] <- pool_i

  required_null_tfs <- unique(
    c(
      required_null_tfs,
      pool_i$source_symbol
    )
  )
}

null_pool <- data.table::rbindlist(
  null_pool_records,
  use.names = TRUE,
  fill = TRUE
)

write_csv_safe(
  null_pool,
  file.path(
    DIRS$tables,
    "09_stage5B_candidate_matched_null_TF_pools.csv"
  ),
  compress = TRUE
)

############################################################
## 10. Precompute full-regulon effects for null TFs
############################################################

log_msg(
  "Precomputing effects for ",
  length(required_null_tfs),
  " unique matched null TFs."
)

null_effect_records <- list()

for (
  null_tf_i in required_null_tfs
) {
  if (
    !null_tf_i %in%
      rownames(stage4_activity)
  ) {
    stop(
      "Matched null TF lacks Stage 4 activity: ",
      null_tf_i
    )
  }

  null_result <- evaluate_tf(
    tf_requested =
      null_tf_i,
    network_table =
      network_table,
    activity_matrix =
      stage4_activity,
    sample_z =
      sample_z,
    sample_meta =
      sample_meta,
    program_definitions =
      program_definitions,
    program_summary =
      program_summary,
    observed_program_scores =
      observed_program_scores
  )

  summary_i <- null_result$summary

  null_effect_records[[length(null_effect_records) + 1L]] <- data.table::data.table(
    null_tf =
      null_tf_i,
    primary_median_gap_reduction =
      summary_i$
        primary_median_gap_reduction,
    primary_positive_fraction =
      summary_i$
        primary_positive_fraction,
    allsize_positive_fraction =
      summary_i$
        allsize_positive_fraction,
    sample_improvement_fraction =
      summary_i$
        sample_improvement_fraction,
    supported_primary_median_gap_reduction =
      summary_i$
        supported_primary_median_gap_reduction,
    inflammation_median_gap_reduction =
      summary_i$
        inflammation_median_gap_reduction,
    global_rms_shift_HFpEF =
      summary_i$
        global_rms_shift_HFpEF,
    specificity_score =
      summary_i$
        specificity_score
  )
}

null_effects <- data.table::rbindlist(
  null_effect_records,
  use.names = TRUE,
  fill = TRUE
)

if (
  data.table::uniqueN(
    null_effects$null_tf
  ) !=
    length(required_null_tfs)
) {
  stop(
    "Not all matched null TF effects were computed."
  )
}

write_csv_safe(
  null_effects,
  file.path(
    DIRS$tables,
    "10_stage5B_precomputed_null_TF_effects.csv"
  )
)

############################################################
## 11. Draw empirical null distributions
############################################################

null_draw_records <- list()
null_summary_records <- list()

for (
  candidate_i in candidate_tfs
) {
  candidate_pool <- null_pool[
    candidate_tf ==
      candidate_i
  ]

  candidate_pool <- merge(
    candidate_pool,
    null_effects,
    by.x = "source_symbol",
    by.y = "null_tf",
    all.x = TRUE
  )

  candidate_pool <- candidate_pool[
    is.finite(
      primary_median_gap_reduction
    ) &
      is.finite(
        sampling_probability
      ) &
      sampling_probability > 0
  ]

  if (
    nrow(candidate_pool) <
      NULL_MIN_POOL_SIZE
  ) {
    stop(
      "Evaluated null pool is too small for ",
      candidate_i,
      "."
    )
  }

  candidate_pool[
    ,
    sampling_probability :=
      sampling_probability /
      sum(sampling_probability)
  ]

  set.seed(
    20260714 +
      match(
        candidate_i,
        candidate_tfs
      ) *
      10000L
  )

  sampled_rows <- sample.int(
    nrow(candidate_pool),
    size =
      NULL_DRAWS_PER_CANDIDATE,
    replace = TRUE,
    prob =
      candidate_pool$
        sampling_probability
  )

  draws_i <- candidate_pool[
    sampled_rows
  ]

  draws_i[
    ,
    `:=`(
      candidate_tf =
        candidate_i,
      draw_id =
        seq_len(.N)
    )
  ]

  candidate_effect <- full_candidate_summary[
    tf_symbol ==
      candidate_i,
    primary_median_gap_reduction
  ]

  if (
    length(candidate_effect) != 1L
  ) {
    stop(
      "Candidate full-regulon effect is not uniquely available for ",
      candidate_i,
      "."
    )
  }

  candidate_percentile <- mean(
    draws_i$
      primary_median_gap_reduction <=
      candidate_effect,
    na.rm = TRUE
  )

  empirical_p <- (
    1 +
      sum(
        draws_i$
          primary_median_gap_reduction >=
          candidate_effect,
        na.rm = TRUE
      )
  ) /
    (
      1 +
        nrow(draws_i)
    )

  null_summary_records[[length(null_summary_records) + 1L]] <- data.table::data.table(
    tf_symbol =
      candidate_i,
    candidate_effect =
      candidate_effect,
    null_pool_size =
      nrow(candidate_pool),
    null_draws =
      nrow(draws_i),
    unique_null_TFs_drawn =
      data.table::uniqueN(
        draws_i$source_symbol
      ),
    null_median =
      safe_median(
        draws_i$
          primary_median_gap_reduction
      ),
    null_q025 =
      safe_quantile(
        draws_i$
          primary_median_gap_reduction,
        0.025
      ),
    null_q975 =
      safe_quantile(
        draws_i$
          primary_median_gap_reduction,
        0.975
      ),
    candidate_percentile =
      candidate_percentile,
    empirical_one_sided_p =
      empirical_p,
    program_overlap_filter_used =
      FALSE
  )

  null_draw_records[[length(null_draw_records) + 1L]] <- draws_i[
    ,
    .(
      candidate_tf,
      draw_id,
      sampled_null_TF =
        source_symbol,
      matching_distance,
      sampling_probability,
      primary_median_gap_reduction,
      primary_positive_fraction,
      allsize_positive_fraction,
      sample_improvement_fraction,
      supported_primary_median_gap_reduction,
      inflammation_median_gap_reduction,
      global_rms_shift_HFpEF,
      specificity_score
    )
  ]
}

null_draws <- data.table::rbindlist(
  null_draw_records,
  use.names = TRUE,
  fill = TRUE
)

null_summary <- data.table::rbindlist(
  null_summary_records,
  use.names = TRUE,
  fill = TRUE
)

if (
  nrow(null_draws) !=
    NULL_DRAWS_PER_CANDIDATE *
      length(candidate_tfs)
) {
  stop(
    "Random matched-TF null draw count is incomplete."
  )
}

if (
  nrow(null_summary) !=
    length(candidate_tfs)
) {
  stop(
    "Random null summary does not contain all candidates."
  )
}

write_csv_safe(
  null_draws,
  file.path(
    DIRS$tables,
    "11_stage5B_random_matched_TF_null_draws.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  null_summary,
  file.path(
    DIRS$tables,
    "12_stage5B_random_matched_TF_null_summary.csv"
  )
)

############################################################
## 12. Final transparent candidate robustness rank
############################################################

final_robustness <- merge(
  full_candidate_summary,
  bootstrap_summary,
  by = "tf_symbol",
  all = TRUE,
  suffixes = c(
    "_full",
    "_bootstrap"
  )
)

final_robustness <- merge(
  final_robustness,
  null_summary[
    ,
    .(
      tf_symbol,
      null_median,
      null_q025,
      null_q975,
      candidate_percentile,
      empirical_one_sided_p
    )
  ],
  by = "tf_symbol",
  all = TRUE
)

final_robustness[
  ,
  rank_full_recovery :=
    rank_metric(
      primary_median_gap_reduction,
      higher_is_better = TRUE
    )
]

final_robustness[
  ,
  rank_bootstrap_recovery :=
    rank_metric(
      median_primary_gap_reduction,
      higher_is_better = TRUE
    )
]

final_robustness[
  ,
  rank_bootstrap_positive :=
    rank_metric(
      positive_recovery_probability,
      higher_is_better = TRUE
    )
]

final_robustness[
  ,
  rank_bootstrap_top3 :=
    rank_metric(
      top3_frequency,
      higher_is_better = TRUE
    )
]

final_robustness[
  ,
  rank_null_percentile :=
    rank_metric(
      candidate_percentile,
      higher_is_better = TRUE
    )
]

final_robustness[
  ,
  rank_empirical_p :=
    rank_metric(
      empirical_one_sided_p,
      higher_is_better = FALSE
    )
]

final_robustness[
  ,
  rank_sample_stability :=
    rank_metric(
      sample_improvement_fraction,
      higher_is_better = TRUE
    )
]

final_robustness[
  ,
  rank_inflammation :=
    rank_metric(
      inflammation_median_gap_reduction,
      higher_is_better = TRUE
    )
]

final_robustness[
  ,
  rank_specificity :=
    rank_metric(
      specificity_score,
      higher_is_better = TRUE
    )
]

rank_columns <- c(
  "rank_full_recovery",
  "rank_bootstrap_recovery",
  "rank_bootstrap_positive",
  "rank_bootstrap_top3",
  "rank_null_percentile",
  "rank_empirical_p",
  "rank_sample_stability",
  "rank_inflammation",
  "rank_specificity"
)

final_robustness[
  ,
  mean_evidence_rank :=
    rowMeans(
      .SD,
      na.rm = TRUE
    ),
  .SDcols =
    rank_columns
]

final_robustness[
  ,
  median_evidence_rank :=
    apply(
      .SD,
      1L,
      stats::median,
      na.rm = TRUE
    ),
  .SDcols =
    rank_columns
]

final_robustness[
  ,
  final_robustness_score :=
    0.5 *
      mean_evidence_rank +
    0.5 *
      median_evidence_rank
]

data.table::setorder(
  final_robustness,
  final_robustness_score,
  rank_null_percentile,
  rank_bootstrap_recovery
)

final_robustness[
  ,
  final_robustness_rank :=
    seq_len(.N)
]

final_robustness[
  ,
  Nfkb1_forced := FALSE
]

write_csv_safe(
  final_robustness,
  file.path(
    DIRS$tables,
    "13_stage5B_final_candidate_robustness_rank.csv"
  )
)

############################################################
## 13. Save checkpoint before figures
############################################################

saveRDS(
  list(
    network_integrity =
      network_integrity,
    candidate_resolution =
      candidate_resolution,
    program_summary =
      program_summary,
    full_candidate_detail =
      full_candidate_detail,
    full_candidate_summary =
      full_candidate_summary,
    stage5_audit =
      stage5_audit,
    bootstrap_results =
      bootstrap_results,
    bootstrap_summary =
      bootstrap_summary,
    tf_covariates =
      tf_covariates,
    null_pool =
      null_pool,
    null_effects =
      null_effects,
    null_draws =
      null_draws,
    null_summary =
      null_summary,
    final_robustness =
      final_robustness
  ),
  file.path(
    DIRS$objects,
    "CHECKPOINT_stage5B_OFFLINE_scientific_results_pre_figures.rds"
  ),
  compress = FALSE
)

log_msg(
  "Stage 5B OFFLINE scientific calculations checkpointed before figures."
)

############################################################
## 14. Figures and source data
############################################################

candidate_order <- final_robustness[
  order(
    final_robustness_rank
  ),
  tf_symbol
]

full_plot_data <- data.table::copy(
  full_candidate_summary
)

full_plot_data[
  ,
  tf_symbol := factor(
    tf_symbol,
    levels =
      candidate_order
  )
]

write_csv_safe(
  full_plot_data,
  file.path(
    DIRS$source,
    "Fig5B_A_full_candidate_recovery_source.csv"
  )
)

plot_full <- ggplot2::ggplot(
  full_plot_data,
  ggplot2::aes(
    x = tf_symbol,
    y =
      primary_median_gap_reduction,
    fill = tf_symbol
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  ggplot2::geom_col(
    width = 0.72,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title =
      "Full-regulon candidate perturbation effects",
    subtitle =
      "Positive values indicate reduced HFpEF-Control separation of primary drug-opposed programs",
    x = NULL,
    y =
      "Median primary-program gap reduction"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
  )

save_plot_bundle(
  plot_full,
  "Fig5B_A_full_candidate_recovery",
  8.5,
  6
)

bootstrap_plot_data <- data.table::copy(
  bootstrap_results
)

bootstrap_plot_data[
  ,
  tf_symbol := factor(
    tf_symbol,
    levels =
      candidate_order
  )
]

write_csv_safe(
  bootstrap_plot_data,
  file.path(
    DIRS$source,
    "Fig5B_B_regulon_bootstrap_source.csv"
  ),
  compress = TRUE
)

plot_bootstrap <- ggplot2::ggplot(
  bootstrap_plot_data,
  ggplot2::aes(
    x = tf_symbol,
    y =
      primary_median_gap_reduction,
    fill = tf_symbol
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  ggplot2::geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title =
      "Regulon-target bootstrap robustness",
    subtitle =
      paste0(
        BOOTSTRAP_REPLICATES,
        " replicates per candidate; ",
        BOOTSTRAP_TARGET_FRACTION * 100,
        "% of regulon targets retained"
      ),
    x = NULL,
    y =
      "Primary-program gap reduction"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
  )

save_plot_bundle(
  plot_bootstrap,
  "Fig5B_B_regulon_bootstrap_distributions",
  9,
  6.5
)

null_plot_data <- data.table::copy(
  null_summary
)

null_plot_data[
  ,
  tf_symbol := factor(
    tf_symbol,
    levels =
      candidate_order
  )
]

write_csv_safe(
  null_plot_data,
  file.path(
    DIRS$source,
    "Fig5B_C_random_null_percentile_source.csv"
  )
)

plot_null <- ggplot2::ggplot(
  null_plot_data,
  ggplot2::aes(
    x = tf_symbol,
    y =
      candidate_percentile,
    fill = tf_symbol
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.95,
    linetype = 2
  ) +
  ggplot2::geom_col(
    width = 0.72,
    show.legend = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  ggplot2::labs(
    title =
      "Candidate percentiles in unbiased matched-TF null distributions",
    subtitle =
      "Matching used regulon size, TF-activity effect, and available TF-expression effect; no program-overlap filter",
    x = NULL,
    y =
      "Candidate percentile"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
  )

save_plot_bundle(
  plot_null,
  "Fig5B_C_random_matched_TF_null_percentiles",
  9,
  6.5
)

rank_plot_data <- data.table::copy(
  final_robustness
)

rank_plot_data[
  ,
  tf_symbol := factor(
    tf_symbol,
    levels = rev(
      candidate_order
    )
  )
]

write_csv_safe(
  rank_plot_data,
  file.path(
    DIRS$source,
    "Fig5B_D_final_robustness_rank_source.csv"
  )
)

plot_rank <- ggplot2::ggplot(
  rank_plot_data,
  ggplot2::aes(
    x =
      -final_robustness_score,
    y = tf_symbol
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x =
        -max(
          final_robustness_score
        ),
      xend =
        -final_robustness_score,
      y = tf_symbol,
      yend = tf_symbol
    ),
    color = "grey65",
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size =
        positive_recovery_probability,
      fill =
        candidate_percentile
    ),
    shape = 21,
    color = "black"
  ) +
  ggplot2::scale_size_continuous(
    limits = c(
      0,
      1
    ),
    range = c(
      3,
      8
    )
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "black"
  ) +
  ggplot2::labs(
    title =
      "Offline Stage 5B candidate robustness rank",
    subtitle =
      "Rank aggregation across full-regulon recovery, bootstrap stability, matched-null performance, inflammation, and specificity",
    x =
      "Higher robustness performance",
    y = NULL,
    size =
      "Bootstrap positive\nprobability",
    fill =
      "Matched-null\npercentile"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  )

save_plot_bundle(
  plot_rank,
  "Fig5B_D_final_candidate_robustness_rank",
  9,
  6.5
)

heatmap_columns <- c(
  "primary_median_gap_reduction",
  "median_primary_gap_reduction",
  "positive_recovery_probability",
  "top3_frequency",
  "candidate_percentile",
  "inflammation_median_gap_reduction",
  "specificity_score"
)

heatmap_table <- final_robustness[
  match(
    candidate_order,
    tf_symbol
  )
]

heatmap_matrix <- as.matrix(
  heatmap_table[
    ,
    ..heatmap_columns
  ]
)

rownames(heatmap_matrix) <-
  heatmap_table$tf_symbol

heatmap_scaled <- scale(
  heatmap_matrix
)

heatmap_scaled[
  !is.finite(
    heatmap_scaled
  )
] <- 0

colnames(heatmap_scaled) <- c(
  "Full recovery",
  "Bootstrap recovery",
  "Bootstrap positive probability",
  "Bootstrap top-3 frequency",
  "Matched-null percentile",
  "Inflammation recovery",
  "Specificity"
)

write_csv_safe(
  data.table::as.data.table(
    heatmap_scaled,
    keep.rownames =
      "tf_symbol"
  ),
  file.path(
    DIRS$source,
    "Fig5B_E_robustness_matrix_source.csv"
  )
)

save_heatmap_bundle(
  heatmap_scaled,
  "Fig5B_E_candidate_robustness_matrix",
  10,
  7,
  title =
    "Full-regulon, bootstrap, and matched-null robustness"
)

############################################################
## 15. Workbook, methods, and parameter records
############################################################

workbook_path <- file.path(
  DIRS$tables,
  "14_stage5B_OFFLINE_bootstrap_null_key_results.xlsx"
)

workbook_sheets <- list(
  Network_integrity =
    as.data.frame(
      network_integrity
    ),
  Candidate_resolution =
    as.data.frame(
      candidate_resolution
    ),
  Program_definitions =
    as.data.frame(
      program_summary
    ),
  Full_candidate_summary =
    as.data.frame(
      full_candidate_summary
    ),
  Stage5_reproduction =
    as.data.frame(
      stage5_audit
    ),
  Bootstrap_summary =
    as.data.frame(
      bootstrap_summary
    ),
  Null_pool =
    as.data.frame(
      null_pool
    ),
  Null_summary =
    as.data.frame(
      null_summary
    ),
  Final_robustness =
    as.data.frame(
      final_robustness
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

xlsx_required_patterns <- c(
  "xl/workbook.xml",
  "xl/worksheets/sheet1.xml",
  "[Content_Types].xml"
)

xlsx_structure_ok <- all(
  vapply(
    xlsx_required_patterns,
    function(pattern_i) {
      any(
        xlsx_contents$Name ==
          pattern_i
      )
    },
    logical(1)
  )
)

if (!xlsx_structure_ok) {
  stop(
    "Generated XLSX did not pass internal structure validation."
  )
}

parameter_table <- data.table::data.table(
  parameter = c(
    "Random seed",
    "Candidate TFs",
    "External network access",
    "Primary perturbation mode",
    "Primary perturbation strength",
    "Perturbation formulation",
    "Bootstrap replicates",
    "Bootstrap target fraction",
    "Random null draws per candidate",
    "Nearest null matching pool",
    "Minimum null pool",
    "Null matching covariates",
    "Program-overlap filter in null matching",
    "Primary signature size",
    "Minimum regulon size",
    "Maximum target shift",
    "Stage 5 reproduction tolerance",
    "Inferential unit",
    "Nfkb1 forced"
  ),
  value = c(
    "20260714",
    paste(
      candidate_tfs,
      collapse = "; "
    ),
    "None; fully offline",
    "disease_normalization",
    "1.00",
    "weighted minimum-norm target adjustment",
    as.character(
      BOOTSTRAP_REPLICATES
    ),
    as.character(
      BOOTSTRAP_TARGET_FRACTION
    ),
    as.character(
      NULL_DRAWS_PER_CANDIDATE
    ),
    as.character(
      NULL_NEAREST_POOL_SIZE
    ),
    as.character(
      NULL_MIN_POOL_SIZE
    ),
    "log regulon size; absolute TF-activity effect; available absolute TF-expression effect",
    "FALSE",
    as.character(
      PRIMARY_SIGNATURE_SIZE
    ),
    as.character(
      MIN_TARGETS_PER_REGULON
    ),
    as.character(
      MAX_ABS_GENE_SHIFT_SD
    ),
    as.character(
      STAGE5_AUDIT_TOLERANCE
    ),
    "Biological sample",
    "FALSE"
  ),
  rationale = c(
    "Reproducibility",
    "Prespecified from completed Stage 4 and Stage 5",
    "Avoid unstable external APIs and preserve an auditable offline workflow",
    "Moves each HFpEF TF activity to the Control reference",
    "Matches the primary Stage 5 v2 comparison",
    "Preserves the signed Stage 4 regulon activity definition",
    "Quantify dependence on candidate target composition",
    "Retain most targets while generating meaningful target-set variation",
    "Generate a stable empirical matched background",
    "Use local covariate matching without program-overlap selection",
    "Avoid unstable small null pools",
    "Control network size and observed TF signal magnitude",
    "Prevent biased selection of low-overlap negative controls",
    "Locked primary Stage 2 program size",
    "Avoid unstable small regulons",
    "Prevent implausibly large predicted expression shifts",
    "Require numerical reproduction of Stage 5 v2",
    "Avoid cell-level pseudoreplication",
    "No candidate is manually promoted"
  )
)

write_csv_safe(
  parameter_table,
  file.path(
    DIRS$methods,
    "stage5B_OFFLINE_parameters_and_rationale.csv"
  )
)

methods_text <- c(
  "HFpEF Stage 5B OFFLINE FIXED v1",
  "DoRothEA regulon bootstrap and unbiased matched-TF null distributions",
  "",
  "Input boundary:",
  "- Stage 5B loaded the completed Stage 4 sample-level macrophage pseudobulk expression, signed DoRothEA A-C TF-target network, regulon activity matrix, and Stage 2/3 program definitions.",
  "- Stage 5B loaded the completed Stage 5 v2 candidate table only for numerical reproduction auditing.",
  "- No raw single-cell processing, clustering, annotation, differential expression, Stage 4 activity inference, or Stage 5 multi-mode perturbation was repeated.",
  "",
  "Full-regulon reproduction:",
  "- The six prespecified candidates were recalculated with the primary Stage 5 v2 disease-normalizing weighted minimum-norm perturbation.",
  "- Recalculated primary program-recovery values were required to match Stage 5 v2 within an absolute tolerance of 1e-6.",
  "",
  "Regulon bootstrap:",
  "- For each candidate, 500 replicates retained 80% of the locked DoRothEA regulon targets without replacement.",
  "- Candidate activity, virtual perturbation, program recovery, and cross-candidate rank were recalculated in every replicate.",
  "- Reported outputs include median recovery, 2.5th and 97.5th percentiles, probability of positive recovery, and top-1/top-3 frequencies.",
  "",
  "Matched random-TF null:",
  "- Null TFs were matched on log regulon size and absolute TF-activity effect.",
  "- Absolute TF-expression effect was also used when the candidate TF itself had detectable expression; missing null-TF expression received a predefined distance penalty rather than being silently excluded.",
  "- Stage 2 or Stage 3 program overlap was not used to select null TFs.",
  "- The 60 nearest eligible TFs formed each local pool and 1,000 distance-weighted draws with replacement generated each empirical null distribution.",
  "- Candidate percentiles and one-sided empirical P values compare candidate primary-program recovery with the matched background.",
  "",
  "Claim boundary:",
  "- Stage 5B is an offline computational robustness analysis, not experimental perturbation.",
  "- Bootstrap intervals quantify regulon-target dependence and are not biological population confidence intervals.",
  "- Empirical P values are candidate-prioritization statistics and do not demonstrate causal regulation.",
  "- Nfkb1 was not forced into the final rank."
)

writeLines(
  methods_text,
  file.path(
    DIRS$methods,
    "stage5B_OFFLINE_methods_and_claim_boundaries.txt"
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
## 16. Completion checks and final status
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
    "15_stage5B_warnings_and_nonfatal_issues.csv"
  )
)

minimum_null_pool_size <- min(
  null_pool[
    ,
    .N,
    by = candidate_tf
  ]$N
)

scientific_checks <- data.table::data.table(
  check = c(
    "Stage 4 completed",
    "Stage 5 v2 completed",
    "Stage 4 failed checks",
    "Stage 5 failed checks",
    "Biological samples",
    "Control samples",
    "HFpEF samples",
    "Offline external network calls",
    "DoRothEA integrity failures",
    "Candidate TFs resolved",
    "Candidate TF activities available",
    "Full candidate summaries",
    "Stage 5 reproduction rows",
    "Stage 5 reproduction failures",
    "Bootstrap rows",
    "Bootstrap replicates per candidate",
    "Minimum matched null pool",
    "Null draw rows",
    "Null summary rows",
    "Program-overlap filtering absent",
    "Final candidate rows",
    "Scientific checkpoint",
    "Workbook",
    "Workbook structure"
  ),
  observed = c(
    as.integer(
      stage4_status$
        overall_status[1L] %in%
        allowed_stage4_status
    ),
    as.integer(
      stage5_status$
        overall_status[1L] ==
        "COMPLETED_STAGE5_READY_FOR_REVIEW"
    ),
    sum(
      stage4_checks$status !=
        "PASS"
    ),
    sum(
      stage5_checks$status !=
        "PASS"
    ),
    nrow(sample_meta),
    sum(
      sample_meta$condition ==
        "Control"
    ),
    sum(
      sample_meta$condition ==
        "HFpEF"
    ),
    0L,
    sum(
      network_integrity$
        target_count_match != TRUE
    ),
    nrow(candidate_resolution),
    sum(
      candidate_resolution$
        activity_available == TRUE
    ),
    nrow(full_candidate_summary),
    nrow(stage5_audit),
    sum(
      stage5_audit$
        within_tolerance != TRUE
    ),
    nrow(bootstrap_results),
    min(
      bootstrap_results[
        ,
        .N,
        by = tf_symbol
      ]$N
    ),
    minimum_null_pool_size,
    nrow(null_draws),
    nrow(null_summary),
    sum(
      null_pool$
        program_overlap_filter_used ==
        TRUE
    ),
    nrow(final_robustness),
    as.integer(
      file.exists(
        file.path(
          DIRS$objects,
          "CHECKPOINT_stage5B_OFFLINE_scientific_results_pre_figures.rds"
        )
      )
    ),
    as.integer(
      file.exists(
        workbook_path
      )
    ),
    as.integer(
      xlsx_structure_ok
    )
  ),
  expected = c(
    1L,
    1L,
    0L,
    0L,
    6L,
    3L,
    3L,
    0L,
    0L,
    length(
      candidate_tfs
    ),
    length(
      candidate_tfs
    ),
    length(
      candidate_tfs
    ),
    length(
      candidate_tfs
    ),
    0L,
    BOOTSTRAP_REPLICATES *
      length(candidate_tfs),
    BOOTSTRAP_REPLICATES,
    NULL_MIN_POOL_SIZE,
    NULL_DRAWS_PER_CANDIDATE *
      length(candidate_tfs),
    length(
      candidate_tfs
    ),
    0L,
    length(
      candidate_tfs
    ),
    1L,
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
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "at_least",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal"
  )
)

scientific_checks[
  ,
  status := data.table::fcase(
    comparison == "equal" &
      observed == expected,
    "PASS",
    comparison == "at_least" &
      observed >= expected,
    "PASS",
    default = "FAIL"
  )
]

write_csv_safe(
  scientific_checks,
  file.path(
    DIRS$tables,
    "16_stage5B_scientific_completion_checks.csv"
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
  "COMPLETED_STAGE5B_OFFLINE_READY_FOR_REVIEW"
} else {
  "COMPLETED_STAGE5B_OFFLINE_REVIEW_REQUIRED"
}

nfkb1_rank <- final_robustness[
  gene_key(tf_symbol) ==
    "NFKB1",
  final_robustness_rank
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
  candidate_TFs =
    paste(
      candidate_tfs,
      collapse = ";"
    ),
  external_network_calls = 0L,
  bootstrap_replicates_per_candidate =
    BOOTSTRAP_REPLICATES,
  bootstrap_rows =
    nrow(bootstrap_results),
  null_draws_per_candidate =
    NULL_DRAWS_PER_CANDIDATE,
  null_draw_rows =
    nrow(null_draws),
  top_robust_candidate =
    final_robustness$
      tf_symbol[1L],
  Nfkb1_rank = if (
    length(nfkb1_rank) == 1L
  ) {
    nfkb1_rank
  } else {
    NA_integer_
  },
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
    "17_stage5B_run_status.csv"
  )
)

readme <- c(
  "HFpEF Reanalysis Project - Stage 5B OFFLINE FIXED v1",
  "DoRothEA regulon bootstrap and unbiased matched-TF null distributions",
  "",
  paste0(
    "Overall status: ",
    overall_status
  ),
  paste0(
    "Top robustness candidate: ",
    final_robustness$
      tf_symbol[1L]
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
  "Primary outputs:",
  "- Exact reproduction audit against Stage 5 v2.",
  "- 500 80%-target regulon bootstrap replicates per candidate.",
  "- 1,000 matched random-TF null draws per candidate.",
  "- Null matching without Stage 2/3 program-overlap filtering.",
  "- Transparent final rank aggregation without forcing Nfkb1.",
  "",
  "This stage is fully offline and does not use CollecTRI, decoupleR, OmnipathR, or internet access.",
  "",
  "Upload the Stage 5B OFFLINE CHECK package before Stage 6."
)

writeLines(
  readme,
  file.path(
    OUT_DIR,
    "README_stage5B_OFFLINE.txt"
  ),
  useBytes = TRUE
)

############################################################
## 17. CHECK package and hashes
############################################################

review_files <- c(
  LOG_FILE,
  file.path(
    DIRS$tables,
    "00_stage5B_network_integrity_audit.csv"
  ),
  file.path(
    DIRS$tables,
    "01_stage5B_candidate_resolution.csv"
  ),
  file.path(
    DIRS$tables,
    "02_stage5B_program_definition_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "04_stage5B_full_candidate_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "05_stage5B_stage5_v2_reproduction_audit.csv"
  ),
  file.path(
    DIRS$tables,
    "07_stage5B_candidate_regulon_bootstrap_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "08_stage5B_all_TF_matching_covariates.csv"
  ),
  file.path(
    DIRS$tables,
    "09_stage5B_candidate_matched_null_TF_pools.csv"
  ),
  file.path(
    DIRS$tables,
    "10_stage5B_precomputed_null_TF_effects.csv"
  ),
  file.path(
    DIRS$tables,
    "12_stage5B_random_matched_TF_null_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "13_stage5B_final_candidate_robustness_rank.csv"
  ),
  workbook_path,
  file.path(
    DIRS$tables,
    "15_stage5B_warnings_and_nonfatal_issues.csv"
  ),
  file.path(
    DIRS$tables,
    "16_stage5B_scientific_completion_checks.csv"
  ),
  file.path(
    DIRS$tables,
    "17_stage5B_run_status.csv"
  ),
  file.path(
    DIRS$methods,
    "stage5B_OFFLINE_parameters_and_rationale.csv"
  ),
  file.path(
    DIRS$methods,
    "stage5B_OFFLINE_methods_and_claim_boundaries.txt"
  ),
  file.path(
    DIRS$methods,
    "sessionInfo.txt"
  ),
  file.path(
    OUT_DIR,
    "README_stage5B_OFFLINE.txt"
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
  "Stage 5B OFFLINE analysis finished."
)

log_msg(
  "Overall status: ",
  overall_status
)

log_msg(
  "Top robustness candidate: ",
  final_robustness$
    tf_symbol[1L]
)

log_msg(
  "Nfkb1 rank: ",
  ifelse(
    length(nfkb1_rank) == 1L,
    nfkb1_rank,
    "not ranked"
  )
)

log_msg(
  "CHECK package: ",
  CHECK_ZIP
)

cat(
  "\n============================================================\n"
)

cat(
  "HFpEF Stage 5B OFFLINE robustness analysis completed\n"
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
  "Top robustness candidate: ",
  final_robustness$
    tf_symbol[1L],
  "\n",
  sep = ""
)

cat(
  "Nfkb1 was not forced.\n"
)

cat(
  "Upload the CHECK package before Stage 6.\n"
)

cat(
  "============================================================\n"
)
