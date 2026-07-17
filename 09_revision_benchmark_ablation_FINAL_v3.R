############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 7 FINAL v2
## Cross-stage constrained sample-level ridge classification
## and exact additive linear-predictor attribution
##
## Project:
##   <HFPEF_PROJECT_DIR>
##
## Required completed inputs:
##   Stage 3:
##     03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH
##   Stage 4:
##     04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1
##   Stage 5B:
##     05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1
##   Stage 6:
##     06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3
##
## Critical FINAL v2 correction:
##   In FINAL v1, ligand_sample was grouped by sample_accession while
##   sample_accession was also returned explicitly inside .(...).
##   data.table therefore created duplicate sample_accession columns,
##   and merge.data.table stopped with check_duplicate_names(y).
##   FINAL v2 removes that duplicated grouping column, validates unique
##   column names and one-row-per-sample keys after grouped summaries,
##   and preserves the locked six-sample order through all merges.
##
## Scientific purpose:
##   Evaluate whether a prespecified compact feature panel spanning
##   macrophage TF activity, drug-opposed transcriptional programs,
##   and TF-linked macrophage-to-vascular/stromal communication carries
##   HFpEF-versus-Control information at the biological-sample level.
##
## Primary validation:
##   - Biological sample is the only modeling unit.
##   - 3 Control and 3 HFpEF samples.
##   - Exhaustive leave-one-Control-plus-one-HFpEF-pair-out validation:
##     3 x 3 = 9 held-out pairs.
##   - Fixed ridge penalty; no outcome-driven feature selection and no
##     hyperparameter optimization on the six samples.
##   - Exact enumeration of all 20 balanced 3-versus-3 label assignments
##     for an empirical sample-label permutation null.
##
## Primary feature panel:
##   1) Stage 4 Bhlhe40 activity
##   2) Stage 4 Nfkb1 activity
##   3) Stage 4 Rela activity
##   4) Stage 2 Top150 drug-opposed macrophage program score
##   5) Stage 6 Nfkb1/Rela communication-burden score
##
## Extended sensitivity panel:
##   - Primary panel
##   - Stage 3-supported Top150 drug-opposed program score
##   - Bhlhe40 communication-burden score
##
## Interpretation boundary:
##   - This is an exploratory internal separability analysis.
##   - Six samples cannot establish a clinical classifier.
##   - Stage 4-6 features were derived from the same biological samples;
##     this is cross-stage internal evidence, not external validation.
##   - Feature attribution equals beta_j * standardized_feature_j on the
##     logistic linear predictor. It is exact additive logit attribution
##     for the fitted ridge model, not a causal effect and not a general
##     nonlinear SHAP analysis.
##   - Cell-level observations are never used as independent model rows.
##
## Output:
##   <HFPEF_PROJECT_DIR>/
##   07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2
##
## CHECK:
##   <HFPEF_PROJECT_DIR>/
##   07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2_CHECK.zip
##
## Run:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage7_CrossStage_Sample_Ridge_Attribution_FINAL_v2.R",
##     encoding = "UTF-8"
##   )
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(warn = 1)
options(encoding = "UTF-8")
options(timeout = 7200)

set.seed(20260714)

############################################################
## 0. Locked paths and prespecified settings
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

STAGE5B_DIR <- file.path(
  PROJECT_DIR,
  "05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1"
)

STAGE6_DIR <- file.path(
  PROJECT_DIR,
  "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3"
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

STAGE6_STATUS_FILE <- file.path(
  STAGE6_DIR,
  "01_tables",
  "23_stage6_run_status.csv"
)

STAGE6_CHECKS_FILE <- file.path(
  STAGE6_DIR,
  "01_tables",
  "22_stage6_scientific_completion_checks.csv"
)

STAGE6_AXES_FILE <- file.path(
  STAGE6_DIR,
  "01_tables",
  "15_stage6_candidate_TF_ligand_receptor_axes.csv"
)

STAGE6_EXPRESSION_FILE <- file.path(
  STAGE6_DIR,
  "01_tables",
  "05_stage6_selected_ligand_receptor_expression_by_sample.csv"
)

STAGE6_CANDIDATE_SUMMARY_FILE <- file.path(
  STAGE6_DIR,
  "01_tables",
  "19_stage6_candidate_TF_communication_summary.csv"
)

STAGE_NAME <-
  "07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2"

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
  "R",
  "07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2.R"
)

REPLACE_EXISTING_STAGE7 <- TRUE

CANDIDATE_TFS <- c(
  "Bhlhe40",
  "Nfkb1",
  "Rela"
)

PRIMARY_SIGNATURE_SIZE <- 150L

PRIMARY_PANEL <- c(
  "TF_Bhlhe40_activity",
  "TF_Nfkb1_activity",
  "TF_Rela_activity",
  "PROGRAM_DrugOpposed_Top150",
  "COMM_NFkB_axis_burden"
)

EXTENDED_PANEL <- c(
  PRIMARY_PANEL,
  "PROGRAM_DrugOpposed_Top150_Stage3Supported",
  "COMM_Bhlhe40_axis_burden"
)

PRIMARY_RIDGE_LAMBDA <- 1.0

RIDGE_LAMBDA_SENSITIVITY <- c(
  0.25,
  1.0,
  4.0
)

MAX_NFKB_AXES <- 30L
MAX_BHLHE40_AXES <- 20L
STRICT_SUPPORT_WEIGHT <- 1.25
RIDGE_MAX_ITERATIONS <- 5000L
RIDGE_CONVERGENCE_TOLERANCE <- 1e-8

MAX_FIGURE_WIDTH_IN <- 12
MAX_FIGURE_HEIGHT_IN <- 9

############################################################
## 1. Preflight, packages, output, and logging
############################################################

detect_script_file <- function() {
  candidates <- character()
  frames <- sys.frames()

  for (frame_index in rev(seq_along(frames))) {
    source_file <- tryCatch(
      frames[[frame_index]]$ofile,
      error = function(e) NULL
    )

    if (
      !is.null(source_file) &&
      length(source_file) == 1L &&
      nzchar(source_file)
    ) {
      candidates <- c(
        candidates,
        source_file
      )
    }
  }

  command_arguments <- commandArgs(
    trailingOnly = FALSE
  )

  file_argument <- grep(
    "^--file=",
    command_arguments,
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
  STAGE4_STATUS_FILE,
  STAGE4_CHECKS_FILE,
  STAGE4_PROGRAM_FILE,
  STAGE4_PSEUDOBULK_RDS,
  STAGE4_ACTIVITY_RDS,
  STAGE5B_STATUS_FILE,
  STAGE5B_CHECKS_FILE,
  STAGE5B_RANK_FILE,
  STAGE6_STATUS_FILE,
  STAGE6_CHECKS_FILE,
  STAGE6_AXES_FILE,
  STAGE6_EXPRESSION_FILE,
  STAGE6_CANDIDATE_SUMMARY_FILE
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs) > 0L) {
  stop(
    "Required Stage 3/4/5B/6 input path(s) are missing:\n",
    paste(
      missing_inputs,
      collapse = "\n"
    )
  )
}

ensure_cran <- function(packages) {
  missing <- packages[
    !vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing) > 0L) {
    install.packages(
      missing,
      repos = "https://cloud.r-project.org",
      dependencies = TRUE
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
    stop(
      "Required CRAN package(s) unavailable: ",
      paste(
        still_missing,
        collapse = ", "
      )
    )
  }
}

ensure_cran(
  c(
    "data.table",
    "ggplot2",
    "writexl",
    "zip",
    "digest"
  )
)

if (REPLACE_EXISTING_STAGE7) {
  if (dir.exists(OUT_DIR)) {
    unlink(
      OUT_DIR,
      recursive = TRUE,
      force = TRUE
    )
  }

  if (file.exists(CHECK_ZIP)) {
    unlink(
      CHECK_ZIP,
      force = TRUE
    )
  }
} else if (
  dir.exists(OUT_DIR) ||
  file.exists(CHECK_ZIP)
) {
  stop(
    "Existing Stage 7 output detected while replacement is disabled."
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

for (directory_i in c(
  OUT_DIR,
  unlist(
    DIRS,
    use.names = FALSE
  )
)) {
  dir.create(
    directory_i,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

START_TIME <- Sys.time()

LOG_FILE <- file.path(
  DIRS$logs,
  "stage7_cross_stage_sample_ridge.log"
)

WARN_FILE <- file.path(
  DIRS$logs,
  "stage7_warnings.log"
)

log_msg <- function(
  ...,
  level = "INFO"
) {
  line <- sprintf(
    "[%s] [%s] %s",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    level,
    paste0(
      ...,
      collapse = ""
    )
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
  "Stage 7 analysis started."
)

############################################################
## 2. General utilities
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

assert_unique_column_names <- function(
  table_object,
  object_name
) {
  column_names <- names(
    table_object
  )

  duplicated_names <- unique(
    column_names[
      duplicated(
        column_names
      )
    ]
  )

  if (length(duplicated_names) > 0L) {
    stop(
      object_name,
      " contains duplicated column name(s): ",
      paste(
        duplicated_names,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}

assert_required_columns <- function(
  table_object,
  required_columns,
  object_name
) {
  missing_columns <- setdiff(
    required_columns,
    names(
      table_object
    )
  )

  if (length(missing_columns) > 0L) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}

assert_unique_rows_by <- function(
  table_object,
  key_columns,
  object_name
) {
  assert_required_columns(
    table_object,
    key_columns,
    object_name
  )

  key_frame <- as.data.frame(
    table_object[
      ,
      key_columns,
      with = FALSE
    ],
    stringsAsFactors = FALSE
  )

  duplicated_key <- duplicated(
    key_frame
  )

  if (any(duplicated_key)) {
    duplicate_preview <- unique(
      key_frame[
        duplicated_key,
        ,
        drop = FALSE
      ]
    )

    preview_text <- paste(
      utils::capture.output(
        print(
          utils::head(
            duplicate_preview,
            5L
          ),
          row.names = FALSE
        )
      ),
      collapse = " | "
    )

    stop(
      object_name,
      " contains duplicated row key(s) for ",
      paste(
        key_columns,
        collapse = "+"
      ),
      ". Preview: ",
      preview_text
    )
  }

  invisible(TRUE)
}

grouped_summary_self_test <- data.table::data.table(
  sample_accession = c(
    "S1",
    "S1",
    "S2"
  ),
  value = c(
    1,
    3,
    2
  )
)[
  ,
  .(
    mean_value =
      mean(value)
  ),
  by = .(
    sample_accession
  )
]

assert_unique_column_names(
  grouped_summary_self_test,
  "grouped_summary_self_test"
)

assert_unique_rows_by(
  grouped_summary_self_test,
  "sample_accession",
  "grouped_summary_self_test"
)

rm(grouped_summary_self_test)

write_csv_safe <- function(
  table_object,
  path,
  compress = FALSE
) {
  table_object <- data.table::as.data.table(
    table_object
  )

  if (ncol(table_object) == 0L) {
    table_object <- data.table::data.table(
      note = "No records generated."
    )
  }

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

safe_auc <- function(
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

roc_curve_dt <- function(
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

  if (length(unique(labels)) < 2L) {
    return(
      data.table::data.table(
        false_positive_rate = c(
          0,
          1
        ),
        true_positive_rate = c(
          0,
          1
        )
      )
    )
  }

  order_index <- order(
    scores,
    decreasing = TRUE
  )

  labels <- labels[
    order_index
  ]

  scores <- scores[
    order_index
  ]

  positive_n <- sum(
    labels == 1L
  )

  negative_n <- sum(
    labels == 0L
  )

  true_positive <- cumsum(
    labels == 1L
  )

  false_positive <- cumsum(
    labels == 0L
  )

  data.table::data.table(
    false_positive_rate = c(
      0,
      false_positive /
        negative_n,
      1
    ),
    true_positive_rate = c(
      0,
      true_positive /
        positive_n,
      1
    ),
    threshold = c(
      Inf,
      scores,
      -Inf
    )
  )
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

map_keys_to_features <- function(
  keys,
  feature_map
) {
  query <- data.table::data.table(
    feature_key = gene_key(keys)
  )

  unique(
    merge(
      query,
      feature_map,
      by = "feature_key",
      all.x = TRUE,
      sort = FALSE
    )[
      !is.na(feature),
      feature
    ]
  )
}

split_semicolon <- function(x) {
  x <- normalize_text(x)

  parts <- unlist(
    strsplit(
      x,
      ";",
      fixed = TRUE
    ),
    use.names = FALSE
  )

  parts <- trimws(parts)

  unique(
    parts[
      nzchar(parts)
    ]
  )
}

safe_scale_training <- function(
  train_matrix,
  test_matrix
) {
  train_matrix <- as.matrix(
    train_matrix
  )

  test_matrix <- as.matrix(
    test_matrix
  )

  training_mean <- colMeans(
    train_matrix,
    na.rm = TRUE
  )

  training_sd <- apply(
    train_matrix,
    2L,
    stats::sd,
    na.rm = TRUE
  )

  training_mean[
    !is.finite(training_mean)
  ] <- 0

  training_sd[
    !is.finite(training_sd) |
      training_sd == 0
  ] <- 1

  for (
    column_index in seq_len(
      ncol(train_matrix)
    )
  ) {
    train_missing <- !is.finite(
      train_matrix[
        ,
        column_index
      ]
    )

    test_missing <- !is.finite(
      test_matrix[
        ,
        column_index
      ]
    )

    train_matrix[
      train_missing,
      column_index
    ] <- training_mean[
      column_index
    ]

    test_matrix[
      test_missing,
      column_index
    ] <- training_mean[
      column_index
    ]
  }

  train_z <- sweep(
    train_matrix,
    2L,
    training_mean,
    "-"
  )

  train_z <- sweep(
    train_z,
    2L,
    training_sd,
    "/"
  )

  test_z <- sweep(
    test_matrix,
    2L,
    training_mean,
    "-"
  )

  test_z <- sweep(
    test_z,
    2L,
    training_sd,
    "/"
  )

  list(
    train_z = train_z,
    test_z = test_z,
    mean = training_mean,
    sd = training_sd
  )
}

expit <- function(x) {
  x <- pmax(
    pmin(
      as.numeric(x),
      35
    ),
    -35
  )

  1 /
    (
      1 +
        exp(-x)
    )
}

fit_ridge_logistic <- function(
  x,
  y,
  lambda
) {
  x <- as.matrix(x)
  y <- as.numeric(y)

  if (
    nrow(x) != length(y) ||
    length(unique(y)) != 2L
  ) {
    stop(
      "Ridge logistic fit requires aligned rows and both classes."
    )
  }

  objective <- function(parameters) {
    intercept <- parameters[1L]
    coefficients <- parameters[-1L]

    linear_predictor <-
      intercept +
      as.numeric(
        x %*% coefficients
      )

    probability <- expit(
      linear_predictor
    )

    negative_log_likelihood <-
      -sum(
        y * log(
          pmax(
            probability,
            1e-15
          )
        ) +
          (
            1 -
              y
          ) *
          log(
            pmax(
              1 -
                probability,
              1e-15
            )
          )
      )

    ridge_penalty <-
      0.5 *
      lambda *
      sum(
        coefficients^2
      )

    negative_log_likelihood +
      ridge_penalty
  }

  gradient <- function(parameters) {
    intercept <- parameters[1L]
    coefficients <- parameters[-1L]

    linear_predictor <-
      intercept +
      as.numeric(
        x %*% coefficients
      )

    probability <- expit(
      linear_predictor
    )

    residual <- probability - y

    c(
      sum(residual),
      as.numeric(
        crossprod(
          x,
          residual
        )
      ) +
        lambda *
        coefficients
    )
  }

  initial_intercept <- stats::qlogis(
    min(
      max(
        mean(y),
        0.05
      ),
      0.95
    )
  )

  initial_parameters <- c(
    initial_intercept,
    rep(
      0,
      ncol(x)
    )
  )

  optimization_attempts <- list()

  optimization_attempts[["BFGS"]] <- tryCatch(
    stats::optim(
      par = initial_parameters,
      fn = objective,
      gr = gradient,
      method = "BFGS",
      control = list(
        maxit =
          RIDGE_MAX_ITERATIONS,
        reltol =
          RIDGE_CONVERGENCE_TOLERANCE
      ),
      hessian = FALSE
    ),
    error = function(e) NULL
  )

  bfgs_valid <- (
    !is.null(
      optimization_attempts[["BFGS"]]
    ) &&
      length(
        optimization_attempts[["BFGS"]]$par
      ) ==
        ncol(x) + 1L &&
      all(
        is.finite(
          optimization_attempts[["BFGS"]]$par
        )
      )
  )

  if (
    !bfgs_valid ||
    optimization_attempts[["BFGS"]]$convergence !=
      0L
  ) {
    start_parameters <- if (
      bfgs_valid
    ) {
      optimization_attempts[["BFGS"]]$par
    } else {
      initial_parameters
    }

    optimization_attempts[["L-BFGS-B"]] <- tryCatch(
      stats::optim(
        par = start_parameters,
        fn = objective,
        gr = gradient,
        method = "L-BFGS-B",
        control = list(
          maxit =
            RIDGE_MAX_ITERATIONS,
          factr = 1e7,
          pgtol =
            RIDGE_CONVERGENCE_TOLERANCE
        ),
        hessian = FALSE
      ),
      error = function(e) NULL
    )
  }

  valid_attempts <- optimization_attempts[
    vapply(
      optimization_attempts,
      function(attempt_i) {
        !is.null(attempt_i) &&
          length(attempt_i$par) ==
            ncol(x) + 1L &&
          all(
            is.finite(
              attempt_i$par
            )
          ) &&
          is.finite(
            attempt_i$value
          )
      },
      logical(1)
    )
  ]

  if (length(valid_attempts) == 0L) {
    return(
      list(
        intercept =
          initial_intercept,
        coefficients =
          rep(
            0,
            ncol(x)
          ),
        converged = FALSE,
        convergence_code = NA_integer_,
        optimization_method =
          "intercept_only_fallback",
        objective = NA_real_
      )
    )
  }

  converged_attempts <- valid_attempts[
    vapply(
      valid_attempts,
      function(attempt_i) {
        attempt_i$convergence == 0L
      },
      logical(1)
    )
  ]

  selected_pool <- if (
    length(converged_attempts) >
      0L
  ) {
    converged_attempts
  } else {
    valid_attempts
  }

  selected_index <- which.min(
    vapply(
      selected_pool,
      function(attempt_i) {
        attempt_i$value
      },
      numeric(1)
    )
  )

  selected_name <- names(
    selected_pool
  )[selected_index]

  optimization <- selected_pool[[selected_index]]

  list(
    intercept =
      optimization$par[1L],
    coefficients =
      optimization$par[-1L],
    converged =
      optimization$convergence == 0L,
    convergence_code =
      optimization$convergence,
    optimization_method =
      selected_name,
    objective =
      optimization$value
  )
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
      "Expected a ggplot object for ",
      stem,
      "."
    )
  }

  width <- min(
    max(
      as.numeric(width),
      3
    ),
    MAX_FIGURE_WIDTH_IN
  )

  height <- min(
    max(
      as.numeric(height),
      3
    ),
    MAX_FIGURE_HEIGHT_IN
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

  if (
    any(
      !file.exists(paths)
    ) ||
    any(
      !is.finite(
        as.numeric(
          file.info(paths)$size
        )
      ) |
        as.numeric(
          file.info(paths)$size
        ) <= 0
    )
  ) {
    stop(
      "Figure export validation failed for ",
      stem,
      "."
    )
  }

  invisible(paths)
}

############################################################
## 3. Lock upstream completion states
############################################################

stage3_status <- data.table::fread(
  STAGE3_STATUS_FILE,
  encoding = "UTF-8"
)

stage4_status <- data.table::fread(
  STAGE4_STATUS_FILE,
  encoding = "UTF-8"
)

stage5b_status <- data.table::fread(
  STAGE5B_STATUS_FILE,
  encoding = "UTF-8"
)

stage6_status <- data.table::fread(
  STAGE6_STATUS_FILE,
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

stage5b_checks <- data.table::fread(
  STAGE5B_CHECKS_FILE,
  encoding = "UTF-8"
)

stage6_checks <- data.table::fread(
  STAGE6_CHECKS_FILE,
  encoding = "UTF-8"
)

expected_status <- c(
  Stage3 =
    "COMPLETED_STAGE3_READY_FOR_REVIEW",
  Stage4 =
    "COMPLETED_STAGE4_READY_FOR_REVIEW",
  Stage5B =
    "COMPLETED_STAGE5B_OFFLINE_READY_FOR_REVIEW",
  Stage6 =
    "COMPLETED_STAGE6_READY_FOR_REVIEW"
)

observed_status <- c(
  Stage3 =
    stage3_status$overall_status[1L],
  Stage4 =
    stage4_status$overall_status[1L],
  Stage5B =
    stage5b_status$overall_status[1L],
  Stage6 =
    stage6_status$overall_status[1L]
)

if (
  !identical(
    unname(observed_status),
    unname(expected_status)
  )
) {
  stop(
    "One or more upstream stages are not in the required completed state:\n",
    paste(
      names(observed_status),
      observed_status,
      sep = "=",
      collapse = "\n"
    )
  )
}

for (check_table in list(
  stage3_checks,
  stage4_checks,
  stage5b_checks,
  stage6_checks
)) {
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
      "At least one upstream scientific completion check is not PASS."
    )
  }
}

upstream_status_audit <- data.table::data.table(
  stage = names(
    observed_status
  ),
  observed_status =
    unname(
      observed_status
    ),
  expected_status =
    unname(
      expected_status
    ),
  status_match = (
    unname(
      observed_status
    ) ==
      unname(
        expected_status
      )
  )
)

write_csv_safe(
  upstream_status_audit,
  file.path(
    DIRS$tables,
    "00_stage7_upstream_status_audit.csv"
  )
)

############################################################
## 4. Load sample-level upstream data
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
    "Stage 7 requires the locked 3-Control plus 3-HFpEF sample design."
  )
}

sample_order <- sample_meta$
  sample_accession

pseudobulk_objects <- readRDS(
  STAGE4_PSEUDOBULK_RDS
)

if (
  !"sample_logcpm" %in%
    names(
      pseudobulk_objects
    )
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

if (
  !all(
    sample_order %in%
      colnames(
        sample_logcpm
      )
  ) ||
  !all(
    sample_order %in%
      colnames(
        stage4_activity
      )
  )
) {
  stop(
    "Stage 4 sample matrices do not contain all six locked samples."
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

program_manifest <- read_table_auto(
  STAGE4_PROGRAM_FILE
)

stage5b_rank <- data.table::fread(
  STAGE5B_RANK_FILE,
  encoding = "UTF-8"
)

stage6_axes <- read_table_auto(
  STAGE6_AXES_FILE
)

stage6_expression <- read_table_auto(
  STAGE6_EXPRESSION_FILE
)

stage6_candidate_summary <- data.table::fread(
  STAGE6_CANDIDATE_SUMMARY_FILE,
  encoding = "UTF-8"
)

for (
  object_name_i in c(
    "sample_meta",
    "program_manifest",
    "stage5b_rank",
    "stage6_axes",
    "stage6_expression",
    "stage6_candidate_summary"
  )
) {
  assert_unique_column_names(
    get(
      object_name_i,
      inherits = FALSE
    ),
    object_name_i
  )
}

############################################################
## 5. Candidate and input-column validation
############################################################

candidate_resolution <- data.table::data.table(
  requested_tf =
    CANDIDATE_TFS,
  requested_key =
    gene_key(
      CANDIDATE_TFS
    ),
  requested_order =
    seq_along(
      CANDIDATE_TFS
    )
)

activity_resolution <- data.table::data.table(
  activity_symbol =
    rownames(
      stage4_activity
    ),
  requested_key =
    gene_key(
      rownames(
        stage4_activity
      )
    )
)

candidate_resolution <- merge(
  candidate_resolution,
  activity_resolution,
  by = "requested_key",
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
    !is.na(activity_symbol) &
      activity_symbol %in%
        rownames(
          stage4_activity
        )
  )
]

assert_unique_column_names(
  candidate_resolution,
  "candidate_resolution"
)

assert_unique_rows_by(
  candidate_resolution,
  "requested_key",
  "candidate_resolution"
)

if (
  nrow(candidate_resolution) !=
    length(
      CANDIDATE_TFS
    ) ||
  any(
    candidate_resolution$
      activity_available != TRUE
  )
) {
  stop(
    "One or more Stage 7 candidate TF activities are unavailable."
  )
}

required_program_columns <- c(
  "program_name",
  "direction",
  "signature_size",
  "symbol_key"
)

missing_program_columns <- setdiff(
  required_program_columns,
  names(program_manifest)
)

if (length(missing_program_columns) > 0L) {
  stop(
    "Stage 4 program manifest is missing column(s): ",
    paste(
      missing_program_columns,
      collapse = ", "
    )
  )
}

required_axis_columns <- c(
  "tf_symbol",
  "final_axis_rank",
  "final_axis_score",
  "sender_ligand_feature",
  "receptor_component_features",
  "receiver",
  "mean_delta_z_HFpEF",
  "predicted_receiver_direction",
  "strict_cross_stage_support",
  "axis_id"
)

missing_axis_columns <- setdiff(
  required_axis_columns,
  names(stage6_axes)
)

if (length(missing_axis_columns) > 0L) {
  stop(
    "Stage 6 axis table is missing column(s): ",
    paste(
      missing_axis_columns,
      collapse = ", "
    )
  )
}

required_expression_columns <- c(
  "feature",
  "feature_key",
  "major_cell_type",
  "sample_accession",
  "condition",
  "log2_cpm",
  "pct_expressed"
)

missing_expression_columns <- setdiff(
  required_expression_columns,
  names(stage6_expression)
)

if (length(missing_expression_columns) > 0L) {
  stop(
    "Stage 6 sample-level expression table is missing column(s): ",
    paste(
      missing_expression_columns,
      collapse = ", "
    )
  )
}

required_stage5b_candidate_columns <- c(
  "tf_symbol",
  "final_robustness_rank",
  "final_robustness_score",
  "positive_recovery_probability",
  "candidate_percentile"
)

missing_stage5b_candidate_columns <- setdiff(
  required_stage5b_candidate_columns,
  names(stage5b_rank)
)

if (length(missing_stage5b_candidate_columns) > 0L) {
  stop(
    "Stage 5B candidate table is missing column(s): ",
    paste(
      missing_stage5b_candidate_columns,
      collapse = ", "
    )
  )
}

required_stage6_candidate_columns <- c(
  "tf_symbol",
  "candidate_role",
  "total_axes",
  "strict_cross_stage_axes",
  "best_axis_rank"
)

missing_stage6_candidate_columns <- setdiff(
  required_stage6_candidate_columns,
  names(stage6_candidate_summary)
)

if (length(missing_stage6_candidate_columns) > 0L) {
  stop(
    "Stage 6 candidate summary is missing column(s): ",
    paste(
      missing_stage6_candidate_columns,
      collapse = ", "
    )
  )
}

candidate_cross_stage_manifest <- merge(
  candidate_resolution,
  stage5b_rank[
    gene_key(tf_symbol) %in%
      gene_key(CANDIDATE_TFS),
    .(
      requested_key =
        gene_key(tf_symbol),
      stage5b_tf_symbol =
        tf_symbol,
      final_robustness_rank,
      final_robustness_score,
      positive_recovery_probability,
      candidate_percentile
    )
  ],
  by = "requested_key",
  all.x = TRUE
)

candidate_cross_stage_manifest <- merge(
  candidate_cross_stage_manifest,
  stage6_candidate_summary[
    gene_key(tf_symbol) %in%
      gene_key(CANDIDATE_TFS),
    .(
      requested_key =
        gene_key(tf_symbol),
      stage6_tf_symbol =
        tf_symbol,
      candidate_role,
      total_axes,
      strict_cross_stage_axes,
      best_axis_rank
    )
  ],
  by = "requested_key",
  all.x = TRUE
)

data.table::setorder(
  candidate_cross_stage_manifest,
  requested_order
)

assert_unique_column_names(
  candidate_cross_stage_manifest,
  "candidate_cross_stage_manifest"
)

assert_unique_rows_by(
  candidate_cross_stage_manifest,
  "requested_key",
  "candidate_cross_stage_manifest"
)

if (
  nrow(candidate_cross_stage_manifest) !=
    length(CANDIDATE_TFS) ||
  any(
    is.na(
      candidate_cross_stage_manifest$
        final_robustness_rank
    )
  ) ||
  any(
    is.na(
      candidate_cross_stage_manifest$
        total_axes
    )
  )
) {
  stop(
    "The prespecified Stage 7 candidates are not fully represented across Stage 5B and Stage 6."
  )
}

write_csv_safe(
  candidate_resolution,
  file.path(
    DIRS$tables,
    "01_stage7_candidate_TF_resolution.csv"
  )
)

write_csv_safe(
  candidate_cross_stage_manifest,
  file.path(
    DIRS$tables,
    "01A_stage7_candidate_cross_stage_manifest.csv"
  )
)

############################################################
## 6. Build prespecified cross-stage features
############################################################

sample_features <- data.table::data.table(
  sample_accession =
    sample_order,
  sample_order_index =
    seq_along(
      sample_order
    )
)

sample_features <- merge(
  sample_features,
  sample_meta[
    ,
    .(
      sample_accession,
      condition
    )
  ],
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

data.table::setorder(
  sample_features,
  sample_order_index
)

sample_features[
  ,
  y := ifelse(
    condition ==
      "HFpEF",
    1L,
    0L
  )
]

assert_unique_column_names(
  sample_features,
  "sample_features_initial"
)

assert_unique_rows_by(
  sample_features,
  "sample_accession",
  "sample_features_initial"
)

feature_definition_records <- list()

for (candidate_i in CANDIDATE_TFS) {
  activity_symbol_i <- candidate_resolution[
    requested_tf ==
      candidate_i,
    activity_symbol
  ]

  feature_name <- paste0(
    "TF_",
    candidate_i,
    "_activity"
  )

  sample_features[
    ,
    (feature_name) :=
      as.numeric(
        stage4_activity[
          activity_symbol_i,
          sample_order
        ]
      )
  ]

  feature_definition_records[[length(feature_definition_records) + 1L]] <- data.table::data.table(
    feature = feature_name,
    biological_layer =
      "Stage4_regulon_activity",
    definition = paste0(
      "Locked Stage 4 weighted regulon activity for ",
      candidate_i,
      "."
    ),
    source_items =
      activity_symbol_i,
    prespecified_primary = (
      feature_name %in%
        PRIMARY_PANEL
    ),
    prespecified_extended = (
      feature_name %in%
        EXTENDED_PANEL
    )
  )
}

feature_map <- make_feature_map(
  rownames(
    sample_logcpm
  )
)

calculate_program_score <- function(
  support_mode
) {
  program_i <- data.table::copy(
    program_manifest[
      signature_size ==
        PRIMARY_SIGNATURE_SIZE
    ]
  )

  program_i[
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

  if (
    support_mode ==
      "Stage3_supported"
  ) {
    program_i <- program_i[
      support_class ==
        "Stage3_supported"
    ]
  } else {
    program_i <- program_i[
      support_class ==
        "Full_Stage2"
    ]
  }

  up_features <- map_keys_to_features(
    program_i[
      direction ==
        "Disease_up_Drug_down",
      symbol_key
    ],
    feature_map
  )

  down_features <- map_keys_to_features(
    program_i[
      direction ==
        "Disease_down_Drug_up",
      symbol_key
    ],
    feature_map
  )

  if (
    length(up_features) +
      length(down_features) <
      5L
  ) {
    stop(
      "Insufficient detected genes for ",
      support_mode,
      " Top150 drug-opposed program."
    )
  }

  up_score <- if (
    length(up_features) > 0L
  ) {
    colMeans(
      sample_logcpm[
        up_features,
        sample_order,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  } else {
    rep(
      0,
      length(sample_order)
    )
  }

  down_score <- if (
    length(down_features) > 0L
  ) {
    colMeans(
      sample_logcpm[
        down_features,
        sample_order,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  } else {
    rep(
      0,
      length(sample_order)
    )
  }

  list(
    score =
      as.numeric(
        up_score -
          down_score
      ),
    up_features =
      up_features,
    down_features =
      down_features
  )
}

program_full <- calculate_program_score(
  "Full_Stage2"
)

program_supported <- calculate_program_score(
  "Stage3_supported"
)

sample_features[
  ,
  PROGRAM_DrugOpposed_Top150 :=
    program_full$score
]

sample_features[
  ,
  PROGRAM_DrugOpposed_Top150_Stage3Supported :=
    program_supported$score
]

feature_definition_records[[length(feature_definition_records) + 1L]] <- data.table::data.table(
  feature =
    "PROGRAM_DrugOpposed_Top150",
  biological_layer =
    "Stage2_Stage4_drug_opposed_program",
  definition =
    "Mean macrophage log2-CPM of Top150 disease-up/drug-down genes minus disease-down/drug-up genes using the Full Stage 2 program.",
  source_items = paste0(
    "up=",
    length(
      program_full$up_features
    ),
    ";down=",
    length(
      program_full$down_features
    )
  ),
  prespecified_primary = TRUE,
  prespecified_extended = TRUE
)

feature_definition_records[[length(feature_definition_records) + 1L]] <- data.table::data.table(
  feature =
    "PROGRAM_DrugOpposed_Top150_Stage3Supported",
  biological_layer =
    "Stage2_Stage3_Stage4_supported_program",
  definition =
    "Mean macrophage log2-CPM of Stage3-supported Top150 disease-up/drug-down genes minus disease-down/drug-up genes.",
  source_items = paste0(
    "up=",
    length(
      program_supported$up_features
    ),
    ";down=",
    length(
      program_supported$down_features
    )
  ),
  prespecified_primary = FALSE,
  prespecified_extended = TRUE
)

build_axis_sample_scores <- function(
  axes_table,
  tf_filter,
  maximum_axes,
  feature_name
) {
  axes_i <- data.table::copy(
    axes_table[
      gene_key(tf_symbol) %in%
        gene_key(tf_filter)
    ]
  )

  data.table::setorder(
    axes_i,
    final_axis_rank,
    final_axis_score,
    axis_id
  )

  axes_i <- unique(
    axes_i,
    by = "axis_id"
  )

  axes_i <- head(
    axes_i,
    min(
      maximum_axes,
      nrow(axes_i)
    )
  )

  if (nrow(axes_i) == 0L) {
    add_warning(
      "COMMUNICATION_FEATURE",
      feature_name,
      "No Stage 6 axes were available; feature values were set to NA."
    )

    return(
      list(
        score = rep(
          NA_real_,
          length(sample_order)
        ),
        axis_sample_table =
          data.table::data.table(),
        axes = axes_i
      )
    )
  }

  axis_sample_records <- list()

  for (
    axis_index in seq_len(
      nrow(axes_i)
    )
  ) {
    axis_row <- axes_i[
      axis_index
    ]

    ligand_feature <-
      axis_row$
        sender_ligand_feature[1L]

    receptor_features <- split_semicolon(
      axis_row$
        receptor_component_features[1L]
    )

    ligand_sign <- if (
      axis_row$
        mean_delta_z_HFpEF[1L] <
        0
    ) {
      1
    } else if (
      axis_row$
        mean_delta_z_HFpEF[1L] >
        0
    ) {
      -1
    } else {
      0
    }

    receiver_sign <- if (
      axis_row$
        predicted_receiver_direction[1L] ==
        "HFpEF_up"
    ) {
      1
    } else if (
      axis_row$
        predicted_receiver_direction[1L] ==
        "HFpEF_down"
    ) {
      -1
    } else {
      0
    }

    axis_weight <- (
      1 /
        sqrt(
          max(
            as.numeric(
              axis_row$
                final_axis_rank[1L]
            ),
            1
          )
        )
    ) *
      if (
        isTRUE(
          axis_row$
            strict_cross_stage_support[1L]
        )
      ) {
        STRICT_SUPPORT_WEIGHT
      } else {
        1
      }

    ligand_sample <- stage6_expression[
      major_cell_type ==
        "Macrophage_Monocyte" &
      sample_accession %in%
        sample_order &
      feature_key ==
        gene_key(
          ligand_feature
        ),
      .(
        ligand_log2_cpm =
          safe_mean(
            log2_cpm
          ),
        ligand_pct =
          safe_mean(
            pct_expressed
          )
      ),
      by = .(
        sample_accession
      )
    ]

    assert_unique_column_names(
      ligand_sample,
      paste0(
        "ligand_sample__",
        axis_row$axis_id[1L]
      )
    )

    assert_unique_rows_by(
      ligand_sample,
      "sample_accession",
      paste0(
        "ligand_sample__",
        axis_row$axis_id[1L]
      )
    )

    receptor_sample <- stage6_expression[
      major_cell_type ==
        axis_row$
          receiver[1L] &
      sample_accession %in%
        sample_order &
      feature_key %in%
        gene_key(
          receptor_features
        ),
      .(
        receptor_components_detected =
          data.table::uniqueN(
            feature_key
          ),
        receptor_log2_cpm =
          safe_mean(
            log2_cpm
          ),
        receptor_pct =
          safe_mean(
            pct_expressed
          )
      ),
      by = .(
        sample_accession
      )
    ]

    assert_unique_column_names(
      receptor_sample,
      paste0(
        "receptor_sample__",
        axis_row$axis_id[1L]
      )
    )

    assert_unique_rows_by(
      receptor_sample,
      "sample_accession",
      paste0(
        "receptor_sample__",
        axis_row$axis_id[1L]
      )
    )

    paired <- merge(
      data.table::data.table(
        sample_accession =
          sample_order,
        sample_order_index =
          seq_along(
            sample_order
          )
      ),
      ligand_sample,
      by = "sample_accession",
      all.x = TRUE,
      sort = FALSE
    )

    assert_unique_column_names(
      paired,
      paste0(
        "paired_after_ligand__",
        axis_row$axis_id[1L]
      )
    )

    assert_unique_rows_by(
      paired,
      "sample_accession",
      paste0(
        "paired_after_ligand__",
        axis_row$axis_id[1L]
      )
    )

    paired <- merge(
      paired,
      receptor_sample,
      by = "sample_accession",
      all.x = TRUE,
      sort = FALSE
    )

    data.table::setorder(
      paired,
      sample_order_index
    )

    assert_unique_column_names(
      paired,
      paste0(
        "paired_after_receptor__",
        axis_row$axis_id[1L]
      )
    )

    assert_unique_rows_by(
      paired,
      "sample_accession",
      paste0(
        "paired_after_receptor__",
        axis_row$axis_id[1L]
      )
    )

    paired[
      ,
      `:=`(
        tf_symbol =
          axis_row$
            tf_symbol[1L],
        axis_id =
          axis_row$
            axis_id[1L],
        final_axis_rank =
          axis_row$
            final_axis_rank[1L],
        receiver =
          axis_row$
            receiver[1L],
        ligand =
          ligand_feature,
        receptor =
          axis_row$
            receptor[1L],
        ligand_sign =
          ligand_sign,
        receiver_sign =
          receiver_sign,
        axis_weight =
          axis_weight,
        receptor_components_required =
          length(
            receptor_features
          )
      )
    ]

    paired[
      ,
      axis_raw_score :=
        0.5 *
        (
          ligand_sign *
            ligand_log2_cpm +
          receiver_sign *
            receptor_log2_cpm
        )
    ]

    paired[
      receptor_components_detected <
        receptor_components_required,
      axis_raw_score := NA_real_
    ]

    axis_sample_records[[length(axis_sample_records) + 1L]] <- paired
  }

  axis_sample_table <- data.table::rbindlist(
    axis_sample_records,
    use.names = TRUE,
    fill = TRUE
  )

  assert_unique_column_names(
    axis_sample_table,
    paste0(
      feature_name,
      "_axis_sample_table"
    )
  )

  assert_unique_rows_by(
    axis_sample_table,
    c(
      "axis_id",
      "sample_accession"
    ),
    paste0(
      feature_name,
      "_axis_sample_table"
    )
  )

  aggregated <- axis_sample_table[
    is.finite(axis_raw_score) &
      is.finite(axis_weight) &
      axis_weight > 0,
    .(
      communication_score =
        sum(
          axis_weight *
            axis_raw_score
        ) /
        sum(axis_weight),
      axes_contributing =
        data.table::uniqueN(
          axis_id
        )
    ),
    by = sample_accession
  ]

  score_table <- merge(
    data.table::data.table(
      sample_accession =
        sample_order,
      sample_order_index =
        seq_along(
          sample_order
        )
    ),
    aggregated,
    by = "sample_accession",
    all.x = TRUE,
    sort = FALSE
  )

  data.table::setorder(
    score_table,
    sample_order_index
  )

  assert_unique_column_names(
    score_table,
    paste0(
      feature_name,
      "_score_table"
    )
  )

  assert_unique_rows_by(
    score_table,
    "sample_accession",
    paste0(
      feature_name,
      "_score_table"
    )
  )

  if (
    !identical(
      as.character(
        score_table$sample_accession
      ),
      as.character(
        sample_order
      )
    )
  ) {
    stop(
      feature_name,
      " score table is not aligned to the locked sample order."
    )
  }

  list(
    score =
      score_table$
        communication_score,
    axes_contributing =
      score_table$
        axes_contributing,
    axis_sample_table =
      axis_sample_table,
    axes =
      axes_i
  )
}

nfkb_communication <- build_axis_sample_scores(
  axes_table =
    stage6_axes,
  tf_filter = c(
    "Nfkb1",
    "Rela"
  ),
  maximum_axes =
    MAX_NFKB_AXES,
  feature_name =
    "COMM_NFkB_axis_burden"
)

bhlhe40_communication <- build_axis_sample_scores(
  axes_table =
    stage6_axes,
  tf_filter =
    "Bhlhe40",
  maximum_axes =
    MAX_BHLHE40_AXES,
  feature_name =
    "COMM_Bhlhe40_axis_burden"
)

sample_features[
  ,
  COMM_NFkB_axis_burden :=
    nfkb_communication$score
]

sample_features[
  ,
  COMM_Bhlhe40_axis_burden :=
    bhlhe40_communication$score
]

feature_definition_records[[length(feature_definition_records) + 1L]] <- data.table::data.table(
  feature =
    "COMM_NFkB_axis_burden",
  biological_layer =
    "Stage6_NFKB1_RELA_communication",
  definition =
    "Inverse-rank-weighted mean of oriented macrophage ligand and receiver receptor-component log2-CPM across the top Nfkb1/Rela communication axes.",
  source_items = paste(
    nfkb_communication$
      axes$axis_id,
    collapse = ";"
  ),
  prespecified_primary = TRUE,
  prespecified_extended = TRUE
)

feature_definition_records[[length(feature_definition_records) + 1L]] <- data.table::data.table(
  feature =
    "COMM_Bhlhe40_axis_burden",
  biological_layer =
    "Stage6_BHLHE40_communication",
  definition =
    "Inverse-rank-weighted mean of oriented macrophage ligand and receiver receptor-component log2-CPM across Bhlhe40 communication axes.",
  source_items = paste(
    bhlhe40_communication$
      axes$axis_id,
    collapse = ";"
  ),
  prespecified_primary = FALSE,
  prespecified_extended = TRUE
)

feature_definitions <- data.table::rbindlist(
  feature_definition_records,
  use.names = TRUE,
  fill = TRUE
)

missing_primary_features <- setdiff(
  PRIMARY_PANEL,
  names(sample_features)
)

missing_extended_features <- setdiff(
  EXTENDED_PANEL,
  names(sample_features)
)

if (
  length(missing_primary_features) > 0L ||
  length(missing_extended_features) > 0L
) {
  stop(
    "Stage 7 feature construction is incomplete. Missing primary: ",
    paste(
      missing_primary_features,
      collapse = ";"
    ),
    " | Missing extended: ",
    paste(
      missing_extended_features,
      collapse = ";"
    )
  )
}

assert_unique_column_names(
  sample_features,
  "sample_features_complete"
)

assert_unique_rows_by(
  sample_features,
  "sample_accession",
  "sample_features_complete"
)

if (
  !identical(
    as.character(
      sample_features$sample_accession
    ),
    as.character(
      sample_order
    )
  )
) {
  stop(
    "Completed Stage 7 sample feature matrix is not aligned to the locked sample order."
  )
}

for (feature_i in EXTENDED_PANEL) {
  values_i <- as.numeric(
    sample_features[[feature_i]]
  )

  finite_n <- sum(
    is.finite(values_i)
  )

  variance_i <- stats::var(
    values_i,
    na.rm = TRUE
  )

  if (
    finite_n < 4L ||
    !is.finite(variance_i) ||
    variance_i <= 0
  ) {
    stop(
      "Stage 7 feature is not model-eligible: ",
      feature_i,
      " | finite samples=",
      finite_n,
      " | variance=",
      variance_i
    )
  }
}

write_csv_safe(
  feature_definitions,
  file.path(
    DIRS$tables,
    "02_stage7_feature_definitions.csv"
  )
)

write_csv_safe(
  sample_features,
  file.path(
    DIRS$tables,
    "03_stage7_sample_level_feature_matrix.csv"
  )
)

write_csv_safe(
  nfkb_communication$
    axis_sample_table,
  file.path(
    DIRS$tables,
    "04_stage7_NFkB_axis_sample_components.csv"
  ),
  compress = TRUE
)

write_csv_safe(
  bhlhe40_communication$
    axis_sample_table,
  file.path(
    DIRS$tables,
    "05_stage7_Bhlhe40_axis_sample_components.csv"
  ),
  compress = TRUE
)

############################################################
## 7. Exhaustive leave-pair-out ridge validation
############################################################

run_leave_pair_out <- function(
  feature_table,
  feature_columns,
  labels,
  lambda,
  panel_name,
  store_details = TRUE
) {
  labels <- as.integer(labels)

  if (
    length(labels) !=
      nrow(feature_table) ||
    sum(labels == 0L) != 3L ||
    sum(labels == 1L) != 3L
  ) {
    stop(
      "Leave-pair-out requires exactly three samples per class."
    )
  }

  control_indices <- which(
    labels == 0L
  )

  hfpef_indices <- which(
    labels == 1L
  )

  pair_grid <- data.table::CJ(
    control_index =
      control_indices,
    hfpef_index =
      hfpef_indices,
    unique = TRUE
  )

  prediction_records <- list()
  coefficient_records <- list()
  contribution_records <- list()
  fold_records <- list()

  feature_matrix <- as.matrix(
    feature_table[
      ,
      feature_columns,
      with = FALSE
    ]
  )

  rownames(feature_matrix) <-
    feature_table$
      sample_accession

  for (
    fold_index in seq_len(
      nrow(pair_grid)
    )
  ) {
    test_indices <- c(
      pair_grid$
        control_index[fold_index],
      pair_grid$
        hfpef_index[fold_index]
    )

    train_indices <- setdiff(
      seq_len(
        nrow(feature_table)
      ),
      test_indices
    )

    scaled <- safe_scale_training(
      train_matrix =
        feature_matrix[
          train_indices,
          ,
          drop = FALSE
        ],
      test_matrix =
        feature_matrix[
          test_indices,
          ,
          drop = FALSE
        ]
    )

    fit <- fit_ridge_logistic(
      x =
        scaled$train_z,
      y =
        labels[
          train_indices
        ],
      lambda =
        lambda
    )

    names(
      fit$coefficients
    ) <- feature_columns

    test_linear_predictor <-
      fit$intercept +
      as.numeric(
        scaled$test_z %*%
          fit$coefficients
      )

    test_probability <- expit(
      test_linear_predictor
    )

    contribution_matrix <- sweep(
      scaled$test_z,
      2L,
      fit$coefficients,
      "*"
    )

    fold_id <- paste0(
      panel_name,
      "__lambda_",
      format(
        lambda,
        scientific = FALSE,
        trim = TRUE
      ),
      "__fold_",
      fold_index
    )

    prediction_i <- data.table::data.table(
      fold_id =
        fold_id,
      panel =
        panel_name,
      lambda =
        lambda,
      sample_accession =
        feature_table$
          sample_accession[
            test_indices
          ],
      true_label =
        labels[
          test_indices
        ],
      true_condition = ifelse(
        labels[
          test_indices
        ] == 1L,
        "HFpEF",
        "Control"
      ),
      predicted_probability =
        test_probability,
      linear_predictor =
        test_linear_predictor,
      heldout_control_sample =
        feature_table$
          sample_accession[
            pair_grid$
              control_index[
                fold_index
              ]
          ],
      heldout_hfpef_sample =
        feature_table$
          sample_accession[
            pair_grid$
              hfpef_index[
                fold_index
              ]
          ]
    )

    prediction_records[[length(prediction_records) + 1L]] <- prediction_i

    pair_control_probability <-
      prediction_i[
        true_label ==
          0L,
        predicted_probability
      ]

    pair_hfpef_probability <-
      prediction_i[
        true_label ==
          1L,
        predicted_probability
      ]

    fold_records[[length(fold_records) + 1L]] <- data.table::data.table(
      fold_id =
        fold_id,
      panel =
        panel_name,
      lambda =
        lambda,
      heldout_control_sample =
        prediction_i$
          heldout_control_sample[1L],
      heldout_hfpef_sample =
        prediction_i$
          heldout_hfpef_sample[1L],
      control_probability =
        pair_control_probability,
      hfpef_probability =
        pair_hfpef_probability,
      pair_margin =
        pair_hfpef_probability -
        pair_control_probability,
      pair_correct = (
        pair_hfpef_probability >
          pair_control_probability
      ),
      model_converged =
        fit$converged,
      convergence_code =
        fit$convergence_code,
      optimization_method =
        fit$optimization_method,
      objective =
        fit$objective
    )

    if (store_details) {
      coefficient_records[[length(coefficient_records) + 1L]] <- data.table::data.table(
        fold_id =
          fold_id,
        panel =
          panel_name,
        lambda =
          lambda,
        feature =
          feature_columns,
        coefficient =
          as.numeric(
            fit$coefficients
          ),
        intercept =
          fit$intercept,
        model_converged =
          fit$converged
      )

      for (
        test_row_index in seq_len(
          nrow(
            contribution_matrix
          )
        )
      ) {
        contribution_records[[length(contribution_records) + 1L]] <- data.table::data.table(
          fold_id =
            fold_id,
          panel =
            panel_name,
          lambda =
            lambda,
          sample_accession =
            feature_table$
              sample_accession[
                test_indices[
                  test_row_index
                ]
              ],
          true_label =
            labels[
              test_indices[
                test_row_index
              ]
            ],
          feature =
            feature_columns,
          standardized_value =
            as.numeric(
              scaled$test_z[
                test_row_index,
                ,
                drop = TRUE
              ]
            ),
          coefficient =
            as.numeric(
              fit$coefficients
            ),
          logit_contribution =
            as.numeric(
              contribution_matrix[
                test_row_index,
                ,
                drop = TRUE
              ]
            ),
          intercept =
            fit$intercept
        )
      }
    }
  }

  predictions <- data.table::rbindlist(
    prediction_records,
    use.names = TRUE,
    fill = TRUE
  )

  folds <- data.table::rbindlist(
    fold_records,
    use.names = TRUE,
    fill = TRUE
  )

  sample_predictions <- predictions[
    ,
    .(
      mean_predicted_probability =
        mean(
          predicted_probability
        ),
      median_predicted_probability =
        stats::median(
          predicted_probability
        ),
      mean_linear_predictor =
        mean(
          linear_predictor
        ),
      heldout_appearances = .N,
      true_label =
        unique(
          true_label
        )[1L],
      true_condition =
        unique(
          true_condition
        )[1L]
    ),
    by = sample_accession
  ]

  pairwise_auc <- mean(
    folds$pair_correct,
    na.rm = TRUE
  )

  sample_auc <- safe_auc(
    sample_predictions$
      true_label,
    sample_predictions$
      mean_predicted_probability
  )

  list(
    predictions =
      predictions,
    folds =
      folds,
    sample_predictions =
      sample_predictions,
    coefficients = if (
      store_details
    ) {
      data.table::rbindlist(
        coefficient_records,
        use.names = TRUE,
        fill = TRUE
      )
    } else {
      data.table::data.table()
    },
    contributions = if (
      store_details
    ) {
      data.table::rbindlist(
        contribution_records,
        use.names = TRUE,
        fill = TRUE
      )
    } else {
      data.table::data.table()
    },
    pairwise_auc =
      pairwise_auc,
    sample_auc =
      sample_auc,
    folds_converged =
      sum(
        folds$model_converged ==
          TRUE
      ),
    total_folds =
      nrow(folds)
  )
}

primary_result <- run_leave_pair_out(
  feature_table =
    sample_features,
  feature_columns =
    PRIMARY_PANEL,
  labels =
    sample_features$y,
  lambda =
    PRIMARY_RIDGE_LAMBDA,
  panel_name =
    "Primary_5_feature_panel",
  store_details = TRUE
)

write_csv_safe(
  primary_result$predictions,
  file.path(
    DIRS$tables,
    "06_stage7_primary_LOPO_heldout_predictions.csv"
  )
)

write_csv_safe(
  primary_result$folds,
  file.path(
    DIRS$tables,
    "07_stage7_primary_LOPO_fold_performance.csv"
  )
)

write_csv_safe(
  primary_result$sample_predictions,
  file.path(
    DIRS$tables,
    "08_stage7_primary_LOPO_sample_predictions.csv"
  )
)

write_csv_safe(
  primary_result$coefficients,
  file.path(
    DIRS$tables,
    "09_stage7_primary_LOPO_coefficients.csv"
  )
)

write_csv_safe(
  primary_result$contributions,
  file.path(
    DIRS$tables,
    "10_stage7_primary_exact_logit_contributions.csv"
  ),
  compress = TRUE
)

############################################################
## 8. Feature attribution and fold stability
############################################################

feature_importance <- primary_result$
  contributions[
    ,
    .(
      mean_logit_contribution =
        mean(
          logit_contribution
        ),
      mean_absolute_logit_contribution =
        mean(
          abs(
            logit_contribution
          )
        ),
      median_absolute_logit_contribution =
        stats::median(
          abs(
            logit_contribution
          )
        ),
      positive_contribution_fraction =
        mean(
          logit_contribution > 0
        ),
      negative_contribution_fraction =
        mean(
          logit_contribution < 0
        ),
      heldout_attributions = .N
    ),
    by = feature
  ]

coefficient_stability <- primary_result$
  coefficients[
    ,
    .(
      mean_coefficient =
        mean(
          coefficient
        ),
      median_coefficient =
        stats::median(
          coefficient
        ),
      mean_absolute_coefficient =
        mean(
          abs(
            coefficient
          )
        ),
      positive_fold_fraction =
        mean(
          coefficient > 0
        ),
      negative_fold_fraction =
        mean(
          coefficient < 0
        ),
      zero_fold_fraction =
        mean(
          coefficient == 0
        ),
      coefficient_folds = .N
    ),
    by = feature
  ]

feature_importance <- merge(
  feature_importance,
  coefficient_stability,
  by = "feature",
  all = TRUE
)

feature_importance <- merge(
  feature_importance,
  feature_definitions[
    ,
    .(
      feature,
      biological_layer,
      definition,
      source_items
    )
  ],
  by = "feature",
  all.x = TRUE
)

data.table::setorder(
  feature_importance,
  -mean_absolute_logit_contribution,
  -mean_absolute_coefficient,
  feature
)

feature_importance[
  ,
  importance_rank :=
    seq_len(.N)
]

write_csv_safe(
  feature_importance,
  file.path(
    DIRS$tables,
    "11_stage7_feature_attribution_and_stability.csv"
  )
)

############################################################
## 9. Exhaustive balanced-label permutation null
############################################################

sample_indices <- seq_len(
  nrow(
    sample_features
  )
)

hfpef_assignments <- utils::combn(
  sample_indices,
  3L,
  simplify = FALSE
)

permutation_records <- list()

observed_hfpef_indices <- which(
  sample_features$y ==
    1L
)

for (
  permutation_index in seq_along(
    hfpef_assignments
  )
) {
  permuted_labels <- rep(
    0L,
    nrow(
      sample_features
    )
  )

  permuted_labels[
    hfpef_assignments[[permutation_index]]
  ] <- 1L

  permutation_result <- run_leave_pair_out(
    feature_table =
      sample_features,
    feature_columns =
      PRIMARY_PANEL,
    labels =
      permuted_labels,
    lambda =
      PRIMARY_RIDGE_LAMBDA,
    panel_name =
      "Primary_5_feature_panel",
    store_details = FALSE
  )

  is_observed <- identical(
    sort(
      hfpef_assignments[[permutation_index]]
    ),
    sort(
      observed_hfpef_indices
    )
  )

  permutation_records[[length(permutation_records) + 1L]] <- data.table::data.table(
    permutation_id =
      permutation_index,
    hfpef_samples = paste(
      sample_features$
        sample_accession[
          hfpef_assignments[[permutation_index]]
        ],
      collapse = ";"
    ),
    is_observed_labeling =
      is_observed,
    pairwise_auc =
      permutation_result$
        pairwise_auc,
    sample_auc =
      permutation_result$
        sample_auc,
    converged_folds =
      permutation_result$
        folds_converged,
    total_folds =
      permutation_result$
        total_folds
  )
}

permutation_null <- data.table::rbindlist(
  permutation_records,
  use.names = TRUE,
  fill = TRUE
)

observed_pairwise_auc <- primary_result$
  pairwise_auc

observed_sample_auc <- primary_result$
  sample_auc

pairwise_empirical_p <- (
  1 +
    sum(
      permutation_null[
        is_observed_labeling ==
          FALSE,
        pairwise_auc
      ] >=
        observed_pairwise_auc
    )
) /
  (
    1 +
      sum(
        permutation_null$
          is_observed_labeling ==
          FALSE
      )
  )

sample_empirical_p <- (
  1 +
    sum(
      permutation_null[
        is_observed_labeling ==
          FALSE,
        sample_auc
      ] >=
        observed_sample_auc
    )
) /
  (
    1 +
      sum(
        permutation_null$
          is_observed_labeling ==
          FALSE
      )
  )

permutation_summary <- data.table::data.table(
  metric = c(
    "pairwise_leave_pair_out_AUC",
    "sample_level_AUC"
  ),
  observed = c(
    observed_pairwise_auc,
    observed_sample_auc
  ),
  null_median = c(
    safe_median(
      permutation_null[
        is_observed_labeling ==
          FALSE,
        pairwise_auc
      ]
    ),
    safe_median(
      permutation_null[
        is_observed_labeling ==
          FALSE,
        sample_auc
      ]
    )
  ),
  null_maximum = c(
    max(
      permutation_null[
        is_observed_labeling ==
          FALSE,
        pairwise_auc
      ],
      na.rm = TRUE
    ),
    max(
      permutation_null[
        is_observed_labeling ==
          FALSE,
        sample_auc
      ],
      na.rm = TRUE
    )
  ),
  empirical_p = c(
    pairwise_empirical_p,
    sample_empirical_p
  ),
  balanced_label_assignments =
    length(
      hfpef_assignments
    ),
  nonobserved_null_assignments =
    sum(
      permutation_null$
        is_observed_labeling ==
        FALSE
    )
)

write_csv_safe(
  permutation_null,
  file.path(
    DIRS$tables,
    "12_stage7_exhaustive_balanced_label_permutation_null.csv"
  )
)

write_csv_safe(
  permutation_summary,
  file.path(
    DIRS$tables,
    "13_stage7_permutation_summary.csv"
  )
)

############################################################
## 10. Prespecified panel and lambda sensitivity
############################################################

sensitivity_scenarios <- data.table::rbindlist(
  list(
    data.table::data.table(
      scenario =
        paste0(
          "Primary_lambda_",
          RIDGE_LAMBDA_SENSITIVITY
        ),
      panel =
        "Primary_5_feature_panel",
      lambda =
        RIDGE_LAMBDA_SENSITIVITY
    ),
    data.table::data.table(
      scenario =
        "Extended_lambda_1",
      panel =
        "Extended_7_feature_panel",
      lambda =
        PRIMARY_RIDGE_LAMBDA
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

sensitivity_records <- list()

for (
  scenario_index in seq_len(
    nrow(
      sensitivity_scenarios
    )
  )
) {
  scenario_i <- sensitivity_scenarios[
    scenario_index
  ]

  feature_columns_i <- if (
    scenario_i$
      panel[1L] ==
      "Primary_5_feature_panel"
  ) {
    PRIMARY_PANEL
  } else {
    EXTENDED_PANEL
  }

  result_i <- run_leave_pair_out(
    feature_table =
      sample_features,
    feature_columns =
      feature_columns_i,
    labels =
      sample_features$y,
    lambda =
      scenario_i$lambda[1L],
    panel_name =
      scenario_i$panel[1L],
    store_details = TRUE
  )

  coefficient_summary_i <- result_i$
    coefficients[
      ,
      .(
        mean_absolute_coefficient =
          mean(
            abs(
              coefficient
            )
          )
      ),
      by = feature
    ]

  top_feature_i <- if (
    nrow(
      coefficient_summary_i
    ) > 0L
  ) {
    coefficient_summary_i[
      order(
        -mean_absolute_coefficient,
        feature
      ),
      feature
    ][1L]
  } else {
    NA_character_
  }

  sensitivity_records[[length(sensitivity_records) + 1L]] <- data.table::data.table(
    scenario =
      scenario_i$
        scenario[1L],
    panel =
      scenario_i$
        panel[1L],
    lambda =
      scenario_i$
        lambda[1L],
    features =
      length(
        feature_columns_i
      ),
    pairwise_auc =
      result_i$
        pairwise_auc,
    sample_auc =
      result_i$
        sample_auc,
    converged_folds =
      result_i$
        folds_converged,
    total_folds =
      result_i$
        total_folds,
    top_absolute_coefficient_feature =
      top_feature_i
  )
}

sensitivity_summary <- data.table::rbindlist(
  sensitivity_records,
  use.names = TRUE,
  fill = TRUE
)

write_csv_safe(
  sensitivity_summary,
  file.path(
    DIRS$tables,
    "14_stage7_panel_and_lambda_sensitivity.csv"
  )
)

############################################################
## 11. Performance summary and checkpoint
############################################################

model_performance <- data.table::data.table(
  metric = c(
    "Primary pairwise leave-pair-out AUC",
    "Primary sample-level AUC",
    "Primary pairwise empirical permutation P",
    "Primary sample-level empirical permutation P",
    "Primary converged folds",
    "Primary total folds",
    "Biological samples",
    "Control samples",
    "HFpEF samples",
    "Primary features",
    "Extended features",
    "Balanced label assignments"
  ),
  value = c(
    primary_result$
      pairwise_auc,
    primary_result$
      sample_auc,
    pairwise_empirical_p,
    sample_empirical_p,
    primary_result$
      folds_converged,
    primary_result$
      total_folds,
    nrow(
      sample_features
    ),
    sum(
      sample_features$y ==
        0L
    ),
    sum(
      sample_features$y ==
        1L
    ),
    length(
      PRIMARY_PANEL
    ),
    length(
      EXTENDED_PANEL
    ),
    nrow(
      permutation_null
    )
  )
)

write_csv_safe(
  model_performance,
  file.path(
    DIRS$tables,
    "15_stage7_model_performance_summary.csv"
  )
)

saveRDS(
  list(
    sample_features =
      sample_features,
    feature_definitions =
      feature_definitions,
    primary_result =
      primary_result,
    feature_importance =
      feature_importance,
    permutation_null =
      permutation_null,
    permutation_summary =
      permutation_summary,
    sensitivity_summary =
      sensitivity_summary,
    nfkb_communication =
      nfkb_communication,
    bhlhe40_communication =
      bhlhe40_communication
  ),
  file.path(
    DIRS$objects,
    "CHECKPOINT_stage7_scientific_results_pre_figures.rds"
  ),
  compress = FALSE
)

log_msg(
  "Stage 7 scientific calculations checkpointed before figures."
)

############################################################
## 12. Figures
############################################################

workflow_data <- data.table::data.table(
  step = factor(
    c(
      "Stage 4\nTF activity",
      "Stage 2/3\nprograms",
      "Stage 6\ncommunication",
      "Sample-level\nfeature panel",
      "9-fold\nleave-pair-out",
      "Exact linear\nattribution"
    ),
    levels = c(
      "Stage 4\nTF activity",
      "Stage 2/3\nprograms",
      "Stage 6\ncommunication",
      "Sample-level\nfeature panel",
      "9-fold\nleave-pair-out",
      "Exact linear\nattribution"
    )
  ),
  x = seq_len(6L),
  y = 1
)

plot_workflow <- ggplot2::ggplot(
  workflow_data,
  ggplot2::aes(
    x = x,
    y = y,
    label = step
  )
) +
  ggplot2::geom_label(
    size = 3.2,
    label.size = 0.3,
    fill = "white"
  ) +
  ggplot2::geom_segment(
    data =
      workflow_data[
        x < 6
      ],
    ggplot2::aes(
      x =
        x + 0.38,
      xend =
        x + 0.62,
      y = y,
      yend = y
    ),
    arrow =
      grid::arrow(
        length =
          grid::unit(
            0.12,
            "inches"
          )
      ),
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0.5,
      6.5
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0.65,
      1.35
    )
  ) +
  ggplot2::labs(
    title =
      "Stage 7 cross-stage constrained sample-level workflow",
    subtitle =
      "Biological samples only; no cell-level pseudoreplication"
  ) +
  ggplot2::theme_void(
    base_size = 11
  )

write_csv_safe(
  workflow_data,
  file.path(
    DIRS$source,
    "Fig7A_workflow_source.csv"
  )
)

save_plot_bundle(
  plot_workflow,
  "Fig7A_cross_stage_sample_level_workflow",
  10.5,
  3.0
)

sample_prediction_plot_data <- merge(
  primary_result$
    sample_predictions,
  sample_meta[
    ,
    .(
      sample_accession,
      original_condition =
        condition
    )
  ],
  by = "sample_accession",
  all.x = TRUE
)

sample_prediction_plot_data[
  ,
  sample_accession := factor(
    sample_accession,
    levels =
      sample_prediction_plot_data[
        order(
          true_label,
          mean_predicted_probability
        ),
        sample_accession
      ]
  )
]

plot_sample_predictions <- ggplot2::ggplot(
  sample_prediction_plot_data,
  ggplot2::aes(
    x = sample_accession,
    y =
      mean_predicted_probability,
    fill =
      true_condition
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = 2
  ) +
  ggplot2::geom_col(
    width = 0.72
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1
    )
  ) +
  ggplot2::labs(
    title =
      "Leave-pair-out sample predictions",
    subtitle = paste0(
      "Pairwise AUC = ",
      round(
        primary_result$
          pairwise_auc,
        3
      ),
      "; sample AUC = ",
      round(
        primary_result$
          sample_auc,
        3
      )
    ),
    x = NULL,
    y =
      "Mean held-out HFpEF probability",
    fill =
      "Observed group"
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

write_csv_safe(
  sample_prediction_plot_data,
  file.path(
    DIRS$source,
    "Fig7B_sample_predictions_source.csv"
  )
)

save_plot_bundle(
  plot_sample_predictions,
  "Fig7B_leave_pair_out_sample_predictions",
  7.5,
  5.2
)

roc_data <- roc_curve_dt(
  primary_result$
    sample_predictions$
      true_label,
  primary_result$
    sample_predictions$
      mean_predicted_probability
)

plot_roc <- ggplot2::ggplot(
  roc_data,
  ggplot2::aes(
    x =
      false_positive_rate,
    y =
      true_positive_rate
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  ggplot2::geom_step(
    linewidth = 0.8
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title =
      "Sample-level ROC curve",
    subtitle = paste0(
      "Six-sample exploratory AUC = ",
      round(
        primary_result$
          sample_auc,
        3
      )
    ),
    x =
      "False-positive rate",
    y =
      "True-positive rate"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  )

write_csv_safe(
  roc_data,
  file.path(
    DIRS$source,
    "Fig7C_sample_ROC_source.csv"
  )
)

save_plot_bundle(
  plot_roc,
  "Fig7C_sample_level_ROC",
  5.2,
  5.2
)

importance_plot_data <- data.table::copy(
  feature_importance
)

importance_plot_data[
  ,
  feature := factor(
    feature,
    levels = rev(
      feature[
        order(
          mean_absolute_logit_contribution
        )
      ]
    )
  )
]

plot_importance <- ggplot2::ggplot(
  importance_plot_data,
  ggplot2::aes(
    x =
      mean_absolute_logit_contribution,
    y = feature,
    size =
      mean_absolute_coefficient,
    fill =
      positive_fold_fraction
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend =
        mean_absolute_logit_contribution,
      y = feature,
      yend = feature
    ),
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    shape = 21
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(
      0,
      1
    )
  ) +
  ggplot2::labs(
    title =
      "Exact additive linear-predictor attribution",
    subtitle =
      "Mean absolute held-out logit contribution; not a causal effect",
    x =
      "Mean absolute logit contribution",
    y = NULL,
    size =
      "Mean absolute\ncoefficient",
    fill =
      "Positive coefficient\nfold fraction"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  )

write_csv_safe(
  importance_plot_data,
  file.path(
    DIRS$source,
    "Fig7D_feature_attribution_source.csv"
  )
)

save_plot_bundle(
  plot_importance,
  "Fig7D_exact_linear_predictor_attribution",
  8.2,
  5.8
)

permutation_plot_data <- data.table::copy(
  permutation_null
)

plot_permutation <- ggplot2::ggplot(
  permutation_plot_data[
    is_observed_labeling ==
      FALSE
  ],
  ggplot2::aes(
    x = pairwise_auc
  )
) +
  ggplot2::geom_histogram(
    binwidth = 1 / 9,
    boundary = 0,
    closed = "left",
    color = "black",
    fill = "white"
  ) +
  ggplot2::geom_vline(
    xintercept =
      primary_result$
        pairwise_auc,
    linewidth = 0.8
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 1 / 9
    )
  ) +
  ggplot2::labs(
    title =
      "Exhaustive balanced-label permutation null",
    subtitle = paste0(
      "Observed pairwise AUC = ",
      round(
        primary_result$
          pairwise_auc,
        3
      ),
      "; empirical P = ",
      signif(
        pairwise_empirical_p,
        3
      )
    ),
    x =
      "Pairwise leave-pair-out AUC",
    y =
      "Balanced label assignments"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  )

write_csv_safe(
  permutation_plot_data,
  file.path(
    DIRS$source,
    "Fig7E_permutation_null_source.csv"
  )
)

save_plot_bundle(
  plot_permutation,
  "Fig7E_exhaustive_balanced_label_permutation",
  7.2,
  5.2
)

heatmap_long <- data.table::melt(
  sample_features[
    ,
    c(
      "sample_accession",
      "condition",
      PRIMARY_PANEL
    ),
    with = FALSE
  ],
  id.vars = c(
    "sample_accession",
    "condition"
  ),
  variable.name =
    "feature",
  value.name =
    "raw_value"
)

heatmap_long[
  ,
  z_value := {
    scaled <- as.numeric(
      scale(raw_value)
    )

    scaled[
      !is.finite(scaled)
    ] <- 0

    scaled
  },
  by = feature
]

heatmap_long[
  ,
  sample_accession := factor(
    sample_accession,
    levels =
      sample_meta[
        order(condition),
        sample_accession
      ]
  )
]

heatmap_long[
  ,
  feature := factor(
    feature,
    levels = rev(
      PRIMARY_PANEL
    )
  )
]

plot_heatmap <- ggplot2::ggplot(
  heatmap_long,
  ggplot2::aes(
    x = sample_accession,
    y = feature,
    fill = z_value
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  ggplot2::scale_fill_gradient2(
    low = "grey20",
    mid = "white",
    high = "grey80",
    midpoint = 0
  ) +
  ggplot2::labs(
    title =
      "Primary cross-stage feature matrix",
    subtitle =
      "Display z-scores are for visualization only; model scaling was training-fold-specific",
    x = NULL,
    y = NULL,
    fill =
      "Display z-score"
  ) +
  ggplot2::theme_bw(
    base_size = 9
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
  )

write_csv_safe(
  heatmap_long,
  file.path(
    DIRS$source,
    "Fig7F_feature_heatmap_source.csv"
  )
)

save_plot_bundle(
  plot_heatmap,
  "Fig7F_primary_cross_stage_feature_heatmap",
  8.5,
  5.5
)

############################################################
## 13. Figure audit, workbook, methods, and limitations
############################################################

expected_figure_stems <- c(
  "Fig7A_cross_stage_sample_level_workflow",
  "Fig7B_leave_pair_out_sample_predictions",
  "Fig7C_sample_level_ROC",
  "Fig7D_exact_linear_predictor_attribution",
  "Fig7E_exhaustive_balanced_label_permutation",
  "Fig7F_primary_cross_stage_feature_heatmap"
)

audit_figure <- function(stem_value) {
  extensions <- c(
    png = ".png",
    pdf = ".pdf",
    tiff = ".tiff"
  )

  paths <- file.path(
    DIRS$figures,
    paste0(
      stem_value,
      extensions
    )
  )

  names(paths) <- names(
    extensions
  )

  exists_vector <- file.exists(
    paths
  )

  size_vector <- rep(
    NA_real_,
    length(paths)
  )

  names(size_vector) <- names(
    paths
  )

  size_vector[
    exists_vector
  ] <- as.numeric(
    file.info(
      paths[
        exists_vector
      ]
    )$size
  )

  data.table::data.table(
    stem =
      stem_value,
    png_exists =
      unname(
        exists_vector["png"]
      ),
    pdf_exists =
      unname(
        exists_vector["pdf"]
      ),
    tiff_exists =
      unname(
        exists_vector["tiff"]
      ),
    png_size_bytes =
      unname(
        size_vector["png"]
      ),
    pdf_size_bytes =
      unname(
        size_vector["pdf"]
      ),
    tiff_size_bytes =
      unname(
        size_vector["tiff"]
      ),
    files_valid = (
      all(exists_vector) &&
        all(
          is.finite(
            size_vector
          ) &
            size_vector > 0
        )
    )
  )
}

figure_audit <- data.table::rbindlist(
  lapply(
    expected_figure_stems,
    audit_figure
  ),
  use.names = TRUE,
  fill = TRUE
)

if (
  nrow(figure_audit) !=
    length(
      expected_figure_stems
    ) ||
  any(
    figure_audit$
      files_valid != TRUE
  )
) {
  stop(
    "At least one Stage 7 figure failed the direct disk audit."
  )
}

write_csv_safe(
  figure_audit,
  file.path(
    DIRS$tables,
    "16_stage7_figure_export_audit.csv"
  )
)

warnings_table <- if (
  length(
    warning_records
  ) > 0L
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
    "17_stage7_warnings_and_nonfatal_issues.csv"
  )
)

stage7_correction_audit <- data.table::data.table(
  item = c(
    "Failed_v1_step",
    "Root_cause",
    "Correction",
    "Additional_runtime_guards",
    "Scientific_definition_changed"
  ),
  value = c(
    "Nfkb1/Rela communication-feature construction",
    "sample_accession appeared both as a grouped key and as an explicit output column in ligand_sample",
    "Removed sample_accession from j and retained it only through by = .(sample_accession)",
    "Unique column-name checks, one-row-per-sample checks, and locked sample-order checks after every grouped summary and merge",
    "FALSE"
  )
)

write_csv_safe(
  stage7_correction_audit,
  file.path(
    DIRS$methods,
    "stage7_FINAL_v2_correction_audit.csv"
  )
)

limitations <- data.table::data.table(
  limitation = c(
    "Sample size",
    "Internal evidence reuse",
    "No clinical-model claim",
    "Attribution interpretation",
    "Permutation resolution",
    "Communication score",
    "Feature panel"
  ),
  statement = c(
    "Only six biological samples were available. Model results are exploratory and cannot establish generalizable predictive performance.",
    "Stage 4-6 features and Stage 7 outcomes derive from the same six samples. Sample-held-out fitting prevents direct row leakage but does not create an independent validation cohort.",
    "The analysis is not a diagnostic model and should not be reported as a clinically validated classifier.",
    "beta multiplied by the held-out standardized feature value is an exact additive contribution to the linear predictor, not a causal effect and not a nonlinear SHAP estimate.",
    "There are only 20 possible balanced 3-versus-3 label assignments, so empirical P-value resolution is coarse.",
    "The communication-burden features summarize Stage 6 ligand-receptor axes and do not prove physical signaling.",
    "The primary and extended feature panels were prespecified; no performance-based feature selection was used."
  )
)

write_csv_safe(
  limitations,
  file.path(
    DIRS$methods,
    "stage7_limitations_and_claim_boundaries.csv"
  )
)

parameter_table <- data.table::data.table(
  parameter = c(
    "Random seed",
    "Biological modeling unit",
    "Validation",
    "Primary panel",
    "Extended panel",
    "Primary ridge lambda",
    "Lambda sensitivity",
    "Model optimization",
    "Balanced label assignments",
    "NFKB communication axes",
    "Bhlhe40 communication axes",
    "Strict-support weight",
    "Feature selection",
    "Clinical classifier claim",
    "Attribution"
  ),
  value = c(
    "20260714",
    "Biological sample",
    "Exhaustive leave-one-Control-plus-one-HFpEF-pair-out; 9 pairs",
    paste(
      PRIMARY_PANEL,
      collapse = ";"
    ),
    paste(
      EXTENDED_PANEL,
      collapse = ";"
    ),
    as.character(
      PRIMARY_RIDGE_LAMBDA
    ),
    paste(
      RIDGE_LAMBDA_SENSITIVITY,
      collapse = ";"
    ),
    "Base-R BFGS ridge logistic regression with analytic gradient",
    as.character(
      nrow(
        permutation_null
      )
    ),
    as.character(
      nrow(
        nfkb_communication$
          axes
      )
    ),
    as.character(
      nrow(
        bhlhe40_communication$
          axes
      )
    ),
    as.character(
      STRICT_SUPPORT_WEIGHT
    ),
    "None; panels prespecified before modeling",
    "FALSE",
    "Exact additive linear-predictor contribution: beta_j multiplied by training-standardized held-out feature_j"
  ),
  rationale = c(
    "Reproducibility",
    "Avoid cell-level pseudoreplication",
    "Every validation fold contains one held-out sample from each class",
    "Compact cross-stage biological hypothesis",
    "Prespecified sensitivity analysis",
    "Avoid unstable tuning with six samples",
    "Penalty sensitivity only",
    "Avoid external ML-package API dependence",
    "Exact enumeration of all balanced sample labels",
    "Summarize prioritized NFKB1-RELA communication",
    "Retain BHLHE40 communication contrast",
    "Upweight strict cross-stage support without excluding other axes",
    "Avoid outcome-driven selection",
    "Six samples are insufficient",
    "Linear-logit decomposition with an explicit interpretation boundary"
  )
)

write_csv_safe(
  parameter_table,
  file.path(
    DIRS$methods,
    "stage7_parameters_and_rationale.csv"
  )
)

methods_text <- c(
  "HFpEF Stage 7 FINAL v2",
  "Cross-stage constrained sample-level ridge classification and exact additive linear-predictor attribution",
  "",
  "Input boundary:",
  "- Stage 4 supplied sample-level macrophage TF activity and macrophage pseudobulk log2-CPM.",
  "- Stage 5B fixed Bhlhe40, Nfkb1, and Rela as the Stage 7 TF candidates.",
  "- Stage 6 supplied prioritized macrophage-to-vascular/stromal ligand-receptor axes and selected ligand/receptor expression by biological sample.",
  "",
  "Feature construction:",
  "- FINAL v2 removed a duplicated sample_accession grouping column from the Stage 6 ligand summary that caused merge.data.table to stop in FINAL v1.",
  "- Every grouped ligand/receptor summary and merge product was required to have unique column names and one row per biological-sample key.",
  "- The locked six-sample order was restored and verified after every communication-feature merge.",
  "- The primary panel contained three candidate TF activities, the Full Stage 2 Top150 drug-opposed macrophage program score, and an Nfkb1/Rela communication-burden score.",
  "- The extended sensitivity panel additionally contained the Stage3-supported Top150 program score and Bhlhe40 communication burden.",
  "- Communication burden was calculated from oriented macrophage ligand and receiver receptor-component log2-CPM, weighted by inverse square-root axis rank and a prespecified strict-support multiplier.",
  "",
  "Validation and model:",
  "- The model used six biological samples only.",
  "- Each of the nine folds held out one Control and one HFpEF sample.",
  "- Feature imputation and standardization were estimated from the four training samples in each fold.",
  "- A fixed-penalty ridge logistic regression was fitted with base-R BFGS optimization and analytic gradients.",
  "- No outcome-driven feature selection or hyperparameter tuning was performed.",
  "",
  "Permutation null:",
  "- All 20 balanced 3-versus-3 sample-label assignments were evaluated with the same primary model and leave-pair-out workflow.",
  "- Empirical P values compare the observed performance with the 19 nonobserved balanced label assignments.",
  "",
  "Attribution:",
  "- For each held-out prediction, feature attribution was beta_j multiplied by the training-standardized held-out feature value.",
  "- Contributions sum with the intercept to the fitted linear predictor.",
  "- These are exact additive logit contributions for the linear model, not causal effects and not a general nonlinear SHAP analysis.",
  "",
  "Claim boundary:",
  "- Stage 7 is an exploratory internal separability analysis.",
  "- It is not external validation and not a clinical prediction model."
)

writeLines(
  methods_text,
  file.path(
    DIRS$methods,
    "stage7_methods_and_claim_boundaries.txt"
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

workbook_path <- file.path(
  DIRS$tables,
  "18_stage7_cross_stage_sample_ridge_key_results.xlsx"
)

workbook_sheets <- list(
  Upstream_status =
    as.data.frame(
      upstream_status_audit
    ),
  Candidate_TFs =
    as.data.frame(
      candidate_cross_stage_manifest
    ),
  Feature_definitions =
    as.data.frame(
      feature_definitions
    ),
  Sample_features =
    as.data.frame(
      sample_features
    ),
  Fold_predictions =
    as.data.frame(
      primary_result$
        predictions
    ),
  Fold_performance =
    as.data.frame(
      primary_result$
        folds
    ),
  Sample_predictions =
    as.data.frame(
      primary_result$
        sample_predictions
    ),
  Attribution =
    as.data.frame(
      feature_importance
    ),
  Permutation_summary =
    as.data.frame(
      permutation_summary
    ),
  Sensitivity =
    as.data.frame(
      sensitivity_summary
    ),
  Figure_audit =
    as.data.frame(
      figure_audit
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

xlsx_structure_ok <- all(
  c(
    "[Content_Types].xml",
    "xl/workbook.xml",
    "xl/worksheets/sheet1.xml"
  ) %in%
    xlsx_contents$Name
)

if (!xlsx_structure_ok) {
  stop(
    "Stage 7 workbook failed internal structure validation."
  )
}

############################################################
## 14. Scientific completion checks and final status
############################################################

scientific_checks <- data.table::data.table(
  check = c(
    "Stage 3 completed",
    "Stage 4 completed",
    "Stage 5B completed",
    "Stage 6 completed",
    "Upstream failed checks",
    "Biological samples",
    "Control samples",
    "HFpEF samples",
    "Candidate TF activities",
    "Primary features",
    "Extended features",
    "Primary leave-pair-out folds",
    "Primary prediction rows",
    "Primary sample prediction rows",
    "Primary converged folds",
    "Balanced label assignments",
    "Permutation observed labeling",
    "Feature attribution rows",
    "Sensitivity scenarios",
    "Required figures",
    "Figure failures",
    "Scientific checkpoint",
    "Workbook",
    "Workbook structure"
  ),
  observed = c(
    as.integer(
      stage3_status$
        overall_status[1L] ==
        expected_status[
          "Stage3"
        ]
    ),
    as.integer(
      stage4_status$
        overall_status[1L] ==
        expected_status[
          "Stage4"
        ]
    ),
    as.integer(
      stage5b_status$
        overall_status[1L] ==
        expected_status[
          "Stage5B"
        ]
    ),
    as.integer(
      stage6_status$
        overall_status[1L] ==
        expected_status[
          "Stage6"
        ]
    ),
    sum(
      stage3_checks$status !=
        "PASS"
    ) +
      sum(
        stage4_checks$status !=
          "PASS"
      ) +
      sum(
        stage5b_checks$status !=
          "PASS"
      ) +
      sum(
        stage6_checks$status !=
          "PASS"
      ),
    nrow(
      sample_features
    ),
    sum(
      sample_features$y ==
        0L
    ),
    sum(
      sample_features$y ==
        1L
    ),
    sum(
      candidate_resolution$
        activity_available ==
        TRUE
    ),
    length(
      intersect(
        PRIMARY_PANEL,
        names(
          sample_features
        )
      )
    ),
    length(
      intersect(
        EXTENDED_PANEL,
        names(
          sample_features
        )
      )
    ),
    nrow(
      primary_result$
        folds
    ),
    nrow(
      primary_result$
        predictions
    ),
    nrow(
      primary_result$
        sample_predictions
    ),
    primary_result$
      folds_converged,
    nrow(
      permutation_null
    ),
    sum(
      permutation_null$
        is_observed_labeling ==
        TRUE
    ),
    nrow(
      feature_importance
    ),
    nrow(
      sensitivity_summary
    ),
    nrow(
      figure_audit
    ),
    sum(
      figure_audit$
        files_valid != TRUE
    ),
    as.integer(
      file.exists(
        file.path(
          DIRS$objects,
          "CHECKPOINT_stage7_scientific_results_pre_figures.rds"
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
    1,
    1,
    1,
    1,
    0,
    6,
    3,
    3,
    3,
    length(
      PRIMARY_PANEL
    ),
    length(
      EXTENDED_PANEL
    ),
    9,
    18,
    6,
    9,
    20,
    1,
    length(
      PRIMARY_PANEL
    ),
    4,
    length(
      expected_figure_stems
    ),
    0,
    1,
    1,
    1
  ),
  comparison = c(
    rep(
      "equal",
      14
    ),
    "at_least",
    rep(
      "equal",
      9
    )
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
    "19_stage7_scientific_completion_checks.csv"
  )
)

if (
  any(
    scientific_checks$status !=
      "PASS"
  )
) {
  failed_checks <- scientific_checks[
    status != "PASS",
    check
  ]

  stop(
    "Stage 7 scientific completion check(s) failed: ",
    paste(
      failed_checks,
      collapse = "; "
    )
  )
}

END_TIME <- Sys.time()

overall_status <-
  "COMPLETED_STAGE7_READY_FOR_REVIEW"

run_status <- data.table::data.table(
  stage =
    STAGE_NAME,
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
    3
  ),
  biological_samples =
    nrow(
      sample_features
    ),
  primary_features =
    length(
      PRIMARY_PANEL
    ),
  extended_features =
    length(
      EXTENDED_PANEL
    ),
  primary_ridge_lambda =
    PRIMARY_RIDGE_LAMBDA,
  leave_pair_out_folds =
    primary_result$
      total_folds,
  converged_folds =
    primary_result$
      folds_converged,
  pairwise_auc =
    primary_result$
      pairwise_auc,
  sample_auc =
    primary_result$
      sample_auc,
  pairwise_empirical_p =
    pairwise_empirical_p,
  sample_empirical_p =
    sample_empirical_p,
  top_attributed_feature =
    feature_importance$
      feature[1L],
  warnings =
    nrow(
      warnings_table
    ),
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
    "20_stage7_run_status.csv"
  )
)

readme <- c(
  "HFpEF Reanalysis Project - Stage 7 FINAL v2",
  "Cross-stage constrained sample-level ridge classification and exact additive linear-predictor attribution",
  "",
  paste0(
    "Status: ",
    overall_status
  ),
  "",
  paste0(
    "Primary pairwise AUC: ",
    round(
      primary_result$
        pairwise_auc,
      4
    )
  ),
  paste0(
    "Primary sample AUC: ",
    round(
      primary_result$
        sample_auc,
      4
    )
  ),
  paste0(
    "Balanced-label pairwise empirical P: ",
    signif(
      pairwise_empirical_p,
      4
    )
  ),
  paste0(
    "Top attributed feature: ",
    feature_importance$
      feature[1L]
  ),
  "",
  "This is an exploratory internal six-sample analysis, not a clinical classifier or external validation.",
  "",
  "Upload the Stage 7 CHECK package before Stage 8."
)

writeLines(
  readme,
  file.path(
    OUT_DIR,
    "README_stage7.txt"
  ),
  useBytes = TRUE
)

############################################################
## 15. CHECK package
############################################################

script_copy_status <-
  "NOT_DETECTED"

if (
  length(SCRIPT_FILE) == 1L &&
  !is.na(SCRIPT_FILE) &&
  file.exists(SCRIPT_FILE)
) {
  methods_copy <- file.copy(
    SCRIPT_FILE,
    file.path(
      DIRS$methods,
      basename(
        EXPECTED_SCRIPT_FILE
      )
    ),
    overwrite = TRUE
  )

  check_copy <- file.copy(
    SCRIPT_FILE,
    file.path(
      DIRS$check,
      basename(
        EXPECTED_SCRIPT_FILE
      )
    ),
    overwrite = TRUE
  )

  script_copy_status <- if (
    isTRUE(
      methods_copy
    ) &&
    isTRUE(
      check_copy
    )
  ) {
    "COPIED"
  } else {
    "COPY_FAILED"
  }
}

review_files <- c(
  LOG_FILE,
  file.path(
    DIRS$tables,
    c(
      "00_stage7_upstream_status_audit.csv",
      "01_stage7_candidate_TF_resolution.csv",
      "01A_stage7_candidate_cross_stage_manifest.csv",
      "02_stage7_feature_definitions.csv",
      "03_stage7_sample_level_feature_matrix.csv",
      "07_stage7_primary_LOPO_fold_performance.csv",
      "08_stage7_primary_LOPO_sample_predictions.csv",
      "11_stage7_feature_attribution_and_stability.csv",
      "13_stage7_permutation_summary.csv",
      "14_stage7_panel_and_lambda_sensitivity.csv",
      "15_stage7_model_performance_summary.csv",
      "16_stage7_figure_export_audit.csv",
      "17_stage7_warnings_and_nonfatal_issues.csv",
      "18_stage7_cross_stage_sample_ridge_key_results.xlsx",
      "19_stage7_scientific_completion_checks.csv",
      "20_stage7_run_status.csv"
    )
  ),
  file.path(
    DIRS$methods,
    c(
      "stage7_parameters_and_rationale.csv",
      "stage7_limitations_and_claim_boundaries.csv",
      "stage7_FINAL_v2_correction_audit.csv",
      "stage7_methods_and_claim_boundaries.txt",
      "sessionInfo.txt"
    )
  ),
  file.path(
    OUT_DIR,
    "README_stage7.txt"
  ),
  list.files(
    DIRS$figures,
    pattern = "\\.png$",
    full.names = TRUE
  )
)

review_files <- unique(
  review_files[
    file.exists(
      review_files
    )
  ]
)

for (source_file in review_files) {
  target_file <- file.path(
    DIRS$check,
    basename(
      source_file
    )
  )

  if (
    normalizePath(
      source_file,
      winslash = "/",
      mustWork = FALSE
    ) !=
      normalizePath(
        target_file,
        winslash = "/",
        mustWork = FALSE
      )
  ) {
    copied <- file.copy(
      source_file,
      target_file,
      overwrite = TRUE
    )

    if (!copied) {
      stop(
        "Failed to copy Stage 7 CHECK file: ",
        source_file
      )
    }
  }
}

check_files <- list.files(
  DIRS$check,
  full.names = TRUE,
  all.files = FALSE
)

check_manifest <- data.table::data.table(
  filename =
    basename(
      check_files
    ),
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
    function(file_path) {
      digest::digest(
        file = file_path,
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
    full.names = TRUE,
    all.files = FALSE
  ),
  root = DIRS$check
)

if (
  !file.exists(CHECK_ZIP) ||
  !is.finite(
    as.numeric(
      file.info(
        CHECK_ZIP
      )$size
    )
  ) ||
  as.numeric(
    file.info(
      CHECK_ZIP
    )$size
  ) <= 0
) {
  stop(
    "Stage 7 CHECK package was not created correctly."
  )
}

log_msg(
  "Stage 7 analysis completed."
)

log_msg(
  "Status: ",
  overall_status
)

log_msg(
  "Pairwise AUC: ",
  round(
    primary_result$
      pairwise_auc,
    4
  ),
  " | empirical P=",
  signif(
    pairwise_empirical_p,
    4
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
  "HFpEF Stage 7 completed\n"
)

cat(
  "Status: ",
  overall_status,
  "\n",
  sep = ""
)

cat(
  "Pairwise leave-pair-out AUC: ",
  round(
    primary_result$
      pairwise_auc,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "Balanced-label empirical P: ",
  signif(
    pairwise_empirical_p,
    4
  ),
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
  "============================================================\n"
)
