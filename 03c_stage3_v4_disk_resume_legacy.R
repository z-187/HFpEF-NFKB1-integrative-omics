############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Major Revision
## Revision Benchmark and Ablation Module FINAL v3
##
## Purpose
##   This is NOT a new biological discovery stage and does not rerun
##   Stages 1-8. It reads the frozen Stage 4-8 outputs and addresses the
##   reviewer's requests for:
##     1) comparison with simpler baselines;
##     2) quantitative ranking metrics;
##     3) leave-one-layer-out ablation;
##     4) feature-importance reporting;
##     5) external-validation alignment.
##
## Clean-start behavior
##   - FINAL_v3 deletes only its own same-version output directory and CHECK zip.
##   - It never deletes or modifies Stage 1-8 outputs.
##   - It recomputes the benchmark/ablation module from frozen Stage 4-8 tables.
##
## Analysis boundaries
##   - External Stage 8 evidence is used only as an evaluation target.
##   - Stage 8 evidence is never used to construct discovery rankings.
##   - Candidate-TF comparisons are descriptive because only three TFs
##     were prospectively frozen for Stage 6-8 validation.
##   - Stage 7 AUC is an internal sample-separability metric, not a
##     clinical diagnostic-performance estimate.
##   - Communication comparisons evaluate ranking/generalization of the
##     frozen axes; they do not prove physical signaling or causality.
##
## Save as
##   <HFPEF_PROJECT_DIR>/
##   HFpEF_Revision_Benchmark_Ablation_FINAL_v3.R
##
## Run from a fresh R session. Do not paste line by line:
##   source(
##     "<HFPEF_ASCII_PROJECT_LINK>/HFpEF_Revision_Benchmark_Ablation_FINAL_v3.R",
##     encoding = "UTF-8",
##     echo = FALSE
##   )
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(warn = 1)
options(encoding = "UTF-8")
options(timeout = 7200)
set.seed(20260715)

############################################################
## 0. Project paths and run settings
############################################################

DIRECT_PROJECT_DIR <- Sys.getenv("HFPEF_PROJECT_DIR", unset = "")
ASCII_PROJECT_LINK <- Sys.getenv(
  "HFPEF_ASCII_PROJECT_LINK",
  unset = file.path(tempdir(), "HFPEF_STAGE8_ASCII_LINK")
)

project_dir_is_valid <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return(FALSE)
  dir.exists(path) &&
    dir.exists(file.path(path, "04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1")) &&
    dir.exists(file.path(path, "05_stage5_multiTF_virtual_perturbation_FIXED_v2")) &&
    dir.exists(file.path(path, "05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1")) &&
    dir.exists(file.path(path, "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3")) &&
    dir.exists(file.path(path, "07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2")) &&
    dir.exists(file.path(path, "08_stage8_multicohort_validation_FINAL_v6"))
}

detect_invoked_script <- function() {
  candidates <- character()
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- tryCatch(frames[[i]]$ofile, error = function(e) NULL)
    if (!is.null(ofile) && length(ofile) == 1L && nzchar(ofile)) {
      candidates <- c(candidates, ofile)
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  candidates <- c(
    candidates,
    sub("^--file=", "", grep("^--file=", args, value = TRUE))
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) == 0L) return(NA_character_)
  gsub("\\\\", "/", path.expand(candidates[1L]))
}

EARLY_SCRIPT_FILE <- detect_invoked_script()
EARLY_SCRIPT_DIR <- if (
  length(EARLY_SCRIPT_FILE) == 1L && !is.na(EARLY_SCRIPT_FILE) && nzchar(EARLY_SCRIPT_FILE)
) dirname(EARLY_SCRIPT_FILE) else NA_character_

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
      paste(paste0("- ", candidates), collapse = "\n")
    )
  }
  gsub("\\\\", "/", path.expand(candidates[which(valid)[1L]]))
})

STAGE4_DIR <- file.path(
  PROJECT_DIR, "04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1"
)
STAGE5_DIR <- file.path(
  PROJECT_DIR, "05_stage5_multiTF_virtual_perturbation_FIXED_v2"
)
STAGE5B_DIR <- file.path(
  PROJECT_DIR, "05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1"
)
STAGE6_DIR <- file.path(
  PROJECT_DIR, "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3"
)
STAGE7_DIR <- file.path(
  PROJECT_DIR, "07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2"
)
STAGE8_DIR <- file.path(
  PROJECT_DIR, "08_stage8_multicohort_validation_FINAL_v6"
)

OUT_NAME <- "REVISION_Benchmark_Ablation_FINAL_v3"
OUT_DIR <- file.path(PROJECT_DIR, OUT_NAME)
CHECK_ZIP <- file.path(PROJECT_DIR, paste0(OUT_NAME, "_CHECK.zip"))
EXPECTED_SCRIPT_FILE <- file.path(
  PROJECT_DIR,
  "R",
  "09_revision_benchmark_ablation_FINAL_v3.R"
)
ANALYSIS_SCHEMA <- "revision_benchmark_ablation_v3_clean_rebuild_20260715"
FORCE_REBUILD <- TRUE
FIGURE_DPI <- 600L
TOPK_VALUES <- c(5L, 10L, 20L)
SENTINEL_AXIS_KEY <- "NFKB1__TNF__TNFRSF1A__ENDOTHELIAL"
EXTERNAL_SUPPORT_THRESHOLD <- 0.50

############################################################
## 1. Packages, folders, and logging
############################################################

ensure_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    install.packages(
      missing,
      repos = "https://cloud.r-project.org",
      dependencies = TRUE
    )
  }
  missing_after <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_after) > 0L) {
    stop("Required package(s) unavailable: ", paste(missing_after, collapse = ", "))
  }
}

ensure_cran(c(
  "data.table", "ggplot2", "openxlsx", "scales", "digest", "zip"
))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

DIRS <- list(
  logs = file.path(OUT_DIR, "00_logs"),
  tables = file.path(OUT_DIR, "01_tables"),
  figures = file.path(OUT_DIR, "02_figures"),
  source = file.path(OUT_DIR, "03_source_data"),
  methods = file.path(OUT_DIR, "04_methods"),
  check = file.path(OUT_DIR, "05_review_check")
)

if (FORCE_REBUILD && dir.exists(OUT_DIR)) {
  unlink(OUT_DIR, recursive = TRUE, force = TRUE)
}
if (FORCE_REBUILD && file.exists(CHECK_ZIP)) {
  unlink(CHECK_ZIP, force = TRUE)
}
for (d in c(OUT_DIR, unlist(DIRS, use.names = FALSE))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(DIRS$logs, "revision_benchmark_ablation.log")
START_TIME <- Sys.time()

log_msg <- function(..., level = "INFO") {
  line <- sprintf(
    "[%s] [%s] %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    level,
    paste0(..., collapse = "")
  )
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
  invisible(line)
}

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, na = "", compress = "auto")
  invisible(path)
}

safe_fread <- function(path) {
  fread(path, encoding = "UTF-8", showProgress = FALSE)
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

log_msg("Revision benchmark/ablation started.")
log_msg("PROJECT_DIR: ", PROJECT_DIR)
log_msg("OUT_DIR: ", OUT_DIR)
log_msg("Analysis schema: ", ANALYSIS_SCHEMA)
log_msg("Clean rebuild enabled: ", FORCE_REBUILD)
log_msg("Script: ", ifelse(is.na(SCRIPT_FILE), "NOT_DETECTED", SCRIPT_FILE))

############################################################
## 2. Input files and upstream validation
############################################################

FILES <- list(
  s4_status = file.path(STAGE4_DIR, "01_tables", "22_stage4_run_status.csv"),
  s4_checks = file.path(STAGE4_DIR, "01_tables", "20_stage4_scientific_completion_checks.csv"),
  s4_weighted = file.path(STAGE4_DIR, "01_tables", "09_stage4_weighted_regulon_activity_HFpEF_vs_Control.csv"),
  s4_aucell = file.path(STAGE4_DIR, "01_tables", "10_stage4_AUCell_regulon_activity_HFpEF_vs_Control.csv"),
  s4_expression = file.path(STAGE4_DIR, "01_tables", "11_stage4_TF_expression_HFpEF_vs_Control.csv"),
  s4_priority = file.path(STAGE4_DIR, "01_tables", "12_stage4_candidate_TF_priority_score.csv"),
  s4_lopo = file.path(STAGE4_DIR, "01_tables", "15_stage4_leave_one_pair_out_TF_robustness_summary.csv"),
  s4_method_comparison = file.path(STAGE4_DIR, "01_tables", "18_stage4_TF_method_comparison_summary.csv"),

  s5_status = file.path(STAGE5_DIR, "01_tables", "21_stage5_run_status.csv"),
  s5_checks = file.path(STAGE5_DIR, "01_tables", "20_stage5_scientific_completion_checks.csv"),
  s5_rank = file.path(STAGE5_DIR, "01_tables", "13_stage5_candidate_TF_rank_aggregation.csv"),
  s5_sensitivity = file.path(STAGE5_DIR, "01_tables", "14_stage5_candidate_ranking_sensitivity_scenarios.csv"),
  s5_mode = file.path(STAGE5_DIR, "01_tables", "16_stage5_normalization_vs_attenuation_results.csv"),

  s5b_status = file.path(STAGE5B_DIR, "01_tables", "17_stage5B_run_status.csv"),
  s5b_checks = file.path(STAGE5B_DIR, "01_tables", "16_stage5B_scientific_completion_checks.csv"),
  s5b_rank = file.path(STAGE5B_DIR, "01_tables", "13_stage5B_final_candidate_robustness_rank.csv"),

  s6_status = file.path(STAGE6_DIR, "01_tables", "23_stage6_run_status.csv"),
  s6_checks = file.path(STAGE6_DIR, "01_tables", "22_stage6_scientific_completion_checks.csv"),
  s6_axes = file.path(STAGE6_DIR, "01_tables", "18_stage6_axis_ranking_stability_summary.csv"),
  s6_candidate_summary = file.path(STAGE6_DIR, "01_tables", "19_stage6_candidate_TF_communication_summary.csv"),
  s6_workbook = file.path(STAGE6_DIR, "01_tables", "20_stage6_TF_dependent_communication_key_results.xlsx"),

  s7_status = file.path(STAGE7_DIR, "01_tables", "20_stage7_run_status.csv"),
  s7_checks = file.path(STAGE7_DIR, "01_tables", "19_stage7_scientific_completion_checks.csv"),
  s7_attribution = file.path(STAGE7_DIR, "01_tables", "11_stage7_feature_attribution_and_stability.csv"),
  s7_performance = file.path(STAGE7_DIR, "01_tables", "15_stage7_model_performance_summary.csv"),

  s8_status = file.path(STAGE8_DIR, "01_tables", "73_stage8_run_status.csv"),
  s8_checks = file.path(STAGE8_DIR, "01_tables", "72_stage8_scientific_completion_checks.csv"),
  s8_tf_evidence = file.path(STAGE8_DIR, "01_tables", "61_multicohort_TF_evidence.csv.gz"),
  s8_tf_summary = file.path(STAGE8_DIR, "01_tables", "64_TF_integrated_summary.csv"),
  s8_axis_evidence = file.path(STAGE8_DIR, "01_tables", "62_multicohort_axis_evidence.csv.gz"),
  s8_axis_summary = file.path(STAGE8_DIR, "01_tables", "65_axis_integrated_summary.csv"),
  s8_roles = file.path(STAGE8_DIR, "01_tables", "66_dataset_roles_and_claim_boundaries.csv")
)

missing_files <- unlist(FILES, use.names = TRUE)[
  !file.exists(unlist(FILES, use.names = FALSE))
]
if (length(missing_files) > 0L) {
  stop(
    "Missing required frozen input(s):\n",
    paste(paste(names(missing_files), missing_files, sep = " = "), collapse = "\n")
  )
}

validate_stage <- function(status_file, checks_file, expected_status, stage_label) {
  status <- safe_fread(status_file)
  checks <- safe_fread(checks_file)
  if (!"overall_status" %in% names(status) || nrow(status) == 0L) {
    stop(stage_label, " status file lacks overall_status.")
  }
  if (!status$overall_status[1L] %in% expected_status) {
    stop(stage_label, " is not ready: ", status$overall_status[1L])
  }
  if (!all(c("status") %in% names(checks)) || any(checks$status != "PASS")) {
    stop(stage_label, " contains a non-PASS scientific check.")
  }
  data.table(
    stage = stage_label,
    overall_status = status$overall_status[1L],
    scientific_checks = nrow(checks),
    failed_checks = sum(checks$status != "PASS")
  )
}

upstream_audit <- rbindlist(list(
  validate_stage(
    FILES$s4_status, FILES$s4_checks,
    c("COMPLETED_STAGE4_READY_FOR_REVIEW", "COMPLETED_STAGE4_READY_WITH_METHOD_CAUTION"),
    "Stage4"
  ),
  validate_stage(
    FILES$s5_status, FILES$s5_checks,
    "COMPLETED_STAGE5_READY_FOR_REVIEW",
    "Stage5"
  ),
  validate_stage(
    FILES$s5b_status, FILES$s5b_checks,
    c("COMPLETED_STAGE5B_OFFLINE_READY_FOR_REVIEW", "COMPLETED_STAGE5B_READY_FOR_REVIEW"),
    "Stage5B"
  ),
  validate_stage(
    FILES$s6_status, FILES$s6_checks,
    "COMPLETED_STAGE6_READY_FOR_REVIEW",
    "Stage6"
  ),
  validate_stage(
    FILES$s7_status, FILES$s7_checks,
    "COMPLETED_STAGE7_READY_FOR_REVIEW",
    "Stage7"
  ),
  validate_stage(
    FILES$s8_status, FILES$s8_checks,
    "COMPLETED_STAGE8_MULTICOHORT_READY_FOR_REVIEW",
    "Stage8"
  )
))
write_csv_safe(upstream_audit, file.path(DIRS$tables, "01_upstream_status_audit.csv"))

############################################################
## 3. General analysis utilities
############################################################

rank_desc <- function(x) {
  frank(-as.numeric(x), ties.method = "average", na.last = "keep")
}

rank_asc <- function(x) {
  frank(as.numeric(x), ties.method = "average", na.last = "keep")
}

rank_mean <- function(...) {
  m <- do.call(cbind, lapply(list(...), as.numeric))
  rowMeans(m, na.rm = TRUE)
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

canonical_gene_key <- function(x) {
  y <- trimws(as.character(x))
  for (cp in c(0x2010, 0x2011, 0x2012, 0x2013, 0x2014, 0x2212)) {
    y <- gsub(intToUtf8(cp), "-", y, fixed = TRUE)
  }
  y <- gsub("[[:space:]]+", "", y)
  toupper(y)
}

canonical_axis_key <- function(tf_symbol, ligand, receptor, receiver) {
  paste(
    canonical_gene_key(tf_symbol),
    canonical_gene_key(ligand),
    canonical_gene_key(receptor),
    toupper(gsub("[^A-Za-z0-9]", "", as.character(receiver))),
    sep = "__"
  )
}

save_plot_all <- function(plot_object, stem, width, height) {
  png_path <- file.path(DIRS$figures, paste0(stem, ".png"))
  pdf_path <- file.path(DIRS$figures, paste0(stem, ".pdf"))
  tif_path <- file.path(DIRS$figures, paste0(stem, ".tiff"))
  ggsave(png_path, plot_object, width = width, height = height, dpi = FIGURE_DPI, bg = "white")
  ggsave(pdf_path, plot_object, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(
    tif_path, plot_object, width = width, height = height,
    dpi = FIGURE_DPI, compression = "lzw", bg = "white"
  )
  data.table(
    stem = stem,
    png = png_path,
    pdf = pdf_path,
    tiff = tif_path,
    png_valid = file.exists(png_path) && file.info(png_path)$size > 1000,
    pdf_valid = file.exists(pdf_path) && file.info(pdf_path)$size > 1000,
    tiff_valid = file.exists(tif_path) && file.info(tif_path)$size > 1000
  )
}

theme_revision <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      strip.background = element_rect(fill = "grey95", colour = "black"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

############################################################
## 4. TF baseline comparison across the full Stage 4 universe
############################################################

s4_priority <- safe_fread(FILES$s4_priority)
s4_weighted <- safe_fread(FILES$s4_weighted)
s4_aucell <- safe_fread(FILES$s4_aucell)
s4_expression <- safe_fread(FILES$s4_expression)
s4_lopo <- safe_fread(FILES$s4_lopo)
s4_method_comparison <- safe_fread(FILES$s4_method_comparison)

required_s4_priority <- c(
  "tf_symbol", "priority_rank", "priority_score", "Nfkb1_forced"
)
if (!all(required_s4_priority %in% names(s4_priority))) {
  stop("Stage 4 priority table lacks required columns.")
}

all_tf <- merge(
  s4_priority,
  s4_expression[, .(
    tf_symbol,
    expression_effect = hfpef_minus_control,
    expression_hedges_g = hedges_g_HFpEF_vs_Control,
    expression_fdr = limma_padj
  )],
  by = "tf_symbol", all.x = TRUE, sort = FALSE,
  suffixes = c("", "_from_expression")
)
all_tf <- merge(
  all_tf,
  s4_weighted[, .(
    tf_symbol,
    regulon_effect = hfpef_minus_control,
    regulon_hedges_g = hedges_g_HFpEF_vs_Control,
    regulon_fdr = limma_padj
  )],
  by = "tf_symbol", all.x = TRUE, sort = FALSE
)
all_tf <- merge(
  all_tf,
  s4_aucell[, .(
    tf_symbol,
    aucell_effect = hfpef_minus_control,
    aucell_hedges_g = hedges_g_HFpEF_vs_Control,
    aucell_fdr = limma_padj
  )],
  by = "tf_symbol", all.x = TRUE, sort = FALSE,
  suffixes = c("", "_from_aucell")
)
all_tf <- merge(
  all_tf,
  s4_lopo,
  by = "tf_symbol", all.x = TRUE, sort = FALSE,
  suffixes = c("", "_from_lopo")
)

all_tf[, rank_TF_expression_only := rank_desc(abs(expression_hedges_g))]
all_tf[, rank_weighted_regulon_only := rank_desc(abs(regulon_hedges_g))]
all_tf[, rank_AUCell_only := rank_desc(abs(aucell_hedges_g))]
all_tf[, rank_regulon_plus_LOPO := rank_mean(
  rank_desc(abs(regulon_hedges_g)),
  rank_desc(sign_stability),
  rank_asc(median_abs_effect_rank)
)]
all_tf[, rank_stage4_multifeature := as.numeric(priority_rank)]
all_tf[, expression_to_integrated_rank_gain :=
         rank_TF_expression_only - rank_stage4_multifeature]
all_tf[, regulon_to_integrated_rank_gain :=
         rank_weighted_regulon_only - rank_stage4_multifeature]
all_tf[, selected_for_virtual_perturbation :=
         tf_symbol %in% c("Bhlhe40", "Runx1", "Spi1", "Rel", "Nfkb1", "Rela")]
all_tf[, frozen_for_stage6_8 := tf_symbol %in% c("Bhlhe40", "Nfkb1", "Rela")]
setorder(all_tf, rank_stage4_multifeature)
write_csv_safe(all_tf, file.path(DIRS$tables, "02_TF_baseline_ranks_all174.csv"))
write_csv_safe(
  s4_method_comparison,
  file.path(DIRS$tables, "02A_stage4_existing_method_correlations.csv")
)

selected_tf_baseline <- all_tf[selected_for_virtual_perturbation == TRUE, .(
  tf_symbol,
  rank_TF_expression_only,
  rank_weighted_regulon_only,
  rank_AUCell_only,
  rank_regulon_plus_LOPO,
  rank_stage4_multifeature,
  expression_to_integrated_rank_gain,
  regulon_to_integrated_rank_gain,
  expression_effect,
  expression_hedges_g,
  regulon_effect,
  regulon_hedges_g,
  priority_score,
  Nfkb1_forced
)]
write_csv_safe(
  selected_tf_baseline,
  file.path(DIRS$tables, "03_selected_TF_baseline_rank_comparison.csv")
)

############################################################
## 5. Candidate-TF method comparison and layer ablation
############################################################

s5_rank <- safe_fread(FILES$s5_rank)
s5_sensitivity <- safe_fread(FILES$s5_sensitivity)
s5_mode <- safe_fread(FILES$s5_mode)
s5b_rank <- safe_fread(FILES$s5b_rank)
s6_candidate <- safe_fread(FILES$s6_candidate_summary)
s8_tf_summary <- safe_fread(FILES$s8_tf_summary)
s8_tf_evidence <- safe_fread(FILES$s8_tf_evidence)

candidate_tfs <- s8_tf_summary[order(integrated_rank), tf_symbol]
if (length(candidate_tfs) != 3L || !setequal(candidate_tfs, c("Bhlhe40", "Nfkb1", "Rela"))) {
  stop(
    "Expected exactly the prospectively frozen Stage 6-8 candidates: ",
    "Bhlhe40, Nfkb1, and Rela."
  )
}

candidate_rank <- all_tf[tf_symbol %in% candidate_tfs, .(
  tf_symbol,
  expression_hedges_g,
  regulon_hedges_g,
  stage4_priority_rank_global = rank_stage4_multifeature,
  Nfkb1_forced
)]
candidate_rank <- merge(
  candidate_rank,
  s5_rank[tf_symbol %in% candidate_tfs],
  by = "tf_symbol", all.x = TRUE, sort = FALSE
)
candidate_rank <- merge(
  candidate_rank,
  s5b_rank[tf_symbol %in% candidate_tfs, .(
    tf_symbol,
    positive_recovery_probability,
    top1_frequency,
    top3_frequency,
    candidate_percentile,
    empirical_one_sided_p,
    final_robustness_rank,
    final_robustness_score,
    Nfkb1_forced_stage5B = Nfkb1_forced
  )],
  by = "tf_symbol", all.x = TRUE, sort = FALSE
)
candidate_rank <- merge(
  candidate_rank,
  s6_candidate[, .(
    tf_symbol,
    Stage5B_rank,
    total_axes,
    strict_cross_stage_axes,
    best_axis_rank,
    median_AUPR_corrected,
    median_NicheNet_pearson
  )],
  by = "tf_symbol", all.x = TRUE, sort = FALSE
)
candidate_rank <- merge(
  candidate_rank,
  s8_tf_summary[, .(
    tf_symbol,
    external_evidence_rows = evidence_rows,
    external_datasets = datasets,
    external_compartments = compartments,
    external_supported_rows = supported_rows,
    external_support_fraction = support_fraction,
    external_median_abs_hedges_g = median_abs_hedges_g,
    external_formal_fdr_rows = formal_fdr_rows,
    external_integrated_rank = integrated_rank
  )],
  by = "tf_symbol", all.x = TRUE, sort = FALSE
)

## Simple baselines and discovery-only integrated rankings.
candidate_rank[, rank_expression_only := rank_desc(abs(expression_hedges_g))]
candidate_rank[, rank_regulon_only := rank_desc(abs(regulon_hedges_g))]
candidate_rank[, rank_stage4_integrated := rank_asc(stage4_priority_rank_global)]

## Perturbation-only score excludes Stage 4 prior and all Stage 6/8 evidence.
candidate_rank[, rank_perturbation_only := rank_asc(rank_mean(
  rank_desc(stage2_primary_median_gap_reduction),
  rank_desc(stage2_primary_positive_fraction),
  rank_desc(biological_sample_improvement_fraction),
  rank_desc(inflammation_median_gap_reduction),
  rank_desc(specificity_score)
))]

## Bootstrap/null robustness is still discovery-only and was frozen before Stage 8.
candidate_rank[, rank_bootstrap_robustness := rank_asc(final_robustness_rank)]

## Communication-only score uses Stage 6 communication evidence and no Stage 8 result.
candidate_rank[, communication_mean_rank := rank_mean(
  rank_desc(strict_cross_stage_axes),
  rank_asc(best_axis_rank),
  rank_desc(total_axes),
  rank_desc(median_AUPR_corrected),
  rank_desc(median_NicheNet_pearson)
)]
candidate_rank[, rank_communication_only := rank_asc(communication_mean_rank)]

## Equal-weight cross-layer rank. External Stage 8 evidence is not included.
candidate_rank[, full_cross_layer_mean_rank := rank_mean(
  rank_stage4_integrated,
  rank_bootstrap_robustness,
  rank_communication_only
)]
candidate_rank[, rank_full_cross_layer := rank_asc(full_cross_layer_mean_rank)]
setorder(candidate_rank, rank_full_cross_layer)
write_csv_safe(
  candidate_rank,
  file.path(DIRS$tables, "04_candidate_TF_method_ranks.csv")
)

ablation_definitions <- list(
  Full_cross_layer = c(
    "rank_stage4_integrated",
    "rank_bootstrap_robustness",
    "rank_communication_only"
  ),
  Without_regulon_layer = c(
    "rank_bootstrap_robustness",
    "rank_communication_only"
  ),
  Without_perturbation_layer = c(
    "rank_stage4_integrated",
    "rank_communication_only"
  ),
  Without_communication_layer = c(
    "rank_stage4_integrated",
    "rank_bootstrap_robustness"
  ),
  TF_expression_only = "rank_expression_only",
  Regulon_only = "rank_regulon_only",
  Perturbation_only = "rank_perturbation_only",
  Communication_only = "rank_communication_only"
)

ablation_rows <- lapply(names(ablation_definitions), function(scenario_name) {
  cols <- ablation_definitions[[scenario_name]]
  dt <- candidate_rank[, c("tf_symbol", cols), with = FALSE]
  rank_matrix <- as.matrix(dt[, ..cols])
  dt[, scenario_score := rowMeans(rank_matrix, na.rm = TRUE)]
  dt[, scenario_rank := rank_asc(scenario_score)]
  dt[, `:=`(
    scenario = scenario_name,
    layers_used = paste(cols, collapse = ";")
  )]
  dt[, .(tf_symbol, scenario, layers_used, scenario_score, scenario_rank)]
})
ablation <- rbindlist(ablation_rows, use.names = TRUE, fill = TRUE)
write_csv_safe(ablation, file.path(DIRS$tables, "05_candidate_TF_ablation_scenarios.csv"))

ablation_stability <- ablation[, .(
  scenarios = uniqueN(scenario),
  median_rank = median(scenario_rank),
  best_rank = min(scenario_rank),
  worst_rank = max(scenario_rank),
  rank_range = max(scenario_rank) - min(scenario_rank),
  top1_frequency = mean(scenario_rank == 1),
  top2_frequency = mean(scenario_rank <= 2)
), by = tf_symbol][order(median_rank, -top1_frequency)]
write_csv_safe(
  ablation_stability,
  file.path(DIRS$tables, "06_candidate_TF_ablation_stability.csv")
)

candidate_method_long <- melt(
  candidate_rank,
  id.vars = c("tf_symbol", "external_integrated_rank"),
  measure.vars = c(
    "rank_expression_only",
    "rank_regulon_only",
    "rank_stage4_integrated",
    "rank_perturbation_only",
    "rank_bootstrap_robustness",
    "rank_communication_only",
    "rank_full_cross_layer"
  ),
  variable.name = "method",
  value.name = "discovery_rank"
)
candidate_method_long[, method := fcase(
  method == "rank_expression_only", "TF expression only",
  method == "rank_regulon_only", "Regulon activity only",
  method == "rank_stage4_integrated", "Stage 4 multifeature",
  method == "rank_perturbation_only", "Perturbation only",
  method == "rank_bootstrap_robustness", "Bootstrap perturbation",
  method == "rank_communication_only", "Communication only",
  method == "rank_full_cross_layer", "Full cross-layer",
  default = as.character(method)
)]

candidate_alignment <- candidate_method_long[, .(
  candidates = .N,
  spearman_vs_external_rank = safe_spearman(discovery_rank, external_integrated_rank),
  mean_absolute_rank_error = mean(abs(discovery_rank - external_integrated_rank)),
  top1_candidate = tf_symbol[which.min(discovery_rank)][1L],
  external_top1_candidate = tf_symbol[which.min(external_integrated_rank)][1L],
  top1_match = tf_symbol[which.min(discovery_rank)][1L] ==
    tf_symbol[which.min(external_integrated_rank)][1L]
), by = method]
write_csv_safe(
  candidate_method_long,
  file.path(DIRS$tables, "07_candidate_TF_method_long.csv")
)
write_csv_safe(
  candidate_alignment,
  file.path(DIRS$tables, "08_candidate_TF_external_alignment.csv")
)

## Preserve the original Stage 5 sensitivity and mode checks as quantitative context.
write_csv_safe(
  s5_sensitivity,
  file.path(DIRS$tables, "08A_stage5_original_ranking_sensitivity.csv")
)
write_csv_safe(
  s5_mode,
  file.path(DIRS$tables, "08B_stage5_normalization_attenuation_sensitivity.csv")
)

############################################################
## 6. Communication-axis method comparison
############################################################

log_msg("Reading Stage 6 Top_axes workbook for simple communication baselines.")
top_axes <- as.data.table(openxlsx::read.xlsx(FILES$s6_workbook, sheet = "Top_axes"))
required_axis_columns <- c(
  "tf_symbol", "nichenet_ligand", "receptor", "receiver", "axis_id",
  "final_axis_rank", "rank_sender_expression", "rank_receptor_expression",
  "rank_NicheNet_AUPR", "rank_NicheNet_pearson", "rank_cross_stage_support"
)
if (!all(required_axis_columns %in% names(top_axes))) {
  stop(
    "Stage 6 Top_axes sheet lacks required columns: ",
    paste(setdiff(required_axis_columns, names(top_axes)), collapse = ", ")
  )
}

for (col in setdiff(required_axis_columns, c(
  "tf_symbol", "nichenet_ligand", "receptor", "receiver", "axis_id"
))) {
  set(top_axes, j = col, value = as.numeric(top_axes[[col]]))
}
top_axes[, axis_key := canonical_axis_key(
  tf_symbol, nichenet_ligand, receptor, receiver
)]

## One record per biological axis. The Stage 6 workbook can contain repeated
## ligand-receptor records from target-level support rows; the minimum/median
## ranks below preserve the locked axis-level prioritization.
axis_discovery <- top_axes[, .(
  axis_id = axis_id[1L],
  tf_symbol = tf_symbol[1L],
  nichenet_ligand = nichenet_ligand[1L],
  receptor = receptor[1L],
  receiver = receiver[1L],
  full_integration_original_rank = min(final_axis_rank, na.rm = TRUE),
  ligand_expression_original_rank = min(rank_sender_expression, na.rm = TRUE),
  receptor_expression_original_rank = min(rank_receptor_expression, na.rm = TRUE),
  nichenet_original_rank = median(
    rowMeans(cbind(rank_NicheNet_AUPR, rank_NicheNet_pearson), na.rm = TRUE),
    na.rm = TRUE
  ),
  cross_stage_support_original_rank = min(rank_cross_stage_support, na.rm = TRUE)
), by = axis_key]

s8_axis_summary <- safe_fread(FILES$s8_axis_summary)
s8_axis_evidence <- safe_fread(FILES$s8_axis_evidence)
if (anyDuplicated(s8_axis_summary$axis_key)) {
  stop("Stage 8 axis summary contains duplicated axis_key values.")
}

axis_eval <- merge(
  axis_discovery,
  s8_axis_summary[, .(
    axis_key,
    external_axis_id = axis_id,
    external_evidence_rows = evidence_rows,
    external_datasets = datasets,
    external_dataset_list = dataset_list,
    external_supported_rows = supported_rows,
    external_support_fraction = support_fraction,
    external_median_abs_hedges_g = median_abs_hedges_g,
    external_formal_fdr_rows = formal_fdr_rows,
    external_integrated_rank = integrated_rank
  )],
  by = "axis_key", all = FALSE, sort = FALSE
)

if (nrow(axis_eval) < 10L) {
  stop("Fewer than 10 frozen Stage 6 axes were recovered in Stage 8 validation.")
}

## Re-rank all methods only within the same frozen externally evaluated axis set.
axis_eval[, rank_full_integration := rank_asc(full_integration_original_rank)]
axis_eval[, rank_ligand_expression_only := rank_asc(ligand_expression_original_rank)]
axis_eval[, rank_receptor_expression_only := rank_asc(receptor_expression_original_rank)]
axis_eval[, rank_LR_expression_only := rank_asc(rank_mean(
  rank_ligand_expression_only,
  rank_receptor_expression_only
))]
axis_eval[, rank_NicheNet_only := rank_asc(nichenet_original_rank)]
axis_eval[, rank_cross_stage_support_only := rank_asc(cross_stage_support_original_rank)]
setorder(axis_eval, external_integrated_rank)
write_csv_safe(axis_eval, file.path(DIRS$tables, "09_axis_method_ranks.csv"))

axis_method_long <- melt(
  axis_eval,
  id.vars = c(
    "axis_key", "tf_symbol", "nichenet_ligand", "receptor", "receiver",
    "external_integrated_rank", "external_support_fraction",
    "external_median_abs_hedges_g", "external_formal_fdr_rows"
  ),
  measure.vars = c(
    "rank_ligand_expression_only",
    "rank_receptor_expression_only",
    "rank_LR_expression_only",
    "rank_NicheNet_only",
    "rank_cross_stage_support_only",
    "rank_full_integration"
  ),
  variable.name = "method",
  value.name = "discovery_rank"
)
axis_method_long[, method := fcase(
  method == "rank_ligand_expression_only", "Ligand expression only",
  method == "rank_receptor_expression_only", "Receptor expression only",
  method == "rank_LR_expression_only", "Ligand-receptor expression",
  method == "rank_NicheNet_only", "NicheNet only",
  method == "rank_cross_stage_support_only", "Cross-stage support only",
  method == "rank_full_integration", "Full integration",
  default = as.character(method)
)]

axis_alignment <- axis_method_long[, .(
  axes = .N,
  spearman_vs_external_rank = safe_spearman(discovery_rank, external_integrated_rank),
  mean_absolute_rank_error = mean(abs(discovery_rank - external_integrated_rank)),
  median_external_support_fraction = median(external_support_fraction),
  externally_supported_axes = sum(
    external_support_fraction >= EXTERNAL_SUPPORT_THRESHOLD,
    na.rm = TRUE
  ),
  top1_axis = axis_key[which.min(discovery_rank)][1L],
  external_top1_axis = axis_key[which.min(external_integrated_rank)][1L],
  top1_match = axis_key[which.min(discovery_rank)][1L] ==
    axis_key[which.min(external_integrated_rank)][1L]
), by = method]
write_csv_safe(
  axis_method_long,
  file.path(DIRS$tables, "10_axis_method_long.csv.gz")
)
write_csv_safe(
  axis_alignment,
  file.path(DIRS$tables, "11_axis_external_alignment.csv")
)

topk_external <- rbindlist(lapply(unique(axis_method_long$method), function(method_i) {
  dt <- axis_method_long[method == method_i][order(discovery_rank)]
  rbindlist(lapply(TOPK_VALUES, function(k) {
    k_use <- min(k, nrow(dt))
    sub <- head(dt, k_use)
    data.table(
      method = method_i,
      top_k = k_use,
      mean_external_support_fraction = mean(sub$external_support_fraction, na.rm = TRUE),
      median_external_support_fraction = median(sub$external_support_fraction, na.rm = TRUE),
      externally_supported_axes = sum(
        sub$external_support_fraction >= EXTERNAL_SUPPORT_THRESHOLD,
        na.rm = TRUE
      ),
      supported_axis_fraction = mean(
        sub$external_support_fraction >= EXTERNAL_SUPPORT_THRESHOLD,
        na.rm = TRUE
      ),
      mean_external_abs_hedges_g = mean(sub$external_median_abs_hedges_g, na.rm = TRUE),
      formal_FDR_supported_axes = sum(sub$external_formal_fdr_rows > 0L, na.rm = TRUE)
    )
  }))
}))
write_csv_safe(
  topk_external,
  file.path(DIRS$tables, "12_axis_topk_external_performance.csv")
)

sentinel_ranks <- axis_method_long[axis_key == SENTINEL_AXIS_KEY, .(
  axis_key,
  method,
  discovery_rank,
  external_integrated_rank,
  external_support_fraction,
  external_median_abs_hedges_g,
  external_formal_fdr_rows
)]
if (nrow(sentinel_ranks) != uniqueN(axis_method_long$method)) {
  stop("The prespecified TNF-TNFRSF1A-Endothelial sentinel axis is incomplete.")
}
write_csv_safe(
  sentinel_ranks,
  file.path(DIRS$tables, "13_sentinel_axis_rank_comparison.csv")
)

############################################################
## 7. Stage 7 quantitative metrics and feature importance
############################################################

stage7_performance <- safe_fread(FILES$s7_performance)
stage7_attribution <- safe_fread(FILES$s7_attribution)
stage7_attribution[, contribution_fraction :=
                     mean_absolute_logit_contribution /
                     sum(mean_absolute_logit_contribution, na.rm = TRUE)]
stage7_attribution[, contribution_percent := 100 * contribution_fraction]
setorder(stage7_attribution, importance_rank)

write_csv_safe(
  stage7_performance,
  file.path(DIRS$tables, "14_stage7_quantitative_metrics.csv")
)
write_csv_safe(
  stage7_attribution,
  file.path(DIRS$tables, "15_stage7_feature_importance.csv")
)

############################################################
## 8. Revision figures
############################################################

figure_audit <- list()

method_order_tf <- c(
  "TF expression only",
  "Regulon activity only",
  "Stage 4 multifeature",
  "Perturbation only",
  "Bootstrap perturbation",
  "Communication only",
  "Full cross-layer",
  "External Stage 8"
)

tf_heat <- rbindlist(list(
  candidate_method_long[, .(
    tf_symbol,
    method,
    rank = discovery_rank
  )],
  candidate_rank[, .(
    tf_symbol,
    method = "External Stage 8",
    rank = external_integrated_rank
  )]
))
tf_heat[, method := factor(method, levels = method_order_tf)]
tf_heat[, tf_symbol := factor(tf_symbol, levels = c("Bhlhe40", "Nfkb1", "Rela"))]

p_tf_heat <- ggplot(tf_heat, aes(x = method, y = tf_symbol, fill = rank)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", rank)), size = 3) +
  scale_fill_gradient(low = "white", high = "grey35", trans = "reverse") +
  labs(
    title = "Candidate TF ranks across simple, integrated, and external analyses",
    x = NULL,
    y = NULL,
    fill = "Rank\n(lower is better)"
  ) +
  theme_revision(9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
figure_audit[[length(figure_audit) + 1L]] <- save_plot_all(
  p_tf_heat, "FigR1A_candidate_TF_method_rank_heatmap", 8.2, 3.3
)
write_csv_safe(tf_heat, file.path(DIRS$source, "FigR1A_source_data.csv"))

ablation_plot <- copy(ablation)
ablation_plot[, scenario := factor(
  scenario,
  levels = c(
    "Full_cross_layer",
    "Without_regulon_layer",
    "Without_perturbation_layer",
    "Without_communication_layer",
    "TF_expression_only",
    "Regulon_only",
    "Perturbation_only",
    "Communication_only"
  )
)]
p_ablation <- ggplot(
  ablation_plot,
  aes(x = scenario, y = scenario_rank, group = tf_symbol, linetype = tf_symbol)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.1) +
  scale_y_reverse(breaks = 1:3) +
  labs(
    title = "Leave-one-layer-out candidate-rank sensitivity",
    x = NULL,
    y = "Candidate rank",
    linetype = "TF"
  ) +
  theme_revision(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
figure_audit[[length(figure_audit) + 1L]] <- save_plot_all(
  p_ablation, "FigR1B_candidate_TF_ablation_stability", 8.2, 3.8
)
write_csv_safe(ablation_plot, file.path(DIRS$source, "FigR1B_source_data.csv"))

candidate_alignment_plot <- copy(candidate_alignment)
candidate_alignment_plot[, method := factor(method, levels = method_order_tf)]
p_tf_align <- ggplot(
  candidate_alignment_plot,
  aes(x = method, y = spearman_vs_external_rank)
) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Discovery-rank alignment with Stage 8 external TF evidence",
    x = NULL,
    y = "Spearman rank correlation"
  ) +
  theme_revision(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
figure_audit[[length(figure_audit) + 1L]] <- save_plot_all(
  p_tf_align, "FigR1C_candidate_TF_external_alignment", 7.2, 3.8
)
write_csv_safe(
  candidate_alignment_plot,
  file.path(DIRS$source, "FigR1C_source_data.csv")
)

method_order_axis <- c(
  "Ligand expression only",
  "Receptor expression only",
  "Ligand-receptor expression",
  "NicheNet only",
  "Cross-stage support only",
  "Full integration"
)
topk_plot <- copy(topk_external)
topk_plot[, method := factor(method, levels = method_order_axis)]
p_topk <- ggplot(
  topk_plot,
  aes(
    x = top_k,
    y = supported_axis_fraction,
    group = method,
    linetype = method
  )
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = TOPK_VALUES) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "External support recovered by simple and integrated axis rankings",
    x = "Top-ranked frozen axes",
    y = paste0("Axes with external support fraction >= ", sprintf("%.2f", EXTERNAL_SUPPORT_THRESHOLD)),
    linetype = "Ranking method"
  ) +
  theme_revision(9)
figure_audit[[length(figure_audit) + 1L]] <- save_plot_all(
  p_topk, "FigR1D_axis_topk_external_performance", 7.2, 4.2
)
write_csv_safe(topk_plot, file.path(DIRS$source, "FigR1D_source_data.csv"))

feature_plot <- copy(stage7_attribution)
feature_plot[, feature_label := fcase(
  feature == "COMM_NFkB_axis_burden", "NF-kB communication-axis burden",
  feature == "PROGRAM_DrugOpposed_Top150", "Drug-opposed Top150 program",
  feature == "TF_Nfkb1_activity", "NFKB1 regulon activity",
  feature == "TF_Bhlhe40_activity", "BHLHE40 regulon activity",
  feature == "TF_Rela_activity", "RELA regulon activity",
  default = as.character(feature)
)]
feature_plot[, feature_label := factor(
  feature_label,
  levels = rev(feature_label[order(importance_rank)])
)]
p_feature <- ggplot(
  feature_plot,
  aes(x = contribution_percent, y = feature_label)
) +
  geom_col(width = 0.7) +
  scale_x_continuous(labels = function(x) paste0(round(x), "%")) +
  labs(
    title = "Stage 7 sample-level feature contribution",
    x = "Share of mean absolute logit contribution",
    y = NULL
  ) +
  theme_revision(9)
figure_audit[[length(figure_audit) + 1L]] <- save_plot_all(
  p_feature, "FigR1E_stage7_feature_importance", 6.6, 3.8
)
write_csv_safe(feature_plot, file.path(DIRS$source, "FigR1E_source_data.csv"))

sentinel_plot <- copy(sentinel_ranks)
sentinel_plot[, method := factor(method, levels = method_order_axis)]
p_sentinel <- ggplot(
  sentinel_plot,
  aes(x = method, y = discovery_rank)
) +
  geom_col(width = 0.72) +
  scale_y_reverse() +
  labs(
    title = "Prespecified TNF-TNFRSF1A-Endothelial sentinel-axis rank",
    subtitle = paste0(
      "Stage 8 external integrated rank = ",
      unique(sentinel_plot$external_integrated_rank)
    ),
    x = NULL,
    y = "Discovery rank (lower is better)"
  ) +
  theme_revision(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
figure_audit[[length(figure_audit) + 1L]] <- save_plot_all(
  p_sentinel, "FigR1F_sentinel_axis_rank_comparison", 7.2, 4.0
)
write_csv_safe(sentinel_plot, file.path(DIRS$source, "FigR1F_source_data.csv"))

figure_audit_dt <- rbindlist(figure_audit, use.names = TRUE, fill = TRUE)
write_csv_safe(figure_audit_dt, file.path(DIRS$tables, "16_figure_export_audit.csv"))

############################################################
## 9. Interpretation boundaries and manuscript-ready summary
############################################################

claim_boundaries <- data.table(
  item = c(
    "TF expression baseline",
    "Regulon-only baseline",
    "Perturbation-only ranking",
    "Communication-only ranking",
    "Full cross-layer ranking",
    "External validation target",
    "Candidate-TF rank correlation",
    "Stage 7 AUC",
    "Sentinel TNF axis"
  ),
  permitted_interpretation = c(
    "Conventional TF differential-expression comparator.",
    "TF activity comparator based on target behavior rather than TF abundance.",
    "Program-recovery comparator excluding Stage 6 and Stage 8 evidence.",
    "Communication-network comparator excluding Stage 8 evidence.",
    "Equal-weight discovery integration of regulon, perturbation robustness, and communication layers.",
    "Independent or orthogonal Stage 8 evidence used only to evaluate frozen discovery rankings.",
    "Descriptive alignment only because the frozen candidate set contains three TFs.",
    "Internal sample-separability metric; not clinical sensitivity, specificity, or diagnostic accuracy.",
    "Prespecified representative inflammatory-endothelial branch; it need not rank first in every external cohort."
  ),
  prohibited_interpretation = c(
    "Do not treat expression rank as a gold-standard truth.",
    "Do not claim direct TF binding or causal regulation.",
    "Do not treat virtual perturbation as experimental knockout evidence.",
    "Do not claim physical ligand-receptor signaling or causal communication.",
    "Do not claim a universally optimal algorithm from one disease application.",
    "Do not feed external evidence back into discovery ranks.",
    "Do not report an inferential P value for a three-item rank correlation.",
    "Do not describe the model as a validated classifier or clinical predictor.",
    "Do not claim universal dominance over PDGFB, IL1B, or other network branches."
  )
)
write_csv_safe(
  claim_boundaries,
  file.path(DIRS$tables, "17_revision_claim_boundaries.csv")
)

summary_lines <- c(
  "Revision benchmark and ablation module",
  paste0("Analysis schema: ", ANALYSIS_SCHEMA),
  "",
  "Primary questions addressed:",
  "1. Are BHLHE40, NFKB1, and RELA recovered by simple TF-expression ranking?",
  "2. What is gained by regulon, perturbation, and communication integration?",
  "3. How stable are candidate ranks when one evidence layer is removed?",
  "4. Do frozen axis rankings recover external Stage 8 support better than simple expression baselines?",
  "5. Which Stage 7 feature layers contribute most to internal sample separation?",
  "",
  paste0(
    "Frozen TF candidates: ",
    paste(candidate_rank[order(rank_full_cross_layer), tf_symbol], collapse = "; ")
  ),
  paste0(
    "Stage 8 external TF order: ",
    paste(candidate_rank[order(external_integrated_rank), tf_symbol], collapse = "; ")
  ),
  paste0(
    "Sentinel axis: ", SENTINEL_AXIS_KEY,
    "; Stage 8 external rank = ",
    unique(sentinel_ranks$external_integrated_rank)
  ),
  "",
  "Interpretation:",
  "- The module tests added information and ranking stability; it does not require NFKB1 or the TNF axis to be universally rank 1.",
  "- BHLHE40 and NFKB1 may occupy complementary program-level and communication-centered roles.",
  "- Stage 8 is an external-evaluation layer, not an ingredient of the discovery score.",
  "- Rank correlations involving only three TFs are descriptive and should be reported without inferential overstatement."
)
writeLines(
  summary_lines,
  con = file.path(DIRS$methods, "revision_benchmark_ablation_summary.txt"),
  useBytes = TRUE
)

methods_lines <- c(
  "Revision benchmark and ablation methods",
  "",
  "TF baselines:",
  "TF-expression-only ranking was defined by the absolute Hedges' g for TF expression in HFpEF versus control macrophage pseudobulk samples.",
  "Regulon-only ranking was defined by the absolute Hedges' g for weighted regulon activity.",
  "The Stage 4 multifeature rank was imported without modification from the frozen Stage 4 priority table.",
  "Perturbation-only ranking was calculated from the mean rank of primary program gap reduction, positive-recovery fraction, biological-sample improvement fraction, inflammatory-program recovery, and perturbation specificity. Stage 4 and Stage 8 evidence were excluded.",
  "Communication-only ranking was calculated from the mean rank of strict cross-stage axes, best axis rank, total axis coverage, NicheNet AUPR, and NicheNet Pearson correlation. Stage 8 evidence was excluded.",
  "The full cross-layer discovery rank was the equal-weight mean of the Stage 4 multifeature rank, Stage 5B bootstrap-robustness rank, and Stage 6 communication-only rank.",
  "",
  "Ablation:",
  "Leave-one-layer-out scenarios omitted the regulon, perturbation, or communication component while retaining the other frozen layers. Single-layer baselines were also reported.",
  "",
  "External evaluation:",
  "Stage 8 integrated ranks, support fractions, effect sizes, and FDR-supported rows were used only after discovery ranks had been fixed. Candidate-TF comparisons were descriptive because three TFs were prospectively frozen for validation.",
  "",
  "Communication-axis benchmarks:",
  "The Stage 6 full integrated axis rank was compared with ligand-expression-only, receptor-expression-only, combined ligand-receptor-expression, NicheNet-only, and cross-stage-support-only rankings. All methods were re-ranked within the same frozen set of axes evaluated in Stage 8.",
  "Top-k external performance was summarized as the fraction of axes with Stage 8 support fraction >= 0.50, the mean external effect size, and the number of FDR-supported axes.",
  "",
  "Quantitative performance:",
  "Stage 7 leave-pair-out AUC, empirical permutation P values, and exact feature-attribution metrics were imported unchanged. These metrics quantify internal sample separability and are not clinical diagnostic-performance estimates."
)
writeLines(
  methods_lines,
  con = file.path(DIRS$methods, "revision_benchmark_ablation_methods.txt"),
  useBytes = TRUE
)

############################################################
## 10. Workbook, completion checks, script archive, and CHECK
############################################################

workbook_path <- file.path(
  DIRS$tables,
  "18_revision_benchmark_ablation_key_results.xlsx"
)
wb <- openxlsx::createWorkbook()
add_sheet <- function(name, x) {
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeDataTable(wb, name, as.data.frame(x))
  openxlsx::freezePane(wb, name, firstRow = TRUE)
  openxlsx::setColWidths(wb, name, cols = seq_len(ncol(x)), widths = "auto")
}
add_sheet("upstream_audit", upstream_audit)
add_sheet("TF_all174", all_tf)
add_sheet("TF_candidates", candidate_rank)
add_sheet("TF_ablation", ablation)
add_sheet("TF_ablation_stability", ablation_stability)
add_sheet("TF_external_alignment", candidate_alignment)
add_sheet("axis_ranks", axis_eval)
add_sheet("axis_external_alignment", axis_alignment)
add_sheet("axis_topk", topk_external)
add_sheet("sentinel_axis", sentinel_ranks)
add_sheet("Stage7_performance", stage7_performance)
add_sheet("Stage7_importance", stage7_attribution)
add_sheet("claim_boundaries", claim_boundaries)
openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)

script_archive_status <- "NOT_DETECTED"
script_archive_path <- file.path(DIRS$methods, basename(EXPECTED_SCRIPT_FILE))
if (!is.na(SCRIPT_FILE) && file.exists(SCRIPT_FILE)) {
  ok <- file.copy(SCRIPT_FILE, script_archive_path, overwrite = TRUE)
  if (isTRUE(ok) && file.exists(script_archive_path)) {
    script_archive_status <- "ARCHIVED"
  } else {
    script_archive_status <- "COPY_FAILED"
  }
}

## Script archiving is a reproducibility/packaging item, not a scientific
## validity criterion. Failure to detect a pasted interactive script is recorded
## transparently but must not invalidate completed benchmark calculations.
session_path <- file.path(DIRS$methods, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_path, useBytes = TRUE)
session_info_status <- if (
  file.exists(session_path) && is.finite(file.info(session_path)$size) &&
    file.info(session_path)$size > 0
) "PASS" else "WARNING"

reproducibility_checks <- data.table(
  item = c(
    "analysis_script_detected",
    "analysis_script_archived",
    "session_info_written",
    "analysis_schema_recorded"
  ),
  status = c(
    ifelse(!is.na(SCRIPT_FILE) && file.exists(SCRIPT_FILE), "PASS", "WARNING"),
    ifelse(script_archive_status == "ARCHIVED", "PASS", "WARNING"),
    session_info_status,
    "PASS"
  ),
  detail = c(
    ifelse(
      !is.na(SCRIPT_FILE) && file.exists(SCRIPT_FILE),
      SCRIPT_FILE,
      paste0(
        "The code was executed interactively or the source file was not ",
        "detectable. This does not alter computed results."
      )
    ),
    ifelse(
      script_archive_status == "ARCHIVED",
      script_archive_path,
      paste0(
        "No script copy was archived during this run. Archive the distributed ",
        "FINAL_v3 script in the public code repository."
      )
    ),
    session_path,
    ANALYSIS_SCHEMA
  )
)
write_csv_safe(
  reproducibility_checks,
  file.path(DIRS$tables, "19A_reproducibility_checks.csv")
)

completion_checks <- data.table(
  check = c(
    "upstream_stages_4_to_8_all_pass",
    "stage8_external_evidence_not_used_in_discovery_rank",
    "stage4_TF_universe_contains_at_least_150_TFs",
    "frozen_candidate_set_is_exactly_three_TFs",
    "Nfkb1_not_forced_in_stage4",
    "Nfkb1_not_forced_in_stage5B",
    "candidate_method_ranks_complete",
    "leave_one_layer_out_scenarios_complete",
    "candidate_external_alignment_complete",
    "stage6_top_axes_sheet_loaded",
    "at_least_10_frozen_axes_recovered_in_stage8",
    "no_duplicate_external_axis_keys",
    "axis_method_ranks_complete",
    "axis_topk_metrics_complete",
    "sentinel_TNF_axis_present_for_all_methods",
    "stage7_AUC_and_permutation_metrics_imported",
    "stage7_feature_importance_complete",
    "all_figure_exports_valid",
    "key_results_workbook_written"
  ),
  passed = c(
    all(upstream_audit$failed_checks == 0L),
    TRUE,
    nrow(all_tf) >= 150L,
    length(candidate_tfs) == 3L && setequal(candidate_tfs, c("Bhlhe40", "Nfkb1", "Rela")),
    all(all_tf[tf_symbol == "Nfkb1", Nfkb1_forced] == FALSE),
    all(candidate_rank[tf_symbol == "Nfkb1", Nfkb1_forced_stage5B] == FALSE),
    nrow(candidate_method_long) == 3L * 7L && all(is.finite(candidate_method_long$discovery_rank)),
    uniqueN(ablation$scenario) == length(ablation_definitions) &&
      nrow(ablation) == 3L * length(ablation_definitions),
    nrow(candidate_alignment) == 7L,
    nrow(top_axes) > 0L,
    nrow(axis_eval) >= 10L,
    !anyDuplicated(s8_axis_summary$axis_key),
    all(is.finite(axis_method_long$discovery_rank)),
    nrow(topk_external) == uniqueN(axis_method_long$method) * length(TOPK_VALUES),
    nrow(sentinel_ranks) == uniqueN(axis_method_long$method),
    all(c(
      "Primary pairwise leave-pair-out AUC",
      "Primary sample-level AUC",
      "Primary pairwise empirical permutation P",
      "Primary sample-level empirical permutation P"
    ) %in% stage7_performance$metric),
    nrow(stage7_attribution) >= 5L && abs(sum(stage7_attribution$contribution_fraction) - 1) < 1e-8,
    nrow(figure_audit_dt) >= 6L &&
      all(figure_audit_dt$png_valid) &&
      all(figure_audit_dt$pdf_valid) &&
      all(figure_audit_dt$tiff_valid),
    file.exists(workbook_path) && file.info(workbook_path)$size > 1000
  )
)
completion_checks[, status := fifelse(passed, "PASS", "FAIL")]
write_csv_safe(
  completion_checks,
  file.path(DIRS$tables, "19_scientific_completion_checks.csv")
)

failed_checks <- completion_checks[status == "FAIL", check]
if (length(failed_checks) > 0L) {
  stop(
    "Revision benchmark calculations finished but completion checks failed:\n",
    paste(failed_checks, collapse = "\n")
  )
}

end_time <- Sys.time()
run_status <- data.table(
  analysis = OUT_NAME,
  analysis_schema = ANALYSIS_SCHEMA,
  start_time = format(START_TIME, "%Y-%m-%d %H:%M:%S"),
  end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
  elapsed_minutes = as.numeric(difftime(end_time, START_TIME, units = "mins")),
  TF_universe = nrow(all_tf),
  frozen_TF_candidates = paste(candidate_tfs, collapse = ";"),
  ablation_scenarios = uniqueN(ablation$scenario),
  frozen_axes_evaluated = nrow(axis_eval),
  sentinel_axis = SENTINEL_AXIS_KEY,
  stage7_pairwise_AUC = stage7_performance[
    metric == "Primary pairwise leave-pair-out AUC", value
  ][1L],
  stage7_sample_AUC = stage7_performance[
    metric == "Primary sample-level AUC", value
  ][1L],
  scientific_checks_failed = 0L,
  script_copy_status = script_archive_status,
  reproducibility_warnings = sum(reproducibility_checks$status == "WARNING"),
  overall_status = "COMPLETED_REVISION_BENCHMARK_ABLATION_READY_FOR_REVIEW"
)
write_csv_safe(run_status, file.path(DIRS$tables, "20_run_status.csv"))

## Build a compact review package.
review_files <- c(
  file.path(DIRS$tables, c(
    "01_upstream_status_audit.csv",
    "02_TF_baseline_ranks_all174.csv",
    "03_selected_TF_baseline_rank_comparison.csv",
    "04_candidate_TF_method_ranks.csv",
    "05_candidate_TF_ablation_scenarios.csv",
    "06_candidate_TF_ablation_stability.csv",
    "08_candidate_TF_external_alignment.csv",
    "09_axis_method_ranks.csv",
    "11_axis_external_alignment.csv",
    "12_axis_topk_external_performance.csv",
    "13_sentinel_axis_rank_comparison.csv",
    "14_stage7_quantitative_metrics.csv",
    "15_stage7_feature_importance.csv",
    "16_figure_export_audit.csv",
    "17_revision_claim_boundaries.csv",
    "18_revision_benchmark_ablation_key_results.xlsx",
    "19_scientific_completion_checks.csv",
    "19A_reproducibility_checks.csv",
    "20_run_status.csv"
  )),
  list.files(DIRS$figures, full.names = TRUE),
  list.files(DIRS$source, full.names = TRUE),
  list.files(DIRS$methods, full.names = TRUE),
  LOG_FILE
)
review_files <- unique(review_files[file.exists(review_files)])
for (src in review_files) {
  target <- file.path(DIRS$check, basename(src))
  file.copy(src, target, overwrite = TRUE)
}

manifest <- data.table(
  filename = list.files(DIRS$check, full.names = FALSE),
  path = list.files(DIRS$check, full.names = TRUE)
)
manifest[, size_bytes := file.info(path)$size]
manifest[, md5 := unname(tools::md5sum(path))]
manifest[, path := NULL]
write_csv_safe(manifest, file.path(DIRS$check, "CHECK_package_file_manifest.csv"))

if (file.exists(CHECK_ZIP)) unlink(CHECK_ZIP, force = TRUE)
zip::zipr(
  zipfile = CHECK_ZIP,
  files = list.files(DIRS$check, full.names = TRUE),
  root = DIRS$check
)
if (!file.exists(CHECK_ZIP) || file.info(CHECK_ZIP)$size <= 1000) {
  stop("CHECK zip was not created correctly: ", CHECK_ZIP)
}

if (any(reproducibility_checks$status == "WARNING")) {
  log_msg(
    "Reproducibility warning(s) recorded: ",
    paste(
      reproducibility_checks[status == "WARNING", item],
      collapse = "; "
    ),
    ". Scientific completion remains PASS.",
    level = "WARN"
  )
}
log_msg("Revision benchmark and ablation completed successfully.")
log_msg("Status: ", run_status$overall_status)
log_msg("CHECK package: ", CHECK_ZIP)

cat("\n============================================================\n")
cat("Revision benchmark and ablation completed.\n")
cat("Status: ", run_status$overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat("============================================================\n")
