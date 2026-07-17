############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 6 FINAL v3 - DELIVERY REPAIR v1
##
## Purpose:
##   Repair only the failed delivery/audit/package section of the
##   completed Stage 6 FINAL v3 run.
##
## This script DOES NOT:
##   - reload the Seurat object;
##   - repeat expression summaries;
##   - repeat receiver pseudobulk analysis;
##   - repeat NicheNet ligand activity;
##   - repeat ligand-receptor-target ranking;
##   - regenerate figures.
##
## It DOES:
##   1) load the Stage 6 scientific checkpoint;
##   2) audit already exported PNG/PDF/TIFF figures directly from disk;
##   3) rebuild the key-results workbook;
##   4) rebuild scientific completion checks and run status;
##   5) create the final CHECK.zip.
##
## Required existing output:
##   <HFPEF_PROJECT_DIR>/
##   06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3
##
## Run:
##   source(
##     "<HFPEF_PROJECT_DIR>/HFpEF_Stage6_FINAL_v3_DELIVERY_REPAIR_v1.R",
##     encoding = "UTF-8"
##   )
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(warn = 1)
options(encoding = "UTF-8")

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

STAGE6_DIR <- file.path(
  PROJECT_DIR,
  "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3"
)

DIRS <- list(
  logs = file.path(
    STAGE6_DIR,
    "00_logs"
  ),
  tables = file.path(
    STAGE6_DIR,
    "01_tables"
  ),
  objects = file.path(
    STAGE6_DIR,
    "02_objects"
  ),
  figures = file.path(
    STAGE6_DIR,
    "03_figures"
  ),
  source = file.path(
    STAGE6_DIR,
    "04_source_data"
  ),
  methods = file.path(
    STAGE6_DIR,
    "05_methods"
  ),
  check = file.path(
    STAGE6_DIR,
    "06_review_check"
  )
)

CHECKPOINT_FILE <- file.path(
  DIRS$objects,
  "CHECKPOINT_stage6_scientific_results_pre_figures.rds"
)

CHECK_ZIP <- file.path(
  PROJECT_DIR,
  "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3_CHECK.zip"
)

EXPECTED_REPAIR_SCRIPT <- file.path(PROJECT_DIR, "06b_stage6_delivery_repair_v1.R")

START_TIME <- Sys.time()

required_directories <- c(
  STAGE6_DIR,
  unlist(
    DIRS[
      c(
        "logs",
        "tables",
        "objects",
        "figures",
        "source",
        "methods"
      )
    ],
    use.names = FALSE
  )
)

missing_directories <- required_directories[
  !dir.exists(required_directories)
]

if (length(missing_directories) > 0L) {
  stop(
    "Required Stage 6 v3 output directory is missing:\n",
    paste(
      missing_directories,
      collapse = "\n"
    )
  )
}

if (!file.exists(CHECKPOINT_FILE)) {
  stop(
    "Stage 6 scientific checkpoint is missing:\n",
    CHECKPOINT_FILE
  )
}

dir.create(
  DIRS$check,
  recursive = TRUE,
  showWarnings = FALSE
)

existing_check_files <- list.files(
  DIRS$check,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)

if (length(existing_check_files) > 0L) {
  unlink(
    existing_check_files,
    recursive = TRUE,
    force = TRUE
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
      "Required package(s) unavailable: ",
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
    "writexl",
    "zip",
    "digest"
  )
)

log_file <- file.path(
  DIRS$logs,
  "stage6_FINAL_v3_delivery_repair.log"
)

log_msg <- function(...) {
  line <- sprintf(
    "[%s] %s",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
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
    file = log_file,
    append = TRUE
  )
}

write_csv_safe <- function(
  table_object,
  path
) {
  data.table::fwrite(
    data.table::as.data.table(
      table_object
    ),
    path
  )
}

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
  file.exists(EXPECTED_REPAIR_SCRIPT)
) {
  SCRIPT_FILE <- normalizePath(
    EXPECTED_REPAIR_SCRIPT,
    winslash = "/",
    mustWork = TRUE
  )
}

log_msg(
  "Loading Stage 6 scientific checkpoint."
)

checkpoint <- readRDS(
  CHECKPOINT_FILE
)

required_checkpoint_objects <- c(
  "resource_audit",
  "candidate_rank",
  "candidate_ligand_coverage",
  "candidate_ligands",
  "cell_count_audit",
  "receiver_gene_set_summary",
  "ligand_activity",
  "ligand_target_links",
  "candidate_axes",
  "axis_ranking_sensitivity",
  "axis_ranking_stability",
  "candidate_summary"
)

missing_checkpoint_objects <- setdiff(
  required_checkpoint_objects,
  names(checkpoint)
)

if (length(missing_checkpoint_objects) > 0L) {
  stop(
    "Scientific checkpoint is missing object(s): ",
    paste(
      missing_checkpoint_objects,
      collapse = ", "
    )
  )
}

for (object_name in required_checkpoint_objects) {
  assign(
    object_name,
    data.table::as.data.table(
      checkpoint[[object_name]]
    )
  )
}

############################################################
## 1. Audit the already exported figures directly from disk
############################################################

base_figure_stems <- c(
  "Fig6A_top_candidate_TF_ligand_receptor_axes",
  "Fig6B_candidate_ligand_NicheNet_activity_heatmap",
  "Fig6C_receiver_receptor_expression_support",
  "Fig6E_candidate_TF_communication_axis_coverage"
)

optional_target_stem <-
  "Fig6D_top_axis_ligand_target_regulatory_potential"

optional_target_paths <- file.path(
  DIRS$figures,
  paste0(
    optional_target_stem,
    c(
      ".png",
      ".pdf",
      ".tiff"
    )
  )
)

expected_figure_stems <- base_figure_stems

if (all(file.exists(optional_target_paths))) {
  expected_figure_stems <- append(
    expected_figure_stems,
    optional_target_stem,
    after = 3L
  )
}

make_figure_audit_row <- function(stem_value) {
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

  exists_vector <- file.exists(paths)

  size_vector <- rep(
    NA_real_,
    length(paths)
  )

  size_vector[exists_vector] <- as.numeric(
    file.info(
      paths[exists_vector]
    )$size
  )

  data.table::data.table(
    stem = as.character(
      stem_value
    ),
    plot_type = if (
      stem_value %in%
        c(
          "Fig6B_candidate_ligand_NicheNet_activity_heatmap",
          "Fig6D_top_axis_ligand_target_regulatory_potential"
        )
    ) {
      "heatmap"
    } else {
      "ggplot"
    },
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
          is.finite(size_vector) &
            size_vector > 0
        )
    ),
    expected_main_figure = TRUE,
    audit_source =
      "DIRECT_DISK_AUDIT_AFTER_V3_DELIVERY_FAILURE"
  )
}

figure_export_audit <- data.table::rbindlist(
  lapply(
    expected_figure_stems,
    make_figure_audit_row
  ),
  use.names = TRUE,
  fill = TRUE
)

if (
  nrow(figure_export_audit) !=
    length(expected_figure_stems) ||
  any(
    figure_export_audit$
      files_valid != TRUE
  )
) {
  failed_figures <- figure_export_audit[
    files_valid != TRUE,
    stem
  ]

  stop(
    "One or more already generated Stage 6 figures are incomplete: ",
    paste(
      failed_figures,
      collapse = "; "
    )
  )
}

write_csv_safe(
  figure_export_audit,
  file.path(
    DIRS$tables,
    "20A_stage6_figure_export_audit.csv"
  )
)

############################################################
## 2. Rebuild workbook
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
    "Rebuilt Stage 6 workbook failed internal structure validation."
  )
}

############################################################
## 3. Repair note, methods, and scientific checks
############################################################

repair_note <- data.table::data.table(
  item = c(
    "Failed_original_step",
    "Root_cause",
    "Scientific_calculations_repeated",
    "Figures_regenerated",
    "Repair_scope",
    "Checkpoint_used"
  ),
  value = c(
    "Figure export audit construction after all figures were exported",
    paste0(
      "record_figure_export() used local list subassignment ",
      "figure_export_records[[...]] <- ... inside a function. ",
      "The outer list remained empty; the missing stem column then ",
      "caused R to resolve stem as graphics::stem and %in%/match failed."
    ),
    "FALSE",
    "FALSE",
    "Audit existing figure files; rebuild workbook, checks, run status, and CHECK package",
    CHECKPOINT_FILE
  )
)

write_csv_safe(
  repair_note,
  file.path(
    DIRS$tables,
    "20B_stage6_delivery_repair_note.csv"
  )
)

methods_text <- c(
  "HFpEF Stage 6 FINAL v3 delivery repair",
  "",
  "The original Stage 6 FINAL v3 run completed and checkpointed all scientific calculations and exported the required figures.",
  "The run stopped during construction of the figure-export audit because record_figure_export() modified a local copy of figure_export_records rather than the outer list.",
  "This repair loaded the existing scientific checkpoint, audited the already exported PNG/PDF/TIFF files directly from disk, rebuilt the key-results workbook, rebuilt completion checks and run status, and generated the CHECK package.",
  "No Seurat processing, pseudobulk analysis, NicheNet activity calculation, ligand-receptor-target ranking, or figure generation was repeated.",
  "The biological interpretation remains ligand-receptor-target prioritization rather than causal validation."
)

writeLines(
  methods_text,
  file.path(
    DIRS$methods,
    "stage6_delivery_repair_methods_and_boundaries.txt"
  ),
  useBytes = TRUE
)

scientific_checks <- data.table::data.table(
  check = c(
    "Scientific checkpoint exists",
    "Required checkpoint objects",
    "Verified NicheNet resources",
    "Candidate TFs",
    "Candidates represented in ligand coverage",
    "Receiver gene-set summaries",
    "NicheNet ligand activities",
    "Positive ligand-target links",
    "Candidate communication axes",
    "Candidate communication summaries",
    "Ranking sensitivity scenarios",
    "Required figures audited",
    "Figure export failures",
    "Workbook exists",
    "Workbook structure"
  ),
  observed = c(
    as.integer(
      file.exists(
        CHECKPOINT_FILE
      )
    ),
    length(
      intersect(
        required_checkpoint_objects,
        names(checkpoint)
      )
    ),
    sum(
      resource_audit$valid ==
        TRUE
    ),
    nrow(
      candidate_rank
    ),
    nrow(
      candidate_ligand_coverage
    ),
    nrow(
      receiver_gene_set_summary
    ),
    nrow(
      ligand_activity
    ),
    nrow(
      ligand_target_links
    ),
    nrow(
      candidate_axes
    ),
    nrow(
      candidate_summary
    ),
    data.table::uniqueN(
      axis_ranking_sensitivity$
        scenario
    ),
    nrow(
      figure_export_audit
    ),
    sum(
      figure_export_audit$
        files_valid != TRUE
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
    length(
      required_checkpoint_objects
    ),
    2,
    3,
    3,
    8,
    1,
    1,
    1,
    3,
    1,
    length(
      expected_figure_stems
    ),
    0,
    1,
    1
  ),
  comparison = c(
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "equal",
    "at_least",
    "at_least",
    "at_least",
    "equal",
    "at_least",
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
    "22_stage6_scientific_completion_checks.csv"
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
    "Delivery repair scientific checks failed: ",
    paste(
      failed_checks,
      collapse = "; "
    )
  )
}

END_TIME <- Sys.time()

run_status <- data.table::data.table(
  stage =
    "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3",
  repair_stage =
    "DELIVERY_REPAIR_v1",
  repair_start_time = format(
    START_TIME,
    "%Y-%m-%d %H:%M:%S"
  ),
  repair_end_time = format(
    END_TIME,
    "%Y-%m-%d %H:%M:%S"
  ),
  repair_elapsed_minutes = round(
    as.numeric(
      difftime(
        END_TIME,
        START_TIME,
        units = "mins"
      )
    ),
    3
  ),
  scientific_calculations_repeated =
    FALSE,
  figures_regenerated =
    FALSE,
  candidate_TFs =
    paste(
      candidate_rank$tf_symbol,
      collapse = ";"
    ),
  candidate_ligands =
    nrow(candidate_ligands),
  NicheNet_activity_rows =
    nrow(ligand_activity),
  ligand_target_links =
    nrow(ligand_target_links),
  candidate_axes =
    nrow(candidate_axes),
  expected_main_figures =
    length(expected_figure_stems),
  validated_main_figures =
    sum(
      figure_export_audit$
        files_valid ==
        TRUE
    ),
  scientific_checks_failed =
    sum(
      scientific_checks$status !=
        "PASS"
    ),
  overall_status =
    "COMPLETED_STAGE6_READY_FOR_REVIEW"
)

write_csv_safe(
  run_status,
  file.path(
    DIRS$tables,
    "23_stage6_run_status.csv"
  )
)

readme <- c(
  "HFpEF Stage 6 FINAL v3",
  "Delivery repaired from the existing scientific checkpoint and exported figures.",
  "",
  "Status: COMPLETED_STAGE6_READY_FOR_REVIEW",
  "",
  "No scientific calculation or figure generation was repeated.",
  "The repair rebuilt the figure audit, workbook, completion checks, run status, and CHECK package.",
  "",
  paste0(
    "Validated figures: ",
    paste(
      expected_figure_stems,
      collapse = "; "
    )
  ),
  "",
  "Upload the FINAL_v3_CHECK.zip for review."
)

writeLines(
  readme,
  file.path(
    STAGE6_DIR,
    "README_stage6.txt"
  ),
  useBytes = TRUE
)

############################################################
## 4. Build final CHECK package
############################################################

script_copy_status <- "NOT_DETECTED"

if (
  length(SCRIPT_FILE) == 1L &&
  !is.na(SCRIPT_FILE) &&
  file.exists(SCRIPT_FILE)
) {
  file.copy(
    SCRIPT_FILE,
    file.path(
      DIRS$methods,
      basename(
        EXPECTED_REPAIR_SCRIPT
      )
    ),
    overwrite = TRUE
  )

  file.copy(
    SCRIPT_FILE,
    file.path(
      DIRS$check,
      basename(
        EXPECTED_REPAIR_SCRIPT
      )
    ),
    overwrite = TRUE
  )

  script_copy_status <- "COPIED"
}

key_table_files <- c(
  "00_stage6_nichenet_resource_audit.csv",
  "01_stage6_candidate_TF_manifest.csv",
  "03_stage6_candidate_ligand_coverage.csv",
  "09_stage6_receiver_gene_set_summary.csv",
  "12_stage6_NicheNet_ligand_activity.csv",
  "16_stage6_integrated_candidate_communication_axes.csv",
  "18_stage6_axis_ranking_stability_summary.csv",
  "19_stage6_candidate_TF_communication_summary.csv",
  "20A_stage6_figure_export_audit.csv",
  "20B_stage6_delivery_repair_note.csv",
  "20_stage6_TF_dependent_communication_key_results.xlsx",
  "22_stage6_scientific_completion_checks.csv",
  "23_stage6_run_status.csv"
)

review_files <- c(
  file.path(
    DIRS$tables,
    key_table_files
  ),
  file.path(
    DIRS$methods,
    "stage6_methods_and_claim_boundaries.txt"
  ),
  file.path(
    DIRS$methods,
    "stage6_delivery_repair_methods_and_boundaries.txt"
  ),
  file.path(
    DIRS$methods,
    "sessionInfo.txt"
  ),
  file.path(
    STAGE6_DIR,
    "README_stage6.txt"
  ),
  log_file,
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
        "Failed to copy CHECK file: ",
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
    "Final Stage 6 CHECK package was not created correctly."
  )
}

log_msg(
  "Stage 6 FINAL v3 delivery repair completed."
)

log_msg(
  "Status: COMPLETED_STAGE6_READY_FOR_REVIEW"
)

log_msg(
  "CHECK package: ",
  CHECK_ZIP
)

cat(
  "\n============================================================\n"
)

cat(
  "HFpEF Stage 6 FINAL v3 delivery repair completed\n"
)

cat(
  "Status: COMPLETED_STAGE6_READY_FOR_REVIEW\n"
)

cat(
  "Scientific calculations repeated: FALSE\n"
)

cat(
  "Figures regenerated: FALSE\n"
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
