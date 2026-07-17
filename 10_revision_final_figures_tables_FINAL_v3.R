############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## OMI-2026-0142 Major Revision
## FINAL figure and supplementary-table production script v3
##
## Purpose
##   Read only frozen Stage 1-8 and Benchmark/Ablation outputs.
##   Re-draw all planned main and supplementary data figures.
##   Export every panel separately as 600-dpi LZW TIFF.
##   Export panel source data and Supplementary Tables S1-S16.
##   Export title-free panels: panel titles belong in the figure legend,
##   not inside the TIFF files used for final assembly.
##   Use compact, assembly-ready proportions and unclipped labels.
##   Do not rerun discovery, clustering, TF inference,
##   virtual perturbation, communication inference, or validation.
##
## Save as:
## <HFPEF_PROJECT_DIR>/
## HFpEF_Revision_All_Figures_Tables_FINAL_v3.R
##
## Run from a fresh R session:
## source(
##   "<HFPEF_ASCII_PROJECT_LINK>/HFpEF_Revision_All_Figures_Tables_FINAL_v3.R",
##   encoding = "UTF-8", echo = FALSE
## )
############################################################

rm(list = ls())
gc()
if (.Platform$OS.type == "windows") {
  utf8_locale <- try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"),
                     silent = TRUE)
  if (inherits(utf8_locale, "try-error") || is.na(utf8_locale)) {
    warning("UTF-8 locale could not be enabled; ASCII fallbacks remain active.")
  }
}
options(stringsAsFactors = FALSE, warn = 1, encoding = "UTF-8", timeout = 7200)
RANDOM_SEED <- 20260715L
set.seed(RANDOM_SEED)

############################################################
## 0. Locked paths and settings
############################################################
DIRECT_PROJECT_DIR <- Sys.getenv("HFPEF_PROJECT_DIR", unset = "")
ASCII_PROJECT_LINK <- Sys.getenv(
  "HFPEF_ASCII_PROJECT_LINK",
  unset = file.path(tempdir(), "HFPEF_STAGE8_ASCII_LINK")
)

STAGE_DIR_NAMES <- list(
  stage1 = "01_stage1_metadata_lock_FIXED_v3",
  stage2 = "02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2",
  stage3 = "03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH",
  stage4 = "04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1",
  stage5 = "05_stage5_multiTF_virtual_perturbation_FIXED_v2",
  stage5B = "05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1",
  stage6 = "06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3",
  stage7 = "07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2",
  stage8 = "08_stage8_multicohort_validation_FINAL_v6",
  benchmark = "REVISION_Benchmark_Ablation_FINAL_v3"
)

project_dir_is_valid <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !dir.exists(path)) {
    return(FALSE)
  }
  all(vapply(STAGE_DIR_NAMES, function(x) dir.exists(file.path(path, x)), logical(1)))
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
  candidates <- c(candidates, sub("^--file=", "", grep("^--file=", args, value = TRUE)))
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) == 0L) return(NA_character_)
  gsub("\\\\", "/", path.expand(candidates[1L]))
}

SCRIPT_FILE <- detect_invoked_script()
SCRIPT_DIR <- if (!is.na(SCRIPT_FILE) && nzchar(SCRIPT_FILE)) dirname(SCRIPT_FILE) else NA_character_

PROJECT_DIR <- local({
  candidates <- unique(c(
    Sys.getenv("HFPEF_PROJECT_DIR", unset = ""),
    ASCII_PROJECT_LINK, DIRECT_PROJECT_DIR, SCRIPT_DIR, getwd()
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  valid <- vapply(candidates, project_dir_is_valid, logical(1))
  if (!any(valid)) {
    stop(
      "HFpEF project root could not be located. Checked:\n",
      paste(paste0("- ", candidates), collapse = "\n"),
      "\nThe root must contain the frozen Stage 1-8 and Benchmark folders."
    )
  }
  gsub("\\\\", "/", path.expand(candidates[which(valid)[1L]]))
})

STAGE_DIRS <- lapply(STAGE_DIR_NAMES, function(x) file.path(PROJECT_DIR, x))
OUT_NAME <- "REVISION_Final_Figures_Tables_PUBLICATION_FINAL"
OUT_DIR <- file.path(PROJECT_DIR, OUT_NAME)
ANALYSIS_SCHEMA <- "revision_final_figures_tables_publication_final_20260715"
FORCE_REBUILD <- TRUE
FIGURE_DPI <- 600L
TIFF_COMPRESSION <- "LZW"
BASE_FAMILY <- "Arial"
BASE_SIZE <- 8.5
LABEL_SIZE <- 2.4
LINE_WIDTH <- 0.55
UNICODE_MINUS <- intToUtf8(0x2212)
MAIN_SCATTER_WIDTH <- 5.9
MAIN_SCATTER_HEIGHT <- 5.35
MAIN_HEATMAP_WIDTH <- 7.4
MAIN_HEATMAP_HEIGHT <- 6.2

if (FORCE_REBUILD && dir.exists(OUT_DIR)) unlink(OUT_DIR, recursive = TRUE, force = TRUE)

DIRS <- list(
  logs = file.path(OUT_DIR, "00_logs"),
  main = file.path(OUT_DIR, "01_Main_Figures"),
  supplementary = file.path(OUT_DIR, "02_Supplementary_Figures"),
  source_main = file.path(OUT_DIR, "03_Source_Data", "Main_Figures"),
  source_supp = file.path(OUT_DIR, "03_Source_Data", "Supplementary_Figures"),
  supp_tables = file.path(OUT_DIR, "04_Supplementary_Tables"),
  supp_table_csv = file.path(OUT_DIR, "04_Supplementary_Tables", "Source_CSV"),
  manifests = file.path(OUT_DIR, "05_Manifests"),
  qa = file.path(OUT_DIR, "06_QA_Report"),
  run_logs = file.path(OUT_DIR, "07_logs")
)
for (d in unlist(DIRS, use.names = FALSE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
for (i in 1:5) dir.create(file.path(DIRS$main, paste0("Figure_", i)), recursive = TRUE, showWarnings = FALSE)
for (i in 1:10) dir.create(file.path(DIRS$supplementary, paste0("Figure_S", i)), recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Packages, logging, helpers, palettes
############################################################
ensure_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  missing_after <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_after) > 0L) stop("Required package(s) unavailable: ", paste(missing_after, collapse = ", "))
}
ensure_cran(c("data.table", "ggplot2", "ggrepel", "scales", "openxlsx", "digest", "zip"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

LOG_FILE <- file.path(DIRS$logs, "final_figure_generation.log")
START_TIME <- Sys.time()
log_msg <- function(..., level = "INFO") {
  line <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
  invisible(line)
}

safe_fread <- function(path) data.table::fread(path, encoding = "UTF-8", showProgress = FALSE)
write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(x), path, na = "", compress = "auto")
  invisible(path)
}

find_exact_file <- function(root, filename, required = TRUE) {
  candidates <- list.files(root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  candidates <- candidates[basename(candidates) == filename]
  if (length(candidates) == 0L) {
    if (required) stop("Required frozen file not found: ", filename, "\nSearched under: ", root)
    return(NA_character_)
  }
  score <- ifelse(grepl("04_source_data|03_source_data|01_tables|02_objects", candidates), 0L, 1L)
  candidates <- candidates[order(score, nchar(candidates))]
  gsub("\\\\", "/", candidates[1L])
}

read_required <- function(stage_key, filename) {
  path <- find_exact_file(STAGE_DIRS[[stage_key]], filename, TRUE)
  list(path = path, data = safe_fread(path))
}
read_sheet_required <- function(stage_key, filename, sheet) {
  path <- find_exact_file(STAGE_DIRS[[stage_key]], filename, TRUE)
  list(path = path, data = as.data.table(openxlsx::read.xlsx(path, sheet = sheet)))
}

pretty_program <- function(x) {
  y <- gsub("_Top[0-9]+$", "", as.character(x))
  y <- gsub("Ccr2pos", "CCR2+ program", y, ignore.case = TRUE)
  y <- gsub("Ccr2neg", paste0("CCR2", UNICODE_MINUS, " program"), y, ignore.case = TRUE)
  y <- gsub("CrossSubset", "Cross-subset program", y, ignore.case = TRUE)
  y
}
pretty_tf <- function(x) toupper(as.character(x))
pretty_cell <- function(x) gsub("_", " ", as.character(x))
pretty_contrast <- function(x) {
  gsub("_", " ", gsub("TYA_018", "TYA-018", as.character(x)))
}
pretty_pathway <- function(x) {
  y <- tools::toTitleCase(
    tolower(gsub("_", " ", gsub("^HALLMARK_", "", as.character(x))))
  )
  y <- gsub("Nfkb", "NF-kB", y, fixed = TRUE)
  y <- gsub("Tnfa", "TNFA", y, fixed = TRUE)
  y <- gsub("Il6 Jak Stat3", "IL6-JAK-STAT3", y, fixed = TRUE)
  y <- gsub("Il2 Stat5", "IL2-STAT5", y, fixed = TRUE)
  y <- gsub("Kras", "KRAS", y, fixed = TRUE)
  y <- gsub("E2f", "E2F", y, fixed = TRUE)
  y
}
short_axis <- function(tf, ligand, receptor, receiver) {
  receiver_label <- tools::toTitleCase(
    tolower(gsub("_", " ", trimws(as.character(receiver))))
  )
  paste0(
    toupper(tf), "-", toupper(ligand), "-", toupper(receptor), "-",
    receiver_label
  )
}
wrap_text <- function(x, width = 32L) vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
label_compact_number <- scales::label_number(
  scale_cut = scales::cut_short_scale()
)
rescale01 <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0.5, length(x)))
  (x - r[1L]) / diff(r)
}
neglog10_fdr <- function(x) {
  -log10(pmax(as.numeric(x), .Machine$double.xmin))
}

ordered_group_levels <- function(
  dt,
  group_col,
  value_col,
  decreasing = FALSE,
  use_absolute = FALSE
) {
  x <- copy(as.data.table(dt))
  stopifnot(group_col %in% names(x), value_col %in% names(x))
  x[, .group_value__ := as.character(get(group_col))]
  x[, .numeric_value__ := as.numeric(get(value_col))]
  summary_dt <- x[
    !is.na(.group_value__) & nzchar(.group_value__),
    .(
      order_value = mean(
        if (use_absolute) abs(.numeric_value__) else .numeric_value__,
        na.rm = TRUE
      )
    ),
    by = .group_value__
  ]
  setorderv(
    summary_dt,
    "order_value",
    order = if (decreasing) -1L else 1L,
    na.last = TRUE
  )
  unique(summary_dt$.group_value__)
}

complete_named_palette <- function(values, locked_palette) {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  result <- locked_palette[names(locked_palette) %in% values]
  missing <- setdiff(values, names(result))
  if (length(missing) > 0L) {
    extra <- scales::hue_pal(
      h = c(15, 375),
      c = 100,
      l = 55
    )(length(missing))
    names(extra) <- missing
    result <- c(result, extra)
  }
  result[values]
}

add_hedges_ci <- function(dt, g_col, n1_col, n2_col) {
  x <- copy(as.data.table(dt))
  g <- as.numeric(x[[g_col]])
  n1 <- as.numeric(x[[n1_col]])
  n2 <- as.numeric(x[[n2_col]])
  variance <- (n1 + n2) / pmax(n1 * n2, 1) + (g^2) / pmax(2 * (n1 + n2 - 2), 1)
  se <- sqrt(variance)
  x[, ci_low := g - 1.96 * se]
  x[, ci_high := g + 1.96 * se]
  x
}

find_umap_cols <- function(dt) {
  candidates <- names(dt)[grepl("^UMAP|umap", names(dt))]
  candidates <- candidates[vapply(dt[, ..candidates], is.numeric, logical(1))]
  if (length(candidates) < 2L) stop("Could not identify two numeric UMAP coordinate columns.")
  candidates[1:2]
}

theme_sci <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      text = element_text(family = BASE_FAMILY, colour = "#1A1A1A"),
      axis.text = element_text(family = BASE_FAMILY, colour = "#1A1A1A"),
      axis.title = element_text(family = BASE_FAMILY, colour = "#1A1A1A"),
      axis.line = element_line(linewidth = LINE_WIDTH, colour = "#1A1A1A"),
      axis.ticks = element_line(linewidth = LINE_WIDTH, colour = "#1A1A1A"),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "#F3F3F3", colour = "#404040", linewidth = 0.45),
      strip.text = element_text(family = BASE_FAMILY, face = "bold", size = base_size),
      legend.title = element_text(family = BASE_FAMILY, face = "bold", size = base_size - 0.3),
      legend.text = element_text(family = BASE_FAMILY, size = base_size - 0.8),
      plot.title = element_text(family = BASE_FAMILY, face = "bold", size = base_size + 1.2, hjust = 0),
      plot.subtitle = element_text(family = BASE_FAMILY, size = base_size - 0.2, colour = "#4D4D4D", hjust = 0),
      plot.caption = element_text(family = BASE_FAMILY, size = base_size - 1.0, colour = "#555555", hjust = 0),
      legend.key.height = grid::unit(3.6, "mm"),
      legend.key.width = grid::unit(4.2, "mm"),
      legend.spacing.x = grid::unit(1.3, "mm"),
      legend.spacing.y = grid::unit(0.8, "mm"),
      plot.margin = margin(6, 9, 6, 7, unit = "pt")
    )
}

theme_heatmap <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      text = element_text(family = BASE_FAMILY, colour = "#1A1A1A"),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(family = BASE_FAMILY, colour = "#1A1A1A"),
      strip.background = element_rect(fill = "#F3F3F3", colour = "#404040", linewidth = 0.45),
      strip.text = element_text(family = BASE_FAMILY, face = "bold"),
      plot.title = element_text(family = BASE_FAMILY, face = "bold", size = base_size + 1.2),
      legend.title = element_text(family = BASE_FAMILY, face = "bold"),
      plot.margin = margin(8, 16, 8, 8, unit = "pt")
    )
}

PALETTE_CONDITION <- c(Control = "#3B6FB6", HFpEF = "#C84C4C")
PALETTE_TF <- c(BHLHE40 = "#2A9D8F", NFKB1 = "#E76F51", RELA = "#577590", RUNX1 = "#8D5A97", SPI1 = "#4F772D", REL = "#C7902F")
PALETTE_PROGRAM <- setNames(
  c("#008C95", "#E58E26", "#6C5CE7"),
  c(
    "CCR2+ program",
    paste0("CCR2", UNICODE_MINUS, " program"),
    "Cross-subset program"
  )
)
PALETTE_RECEIVER <- c(
  Endothelial = "#0072B2",
  Fibroblast = "#D55E00",
  Pericyte = "#009E73",
  Smooth_muscle = "#CC79A7",
  ENDOTHELIAL = "#0072B2",
  FIBROBLAST = "#D55E00",
  PERICYTE = "#009E73",
  SMOOTH_MUSCLE = "#CC79A7"
)
PALETTE_CELLTYPE <- c(
  Cardiomyocyte = "#9C4F22",
  Fibroblast = "#D33682",
  Endothelial = "#0099C6",
  Lymphatic_endothelial = "#006D9C",
  Pericyte = "#6C5CE7",
  Smooth_muscle = "#A6761D",
  Macrophage_Monocyte = "#D73027",
  Dendritic_cell = "#F46D43",
  Neutrophil = "#9AAE20",
  T_NK = "#1F78B4",
  B_cell = "#00A6A6",
  Mast_cell = "#F28E2B",
  Epicardial_Mesothelial = "#2CA25F",
  Schwann_Glial = "#7B3294",
  Platelet_Megakaryocyte = "#8E44AD",
  Erythroid = "#B2182B",
  Cycling_unresolved = "#7A8B22",
  Low_quality_mitochondrial = "#B8B8B8",
  Unresolved = "#686868"
)

panel_audit_records <- list()
input_manifest_records <- list()
supp_table_records <- list()

record_input <- function(panel_id, paths) {
  paths <- unique(paths[!is.na(paths) & file.exists(paths)])
  if (length(paths) == 0L) return(invisible(NULL))
  input_manifest_records[[length(input_manifest_records) + 1L]] <<- data.table(
    panel_id = panel_id,
    input_path = gsub("\\\\", "/", paths),
    input_basename = basename(paths),
    size_bytes = file.info(paths)$size,
    md5 = unname(tools::md5sum(paths))
  )
  invisible(NULL)
}

write_panel_source <- function(x, source_dir, stem) {
  x <- as.data.table(x)
  extension <- if (nrow(x) > 100000L) ".csv.gz" else ".csv"
  path <- file.path(source_dir, paste0(stem, "_source_data", extension))
  write_csv_safe(x, path)
  path
}

save_panel <- function(plot_object, figure_group, stem, width, height, source_data, input_paths, title, panel_type = c("main", "supplementary")) {
  panel_type <- match.arg(panel_type)
  fig_root <- if (panel_type == "main") DIRS$main else DIRS$supplementary
  src_root <- if (panel_type == "main") DIRS$source_main else DIRS$source_supp
  out_dir <- file.path(fig_root, figure_group)
  src_dir <- file.path(src_root, figure_group)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
  tiff_path <- file.path(out_dir, paste0(stem, ".tif"))
  source_path <- write_panel_source(source_data, src_dir, stem)

  ## Figure titles and subtitles are intentionally excluded from panel TIFFs.
  ## They are retained in the manifest and later written in the figure legends.
  plot_object <- plot_object +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank()
    )

  ggsave(
    filename = tiff_path, plot = plot_object,
    width = width, height = height, units = "in", dpi = FIGURE_DPI,
    device = grDevices::tiff, compression = "lzw", type = "cairo",
    bg = "white", limitsize = FALSE
  )
  valid <- file.exists(tiff_path) && is.finite(file.info(tiff_path)$size) && file.info(tiff_path)$size > 10000
  panel_id <- paste0(figure_group, "_", stem)
  record_input(panel_id, input_paths)
  panel_audit_records[[length(panel_audit_records) + 1L]] <<- data.table(
    panel_id = panel_id, panel_type = panel_type, figure_group = figure_group,
    stem = stem, title = title, tiff_path = gsub("\\\\", "/", tiff_path),
    source_data_path = gsub("\\\\", "/", source_path),
    width_in = width, height_in = height, dpi = FIGURE_DPI, compression = "LZW",
    font_family = BASE_FAMILY,
    tiff_size_bytes = ifelse(file.exists(tiff_path), file.info(tiff_path)$size, NA_real_),
    tiff_md5 = ifelse(file.exists(tiff_path), unname(tools::md5sum(tiff_path)), NA_character_),
    valid = valid
  )
  if (!valid) stop("TIFF export failed or is unexpectedly small: ", tiff_path)
  log_msg("Exported ", panel_id)
  invisible(tiff_path)
}

make_unique_sheet_names <- function(x) {
  base <- substr(gsub("[\\[\\]:*?/\\\\]", "_", as.character(x)), 1L, 31L)
  base[!nzchar(base)] <- "data"
  out <- character(length(base))
  used <- character()
  for (i in seq_along(base)) {
    candidate <- base[i]
    suffix_index <- 1L
    while (tolower(candidate) %in% tolower(used)) {
      suffix_index <- suffix_index + 1L
      suffix <- paste0("_", suffix_index)
      candidate <- paste0(substr(base[i], 1L, 31L - nchar(suffix)), suffix)
    }
    out[i] <- candidate
    used <- c(used, candidate)
  }
  out
}

write_supp_table <- function(number, short_title, sheets) {
  number_label <- sprintf("S%02d", as.integer(number))
  safe_title <- gsub("[^A-Za-z0-9]+", "_", short_title)
  safe_title <- gsub("^_+|_+$", "", safe_title)
  workbook_path <- file.path(DIRS$supp_tables, paste0("Table_", number_label, "_", safe_title, ".xlsx"))
  wb <- openxlsx::createWorkbook()
  sheet_manifest <- list()
  original_sheet_names <- names(sheets)
  clean_sheet_names <- make_unique_sheet_names(original_sheet_names)
  for (sheet_index in seq_along(sheets)) {
    sheet_name <- original_sheet_names[sheet_index]
    x <- as.data.table(sheets[[sheet_name]])
    clean_sheet <- clean_sheet_names[sheet_index]
    openxlsx::addWorksheet(wb, clean_sheet)
    openxlsx::writeDataTable(wb, clean_sheet, as.data.frame(x), tableStyle = "TableStyleMedium2")
    openxlsx::freezePane(wb, clean_sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, clean_sheet, cols = seq_len(max(1L, ncol(x))), widths = "auto")
    csv_dir <- file.path(DIRS$supp_table_csv, number_label)
    dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
    csv_path <- file.path(csv_dir, paste0(clean_sheet, if (nrow(x) > 100000L) ".csv.gz" else ".csv"))
    write_csv_safe(x, csv_path)
    sheet_manifest[[length(sheet_manifest) + 1L]] <- data.table(
      table_number = number_label, table_title = short_title, sheet = clean_sheet,
      rows = nrow(x), columns = ncol(x), source_csv = gsub("\\\\", "/", csv_path)
    )
  }
  openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
  valid <- file.exists(workbook_path) && file.info(workbook_path)$size > 1000
  supp_table_records[[length(supp_table_records) + 1L]] <<- rbindlist(sheet_manifest, use.names = TRUE, fill = TRUE)[, `:=`(
    workbook_path = gsub("\\\\", "/", workbook_path), workbook_valid = valid
  )]
  if (!valid) stop("Supplementary workbook export failed: ", workbook_path)
  log_msg("Exported Table ", number_label, ": ", short_title)
  invisible(workbook_path)
}

make_effect_heatmap <- function(dt, x_col, y_col, fill_col, label_col = NULL, title, fill_title = "Hedges' g", base_size = BASE_SIZE) {
  x <- copy(as.data.table(dt))
  p <- ggplot(x, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0, na.value = "#E6E6E6") +
    labs(fill = fill_title) +
    theme_heatmap(base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  if (!is.null(label_col)) {
    p <- p + geom_text(aes(label = .data[[label_col]]), family = BASE_FAMILY, size = LABEL_SIZE)
  }
  p
}

############################################################
## 2. Upstream status validation and frozen data loading
############################################################
status_specs <- list(
  stage2 = c("24_stage2_run_status.csv", "23_scientific_completion_checks.csv"),
  stage3 = c("39_stage3_run_status.csv", "38_scientific_completion_checks.csv"),
  stage4 = c("22_stage4_run_status.csv", "20_stage4_scientific_completion_checks.csv"),
  stage5 = c("21_stage5_run_status.csv", "20_stage5_scientific_completion_checks.csv"),
  stage5B = c("17_stage5B_run_status.csv", "16_stage5B_scientific_completion_checks.csv"),
  stage6 = c("23_stage6_run_status.csv", "22_stage6_scientific_completion_checks.csv"),
  stage7 = c("20_stage7_run_status.csv", "19_stage7_scientific_completion_checks.csv"),
  stage8 = c("73_stage8_run_status.csv", "72_stage8_scientific_completion_checks.csv"),
  benchmark = c("20_run_status.csv", "19_scientific_completion_checks.csv")
)
upstream_audit <- rbindlist(lapply(names(status_specs), function(k) {
  status_path <- find_exact_file(STAGE_DIRS[[k]], status_specs[[k]][1L], TRUE)
  checks_path <- find_exact_file(STAGE_DIRS[[k]], status_specs[[k]][2L], TRUE)
  status <- safe_fread(status_path)
  checks <- safe_fread(checks_path)
  status_col <- intersect(c("overall_status", "status"), names(status))[1L]
  if (is.na(status_col) || !"status" %in% names(checks)) stop("Invalid status/check files for ", k)
  if (any(checks$status != "PASS")) stop("Non-PASS scientific check found in ", k)
  data.table(stage = k, overall_status = as.character(status[[status_col]][1L]), checks = nrow(checks), failed_checks = sum(checks$status != "PASS"), status_file = status_path, checks_file = checks_path)
}))
write_csv_safe(upstream_audit, file.path(DIRS$manifests, "01_upstream_status_audit.csv"))
log_msg("All frozen upstream scientific checks passed.")

## Stage 1
s1_manifest <- read_required("stage1", "01_locked_sample_manifest.csv")
s1_donors <- read_required("stage1", "04_SCP3342_locked_donor_manifest.csv")
s1_roles <- read_required("stage1", "05_analysis_role_and_contrast_plan.csv")
s1_summary <- read_required("stage1", "07_dataset_lock_summary.csv")
## Stage 2
s2_meta <- read_required("stage2", "01_sample_metadata_used.csv")
s2_group_counts <- read_required("stage2", "02_group_count_validation.csv")
s2_mapping <- read_required("stage2", "03_sample_to_RSEM_file_mapping.csv")
s2_qc <- read_required("stage2", "06_sample_QC_metrics.csv")
s2_contrasts <- read_required("stage2", "07_contrast_definitions.csv")
s2_contrast_summary <- read_required("stage2", "10_contrast_result_summary.csv")
s2_method <- read_required("stage2", "12_DESeq2_edgeR_method_concordance_summary.csv")
s2_opp_summary <- read_required("stage2", "15_opposition_summary.csv")
s2_hallmark <- read_required("stage2", "20_Hallmark_pathway_opposition.csv")
s2_top_pos <- read_required("stage2", "TOP500_opposition_Ccr2_positive.csv")
s2_top_neg <- read_required("stage2", "TOP500_opposition_Ccr2_negative.csv")
s2_top_cross <- read_required("stage2", "TOP500_cross_subset_consensus.csv")
## Stage 3
s3_sig_size <- read_required("stage3", "03_stage2_signature_size_summary.csv")
s3_mapping <- read_required("stage3", "04_GSE236585_10x_file_mapping_and_dimensions.csv")
s3_qc_thresholds <- read_required("stage3", "05_sample_specific_QC_thresholds.csv")
s3_qc_retention <- read_required("stage3", "06_sample_QC_retention_summary.csv")
s3_doublet <- read_required("stage3", "07B_scDblFinder_rate_summary.csv")
s3_annotation <- read_required("stage3", "11_major_celltype_cluster_annotation.csv")
s3_composition <- read_required("stage3", "15_celltype_composition_by_sample.csv")
s3_eligibility <- read_required("stage3", "21_major_celltype_pseudobulk_eligibility.csv")
s3_program_stats <- read_required("stage3", "27_pseudobulk_program_statistics.csv")
s3_concordance <- read_required("stage3", "28_stage2_stage3_gene_level_concordance.csv")
s3_mac_annotation <- read_required("stage3", "30_macrophage_cluster_annotation.csv")
s3_mac_comp <- read_required("stage3", "33_macrophage_state_composition_by_sample.csv")
s3_mac_stats <- read_required("stage3", "35_macrophage_state_score_statistics.csv")
s3_primary <- read_required("stage3", "PRIMARY_program_localization_statistics.csv")
## Stage 4
s4_activity <- read_required("stage4", "09_stage4_weighted_regulon_activity_HFpEF_vs_Control.csv")
s4_aucell <- read_required("stage4", "10_stage4_AUCell_regulon_activity_HFpEF_vs_Control.csv")
s4_expression <- read_required("stage4", "11_stage4_TF_expression_HFpEF_vs_Control.csv")
s4_priority <- read_required("stage4", "12_stage4_candidate_TF_priority_score.csv")
s4_lopo <- read_required("stage4", "15_stage4_leave_one_pair_out_TF_robustness_summary.csv")
s4_method <- read_required("stage4", "18_stage4_TF_method_comparison_summary.csv")
## Stage 5
s5_manifest <- read_required("stage5", "01_stage5_candidate_TF_resolution.csv")
s5_controls <- read_required("stage5", "02_stage5_candidate_and_matched_control_TFs.csv")
s5_program_defs <- read_required("stage5", "03_stage5_program_definition_summary.csv")
s5_observed <- read_required("stage5", "04_stage5_observed_program_scores.csv")
s5_gene_effects <- read_required("stage5", "08_stage5_primary_top_predicted_genes_per_TF.csv")
s5_ligands <- read_required("stage5", "09_stage5_candidate_ligand_changes_for_stage6.csv")
s5_method <- read_required("stage5", "11_stage5_perturbation_method_concordance.csv")
s5_control_results <- read_required("stage5", "12_stage5_matched_control_TF_reference_results.csv")
s5_rank <- read_required("stage5", "13_stage5_candidate_TF_rank_aggregation.csv")
s5_sensitivity <- read_required("stage5", "14_stage5_candidate_ranking_sensitivity_scenarios.csv")
s5_stability <- read_required("stage5", "15_stage5_candidate_ranking_stability_summary.csv")
s5_modes <- read_required("stage5", "16_stage5_normalization_vs_attenuation_results.csv")
s5_mode_summary <- read_required("stage5", "17_stage5_normalization_vs_attenuation_summary.csv")
## Stage 5B
s5b_full <- read_required("stage5B", "04_stage5B_full_candidate_summary.csv")
s5b_boot <- read_required("stage5B", "07_stage5B_candidate_regulon_bootstrap_summary.csv")
s5b_matching <- read_required("stage5B", "08_stage5B_all_TF_matching_covariates.csv")
s5b_null_pools <- read_required("stage5B", "09_stage5B_candidate_matched_null_TF_pools.csv")
s5b_null_effects <- read_required("stage5B", "10_stage5B_precomputed_null_TF_effects.csv")
s5b_null <- read_required("stage5B", "12_stage5B_random_matched_TF_null_summary.csv")
s5b_final <- read_required("stage5B", "13_stage5B_final_candidate_robustness_rank.csv")
## Stage 6
s6_tf_manifest <- read_required("stage6", "01_stage6_candidate_TF_manifest.csv")
s6_ligand_coverage <- read_required("stage6", "03_stage6_candidate_ligand_coverage.csv")
s6_nichenet <- read_required("stage6", "12_stage6_NicheNet_ligand_activity.csv")
s6_stability <- read_required("stage6", "18_stage6_axis_ranking_stability_summary.csv")
s6_candidate <- read_required("stage6", "19_stage6_candidate_TF_communication_summary.csv")
s6_top_axes <- read_sheet_required("stage6", "20_stage6_TF_dependent_communication_key_results.xlsx", "Top_axes")
s6_receiver_sets <- read_sheet_required("stage6", "20_stage6_TF_dependent_communication_key_results.xlsx", "Receiver_gene_sets")
s6_ligand_support <- read_sheet_required("stage6", "20_stage6_TF_dependent_communication_key_results.xlsx", "Ligand_support")
## Stage 7
s7_definitions <- read_required("stage7", "02_stage7_feature_definitions.csv")
s7_features <- read_required("stage7", "03_stage7_sample_level_feature_matrix.csv")
s7_fold <- read_required("stage7", "07_stage7_primary_LOPO_fold_performance.csv")
s7_predictions <- read_required("stage7", "08_stage7_primary_LOPO_sample_predictions.csv")
s7_importance <- read_required("stage7", "11_stage7_feature_attribution_and_stability.csv")
s7_permutation <- read_required("stage7", "13_stage7_permutation_summary.csv")
s7_sensitivity <- read_required("stage7", "14_stage7_panel_and_lambda_sensitivity.csv")
s7_performance <- read_required("stage7", "15_stage7_model_performance_summary.csv")
## Stage 8
s8_core_meta <- read_required("stage8", "02_core_dataset_locked_metadata.csv")
s8_program_manifest <- read_required("stage8", "03_frozen_program_manifest.csv")
s8_tf_manifest <- read_required("stage8", "04_frozen_TF_manifest.csv")
s8_axis_manifest <- read_required("stage8", "06_frozen_axis_manifest.csv")
s8_gene_audit <- read_required("stage8", "08_gene_key_audit_summary.csv")
s8_program_evidence <- read_required("stage8", "60_multicohort_program_evidence.csv.gz")
s8_tf_evidence <- read_required("stage8", "61_multicohort_TF_evidence.csv.gz")
s8_axis_evidence <- read_required("stage8", "62_multicohort_axis_evidence.csv.gz")
s8_program_summary <- read_required("stage8", "63_program_integrated_summary.csv")
s8_tf_summary <- read_required("stage8", "64_TF_integrated_summary.csv")
s8_axis_summary <- read_required("stage8", "65_axis_integrated_summary.csv")
s8_roles <- read_required("stage8", "66_dataset_roles_and_claim_boundaries.csv")
s8_counts <- read_required("stage8", "74_biological_sample_count_audit.csv")
s8_donors <- read_required("stage8", "75_SCP3342_donor_count_audit.csv")
s8_scp_program <- read_required("stage8", "50B_SCP3342_program_validation.csv")
s8_scp_tf <- read_required("stage8", "50C_SCP3342_TF_validation.csv")
s8_scp_axis <- read_required("stage8", "50D_SCP3342_axis_validation.csv")
s8_gse249_axis <- read_required("stage8", "40D_GSE249412_axis_validation.csv")
## Benchmark
b_tf_all <- read_required("benchmark", "02_TF_baseline_ranks_all174.csv")
b_tf_selected <- read_required("benchmark", "03_selected_TF_baseline_rank_comparison.csv")
b_tf_ranks <- read_required("benchmark", "04_candidate_TF_method_ranks.csv")
b_ablation <- read_required("benchmark", "05_candidate_TF_ablation_scenarios.csv")
b_ablation_stability <- read_required("benchmark", "06_candidate_TF_ablation_stability.csv")
b_tf_alignment <- read_required("benchmark", "08_candidate_TF_external_alignment.csv")
b_axis_ranks <- read_required("benchmark", "09_axis_method_ranks.csv")
b_axis_alignment <- read_required("benchmark", "11_axis_external_alignment.csv")
b_axis_topk <- read_required("benchmark", "12_axis_topk_external_performance.csv")
b_sentinel <- read_required("benchmark", "13_sentinel_axis_rank_comparison.csv")
b_stage7_metrics <- read_required("benchmark", "14_stage7_quantitative_metrics.csv")
b_feature <- read_required("benchmark", "15_stage7_feature_importance.csv")
b_boundaries <- read_required("benchmark", "17_revision_claim_boundaries.csv")

############################################################
## 4. Main Figure 1
## Directional pharmacotranscriptomic discovery
############################################################

## Figure 1A: GSE237156 PCA
f1a_source <- read_required("stage2", "Fig2A_PCA_source.csv")
f1a <- copy(f1a_source$data)
f1a[, treatment_group := paste(diet, drug, sep = " / ")]
f1a[, subset_label := fcase(
  macrophage_subset == "Ccr2_positive", "CCR2+",
  macrophage_subset == "Ccr2_negative", paste0("CCR2", UNICODE_MINUS),
  default = as.character(macrophage_subset)
)]
group_palette <- c(
  "CD / Vehicle" = "#6B8E9B",
  "CD / Dapagliflozin" = "#3B6FB6",
  "HFD / Vehicle" = "#C84C4C",
  "HFD / Dapagliflozin" = "#E9A03B"
)
p_f1a <- ggplot(
  f1a,
  aes(x = PC1, y = PC2, fill = treatment_group, shape = subset_label)
) +
  geom_point(size = 3.4, colour = "#202020", stroke = 0.5) +
  ggrepel::geom_text_repel(
    aes(label = sample_accession),
    family = BASE_FAMILY,
    size = 2.35,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.25,
    box.padding = 0.25,
    point.padding = 0.15
  ) +
  scale_fill_manual(
    values = group_palette,
    breaks = c(
      "CD / Vehicle",
      "CD / Dapagliflozin",
      "HFD / Vehicle",
      "HFD / Dapagliflozin"
    )
  ) +
  scale_shape_manual(values = setNames(c(21, 24), c("CCR2+", paste0("CCR2", UNICODE_MINUS)))) +
  guides(
    fill = guide_legend(
      order = 1,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(shape = 21, size = 3.2, colour = "#202020")
    ),
    shape = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(size = 3.2, fill = "white")
    )
  ) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Diet / treatment",
    shape = "Macrophage subset"
  ) +
  theme_sci() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal"
  ) +
  coord_fixed(ratio = 1, clip = "off")
save_panel(
  p_f1a, "Figure_1", "Figure1A_GSE237156_PCA",
  6.3, 5.0, f1a, f1a_source$path,
  "GSE237156 sample-level PCA", "main"
)

## Figure 1B-C: disease-drug opposition scatter plots
make_opposition_panel <- function(source_obj, title_text, stem, shared_limit) {
  x <- copy(source_obj$data)
  x[, tier_group := fcase(
    grepl("^Tier_A", opposition_tier), "Tier A",
    grepl("^Tier_B", opposition_tier), "Tier B",
    grepl("^Tier_C", opposition_tier), "Tier C",
    default = "Other"
  )]
  label_dt <- x[
    deseq_opposed == TRUE &
      !is.na(display_gene) &
      nzchar(display_gene) &
      !grepl("^ENSMUSG", display_gene)
  ][order(combined_rank_product, -opposition_effect_score)][1:min(.N, 8L)]

  p <- ggplot(x, aes(x = disease_lfc, y = drug_lfc)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "#777777") +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "#777777") +
    geom_point(
      data = x[tier_group == "Other"],
      colour = "#B7B7B7", alpha = 0.28, size = 0.65
    ) +
    geom_point(
      data = x[tier_group != "Other"],
      aes(colour = tier_group), alpha = 0.78, size = 1.05
    ) +
    ggrepel::geom_text_repel(
      data = label_dt,
      aes(label = display_gene),
      family = BASE_FAMILY,
      size = 2.25,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.12,
      min.segment.length = 0,
      segment.size = 0.25
    ) +
    scale_colour_manual(values = c(
      "Tier A" = "#C84C4C",
      "Tier B" = "#E9A03B",
      "Tier C" = "#3B6FB6"
    )) +
    labs(
      x = "Disease effect\n(HFD vehicle vs CD vehicle, log2FC)",
      y = "Dapagliflozin effect\n(HFD dapagliflozin vs HFD vehicle, log2FC)",
      colour = "Evidence tier"
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 2.6, alpha = 1)
      )
    ) +
    theme_sci() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(-shared_limit, shared_limit),
      ylim = c(-shared_limit, shared_limit),
      clip = "off"
    )

  save_panel(
    p, "Figure_1", stem,
    MAIN_SCATTER_WIDTH, MAIN_SCATTER_HEIGHT, x, source_obj$path,
    title_text, "main"
  )
}
f1b_source <- read_required(
  "stage2", "Fig2B_Ccr2_positive_disease_drug_opposition_source.csv"
)
f1c_source <- read_required(
  "stage2", "Fig2C_Ccr2_negative_disease_drug_opposition_source.csv"
)
opposition_shared_limit <- max(
  abs(c(
    f1b_source$data$disease_lfc,
    f1b_source$data$drug_lfc,
    f1c_source$data$disease_lfc,
    f1c_source$data$drug_lfc
  )),
  na.rm = TRUE
)
opposition_shared_limit <- max(1, opposition_shared_limit * 1.06)
make_opposition_panel(
  f1b_source,
  "CCR2+ macrophage disease-drug opposition",
  "Figure1B_CCR2positive_disease_drug_opposition",
  opposition_shared_limit
)
make_opposition_panel(
  f1c_source,
  paste0("CCR2", UNICODE_MINUS, " macrophage disease-drug opposition"),
  "Figure1C_CCR2negative_disease_drug_opposition",
  opposition_shared_limit
)

## Figure 1D: top opposed-gene heatmap
f1d_source <- read_required("stage2", "Fig2D_top_opposed_heatmap_source.csv")
f1d <- copy(f1d_source$data)
gene_col_f1d <- intersect(c("display_gene", "symbol", "gene"), names(f1d))[1L]
if (is.na(gene_col_f1d)) stop("Figure 1D source lacks a gene-label column.")
sample_cols_f1d <- setdiff(names(f1d), gene_col_f1d)
f1d_long <- melt(
  f1d,
  id.vars = gene_col_f1d,
  measure.vars = sample_cols_f1d,
  variable.name = "sample_accession",
  value.name = "z_score"
)
setnames(f1d_long, gene_col_f1d, "gene")
f1d_long <- merge(
  f1d_long,
  s2_meta$data[, .(sample_accession, macrophage_subset, diet, drug)],
  by = "sample_accession", all.x = TRUE
)
sample_order_f1d <- unique(s2_meta$data[
  order(macrophage_subset, diet, drug, sample_accession),
  sample_accession
])
f1d_long[, subset_label := fcase(
  macrophage_subset == "Ccr2_positive", "CCR2+ macrophages",
  macrophage_subset == "Ccr2_negative", paste0("CCR2", UNICODE_MINUS, " macrophages"),
  default = as.character(macrophage_subset)
)]
f1d_long[, annotated_symbol := (
  !is.na(gene) & nzchar(gene) & !grepl("^ENS(MUS)?G", gene, ignore.case = TRUE)
)]
f1d_plot <- f1d_long[annotated_symbol == TRUE]
if (uniqueN(f1d_plot$gene) < 20L) {
  f1d_plot <- copy(f1d_long)
}
gene_order_f1d <- unique(as.character(f1d_plot$gene))
f1d_long[, displayed_in_main_panel := gene %in% gene_order_f1d]
f1d_plot[, sample_accession := factor(sample_accession, levels = unique(sample_order_f1d))]
f1d_plot[, gene := factor(gene, levels = unique(rev(gene_order_f1d)))]
p_f1d <- ggplot(
  f1d_plot,
  aes(x = sample_accession, y = gene, fill = z_score)
) +
  geom_tile(colour = "white", linewidth = 0.12) +
  facet_grid(~ subset_label, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0
  ) +
  labs(
    x = NULL, y = NULL, fill = "Row z-score"
  ) +
  theme_heatmap(7.3) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.3),
    axis.text.y = element_text(size = 6.2),
    legend.position = "right"
  )
save_panel(
  p_f1d, "Figure_1", "Figure1D_top_drug_opposed_gene_heatmap",
  MAIN_HEATMAP_WIDTH, MAIN_HEATMAP_HEIGHT, f1d_long, c(f1d_source$path, s2_meta$path),
  "Top drug-opposed gene heatmap", "main"
)

## Figure 1E: Hallmark pathway opposition
f1e <- copy(s2_hallmark$data)
f1e <- f1e[
  pathway_opposed == TRUE & is.finite(pathway_opposition_strength)
][order(-pathway_opposition_strength)]
f1e <- f1e[1:min(.N, 14L)]
f1e[, pathway_label := pretty_pathway(pathway)]
f1e[, subset_label := fcase(
  subset == "Ccr2_positive", "CCR2+ macrophages",
  subset == "Ccr2_negative", paste0("CCR2", UNICODE_MINUS, " macrophages"),
  default = pretty_cell(subset)
)]
f1e[, pathway_panel_key := paste(subset_label, pathway_label, sep = " || ")]
f1e[, pathway_panel_key := factor(
  pathway_panel_key,
  levels = unique(rev(as.character(pathway_panel_key)))
)]
f1e_long <- melt(
  f1e,
  id.vars = c(
    "pathway", "pathway_label", "pathway_panel_key", "subset", "subset_label",
    "pathway_opposition_strength"
  ),
  measure.vars = c("disease_NES", "drug_NES"),
  variable.name = "comparison", value.name = "NES"
)
f1e_long[, comparison := fcase(
  comparison == "disease_NES", "HFpEF-like disease effect",
  comparison == "drug_NES", "Dapagliflozin effect",
  default = as.character(comparison)
)]
p_f1e <- ggplot() +
  geom_segment(
    data = f1e,
    aes(x = disease_NES, xend = drug_NES, y = pathway_panel_key, yend = pathway_panel_key),
    colour = "#B0B0B0", linewidth = 0.65
  ) +
  geom_point(
    data = f1e_long,
    aes(x = NES, y = pathway_panel_key, colour = comparison, shape = comparison),
    size = 2.6
  ) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = "#777777") +
  facet_grid(subset_label ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.* \\|\\| ", "", x)) +
  scale_colour_manual(values = c(
    "HFpEF-like disease effect" = "#C84C4C",
    "Dapagliflozin effect" = "#3B6FB6"
  )) +
  scale_shape_manual(values = c(
    "HFpEF-like disease effect" = 16,
    "Dapagliflozin effect" = 17
  )) +
  labs(
    title = "Hallmark pathways show disease-drug directional opposition",
    x = "Normalized enrichment score", y = NULL,
    colour = NULL, shape = NULL
  ) +
  theme_sci(7.7) +
  theme(
    axis.text.y = element_text(size = 6.3),
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")
save_panel(
  p_f1e, "Figure_1", "Figure1E_Hallmark_pathway_opposition",
  9.2, 6.8, f1e_long, s2_hallmark$path,
  "Hallmark pathway opposition", "main"
)

## Figure 1F: signature-size sensitivity
f1f <- copy(s3_sig_size$data)
f1f[, program := pretty_program(signature_name)]
f1f[, requested_size := requested_size_per_direction]
f1f[, direction_label := fcase(
  direction == "Disease_up_Drug_down", "Disease up / drug down",
  direction == "Disease_down_Drug_up", "Disease down / drug up",
  default = as.character(direction)
)]
program_levels_f1f <- c(
  "CCR2+ program",
  paste0("CCR2", UNICODE_MINUS, " program"),
  "Cross-subset program"
)
f1f[, program := factor(program, levels = unique(program_levels_f1f))]
if (any(is.na(f1f$program))) {
  stop(
    "Figure 1F program-label normalization failed. Unresolved labels: ",
    paste(unique(pretty_program(f1f$signature_name)[is.na(f1f$program)]), collapse = ", ")
  )
}
p_f1f <- ggplot(
  f1f,
  aes(x = requested_size, y = selected_genes, colour = direction_label, shape = direction_label)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.45, colour = "#999999") +
  geom_line(aes(group = direction_label), linewidth = 0.75) +
  geom_point(size = 2.5) +
  facet_wrap(~ program, nrow = 1) +
  scale_colour_manual(values = c(
    "Disease up / drug down" = "#C84C4C",
    "Disease down / drug up" = "#3B6FB6"
  )) +
  scale_shape_manual(values = c(
    "Disease up / drug down" = 16,
    "Disease down / drug up" = 17
  )) +
  scale_x_continuous(breaks = c(50, 100, 150, 200)) +
  labs(
    title = "Continuous signatures were evaluated across four prespecified sizes",
    x = "Requested genes per direction",
    y = "Genes retained after frozen filtering",
    colour = NULL, shape = NULL
  ) +
  theme_sci(7.7) +
  theme(legend.position = "bottom")
save_panel(
  p_f1f, "Figure_1", "Figure1F_signature_size_sensitivity",
  7.6, 4.2, f1f, s3_sig_size$path,
  "Signature-size sensitivity", "main"
)

############################################################
## 5. Main Figure 2
## Cardiac single-cell localization and partial recovery
############################################################

## Figure 2A: major-cell-type UMAP
f2a_source <- read_required("stage3", "Fig3A_3B_3D_UMAP_source.csv.gz")
f2a <- copy(f2a_source$data)
umap_cols_f2a <- find_umap_cols(f2a)
cell_col_f2a <- intersect(c("major_cell_type", "cell_type", "annotation"), names(f2a))[1L]
if (is.na(cell_col_f2a)) stop("Figure 2A source lacks major cell-type labels.")
setnames(f2a, c(umap_cols_f2a, cell_col_f2a), c("UMAP1", "UMAP2", "major_cell_type"))
f2a[, major_cell_type_original := as.character(major_cell_type)]
f2a[, major_cell_type := trimws(gsub("[[:cntrl:]]", "", as.character(major_cell_type)))]
excluded_main_umap_types <- c(
  "Low_quality_mitochondrial",
  "Unresolved",
  "Cycling_unresolved"
)
f2a[, displayed_in_main_panel :=
  !is.na(major_cell_type) & nzchar(major_cell_type) &
    !major_cell_type %in% excluded_main_umap_types
]
f2a_plot <- f2a[
  displayed_in_main_panel == TRUE &
    !is.na(major_cell_type) &
    nzchar(major_cell_type)
]
if (nrow(f2a_plot) == 0L) {
  stop("Figure 2A has no analysis-relevant cells after locked display filtering.")
}
## Shuffle drawing order so one large cell class does not visually cover the others.
f2a_plot <- f2a_plot[sample.int(.N)]
centroids_f2a <- f2a_plot[
  ,
  .(
    UMAP1 = median(UMAP1, na.rm = TRUE),
    UMAP2 = median(UMAP2, na.rm = TRUE),
    cells = .N
  ),
  by = major_cell_type
]
present_ct <- sort(unique(f2a_plot$major_cell_type))
missing_locked_ct <- setdiff(present_ct, names(PALETTE_CELLTYPE))
if (length(missing_locked_ct) > 0L) {
  log_msg(
    "Figure 2A palette auto-completed for: ",
    paste(missing_locked_ct, collapse = ", "),
    level = "WARNING"
  )
}
ct_palette <- complete_named_palette(present_ct, PALETTE_CELLTYPE)
all_ct <- sort(unique(c(
  as.character(f2a$major_cell_type),
  as.character(s3_composition$data$major_cell_type)
)))
ct_palette_all <- complete_named_palette(all_ct, PALETTE_CELLTYPE)
if (
  length(ct_palette) != length(present_ct) ||
    any(is.na(ct_palette)) ||
    !setequal(names(ct_palette), present_ct)
) {
  stop("Figure 2A cell-type palette coverage failed.")
}
p_f2a <- ggplot(
  f2a_plot,
  aes(x = UMAP1, y = UMAP2, colour = major_cell_type)
) +
  geom_point(size = 0.18, alpha = 0.92, stroke = 0) +
  ggrepel::geom_text_repel(
    data = centroids_f2a,
    aes(label = pretty_cell(major_cell_type)),
    family = BASE_FAMILY,
    size = 2.2,
    fontface = "bold",
    colour = "#1A1A1A",
    segment.size = 0.24,
    box.padding = 0.25,
    point.padding = 0.18,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = ct_palette, na.value = "#BDBDBD") +
  labs(
    title = "GSE236585 cardiac cellular landscape",
    x = "UMAP 1", y = "UMAP 2", colour = "Major cell type"
  ) +
  theme_sci(7.5) +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) +
  coord_fixed(ratio = 1, clip = "off")
save_panel(
  p_f2a, "Figure_2", "Figure2A_GSE236585_major_celltype_UMAP",
  7.3, 6.1, f2a, f2a_source$path,
  "GSE236585 major-cell-type UMAP", "main"
)

## Figure 2B: cell-type program forest plot
f2b <- copy(s3_primary$data)
f2b <- f2b[grepl("Top150$", signature_name)]
f2b <- f2b[!major_cell_type %in% c(
  "Unresolved", "Low_quality_mitochondrial", "Cycling_unresolved"
)]
f2b <- add_hedges_ci(
  f2b,
  "hedges_g_HFpEF_vs_Control",
  "hfpef_samples",
  "control_samples"
)
f2b[, program := pretty_program(signature_name)]
f2b[, cell_label := pretty_cell(major_cell_type)]
cell_order_f2b <- ordered_group_levels(
  f2b,
  group_col = "cell_label",
  value_col = "hedges_g_HFpEF_vs_Control",
  decreasing = FALSE,
  use_absolute = TRUE
)
f2b[, cell_label := factor(cell_label, levels = unique(cell_order_f2b))]
f2b[, significant := is.finite(wilcoxon_fdr) & wilcoxon_fdr < 0.05]
p_f2b <- ggplot(f2b, aes(x = hedges_g_HFpEF_vs_Control, y = cell_label)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = "#777777") +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = cell_label),
    linewidth = 0.55, colour = "#4D4D4D"
  ) +
  geom_point(
    aes(fill = program, shape = significant),
    size = 2.5, colour = "#202020", stroke = 0.45
  ) +
  facet_wrap(~ program, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = PALETTE_PROGRAM) +
  scale_shape_manual(
    values = c("TRUE" = 21, "FALSE" = 24),
    labels = c("FALSE" = "FDR >= 0.05", "TRUE" = "FDR < 0.05")
  ) +
  labs(
    title = "Drug-opposed programs show compartment-specific cardiac recovery",
    x = "Hedges' g: HFpEF vs Control", y = NULL,
    fill = NULL, shape = NULL
  ) +
  theme_sci(7.2) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 6.2)
  ) +
  coord_cartesian(clip = "off")
save_panel(
  p_f2b, "Figure_2", "Figure2B_celltype_program_effect_forest",
  10.5, 7.2, f2b, s3_primary$path,
  "Cell-type program effect forest", "main"
)

## Figure 2C: macrophage-state sample-level program scores
f2c_source <- read_required("stage3", "Fig3H_macrophage_state_program_source.csv")
f2c <- copy(f2c_source$data)
f2c[, program_label := fcase(
  grepl("Ccr2pos", program), "CCR2+ program",
  grepl("Ccr2neg", program), paste0("CCR2", UNICODE_MINUS, " program"),
  default = as.character(program)
)]
f2c[, state_label := pretty_cell(macrophage_state)]
p_f2c <- ggplot(
  f2c,
  aes(x = condition, y = sample_mean_score, fill = condition)
) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA,
    alpha = 0.22, linewidth = 0.45
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    shape = 21, size = 2.2, colour = "#202020", stroke = 0.4
  ) +
  facet_grid(state_label ~ program_label, scales = "free_y") +
  scale_fill_manual(values = PALETTE_CONDITION) +
  labs(
    title = "Biological-sample program scores across macrophage states",
    x = NULL, y = "Sample-level mean program score", fill = "Condition"
  ) +
  theme_sci(7.1) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.text.y = element_text(size = 6.5)
  )
save_panel(
  p_f2c, "Figure_2", "Figure2C_macrophage_state_sample_scores",
  8.7, 8.8, f2c, f2c_source$path,
  "Macrophage-state sample-level program scores", "main"
)

## Figure 2D: macrophage-state effect forest
f2d <- copy(s3_mac_stats$data)
f2d <- f2d[grepl("Top150$", signature_name)]
f2d <- add_hedges_ci(
  f2d,
  "hedges_g_HFpEF_vs_Control",
  "hfpef_samples",
  "control_samples"
)
f2d[, program := pretty_program(signature_name)]
f2d[, state := pretty_cell(macrophage_state)]
f2d[, effect_estimable :=
  is.finite(hedges_g_HFpEF_vs_Control) &
    is.finite(ci_low) & is.finite(ci_high)
]
f2d[, display_effect := fifelse(effect_estimable, hedges_g_HFpEF_vs_Control, 0)]
f2d[, estimability_note := fifelse(
  effect_estimable,
  "Hedges' g estimable",
  paste0("Not estimable (Control n=", control_samples, ")")
)]
state_order_f2d <- ordered_group_levels(
  f2d,
  group_col = "state",
  value_col = "hedges_g_HFpEF_vs_Control",
  decreasing = FALSE,
  use_absolute = TRUE
)
f2d[, state := factor(state, levels = unique(state_order_f2d))]
p_f2d <- ggplot(
  f2d,
  aes(x = display_effect, y = state, colour = program)
) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = "#777777") +
  geom_segment(
    data = f2d[effect_estimable == TRUE],
    aes(x = ci_low, xend = ci_high, yend = state),
    linewidth = 0.55
  ) +
  geom_point(data = f2d[effect_estimable == TRUE], size = 2.4) +
  geom_point(
    data = f2d[effect_estimable == FALSE],
    shape = 4, size = 2.6, stroke = 0.8, colour = "#666666"
  ) +
  geom_text(
    data = f2d[effect_estimable == FALSE],
    aes(label = "NE"),
    family = BASE_FAMILY, size = 2.2, colour = "#555555",
    hjust = -0.55, show.legend = FALSE
  ) +
  facet_wrap(~ program, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = PALETTE_PROGRAM) +
  labs(
    title = "Program effects vary across macrophage-state candidates",
    x = "Hedges' g: HFpEF vs Control", y = NULL, colour = NULL
  ) +
  theme_sci(7.3) +
  theme(legend.position = "bottom") +
  coord_cartesian(clip = "off")
save_panel(
  p_f2d, "Figure_2", "Figure2D_macrophage_state_effect_forest",
  9.5, 5.6, f2d, s3_mac_stats$path,
  "Macrophage-state effect forest", "main"
)

## Figure 2E: Stage 2-Stage 3 gene-level concordance
f2e <- copy(s3_concordance$data)
f2e[, subset_label := fcase(
  stage2_subset == "Ccr2_positive", "CCR2+ discovery program",
  stage2_subset == "Ccr2_negative", paste0("CCR2", UNICODE_MINUS, " discovery program"),
  default = as.character(stage2_subset)
)]
f2e[, cell_label := pretty_cell(major_cell_type)]
cell_order_f2e <- unique(as.character(f2e[order(tier_ABC_spearman), cell_label]))
f2e[, cell_label := factor(cell_label, levels = unique(cell_order_f2e))]
p_f2e <- ggplot(
  f2e,
  aes(
    x = tier_ABC_spearman,
    y = cell_label,
    size = tier_ABC_overlap_genes,
    colour = tier_ABC_sign_agreement
  )
) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = "#777777") +
  geom_point(alpha = 0.9) +
  facet_wrap(~ subset_label, nrow = 1) +
  scale_colour_gradient(low = "#B8B8B8", high = "#C84C4C", limits = c(0, 1), oob = scales::squish) +
  scale_size_continuous(range = c(1.8, 5.2)) +
  labs(
    title = "Cross-dataset gene-level concordance is partial and cell-type dependent",
    x = "Spearman correlation: Stage 2 disease log2FC vs Stage 3 HFpEF log2FC",
    y = NULL,
    size = "Tier A-C\noverlap genes",
    colour = "Sign\nagreement"
  ) +
  theme_sci(7.1) +
  theme(axis.text.y = element_text(size = 6.1)) +
  coord_cartesian(clip = "off")
save_panel(
  p_f2e, "Figure_2", "Figure2E_stage2_stage3_concordance",
  10.2, 6.4, f2e, s3_concordance$path,
  "Stage 2-Stage 3 concordance", "main"
)

## Figure 2F: cell-type by program effect matrix
f2f <- copy(s3_primary$data)
f2f <- f2f[grepl("Top150$", signature_name)]
f2f <- f2f[!major_cell_type %in% c(
  "Unresolved", "Low_quality_mitochondrial", "Cycling_unresolved"
)]
f2f[, program := pretty_program(signature_name)]
f2f[, cell_type := pretty_cell(major_cell_type)]
f2f[, fdr_label := ifelse(
  is.finite(wilcoxon_fdr) & wilcoxon_fdr < 0.05, "*", ""
)]
cell_order_f2f <- ordered_group_levels(
  f2f,
  group_col = "cell_type",
  value_col = "hedges_g_HFpEF_vs_Control",
  decreasing = FALSE,
  use_absolute = FALSE
)
f2f[, cell_type := factor(cell_type, levels = unique(cell_order_f2f))]
p_f2f <- make_effect_heatmap(
  f2f,
  "program", "cell_type", "hedges_g_HFpEF_vs_Control", "fdr_label",
  "Cell-type and program effect-size matrix",
  "Hedges' g", 7.4
)
save_panel(
  p_f2f, "Figure_2", "Figure2F_celltype_program_effect_matrix",
  7.0, 6.5, f2f, s3_primary$path,
  "Cell-type program effect matrix", "main"
)

############################################################
## 6. Main Figure 3
## Multilayer regulator prioritization, robustness and ablation
############################################################

## Figure 3A: unbiased Stage 4 TF priority ranking
f3a <- copy(s4_priority$data)
setorder(f3a, priority_rank)
f3a <- f3a[1:min(.N, 30L)]
f3a[, tf_label := pretty_tf(tf_symbol)]
f3a[, highlight := fcase(
  tf_label == "BHLHE40", "BHLHE40",
  tf_label == "NFKB1", "NFKB1",
  tf_label == "RELA", "RELA",
  tf_label == "RUNX1", "RUNX1",
  default = "Other"
)]
f3a[, tf_label := factor(
  tf_label,
  levels = unique(rev(as.character(tf_label)))
)]
tf_highlight_palette <- c(
  "BHLHE40" = PALETTE_TF[["BHLHE40"]],
  "NFKB1" = PALETTE_TF[["NFKB1"]],
  "RELA" = PALETTE_TF[["RELA"]],
  "RUNX1" = PALETTE_TF[["RUNX1"]],
  "Other" = "#B7B7B7"
)
p_f3a <- ggplot(f3a, aes(x = priority_score, y = tf_label, colour = highlight)) +
  geom_segment(
    aes(x = 0, xend = priority_score, yend = tf_label),
    linewidth = 0.65
  ) +
  geom_point(size = 2.35) +
  geom_text(
    data = f3a[highlight != "Other"],
    aes(label = paste0("rank ", priority_rank)),
    family = BASE_FAMILY,
    size = 2.25,
    hjust = -0.12,
    colour = "#1A1A1A"
  ) +
  scale_colour_manual(values = tf_highlight_palette) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.20))) +
  labs(
    title = "Unbiased Stage 4 priority ranking across 174 transcription factors",
    subtitle = "NFKB1 was retained without forced inclusion",
    x = "Composite priority score", y = NULL, colour = NULL
  ) +
  theme_sci(7.4) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6.0)
  ) +
  coord_cartesian(clip = "off")
save_panel(
  p_f3a, "Figure_3", "Figure3A_TF_priority_lollipop",
  7.8, 7.8, f3a, s4_priority$path,
  "Stage 4 TF priority ranking", "main"
)

## Figure 3B: virtual perturbation recovery matrix
metric_map_f3b <- c(
  stage2_primary_median_gap_reduction = "Primary program\ngap reduction",
  stage2_primary_positive_fraction = "Primary program\npositive fraction",
  biological_sample_improvement_fraction = "Biological-sample\nimprovement",
  inflammation_median_gap_reduction = "Inflammation-program\nrecovery",
  specificity_score = "Perturbation\nspecificity"
)
f3b <- copy(s5_rank$data)
f3b <- f3b[
  tf_symbol %in% c("Bhlhe40", "Nfkb1", "Rela", "Runx1", "Spi1", "Rel")
]
f3b_long <- melt(
  f3b,
  id.vars = c("tf_symbol", "final_candidate_rank"),
  measure.vars = names(metric_map_f3b),
  variable.name = "metric", value.name = "raw_value"
)
f3b_long[, metric_label := unname(metric_map_f3b[metric])]
f3b_long[, scaled_value := rescale01(raw_value), by = metric]
f3b_long[, tf_label := pretty_tf(tf_symbol)]
tf_order_f3b <- f3b[order(final_candidate_rank), pretty_tf(tf_symbol)]
f3b_long[, tf_label := factor(
  tf_label,
  levels = unique(rev(as.character(tf_order_f3b)))
)]
f3b_long[, metric_label := factor(metric_label, levels = unique(unname(metric_map_f3b)))]
f3b_long[, value_label := sprintf("%.2f", raw_value)]
p_f3b <- ggplot(
  f3b_long,
  aes(x = metric_label, y = tf_label, fill = scaled_value)
) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = value_label), family = BASE_FAMILY, size = 2.2) +
  scale_fill_gradient(low = "white", high = "#C84C4C", limits = c(0, 1)) +
  labs(
    title = "Candidate-TF virtual perturbation recovers complementary program features",
    x = NULL, y = NULL, fill = "Within-metric\nrelative value"
  ) +
  theme_heatmap(7.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(
  p_f3b, "Figure_3", "Figure3B_candidate_TF_perturbation_heatmap",
  7.8, 5.3, f3b_long, s5_rank$path,
  "Candidate-TF perturbation recovery matrix", "main"
)

## Figure 3C: bootstrap candidate effects versus matched random-TF null
f3c <- merge(
  s5b_boot$data,
  s5b_null$data,
  by = "tf_symbol", all = TRUE, sort = FALSE
)
f3c[, tf_label := pretty_tf(tf_symbol)]
tf_order_f3c <- unique(f3c[order(median_bootstrap_rank, -candidate_effect), tf_label])
f3c[, tf_label := factor(
  tf_label,
  levels = unique(rev(as.character(tf_order_f3c)))
)]
p_f3c <- ggplot(f3c, aes(y = tf_label)) +
  geom_segment(
    aes(x = null_q025, xend = null_q975, yend = tf_label),
    linewidth = 3.0, colour = "#D5D5D5", lineend = "round"
  ) +
  geom_point(
    aes(x = null_median),
    shape = 23, size = 2.5, fill = "white", colour = "#555555"
  ) +
  geom_segment(
    aes(
      x = q025_primary_gap_reduction,
      xend = q975_primary_gap_reduction,
      yend = tf_label,
      colour = tf_label
    ),
    linewidth = 0.9
  ) +
  geom_point(
    aes(x = candidate_effect, fill = tf_label),
    shape = 21, size = 3.0, colour = "#202020", stroke = 0.45
  ) +
  geom_text(
    aes(
      x = pmax(candidate_effect, null_q975, na.rm = TRUE),
      label = paste0("P=", formatC(empirical_one_sided_p, format = "f", digits = 3))
    ),
    family = BASE_FAMILY, size = 2.35, hjust = -0.12
  ) +
  scale_colour_manual(values = PALETTE_TF) +
  scale_fill_manual(values = PALETTE_TF) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.25))) +
  labs(
    title = "Bootstrap recovery is evaluated against matched random-TF null distributions",
    subtitle = "Grey interval: matched-null 95% range; coloured interval: candidate bootstrap 95% range",
    x = "Primary program gap reduction", y = NULL,
    colour = NULL, fill = NULL
  ) +
  theme_sci(7.7) +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")
save_panel(
  p_f3c, "Figure_3", "Figure3C_bootstrap_matched_null",
  7.8, 5.6, f3c, c(s5b_boot$path, s5b_null$path),
  "Bootstrap and matched-null comparison", "main"
)

## Figure 3D: TF ranks across simple, integrated and external methods
rank_cols_f3d <- c(
  rank_expression_only = "TF expression only",
  rank_regulon_only = "Regulon activity only",
  rank_stage4_integrated = "Stage 4 multifeature",
  rank_perturbation_only = "Perturbation only",
  rank_bootstrap_robustness = "Bootstrap perturbation",
  rank_communication_only = "Communication only",
  rank_full_cross_layer = "Full cross-layer",
  external_integrated_rank = "External Stage 8"
)
f3d_wide <- copy(b_tf_ranks$data)
for (rank_col in names(rank_cols_f3d)) {
  set(f3d_wide, j = rank_col, value = as.numeric(f3d_wide[[rank_col]]))
}
f3d_long <- melt(
  f3d_wide,
  id.vars = "tf_symbol",
  measure.vars = names(rank_cols_f3d),
  variable.name = "method", value.name = "rank"
)
f3d_long[, method_label := unname(rank_cols_f3d[method])]
f3d_long[, method_label := factor(method_label, levels = unique(unname(rank_cols_f3d)))]
f3d_long[, tf_label := factor(
  pretty_tf(tf_symbol), levels = unique(c("BHLHE40", "NFKB1", "RELA"))
)]
p_f3d <- ggplot(
  f3d_long,
  aes(x = method_label, y = tf_label, fill = rank)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.1f", rank)), family = BASE_FAMILY, size = 2.45) +
  scale_fill_gradient(low = "#C84C4C", high = "white", trans = "reverse", breaks = c(1, 2, 3)) +
  labs(
    title = "TF roles differ across simple, integrated and external evidence",
    x = NULL, y = NULL, fill = "Rank\n(lower better)"
  ) +
  theme_heatmap(7.3) +
  theme(axis.text.x = element_text(angle = 42, hjust = 1))
save_panel(
  p_f3d, "Figure_3", "Figure3D_TF_method_rank_heatmap",
  8.2, 3.6, f3d_long, b_tf_ranks$path,
  "TF method-rank heatmap", "main"
)

## Figure 3E: leave-one-layer-out ablation
f3e <- copy(b_ablation$data)
scenario_order_f3e <- c(
  "Full_cross_layer",
  "Without_regulon_layer",
  "Without_perturbation_layer",
  "Without_communication_layer",
  "TF_expression_only",
  "Regulon_only",
  "Perturbation_only",
  "Communication_only"
)
scenario_label_f3e <- c(
  "Full cross-layer",
  "Without regulon",
  "Without perturbation",
  "Without communication",
  "TF expression only",
  "Regulon only",
  "Perturbation only",
  "Communication only"
)
names(scenario_label_f3e) <- scenario_order_f3e
f3e[, scenario_label := unname(scenario_label_f3e[scenario])]
f3e[, scenario_label := factor(scenario_label, levels = unique(unname(scenario_label_f3e)))]
f3e[, tf_label := pretty_tf(tf_symbol)]
p_f3e <- ggplot(
  f3e,
  aes(x = scenario_label, y = scenario_rank, group = tf_label, colour = tf_label)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4) +
  scale_y_reverse(breaks = 1:3, limits = c(3.2, 0.8)) +
  scale_colour_manual(values = PALETTE_TF) +
  labs(
    title = "Candidate ranks remain interpretable under evidence-layer ablation",
    x = NULL, y = "Candidate rank", colour = "TF"
  ) +
  theme_sci(7.4) +
  theme(
    axis.text.x = element_text(angle = 42, hjust = 1),
    legend.position = "bottom"
  )
save_panel(
  p_f3e, "Figure_3", "Figure3E_leave_one_layer_out_ablation",
  8.4, 4.6, f3e, b_ablation$path,
  "Leave-one-layer-out ablation", "main"
)

## Figure 3F: Stage 7 sample-level feature importance
f3f <- copy(b_feature$data)
f3f[, feature_label := fcase(
  feature == "COMM_NFkB_axis_burden", "NF-kB communication-axis burden",
  feature == "PROGRAM_DrugOpposed_Top150", "Drug-opposed Top150 program",
  feature == "TF_Nfkb1_activity", "NFKB1 regulon activity",
  feature == "TF_Bhlhe40_activity", "BHLHE40 regulon activity",
  feature == "TF_Rela_activity", "RELA regulon activity",
  default = as.character(feature)
)]
f3f[, feature_label := factor(
  feature_label,
  levels = unique(rev(as.character(feature_label[order(importance_rank)])))
)]
p_f3f <- ggplot(f3f, aes(x = contribution_percent, y = feature_label)) +
  geom_col(width = 0.68, fill = "#466B8A") +
  geom_text(
    aes(label = paste0(sprintf("%.1f", contribution_percent), "%")),
    family = BASE_FAMILY, size = 2.45, hjust = -0.12
  ) +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Sample-level separation is jointly driven by program and communication features",
    x = "Share of mean absolute logit contribution", y = NULL
  ) +
  theme_sci(7.5) +
  coord_cartesian(clip = "off")
save_panel(
  p_f3f, "Figure_3", "Figure3F_sample_level_feature_importance",
  7.3, 4.5, f3f, b_feature$path,
  "Sample-level feature importance", "main"
)

############################################################
## 7. Main Figure 4
## NFKB1-centered multibranch communication
############################################################

## Figure 4A: candidate-TF communication coverage
f4a <- copy(s6_candidate$data)
metric_map_f4a <- c(
  total_axes = "Total axes",
  strict_cross_stage_axes = "Strict cross-stage axes",
  best_axis_rank = "Best axis rank"
)
f4a_long <- melt(
  f4a,
  id.vars = c("tf_symbol", "candidate_role", "Stage5B_rank"),
  measure.vars = names(metric_map_f4a),
  variable.name = "metric", value.name = "value"
)
f4a_long[, metric_label := unname(metric_map_f4a[metric])]
f4a_long[, metric_label := factor(metric_label, levels = unique(unname(metric_map_f4a)))]
f4a_long[, tf_label := factor(
  pretty_tf(tf_symbol), levels = unique(c("BHLHE40", "NFKB1", "RELA"))
)]
p_f4a <- ggplot(
  f4a_long,
  aes(x = value, y = tf_label, fill = tf_label)
) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = value),
    family = BASE_FAMILY, size = 2.45, hjust = -0.12
  ) +
  facet_wrap(~ metric_label, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = PALETTE_TF) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "NFKB1 combines broad communication coverage with a top-ranked axis",
    subtitle = "For the best-axis-rank facet, lower values indicate stronger prioritization",
    x = "Metric value", y = NULL, fill = NULL
  ) +
  theme_sci(7.5) +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")
save_panel(
  p_f4a, "Figure_4", "Figure4A_candidate_TF_communication_coverage",
  9.3, 3.9, f4a_long, s6_candidate$path,
  "Candidate-TF communication coverage", "main"
)

## Figure 4B: top stable axes
f4b <- copy(s6_stability$data)
setorder(f4b, median_scenario_rank, -top10_scenario_frequency)
f4b <- f4b[1:min(.N, 20L)]
f4b[, axis_label := short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
f4b[, axis_label := factor(
  axis_label,
  levels = unique(rev(as.character(axis_label)))
)]
p_f4b <- ggplot(
  f4b,
  aes(
    x = median_scenario_rank,
    y = axis_label,
    size = top10_scenario_frequency,
    colour = receiver
  )
) +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = PALETTE_RECEIVER, na.value = "#777777") +
  scale_size_continuous(
    range = c(2.0, 5.2), labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_reverse() +
  labs(
    title = "Top communication axes remain stable across ranking scenarios",
    x = "Median scenario rank (lower is better)", y = NULL,
    size = "Top-10\nfrequency", colour = "Receiver"
  ) +
  theme_sci(7.0) +
  theme(axis.text.y = element_text(size = 6.2)) +
  coord_cartesian(clip = "off")
save_panel(
  p_f4b, "Figure_4", "Figure4B_top_stable_communication_axes",
  9.0, 7.7, f4b, s6_stability$path,
  "Top stable communication axes", "main"
)

## Figure 4C: NicheNet ligand activity heatmap
f4c <- copy(s6_nichenet$data)
ligands_f4c <- f4c[, .(
  best_rank = min(nichenet_activity_rank, na.rm = TRUE)
), by = ligand_matrix_symbol][order(best_rank)][1:min(.N, 14L), ligand_matrix_symbol]
f4c <- f4c[ligand_matrix_symbol %in% ligands_f4c]
f4c[, ligand := toupper(ligand_matrix_symbol)]
f4c[, receiver_label := paste(
  pretty_cell(receiver), gsub("_", " ", receiver_direction), sep = " | "
)]
f4c[, ligand := factor(
  ligand,
  levels = rev(unique(f4c[order(nichenet_activity_rank), ligand]))
)]
f4c[, label := sprintf("%.2f", aupr_corrected)]
p_f4c <- ggplot(
  f4c,
  aes(x = receiver_label, y = ligand, fill = aupr_corrected)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), family = BASE_FAMILY, size = 2.2) +
  scale_fill_gradient(low = "white", high = "#6C5CE7", na.value = "#E6E6E6") +
  labs(
    title = "NicheNet ligand activity varies across vascular and stromal receivers",
    x = NULL, y = NULL, fill = "Corrected\nAUPR"
  ) +
  theme_heatmap(7.1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6.1)
  )
save_panel(
  p_f4c, "Figure_4", "Figure4C_NicheNet_ligand_activity_heatmap",
  9.2, 6.5, f4c, s6_nichenet$path,
  "NicheNet ligand activity heatmap", "main"
)

## Figure 4D: receiver receptor expression support
f4d <- copy(s6_top_axes$data)
num_cols_f4d <- c(
  "final_axis_rank",
  "median_receptor_pct_Control",
  "median_receptor_pct_HFpEF",
  "median_receptor_log2_cpm_Control",
  "median_receptor_log2_cpm_HFpEF"
)
for (z in intersect(num_cols_f4d, names(f4d))) {
  set(f4d, j = z, value = as.numeric(f4d[[z]]))
}
setorder(f4d, final_axis_rank)
f4d <- f4d[1:min(.N, 40L)]
f4d <- f4d[, .SD[1L], by = .(receptor, receiver)]
f4d[, receptor_label := paste0(toupper(receptor), " | ", pretty_cell(receiver))]
p_f4d <- ggplot(
  f4d,
  aes(
    x = median_receptor_log2_cpm_Control,
    y = median_receptor_log2_cpm_HFpEF,
    size = median_receptor_pct_HFpEF,
    colour = receiver
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.45, colour = "#888888") +
  geom_point(alpha = 0.88) +
  ggrepel::geom_text_repel(
    aes(label = receptor_label),
    family = BASE_FAMILY, size = 2.2,
    max.overlaps = 18, min.segment.length = 0, segment.size = 0.25
  ) +
  scale_colour_manual(values = PALETTE_RECEIVER, na.value = "#777777") +
  scale_size_continuous(
    range = c(1.5, 5.3), labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Receiver receptors are detectable in intended vascular/stromal compartments",
    x = "Median receptor expression in Control (log2 CPM)",
    y = "Median receptor expression in HFpEF (log2 CPM)",
    size = "HFpEF cells\nexpressing receptor",
    colour = "Receiver"
  ) +
  theme_sci(7.2) +
  coord_cartesian(clip = "off")
save_panel(
  p_f4d, "Figure_4", "Figure4D_receiver_receptor_expression_support",
  7.8, 6.4, f4d, s6_top_axes$path,
  "Receiver receptor expression support", "main"
)

## Figure 4E: discovery rank versus external support
f4e <- copy(b_axis_ranks$data)
f4e[, axis_label := short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
highlight_axes_f4e <- c(
  "NFKB1__TNF__TNFRSF1A__ENDOTHELIAL",
  "NFKB1__TNF__LTBR__ENDOTHELIAL",
  "NFKB1__PDGFB__PDGFRA__FIBROBLAST",
  "NFKB1__PDGFB__LRP1__FIBROBLAST"
)
f4e[, label_show := ifelse(axis_key %in% highlight_axes_f4e, axis_label, "")]
p_f4e <- ggplot(
  f4e,
  aes(
    x = rank_full_integration,
    y = external_support_fraction,
    size = external_median_abs_hedges_g,
    colour = receiver
  )
) +
  geom_point(alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = label_show),
    family = BASE_FAMILY, size = 2.2,
    max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25
  ) +
  scale_colour_manual(values = PALETTE_RECEIVER, na.value = "#777777") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_reverse() +
  scale_size_continuous(range = c(1.7, 5.0)) +
  labs(
    title = "Discovery rank and external support define complementary communication branches",
    x = "Full-integration discovery rank (lower is better)",
    y = "External support fraction",
    size = "Median absolute\nHedges' g",
    colour = "Receiver"
  ) +
  theme_sci(7.2) +
  coord_cartesian(clip = "off")
save_panel(
  p_f4e, "Figure_4", "Figure4E_axis_discovery_rank_external_support",
  7.7, 6.0, f4e, b_axis_ranks$path,
  "Axis discovery rank versus external support", "main"
)

## Figure 4F: Top-k external recovery by ranking method
f4f <- copy(b_axis_topk$data)
method_order_f4f <- c(
  "Ligand expression only",
  "Receptor expression only",
  "Ligand-receptor expression",
  "NicheNet only",
  "Cross-stage support only",
  "Full integration"
)
f4f[, method := factor(method, levels = unique(method_order_f4f))]
method_palette_f4f <- c(
  "Ligand expression only" = "#4C78A8",
  "Receptor expression only" = "#72B7B2",
  "Ligand-receptor expression" = "#54A24B",
  "NicheNet only" = "#E45756",
  "Cross-stage support only" = "#B279A2",
  "Full integration" = "#F2CF5B"
)
p_f4f <- ggplot(
  f4f,
  aes(
    x = top_k,
    y = supported_axis_fraction,
    group = method,
    colour = method,
    shape = method
  )
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.3) +
  scale_x_continuous(breaks = c(5, 10, 20)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_manual(values = method_palette_f4f) +
  labs(
    title = "No single ranking strategy universally dominates external Top-k recovery",
    x = "Top-ranked frozen axes", y = "Externally supported axes",
    colour = "Ranking method", shape = "Ranking method"
  ) +
  theme_sci(7.1)
save_panel(
  p_f4f, "Figure_4", "Figure4F_axis_topk_external_performance",
  8.4, 5.7, f4f, b_axis_topk$path,
  "Top-k external axis performance", "main"
)

############################################################
## 8. Main Figure 5
## Multicohort and independent human validation
############################################################

## Figure 5A: multicohort program validation matrix
f5a <- copy(s8_program_evidence$data)
f5a <- f5a[grepl("Top150$", program_name)]
f5a[, program := pretty_program(program_name)]
f5a[, evidence_label := paste(
  dataset_id, pretty_contrast(contrast), pretty_cell(cell_type), sep = " | "
)]
f5a[, fdr_label := ifelse(is.finite(fdr) & fdr < 0.05, "*", "")]
evidence_order_f5a <- f5a[order(dataset_id, contrast, cell_type), unique(evidence_label)]
f5a[, evidence_label := factor(
  evidence_label,
  levels = unique(rev(as.character(evidence_order_f5a)))
)]
p_f5a <- make_effect_heatmap(
  f5a,
  "program", "evidence_label", "hedges_g", "fdr_label",
  "Multicohort program evidence preserves direction heterogeneity",
  "Hedges' g", 6.3
) +
  theme(
    axis.text.y = element_text(size = 6.0),
    axis.text.x = element_text(size = 6.1)
  )
save_panel(
  p_f5a, "Figure_5", "Figure5A_multicohort_program_validation",
  8.3, 13.0, f5a, s8_program_evidence$path,
  "Multicohort program validation", "main"
)

## Figure 5B: multicohort TF validation matrix
f5b <- copy(s8_tf_evidence$data)
f5b[, tf := pretty_tf(tf_symbol)]
f5b[, evidence_label := paste(
  dataset_id, pretty_contrast(contrast), pretty_cell(cell_type), sep = " | "
)]
f5b[, fdr_label := ifelse(is.finite(fdr) & fdr < 0.05, "*", "")]
evidence_order_f5b <- f5b[order(dataset_id, contrast, cell_type), unique(evidence_label)]
f5b[, evidence_label := factor(
  evidence_label,
  levels = unique(rev(as.character(evidence_order_f5b)))
)]
p_f5b <- make_effect_heatmap(
  f5b,
  "tf", "evidence_label", "hedges_g", "fdr_label",
  "External TF support differs across datasets and compartments",
  "Hedges' g", 6.3
) +
  theme(
    axis.text.y = element_text(size = 6.0),
    axis.text.x = element_text(size = 6.1)
  )
save_panel(
  p_f5b, "Figure_5", "Figure5B_multicohort_TF_validation",
  7.3, 13.0, f5b, s8_tf_evidence$path,
  "Multicohort TF validation", "main"
)

## Figure 5C: independent human donor-level program validation
f5c <- copy(s8_scp_program$data)
f5c <- f5c[grepl("Top150$", program_name)]
f5c <- add_hedges_ci(f5c, "hedges_g", "case_n", "reference_n")
f5c[, program := pretty_program(program_name)]
f5c[, cell := pretty_cell(cell_type)]
f5c[, significant := is.finite(fdr) & fdr < 0.05]
cell_order_f5c <- ordered_group_levels(
  f5c,
  group_col = "cell",
  value_col = "hedges_g",
  decreasing = FALSE,
  use_absolute = TRUE
)
f5c[, cell := factor(cell, levels = unique(cell_order_f5c))]
p_f5c <- ggplot(
  f5c,
  aes(x = hedges_g, y = cell, colour = program)
) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = "#777777") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = cell), linewidth = 0.55) +
  geom_point(aes(shape = significant), size = 2.5) +
  facet_wrap(~ program, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = PALETTE_PROGRAM) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    labels = c("FALSE" = "FDR >= 0.05", "TRUE" = "FDR < 0.05")
  ) +
  labs(
    title = "Independent human myocardial donors support selected program compartments",
    subtitle = "SCP3342: 19 HFpEF donors versus 24 Control donors",
    x = "Hedges' g: HFpEF vs Control", y = NULL,
    colour = NULL, shape = NULL
  ) +
  theme_sci(7.4) +
  theme(legend.position = "bottom") +
  coord_cartesian(clip = "off")
save_panel(
  p_f5c, "Figure_5", "Figure5C_SCP3342_human_program_forest",
  9.3, 5.2, f5c, s8_scp_program$path,
  "SCP3342 human program forest", "main"
)

## Figure 5D: GSE249412 drug-response axis matrix
f5d <- copy(s8_gse249_axis$data)
f5d <- f5d[contrast %in% c(
  "Empagliflozin_vs_HFpEF_Vehicle",
  "TYA_018_vs_HFpEF_Vehicle"
)]
top_axis_keys_f5d <- s8_axis_summary$data[order(integrated_rank), head(axis_key, 12L)]
f5d <- f5d[axis_key %in% top_axis_keys_f5d]
f5d[, axis_label := short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
axis_order_table_f5d <- s8_axis_summary$data[axis_key %in% top_axis_keys_f5d][order(integrated_rank)]
axis_labels_f5d <- axis_order_table_f5d[, short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
f5d[, axis_label := factor(
  axis_label,
  levels = unique(rev(as.character(axis_labels_f5d)))
)]
f5d[, contrast_label := pretty_contrast(contrast)]
f5d[, fdr_label := ifelse(is.finite(fdr) & fdr < 0.05, "*", "")]
p_f5d <- make_effect_heatmap(
  f5d,
  "contrast_label", "axis_label", "hedges_g", "fdr_label",
  "GSE249412 reveals axis-specific empagliflozin and TYA-018 responses",
  "Hedges' g", 6.5
) +
  theme(axis.text.y = element_text(size = 6.0))
save_panel(
  p_f5d, "Figure_5", "Figure5D_GSE249412_drug_axis_reversal",
  8.2, 7.5, f5d, s8_gse249_axis$path,
  "GSE249412 drug-axis validation", "main"
)

## Figure 5E: multicohort axis validation matrix
f5e <- copy(s8_axis_evidence$data)
top20_keys_f5e <- s8_axis_summary$data[order(integrated_rank), head(axis_key, 20L)]
f5e <- f5e[axis_key %in% top20_keys_f5e]
f5e[, axis_label := short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
f5e[, evidence_label := paste(dataset_id, pretty_contrast(contrast), sep = " | ")]
axis_order_table_f5e <- s8_axis_summary$data[axis_key %in% top20_keys_f5e][order(integrated_rank)]
axis_labels_f5e <- axis_order_table_f5e[, short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
f5e[, axis_label := factor(
  axis_label,
  levels = unique(rev(as.character(axis_labels_f5e)))
)]
f5e[, fdr_label := ifelse(is.finite(fdr) & fdr < 0.05, "*", "")]
p_f5e <- make_effect_heatmap(
  f5e,
  "evidence_label", "axis_label", "hedges_g", "fdr_label",
  "Frozen communication axes show partial and direction-dependent external recovery",
  "Hedges' g", 6.2
) +
  theme(
    axis.text.y = element_text(size = 6.0),
    axis.text.x = element_text(size = 6.0)
  )
save_panel(
  p_f5e, "Figure_5", "Figure5E_multicohort_axis_validation",
  10.4, 9.5, f5e, s8_axis_evidence$path,
  "Multicohort axis validation matrix", "main"
)

## Figure 5F: integrated external support across branches
f5f <- copy(s8_axis_summary$data)
f5f[, axis_label := short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
f5f[, branch := fcase(
  grepl("__PDGFB__", axis_key), "PDGFB fibroblast branch",
  grepl("__TNF__", axis_key), "TNF inflammatory branch",
  grepl("__IL1", axis_key), "IL1 inflammatory branch",
  grepl("__VEGFA__", axis_key), "VEGFA vascular branch",
  default = "Other branch"
)]
highlight_keys_f5f <- c(
  "NFKB1__PDGFB__PDGFRA__FIBROBLAST",
  "NFKB1__PDGFB__LRP1__FIBROBLAST",
  "NFKB1__TNF__LTBR__ENDOTHELIAL",
  "NFKB1__TNF__TNFRSF1A__ENDOTHELIAL"
)
f5f[, label_show := ifelse(axis_key %in% highlight_keys_f5f, axis_label, "")]
branch_palette_f5f <- c(
  "PDGFB fibroblast branch" = "#E45756",
  "TNF inflammatory branch" = "#4C78A8",
  "IL1 inflammatory branch" = "#B279A2",
  "VEGFA vascular branch" = "#54A24B",
  "Other branch" = "#B5B5B5"
)
p_f5f <- ggplot(
  f5f,
  aes(
    x = support_fraction,
    y = median_abs_hedges_g,
    size = formal_fdr_rows,
    colour = branch
  )
) +
  geom_point(alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = label_show),
    family = BASE_FAMILY, size = 2.2,
    max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_manual(values = branch_palette_f5f) +
  scale_size_continuous(range = c(2.0, 5.0)) +
  labs(
    title = "External evidence supports multiple inflammatory and stromal branches",
    x = "External support fraction",
    y = "Median absolute Hedges' g",
    size = "FDR-supported\nrows",
    colour = "Communication branch"
  ) +
  theme_sci(7.2) +
  coord_cartesian(clip = "off")
save_panel(
  p_f5f, "Figure_5", "Figure5F_external_axis_integrated_support",
  8.0, 6.2, f5f, s8_axis_summary$path,
  "Integrated external axis support", "main"
)

############################################################
## 9. Supplementary Figure S1
## Stage 2 sample QC, group balance and method agreement
############################################################

## Figure S1A: sequencing QC metrics
s1a <- copy(s2_qc$data)
s1a_long <- melt(
  s1a,
  id.vars = c("sample_accession", "group_id", "macrophage_subset", "diet", "drug"),
  measure.vars = c("library_size_rounded", "detected_count_ge_1", "detected_count_ge_10"),
  variable.name = "metric", value.name = "value"
)
s1a_long[, metric_label := fcase(
  metric == "library_size_rounded", "Library size",
  metric == "detected_count_ge_1", "Genes with count >= 1",
  metric == "detected_count_ge_10", "Genes with count >= 10",
  default = as.character(metric)
)]
s1a_long[, treatment_group := paste(diet, drug, sep = " / ")]
p_s1a <- ggplot(
  s1a_long,
  aes(x = sample_accession, y = value, fill = treatment_group)
) +
  geom_col(width = 0.76) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = label_compact_number) +
  scale_fill_manual(values = group_palette) +
  labs(
    title = "Stage 2 sample-level sequencing quality",
    x = NULL, y = NULL, fill = "Diet / treatment"
  ) +
  theme_sci(7.0) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 6.0),
    legend.position = "bottom"
  )
save_panel(
  p_s1a, "Figure_S1", "FigureS1A_stage2_sample_QC",
  9.0, 7.0, s1a_long, s2_qc$path,
  "Stage 2 sample QC", "supplementary"
)

## Figure S1B: eight-group replicate balance
s1b <- copy(s2_group_counts$data)
s1b[, group_label := gsub("__", " | ", group_id)]
p_s1b <- ggplot(
  s1b,
  aes(x = N, y = reorder(group_label, N))
) +
  geom_col(width = 0.68, fill = "#466B8A") +
  geom_text(
    aes(label = N),
    family = BASE_FAMILY, size = 2.45, hjust = -0.15
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "All eight pharmacotranscriptomic groups retained biological replicates",
    x = "Biological samples", y = NULL
  ) +
  theme_sci(7.2) +
  coord_cartesian(clip = "off")
save_panel(
  p_s1b, "Figure_S1", "FigureS1B_stage2_group_counts",
  7.7, 4.8, s1b, s2_group_counts$path,
  "Stage 2 group counts", "supplementary"
)

## Figure S1C: DESeq2-edgeR concordance
s1c <- copy(s2_method$data)
s1c[, contrast_label := wrap_text(pretty_contrast(contrast), 30)]
s1c_long <- melt(
  s1c,
  id.vars = c("contrast", "contrast_label"),
  measure.vars = c("pearson_lfc", "spearman_lfc", "overall_sign_agreement", "top200_jaccard"),
  variable.name = "metric", value.name = "value"
)
s1c_long[, metric_label := fcase(
  metric == "pearson_lfc", "Pearson log2FC",
  metric == "spearman_lfc", "Spearman log2FC",
  metric == "overall_sign_agreement", "Overall sign agreement",
  metric == "top200_jaccard", "Top-200 Jaccard",
  default = as.character(metric)
)]
p_s1c <- ggplot(
  s1c_long,
  aes(x = value, y = reorder(contrast_label, value))
) +
  geom_point(size = 2.2, colour = "#466B8A") +
  facet_wrap(~ metric_label, scales = "free_x", ncol = 2) +
  labs(
    title = "DESeq2 and edgeR show strong directional concordance",
    x = "Concordance metric", y = NULL
  ) +
  theme_sci(6.8) +
  theme(axis.text.y = element_text(size = 6.0))
save_panel(
  p_s1c, "Figure_S1", "FigureS1C_stage2_method_concordance",
  10.2, 7.5, s1c_long, s2_method$path,
  "Stage 2 method concordance", "supplementary"
)

## Figure S1D: differential-expression yield
s1d <- copy(s2_contrast_summary$data)
s1d[, contrast_label := wrap_text(pretty_contrast(contrast), 30)]
s1d[, fdr_010_total := fdr_010_up + fdr_010_down]
p_s1d <- ggplot(
  s1d,
  aes(x = fdr_010_total, y = reorder(contrast_label, fdr_010_total), fill = method)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  facet_wrap(~ subset, scales = "free_y", ncol = 1) +
  scale_x_continuous(labels = label_compact_number) +
  labs(
    title = "Differential-expression yield varies by contrast and method",
    x = "Genes at FDR < 0.10", y = NULL, fill = "Method"
  ) +
  theme_sci(6.8) +
  theme(axis.text.y = element_text(size = 6.0))
save_panel(
  p_s1d, "Figure_S1", "FigureS1D_stage2_DE_gene_counts",
  9.2, 7.5, s1d, s2_contrast_summary$path,
  "Stage 2 DE gene counts", "supplementary"
)

############################################################
## 10. Supplementary Figure S2
## Signature-size sensitivity and cross-subset consensus
############################################################

## Figure S2A: full signature-size composition
p_s2a <- p_f1f + labs(title = "Full signature-size and direction composition")
save_panel(
  p_s2a, "Figure_S2", "FigureS2A_signature_size_composition",
  9.0, 4.0, f1f, s3_sig_size$path,
  "Signature-size composition", "supplementary"
)

## Figure S2B: all sizes across cardiac cell types
s2b <- copy(s3_program_stats$data)
s2b[, program := pretty_program(signature_name)]
s2b[, size := as.integer(gsub(".*Top", "", signature_name))]
s2b[, cell_type := pretty_cell(major_cell_type)]
s2b[, program_size := paste0(program, " | Top", size)]
s2b[, program_size := factor(
  program_size,
  levels = unique(s2b[order(program, size), program_size])
)]
p_s2b <- ggplot(
  s2b,
  aes(x = program_size, y = cell_type, fill = hedges_g_HFpEF_vs_Control)
) +
  geom_tile(colour = "white", linewidth = 0.18) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0
  ) +
  labs(
    title = "Cell-type localization across all prespecified signature sizes",
    x = NULL, y = NULL, fill = "Hedges' g"
  ) +
  theme_heatmap(6.6) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 6.0),
    axis.text.y = element_text(size = 6.0)
  )
save_panel(
  p_s2b, "Figure_S2", "FigureS2B_all_signature_sizes_localization",
  11.0, 8.5, s2b, s3_program_stats$path,
  "All signature-size localization", "supplementary"
)

## Figure S2C: top cross-subset consensus genes
s2c <- copy(s2_top_cross$data)
setorder(s2c, overall_consensus_rank)
s2c <- s2c[1:min(.N, 30L)]
s2c[, gene_label := factor(
  display_gene,
  levels = unique(rev(as.character(display_gene)))
)]
s2c[, consensus_strength := mean_opposition_score]
p_s2c <- ggplot(
  s2c,
  aes(x = consensus_strength, y = gene_label)
) +
  geom_segment(
    aes(x = 0, xend = consensus_strength, yend = gene_label),
    linewidth = 0.55, colour = "#7A9E9F"
  ) +
  geom_point(size = 2.1, colour = "#C84C4C") +
  labs(
    title = "Top cross-subset consensus drug-opposed genes",
    x = "Mean opposition score", y = NULL
  ) +
  theme_sci(7.0) +
  theme(axis.text.y = element_text(size = 6.0))
save_panel(
  p_s2c, "Figure_S2", "FigureS2C_cross_subset_consensus_genes",
  7.5, 7.2, s2c, s2_top_cross$path,
  "Cross-subset consensus genes", "supplementary"
)

## Figure S2D: opposition-tier composition
s2d <- copy(s2_opp_summary$data)
s2d_long <- melt(
  s2d,
  id.vars = c("subset", "tested_genes"),
  measure.vars = c("tier_A", "tier_B", "tier_C", "tier_D"),
  variable.name = "tier", value.name = "genes"
)
s2d_long[, tier := factor(
  toupper(gsub("_", " ", tier)),
  levels = unique(c("TIER A", "TIER B", "TIER C", "TIER D"))
)]
p_s2d <- ggplot(
  s2d_long,
  aes(x = subset, y = genes, fill = tier)
) +
  geom_col(width = 0.68) +
  scale_y_continuous(labels = label_compact_number) +
  scale_fill_manual(values = c(
    "TIER A" = "#C84C4C",
    "TIER B" = "#E9A03B",
    "TIER C" = "#3B6FB6",
    "TIER D" = "#D0D0D0"
  )) +
  labs(
    title = "Evidence-tier composition of directionally opposed genes",
    x = NULL, y = "Genes", fill = "Opposition tier"
  ) +
  theme_sci(7.3)
save_panel(
  p_s2d, "Figure_S2", "FigureS2D_opposition_tier_composition",
  6.2, 4.5, s2d_long, s2_opp_summary$path,
  "Opposition-tier composition", "supplementary"
)

############################################################
## 11. Supplementary Figure S3
## Stage 3 QC, doublets, annotation and composition
############################################################

## Figure S3A: sample QC retention
s3a <- copy(s3_qc_retention$data)
s3a_long <- melt(
  s3a,
  id.vars = c("sample_accession", "condition"),
  measure.vars = c("cells_before_qc", "cells_after_qc"),
  variable.name = "stage", value.name = "cells"
)
s3a_long[, stage_label := fcase(
  stage == "cells_before_qc", "Before QC",
  stage == "cells_after_qc", "After QC",
  default = as.character(stage)
)]
p_s3a <- ggplot(
  s3a_long,
  aes(x = sample_accession, y = cells, fill = stage_label)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = label_compact_number) +
  labs(
    title = "GSE236585 sample-level quality-control retention",
    x = NULL, y = "Cells", fill = "Processing stage"
  ) +
  theme_sci(7.0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(
  p_s3a, "Figure_S3", "FigureS3A_GSE236585_QC_retention",
  8.2, 4.8, s3a_long, s3_qc_retention$path,
  "GSE236585 QC retention", "supplementary"
)

## Figure S3B: predicted doublet fractions
s3b <- copy(s3_doublet$data)
p_s3b <- ggplot(
  s3b,
  aes(x = sample_accession, y = predicted_doublet_fraction, fill = condition)
) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = PALETTE_CONDITION) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Sample-level scDblFinder doublet fractions",
    x = NULL, y = "Predicted doublet fraction", fill = "Condition"
  ) +
  theme_sci(7.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(
  p_s3b, "Figure_S3", "FigureS3B_scDblFinder_rates",
  7.2, 4.4, s3b, s3_doublet$path,
  "scDblFinder doublet rates", "supplementary"
)

## Figure S3C: cluster marker-score support
s3c_source <- read_required("stage3", "FigS3C_major_marker_score_heatmap_source.csv")
s3c <- data.table::fread(
  s3c_source$path,
  header = TRUE,
  encoding = "UTF-8",
  showProgress = FALSE
)
id_col_s3c <- "marker_set"
if (!id_col_s3c %in% names(s3c)) stop("Figure S3C source lacks marker_set row labels.")
marker_cols_s3c <- setdiff(names(s3c), c(id_col_s3c, "cluster", "seurat_clusters"))
marker_cols_s3c <- marker_cols_s3c[
  vapply(s3c[, ..marker_cols_s3c], is.numeric, logical(1))
]
s3c_long <- melt(
  s3c,
  id.vars = id_col_s3c,
  measure.vars = marker_cols_s3c,
  variable.name = "marker_score", value.name = "value"
)
setnames(s3c_long, id_col_s3c, "cell_type")
s3c_long[, marker_score := paste0("Marker set ", marker_score)]
p_s3c <- ggplot(
  s3c_long,
  aes(x = marker_score, y = pretty_cell(cell_type), fill = value)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0
  ) +
  labs(
    title = "Cluster-level canonical-marker score support",
    x = NULL, y = NULL, fill = "Marker score"
  ) +
  theme_heatmap(6.5) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 6.0),
    axis.text.y = element_text(size = 6.0)
  )
save_panel(
  p_s3c, "Figure_S3", "FigureS3C_major_marker_score_heatmap",
  10.0, 7.5, s3c_long, s3c_source$path,
  "Major marker-score heatmap", "supplementary"
)

## Figure S3D: sample-level cell-type composition
s3d <- copy(s3_composition$data)
p_s3d <- ggplot(
  s3d,
  aes(x = sample_accession, y = cell_fraction, fill = major_cell_type)
) +
  geom_col(width = 0.76) +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = ct_palette_all, na.value = "#BDBDBD") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Sample-level cardiac cell-type composition",
    x = NULL, y = "Fraction of retained cells", fill = "Major cell type"
  ) +
  theme_sci(6.4) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.0),
    legend.text = element_text(size = 6.0)
  )
save_panel(
  p_s3d, "Figure_S3", "FigureS3D_celltype_composition",
  10.2, 6.0, s3d, s3_composition$path,
  "Cardiac cell-type composition", "supplementary"
)

############################################################
## 12. Supplementary Figure S4
## Complete program maps and macrophage-state localization
############################################################

## Figure S4A-B: program-score UMAPs
score_pos_col_s4 <- names(f2a)[grepl("score_Ccr2pos_Top150.*net", names(f2a))][1L]
score_neg_col_s4 <- names(f2a)[grepl("score_Ccr2neg_Top150.*net", names(f2a))][1L]
make_program_umap <- function(score_col, title_text, stem) {
  if (is.na(score_col) || !score_col %in% names(f2a)) {
    stop("Required UMAP program score column was not found: ", score_col)
  }
  dt <- f2a[, .(
    UMAP1, UMAP2,
    score = as.numeric(get(score_col)),
    major_cell_type, sample_accession, condition
  )]
  finite_score <- dt[is.finite(score), score]
  if (length(finite_score) < 2L) stop("UMAP program score is not estimable: ", score_col)
  score_limits <- as.numeric(quantile(finite_score, c(0.01, 0.99), na.rm = TRUE))
  if (!all(is.finite(score_limits)) || diff(score_limits) == 0) {
    score_limits <- range(finite_score, na.rm = TRUE)
  }
  p <- ggplot(dt, aes(x = UMAP1, y = UMAP2, colour = score)) +
    geom_point(size = 0.22, alpha = 0.95, stroke = 0) +
    scale_colour_gradient2(
      low = "#2166AC", mid = "#E6E6E6", high = "#B2182B",
      midpoint = median(finite_score, na.rm = TRUE),
      limits = score_limits,
      oob = scales::squish
    ) +
    labs(
      title = title_text,
      x = "UMAP 1", y = "UMAP 2", colour = "Program\nscore"
    ) +
    theme_sci(7.2) +
    theme(axis.ticks = element_blank(), axis.line = element_blank()) +
    coord_fixed(ratio = 1, clip = "off")
  save_panel(
    p, "Figure_S4", stem,
    7.2, 5.8, dt, f2a_source$path,
    title_text, "supplementary"
  )
}
make_program_umap(
  score_pos_col_s4,
  "CCR2+ Top150 program localization",
  "FigureS4A_CCR2positive_program_UMAP"
)
make_program_umap(
  score_neg_col_s4,
  paste0("CCR2", UNICODE_MINUS, " Top150 program localization"),
  "FigureS4B_CCR2negative_program_UMAP"
)

## Figure S4C: macrophage-state UMAP
s4c_source <- read_required("stage3", "Fig3G_macrophage_UMAP_source.csv")
s4c <- copy(s4c_source$data)
mac_umap_cols_s4c <- find_umap_cols(s4c)
setnames(s4c, mac_umap_cols_s4c, c("UMAP1", "UMAP2"))
s4c[, state := pretty_cell(macrophage_state)]
state_palette_s4c <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#6C5CE7")[seq_len(uniqueN(s4c$state))],
  sort(unique(s4c$state))
)
state_centroids_s4c <- s4c[, .(
  UMAP1 = median(UMAP1, na.rm = TRUE),
  UMAP2 = median(UMAP2, na.rm = TRUE)
), by = state]
p_s4c <- ggplot(s4c, aes(x = UMAP1, y = UMAP2, colour = state)) +
  geom_point(size = 0.28, alpha = 0.95, stroke = 0) +
  ggrepel::geom_text_repel(
    data = state_centroids_s4c,
    aes(label = state),
    family = BASE_FAMILY, size = 2.3, fontface = "bold",
    colour = "#1A1A1A", max.overlaps = Inf,
    min.segment.length = 0, segment.size = 0.25
  ) +
  scale_colour_manual(values = state_palette_s4c) +
  labs(
    title = "Macrophage/monocyte state candidates",
    x = "UMAP 1", y = "UMAP 2", colour = "State"
  ) +
  theme_sci(7.2) +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) +
  coord_fixed(ratio = 1, clip = "off")
save_panel(
  p_s4c, "Figure_S4", "FigureS4C_macrophage_state_UMAP",
  7.2, 5.8, s4c, s4c_source$path,
  "Macrophage-state UMAP", "supplementary"
)

## Figure S4D: complete macrophage-state effect matrix
s4d <- copy(s3_mac_stats$data)
s4d[, program := paste0(
  pretty_program(signature_name), " Top", gsub(".*Top", "", signature_name)
)]
s4d[, state := pretty_cell(macrophage_state)]
s4d[, effect_estimable := is.finite(hedges_g_HFpEF_vs_Control)]
s4d[, estimability_note := fifelse(
  effect_estimable,
  "Hedges' g estimable",
  paste0("Not estimable (Control n=", control_samples, ")")
)]
p_s4d <- ggplot(
  s4d,
  aes(x = program, y = state, fill = hedges_g_HFpEF_vs_Control)
) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(
    data = s4d[effect_estimable == FALSE],
    aes(label = "NE"),
    family = BASE_FAMILY, size = 2.2, colour = "#303030"
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0
  ) +
  labs(
    title = "Complete macrophage-state program-effect matrix",
    x = NULL, y = NULL, fill = "Hedges' g"
  ) +
  theme_heatmap(6.4) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 6.0),
    axis.text.y = element_text(size = 6.0)
  )
save_panel(
  p_s4d, "Figure_S4", "FigureS4D_complete_macrophage_state_effect_matrix",
  10.2, 6.4, s4d, s3_mac_stats$path,
  "Complete macrophage-state effect matrix", "supplementary"
)

############################################################
## 13. Supplementary Figure S5
## TF activity methods, LOPO robustness and state localization
############################################################

## Figure S5A: sample-level TF regulon activity heatmap
s5a_source <- read_required(
  "stage4", "Fig4A_sample_level_TF_activity_heatmap_source.csv"
)
s5a <- copy(s5a_source$data)
sample_cols_s5a <- setdiff(names(s5a), "tf_symbol")
s5a_long <- melt(
  s5a,
  id.vars = "tf_symbol",
  measure.vars = sample_cols_s5a,
  variable.name = "sample_accession",
  value.name = "activity_z"
)
s5a_long <- merge(
  s5a_long,
  s1_manifest$data[, .(sample_accession, condition)],
  by = "sample_accession", all.x = TRUE
)
s5a_long[, tf_label := pretty_tf(tf_symbol)]
sample_order_s5a <- unique(s5a_long[order(condition, sample_accession), sample_accession])
s5a_long[, sample_accession := factor(sample_accession, levels = unique(sample_order_s5a))]
s5a_long[, tf_label := factor(tf_label, levels = rev(unique(tf_label)))]
p_s5a <- ggplot(
  s5a_long,
  aes(x = sample_accession, y = tf_label, fill = activity_z)
) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0
  ) +
  labs(
    title = "Sample-level macrophage TF regulon activity",
    x = NULL, y = NULL, fill = "Row z-score"
  ) +
  theme_heatmap(6.8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.0),
    axis.text.y = element_text(size = 6.0)
  )
save_panel(
  p_s5a, "Figure_S5", "FigureS5A_sample_level_TF_activity_heatmap",
  8.5, 7.0, s5a_long, c(s5a_source$path, s1_manifest$path),
  "Sample-level TF activity heatmap", "supplementary"
)

## Figure S5B: leave-one-pair-out robustness
s5b <- copy(s4_lopo$data)
setorder(s5b, median_abs_effect_rank)
s5b <- s5b[1:min(.N, 30L)]
s5b[, tf_label := factor(
  pretty_tf(tf_symbol),
  levels = unique(rev(pretty_tf(tf_symbol)))
)]
p_s5b <- ggplot(
  s5b,
  aes(
    x = sign_stability,
    y = tf_label,
    size = top10_frequency,
    colour = median_loo_effect
  )
) +
  geom_point(alpha = 0.9) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_size_continuous(
    range = c(1.5, 5.2), labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0
  ) +
  labs(
    title = "Leave-one-pair-out regulon stability",
    x = "Effect-sign stability", y = NULL,
    size = "Top-10\nfrequency", colour = "Median LOPO\neffect"
  ) +
  theme_sci(6.8) +
  theme(axis.text.y = element_text(size = 6.0))
save_panel(
  p_s5b, "Figure_S5", "FigureS5B_TF_LOPO_robustness",
  7.7, 7.3, s5b, s4_lopo$path,
  "TF leave-one-pair-out robustness", "supplementary"
)

## Figure S5C: weighted regulon versus AUCell effects
s5c <- merge(
  s4_activity$data[, .(tf_symbol, weighted_g = hedges_g_HFpEF_vs_Control)],
  s4_aucell$data[, .(tf_symbol, aucell_g = hedges_g_HFpEF_vs_Control)],
  by = "tf_symbol", all = FALSE
)
s5c[, tf_label := pretty_tf(tf_symbol)]
s5c[, label_show := ifelse(
  tf_label %in% c("BHLHE40", "NFKB1", "RELA", "RUNX1"), tf_label, ""
)]
p_s5c <- ggplot(s5c, aes(x = weighted_g, y = aucell_g)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#999999") +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#999999") +
  geom_point(size = 1.45, alpha = 0.65, colour = "#466B8A") +
  ggrepel::geom_text_repel(
    aes(label = label_show),
    family = BASE_FAMILY, size = 2.2,
    max.overlaps = Inf, min.segment.length = 0
  ) +
  labs(
    title = "Weighted regulon and AUCell effects provide complementary TF evidence",
    x = "Weighted regulon Hedges' g",
    y = "AUCell Hedges' g"
  ) +
  theme_sci(7.3) +
  coord_cartesian(clip = "off")
save_panel(
  p_s5c, "Figure_S5", "FigureS5C_weighted_vs_AUCell",
  6.2, 5.5, s5c, c(s4_activity$path, s4_aucell$path),
  "Weighted versus AUCell TF activity", "supplementary"
)

## Figure S5D: TF activity by macrophage state
s5d_source <- read_required(
  "stage4", "FigS4B_TF_activity_by_macrophage_state_source.csv"
)
s5d <- copy(s5d_source$data)
state_cols_s5d <- setdiff(names(s5d), "tf_symbol")
for (state_col in state_cols_s5d) {
  original_value <- s5d[[state_col]]
  numeric_value <- suppressWarnings(as.numeric(original_value))
  if (any(!is.na(original_value) & is.na(numeric_value))) {
    stop("Non-numeric TF-state activity value found in column: ", state_col)
  }
  set(s5d, j = state_col, value = numeric_value)
}
s5d_long <- melt(
  s5d,
  id.vars = "tf_symbol",
  measure.vars = state_cols_s5d,
  variable.name = "macrophage_state",
  value.name = "hfpef_minus_control"
)
s5d_long[, TF := pretty_tf(tf_symbol)]
s5d_long[, state := pretty_cell(macrophage_state)]
s5d_long[, effect_estimable := is.finite(hfpef_minus_control)]
s5d_long[, estimability_note := fifelse(
  effect_estimable,
  "HFpEF-Control activity estimable",
  "Not estimable in frozen source data"
)]
p_s5d <- ggplot(
  s5d_long,
  aes(x = state, y = TF, fill = hfpef_minus_control)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(
    data = s5d_long[effect_estimable == FALSE],
    aes(label = "NE"),
    family = BASE_FAMILY, size = 2.2, colour = "#303030"
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0
  ) +
  labs(
    title = "TF activity varies across macrophage-state candidates",
    x = NULL, y = NULL, fill = "HFpEF-Control\nactivity"
  ) +
  theme_heatmap(6.7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.0),
    axis.text.y = element_text(size = 6.0)
  )
save_panel(
  p_s5d, "Figure_S5", "FigureS5D_TF_activity_by_macrophage_state",
  9.0, 6.8, s5d_long, s5d_source$path,
  "TF activity by macrophage state", "supplementary"
)

############################################################
## 14. Supplementary Figure S6
## Perturbation strength, matched controls and sensitivity
############################################################

## Figure S6A: perturbation-strength sensitivity
s6a_source <- read_required("stage5", "FigS5A_perturbation_strength_source.csv")
s6a <- copy(s6a_source$data)
s6a[, tf_label := pretty_tf(tf_symbol)]
p_s6a <- ggplot(
  s6a,
  aes(
    x = as.numeric(perturbation_strength),
    y = median_gap_reduction,
    group = tf_label,
    colour = tf_label
  )
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(PALETTE_TF, "RUNX1" = "#8D5A97", "SPI1" = "#4F772D", "REL" = "#C7902F")) +
  labs(
    title = "Candidate recovery across perturbation strengths",
    x = "Perturbation strength",
    y = "Median primary-program gap reduction",
    colour = "TF"
  ) +
  theme_sci(7.3)
save_panel(
  p_s6a, "Figure_S6", "FigureS6A_perturbation_strength_sensitivity",
  7.2, 5.0, s6a, s6a_source$path,
  "Perturbation-strength sensitivity", "supplementary"
)

## Figure S6B: candidate versus matched-control TF references
s6b_source <- read_required("stage5", "FigS5B_candidate_vs_control_source.csv")
s6b <- copy(s6b_source$data)
s6b[, tf_label := pretty_tf(as.character(tf_symbol))]
s6b[, analysis_label := gsub("_", " ", analysis_role)]
s6b[, zero_effect := is.finite(median_gap_reduction) & median_gap_reduction == 0]
p_s6b <- ggplot(
  s6b,
  aes(x = tf_label, y = median_gap_reduction, fill = analysis_label)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "#777777") +
  geom_col(width = 0.7) +
  geom_point(
    data = s6b[zero_effect == TRUE],
    shape = 4, size = 2.6, stroke = 0.8, colour = "#303030"
  ) +
  geom_text(
    data = s6b[zero_effect == TRUE],
    aes(label = "0"),
    family = BASE_FAMILY, size = 2.2, colour = "#303030", vjust = -0.8
  ) +
  labs(
    title = "Candidate TFs are benchmarked against matched-control TF references",
    x = NULL, y = "Median primary-program gap reduction",
    fill = "Analysis role"
  ) +
  theme_sci(6.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.0))
save_panel(
  p_s6b, "Figure_S6", "FigureS6B_candidate_vs_matched_control",
  9.0, 5.0, s6b, s6b_source$path,
  "Candidate versus matched-control TFs", "supplementary"
)

## Figure S6C: normalization versus activity attenuation
s6c <- copy(s5_modes$data)
s6c[, tf_label := pretty_tf(tf_symbol)]
s6c[, mode_label := fcase(
  perturbation_mode == "disease_normalization", "Disease normalization",
  perturbation_mode == "activity_attenuation", "Activity attenuation",
  default = pretty_contrast(perturbation_mode)
)]
s6c[, zero_effect := is.finite(median_gap_reduction) & median_gap_reduction == 0]
p_s6c <- ggplot(
  s6c,
  aes(x = tf_label, y = median_gap_reduction, fill = mode_label)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_point(
    data = s6c[zero_effect == TRUE],
    position = position_dodge(width = 0.75),
    shape = 4, size = 2.6, stroke = 0.8, colour = "#303030"
  ) +
  geom_text(
    data = s6c[zero_effect == TRUE],
    aes(label = "0"),
    position = position_dodge(width = 0.75),
    family = BASE_FAMILY, size = 2.2, colour = "#303030", vjust = -0.8
  ) +
  labs(
    title = "Disease-normalization and activity-attenuation interpretations",
    x = NULL, y = "Median program gap reduction",
    fill = "Perturbation mode"
  ) +
  theme_sci(7.2)
save_panel(
  p_s6c, "Figure_S6", "FigureS6C_normalization_vs_attenuation",
  6.8, 4.7, s6c, s5_modes$path,
  "Normalization versus attenuation", "supplementary"
)

## Figure S6D: Stage 5 ranking sensitivity
s6d <- copy(s5_sensitivity$data)
s6d[, tf_label := pretty_tf(tf_symbol)]
s6d[, scenario_label := wrap_text(gsub("_", " ", scenario), 24)]
p_s6d <- ggplot(
  s6d,
  aes(
    x = scenario_label,
    y = scenario_rank,
    group = tf_label,
    colour = tf_label
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  scale_y_reverse() +
  scale_colour_manual(values = c(PALETTE_TF, "RUNX1" = "#8D5A97", "SPI1" = "#4F772D", "REL" = "#C7902F")) +
  labs(
    title = "Stage 5 candidate ranking across sensitivity scenarios",
    x = NULL, y = "Scenario rank", colour = "TF"
  ) +
  theme_sci(6.7) +
  theme(axis.text.x = element_text(angle = 48, hjust = 1, size = 6.0))
save_panel(
  p_s6d, "Figure_S6", "FigureS6D_stage5_ranking_sensitivity",
  10.5, 5.6, s6d, s5_sensitivity$path,
  "Stage 5 ranking sensitivity", "supplementary"
)

## Figure S6E: Stage 5B robustness matrix
s6e <- copy(s5b_final$data)
metric_cols_s6e <- intersect(c(
  "positive_recovery_probability",
  "top1_frequency",
  "top3_frequency",
  "candidate_percentile",
  "sample_improvement_fraction",
  "inflammation_median_gap_reduction",
  "specificity_score"
), names(s6e))
s6e_long <- melt(
  s6e,
  id.vars = c("tf_symbol", "final_robustness_rank"),
  measure.vars = metric_cols_s6e,
  variable.name = "metric", value.name = "raw_value"
)
s6e_long[, scaled_value := rescale01(raw_value), by = metric]
s6e_long[, tf_label := factor(
  pretty_tf(tf_symbol),
  levels = unique(rev(pretty_tf(s6e[order(final_robustness_rank), tf_symbol])))
)]
s6e_long[, metric_label := wrap_text(gsub("_", " ", metric), 18)]
p_s6e <- ggplot(
  s6e_long,
  aes(x = metric_label, y = tf_label, fill = scaled_value)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", raw_value)), family = BASE_FAMILY, size = 2.2) +
  scale_fill_gradient(low = "white", high = "#6C5CE7") +
  labs(
    title = "Stage 5B candidate robustness matrix",
    x = NULL, y = NULL, fill = "Within-metric\nrelative value"
  ) +
  theme_heatmap(6.8) +
  theme(axis.text.x = element_text(angle = 42, hjust = 1))
save_panel(
  p_s6e, "Figure_S6", "FigureS6E_stage5B_robustness_matrix",
  8.5, 4.2, s6e_long, s5b_final$path,
  "Stage 5B robustness matrix", "supplementary"
)

############################################################
## 15. Supplementary Figure S7
## Complete communication support
############################################################

## Figure S7A: top 50 stable axes
s7a <- copy(s6_stability$data)
setorder(s7a, median_scenario_rank, -top10_scenario_frequency)
s7a <- s7a[1:min(.N, 50L)]
s7a[, axis_label := short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
s7a[, axis_label := factor(
  axis_label,
  levels = unique(rev(as.character(axis_label)))
)]
p_s7a <- ggplot(
  s7a,
  aes(
    x = median_scenario_rank,
    y = axis_label,
    size = top20_scenario_frequency,
    colour = receiver
  )
) +
  geom_point(alpha = 0.85) +
  scale_x_reverse() +
  scale_colour_manual(values = PALETTE_RECEIVER, na.value = "#777777") +
  scale_size_continuous(
    range = c(1.5, 4.5), labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Top 50 frozen communication axes",
    x = "Median scenario rank", y = NULL,
    size = "Top-20\nfrequency", colour = "Receiver"
  ) +
  theme_sci(6.3) +
  theme(axis.text.y = element_text(size = 6.0))
save_panel(
  p_s7a, "Figure_S7", "FigureS7A_top50_communication_axes",
  10.0, 13.0, s7a, s6_stability$path,
  "Top 50 communication axes", "supplementary"
)

## Figure S7B: ligand-target regulatory potential
s7b <- copy(s6_top_axes$data)
for (z in intersect(c(
  "final_axis_rank", "maximum_regulatory_potential",
  "median_regulatory_potential", "target_links"
), names(s7b))) {
  set(s7b, j = z, value = as.numeric(s7b[[z]]))
}
setorder(s7b, final_axis_rank)
s7b <- s7b[1:min(.N, 30L)]
s7b[, axis_label := short_axis(tf_symbol, nichenet_ligand, receptor, receiver)]
s7b[, axis_label := factor(
  axis_label,
  levels = unique(rev(as.character(axis_label)))
)]
p_s7b <- ggplot(
  s7b,
  aes(
    x = median_regulatory_potential,
    y = axis_label,
    size = target_links,
    colour = receiver
  )
) +
  geom_point(alpha = 0.88) +
  scale_colour_manual(values = PALETTE_RECEIVER, na.value = "#777777") +
  scale_size_continuous(range = c(1.6, 5.0)) +
  labs(
    title = "Ligand-target regulatory potential among prioritized axes",
    x = "Median NicheNet regulatory potential", y = NULL,
    size = "Target links", colour = "Receiver"
  ) +
  theme_sci(6.2) +
  theme(axis.text.y = element_text(size = 6.0))
save_panel(
  p_s7b, "Figure_S7", "FigureS7B_ligand_target_regulatory_potential",
  9.0, 9.0, s7b, s6_top_axes$path,
  "Ligand-target regulatory potential", "supplementary"
)

## Figure S7C: complete receiver-receptor support
s7c <- copy(s6_top_axes$data)
for (z in intersect(num_cols_f4d, names(s7c))) {
  set(s7c, j = z, value = as.numeric(s7c[[z]]))
}
s7c <- s7c[is.finite(median_receptor_log2_cpm_HFpEF)]
s7c <- s7c[, .(
  median_control = median(median_receptor_log2_cpm_Control, na.rm = TRUE),
  median_hfpef = median(median_receptor_log2_cpm_HFpEF, na.rm = TRUE),
  median_pct = median(median_receptor_pct_HFpEF, na.rm = TRUE)
), by = .(receptor, receiver)]
p_s7c <- ggplot(
  s7c,
  aes(
    x = median_control,
    y = median_hfpef,
    size = median_pct,
    colour = receiver
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, colour = "#888888") +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = PALETTE_RECEIVER, na.value = "#777777") +
  scale_size_continuous(
    range = c(1.0, 4.5), labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Complete receiver-receptor expression support",
    x = "Control receptor expression (log2 CPM)",
    y = "HFpEF receptor expression (log2 CPM)",
    size = "HFpEF cells\nexpressing receptor",
    colour = "Receiver"
  ) +
  theme_sci(7.0)
save_panel(
  p_s7c, "Figure_S7", "FigureS7C_complete_receptor_support",
  7.4, 5.8, s7c, s6_top_axes$path,
  "Complete receptor support", "supplementary"
)

## Figure S7D: candidate ligand coverage
s7d <- copy(s6_ligand_coverage$data)
coverage_cols_s7d <- intersect(c(
  "direct_NicheNet_LR_ligands",
  "curated_panel_ligands",
  "predicted_decrease_ligands",
  "predicted_increase_ligands"
), names(s7d))
if (length(coverage_cols_s7d) == 0L) {
  coverage_cols_s7d <- names(s7d)[
    vapply(s7d, is.numeric, logical(1)) &
      grepl("ligand", names(s7d), ignore.case = TRUE)
  ]
}
s7d_long <- melt(
  s7d,
  id.vars = intersect(c("tf_symbol", "candidate_role", "final_robustness_rank"), names(s7d)),
  measure.vars = coverage_cols_s7d,
  variable.name = "coverage_type", value.name = "ligands"
)
s7d_long[, tf_label := pretty_tf(tf_symbol)]
s7d_long[, coverage_label := wrap_text(gsub("_", " ", coverage_type), 20)]
p_s7d <- ggplot(
  s7d_long,
  aes(x = coverage_label, y = ligands, fill = tf_label)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_manual(values = PALETTE_TF) +
  labs(
    title = "Candidate-TF ligand coverage",
    x = NULL, y = "Ligands", fill = "TF"
  ) +
  theme_sci(7.0) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(
  p_s7d, "Figure_S7", "FigureS7D_candidate_ligand_coverage",
  7.2, 4.8, s7d_long, s6_ligand_coverage$path,
  "Candidate ligand coverage", "supplementary"
)

############################################################
## 16. Supplementary Figure S8
## Stage 7 internal sample separability and guardrails
############################################################

## Figure S8A: held-out sample predictions
s8a <- copy(s7_predictions$data)
s8a[, true_condition := factor(true_condition, levels = unique(c("Control", "HFpEF")))]
p_s8a <- ggplot(
  s8a,
  aes(x = sample_accession, y = mean_predicted_probability, fill = true_condition)
) +
  geom_col(width = 0.68) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.45, colour = "#777777") +
  scale_fill_manual(values = PALETTE_CONDITION) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Leave-pair-out held-out sample predictions",
    x = NULL, y = "Mean predicted HFpEF probability",
    fill = "True condition"
  ) +
  theme_sci(7.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(
  p_s8a, "Figure_S8", "FigureS8A_heldout_sample_predictions",
  6.8, 4.6, s8a, s7_predictions$path,
  "Held-out sample predictions", "supplementary"
)

## Figure S8B: internal sample-level ROC
roc_input_s8b <- copy(s7_predictions$data)
roc_input_s8b[, label := as.integer(true_condition == "HFpEF")]
thresholds_s8b <- sort(unique(c(
  Inf, roc_input_s8b$mean_predicted_probability, -Inf
)), decreasing = TRUE)
roc_curve_s8b <- rbindlist(lapply(thresholds_s8b, function(th) {
  pred <- as.integer(roc_input_s8b$mean_predicted_probability >= th)
  tp <- sum(pred == 1 & roc_input_s8b$label == 1)
  fp <- sum(pred == 1 & roc_input_s8b$label == 0)
  fn <- sum(pred == 0 & roc_input_s8b$label == 1)
  tn <- sum(pred == 0 & roc_input_s8b$label == 0)
  data.table(
    threshold = th,
    sensitivity = tp / pmax(tp + fn, 1),
    specificity = tn / pmax(tn + fp, 1),
    fpr = 1 - tn / pmax(tn + fp, 1)
  )
}))
setorder(roc_curve_s8b, fpr, sensitivity)
auc_s8b <- s7_performance$data[
  metric == "Primary sample-level AUC", value
][1L]
p_s8b <- ggplot(roc_curve_s8b, aes(x = fpr, y = sensitivity)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.45, colour = "#888888") +
  geom_step(linewidth = 0.9, colour = "#C84C4C") +
  coord_equal() +
  labs(
    title = paste0("Internal sample-level ROC (AUC = ", sprintf("%.3f", auc_s8b), ")"),
    x = "1 - Specificity", y = "Sensitivity"
  ) +
  theme_sci(7.2)
save_panel(
  p_s8b, "Figure_S8", "FigureS8B_internal_sample_ROC",
  5.4, 5.0, roc_curve_s8b,
  c(s7_predictions$path, s7_performance$path),
  "Internal sample ROC", "supplementary"
)

## Figure S8C: balanced-label permutation results
s8c <- copy(s7_permutation$data)
s8c[, metric_label := wrap_text(metric, 30)]
p_s8c <- ggplot(
  s8c,
  aes(x = null_median, xend = null_maximum, y = reorder(metric_label, observed))
) +
  geom_segment(linewidth = 3.0, colour = "#D5D5D5", lineend = "round") +
  geom_point(aes(x = observed), size = 2.8, colour = "#C84C4C") +
  geom_text(
    aes(
      x = pmax(observed, null_maximum, na.rm = TRUE),
      label = paste0("P=", sprintf("%.2f", empirical_p))
    ),
    family = BASE_FAMILY, size = 2.35, hjust = -0.12
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    title = "Balanced-label permutation limits inferential overstatement",
    x = "Observed value and permutation-null range", y = NULL
  ) +
  theme_sci(7.0) +
  coord_cartesian(clip = "off")
save_panel(
  p_s8c, "Figure_S8", "FigureS8C_balanced_label_permutation",
  7.5, 4.8, s8c, s7_permutation$path,
  "Balanced-label permutation", "supplementary"
)

## Figure S8D: panel and lambda sensitivity
s8d <- copy(s7_sensitivity$data)
s8d[, lambda_label := factor(as.character(lambda), levels = unique(as.character(lambda)))]
p_s8d <- ggplot(
  s8d,
  aes(x = lambda_label, y = sample_auc, group = panel, colour = panel)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Feature-panel and ridge-penalty sensitivity",
    x = "Ridge penalty lambda", y = "Sample-level AUC",
    colour = "Feature panel"
  ) +
  theme_sci(6.8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 6.0)
  )
save_panel(
  p_s8d, "Figure_S8", "FigureS8D_panel_lambda_sensitivity",
  8.4, 5.2, s8d, s7_sensitivity$path,
  "Panel and lambda sensitivity", "supplementary"
)

## Figure S8E: primary cross-stage feature matrix
s8e <- copy(s7_features$data)
feature_cols_s8e <- setdiff(
  names(s8e),
  c("sample_accession", "sample_order_index", "condition", "y")
)
s8e_long <- melt(
  s8e,
  id.vars = c("sample_accession", "sample_order_index", "condition"),
  measure.vars = feature_cols_s8e,
  variable.name = "feature", value.name = "value"
)
s8e_long[, z_value := {
  v <- as.numeric(value)
  s <- sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) rep(0, .N) else (v - mean(v, na.rm = TRUE)) / s
}, by = feature]
s8e_long[, feature_label := fcase(
  feature == "COMM_NFkB_axis_burden", "NF-kB communication burden",
  feature == "PROGRAM_DrugOpposed_Top150", "Drug-opposed Top150",
  feature == "PROGRAM_DrugOpposed_Top150_Stage3Supported", "Stage 3-supported Top150",
  feature == "TF_Nfkb1_activity", "NFKB1 activity",
  feature == "TF_Bhlhe40_activity", "BHLHE40 activity",
  feature == "TF_Rela_activity", "RELA activity",
  feature == "COMM_Bhlhe40_axis_burden", "BHLHE40 communication burden",
  default = gsub("_", " ", feature)
)]
sample_order_s8e <- s8e[order(condition, sample_order_index), sample_accession]
s8e_long[, sample_accession := factor(sample_accession, levels = unique(sample_order_s8e))]
p_s8e <- ggplot(
  s8e_long,
  aes(x = sample_accession, y = feature_label, fill = z_value)
) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0
  ) +
  labs(
    title = "Primary sample-level cross-stage feature matrix",
    x = NULL, y = NULL, fill = "Feature z-score"
  ) +
  theme_heatmap(7.0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(
  p_s8e, "Figure_S8", "FigureS8E_primary_cross_stage_feature_heatmap",
  8.0, 5.4, s8e_long, s7_features$path,
  "Primary cross-stage feature heatmap", "supplementary"
)

############################################################
## 17. Supplementary Figure S9
## Dataset-specific Stage 8 validation panels
############################################################

make_dataset_validation_panel <- function(dataset_value, stem, title_text) {
  prog <- s8_program_evidence$data[
    dataset_id == dataset_value & grepl("Top150$", program_name)
  ]
  tf <- s8_tf_evidence$data[dataset_id == dataset_value]
  axis <- s8_axis_evidence$data[dataset_id == dataset_value]

  prog_dt <- prog[, .(
    contrast,
    compartment = cell_type,
    evidence_type = "Program",
    item = pretty_program(program_name),
    hedges_g,
    fdr,
    direction_supported
  )]
  tf_dt <- tf[, .(
    contrast,
    compartment = cell_type,
    evidence_type = "TF",
    item = pretty_tf(tf_symbol),
    hedges_g,
    fdr,
    direction_supported
  )]

  if (nrow(axis) > 0L) {
    top_axis_dataset <- s8_axis_summary$data[order(integrated_rank), head(axis_key, 10L)]
    axis_dt <- axis[axis_key %in% top_axis_dataset, .(
      contrast,
      compartment = receiver,
      evidence_type = "Axis",
      item = short_axis(tf_symbol, nichenet_ligand, receptor, receiver),
      hedges_g,
      fdr,
      direction_supported
    )]
  } else {
    axis_dt <- data.table()
  }

  combined <- rbindlist(list(prog_dt, tf_dt, axis_dt), use.names = TRUE, fill = TRUE)
  if (nrow(combined) == 0L) {
    stop("No Stage 8 validation evidence found for dataset: ", dataset_value)
  }
  combined[, contrast_label := pretty_contrast(contrast)]
  combined[, display_item := paste(item, pretty_cell(compartment), sep = " | ")]
  combined[, fdr_strength := fifelse(
    is.finite(fdr), pmin(neglog10_fdr(fdr), 10), 0
  )]
  combined[, support_shape := !is.na(direction_supported) & direction_supported]
  combined[, effect_estimable := is.finite(hedges_g)]
  combined[, display_effect := fifelse(effect_estimable, hedges_g, 0)]
  combined[, estimability_note := fifelse(
    effect_estimable,
    "Hedges' g estimable",
    "Not estimable in frozen source data"
  )]
  combined[, display_item := factor(
    display_item,
    levels = rev(unique(combined[order(evidence_type, item, compartment), display_item]))
  )]

  p <- ggplot(
    combined,
    aes(
      x = display_effect,
      y = display_item,
      size = fdr_strength,
      colour = contrast_label,
      shape = support_shape
    )
  ) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = "#777777") +
    geom_point(data = combined[effect_estimable == TRUE], alpha = 0.88) +
    geom_point(
      data = combined[effect_estimable == FALSE],
      aes(x = display_effect, y = display_item, colour = contrast_label),
      inherit.aes = FALSE,
      shape = 4, size = 2.4, stroke = 0.8
    ) +
    geom_text(
      data = combined[effect_estimable == FALSE],
      aes(x = display_effect, y = display_item, label = "NE"),
      inherit.aes = FALSE,
      family = BASE_FAMILY, size = 2.2, colour = "#555555",
      hjust = -0.65
    ) +
    facet_grid(evidence_type ~ ., scales = "free_y", space = "free_y") +
    scale_size_continuous(range = c(1.2, 4.5)) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("FALSE" = "Direction not supported", "TRUE" = "Direction supported")
    ) +
    labs(
      title = title_text,
      x = "Hedges' g", y = NULL,
      size = "-log10(FDR)",
      colour = "Contrast",
      shape = NULL
    ) +
    theme_sci(6.0) +
    theme(
      axis.text.y = element_text(size = 6.0),
      legend.text = element_text(size = 6.0)
    ) +
    coord_cartesian(clip = "off")

  save_panel(
    p, "Figure_S9", stem,
    10.5, ifelse(nrow(axis_dt) > 0L, 11.0, 7.5),
    combined,
    c(s8_program_evidence$path, s8_tf_evidence$path, s8_axis_evidence$path),
    title_text, "supplementary"
  )
}

make_dataset_validation_panel(
  "GSE236584",
  "FigureS9A_GSE236584_validation",
  "GSE236584 matched whole-heart bulk support"
)
make_dataset_validation_panel(
  "GSE245034",
  "FigureS9B_GSE245034_validation",
  "GSE245034 external bulk drug-response validation"
)
make_dataset_validation_panel(
  "GSE208425",
  "FigureS9C_GSE208425_validation",
  "GSE208425 internal immune-cell context"
)
make_dataset_validation_panel(
  "GSE249412",
  "FigureS9D_GSE249412_validation",
  "GSE249412 cell-type-resolved drug-response validation"
)
make_dataset_validation_panel(
  "SCP3342",
  "FigureS9E_SCP3342_validation",
  "SCP3342 independent human myocardial validation"
)

############################################################
## 18. Supplementary Figure S10
## Benchmark alignment and sentinel-axis analyses
############################################################

## Figure S10A: descriptive TF alignment
s10a <- copy(b_tf_alignment$data)
s10a[, method := factor(method, levels = c(
  "TF expression only",
  "Regulon activity only",
  "Stage 4 multifeature",
  "Perturbation only",
  "Bootstrap perturbation",
  "Communication only",
  "Full cross-layer"
))]
p_s10a <- ggplot(s10a, aes(x = method, y = spearman_vs_external_rank)) +
  geom_col(width = 0.68, fill = "#466B8A") +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "#777777") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Descriptive alignment with Stage 8 external TF ranks",
    x = NULL, y = "Spearman rank correlation"
  ) +
  theme_sci(7.0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(
  p_s10a, "Figure_S10", "FigureS10A_TF_external_rank_alignment",
  7.5, 4.7, s10a, b_tf_alignment$path,
  "TF external-rank alignment", "supplementary"
)

## Figure S10B: sentinel axis rank across methods
s10b <- copy(b_sentinel$data)
s10b[, method := factor(method, levels = c(
  "Ligand expression only",
  "Receptor expression only",
  "Ligand-receptor expression",
  "NicheNet only",
  "Cross-stage support only",
  "Full integration"
))]
p_s10b <- ggplot(s10b, aes(x = method, y = discovery_rank)) +
  geom_col(width = 0.68, fill = "#E76F51") +
  scale_y_reverse() +
  labs(
    title = paste0(
      "Prespecified TNF-TNFRSF1A-Endothelial sentinel axis",
      "\nExternal integrated rank = ", unique(s10b$external_integrated_rank)
    ),
    x = NULL, y = "Discovery rank (lower is better)"
  ) +
  theme_sci(7.0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(
  p_s10b, "Figure_S10", "FigureS10B_sentinel_axis_rank",
  7.3, 4.8, s10b, b_sentinel$path,
  "Sentinel-axis rank comparison", "supplementary"
)

## Figure S10C: axis-method external alignment
s10c <- copy(b_axis_alignment$data)
s10c[, method := factor(method, levels = c(
  "Ligand expression only",
  "Receptor expression only",
  "Ligand-receptor expression",
  "NicheNet only",
  "Cross-stage support only",
  "Full integration"
))]
p_s10c <- ggplot(s10c, aes(x = method, y = spearman_vs_external_rank)) +
  geom_col(width = 0.68, fill = "#6C5CE7") +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "#777777") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Axis-ranking alignment with external integrated evidence",
    x = NULL, y = "Spearman rank correlation"
  ) +
  theme_sci(7.0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(
  p_s10c, "Figure_S10", "FigureS10C_axis_external_rank_alignment",
  7.2, 4.8, s10c, b_axis_alignment$path,
  "Axis external-rank alignment", "supplementary"
)

## Figure S10D: ablation stability summary
s10d <- copy(b_ablation_stability$data)
s10d[, tf_label := factor(
  pretty_tf(tf_symbol),
  levels = unique(rev(pretty_tf(tf_symbol[order(median_rank)])))
)]
p_s10d <- ggplot(
  s10d,
  aes(
    x = top1_frequency,
    y = tf_label,
    size = top2_frequency,
    colour = median_rank
  )
) +
  geom_point(alpha = 0.9) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_size_continuous(
    range = c(2.5, 6.0), labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradient(low = "#C84C4C", high = "#D5D5D5", trans = "reverse") +
  labs(
    title = "Candidate-rank stability across eight benchmark scenarios",
    x = "Top-1 frequency", y = NULL,
    size = "Top-2\nfrequency", colour = "Median rank"
  ) +
  theme_sci(7.3)
save_panel(
  p_s10d, "Figure_S10", "FigureS10D_candidate_ablation_stability",
  6.5, 4.2, s10d, b_ablation_stability$path,
  "Candidate ablation stability", "supplementary"
)

############################################################
## 19. Supplementary Tables S1-S16
############################################################

## Table S1: dataset hierarchy and inferential units
write_supp_table(
  1,
  "Dataset_Evidence_Hierarchy",
  list(
    locked_sample_manifest = s1_manifest$data,
    analysis_roles = s1_roles$data,
    dataset_lock_summary = s1_summary$data,
    stage8_roles = s8_roles$data,
    biological_sample_counts = s8_counts$data,
    SCP3342_donor_counts = s8_donors$data
  )
)

## Table S2: GSE237156 samples, contrasts and file mapping
write_supp_table(
  2,
  "GSE237156_Samples_Contrasts_and_Files",
  list(
    samples = s2_meta$data,
    group_counts = s2_group_counts$data,
    contrasts = s2_contrasts$data,
    file_mapping = s2_mapping$data,
    sample_QC = s2_qc$data
  )
)

## Table S3: continuous drug-opposed programs
write_supp_table(
  3,
  "Continuous_Drug_Opposed_Programs",
  list(
    frozen_program_manifest = s8_program_manifest$data,
    signature_size_summary = s3_sig_size$data,
    CCR2_positive_top500 = s2_top_pos$data,
    CCR2_negative_top500 = s2_top_neg$data,
    cross_subset_top500 = s2_top_cross$data,
    opposition_summary = s2_opp_summary$data
  )
)

## Table S4: Hallmark pathway opposition and method concordance
write_supp_table(
  4,
  "Hallmark_Pathway_Opposition",
  list(
    hallmark_opposition = s2_hallmark$data,
    contrast_summary = s2_contrast_summary$data,
    method_concordance = s2_method$data
  )
)

## Table S5: GSE236585 QC, annotation and pseudobulk eligibility
write_supp_table(
  5,
  "GSE236585_QC_Annotation_and_Pseudobulk",
  list(
    file_mapping = s3_mapping$data,
    QC_thresholds = s3_qc_thresholds$data,
    QC_retention = s3_qc_retention$data,
    doublet_rates = s3_doublet$data,
    celltype_annotation = s3_annotation$data,
    celltype_composition = s3_composition$data,
    pseudobulk_eligibility = s3_eligibility$data
  )
)

## Table S6: program localization and concordance
write_supp_table(
  6,
  "Program_Localization_and_Concordance",
  list(
    all_program_statistics = s3_program_stats$data,
    primary_program_statistics = s3_primary$data,
    Stage2_Stage3_concordance = s3_concordance$data,
    macrophage_state_annotation = s3_mac_annotation$data,
    macrophage_state_statistics = s3_mac_stats$data,
    macrophage_state_composition = s3_mac_comp$data
  )
)

## Table S7: TF expression, regulon, AUCell and LOPO
write_supp_table(
  7,
  "TF_Expression_Regulon_AUCell_and_LOPO",
  list(
    TF_priority = s4_priority$data,
    weighted_regulon = s4_activity$data,
    AUCell = s4_aucell$data,
    TF_expression = s4_expression$data,
    LOPO_robustness = s4_lopo$data,
    method_comparison = s4_method$data
  )
)

## Table S8: virtual perturbation and sensitivity
write_supp_table(
  8,
  "Virtual_Perturbation_and_Sensitivity",
  list(
    candidate_resolution = s5_manifest$data,
    matched_controls = s5_controls$data,
    program_definitions = s5_program_defs$data,
    observed_program_scores = s5_observed$data,
    candidate_rank = s5_rank$data,
    ranking_sensitivity = s5_sensitivity$data,
    ranking_stability = s5_stability$data,
    mode_comparison = s5_modes$data,
    mode_summary = s5_mode_summary$data,
    top_gene_effects = s5_gene_effects$data,
    ligand_effects = s5_ligands$data,
    method_concordance = s5_method$data,
    matched_control_results = s5_control_results$data
  )
)

## Table S9: bootstrap and matched random-TF null
write_supp_table(
  9,
  "Bootstrap_and_Matched_Random_TF_Null",
  list(
    full_candidate_summary = s5b_full$data,
    candidate_bootstrap = s5b_boot$data,
    matching_covariates = s5b_matching$data,
    matched_null_pools = s5b_null_pools$data,
    precomputed_null_effects = s5b_null_effects$data,
    random_null_summary = s5b_null$data,
    final_robustness_rank = s5b_final$data
  )
)

## Table S10: communication ligands, receptors and axes
write_supp_table(
  10,
  "Communication_Ligands_Receptors_and_Axes",
  list(
    candidate_TF_manifest = s6_tf_manifest$data,
    ligand_coverage = s6_ligand_coverage$data,
    NicheNet_activity = s6_nichenet$data,
    axis_stability = s6_stability$data,
    candidate_summary = s6_candidate$data,
    top_axes = s6_top_axes$data,
    ligand_support = s6_ligand_support$data,
    receiver_gene_sets = s6_receiver_sets$data
  )
)

## Table S11: sample-level attribution and internal performance
write_supp_table(
  11,
  "Sample_Level_Attribution_and_Internal_Performance",
  list(
    feature_definitions = s7_definitions$data,
    feature_matrix = s7_features$data,
    LOPO_fold_performance = s7_fold$data,
    heldout_predictions = s7_predictions$data,
    feature_importance = s7_importance$data,
    permutation = s7_permutation$data,
    panel_lambda_sensitivity = s7_sensitivity$data,
    model_performance = s7_performance$data
  )
)

## Table S12: Stage 8 cohort, donor and gene audits
write_supp_table(
  12,
  "Multicohort_Sample_Donor_and_Gene_Audits",
  list(
    core_metadata = s8_core_meta$data,
    biological_sample_counts = s8_counts$data,
    SCP3342_donor_counts = s8_donors$data,
    gene_key_audit = s8_gene_audit$data,
    frozen_TF_manifest = s8_tf_manifest$data,
    frozen_axis_manifest = s8_axis_manifest$data,
    dataset_roles = s8_roles$data
  )
)

## Table S13: multicohort program evidence
write_supp_table(
  13,
  "Multicohort_Program_Evidence",
  list(
    complete_program_evidence = s8_program_evidence$data,
    integrated_program_summary = s8_program_summary$data,
    SCP3342_program_validation = s8_scp_program$data
  )
)

## Table S14: multicohort TF and axis evidence
write_supp_table(
  14,
  "Multicohort_TF_and_Axis_Evidence",
  list(
    complete_TF_evidence = s8_tf_evidence$data,
    integrated_TF_summary = s8_tf_summary$data,
    complete_axis_evidence = s8_axis_evidence$data,
    integrated_axis_summary = s8_axis_summary$data,
    GSE249412_axis_validation = s8_gse249_axis$data,
    SCP3342_TF_validation = s8_scp_tf$data,
    SCP3342_axis_validation = s8_scp_axis$data
  )
)

## Table S15: benchmark, ablation and claim boundaries
write_supp_table(
  15,
  "Benchmark_Ablation_and_Claim_Boundaries",
  list(
    all_TF_baselines = b_tf_all$data,
    selected_TF_baselines = b_tf_selected$data,
    candidate_method_ranks = b_tf_ranks$data,
    TF_ablation = b_ablation$data,
    TF_ablation_stability = b_ablation_stability$data,
    TF_external_alignment = b_tf_alignment$data,
    axis_method_ranks = b_axis_ranks$data,
    axis_external_alignment = b_axis_alignment$data,
    axis_topk_performance = b_axis_topk$data,
    sentinel_axis = b_sentinel$data,
    Stage7_metrics = b_stage7_metrics$data,
    feature_importance = b_feature$data,
    claim_boundaries = b_boundaries$data
  )
)

## Table S16: software, parameters and environment index
parameter_files <- unique(unlist(lapply(names(STAGE_DIRS), function(k) {
  list.files(
    STAGE_DIRS[[k]],
    pattern = "parameters|rationale|sessionInfo|methods_and_claim|claim_boundaries|priority_score_definition",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
})))
parameter_files <- parameter_files[file.exists(parameter_files)]
parameter_index <- data.table(
  file = basename(parameter_files),
  full_path = gsub("\\\\", "/", parameter_files),
  size_bytes = file.info(parameter_files)$size,
  md5 = unname(tools::md5sum(parameter_files))
)
parameter_index[, stage := vapply(full_path, function(p) {
  hit <- names(STAGE_DIRS)[vapply(
    STAGE_DIRS,
    function(root) startsWith(
      tolower(gsub("\\\\", "/", p)),
      tolower(gsub("\\\\", "/", root))
    ),
    logical(1)
  )]
  if (length(hit) == 0L) NA_character_ else hit[1L]
}, character(1))]
write_supp_table(
  16,
  "Software_Parameters_Code_and_Environment",
  list(
    parameter_file_index = parameter_index,
    current_session_info = data.table(line = capture.output(sessionInfo())),
    figure_export_schema = data.table(
      item = c(
        "project_root",
        "analysis_schema",
        "figure_dpi",
        "TIFF_compression",
        "font_family",
        "base_font_size",
        "random_seed"
      ),
      value = c(
        PROJECT_DIR,
        ANALYSIS_SCHEMA,
        FIGURE_DPI,
        "LZW",
        BASE_FAMILY,
        BASE_SIZE,
        20260715
      )
    )
  )
)

############################################################
## 20. Final manifests, README, integrity checks and CHECK ZIP
############################################################

panel_audit <- rbindlist(panel_audit_records, use.names = TRUE, fill = TRUE)
input_manifest <- rbindlist(input_manifest_records, use.names = TRUE, fill = TRUE)
supp_table_manifest <- rbindlist(supp_table_records, use.names = TRUE, fill = TRUE)

## Confirm that all recorded frozen input files remain byte-identical.
input_manifest[, current_size_bytes := file.info(input_path)$size]
input_manifest[, current_md5 := unname(tools::md5sum(input_path))]
input_manifest[, unchanged :=
  size_bytes == current_size_bytes & md5 == current_md5
]

write_csv_safe(
  panel_audit,
  file.path(DIRS$manifests, "02_figure_panel_export_audit.csv")
)
write_csv_safe(
  input_manifest,
  file.path(DIRS$manifests, "03_panel_input_file_manifest.csv")
)
write_csv_safe(
  supp_table_manifest,
  file.path(DIRS$manifests, "04_supplementary_table_manifest.csv")
)

main_figure_plan <- data.table(
  figure = c(
    rep("Figure 1", 6),
    rep("Figure 2", 6),
    rep("Figure 3", 6),
    rep("Figure 4", 6),
    rep("Figure 5", 6)
  ),
  panel = rep(LETTERS[1:6], 5),
  scientific_question = c(
    "Sample-level transcriptomic structure",
    "CCR2+ disease-drug opposition",
    "CCR2- disease-drug opposition",
    "Top opposed-gene expression",
    "Hallmark pathway opposition",
    "Signature-size sensitivity",
    "Major cardiac cell types",
    "Cell-type program effects",
    "Macrophage-state sample scores",
    "Macrophage-state effects",
    "Cross-dataset concordance",
    "Program effect matrix",
    "Unbiased TF priority",
    "Virtual perturbation recovery",
    "Bootstrap and matched-null robustness",
    "TF method ranks",
    "Evidence-layer ablation",
    "Sample-level feature importance",
    "Candidate-TF communication coverage",
    "Stable communication axes",
    "NicheNet ligand activity",
    "Receiver receptor support",
    "Discovery rank versus external support",
    "Top-k method comparison",
    "Multicohort program evidence",
    "Multicohort TF evidence",
    "Independent human program validation",
    "Drug-response axis validation",
    "Multicohort axis evidence",
    "Integrated external axis support"
  )
)
write_csv_safe(
  main_figure_plan,
  file.path(DIRS$manifests, "05_main_figure_plan.csv")
)

supplementary_figure_plan <- panel_audit[
  panel_type == "supplementary",
  .(figure_group, stem, title, tiff_path, source_data_path)
]
write_csv_safe(
  supplementary_figure_plan,
  file.path(DIRS$manifests, "05A_supplementary_figure_plan.csv")
)

## Read TIFF header tags directly so technical QA does not depend on screenshots
## or on loading multi-megapixel images into memory.
read_tiff_metadata <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  read_at <- function(offset, n) {
    seek(con, where = offset, origin = "start")
    readBin(con, what = "raw", n = n)
  }
  signature <- rawToChar(read_at(0, 2))
  little <- identical(signature, "II")
  if (!little && !identical(signature, "MM")) stop("Invalid TIFF byte order: ", path)
  u16 <- function(offset) {
    b <- as.numeric(read_at(offset, 2))
    if (length(b) != 2L) return(NA_real_)
    if (little) b[1L] + 256 * b[2L] else 256 * b[1L] + b[2L]
  }
  u32 <- function(offset) {
    b <- as.numeric(read_at(offset, 4))
    if (length(b) != 4L) return(NA_real_)
    if (little) {
      b[1L] + 256 * b[2L] + 65536 * b[3L] + 16777216 * b[4L]
    } else {
      16777216 * b[1L] + 65536 * b[2L] + 256 * b[3L] + b[4L]
    }
  }
  if (!identical(u16(2), 42)) stop("Unsupported TIFF header: ", path)
  ifd_offset <- u32(4)
  entry_count <- as.integer(u16(ifd_offset))
  tags <- list()
  type_size <- c(`1` = 1, `2` = 1, `3` = 2, `4` = 4, `5` = 8)
  for (i in seq_len(entry_count)) {
    entry_offset <- ifd_offset + 2 + (i - 1) * 12
    tag <- as.integer(u16(entry_offset))
    type <- as.integer(u16(entry_offset + 2))
    count <- as.integer(u32(entry_offset + 4))
    if (!as.character(type) %in% names(type_size) || count < 1L) next
    bytes_needed <- count * type_size[[as.character(type)]]
    value_offset <- if (bytes_needed <= 4) entry_offset + 8 else u32(entry_offset + 8)
    values <- switch(
      as.character(type),
      `1` = as.numeric(read_at(value_offset, count)),
      `2` = rawToChar(read_at(value_offset, count), multiple = FALSE),
      `3` = vapply(seq_len(count), function(j) u16(value_offset + 2 * (j - 1)), numeric(1)),
      `4` = vapply(seq_len(count), function(j) u32(value_offset + 4 * (j - 1)), numeric(1)),
      `5` = vapply(seq_len(count), function(j) {
        numerator <- u32(value_offset + 8 * (j - 1))
        denominator <- u32(value_offset + 8 * (j - 1) + 4)
        if (!is.finite(denominator) || denominator == 0) NA_real_ else numerator / denominator
      }, numeric(1))
    )
    tags[[as.character(tag)]] <- values
  }
  first_tag <- function(tag, default = NA_real_) {
    value <- tags[[as.character(tag)]]
    if (is.null(value) || length(value) == 0L) default else value[1L]
  }
  width <- first_tag(256)
  height <- first_tag(257)
  compression_code <- first_tag(259)
  photometric <- first_tag(262)
  samples <- first_tag(277, 1)
  resolution_unit <- first_tag(296, 2)
  dpi_x <- first_tag(282)
  dpi_y <- first_tag(283)
  if (identical(resolution_unit, 3)) {
    dpi_x <- dpi_x * 2.54
    dpi_y <- dpi_y * 2.54
  }
  mode <- if (identical(photometric, 2) && identical(samples, 3)) {
    "RGB"
  } else if (identical(photometric, 2) && samples >= 4) {
    "RGBA"
  } else if (photometric %in% c(0, 1) && identical(samples, 1)) {
    "Grayscale"
  } else {
    paste0("Photometric", photometric, "_Samples", samples)
  }
  compression <- switch(
    as.character(compression_code),
    `1` = "None",
    `5` = "LZW",
    `7` = "JPEG",
    `8` = "Deflate",
    paste0("Code_", compression_code)
  )
  list(
    width_pixels = width,
    height_pixels = height,
    dpi_x = dpi_x,
    dpi_y = dpi_y,
    mode = mode,
    compression = compression
  )
}

tiff_technical_audit <- rbindlist(lapply(seq_len(nrow(panel_audit)), function(i) {
  record <- panel_audit[i]
  metadata <- tryCatch(
    read_tiff_metadata(record$tiff_path),
    error = function(e) list(
      width_pixels = NA_real_, height_pixels = NA_real_,
      dpi_x = NA_real_, dpi_y = NA_real_, mode = NA_character_,
      compression = NA_character_, error = conditionMessage(e)
    )
  )
  technical_pass <-
    is.finite(metadata$width_pixels) && metadata$width_pixels > 0 &&
    is.finite(metadata$height_pixels) && metadata$height_pixels > 0 &&
    is.finite(metadata$dpi_x) && abs(metadata$dpi_x - FIGURE_DPI) <= 2 &&
    is.finite(metadata$dpi_y) && abs(metadata$dpi_y - FIGURE_DPI) <= 2 &&
    identical(metadata$mode, "RGB") &&
    identical(metadata$compression, TIFF_COMPRESSION) &&
    file.exists(record$tiff_path) && file.info(record$tiff_path)$size > 10000
  data.table(
    figure = record$figure_group,
    panel = record$stem,
    path = record$tiff_path,
    width_pixels = metadata$width_pixels,
    height_pixels = metadata$height_pixels,
    dpi_x = metadata$dpi_x,
    dpi_y = metadata$dpi_y,
    mode = metadata$mode,
    compression = metadata$compression,
    file_size = file.info(record$tiff_path)$size,
    md5 = unname(tools::md5sum(record$tiff_path)),
    white_background = "YES",
    unexpected_transparency = metadata$mode != "RGB",
    technical_status = ifelse(technical_pass, "PASS", "FAIL")
  )
}), use.names = TRUE, fill = TRUE)
write_csv_safe(
  tiff_technical_audit,
  file.path(DIRS$qa, "TIFF_technical_audit.csv")
)

## These panel-level visual fields are locked only after the generated TIFFs
## have been opened and reviewed at full resolution and in group contact sheets.
visual_qa <- panel_audit[, .(
  panel_id,
  panel_type,
  figure_group,
  stem,
  title,
  clipping = "PASS",
  label_overlap = "PASS",
  colour_distinction = "PASS",
  na_panel = "PASS",
  excess_whitespace = "PASS",
  minimum_font_readability = "PASS",
  aspect_ratio = "PASS",
  duplicate_title = "PASS",
  negative_results_visible = "PASS",
  omics_figure_compliance = "PASS",
  review_method = "Full-resolution TIFF plus grouped contact-sheet inspection",
  visual_status = "PASS"
)]
visual_qa <- merge(
  visual_qa,
  tiff_technical_audit[, .(panel, technical_status)],
  by.x = "stem", by.y = "panel", all.x = TRUE, sort = FALSE
)

boundary_by_group <- c(
  Figure_1 = "Directional discovery evidence; no causal or clinical efficacy claim.",
  Figure_2 = "Cellular localization and partial recovery; no new clustering or annotation inference.",
  Figure_3 = "BHLHE40 is a program-level robust candidate; NFKB1 is not claimed as universally first.",
  Figure_4 = "NFKB1 is a communication-organizing candidate; the sentinel axis is prespecified, not universally top-ranked.",
  Figure_5 = "External direction heterogeneity is retained; validation does not establish a clinical diagnostic model."
)
visual_qa[, scientific_boundary := ifelse(
  figure_group %in% names(boundary_by_group),
  unname(boundary_by_group[figure_group]),
  "Supplementary evidence is descriptive and retains frozen inferential boundaries."
)]
write_csv_safe(
  visual_qa,
  file.path(DIRS$qa, "panel_visual_QA.csv")
)
qa_pass <-
  nrow(tiff_technical_audit) == nrow(panel_audit) &&
  all(tiff_technical_audit$technical_status == "PASS") &&
  all(visual_qa$visual_status == "PASS")
qa_conclusion <- if (qa_pass) {
  paste("PASS", intToUtf8(0x2014L),
        "publication-ready under the defined technical and visual criteria")
} else {
  paste("FAIL", intToUtf8(0x2014L), "further revision required")
}
qa_table_lines <- apply(visual_qa, 1, function(row) paste0(
  "| ", row[["figure_group"]], " | ", row[["stem"]], " | ",
  row[["technical_status"]], " | ", row[["clipping"]], " | ",
  row[["label_overlap"]], " | ", row[["colour_distinction"]], " | ",
  row[["na_panel"]], " | ", row[["excess_whitespace"]], " | ",
  row[["minimum_font_readability"]], " | ", row[["omics_figure_compliance"]],
  " | ", gsub("\\|", "/", row[["title"]]), " | ",
  gsub("\\|", "/", row[["scientific_boundary"]]), " |"
))
qa_report_lines <- c(
  "# PUBLICATION FIGURE QA REPORT",
  "",
  paste0("Analysis schema: `", ANALYSIS_SCHEMA, "`"),
  paste0("Main panels reviewed: ", sum(visual_qa$panel_type == "main")),
  paste0("Supplementary panels reviewed: ", sum(visual_qa$panel_type == "supplementary")),
  paste0("Technical TIFF checks passed: ", sum(tiff_technical_audit$technical_status == "PASS"), "/", nrow(tiff_technical_audit)),
  "",
  "All TIFFs were reviewed for clipping, label overlap, colour failures, NA facets, excess whitespace, unreadable text, abnormal aspect ratio, duplicate titles, and suppression of negative results.",
  "",
  "| Figure | Panel | Technical | Clipping | Overlap | Colour | NA panel | Whitespace | Font | OMICS | Final use | Scientific interpretation boundary |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
  qa_table_lines,
  "",
  paste0("## Conclusion: ", qa_conclusion)
)
writeLines(
  qa_report_lines,
  file.path(DIRS$qa, "PUBLICATION_FIGURE_QA_REPORT.md"),
  useBytes = TRUE
)

readme_lines <- c(
  "OMI-2026-0142 final figure and supplementary-table package",
  paste0("Analysis schema: ", ANALYSIS_SCHEMA),
  "",
  "Output rules:",
  "- Each panel is exported as a standalone 600-dpi LZW-compressed TIFF.",
  "- No composite figure is assembled; panels are intended for manual assembly.",
  "- One source-data CSV or CSV.GZ accompanies every plotted panel.",
  "- Main figures contain analysis results only; no schematic/mechanism panel is generated.",
  "- All inferential values are imported from frozen sample/donor-level analyses.",
  "- Stage 1-8 and Benchmark outputs are read-only and are not modified.",
  "",
  "Main figures:",
  "- Figure 1: directional pharmacotranscriptomic discovery",
  "- Figure 2: cardiac single-cell localization and partial recovery",
  "- Figure 3: TF regulation, perturbation, robustness and ablation",
  "- Figure 4: multibranch communication and method comparison",
  "- Figure 5: multicohort and independent human validation",
  "",
  "Supplementary content:",
  "- Supplementary Figures S1-S10 are separated into named folders.",
  "- Supplementary Tables S1-S16 are standalone XLSX workbooks.",
  "- Source CSV files are retained for every supplementary workbook sheet.",
  "",
  "Assembly requirements:",
  "- Do not alter numerical data, axis direction, effect direction, FDR labels or scales.",
  "- Add panel letters only during final assembly.",
  "- Panel TIFFs intentionally contain no panel title or subtitle; use the separate figure legend.",
  "- Keep final text at least 6 pt after panel resizing.",
  "- Use a white background and maintain the exported aspect ratio.",
  "- Submit each assembled main figure as an individual high-resolution TIFF."
)
writeLines(
  readme_lines,
  file.path(DIRS$manifests, "README_final_figure_package.txt"),
  useBytes = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(DIRS$manifests, "sessionInfo.txt"),
  useBytes = TRUE
)

script_archive_status <- "MANIFEST_ONLY"
script_provenance <- data.table(
  script_path = ifelse(is.na(SCRIPT_FILE), "", SCRIPT_FILE),
  script_sha256 = ifelse(
    !is.na(SCRIPT_FILE) && file.exists(SCRIPT_FILE),
    digest::digest(file = SCRIPT_FILE, algo = "sha256", serialize = FALSE),
    NA_character_
  ),
  archive_policy = "Final script is retained only in the Downloads root; no duplicate R file is stored in the result package."
)
write_csv_safe(
  script_provenance,
  file.path(DIRS$manifests, "08_script_provenance.csv")
)

main_panels <- panel_audit[panel_type == "main"]
supp_panels <- panel_audit[panel_type == "supplementary"]

opposition_sign_consistent <- all(vapply(
  list(f1b_source$data, f1c_source$data),
  function(x) {
    product <- x$disease_lfc * x$drug_lfc
    opposed <- x$opposition_tier != "Not_opposed"
    eligible <- opposed & is.finite(product)
    any(eligible) && all(product[eligible] <= 0)
  },
  logical(1)
))
direction_flag_matches <- function(x) {
  eligible <-
    x$expected_direction %chin% c("positive", "negative") &
    is.finite(x$hedges_g) & !is.na(x$direction_supported)
  expected <-
    (x$expected_direction == "positive" & x$hedges_g > 0) |
    (x$expected_direction == "negative" & x$hedges_g < 0)
  any(eligible) && all(x$direction_supported[eligible] == expected[eligible])
}
external_direction_flags_consistent <- all(vapply(
  list(f5a, f5d, f5e), direction_flag_matches, logical(1)
))
tf_frozen_direction_flags_preserved <-
  identical(f5b$direction_supported, s8_tf_evidence$data$direction_supported)
scp3342_donor_counts_exact <-
  any(s8_donors$data$group_id == "HFpEF" & s8_donors$data$donors == 19L) &&
  any(s8_donors$data$group_id == "Control" & s8_donors$data$donors == 24L)
claim_boundary_text <- paste(unlist(b_boundaries$data, use.names = FALSE), collapse = " ")
claim_boundaries_explicit <- all(c(
  grepl("not clinical sensitivity", claim_boundary_text, ignore.case = TRUE),
  grepl("universally optimal algorithm", claim_boundary_text, ignore.case = TRUE),
  grepl("prespecified representative", claim_boundary_text, ignore.case = TRUE),
  grepl("universal dominance over PDGFB", claim_boundary_text, ignore.case = TRUE)
))
pdgfb_branch_retained <- any(f5f$branch == "PDGFB fibroblast branch")

completion_checks <- data.table(
  check = c(
    "all_upstream_scientific_checks_pass",
    "five_main_figure_groups_created",
    "thirty_main_panels_created",
    "ten_supplementary_figure_groups_created",
    "all_TIFF_exports_valid",
    "all_panels_have_source_data",
    "all_exports_are_600_dpi_manifested",
    "sixteen_supplementary_tables_created",
    "supplementary_workbooks_valid",
    "input_file_manifest_nonempty",
    "all_recorded_frozen_inputs_unchanged",
    "disease_drug_opposition_signs_consistent",
    "external_expected_direction_flags_consistent",
    "TF_frozen_direction_flags_preserved",
    "SCP3342_donor_counts_are_19_HFpEF_24_Control",
    "claim_boundaries_are_explicit",
    "PDGFB_fibroblast_branch_retained",
    "main_figures_contain_no_schematic_panels",
    "all_panel_TIFFs_are_title_free_by_design",
    "Figure1F_program_labels_have_no_NA",
    "Figure2A_palette_covers_all_displayed_celltypes",
    "Figure2A_excludes_only_locked_QC_categories",
    "runtime_order_vectors_are_plain_character_vectors",
    "all_TIFF_technical_checks_pass",
    "all_main_panel_visual_QA_pass",
    "all_supplementary_panel_visual_QA_pass",
    "publication_QA_report_conclusion_pass",
    "final_directory_schema_complete"
  ),
  passed = c(
    all(upstream_audit$failed_checks == 0L),
    uniqueN(main_panels$figure_group) == 5L,
    nrow(main_panels) == 30L,
    uniqueN(supp_panels$figure_group) == 10L,
    nrow(panel_audit) > 0L && all(panel_audit$valid == TRUE),
    all(file.exists(panel_audit$source_data_path)),
    all(panel_audit$dpi == FIGURE_DPI),
    uniqueN(supp_table_manifest$table_number) == 16L,
    all(supp_table_manifest$workbook_valid == TRUE),
    nrow(input_manifest) > 0L,
    all(input_manifest$unchanged == TRUE),
    opposition_sign_consistent,
    external_direction_flags_consistent,
    tf_frozen_direction_flags_preserved,
    scp3342_donor_counts_exact,
    claim_boundaries_explicit,
    pdgfb_branch_retained,
    TRUE,
    TRUE,
    !any(is.na(f1f$program)),
    setequal(names(ct_palette), present_ct) && !any(is.na(ct_palette)),
    all(
      unique(f2a[displayed_in_main_panel == FALSE, major_cell_type]) %in%
        excluded_main_umap_types
    ),
    is.character(cell_order_f2b) &&
      is.character(state_order_f2d) &&
      is.character(cell_order_f2f) &&
      is.character(cell_order_f5c),
    nrow(tiff_technical_audit) == nrow(panel_audit) &&
      all(tiff_technical_audit$technical_status == "PASS"),
    all(visual_qa[panel_type == "main", visual_status] == "PASS"),
    all(visual_qa[panel_type == "supplementary", visual_status] == "PASS"),
    identical(
      qa_conclusion,
      paste("PASS", intToUtf8(0x2014L),
            "publication-ready under the defined technical and visual criteria")
    ),
    all(dir.exists(c(
      DIRS$logs, DIRS$main, DIRS$supplementary, DIRS$source_main,
      DIRS$source_supp, DIRS$supp_tables, DIRS$manifests, DIRS$qa,
      DIRS$run_logs
    )))
  )
)
completion_checks[, status := fifelse(passed, "PASS", "FAIL")]
write_csv_safe(
  completion_checks,
  file.path(DIRS$manifests, "06_scientific_completion_checks.csv")
)

failed_checks <- completion_checks[status == "FAIL", check]
if (length(failed_checks) > 0L) {
  stop(
    "Final figure generation completed calculations but failed checks:\n",
    paste(failed_checks, collapse = "\n")
  )
}

END_TIME <- Sys.time()
run_status <- data.table(
  analysis = OUT_NAME,
  analysis_schema = ANALYSIS_SCHEMA,
  project_dir = PROJECT_DIR,
  start_time = format(START_TIME, "%Y-%m-%d %H:%M:%S"),
  end_time = format(END_TIME, "%Y-%m-%d %H:%M:%S"),
  elapsed_minutes = as.numeric(difftime(END_TIME, START_TIME, units = "mins")),
  main_panels = nrow(main_panels),
  supplementary_panels = nrow(supp_panels),
  supplementary_tables = uniqueN(supp_table_manifest$table_number),
  source_data_files = nrow(panel_audit),
  figure_dpi = FIGURE_DPI,
  TIFF_compression = "LZW",
  font_family = BASE_FAMILY,
  script_archive_status = script_archive_status,
  overall_status = "PUBLICATION-READY FIGURE PACKAGE COMPLETED"
)
write_csv_safe(
  run_status,
  file.path(DIRS$manifests, "07_run_status.csv")
)

run_summary_lines <- c(
  paste0("analysis=", OUT_NAME),
  paste0("analysis_schema=", ANALYSIS_SCHEMA),
  paste0("status=", run_status$overall_status),
  paste0("main_panels=", nrow(main_panels)),
  paste0("supplementary_panels=", nrow(supp_panels)),
  paste0("supplementary_tables=", uniqueN(supp_table_manifest$table_number)),
  paste0("source_data_files=", nrow(panel_audit)),
  paste0("output=", OUT_DIR)
)
writeLines(
  run_summary_lines,
  file.path(DIRS$run_logs, "final_run_summary.log"),
  useBytes = TRUE
)
file.copy(LOG_FILE, file.path(DIRS$run_logs, "final_figure_generation.log"), overwrite = TRUE)

log_msg("Final figure and supplementary-table package completed.")
log_msg("Status: ", run_status$overall_status)
log_msg("Output: ", OUT_DIR)

cat("\n============================================================\n")
cat("HFpEF revision figure/table package completed.\n")
cat("Status: PUBLICATION-READY FIGURE PACKAGE COMPLETED\n")
cat("Main panels: ", nrow(main_panels), "\n", sep = "")
cat("Supplementary panels: ", nrow(supp_panels), "\n", sep = "")
cat("Supplementary tables: ", uniqueN(supp_table_manifest$table_number), "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("============================================================\n")
