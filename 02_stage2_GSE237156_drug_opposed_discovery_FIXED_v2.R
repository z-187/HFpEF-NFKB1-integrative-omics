############################################################
## Public repository configuration
##   Set the project root before execution, for example:
##   Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
##   Raw public data are not bundled with this repository.
##
## HFpEF Reanalysis Project
## Stage 2 FIXED v2 REPLACE
## GSE237156 drug-opposed macrophage transcriptomic discovery (clean replacement run)
##
## FIXED v2 changes:
##   - Deletes the incomplete FIXED_v1 output directory and CHECK archive.
##   - Fixes the GSEA ordering error caused by setorder(..., -abs(stat)).
##   - Uses an explicit abs_stat column before ordering ranked genes.
##   - Prefers mouse-native MSigDB MM/MH Hallmark gene sets.
##   - Makes optional GSEA nonfatal to the core transcriptomic analysis.
##   - Uses DESeq2 normal-prior shrinkage without optional ashr compilation.
##
## Project root:
##   <HFPEF_PROJECT_DIR>
##
## Read-only input:
##   <HFPEF_PROJECT_DIR>/0.GEO/GSE237156_RAW.tar
##
## Locked metadata:
##   <HFPEF_PROJECT_DIR>/
##   01_stage1_metadata_lock_FIXED_v3/01_tables/
##   01_locked_sample_manifest.csv
##
## Output:
##   <HFPEF_PROJECT_DIR>/
##   02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2
##
## CHECK package:
##   <HFPEF_PROJECT_DIR>/
##   02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2_CHECK.zip
##
## Primary objectives:
##   1) Reconstruct the 16-sample RSEM expected-count matrix.
##   2) Validate the locked 2 x 2 x 2 design:
##      Ccr2 status x diet x dapagliflozin.
##   3) Perform DESeq2 primary differential-expression analysis.
##   4) Perform edgeR quasi-likelihood sensitivity analysis.
##   5) Rank all genes by continuous drug-opposition evidence without
##      preselecting Nfkb1 or a fixed candidate-gene list.
##   6) Quantify method agreement and cross-subset reproducibility.
##   7) Run transparent factorial LRT sensitivity analyses.
##   8) Optionally evaluate Hallmark pathway opposition by GSEA.
##
## Scientific boundary:
##   - n = 2 biological samples per experimental cell.
##   - Results are hypothesis-generating and effect-size oriented.
##   - Drug-opposed expression does not prove direct drug-target binding,
##     transcription-factor causality, or therapeutic efficacy.
##   - No candidate gene is forced into the ranked results.
##
## Recommended run:
##   source(
##     "<HFPEF_PROJECT_DIR>/
##      HFpEF_Stage2_GSE237156_Drug_Opposed_Discovery_FIXED_v2_REPLACE.R",
##     encoding = "UTF-8"
##   )
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(warn = 1)
options(encoding = "UTF-8")
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

STAGE1_LOCK_DIR <- file.path(
  PROJECT_DIR,
  "01_stage1_metadata_lock_FIXED_v3"
)
LOCKED_MANIFEST_FILE <- file.path(
  STAGE1_LOCK_DIR,
  "01_tables",
  "01_locked_sample_manifest.csv"
)
RAW_TAR_FILE <- file.path(DATA_DIR, "GSE237156_RAW.tar")

## The prior incomplete Stage 2 output is intentionally deleted and replaced.
OLD_STAGE_NAME <- "02_stage2_GSE237156_drug_opposed_discovery_FIXED_v1"
OLD_OUT_DIR <- file.path(PROJECT_DIR, OLD_STAGE_NAME)
OLD_CHECK_ZIP <- file.path(
  PROJECT_DIR,
  paste0(OLD_STAGE_NAME, "_CHECK.zip")
)

STAGE_NAME <- "02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2"
OUT_DIR <- file.path(PROJECT_DIR, STAGE_NAME)
CHECK_ZIP <- file.path(
  PROJECT_DIR,
  paste0(STAGE_NAME, "_CHECK.zip")
)

REPLACE_PREVIOUS_STAGE2 <- TRUE
RUN_OPTIONAL_GSEA <- TRUE

## Filtering and reporting parameters.
MIN_COUNT <- 10L
MIN_SAMPLES_WITH_MIN_COUNT <- 2L
DESEQ_ALPHA <- 0.10
FORMAL_FDR <- 0.05
EXPLORATORY_FDR <- 0.10
MIN_ABS_LFC_FOR_TIER <- 0.50
TOP_HEATMAP_GENES <- 40L
TOP_LABEL_GENES <- 12L
TOP_CHECK_GENES <- 500L

## The script snapshot is useful but is not a scientific pass criterion.
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

  if (length(candidates) == 0L) return(NA_character_)

  normalizePath(
    candidates[1L],
    winslash = "/",
    mustWork = TRUE
  )
}
SCRIPT_FILE <- detect_script_file()

############################################################
## 1. Output setup and logging
############################################################

required_inputs <- c(
  PROJECT_DIR,
  DATA_DIR,
  LOCKED_MANIFEST_FILE,
  RAW_TAR_FILE
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Required input path(s) are missing:\n",
    paste(missing_inputs, collapse = "\n")
  )
}

replacement_records <- data.frame(
  path = c(
    OLD_OUT_DIR,
    OLD_CHECK_ZIP,
    OUT_DIR,
    CHECK_ZIP
  ),
  path_type = c(
    "prior_incomplete_output_directory",
    "prior_incomplete_check_zip",
    "current_v2_output_directory",
    "current_v2_check_zip"
  ),
  existed_before = FALSE,
  deletion_attempted = FALSE,
  deletion_succeeded = FALSE,
  stringsAsFactors = FALSE
)

if (REPLACE_PREVIOUS_STAGE2) {
  for (i in seq_len(nrow(replacement_records))) {
    target <- replacement_records$path[i]
    existed <- dir.exists(target) || file.exists(target)
    replacement_records$existed_before[i] <- existed

    if (existed) {
      replacement_records$deletion_attempted[i] <- TRUE
      unlink(
        target,
        recursive = dir.exists(target),
        force = TRUE
      )
      replacement_records$deletion_succeeded[i] <- !(
        dir.exists(target) || file.exists(target)
      )
      if (!replacement_records$deletion_succeeded[i]) {
        stop(
          "Failed to delete previous Stage 2 path:\n",
          target
        )
      }
    } else {
      replacement_records$deletion_succeeded[i] <- TRUE
    }
  }
} else {
  if (
    dir.exists(OLD_OUT_DIR) ||
    file.exists(OLD_CHECK_ZIP) ||
    dir.exists(OUT_DIR) ||
    file.exists(CHECK_ZIP)
  ) {
    stop(
      "Existing Stage 2 output was detected while replacement is disabled."
    )
  }
}

DIRS <- list(
  logs = file.path(OUT_DIR, "00_logs"),
  tables = file.path(OUT_DIR, "01_tables"),
  objects = file.path(OUT_DIR, "02_objects"),
  figures = file.path(OUT_DIR, "03_figures"),
  source = file.path(OUT_DIR, "04_source_data"),
  methods = file.path(OUT_DIR, "05_methods"),
  check = file.path(OUT_DIR, "06_review_check"),
  extracted = file.path(OUT_DIR, "07_extracted_RSEM")
)
for (d in c(OUT_DIR, unlist(DIRS, use.names = FALSE))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(DIRS$logs, "stage2_GSE237156.log")
WARN_FILE <- file.path(DIRS$logs, "stage2_warnings.log")
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

log_msg("Stage 2 GSE237156 analysis started.")
log_msg(
  "Previous incomplete Stage 2 output replacement enabled: ",
  REPLACE_PREVIOUS_STAGE2
)
for (i in seq_len(nrow(replacement_records))) {
  log_msg(
    "Replacement audit | ",
    replacement_records$path_type[i],
    " | existed_before=",
    replacement_records$existed_before[i],
    " | deletion_succeeded=",
    replacement_records$deletion_succeeded[i],
    " | path=",
    replacement_records$path[i]
  )
}
write.csv(
  replacement_records,
  file.path(
    DIRS$logs,
    "previous_stage2_replacement_audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
log_msg("RAW_TAR_FILE: ", RAW_TAR_FILE)
log_msg("LOCKED_MANIFEST_FILE: ", LOCKED_MANIFEST_FILE)
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
## 2. Package installation and loading
############################################################

ensure_cran <- function(pkgs, required = TRUE) {
  missing <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
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
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
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
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
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
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
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
    "data.table",
    "ggplot2",
    "openxlsx",
    "pheatmap",
    "ggrepel",
    "zip",
    "digest",
    "matrixStats"
  ),
  required = TRUE
)
ensure_bioc(
  c(
    "DESeq2",
    "edgeR",
    "limma",
    "AnnotationDbi",
    "org.Mm.eg.db"
  ),
  required = TRUE
)

## DESeq2 normal-prior shrinkage is used for reproducibility.
## Optional GSEA packages are nonfatal.
if (RUN_OPTIONAL_GSEA) {
  ensure_cran("msigdbr", required = FALSE)
  ensure_bioc("fgsea", required = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(openxlsx)
  library(pheatmap)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
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

write_csv_safe <- function(x, path, compress = FALSE) {
  if (is.null(x) || ncol(x) == 0L) {
    fwrite(
      data.table(note = "No records generated."),
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
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
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

read_rsem_gene_results <- function(path, sample_accession) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  x <- tryCatch(
    read.delim(
      con,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      stop(
        "Failed to read RSEM file ",
        basename(path),
        ": ",
        conditionMessage(e)
      )
    }
  )
  setDT(x)

  original_names <- names(x)
  normalized_names <- tolower(
    gsub(
      "[^A-Za-z0-9]+",
      "_",
      original_names
    )
  )
  setnames(x, original_names, normalized_names)

  required <- c("gene_id", "expected_count")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      "RSEM file ",
      basename(path),
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  if (!"tpm" %in% names(x)) {
    x[, tpm := NA_real_]
  }
  if (!"fpkm" %in% names(x)) {
    x[, fpkm := NA_real_]
  }

  x[, gene_id_raw := as.character(gene_id)]
  x[, gene_id := sub("\\.[0-9]+$", "", gene_id_raw)]
  x[, expected_count := suppressWarnings(
    as.numeric(expected_count)
  )]
  x[, tpm := suppressWarnings(as.numeric(tpm))]
  x[, fpkm := suppressWarnings(as.numeric(fpkm))]

  if (any(!is.finite(x$expected_count))) {
    stop(
      "Non-finite expected_count values were detected in ",
      basename(path)
    )
  }

  ## Aggregate only if version stripping creates duplicate identifiers.
  out <- x[, .(
    expected_count = sum(
      expected_count,
      na.rm = TRUE
    ),
    tpm = sum(tpm, na.rm = TRUE),
    fpkm = sum(fpkm, na.rm = TRUE),
    collapsed_source_rows = .N
  ), by = gene_id]

  out[, sample_accession := sample_accession]
  out
}

safe_rank <- function(x, decreasing = TRUE) {
  x2 <- x
  x2[!is.finite(x2)] <- NA_real_
  rank(
    if (decreasing) -x2 else x2,
    na.last = "keep",
    ties.method = "average"
  )
}

safe_neglog10 <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 1
  -log10(pmax(x, .Machine$double.xmin))
}

get_gene_annotation <- function(gene_ids) {
  gene_ids <- unique(as.character(gene_ids))
  gene_ids <- gene_ids[nzchar(gene_ids)]

  ensembl_fraction <- mean(
    grepl("^ENSMUSG", gene_ids),
    na.rm = TRUE
  )
  entrez_fraction <- mean(
    grepl("^[0-9]+$", gene_ids),
    na.rm = TRUE
  )

  if (
    is.finite(ensembl_fraction) &&
    ensembl_fraction >= 0.50
  ) {
    keytype <- "ENSEMBL"
  } else if (
    is.finite(entrez_fraction) &&
    entrez_fraction >= 0.50
  ) {
    keytype <- "ENTREZID"
  } else {
    keytype <- "SYMBOL"
  }

  ann <- tryCatch(
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = gene_ids,
      keytype = keytype,
      columns = unique(
        c(keytype, "SYMBOL", "ENTREZID", "GENENAME")
      )
    ),
    error = function(e) {
      add_warning(
        "ANNOTATION",
        keytype,
        conditionMessage(e)
      )
      data.frame()
    }
  )
  setDT(ann)

  if (nrow(ann) == 0L) {
    return(
      data.table(
        gene_id = gene_ids,
        symbol = NA_character_,
        entrez_id = NA_character_,
        gene_name = NA_character_,
        annotation_keytype = keytype
      )
    )
  }

  id_col <- keytype
  setnames(
    ann,
    old = intersect(
      c(id_col, "SYMBOL", "ENTREZID", "GENENAME"),
      names(ann)
    ),
    new = c(
      "gene_id",
      "symbol",
      "entrez_id",
      "gene_name"
    )[
      match(
        intersect(
          c(id_col, "SYMBOL", "ENTREZID", "GENENAME"),
          names(ann)
        ),
        c(id_col, "SYMBOL", "ENTREZID", "GENENAME")
      )
    ]
  )

  for (nm in c(
    "gene_id",
    "symbol",
    "entrez_id",
    "gene_name"
  )) {
    if (!nm %in% names(ann)) ann[, (nm) := NA_character_]
  }

  ann <- ann[
    !is.na(gene_id) &
      nzchar(gene_id)
  ]
  ann[, symbol_present := (
    !is.na(symbol) &
      nzchar(symbol)
  )]
  setorder(ann, gene_id, -symbol_present)
  ann <- ann[, .SD[1L], by = gene_id]
  ann[, symbol_present := NULL]
  ann[, annotation_keytype := keytype]

  merge(
    data.table(gene_id = gene_ids),
    ann,
    by = "gene_id",
    all.x = TRUE
  )
}

make_deseq_result <- function(
  dds,
  contrast_name,
  numerator,
  denominator,
  alpha = DESEQ_ALPHA
) {
  res_raw <- DESeq2::results(
    dds,
    contrast = c(
      "group",
      numerator,
      denominator
    ),
    alpha = alpha,
    independentFiltering = TRUE
  )

  ## Use DESeq2's normal-prior shrinkage to avoid optional compiled
  ## dependencies and to keep the rerun stable on Windows.
  shrink_type <- "normal"
  res_shrunk <- tryCatch(
    DESeq2::lfcShrink(
      dds,
      contrast = c(
        "group",
        numerator,
        denominator
      ),
      res = res_raw,
      type = "normal",
      quiet = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(res_shrunk)) {
    res_shrunk <- res_raw
    shrink_type <- "unshrunken_fallback"
    add_warning(
      "LFC_SHRINK",
      contrast_name,
      "Normal-prior shrinkage failed; unshrunken DESeq2 log2FC retained."
    )
  }

  out <- data.table(
    gene_id = rownames(res_raw),
    base_mean = res_raw$baseMean,
    log2fc_raw = res_raw$log2FoldChange,
    lfc_se_raw = res_raw$lfcSE,
    stat = res_raw$stat,
    pvalue = res_raw$pvalue,
    padj = res_raw$padj,
    log2fc_shrunk = res_shrunk$log2FoldChange,
    lfc_se_shrunk = res_shrunk$lfcSE,
    contrast = contrast_name,
    numerator = numerator,
    denominator = denominator,
    lfc_shrink_type = shrink_type,
    method = "DESeq2"
  )

  out[, formal_fdr_005 := (
    !is.na(padj) &
      padj <= FORMAL_FDR
  )]
  out[, exploratory_fdr_010 := (
    !is.na(padj) &
      padj <= EXPLORATORY_FDR
  )]
  out
}

make_edger_result <- function(
  fit,
  contrast_name,
  numerator,
  denominator,
  design_columns
) {
  contrast_vector <- rep(
    0,
    length(design_columns)
  )
  names(contrast_vector) <- design_columns

  if (
    !numerator %in% design_columns ||
    !denominator %in% design_columns
  ) {
    stop(
      "edgeR contrast level not found: ",
      numerator,
      " vs ",
      denominator
    )
  }

  contrast_vector[numerator] <- 1
  contrast_vector[denominator] <- -1

  test <- edgeR::glmQLFTest(
    fit,
    contrast = contrast_vector
  )
  tab <- edgeR::topTags(
    test,
    n = Inf,
    sort.by = "none"
  )$table
  setDT(tab, keep.rownames = "gene_id")

  out <- tab[, .(
    gene_id,
    log2fc = logFC,
    log_cpm = logCPM,
    qlf = F,
    pvalue = PValue,
    padj = FDR
  )]
  out[, contrast := contrast_name]
  out[, numerator := numerator]
  out[, denominator := denominator]
  out[, method := "edgeR_QL"]
  out[, formal_fdr_005 := (
    !is.na(padj) &
      padj <= FORMAL_FDR
  )]
  out[, exploratory_fdr_010 := (
    !is.na(padj) &
      padj <= EXPLORATORY_FDR
  )]
  out
}

summarize_contrast <- function(x, method_name) {
  data.table(
    contrast = unique(x$contrast)[1L],
    method = method_name,
    tested_genes = nrow(x),
    fdr_005_up = sum(
      x$padj <= FORMAL_FDR &
        if (
          "log2fc_shrunk" %in% names(x)
        ) {
          x$log2fc_shrunk > 0
        } else {
          x$log2fc > 0
        },
      na.rm = TRUE
    ),
    fdr_005_down = sum(
      x$padj <= FORMAL_FDR &
        if (
          "log2fc_shrunk" %in% names(x)
        ) {
          x$log2fc_shrunk < 0
        } else {
          x$log2fc < 0
        },
      na.rm = TRUE
    ),
    fdr_010_up = sum(
      x$padj <= EXPLORATORY_FDR &
        if (
          "log2fc_shrunk" %in% names(x)
        ) {
          x$log2fc_shrunk > 0
        } else {
          x$log2fc > 0
        },
      na.rm = TRUE
    ),
    fdr_010_down = sum(
      x$padj <= EXPLORATORY_FDR &
        if (
          "log2fc_shrunk" %in% names(x)
        ) {
          x$log2fc_shrunk < 0
        } else {
          x$log2fc < 0
        },
      na.rm = TRUE
    )
  )
}

############################################################
## 4. Read and validate locked GSE237156 metadata
############################################################

locked_manifest <- fread(
  LOCKED_MANIFEST_FILE,
  encoding = "UTF-8",
  na.strings = c("", "NA", "NaN")
)

required_manifest_columns <- c(
  "dataset_id",
  "sample_accession",
  "original_title",
  "group_id",
  "condition",
  "diet_or_model",
  "drug",
  "macrophage_subset",
  "lock_status"
)
missing_manifest_columns <- setdiff(
  required_manifest_columns,
  names(locked_manifest)
)
if (length(missing_manifest_columns) > 0L) {
  stop(
    "Locked manifest is missing column(s): ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}

sample_meta <- locked_manifest[
  dataset_id == "GSE237156"
]
if (nrow(sample_meta) != 16L) {
  stop(
    "Expected 16 locked GSE237156 samples, found ",
    nrow(sample_meta),
    "."
  )
}
if (uniqueN(sample_meta$sample_accession) != 16L) {
  stop("GSE237156 sample accessions are not unique.")
}
if (any(!grepl("^LOCKED", sample_meta$lock_status))) {
  stop("At least one GSE237156 sample is not metadata-locked.")
}

sample_meta[, macrophage_subset := factor(
  macrophage_subset,
  levels = c(
    "Ccr2_negative",
    "Ccr2_positive"
  )
)]
sample_meta[, diet := factor(
  diet_or_model,
  levels = c("CD", "HFD")
)]
sample_meta[, drug_factor := factor(
  drug,
  levels = c(
    "Vehicle",
    "Dapagliflozin"
  )
)]
sample_meta[, group := factor(group_id)]

group_check <- sample_meta[, .N, by = group_id]
if (
  nrow(group_check) != 8L ||
  any(group_check$N != 2L)
) {
  fwrite(
    group_check,
    file.path(
      DIRS$tables,
      "FATAL_GSE237156_group_count_check.csv"
    )
  )
  stop(
    "The locked GSE237156 design is not 8 groups x 2 samples."
  )
}

setorder(
  sample_meta,
  macrophage_subset,
  diet,
  drug_factor,
  sample_accession
)
fwrite(
  sample_meta,
  file.path(
    DIRS$tables,
    "01_sample_metadata_used.csv"
  )
)
fwrite(
  group_check,
  file.path(
    DIRS$tables,
    "02_group_count_validation.csv"
  )
)
log_msg(
  "Validated GSE237156 design: 16 samples, 8 groups, 2 samples per group."
)

############################################################
## 5. Map and extract RSEM files from GEO RAW tar
############################################################

tar_members <- tryCatch(
  utils::untar(
    RAW_TAR_FILE,
    list = TRUE
  ),
  error = function(e) {
    stop(
      "Unable to list GSE237156_RAW.tar: ",
      conditionMessage(e)
    )
  }
)

if (length(tar_members) == 0L) {
  stop("GSE237156_RAW.tar contains no listed members.")
}

mapping_records <- list()
for (gsm in sample_meta$sample_accession) {
  hits <- tar_members[
    grepl(
      gsm,
      tar_members,
      fixed = TRUE
    ) &
      grepl(
        "\\.genes\\.results(\\.gz)?$",
        tar_members,
        ignore.case = TRUE
      )
  ]

  if (length(hits) != 1L) {
    mapping_records[[length(mapping_records) + 1L]] <- data.table(
      sample_accession = gsm,
      archive_member = paste(hits, collapse = "; "),
      match_count = length(hits),
      mapping_status = "FAIL"
    )
  } else {
    mapping_records[[length(mapping_records) + 1L]] <- data.table(
      sample_accession = gsm,
      archive_member = hits,
      match_count = 1L,
      mapping_status = "PASS"
    )
  }
}

sample_file_map <- rbindlist(mapping_records)
sample_file_map <- merge(
  sample_meta[, .(
    sample_accession,
    original_title,
    group_id,
    macrophage_subset,
    diet,
    drug = drug_factor
  )],
  sample_file_map,
  by = "sample_accession",
  all.x = TRUE
)

fwrite(
  sample_file_map,
  file.path(
    DIRS$tables,
    "03_sample_to_RSEM_file_mapping.csv"
  )
)

if (any(sample_file_map$mapping_status != "PASS")) {
  stop(
    "At least one GSM accession did not map uniquely to an RSEM gene-results file. ",
    "See 03_sample_to_RSEM_file_mapping.csv."
  )
}

members_to_extract <- unique(
  sample_file_map$archive_member
)
log_msg(
  "Extracting ",
  length(members_to_extract),
  " RSEM files."
)

utils::untar(
  RAW_TAR_FILE,
  files = members_to_extract,
  exdir = DIRS$extracted
)

sample_file_map[, local_file := file.path(
  DIRS$extracted,
  archive_member
)]
sample_file_map[, extracted_exists := file.exists(local_file)]

if (any(!sample_file_map$extracted_exists)) {
  stop("One or more selected RSEM files were not extracted successfully.")
}

fwrite(
  sample_file_map,
  file.path(
    DIRS$tables,
    "03_sample_to_RSEM_file_mapping.csv"
  )
)

############################################################
## 6. Build expected-count, TPM, and FPKM matrices
############################################################

rsem_list <- list()
for (i in seq_len(nrow(sample_file_map))) {
  gsm <- sample_file_map$sample_accession[i]
  f <- sample_file_map$local_file[i]
  log_msg(
    "Reading RSEM file ",
    i,
    "/",
    nrow(sample_file_map),
    ": ",
    basename(f)
  )
  rsem_list[[gsm]] <- read_rsem_gene_results(
    f,
    gsm
  )
}

rsem_long <- rbindlist(
  rsem_list,
  use.names = TRUE,
  fill = TRUE
)

count_wide <- dcast(
  rsem_long,
  gene_id ~ sample_accession,
  value.var = "expected_count",
  fill = 0
)
tpm_wide <- dcast(
  rsem_long,
  gene_id ~ sample_accession,
  value.var = "tpm",
  fill = 0
)
fpkm_wide <- dcast(
  rsem_long,
  gene_id ~ sample_accession,
  value.var = "fpkm",
  fill = 0
)

sample_order <- sample_meta$sample_accession
missing_count_columns <- setdiff(
  sample_order,
  names(count_wide)
)
if (length(missing_count_columns) > 0L) {
  stop(
    "Count matrix is missing sample column(s): ",
    paste(missing_count_columns, collapse = ", ")
  )
}

count_matrix_fractional <- as.matrix(
  count_wide[, ..sample_order]
)
rownames(count_matrix_fractional) <- count_wide$gene_id
storage.mode(count_matrix_fractional) <- "numeric"

tpm_matrix <- as.matrix(
  tpm_wide[, ..sample_order]
)
rownames(tpm_matrix) <- tpm_wide$gene_id
storage.mode(tpm_matrix) <- "numeric"

fpkm_matrix <- as.matrix(
  fpkm_wide[, ..sample_order]
)
rownames(fpkm_matrix) <- fpkm_wide$gene_id
storage.mode(fpkm_matrix) <- "numeric"

if (
  any(!is.finite(count_matrix_fractional)) ||
  any(count_matrix_fractional < 0)
) {
  stop("Count matrix contains non-finite or negative values.")
}

count_matrix <- round(count_matrix_fractional)
storage.mode(count_matrix) <- "integer"

gene_filter <- (
  rowSums(count_matrix >= MIN_COUNT) >=
    MIN_SAMPLES_WITH_MIN_COUNT
)
if (sum(gene_filter) < 5000L) {
  stop(
    "Only ",
    sum(gene_filter),
    " genes passed count filtering; inspect input reconstruction."
  )
}
count_matrix_filtered <- count_matrix[
  gene_filter,
  ,
  drop = FALSE
]

gene_filter_table <- data.table(
  gene_id = rownames(count_matrix),
  total_count = rowSums(count_matrix),
  samples_count_ge_min = rowSums(
    count_matrix >= MIN_COUNT
  ),
  retained = gene_filter
)
fwrite(
  gene_filter_table,
  file.path(
    DIRS$tables,
    "04_gene_filter_status.csv.gz"
  ),
  compress = "gzip"
)

saveRDS(
  count_matrix_fractional,
  file.path(
    DIRS$objects,
    "GSE237156_fractional_expected_counts.rds"
  )
)
saveRDS(
  count_matrix,
  file.path(
    DIRS$objects,
    "GSE237156_rounded_expected_counts.rds"
  )
)
saveRDS(
  count_matrix_filtered,
  file.path(
    DIRS$objects,
    "GSE237156_filtered_integer_counts.rds"
  )
)
saveRDS(
  tpm_matrix,
  file.path(
    DIRS$objects,
    "GSE237156_TPM_matrix.rds"
  )
)
saveRDS(
  fpkm_matrix,
  file.path(
    DIRS$objects,
    "GSE237156_FPKM_matrix.rds"
  )
)

log_msg(
  "Reconstructed ",
  nrow(count_matrix),
  " genes x ",
  ncol(count_matrix),
  " samples; ",
  sum(gene_filter),
  " genes retained."
)

############################################################
## 7. Gene annotation
############################################################

gene_annotation <- get_gene_annotation(
  rownames(count_matrix)
)
gene_annotation[, display_gene := fifelse(
  !is.na(symbol) &
    nzchar(symbol),
  symbol,
  gene_id
)]
fwrite(
  gene_annotation,
  file.path(
    DIRS$tables,
    "05_gene_annotation.csv"
  )
)

############################################################
## 8. Sample-level quality control
############################################################

sample_qc <- data.table(
  sample_accession = colnames(count_matrix),
  library_size_fractional = colSums(
    count_matrix_fractional
  ),
  library_size_rounded = colSums(count_matrix),
  detected_count_ge_1 = colSums(
    count_matrix >= 1L
  ),
  detected_count_ge_10 = colSums(
    count_matrix >= 10L
  ),
  detected_tpm_gt_0 = colSums(
    tpm_matrix > 0
  )
)
sample_qc <- merge(
  sample_meta[, .(
    sample_accession,
    original_title,
    group_id,
    macrophage_subset,
    diet,
    drug = drug_factor
  )],
  sample_qc,
  by = "sample_accession",
  all.x = TRUE
)
setorder(
  sample_qc,
  macrophage_subset,
  diet,
  drug,
  sample_accession
)
fwrite(
  sample_qc,
  file.path(
    DIRS$tables,
    "06_sample_QC_metrics.csv"
  )
)

coldata <- as.data.frame(
  sample_meta[, .(
    sample_accession,
    macrophage_subset,
    diet,
    drug = drug_factor,
    group
  )]
)
rownames(coldata) <- coldata$sample_accession
coldata <- coldata[
  colnames(count_matrix_filtered),
  ,
  drop = FALSE
]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = coldata,
  design = ~ group
)
dds <- DESeq(
  dds,
  minReplicatesForReplace = Inf,
  quiet = TRUE
)

vsd <- varianceStabilizingTransformation(
  dds,
  blind = FALSE
)
vst_matrix <- assay(vsd)
saveRDS(
  dds,
  file.path(
    DIRS$objects,
    "GSE237156_DESeq2_dds.rds"
  )
)
saveRDS(
  vsd,
  file.path(
    DIRS$objects,
    "GSE237156_DESeq2_vsd.rds"
  )
)

## PCA.
variable_genes <- head(
  order(
    matrixStats::rowVars(vst_matrix),
    decreasing = TRUE
  ),
  min(1000L, nrow(vst_matrix))
)
pca_fit <- prcomp(
  t(vst_matrix[variable_genes, , drop = FALSE]),
  center = TRUE,
  scale. = FALSE
)
pca_percent <- (
  100 *
    pca_fit$sdev^2 /
    sum(pca_fit$sdev^2)
)
pca_dt <- data.table(
  sample_accession = rownames(pca_fit$x),
  PC1 = pca_fit$x[, 1L],
  PC2 = pca_fit$x[, 2L]
)
pca_dt <- merge(
  sample_meta[, .(
    sample_accession,
    original_title,
    group_id,
    macrophage_subset,
    diet,
    drug = drug_factor
  )],
  pca_dt,
  by = "sample_accession",
  all.x = TRUE
)
fwrite(
  pca_dt,
  file.path(
    DIRS$source,
    "Fig2A_PCA_source.csv"
  )
)

p_pca <- ggplot(
  pca_dt,
  aes(
    x = PC1,
    y = PC2,
    shape = macrophage_subset
  )
) +
  geom_point(
    aes(
      fill = interaction(diet, drug)
    ),
    size = 3.5
  ) +
  geom_text_repel(
    aes(label = sample_accession),
    size = 2.7,
    max.overlaps = Inf
  ) +
  labs(
    x = sprintf(
      "PC1 (%.1f%%)",
      pca_percent[1L]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      pca_percent[2L]
    ),
    fill = "Diet × drug",
    shape = "Macrophage subset",
    title = "GSE237156 sample-level variance-stabilized PCA"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(
    DIRS$figures,
    "Fig2A_GSE237156_PCA.png"
  ),
  p_pca,
  width = 9,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(
    DIRS$figures,
    "Fig2A_GSE237156_PCA.pdf"
  ),
  p_pca,
  width = 9,
  height = 7
)

## Library-size and detection plot.
qc_long <- melt(
  sample_qc,
  id.vars = c(
    "sample_accession",
    "original_title",
    "group_id",
    "macrophage_subset",
    "diet",
    "drug"
  ),
  measure.vars = c(
    "library_size_rounded",
    "detected_count_ge_10"
  ),
  variable.name = "metric",
  value.name = "value"
)
fwrite(
  qc_long,
  file.path(
    DIRS$source,
    "FigS2A_sample_QC_source.csv"
  )
)
p_qc <- ggplot(
  qc_long,
  aes(
    x = sample_accession,
    y = value,
    fill = interaction(diet, drug)
  )
) +
  geom_col() +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 1
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Diet × drug",
    title = "GSE237156 library and detection metrics"
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(
    DIRS$figures,
    "FigS2A_GSE237156_sample_QC.png"
  ),
  p_qc,
  width = 9,
  height = 9,
  dpi = 300
)

## Correlation heatmap.
sample_cor <- cor(
  vst_matrix,
  method = "pearson"
)
fwrite(
  as.data.table(
    sample_cor,
    keep.rownames = "sample_accession"
  ),
  file.path(
    DIRS$source,
    "FigS2B_sample_correlation_source.csv"
  )
)

cor_annotation <- data.frame(
  subset = coldata$macrophage_subset,
  diet = coldata$diet,
  drug = coldata$drug
)
rownames(cor_annotation) <- rownames(coldata)

png(
  file.path(
    DIRS$figures,
    "FigS2B_GSE237156_sample_correlation.png"
  ),
  width = 2800,
  height = 2500,
  res = 300
)
pheatmap(
  sample_cor,
  annotation_col = cor_annotation,
  annotation_row = cor_annotation,
  border_color = NA,
  main = "Variance-stabilized sample correlation"
)
dev.off()

pdf(
  file.path(
    DIRS$figures,
    "FigS2B_GSE237156_sample_correlation.pdf"
  ),
  width = 10,
  height = 9
)
pheatmap(
  sample_cor,
  annotation_col = cor_annotation,
  annotation_row = cor_annotation,
  border_color = NA,
  main = "Variance-stabilized sample correlation"
)
dev.off()

############################################################
## 9. Define primary and supporting contrasts
############################################################

contrast_definitions <- rbindlist(list(
  data.table(
    subset = "Ccr2_positive",
    effect = "Disease_under_vehicle",
    contrast = "Ccr2_positive_HFD_Vehicle_vs_CD_Vehicle",
    numerator = "Ccr2_positive__HFD__Vehicle",
    denominator = "Ccr2_positive__CD__Vehicle",
    role = "Primary disease-associated effect"
  ),
  data.table(
    subset = "Ccr2_positive",
    effect = "Drug_under_HFD",
    contrast = "Ccr2_positive_HFD_Dapagliflozin_vs_HFD_Vehicle",
    numerator = "Ccr2_positive__HFD__Dapagliflozin",
    denominator = "Ccr2_positive__HFD__Vehicle",
    role = "Primary drug-reversal effect"
  ),
  data.table(
    subset = "Ccr2_positive",
    effect = "Drug_under_CD",
    contrast = "Ccr2_positive_CD_Dapagliflozin_vs_CD_Vehicle",
    numerator = "Ccr2_positive__CD__Dapagliflozin",
    denominator = "Ccr2_positive__CD__Vehicle",
    role = "Drug effect in control-diet context"
  ),
  data.table(
    subset = "Ccr2_positive",
    effect = "Disease_under_drug",
    contrast = "Ccr2_positive_HFD_Dapagliflozin_vs_CD_Dapagliflozin",
    numerator = "Ccr2_positive__HFD__Dapagliflozin",
    denominator = "Ccr2_positive__CD__Dapagliflozin",
    role = "Residual disease effect under drug"
  ),
  data.table(
    subset = "Ccr2_negative",
    effect = "Disease_under_vehicle",
    contrast = "Ccr2_negative_HFD_Vehicle_vs_CD_Vehicle",
    numerator = "Ccr2_negative__HFD__Vehicle",
    denominator = "Ccr2_negative__CD__Vehicle",
    role = "Primary disease-associated effect"
  ),
  data.table(
    subset = "Ccr2_negative",
    effect = "Drug_under_HFD",
    contrast = "Ccr2_negative_HFD_Dapagliflozin_vs_HFD_Vehicle",
    numerator = "Ccr2_negative__HFD__Dapagliflozin",
    denominator = "Ccr2_negative__HFD__Vehicle",
    role = "Primary drug-reversal effect"
  ),
  data.table(
    subset = "Ccr2_negative",
    effect = "Drug_under_CD",
    contrast = "Ccr2_negative_CD_Dapagliflozin_vs_CD_Vehicle",
    numerator = "Ccr2_negative__CD__Dapagliflozin",
    denominator = "Ccr2_negative__CD__Vehicle",
    role = "Drug effect in control-diet context"
  ),
  data.table(
    subset = "Ccr2_negative",
    effect = "Disease_under_drug",
    contrast = "Ccr2_negative_HFD_Dapagliflozin_vs_CD_Dapagliflozin",
    numerator = "Ccr2_negative__HFD__Dapagliflozin",
    denominator = "Ccr2_negative__CD__Dapagliflozin",
    role = "Residual disease effect under drug"
  )
))
fwrite(
  contrast_definitions,
  file.path(
    DIRS$tables,
    "07_contrast_definitions.csv"
  )
)

missing_levels <- setdiff(
  unique(
    c(
      contrast_definitions$numerator,
      contrast_definitions$denominator
    )
  ),
  levels(coldata$group)
)
if (length(missing_levels) > 0L) {
  stop(
    "Contrast group level(s) absent from DESeq2 design: ",
    paste(missing_levels, collapse = ", ")
  )
}

############################################################
## 10. DESeq2 primary analysis
############################################################

deseq_results_list <- list()
for (i in seq_len(nrow(contrast_definitions))) {
  cd <- contrast_definitions[i]
  log_msg(
    "DESeq2 contrast ",
    i,
    "/",
    nrow(contrast_definitions),
    ": ",
    cd$contrast
  )
  deseq_results_list[[cd$contrast]] <- make_deseq_result(
    dds = dds,
    contrast_name = cd$contrast,
    numerator = cd$numerator,
    denominator = cd$denominator
  )
}

deseq_all <- rbindlist(
  deseq_results_list,
  use.names = TRUE,
  fill = TRUE
)
deseq_all <- merge(
  deseq_all,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)

fwrite(
  deseq_all,
  file.path(
    DIRS$tables,
    "08_DESeq2_all_contrasts.csv.gz"
  ),
  compress = "gzip"
)

deseq_summary <- rbindlist(
  lapply(
    deseq_results_list,
    summarize_contrast,
    method_name = "DESeq2"
  )
)

############################################################
## 11. edgeR quasi-likelihood sensitivity analysis
############################################################

edge_group <- factor(
  coldata$group,
  levels = levels(coldata$group)
)
y <- DGEList(
  counts = count_matrix,
  group = edge_group
)
keep_edge <- filterByExpr(
  y,
  group = edge_group,
  min.count = MIN_COUNT
)
y <- y[keep_edge, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design_edge <- model.matrix(
  ~ 0 + edge_group
)
colnames(design_edge) <- levels(edge_group)
rownames(design_edge) <- rownames(coldata)

y <- estimateDisp(
  y,
  design_edge,
  robust = TRUE
)
fit_edge <- glmQLFit(
  y,
  design_edge,
  robust = TRUE
)
saveRDS(
  y,
  file.path(
    DIRS$objects,
    "GSE237156_edgeR_DGEList.rds"
  )
)
saveRDS(
  fit_edge,
  file.path(
    DIRS$objects,
    "GSE237156_edgeR_QL_fit.rds"
  )
)

edger_results_list <- list()
for (i in seq_len(nrow(contrast_definitions))) {
  cd <- contrast_definitions[i]
  log_msg(
    "edgeR QL contrast ",
    i,
    "/",
    nrow(contrast_definitions),
    ": ",
    cd$contrast
  )
  edger_results_list[[cd$contrast]] <- make_edger_result(
    fit = fit_edge,
    contrast_name = cd$contrast,
    numerator = cd$numerator,
    denominator = cd$denominator,
    design_columns = colnames(design_edge)
  )
}

edger_all <- rbindlist(
  edger_results_list,
  use.names = TRUE,
  fill = TRUE
)
edger_all <- merge(
  edger_all,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)
fwrite(
  edger_all,
  file.path(
    DIRS$tables,
    "09_edgeR_QL_all_contrasts.csv.gz"
  ),
  compress = "gzip"
)

edger_summary <- rbindlist(
  lapply(
    edger_results_list,
    summarize_contrast,
    method_name = "edgeR_QL"
  )
)

contrast_summary <- rbind(
  deseq_summary,
  edger_summary,
  fill = TRUE
)
contrast_summary <- merge(
  contrast_summary,
  contrast_definitions[, .(
    contrast,
    subset,
    effect,
    role
  )],
  by = "contrast",
  all.x = TRUE
)
setorder(
  contrast_summary,
  subset,
  effect,
  method
)
fwrite(
  contrast_summary,
  file.path(
    DIRS$tables,
    "10_contrast_result_summary.csv"
  )
)

############################################################
## 12. DESeq2-edgeR method concordance
############################################################

method_comparison_list <- list()
method_summary_list <- list()

for (contrast_name in contrast_definitions$contrast) {
  d <- deseq_all[
    contrast == contrast_name,
    .(
      gene_id,
      deseq_log2fc = log2fc_shrunk,
      deseq_pvalue = pvalue,
      deseq_padj = padj
    )
  ]
  e <- edger_all[
    contrast == contrast_name,
    .(
      gene_id,
      edger_log2fc = log2fc,
      edger_pvalue = pvalue,
      edger_padj = padj
    )
  ]

  m <- merge(
    d,
    e,
    by = "gene_id",
    all = FALSE
  )
  m[, contrast := contrast_name]
  m[, sign_agreement := (
    sign(deseq_log2fc) ==
      sign(edger_log2fc)
  )]
  m[, both_fdr_010 := (
    !is.na(deseq_padj) &
      deseq_padj <= EXPLORATORY_FDR &
      !is.na(edger_padj) &
      edger_padj <= EXPLORATORY_FDR
  )]
  method_comparison_list[[contrast_name]] <- m

  complete <- m[
    is.finite(deseq_log2fc) &
      is.finite(edger_log2fc)
  ]
  top_d <- head(
    complete[
      order(
        deseq_padj,
        -abs(deseq_log2fc),
        na.last = TRUE
      ),
      gene_id
    ],
    200L
  )
  top_e <- head(
    complete[
      order(
        edger_padj,
        -abs(edger_log2fc),
        na.last = TRUE
      ),
      gene_id
    ],
    200L
  )

  method_summary_list[[contrast_name]] <- data.table(
    contrast = contrast_name,
    common_tested_genes = nrow(complete),
    pearson_lfc = suppressWarnings(
      cor(
        complete$deseq_log2fc,
        complete$edger_log2fc,
        method = "pearson"
      )
    ),
    spearman_lfc = suppressWarnings(
      cor(
        complete$deseq_log2fc,
        complete$edger_log2fc,
        method = "spearman"
      )
    ),
    overall_sign_agreement = mean(
      complete$sign_agreement,
      na.rm = TRUE
    ),
    sign_agreement_among_any_fdr_010 = mean(
      complete[
        (
          deseq_padj <= EXPLORATORY_FDR |
            edger_padj <= EXPLORATORY_FDR
        ),
        sign_agreement
      ],
      na.rm = TRUE
    ),
    top200_jaccard = length(
      intersect(top_d, top_e)
    ) /
      length(
        union(top_d, top_e)
      )
  )
}

method_comparison <- rbindlist(
  method_comparison_list,
  use.names = TRUE,
  fill = TRUE
)
method_summary <- rbindlist(
  method_summary_list,
  use.names = TRUE,
  fill = TRUE
)
fwrite(
  method_comparison,
  file.path(
    DIRS$tables,
    "11_DESeq2_edgeR_gene_level_comparison.csv.gz"
  ),
  compress = "gzip"
)
fwrite(
  method_summary,
  file.path(
    DIRS$tables,
    "12_DESeq2_edgeR_method_concordance_summary.csv"
  )
)

## Method agreement plots for primary disease and drug contrasts.
primary_contrasts <- contrast_definitions[
  effect %in% c(
    "Disease_under_vehicle",
    "Drug_under_HFD"
  )
]
for (contrast_name in primary_contrasts$contrast) {
  plot_dt <- method_comparison[
    contrast == contrast_name
  ]
  plot_dt <- merge(
    plot_dt,
    gene_annotation[, .(
      gene_id,
      display_gene
    )],
    by = "gene_id",
    all.x = TRUE
  )
  label_dt <- plot_dt[
    order(
      pmin(
        deseq_padj,
        edger_padj,
        na.rm = TRUE
      ),
      -abs(deseq_log2fc)
    )
  ][seq_len(min(TOP_LABEL_GENES, .N))]

  p_method <- ggplot(
    plot_dt,
    aes(
      x = deseq_log2fc,
      y = edger_log2fc
    )
  ) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(
      aes(shape = sign_agreement),
      alpha = 0.5,
      size = 1.5
    ) +
    geom_text_repel(
      data = label_dt,
      aes(label = display_gene),
      size = 2.7,
      max.overlaps = Inf
    ) +
    labs(
      x = "DESeq2 shrunken log2FC",
      y = "edgeR QL log2FC",
      shape = "Sign agreement",
      title = paste0(
        "Method agreement: ",
        contrast_name
      )
    ) +
    theme_bw(base_size = 10)

  safe_name <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    contrast_name
  )
  ggsave(
    file.path(
      DIRS$figures,
      paste0(
        "FigS2C_method_agreement_",
        safe_name,
        ".png"
      )
    ),
    p_method,
    width = 7,
    height = 6,
    dpi = 300
  )
}

############################################################
## 13. Continuous drug-opposition ranking by macrophage subset
############################################################

build_opposition_table <- function(subset_name) {
  disease_name <- contrast_definitions[
    subset == subset_name &
      effect == "Disease_under_vehicle",
    contrast
  ]
  drug_name <- contrast_definitions[
    subset == subset_name &
      effect == "Drug_under_HFD",
    contrast
  ]

  disease_d <- deseq_all[
    contrast == disease_name,
    .(
      gene_id,
      disease_base_mean = base_mean,
      disease_lfc = log2fc_shrunk,
      disease_lfc_raw = log2fc_raw,
      disease_pvalue = pvalue,
      disease_padj = padj,
      disease_stat = stat
    )
  ]
  drug_d <- deseq_all[
    contrast == drug_name,
    .(
      gene_id,
      drug_base_mean = base_mean,
      drug_lfc = log2fc_shrunk,
      drug_lfc_raw = log2fc_raw,
      drug_pvalue = pvalue,
      drug_padj = padj,
      drug_stat = stat
    )
  ]

  disease_e <- edger_all[
    contrast == disease_name,
    .(
      gene_id,
      disease_edger_lfc = log2fc,
      disease_edger_pvalue = pvalue,
      disease_edger_padj = padj
    )
  ]
  drug_e <- edger_all[
    contrast == drug_name,
    .(
      gene_id,
      drug_edger_lfc = log2fc,
      drug_edger_pvalue = pvalue,
      drug_edger_padj = padj
    )
  ]

  out <- Reduce(
    function(x, y) merge(
      x,
      y,
      by = "gene_id",
      all = FALSE
    ),
    list(
      disease_d,
      drug_d,
      disease_e,
      drug_e
    )
  )
  out <- merge(
    out,
    gene_annotation,
    by = "gene_id",
    all.x = TRUE
  )

  out[, subset := subset_name]
  out[, deseq_opposed := (
    is.finite(disease_lfc) &
      is.finite(drug_lfc) &
      disease_lfc * drug_lfc < 0
  )]
  out[, edger_opposed := (
    is.finite(disease_edger_lfc) &
      is.finite(drug_edger_lfc) &
      disease_edger_lfc *
        drug_edger_lfc < 0
  )]
  out[, disease_method_sign_agreement := (
    sign(disease_lfc) ==
      sign(disease_edger_lfc)
  )]
  out[, drug_method_sign_agreement := (
    sign(drug_lfc) ==
      sign(drug_edger_lfc)
  )]
  out[, four_effect_signs_consistent := (
    disease_method_sign_agreement &
      drug_method_sign_agreement
  )]

  ## Positive values indicate opposition; negative values indicate
  ## disease and drug effects in the same direction.
  out[, opposition_effect_score := (
    -sign(disease_lfc * drug_lfc) *
      sqrt(
        abs(disease_lfc * drug_lfc)
      )
  )]

  out[, disease_abs_lfc_rank := safe_rank(
    abs(disease_lfc),
    decreasing = TRUE
  )]
  out[, drug_abs_lfc_rank := safe_rank(
    abs(drug_lfc),
    decreasing = TRUE
  )]
  out[, rank_product_effect := sqrt(
    disease_abs_lfc_rank *
      drug_abs_lfc_rank
  )]

  out[, disease_evidence_rank := safe_rank(
    safe_neglog10(disease_padj),
    decreasing = TRUE
  )]
  out[, drug_evidence_rank := safe_rank(
    safe_neglog10(drug_padj),
    decreasing = TRUE
  )]
  out[, rank_product_statistical := sqrt(
    disease_evidence_rank *
      drug_evidence_rank
  )]

  out[, combined_rank_product := (
    rank_product_effect *
      rank_product_statistical
  )^(1 / 2)]

  out[, opposition_tier := fcase(
    deseq_opposed &
      edger_opposed &
      four_effect_signs_consistent &
      disease_padj <= EXPLORATORY_FDR &
      drug_padj <= EXPLORATORY_FDR,
    "Tier_A_both_DESeq2_FDR_and_edgeR_direction",

    deseq_opposed &
      edger_opposed &
      four_effect_signs_consistent &
      (
        disease_padj <= EXPLORATORY_FDR |
          drug_padj <= EXPLORATORY_FDR
      ) &
      abs(disease_lfc) >= MIN_ABS_LFC_FOR_TIER &
      abs(drug_lfc) >= MIN_ABS_LFC_FOR_TIER,
    "Tier_B_one_DESeq2_FDR_effect_supported",

    deseq_opposed &
      edger_opposed &
      four_effect_signs_consistent &
      abs(disease_lfc) >= MIN_ABS_LFC_FOR_TIER &
      abs(drug_lfc) >= MIN_ABS_LFC_FOR_TIER,
    "Tier_C_effect_and_method_supported",

    deseq_opposed,
    "Tier_D_directional_opposition_only",

    default = "Not_opposed"
  )]

  tier_order <- c(
    "Tier_A_both_DESeq2_FDR_and_edgeR_direction",
    "Tier_B_one_DESeq2_FDR_effect_supported",
    "Tier_C_effect_and_method_supported",
    "Tier_D_directional_opposition_only",
    "Not_opposed"
  )
  out[, tier_order := match(
    opposition_tier,
    tier_order
  )]
  setorder(
    out,
    tier_order,
    combined_rank_product,
    -opposition_effect_score,
    gene_id
  )
  out[, within_subset_rank := seq_len(.N)]
  out[, tier_order := NULL]
  out
}

opposition_positive <- build_opposition_table(
  "Ccr2_positive"
)
opposition_negative <- build_opposition_table(
  "Ccr2_negative"
)

fwrite(
  opposition_positive,
  file.path(
    DIRS$tables,
    "13_opposition_rank_Ccr2_positive.csv.gz"
  ),
  compress = "gzip"
)
fwrite(
  opposition_negative,
  file.path(
    DIRS$tables,
    "14_opposition_rank_Ccr2_negative.csv.gz"
  ),
  compress = "gzip"
)

opposition_summary <- rbindlist(list(
  opposition_positive[, .(
    subset = "Ccr2_positive",
    tested_genes = .N,
    DESeq2_opposed = sum(deseq_opposed),
    edgeR_opposed = sum(edger_opposed),
    both_methods_opposed = sum(
      deseq_opposed &
        edger_opposed &
        four_effect_signs_consistent
    ),
    tier_A = sum(
      grepl("^Tier_A", opposition_tier)
    ),
    tier_B = sum(
      grepl("^Tier_B", opposition_tier)
    ),
    tier_C = sum(
      grepl("^Tier_C", opposition_tier)
    ),
    tier_D = sum(
      grepl("^Tier_D", opposition_tier)
    )
  )],
  opposition_negative[, .(
    subset = "Ccr2_negative",
    tested_genes = .N,
    DESeq2_opposed = sum(deseq_opposed),
    edgeR_opposed = sum(edger_opposed),
    both_methods_opposed = sum(
      deseq_opposed &
        edger_opposed &
        four_effect_signs_consistent
    ),
    tier_A = sum(
      grepl("^Tier_A", opposition_tier)
    ),
    tier_B = sum(
      grepl("^Tier_B", opposition_tier)
    ),
    tier_C = sum(
      grepl("^Tier_C", opposition_tier)
    ),
    tier_D = sum(
      grepl("^Tier_D", opposition_tier)
    )
  )]
))
fwrite(
  opposition_summary,
  file.path(
    DIRS$tables,
    "15_opposition_summary.csv"
  )
)

############################################################
## 14. Cross-subset consensus ranking
############################################################

pos_for_merge <- opposition_positive[, .(
  gene_id,
  symbol,
  display_gene,
  pos_disease_lfc = disease_lfc,
  pos_drug_lfc = drug_lfc,
  pos_disease_padj = disease_padj,
  pos_drug_padj = drug_padj,
  pos_deseq_opposed = deseq_opposed,
  pos_edger_opposed = edger_opposed,
  pos_method_sign_consistent = four_effect_signs_consistent,
  pos_opposition_score = opposition_effect_score,
  pos_tier = opposition_tier,
  pos_rank = within_subset_rank
)]
neg_for_merge <- opposition_negative[, .(
  gene_id,
  neg_disease_lfc = disease_lfc,
  neg_drug_lfc = drug_lfc,
  neg_disease_padj = disease_padj,
  neg_drug_padj = drug_padj,
  neg_deseq_opposed = deseq_opposed,
  neg_edger_opposed = edger_opposed,
  neg_method_sign_consistent = four_effect_signs_consistent,
  neg_opposition_score = opposition_effect_score,
  neg_tier = opposition_tier,
  neg_rank = within_subset_rank
)]

cross_subset <- merge(
  pos_for_merge,
  neg_for_merge,
  by = "gene_id",
  all = TRUE
)
cross_subset[, opposed_in_both_subsets := (
  pos_deseq_opposed &
    neg_deseq_opposed
)]
cross_subset[, disease_direction_concordant := (
  sign(pos_disease_lfc) ==
    sign(neg_disease_lfc)
)]
cross_subset[, drug_direction_concordant := (
  sign(pos_drug_lfc) ==
    sign(neg_drug_lfc)
)]
cross_subset[, both_methods_and_subsets_supported := (
  pos_deseq_opposed &
    neg_deseq_opposed &
    pos_edger_opposed &
    neg_edger_opposed &
    pos_method_sign_consistent &
    neg_method_sign_consistent
)]
cross_subset[, consensus_rank_product := sqrt(
  pos_rank * neg_rank
)]
cross_subset[, mean_opposition_score := rowMeans(
  cbind(
    pos_opposition_score,
    neg_opposition_score
  ),
  na.rm = TRUE
)]
cross_subset[, consensus_category := fcase(
  both_methods_and_subsets_supported &
    disease_direction_concordant &
    drug_direction_concordant,
  "Cross_subset_full_directional_consensus",

  opposed_in_both_subsets,
  "Opposed_in_both_subsets",

  pos_deseq_opposed |
    neg_deseq_opposed,
  "Subset_specific_opposition",

  default = "Not_opposed"
)]
category_order <- c(
  "Cross_subset_full_directional_consensus",
  "Opposed_in_both_subsets",
  "Subset_specific_opposition",
  "Not_opposed"
)
cross_subset[, category_order := match(
  consensus_category,
  category_order
)]
setorder(
  cross_subset,
  category_order,
  consensus_rank_product,
  -mean_opposition_score,
  gene_id
)
cross_subset[, overall_consensus_rank := seq_len(.N)]
cross_subset[, category_order := NULL]

fwrite(
  cross_subset,
  file.path(
    DIRS$tables,
    "16_cross_subset_consensus_ranking.csv.gz"
  ),
  compress = "gzip"
)

############################################################
## 15. Factorial LRT sensitivity analyses
############################################################

factorial_coldata <- as.data.frame(
  sample_meta[, .(
    sample_accession,
    subset = macrophage_subset,
    diet,
    drug = drug_factor
  )]
)
rownames(factorial_coldata) <- (
  factorial_coldata$sample_accession
)
factorial_coldata <- factorial_coldata[
  colnames(count_matrix_filtered),
  ,
  drop = FALSE
]

## Global evidence for any drug-associated term.
dds_lrt_drug <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = factorial_coldata,
  design = ~ subset * diet * drug
)
dds_lrt_drug <- DESeq(
  dds_lrt_drug,
  test = "LRT",
  reduced = ~ subset * diet,
  minReplicatesForReplace = Inf,
  quiet = TRUE
)
lrt_drug_res <- results(dds_lrt_drug)
lrt_drug <- data.table(
  gene_id = rownames(lrt_drug_res),
  base_mean = lrt_drug_res$baseMean,
  log2fc_display_coefficient = lrt_drug_res$log2FoldChange,
  stat = lrt_drug_res$stat,
  pvalue = lrt_drug_res$pvalue,
  padj = lrt_drug_res$padj,
  test = "LRT_any_drug_related_term",
  full_model = "subset*diet*drug",
  reduced_model = "subset*diet"
)
lrt_drug <- merge(
  lrt_drug,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)

## Heterogeneity of drug-diet interaction across macrophage subsets.
dds_lrt_threeway <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = factorial_coldata,
  design = ~ subset * diet * drug
)
dds_lrt_threeway <- DESeq(
  dds_lrt_threeway,
  test = "LRT",
  reduced = ~ subset * diet +
    subset * drug +
    diet * drug,
  minReplicatesForReplace = Inf,
  quiet = TRUE
)
lrt_threeway_res <- results(dds_lrt_threeway)
lrt_threeway <- data.table(
  gene_id = rownames(lrt_threeway_res),
  base_mean = lrt_threeway_res$baseMean,
  log2fc_display_coefficient = lrt_threeway_res$log2FoldChange,
  stat = lrt_threeway_res$stat,
  pvalue = lrt_threeway_res$pvalue,
  padj = lrt_threeway_res$padj,
  test = "LRT_three_way_interaction",
  full_model = "subset*diet*drug",
  reduced_model = "subset*diet + subset*drug + diet*drug"
)
lrt_threeway <- merge(
  lrt_threeway,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)

fwrite(
  lrt_drug,
  file.path(
    DIRS$tables,
    "17_factorial_LRT_any_drug_term.csv.gz"
  ),
  compress = "gzip"
)
fwrite(
  lrt_threeway,
  file.path(
    DIRS$tables,
    "18_factorial_LRT_three_way_interaction.csv.gz"
  ),
  compress = "gzip"
)
saveRDS(
  dds_lrt_drug,
  file.path(
    DIRS$objects,
    "GSE237156_factorial_LRT_drug_dds.rds"
  )
)
saveRDS(
  dds_lrt_threeway,
  file.path(
    DIRS$objects,
    "GSE237156_factorial_LRT_threeway_dds.rds"
  )
)

############################################################
## 16. Optional mouse-native Hallmark GSEA and pathway opposition
############################################################

gsea_all <- data.table()
pathway_opposition <- data.table()
gsea_status <- "NOT_REQUESTED"
gsea_database_mode <- "NOT_RUN"
gsea_section_error <- NULL

if (RUN_OPTIONAL_GSEA) {
  tryCatch({
    if (
      requireNamespace("fgsea", quietly = TRUE) &&
      requireNamespace("msigdbr", quietly = TRUE)
    ) {
      log_msg("Running optional mouse-native Hallmark GSEA.")
      gsea_status <- "RUNNING"

      msig_formals <- names(
        formals(msigdbr::msigdbr)
      )

      ## Prefer the native mouse MSigDB division and MH Hallmark collection.
      ## Fall back to human Hallmark ortholog mapping only for older msigdbr
      ## releases that do not expose db_species.
      if (
        "db_species" %in% msig_formals &&
        "collection" %in% msig_formals
      ) {
        hallmark <- msigdbr::msigdbr(
          db_species = "MM",
          species = "Mus musculus",
          collection = "MH"
        )
        gsea_database_mode <- "Mouse_MSigDB_MM_MH"
      } else if (
        "db_species" %in% msig_formals &&
        "category" %in% msig_formals
      ) {
        hallmark <- msigdbr::msigdbr(
          db_species = "MM",
          species = "Mus musculus",
          category = "MH"
        )
        gsea_database_mode <- "Mouse_MSigDB_MM_MH_legacy_category"
      } else if ("collection" %in% msig_formals) {
        hallmark <- msigdbr::msigdbr(
          species = "Mus musculus",
          collection = "H"
        )
        gsea_database_mode <- "Human_Hallmark_mouse_ortholog_fallback"
      } else {
        hallmark <- msigdbr::msigdbr(
          species = "Mus musculus",
          category = "H"
        )
        gsea_database_mode <- "Human_Hallmark_mouse_ortholog_legacy_fallback"
      }

      hallmark <- as.data.table(hallmark)
      symbol_candidates <- c(
        "gene_symbol",
        "target_gene",
        "db_gene_symbol"
      )
      symbol_column <- intersect(
        symbol_candidates,
        names(hallmark)
      )[1L]

      if (
        length(symbol_column) == 0L ||
        is.na(symbol_column) ||
        !"gs_name" %in% names(hallmark)
      ) {
        stop(
          "Expected gs_name and mouse gene-symbol columns were not found in msigdbr output."
        )
      }

      pathways <- split(
        hallmark[[symbol_column]],
        hallmark$gs_name
      )
      pathways <- lapply(
        pathways,
        function(x) unique(
          x[
            !is.na(x) &
              nzchar(x)
          ]
        )
      )
      pathways <- pathways[
        lengths(pathways) >= 15L
      ]

      if (length(pathways) == 0L) {
        stop("No eligible Hallmark pathways remained after gene-set filtering.")
      }

      gsea_records <- list()
      for (i in seq_len(nrow(primary_contrasts))) {
        contrast_name <- primary_contrasts$contrast[i]

        rank_dt <- deseq_all[
          contrast == contrast_name &
            !is.na(symbol) &
            nzchar(symbol) &
            is.finite(stat),
          .(
            symbol,
            stat
          )
        ]

        ## FIXED v2: setorder() accepts column names, not expressions such
        ## as -abs(stat). Create an explicit absolute-statistic column.
        rank_dt[, abs_stat := abs(stat)]
        setorder(
          rank_dt,
          symbol,
          -abs_stat
        )
        rank_dt <- rank_dt[, .SD[1L], by = symbol]

        ranks <- rank_dt$stat
        names(ranks) <- rank_dt$symbol
        ranks <- sort(
          ranks,
          decreasing = TRUE
        )

        if (length(ranks) < 1000L) {
          add_warning(
            "GSEA",
            contrast_name,
            paste0(
              "Only ",
              length(ranks),
              " unique ranked symbols were available; contrast skipped."
            )
          )
          next
        }

        fg <- tryCatch(
          fgsea::fgseaMultilevel(
            pathways = pathways,
            stats = ranks,
            minSize = 15,
            maxSize = 500,
            eps = 0
          ),
          error = function(e) {
            add_warning(
              "GSEA_CONTRAST",
              contrast_name,
              conditionMessage(e)
            )
            NULL
          }
        )

        if (!is.null(fg) && nrow(fg) > 0L) {
          fg <- as.data.table(fg)
          fg[, leadingEdge := vapply(
            leadingEdge,
            paste,
            collapse = ";",
            character(1)
          )]
          fg[, contrast := contrast_name]
          fg[, subset := primary_contrasts[
            contrast == contrast_name,
            subset
          ]]
          fg[, effect := primary_contrasts[
            contrast == contrast_name,
            effect
          ]]
          fg[, gene_set_database_mode := gsea_database_mode]
          gsea_records[[contrast_name]] <- fg
        }
      }

      if (length(gsea_records) > 0L) {
        gsea_all <- rbindlist(
          gsea_records,
          use.names = TRUE,
          fill = TRUE
        )
        fwrite(
          gsea_all,
          file.path(
            DIRS$tables,
            "19_Hallmark_GSEA_primary_contrasts.csv"
          )
        )

        pathway_list <- list()
        for (subset_name in unique(
          primary_contrasts$subset
        )) {
          disease_contrast <- primary_contrasts[
            subset == subset_name &
              effect == "Disease_under_vehicle",
            contrast
          ]
          drug_contrast <- primary_contrasts[
            subset == subset_name &
              effect == "Drug_under_HFD",
            contrast
          ]

          gd <- gsea_all[
            contrast == disease_contrast,
            .(
              pathway,
              disease_NES = NES,
              disease_padj = padj
            )
          ]
          gt <- gsea_all[
            contrast == drug_contrast,
            .(
              pathway,
              drug_NES = NES,
              drug_padj = padj
            )
          ]

          gp <- merge(
            gd,
            gt,
            by = "pathway",
            all = FALSE
          )
          gp[, subset := subset_name]
          gp[, pathway_opposed := (
            disease_NES * drug_NES < 0
          )]
          gp[, pathway_opposition_strength := (
            -sign(disease_NES * drug_NES) *
              sqrt(
                abs(disease_NES * drug_NES)
              )
          )]
          setorder(
            gp,
            -pathway_opposed,
            -pathway_opposition_strength,
            disease_padj,
            drug_padj
          )
          pathway_list[[subset_name]] <- gp
        }

        pathway_opposition <- rbindlist(
          pathway_list,
          use.names = TRUE,
          fill = TRUE
        )
        fwrite(
          pathway_opposition,
          file.path(
            DIRS$tables,
            "20_Hallmark_pathway_opposition.csv"
          )
        )
        gsea_status <- "COMPLETED"
      } else {
        gsea_status <- "NO_GSEA_RESULTS"
      }
    } else {
      gsea_status <- "OPTIONAL_PACKAGES_UNAVAILABLE"
      add_warning(
        "GSEA",
        "packages",
        "fgsea and/or msigdbr unavailable; optional GSEA was skipped."
      )
    }
  }, error = function(e) {
    gsea_section_error <<- conditionMessage(e)
  })

  ## GSEA is optional and must never terminate the core Stage 2 analysis.
  if (!is.null(gsea_section_error)) {
    gsea_status <- "FAILED_NONFATAL"
    add_warning(
      "GSEA_SECTION",
      "optional_mouse_hallmark",
      paste0(
        gsea_section_error,
        " Core DESeq2, edgeR, opposition ranking, and factorial analyses remain valid."
      )
    )
  }
}

############################################################
## 17. Core drug-opposition figures
############################################################

plot_opposition_scatter <- function(
  opposition_dt,
  subset_name,
  figure_stem
) {
  plot_dt <- copy(opposition_dt)
  label_dt <- plot_dt[
    deseq_opposed == TRUE
  ][
    order(
      combined_rank_product,
      -opposition_effect_score
    )
  ][
    seq_len(min(TOP_LABEL_GENES, .N))
  ]

  fwrite(
    plot_dt[, .(
      gene_id,
      display_gene,
      disease_lfc,
      drug_lfc,
      disease_padj,
      drug_padj,
      deseq_opposed,
      edger_opposed,
      opposition_tier,
      opposition_effect_score,
      combined_rank_product
    )],
    file.path(
      DIRS$source,
      paste0(
        figure_stem,
        "_source.csv"
      )
    )
  )

  p <- ggplot(
    plot_dt,
    aes(
      x = disease_lfc,
      y = drug_lfc
    )
  ) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(
      aes(shape = deseq_opposed),
      alpha = 0.55,
      size = 1.6
    ) +
    geom_text_repel(
      data = label_dt,
      aes(label = display_gene),
      size = 2.8,
      max.overlaps = Inf
    ) +
    labs(
      x = "Disease effect: HFD vehicle vs CD vehicle\nDESeq2 shrunken log2FC",
      y = "Drug effect: HFD dapagliflozin vs HFD vehicle\nDESeq2 shrunken log2FC",
      shape = "Opposed direction",
      title = paste0(
        subset_name,
        ": disease–drug directional relationship"
      ),
      subtitle = "Upper-left and lower-right quadrants indicate transcriptional opposition"
    ) +
    theme_bw(base_size = 10)

  ggsave(
    file.path(
      DIRS$figures,
      paste0(figure_stem, ".png")
    ),
    p,
    width = 8,
    height = 7,
    dpi = 300
  )
  ggsave(
    file.path(
      DIRS$figures,
      paste0(figure_stem, ".pdf")
    ),
    p,
    width = 8,
    height = 7
  )
}

plot_opposition_scatter(
  opposition_positive,
  "Ccr2-positive macrophages",
  "Fig2B_Ccr2_positive_disease_drug_opposition"
)
plot_opposition_scatter(
  opposition_negative,
  "Ccr2-negative macrophages",
  "Fig2C_Ccr2_negative_disease_drug_opposition"
)

## Heatmap of top consensus/opposed genes.
heatmap_candidates <- unique(c(
  cross_subset[
    consensus_category ==
      "Cross_subset_full_directional_consensus",
    gene_id
  ],
  opposition_positive[
    deseq_opposed == TRUE,
    gene_id
  ],
  opposition_negative[
    deseq_opposed == TRUE,
    gene_id
  ]
))
rank_lookup <- cross_subset[
  gene_id %in% heatmap_candidates
]
setorder(
  rank_lookup,
  overall_consensus_rank
)
heatmap_genes <- head(
  rank_lookup$gene_id,
  TOP_HEATMAP_GENES
)
heatmap_genes <- heatmap_genes[
  heatmap_genes %in% rownames(vst_matrix)
]

if (length(heatmap_genes) >= 2L) {
  heat_mat <- vst_matrix[
    heatmap_genes,
    ,
    drop = FALSE
  ]
  heat_z <- t(scale(t(heat_mat)))
  heat_z[!is.finite(heat_z)] <- 0

  label_map <- gene_annotation[
    match(
      rownames(heat_z),
      gene_id
    ),
    display_gene
  ]
  label_map[
    is.na(label_map) |
      !nzchar(label_map)
  ] <- rownames(heat_z)[
    is.na(label_map) |
      !nzchar(label_map)
  ]
  rownames(heat_z) <- make.unique(label_map)

  heat_annotation <- data.frame(
    subset = coldata$macrophage_subset,
    diet = coldata$diet,
    drug = coldata$drug
  )
  rownames(heat_annotation) <- rownames(coldata)

  heat_source <- as.data.table(
    heat_z,
    keep.rownames = "display_gene"
  )
  fwrite(
    heat_source,
    file.path(
      DIRS$source,
      "Fig2D_top_opposed_heatmap_source.csv"
    )
  )

  png(
    file.path(
      DIRS$figures,
      "Fig2D_top_drug_opposed_genes_heatmap.png"
    ),
    width = 3000,
    height = 3300,
    res = 300
  )
  pheatmap(
    heat_z,
    annotation_col = heat_annotation,
    show_colnames = FALSE,
    border_color = NA,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    main = "Top cross-subset and subset-specific drug-opposed genes"
  )
  dev.off()

  pdf(
    file.path(
      DIRS$figures,
      "Fig2D_top_drug_opposed_genes_heatmap.pdf"
    ),
    width = 10,
    height = 12
  )
  pheatmap(
    heat_z,
    annotation_col = heat_annotation,
    show_colnames = FALSE,
    border_color = NA,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    main = "Top cross-subset and subset-specific drug-opposed genes"
  )
  dev.off()
} else {
  add_warning(
    "FIGURE",
    "opposition_heatmap",
    "Fewer than two eligible genes were available; heatmap was skipped."
  )
}

############################################################
## 18. Compact workbook and parameter documentation
############################################################

workbook_path <- file.path(
  DIRS$tables,
  "21_GSE237156_stage2_key_results.xlsx"
)
wb <- createWorkbook()
write_sheet_safe(
  wb,
  "Sample_metadata",
  sample_meta
)
write_sheet_safe(
  wb,
  "Sample_file_mapping",
  sample_file_map
)
write_sheet_safe(
  wb,
  "Sample_QC",
  sample_qc
)
write_sheet_safe(
  wb,
  "Contrast_definitions",
  contrast_definitions
)
write_sheet_safe(
  wb,
  "Contrast_summary",
  contrast_summary
)
write_sheet_safe(
  wb,
  "Method_concordance",
  method_summary
)
write_sheet_safe(
  wb,
  "Opposition_summary",
  opposition_summary
)
write_sheet_safe(
  wb,
  "Top_Ccr2_positive",
  head(
    opposition_positive,
    1000L
  )
)
write_sheet_safe(
  wb,
  "Top_Ccr2_negative",
  head(
    opposition_negative,
    1000L
  )
)
write_sheet_safe(
  wb,
  "Cross_subset_top",
  head(
    cross_subset,
    1000L
  )
)
write_sheet_safe(
  wb,
  "Pathway_opposition",
  pathway_opposition
)
saveWorkbook(
  wb,
  workbook_path,
  overwrite = TRUE
)

parameter_table <- data.table(
  parameter = c(
    "Random seed",
    "Minimum count",
    "Minimum samples meeting minimum count",
    "DESeq2 alpha",
    "Formal FDR",
    "Exploratory FDR",
    "Tier minimum absolute log2FC",
    "DESeq2 design",
    "edgeR design",
    "Factorial full model",
    "Factorial reduced model for any drug term",
    "Factorial reduced model for three-way interaction",
    "Primary analysis unit",
    "Candidate preselection",
    "DESeq2 LFC shrinkage",
    "Optional GSEA status",
    "GSEA database mode",
    "Previous Stage 2 replacement"
  ),
  value = c(
    "20260714",
    as.character(MIN_COUNT),
    as.character(MIN_SAMPLES_WITH_MIN_COUNT),
    as.character(DESEQ_ALPHA),
    as.character(FORMAL_FDR),
    as.character(EXPLORATORY_FDR),
    as.character(MIN_ABS_LFC_FOR_TIER),
    "~ group",
    "~ 0 + group",
    "~ subset * diet * drug",
    "~ subset * diet",
    "~ subset*diet + subset*drug + diet*drug",
    "Biological sample",
    "None; all retained genes ranked",
    "DESeq2 normal prior; unshrunken only as explicit fallback",
    gsea_status,
    gsea_database_mode,
    as.character(REPLACE_PREVIOUS_STAGE2)
  ),
  rationale = c(
    "Reproducibility",
    "Remove extremely sparse genes",
    "Retain genes measurable in at least one experimental cell",
    "Exploratory small-n reporting threshold",
    "Conventional formal threshold",
    "Small-n hypothesis-generating threshold",
    "Used only for transparent evidence tiers, not hard candidate selection",
    "Direct estimable pairwise contrasts",
    "Independent count-based quasi-likelihood sensitivity method",
    "Tests full 2 x 2 x 2 design",
    "Global test of drug-associated terms",
    "Tests subset-specific diet-by-drug heterogeneity",
    "Avoid pseudoreplication",
    "Prevents circular prioritization of Nfkb1 or any fixed panel",
    "Stable shrinkage without requiring optional compiled ashr dependencies",
    "Optional pathway-level contextualization; nonfatal if unavailable",
    "Records whether mouse-native MH or an older-version fallback was used",
    "Deletes the incomplete FIXED_v1 outputs and regenerates a clean FIXED_v2 result set"
  )
)
fwrite(
  parameter_table,
  file.path(
    DIRS$methods,
    "stage2_parameters_and_rationale.csv"
  )
)

methods_text <- c(
  "HFpEF Stage 2: GSE237156 drug-opposed macrophage transcriptomic discovery",
  "",
  "Input:",
  "- Sixteen RSEM gene-results files reconstructed from GSE237156_RAW.tar.",
  "- Metadata derived from the Stage 1 locked sample manifest.",
  "",
  "Experimental design:",
  "- Ccr2-positive versus Ccr2-negative sorted cardiac macrophages.",
  "- Control diet versus high-fat diet.",
  "- Vehicle versus dapagliflozin.",
  "- Two biological samples per experimental cell.",
  "",
  "Primary differential-expression method:",
  "- DESeq2 on rounded RSEM expected counts.",
  "- Design: group factor containing the eight experimental cells.",
  "- Log2 fold changes were shrunk with the DESeq2 normal prior; unshrunken estimates were retained only as an explicit fallback.",
  "",
  "Sensitivity method:",
  "- edgeR quasi-likelihood generalized linear models with robust dispersion estimation.",
  "- Method agreement is reported by fold-change correlation, sign agreement, and top-rank overlap.",
  "",
  "Drug-opposition definition:",
  "- Disease effect: HFD vehicle versus CD vehicle within each Ccr2 subset.",
  "- Drug effect: HFD dapagliflozin versus HFD vehicle within the same subset.",
  "- Opposition requires opposite signs of the disease and drug effects.",
  "- Continuous opposition effect score: -sign(disease LFC x drug LFC) x sqrt(abs(product)).",
  "- All retained genes are ranked; evidence tiers are descriptive and do not define a fixed candidate list.",
  "",
  "Factorial sensitivity analyses:",
  "- Full model: subset x diet x drug.",
  "- Any drug-related term LRT: reduced model subset x diet.",
  "- Three-way interaction LRT: reduced model containing all pairwise interactions.",
  "",
  "Claim boundary:",
  "- With two samples per experimental cell, findings are hypothesis-generating.",
  "- Transcriptional opposition does not establish direct pharmacologic targeting or causal regulation.",
  "- No gene, including Nfkb1, was forced into the ranked output.",
  "",
  "Replacement behavior:",
  "- The incomplete FIXED_v1 output directory and CHECK archive are deleted before the clean FIXED_v2 run.",
  "- Optional Hallmark GSEA uses mouse-native MM/MH resources when supported and cannot abort the core analysis."
)
writeLines(
  methods_text,
  file.path(
    DIRS$methods,
    "stage2_methods_and_claim_boundaries.txt"
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
## 19. Run status and review package
############################################################

warnings_dt <- if (length(warning_records) > 0L) {
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
    "22_warnings_and_nonfatal_issues.csv"
  )
)

scientific_checks <- data.table(
  check = c(
    "Locked sample count",
    "Unique GSM count",
    "Eight groups",
    "Two samples per group",
    "Unique RSEM mapping",
    "Extracted RSEM files",
    "Filtered genes >= 5000",
    "Eight DESeq2 contrasts",
    "Eight edgeR contrasts",
    "Two opposition tables",
    "Cross-subset consensus nonempty",
    "Factorial drug LRT nonempty",
    "Factorial three-way LRT nonempty",
    "Prior incomplete Stage 2 paths removed"
  ),
  observed = c(
    nrow(sample_meta),
    uniqueN(sample_meta$sample_accession),
    nrow(group_check),
    ifelse(all(group_check$N == 2L), 2L, NA_integer_),
    sum(sample_file_map$mapping_status == "PASS"),
    sum(sample_file_map$extracted_exists),
    sum(gene_filter),
    uniqueN(deseq_all$contrast),
    uniqueN(edger_all$contrast),
    as.integer(
      nrow(opposition_positive) > 0L
    ) +
      as.integer(
        nrow(opposition_negative) > 0L
      ),
    nrow(cross_subset),
    nrow(lrt_drug),
    nrow(lrt_threeway),
    as.integer(
      all(
        replacement_records$deletion_succeeded
      )
    )
  ),
  expected = c(
    16L,
    16L,
    8L,
    2L,
    16L,
    16L,
    5000L,
    8L,
    8L,
    2L,
    1L,
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
    "at_least",
    "equal",
    "equal",
    "equal",
    "at_least",
    "at_least",
    "at_least",
    "equal"
  )
)
scientific_checks[, status := fifelse(
  comparison == "equal" &
    observed == expected,
  "PASS",
  fifelse(
    comparison == "at_least" &
      observed >= expected,
    "PASS",
    "FAIL"
  )
)]
fwrite(
  scientific_checks,
  file.path(
    DIRS$tables,
    "23_scientific_completion_checks.csv"
  )
)

script_copy_status <- "NOT_DETECTED"
if (
  length(SCRIPT_FILE) == 1L &&
  !is.na(SCRIPT_FILE) &&
  file.exists(SCRIPT_FILE)
) {
  script_methods <- file.path(
    DIRS$methods,
    "HFpEF_Stage2_GSE237156_Drug_Opposed_Discovery_FIXED_v2_REPLACE.R"
  )
  script_check <- file.path(
    DIRS$check,
    "HFpEF_Stage2_GSE237156_Drug_Opposed_Discovery_FIXED_v2_REPLACE.R"
  )
  copy1 <- file.copy(
    SCRIPT_FILE,
    script_methods,
    overwrite = TRUE
  )
  copy2 <- file.copy(
    SCRIPT_FILE,
    script_check,
    overwrite = TRUE
  )
  script_copy_status <- if (
    isTRUE(copy1) &&
      isTRUE(copy2)
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
  "COMPLETED_STAGE2_READY_FOR_REVIEW"
} else {
  "COMPLETED_STAGE2_REVIEW_REQUIRED"
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
  samples = nrow(sample_meta),
  raw_genes = nrow(count_matrix),
  retained_genes = sum(gene_filter),
  deseq_contrasts = uniqueN(deseq_all$contrast),
  edger_contrasts = uniqueN(edger_all$contrast),
  positive_subset_opposed_genes = sum(
    opposition_positive$deseq_opposed,
    na.rm = TRUE
  ),
  negative_subset_opposed_genes = sum(
    opposition_negative$deseq_opposed,
    na.rm = TRUE
  ),
  cross_subset_full_consensus_genes = sum(
    cross_subset$consensus_category ==
      "Cross_subset_full_directional_consensus",
    na.rm = TRUE
  ),
  optional_gsea_status = gsea_status,
  gsea_database_mode = gsea_database_mode,
  replaced_previous_stage2 = REPLACE_PREVIOUS_STAGE2,
  previous_v1_output_existed = any(
    replacement_records$path_type %in% c(
      "prior_incomplete_output_directory",
      "prior_incomplete_check_zip"
    ) &
      replacement_records$existed_before
  ),
  warnings = nrow(warnings_dt),
  script_copy_status = script_copy_status,
  scientific_checks_failed = sum(
    scientific_checks$status != "PASS"
  ),
  overall_status = overall_status
)
fwrite(
  run_status,
  file.path(
    DIRS$tables,
    "24_stage2_run_status.csv"
  )
)

readme <- c(
  "HFpEF Reanalysis Project - Stage 2",
  "GSE237156 drug-opposed macrophage transcriptomic discovery",
  "",
  paste0("Overall status: ", overall_status),
  paste0("Samples: ", nrow(sample_meta)),
  paste0("Raw genes: ", nrow(count_matrix)),
  paste0("Retained genes: ", sum(gene_filter)),
  paste0("DESeq2 contrasts: ", uniqueN(deseq_all$contrast)),
  paste0("edgeR contrasts: ", uniqueN(edger_all$contrast)),
  paste0("Optional GSEA status: ", gsea_status),
  paste0("GSEA database mode: ", gsea_database_mode),
  paste0("Previous incomplete Stage 2 replaced: ", REPLACE_PREVIOUS_STAGE2),
  paste0("Script snapshot status: ", script_copy_status),
  "",
  "Primary interpretation:",
  "- Disease-associated and dapagliflozin-associated effects were estimated separately within Ccr2-positive and Ccr2-negative macrophages.",
  "- All retained genes were ranked by continuous drug-opposition evidence.",
  "- No Nfkb1-centered or fixed-gene assumption was imposed.",
  "",
  "Do not interpret Stage 2 alone as causal or independently validated.",
  "Proceed to Stage 3 only after reviewing the CHECK package."
)
writeLines(
  readme,
  file.path(
    OUT_DIR,
    "README_stage2.txt"
  ),
  useBytes = TRUE
)

## Compact review tables.
fwrite(
  head(
    opposition_positive,
    TOP_CHECK_GENES
  ),
  file.path(
    DIRS$check,
    "TOP500_opposition_Ccr2_positive.csv"
  )
)
fwrite(
  head(
    opposition_negative,
    TOP_CHECK_GENES
  ),
  file.path(
    DIRS$check,
    "TOP500_opposition_Ccr2_negative.csv"
  )
)
fwrite(
  head(
    cross_subset,
    TOP_CHECK_GENES
  ),
  file.path(
    DIRS$check,
    "TOP500_cross_subset_consensus.csv"
  )
)

review_files <- c(
  file.path(
    DIRS$logs,
    "previous_stage2_replacement_audit.csv"
  ),
  file.path(
    DIRS$tables,
    "01_sample_metadata_used.csv"
  ),
  file.path(
    DIRS$tables,
    "02_group_count_validation.csv"
  ),
  file.path(
    DIRS$tables,
    "03_sample_to_RSEM_file_mapping.csv"
  ),
  file.path(
    DIRS$tables,
    "06_sample_QC_metrics.csv"
  ),
  file.path(
    DIRS$tables,
    "07_contrast_definitions.csv"
  ),
  file.path(
    DIRS$tables,
    "10_contrast_result_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "12_DESeq2_edgeR_method_concordance_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "15_opposition_summary.csv"
  ),
  file.path(
    DIRS$tables,
    "20_Hallmark_pathway_opposition.csv"
  ),
  workbook_path,
  file.path(
    DIRS$tables,
    "22_warnings_and_nonfatal_issues.csv"
  ),
  file.path(
    DIRS$tables,
    "23_scientific_completion_checks.csv"
  ),
  file.path(
    DIRS$tables,
    "24_stage2_run_status.csv"
  ),
  file.path(
    DIRS$methods,
    "stage2_parameters_and_rationale.csv"
  ),
  file.path(
    DIRS$methods,
    "stage2_methods_and_claim_boundaries.txt"
  ),
  file.path(
    DIRS$methods,
    "sessionInfo.txt"
  ),
  file.path(
    OUT_DIR,
    "README_stage2.txt"
  ),
  LOG_FILE,
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
check_manifest <- data.table(
  filename = basename(check_files),
  size_bytes = as.numeric(
    file.info(check_files)$size
  )
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
fwrite(
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

log_msg("Stage 2 analysis finished.")
log_msg("Overall status: ", overall_status)
log_msg("Retained genes: ", sum(gene_filter))
log_msg(
  "Ccr2-positive DESeq2-opposed genes: ",
  sum(
    opposition_positive$deseq_opposed,
    na.rm = TRUE
  )
)
log_msg(
  "Ccr2-negative DESeq2-opposed genes: ",
  sum(
    opposition_negative$deseq_opposed,
    na.rm = TRUE
  )
)
log_msg(
  "Cross-subset full-consensus genes: ",
  sum(
    cross_subset$consensus_category ==
      "Cross_subset_full_directional_consensus",
    na.rm = TRUE
  )
)
log_msg("CHECK package: ", CHECK_ZIP)

cat("\n============================================================\n")
cat("HFpEF Stage 2 GSE237156 analysis completed\n")
cat("Status: ", overall_status, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("CHECK: ", CHECK_ZIP, "\n", sep = "")
cat("Upload the CHECK package for review before Stage 3.\n")
cat("============================================================\n")
