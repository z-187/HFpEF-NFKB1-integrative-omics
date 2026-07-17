output_file <- file.path("environment", "sessionInfo.txt")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

sink(output_file)
cat("Captured:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "

")
print(sessionInfo())
cat("
Installed package versions:
")
packages <- sort(unique(c(
  "AnnotationDbi", "AUCell", "BiocManager", "BiocParallel", "data.table",
  "decoupleR", "DESeq2", "digest", "dorothea", "edgeR", "fgsea",
  "ggplot2", "ggrepel", "hdf5r", "limma", "Matrix", "matrixStats",
  "msigdbr", "openxlsx", "org.Mm.eg.db", "patchwork", "pheatmap",
  "scales", "scDblFinder", "Seurat", "SeuratObject",
  "SingleCellExperiment", "SummarizedExperiment", "writexl", "zip"
)))
for (package in packages) {
  version <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    "NOT_INSTALLED"
  }
  cat(sprintf("%-28s %s
", package, version))
}
sink()
message("Wrote ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
