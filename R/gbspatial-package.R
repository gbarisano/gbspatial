#' gbspatial: functions for spatial transcriptomic data
#'
#' @keywords internal
"_PACKAGE"

# ------------------------------------------------------------------------------
# Quiet R CMD check's "checking R code for possible problems" NOTE.
#
# 1) Base functions used across the package (the check's "Consider adding
#    importFrom(...)" suggestion). Importing them removes the
#    "no visible global function definition" lines.
# 2) `.data` is the tidy-eval pronoun used in ggplot2 aes(); ggplot2 re-exports
#    it, so no new dependency is needed.
# 3) Column names referenced inside dplyr verbs are non-standard-evaluation
#    symbols; declaring them as global variables removes the
#    "no visible binding for global variable" lines.
# ------------------------------------------------------------------------------

#' @importFrom graphics legend par plot.new rect text
#' @importFrom grDevices colorRampPalette dev.control dev.off pdf png recordPlot
#' @importFrom stats aggregate ave lm median na.omit setNames
#' @importFrom utils capture.output installed.packages read.csv setTxtProgressBar txtProgressBar
#' @importFrom ggplot2 .data
NULL

utils::globalVariables(c(
  # run_spatial_qc.R
  "Category", "n", "Percentage", "Label", "nCount_RNA", "Area",
  "SplitRatioToLocal", "x_slide_mm", "y_slide_mm", "SBR",
  # other pre-existing functions
  "geometry", "grp", "vertex_order", "x", "y", "cell",
  # SpatialLRContrast.R volcano (only if you switch from .data$ to bare names)
  "log2FC", "padj", "contrast"
))
