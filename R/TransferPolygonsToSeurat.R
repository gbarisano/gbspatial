#!/usr/bin/env Rscript

# Load required libraries if needed
suppressPackageStartupMessages({
  # CreateSegmentation() / CreateFOV() live in SeuratObject (re-exported by Seurat)
  if (requireNamespace("SeuratObject", quietly = TRUE)) library(SeuratObject)
})

#' Transfer FastReSeg polygon coordinates into a Seurat spatial object
#'
#' Converts a nested list of per-FOV, per-cell polygon vertex matrices (as produced by
#' FastReSeg) into Seurat \code{Segmentation}/\code{FOV} objects and attaches them to a
#' Seurat object, creating one FOV slot per slide.
#'
#' @param polys A list with one element per FOV. Each FOV element is itself a *named* list
#'   with one element per cell (the names must be the cell IDs that match
#'   \code{colnames(seurat_object)}), and each cell element is a matrix or data.frame whose
#'   first two columns are the x and y vertex coordinates of that cell's polygon.
#' @param seurat_object A Seurat object to attach the segmentation FOV(s) to.
#' @param slide_col Name of the metadata column that groups cells into slides. One FOV slot
#'   is created per unique value of this column, named by that value (this mirrors the
#'   original per-slide loop). If \code{NULL}, all cells are placed into a single FOV named
#'   \code{fov_name}. Defaults to "slide".
#' @param assay Assay name passed to \code{CreateFOV}. Defaults to "RNA".
#' @param invert_x Logical; if TRUE the x coordinates are multiplied by -1. Defaults to FALSE.
#' @param invert_y Logical; if TRUE the y coordinates are multiplied by -1. FastReSeg's y-axis
#'   usually needs inverting to match the Seurat object's orientation, so defaults to TRUE.
#' @param divide_by_1000 Logical; if TRUE the coordinates are divided by 1000. Use this when
#'   the polygon coordinates are in microns while the Seurat object stores coordinates in
#'   millimeters (the usual CosMx case). Defaults to TRUE.
#' @param boundary Name of the segmentation boundary / FOV type. Defaults to "segmentation".
#' @param fov_name FOV slot name used only when \code{slide_col} is NULL. Defaults to "fov".
#' @param verbose Logical; print progress / coverage messages. Defaults to TRUE.
#' @return The Seurat object with one segmentation FOV added per slide (returned invisibly-free
#'   so it can be reassigned, e.g. \code{obj <- TransferPolygonsToSeurat(polys, obj)}).
#' @export
TransferPolygonsToSeurat <- function(polys,
                                     seurat_object,
                                     slide_col = "slide",
                                     assay = "RNA",
                                     invert_x = FALSE,
                                     invert_y = TRUE,
                                     divide_by_1000 = TRUE,
                                     boundary = "segmentation",
                                     fov_name = "fov",
                                     verbose = TRUE) {

  # Resolve the Seurat spatial constructors from whichever namespace provides them,
  # so the function works whether the user has them exported via SeuratObject or Seurat.
  .resolve_fn <- function(fn) {
    for (pkg in c("SeuratObject", "Seurat")) {
      if (requireNamespace(pkg, quietly = TRUE) &&
          exists(fn, where = asNamespace(pkg), inherits = FALSE)) {
        return(get(fn, envir = asNamespace(pkg)))
      }
    }
    stop("Could not find '", fn, "()' in SeuratObject or Seurat. Please install/load Seurat.")
  }
  CreateSegmentation_fn <- .resolve_fn("CreateSegmentation")
  CreateFOV_fn          <- .resolve_fn("CreateFOV")

  if (!is.list(polys) || length(polys) == 0) {
    stop("`polys` must be a non-empty list with one element per FOV.")
  }

  # Sign / scale factors derived from the options
  sx <- if (isTRUE(invert_x)) -1 else 1
  sy <- if (isTRUE(invert_y)) -1 else 1
  divisor <- if (isTRUE(divide_by_1000)) 1000 else 1

  # ---- 1. Flatten the nested list into a single x / y / cell data.frame ----
  # (Vectorized per FOV: one rbind of numeric matrices instead of one data.frame per cell.)
  build_fov_df <- function(sublist) {
    if (is.null(sublist) || length(sublist) == 0) return(NULL)
    cell_ids <- names(sublist)
    if (is.null(cell_ids) || any(cell_ids == "" | is.na(cell_ids))) {
      stop("Each per-FOV element of `polys` must be a *named* list (names = cell IDs).")
    }
    mats <- lapply(sublist, function(m) {
      m <- as.matrix(m)
      if (ncol(m) < 2) stop("Each cell's coordinates must have at least 2 columns (x, y).")
      m[, 1:2, drop = FALSE]
    })
    counts <- vapply(mats, nrow, integer(1))
    keep <- counts > 0
    if (!any(keep)) return(NULL)
    mats <- mats[keep]; cell_ids <- cell_ids[keep]; counts <- counts[keep]
    allmat <- do.call(rbind, mats)
    data.frame(
      x = sx * as.numeric(allmat[, 1]) / divisor,
      y = sy * as.numeric(allmat[, 2]) / divisor,
      cell = rep(cell_ids, times = counts),
      stringsAsFactors = FALSE
    )
  }

  poly_df <- do.call(rbind, lapply(polys, build_fov_df))
  if (is.null(poly_df) || nrow(poly_df) == 0) {
    stop("No polygon coordinates could be extracted from `polys`.")
  }
  rownames(poly_df) <- NULL

  # ---- 2. Coverage diagnostics ----
  obj_cells  <- colnames(seurat_object)
  poly_cells <- unique(poly_df$cell)
  if (isTRUE(verbose)) {
    n_no_poly <- sum(!obj_cells %in% poly_cells)
    n_no_cell <- sum(!poly_cells %in% obj_cells)
    message("Extracted polygons for ", length(poly_cells), " cells (",
            nrow(poly_df), " vertices).")
    if (n_no_poly > 0) message("  ", n_no_poly, " object cell(s) have no polygon and will lack segmentation.")
    if (n_no_cell > 0) message("  ", n_no_cell, " polygon cell(s) are not present in the object and will be ignored.")
  }

  # ---- 3. Determine the slide -> cells grouping ----
  if (is.null(slide_col)) {
    groups <- stats::setNames(list(obj_cells), fov_name)
  } else {
    meta <- seurat_object@meta.data
    if (!slide_col %in% colnames(meta)) {
      stop("slide_col '", slide_col, "' not found in seurat_object metadata.")
    }
    sv <- as.character(meta[[slide_col]])
    uniq <- unique(sv[!is.na(sv)])
    groups <- stats::setNames(lapply(uniq, function(u) obj_cells[which(sv == u)]), uniq)
  }

  # ---- 4. Build one Segmentation / FOV per slide and attach it ----
  added <- 0L
  for (slot_name in names(groups)) {
    cells_i <- groups[[slot_name]]
    sub_df <- poly_df[poly_df$cell %in% cells_i, , drop = FALSE]

    if (nrow(sub_df) == 0) {
      warning("No polygon coordinates for slide '", slot_name, "' - skipping.")
      next
    }
    if (isTRUE(verbose)) {
      message("Adding FOV '", slot_name, "' (", length(unique(sub_df$cell)), " cells).")
    }

    segmentation_obj <- CreateSegmentation_fn(coords = sub_df)
    fov_obj <- CreateFOV_fn(
      coords = stats::setNames(list(segmentation_obj), boundary),
      type   = boundary,
      assay  = assay
    )

    ok <- tryCatch({
      seurat_object[[slot_name]] <- fov_obj
      TRUE
    }, error = function(e) {
      warning("Failed to attach FOV '", slot_name, "': ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) added <- added + 1L
  }

  if (isTRUE(verbose)) message("Done. Added ", added, " FOV slot(s).")
  seurat_object
}
