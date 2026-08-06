#!/usr/bin/env Rscript

# Load required libraries if needed
suppressPackageStartupMessages({
  # CreateSegmentation() / CreateFOV() live in SeuratObject (re-exported by Seurat)
  if (requireNamespace("SeuratObject", quietly = TRUE)) library(SeuratObject)
})

#' Transfer FastReSeg polygon coordinates into a Seurat spatial object
#'
#' Converts polygon vertex coordinates into Seurat \code{Segmentation}/\code{FOV} objects and
#' attaches them to a Seurat object, creating one FOV slot per slide.
#'
#' The \code{polys} argument now accepts several forms (all describing the same underlying
#' polygon vertices):
#' \enumerate{
#'   \item \strong{Nested list} (the original format): a list with one element per FOV; each FOV
#'     element is a *named* list (names = cell IDs) whose elements are matrices/data.frames whose
#'     first two columns are the x and y vertex coordinates of that cell's polygon.
#'   \item \strong{Flat table}: a \code{data.frame} / \code{data.table} in long format, one vertex
#'     per row, with columns for the cell ID and the x / y coordinates (e.g. a CosMx
#'     \code{*-polygons} table with \code{cell}, \code{x_global_px}, \code{y_global_px}, ...).
#'   \item \strong{File path(s)}: a character vector of one \emph{or more} paths, each pointing to
#'     \itemize{
#'       \item an \code{.rds}/\code{.RDS} file (read with \code{readRDS()}; may contain either a
#'             nested list or a flat table), or
#'       \item a delimited text file - \code{.csv}, \code{.csv.gz}, \code{.tsv}, \code{.txt},
#'             or \code{.gz} - read as a flat table.
#'     }
#'     When more than one path is given (e.g. one polygon file per slide), every file is read and
#'     the vertices are pooled; each slide's polygons are then routed to the correct FOV by matching
#'     cell IDs, so you do \emph{not} need to state which file belongs to which slide. Files may mix
#'     formats (RDS + csv.gz) but should share the same column schema / coordinate units.
#'   \item \strong{List of sources}: a list whose elements are any mix of the above file paths and/or
#'     in-memory flat tables (equivalent to passing a character vector of paths).
#' }
#'
#' @param polys Polygon source(s). One of: a nested per-FOV list; a flat
#'   \code{data.frame}/\code{data.table}; a character vector of one or more file paths
#'   (\code{.rds}/\code{.RDS} or \code{.csv(.gz)}/\code{.tsv}); or a list of such sources. When
#'   several sources are supplied their vertices are combined. See Details above.
#' @param seurat_object A Seurat object to attach the segmentation FOV(s) to.
#' @param slide_col Name of the \emph{Seurat metadata} column that groups cells into slides. One
#'   FOV slot is created per unique value of this column, named by that value. If \code{NULL}, all
#'   cells are placed into a single FOV named \code{fov_name}. Defaults to "slide". (This concerns
#'   the Seurat object's metadata and is unrelated to any column in a flat polygon table.)
#' @param cell_col For \strong{flat-table} input only: name of the column holding cell IDs (these
#'   must match \code{colnames(seurat_object)}). If \code{NULL} (default) it is auto-detected from
#'   common CosMx names (\code{cell}, \code{cell_ID}, \code{cellID}, ...) by choosing whichever
#'   column overlaps most with \code{colnames(seurat_object)}.
#' @param x_col,y_col For \strong{flat-table} input only: names of the x / y coordinate columns. If
#'   \code{NULL} (default) they are auto-detected, preferring a continuous whole-slide coordinate
#'   pair in this order: \code{x_global_px}/\code{y_global_px}, then \code{x_slide_mm}/\code{y_slide_mm},
#'   then \code{x_local_px}/\code{y_local_px}, then \code{x}/\code{y}. Note that local (per-FOV) pixel
#'   coordinates reset in every FOV and will overlap if combined, so a global/slide pair is preferred.
#' @param assay Assay name passed to \code{CreateFOV}. Defaults to "RNA".
#' @param invert_x Logical; if TRUE the x coordinates are multiplied by -1. Defaults to FALSE.
#' @param invert_y Logical; if TRUE the y coordinates are multiplied by -1. FastReSeg's y-axis
#'   usually needs inverting to match the Seurat object's orientation, so defaults to TRUE.
#' @param divide_by_1000 Logical; if TRUE the coordinates are divided by 1000. Use this when the
#'   polygon coordinates are in microns while the Seurat object stores coordinates in millimeters
#'   (the usual CosMx case). Defaults to TRUE. \strong{Units caveat for flat tables:} if you select
#'   the \code{*_slide_mm} columns the values are already in mm (set \code{divide_by_1000 = FALSE});
#'   raw pixel columns are \emph{not} microns, so \code{/1000} will not convert them to mm - check
#'   the coordinate range printed by \code{verbose} and set the transform flags accordingly.
#' @param boundary Name of the segmentation boundary / FOV type. Defaults to "segmentation".
#' @param fov_name FOV slot name used only when \code{slide_col} is NULL. Defaults to "fov".
#' @param verbose Logical; print progress / coverage messages. Defaults to TRUE.
#' @return The Seurat object with one segmentation FOV added per slide.
#' @export
TransferPolygonsToSeurat <- function(polys,
                                     seurat_object,
                                     slide_col = "slide",
                                     cell_col = NULL,
                                     x_col = NULL,
                                     y_col = NULL,
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

  obj_cells <- colnames(seurat_object)

  # ---- 0. Normalize the `polys` argument into RAW x / y / cell data.frame -------------------
  # `polys` may be a file path, a flat data.frame/data.table, or a nested per-FOV list.
  # Whatever the form, produce a data.frame with columns x, y, cell holding the *untransformed*
  # coordinates; the sign/scale transforms are applied once, centrally, further below.

  # If a file path was given, read it into an R object first.
  # `polys` may describe ONE source or SEVERAL. Build a list of individual sources, convert
  # each to a RAW (x, y, cell) data.frame, then row-bind them. Pooling all vertices and then
  # routing them to slides by cell-ID membership (step 3) means each slide's file lands in the
  # right FOV automatically - no explicit file->slide mapping is required.
  norm_args <- list(cell_col = cell_col, x_col = x_col, y_col = y_col,
                    obj_cells = obj_cells, verbose = verbose)

  if (is.character(polys)) {
    # One or more file paths (.rds / .csv(.gz) / .tsv). Any length >= 1 is allowed.
    if (length(polys) == 0) stop("`polys` is an empty character vector - no file paths given.")
    sources <- as.list(polys)
  } else if (is.data.frame(polys)) {
    # A single already-loaded flat table (data.frame / data.table).
    sources <- list(polys)
  } else if (is.list(polys)) {
    # Either (a) a list of sources - paths and/or data.frames, one per slide/file - or
    # (b) the original single nested per-FOV list (elements are named lists of matrices).
    if (length(polys) == 0) stop("`polys` is an empty list.")
    is_source_list <- all(vapply(polys, function(e)
      (is.character(e) && length(e) == 1L) || is.data.frame(e), logical(1)))
    sources <- if (is_source_list) polys else list(polys)
  } else {
    stop("`polys` must be a nested list, a data.frame/data.table, a file path, ",
         "or a vector/list of such sources; got an object of class '",
         paste(class(polys), collapse = "/"), "'.")
  }

  multi <- length(sources) > 1L
  parts <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    if (isTRUE(verbose) && multi) {
      lbl <- if (is.character(sources[[i]])) sources[[i]] else class(sources[[i]])[1]
      message("--- Source ", i, "/", length(sources), ": ", lbl, " ---")
    }
    parts[[i]] <- do.call(.source_to_raw_df, c(list(src = sources[[i]]), norm_args))
  }
  poly_df <- do.call(rbind, parts)

  if (is.null(poly_df) || nrow(poly_df) == 0) {
    stop("No polygon coordinates could be extracted from `polys`.")
  }
  rownames(poly_df) <- NULL

  # ---- 1. Apply sign / scale transforms once (identical maths for every input form) ----
  sx <- if (isTRUE(invert_x)) -1 else 1
  sy <- if (isTRUE(invert_y)) -1 else 1
  divisor <- if (isTRUE(divide_by_1000)) 1000 else 1
  poly_df$x <- sx * poly_df$x / divisor
  poly_df$y <- sy * poly_df$y / divisor

  # ---- 2. Coverage diagnostics ----
  poly_cells <- unique(poly_df$cell)
  if (isTRUE(verbose)) {
    n_no_poly <- sum(!obj_cells %in% poly_cells)
    n_no_cell <- sum(!poly_cells %in% obj_cells)
    message("Extracted polygons for ", length(poly_cells), " cells (",
            nrow(poly_df), " vertices).")
    message("  x range [", signif(min(poly_df$x), 4), ", ", signif(max(poly_df$x), 4), "]",
            "  y range [", signif(min(poly_df$y), 4), ", ", signif(max(poly_df$y), 4), "]",
            "  (after invert/scale transforms)")
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


# ============================================================================================
# Internal helpers
# ============================================================================================

#' Convert ONE polygon source into a RAW (untransformed) x / y / cell data.frame.
#' A source is a file path (.rds/.csv(.gz)/.tsv), an in-memory flat data.frame/data.table,
#' or a nested per-FOV list. Column detection is performed per source, so when several files
#' are combined they should share the same schema / units (a global/slide coordinate pair is
#' recommended so coordinates are comparable across files).
#' @keywords internal
#' @noRd
.source_to_raw_df <- function(src, cell_col, x_col, y_col, obj_cells, verbose = TRUE) {
  if (is.character(src)) {
    if (length(src) != 1L) stop("Each file-path source must be a single string.")
    src <- .read_polys_file(src, verbose = verbose)
  }
  if (is.data.frame(src)) {
    return(.flat_table_to_raw_df(src, cell_col = cell_col, x_col = x_col, y_col = y_col,
                                 obj_cells = obj_cells, verbose = verbose))
  }
  if (is.list(src)) {
    if (length(src) == 0) stop("A nested-list source is empty.")
    return(do.call(rbind, lapply(src, .build_fov_raw_df)))
  }
  stop("Unsupported polygon source of class '", paste(class(src), collapse = "/"), "'.")
}

#' Read a polygon source file (.rds/.RDS or delimited text) into an R object.
#' @keywords internal
#' @noRd
.read_polys_file <- function(path, verbose = TRUE) {
  if (!file.exists(path)) stop("File not found: ", path)
  lower <- tolower(path)

  if (grepl("\\.rds$", lower)) {
    if (isTRUE(verbose)) message("Reading RDS: ", path)
    return(readRDS(path))
  }

  # Otherwise treat as a delimited text table (optionally gzipped).
  # Strip a trailing .gz to detect the real extension for separator inference.
  base_lower <- sub("\\.gz$", "", lower)
  sep <- if (grepl("\\.(tsv|tab)$", base_lower)) "\t" else ","

  if (isTRUE(verbose)) message("Reading table: ", path)

  # Prefer data.table::fread (fast, handles .gz natively); fall back to base read.table.
  if (requireNamespace("data.table", quietly = TRUE)) {
    df <- tryCatch(
      data.table::fread(path, sep = sep, showProgress = isTRUE(verbose),
                        data.table = FALSE),
      error = function(e) NULL
    )
    if (!is.null(df)) return(df)
    if (isTRUE(verbose)) message("  data.table::fread failed; falling back to read.table().")
  }

  con <- if (grepl("\\.gz$", lower)) gzfile(path) else path
  utils::read.table(con, header = TRUE, sep = sep,
                    stringsAsFactors = FALSE, check.names = FALSE)
}


#' Choose the cell-ID column of a flat polygon table.
#' Auto-detects by maximal overlap with the Seurat object's cell names.
#' @keywords internal
#' @noRd
.pick_cell_col <- function(df, cell_col, obj_cells) {
  cols <- colnames(df)
  if (!is.null(cell_col)) {
    if (!cell_col %in% cols) {
      stop("cell_col '", cell_col, "' not found. Available columns: ",
           paste(cols, collapse = ", "))
    }
    return(cell_col)
  }
  candidates <- c("cell", "cell_ID", "cellID", "cell_id", "CellId", "cellId", "cell_id_int")
  present <- candidates[candidates %in% cols]
  if (length(present) == 0) {
    stop("Could not auto-detect a cell-ID column. Please set `cell_col`. ",
         "Available columns: ", paste(cols, collapse = ", "))
  }
  # Prefer the candidate whose values overlap the object's cell names the most.
  overlaps <- vapply(present, function(cc) {
    sum(unique(as.character(df[[cc]])) %in% obj_cells)
  }, numeric(1))
  if (max(overlaps) > 0) {
    return(present[which.max(overlaps)])
  }
  # No overlap with the object at all - fall back to the first candidate and warn.
  warning("None of the candidate cell-ID columns (", paste(present, collapse = ", "),
          ") overlap colnames(seurat_object). Using '", present[1],
          "'; check that IDs match the object.")
  present[1]
}


#' Choose the x / y coordinate columns of a flat polygon table.
#' Prefers a continuous whole-slide coordinate pair over per-FOV local pixels.
#' @keywords internal
#' @noRd
.pick_coord_cols <- function(df, x_col, y_col, verbose = TRUE) {
  cols <- colnames(df)

  if (!is.null(x_col) || !is.null(y_col)) {
    if (is.null(x_col) || is.null(y_col)) {
      stop("Set both `x_col` and `y_col`, or neither (for auto-detection).")
    }
    for (cc in c(x_col, y_col)) {
      if (!cc %in% cols) {
        stop("coordinate column '", cc, "' not found. Available columns: ",
             paste(cols, collapse = ", "))
      }
    }
    return(c(x_col, y_col))
  }

  # Priority: global px -> slide mm -> local px -> generic x/y.
  pairs <- list(
    c("x_global_px", "y_global_px"),
    c("x_slide_mm",  "y_slide_mm"),
    c("x_local_px",  "y_local_px"),
    c("x",           "y")
  )
  for (p in pairs) {
    if (all(p %in% cols)) {
      if (identical(p, c("x_local_px", "y_local_px")) && isTRUE(verbose)) {
        warning("Using local (per-FOV) pixel columns; these reset per FOV and will overlap ",
                "across FOVs in a combined coordinate space. Prefer a global/slide pair if available.")
      }
      return(p)
    }
  }

  # Fallback: first two numeric columns.
  num <- cols[vapply(df, is.numeric, logical(1))]
  if (length(num) >= 2) {
    if (isTRUE(verbose)) {
      message("No standard coordinate columns found; using first two numeric columns: ",
              num[1], ", ", num[2], ".")
    }
    return(num[1:2])
  }
  stop("Could not auto-detect x / y coordinate columns. Please set `x_col` and `y_col`. ",
       "Available columns: ", paste(cols, collapse = ", "))
}


#' Convert a flat long-format polygon table into a RAW (untransformed) x / y / cell data.frame.
#' @keywords internal
#' @noRd
.flat_table_to_raw_df <- function(df, cell_col, x_col, y_col, obj_cells, verbose = TRUE) {
  if (nrow(df) == 0) stop("The flat polygon table has no rows.")

  cc <- .pick_cell_col(df, cell_col, obj_cells)
  xy <- .pick_coord_cols(df, x_col, y_col, verbose = verbose)
  if (isTRUE(verbose)) {
    message("Flat table: using cell column '", cc, "', coordinates '",
            xy[1], "' / '", xy[2], "'.")
  }

  # Pull out just the three vectors we need (avoids copying the whole - possibly huge - table).
  as_num <- function(v) if (is.numeric(v)) v else as.numeric(as.character(v))
  x <- as_num(df[[xy[1]]])
  y <- as_num(df[[xy[2]]])
  cell <- as.character(df[[cc]])

  keep <- !(is.na(x) | is.na(y) | is.na(cell) | cell == "")
  if (!all(keep)) {
    x <- x[keep]; y <- y[keep]; cell <- cell[keep]
  }
  if (length(x) == 0) stop("No usable rows (all coordinates/cell IDs were NA or empty).")

  data.frame(x = x, y = y, cell = cell, stringsAsFactors = FALSE)
}


#' Convert one FOV element of a nested per-FOV list into a RAW x / y / cell data.frame.
#' (Vectorized per FOV: one rbind of numeric matrices instead of one data.frame per cell.)
#' @keywords internal
#' @noRd
.build_fov_raw_df <- function(sublist) {
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
    x = as.numeric(allmat[, 1]),
    y = as.numeric(allmat[, 2]),
    cell = rep(cell_ids, times = counts),
    stringsAsFactors = FALSE
  )
}
