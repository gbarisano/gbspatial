#!/usr/bin/env Rscript

#' Read manually traced QuPath ROIs (GeoJSON) into a coordinate data.frame (one slide)
#'
#' Reads a QuPath \code{.geojson}/\code{.json} export of manually drawn regions of interest and returns a tidy
#' data.frame of polygon vertices (columns \code{aoi_id, x, y}). QuPath coordinates are in pixels of the (possibly
#' down/up-sampled) image the ROIs were drawn on; they are converted to microns via \code{pixel_size_um} and
#' \code{downsample_factor}. If a FOV positions file is supplied, the coordinates are additionally mapped into the
#' same coordinate frame as the spatial object (see details), which is the frame you need for overlaying on cells.
#'
#' @details
#' Mapping happens in two steps. First, QuPath pixels are converted to GLOBAL microns using the FOV positions file:
#' x scales directly, and y is flipped because QuPath's y grows top->bottom while the metadata frame grows
#' bottom->top (\code{Yg = y_top - y_px * pixel_size_um * downsample_factor}, with
#' \code{y_top = min(y_global_mm)*1000 + fov_size_px*pixel_size_um}). Second, the coordinates are negated according to
#' \code{x_invert}/\code{y_invert} to match how the cell polygons were imported, and scaled to the output units.
#'
#' IMPORTANT: \code{y_invert} (and \code{x_invert}) must match the \code{invert_y}/\code{invert_x} you used in
#' \code{TransferPolygonsToSeurat} (or however the cell polygons were put into the object). The defaults here
#' (\code{y_invert = TRUE}, \code{x_invert = FALSE}) match \code{TransferPolygonsToSeurat}'s defaults, so if you
#' imported cells with those defaults you should NOT override them. Setting the wrong \code{y_invert} shifts the ROIs
#' by roughly \code{2 * min(y_global_mm)} along y.
#'
#' The returned coordinates are in the units set by \code{coord_scaling} (default "divide_by_1000" = millimeters,
#' matching the Seurat object). To overlay them with \code{PlotFOVSegmentation(..., extra_polygons = )}, set
#' \code{extra_polygons_invert_y = FALSE}, \code{extra_polygons_invert_x = FALSE} (inversion already applied here),
#' and \code{extra_polygons_divide_by_1000 = FALSE} when the output is already mm (default) or TRUE for microns.
#'
#' @param json_file Path to the QuPath GeoJSON file for this slide.
#' @param fov_pos_file Path to (or an already-read data.frame of) the FOV positions file for the SAME slide. If NULL,
#'   coordinates are only converted to image microns and NOT mapped to the metadata frame.
#' @param downsample_factor Numeric. Factor by which the ROI image was down-sampled relative to the original image the
#'   FOV positions refer to (e.g. 8 for an 8x-downsampled image). Use a fractional value (e.g. 0.5) for up-sampling.
#'   Defaults to 1.
#' @param fov_size_px Numeric. Size (edge length) of one FOV in pixels of the original image. Defaults to 4256.
#' @param pixel_size_um Numeric. Physical size of one original-image pixel, in microns. Defaults to 0.12028.
#' @param coord_scaling One of "divide_by_1000" (default; output in mm, matching the Seurat object), "none"
#'   (output in microns), or "multiply_by_1000". Controls the final unit of all output coordinates AND the offset,
#'   consistently.
#' @param y_invert,x_invert Logical. Negate the (global) y / x coordinate to match how the cell polygons were
#'   imported into the object. Must match \code{TransferPolygonsToSeurat}'s \code{invert_y}/\code{invert_x}. Defaults
#'   \code{y_invert = TRUE}, \code{x_invert = FALSE}. (The QuPath top->bottom y-flip is applied separately and always.)
#' @param x_offset,y_offset Numeric, in GLOBAL microns. \code{x_offset} shifts the x origin (default 0). \code{y_offset}
#'   overrides the computed \code{y_top} (global-micron y that QuPath y = 0 maps to); if NULL it is computed as in
#'   Details.
#' @param fov_x_col,fov_y_col Column names in the FOV positions file. Default "x_global_mm" / "y_global_mm".
#' @param slide Optional slide identifier; if given, added as a \code{slide} column (handy for rbind across slides).
#' @param cell_polygons Optional. If supplied (a data.frame of cell polygon vertices or a Seurat object), the function
#'   also computes, for each cell, the ROI it mostly falls in, and returns a list with elements \code{aoi_polygons}
#'   and \code{aoi_cells}. See \code{AssignCellsToAOIs}. Extra arguments are forwarded to it via \code{...}.
#' @param ... Passed to \code{AssignCellsToAOIs} when \code{cell_polygons} is supplied.
#' @return A data.frame with columns \code{aoi_id, x, y} (and \code{slide} if provided). If \code{cell_polygons} is
#'   supplied, a list with elements \code{aoi_polygons} (that data.frame) and \code{aoi_cells} (the per-cell
#'   assignment).
#' @export
ReadQuPathROIs <- function(json_file,
                           fov_pos_file = NULL,
                           downsample_factor = 1,
                           fov_size_px = 4256,
                           pixel_size_um = 0.12028,
                           coord_scaling = c("divide_by_1000", "none", "multiply_by_1000"),
                           y_invert = TRUE,
                           x_invert = FALSE,
                           x_offset = 0,
                           y_offset = NULL,
                           fov_x_col = "x_global_mm",
                           fov_y_col = "y_global_mm",
                           slide = NULL,
                           cell_polygons = NULL,
                           ...) {

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read GeoJSON files.")
  }

  gj <- jsonlite::fromJSON(json_file, simplifyVector = FALSE)

  # Accept a FeatureCollection, a single Feature, or a bare list of features
  if (!is.null(gj$type) && identical(gj$type, "FeatureCollection")) {
    features <- gj$features
  } else if (!is.null(gj$type) && identical(gj$type, "Feature")) {
    features <- list(gj)
  } else {
    features <- gj
  }
  if (length(features) == 0) stop("No features found in ", json_file)

  # ---- 1. Flatten features into aoi_id / x / y rows ----
  rows <- list()
  for (k in seq_along(features)) {
    feature <- features[[k]]
    geom <- feature$geometry
    if (is.null(geom) || is.null(geom$coordinates)) next

    poly_id <- feature$id
    if (is.null(poly_id)) poly_id <- feature$properties$name
    if (is.null(poly_id)) poly_id <- paste0("roi_", k)
    poly_id <- as.character(poly_id)

    gtype <- geom$type
    rings <- list()
    if (identical(gtype, "Polygon")) {
      rings[[1]] <- geom$coordinates[[1]]                 # outer ring
    } else if (identical(gtype, "MultiPolygon")) {
      for (pp in seq_along(geom$coordinates)) {
        rings[[length(rings) + 1L]] <- geom$coordinates[[pp]][[1]]
      }
    } else {
      warning("Skipping feature ", k, ": unsupported geometry type '", gtype, "'.")
      next
    }

    for (r in seq_along(rings)) {
      ring <- rings[[r]]
      if (length(ring) < 3) next
      xy <- do.call(rbind, lapply(ring, function(pt) c(as.numeric(pt[[1]]), as.numeric(pt[[2]]))))
      gid <- if (length(rings) > 1) paste0(poly_id, "-", r) else poly_id
      rows[[length(rows) + 1L]] <- data.frame(aoi_id = gid, x = xy[, 1], y = xy[, 2],
                                              stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) stop("No usable polygon coordinates extracted from ", json_file)
  aoi <- do.call(rbind, rows)

  # ---- 2. Map QuPath pixels to the Seurat object's coordinate frame ----
  # Step A (fixed geometry): pixels -> GLOBAL microns. QuPath y grows top->bottom while the global/metadata
  #   frame grows bottom->top, so y is always flipped here: Yg = y_top - y_px*scale_um. x has no flip.
  #   y_top (global microns of QuPath y = 0) = min(y_global_mm)*1000 + one FOV height (fov_size_px*pixel_size_um).
  # Step B (match the object): negate x / y to match how the cell polygons were imported (invert_x / invert_y,
  #   e.g. TransferPolygonsToSeurat's defaults), then scale to the output units (mm by default).
  coord_scaling <- match.arg(coord_scaling)
  div <- switch(coord_scaling, divide_by_1000 = 1000, none = 1, multiply_by_1000 = 1 / 1000)
  scale_um <- pixel_size_um * downsample_factor       # pixels -> microns
  sgn_x <- if (isTRUE(x_invert)) -1 else 1
  sgn_y <- if (isTRUE(y_invert)) -1 else 1

  if (!is.null(fov_pos_file)) {
    fov <- if (is.data.frame(fov_pos_file)) fov_pos_file else utils::read.csv(fov_pos_file)
    if (!fov_y_col %in% names(fov)) {
      stop("Column '", fov_y_col, "' not found in the FOV positions file.")
    }
    if (is.null(y_offset)) {
      y_offset <- min(fov[[fov_y_col]], na.rm = TRUE) * 1000 + fov_size_px * pixel_size_um  # global um
    }
    x_global_um <- aoi$x * scale_um + x_offset          # x_offset: global-micron origin shift (default 0)
    y_global_um <- y_offset - aoi$y * scale_um          # QuPath y-down -> global y-up
    aoi$x <- sgn_x * x_global_um / div
    aoi$y <- sgn_y * y_global_um / div
  } else {
    aoi$x <- sgn_x * (aoi$x * scale_um + x_offset) / div
    aoi$y <- sgn_y * (aoi$y * scale_um) / div
    message("No fov_pos_file supplied: coordinates are in image units and are NOT mapped to the metadata frame.")
  }

  if (!is.null(slide)) aoi$slide <- slide
  rownames(aoi) <- NULL

  # ---- 3. Optionally assign cells to ROIs ----
  if (is.null(cell_polygons)) {
    return(aoi)
  }
  aoi_cells <- AssignCellsToAOIs(aoi = aoi, cell_polygons = cell_polygons, slide = slide, ...)
  list(aoi_polygons = aoi, aoi_cells = aoi_cells)
}


#' Assign each cell to the ROI it mostly overlaps (one slide)
#'
#' Given AOI polygons (e.g. from \code{ReadQuPathROIs}) and cell polygon vertices, returns one row per cell with the
#' \code{aoi_id} of the ROI that covers the majority of the cell's area (or NA if none exceeds \code{min_overlap}).
#' Uses \pkg{sf} for the geometry and reproduces the majority-overlap logic of the reference code.
#'
#' @param aoi data.frame of ROI vertices with columns \code{aoi_id}, x, y (as returned by \code{ReadQuPathROIs}).
#' @param cell_polygons Either a data.frame of cell polygon vertices (one row per vertex) with a cell-ID column and
#'   x/y columns, or a Seurat object. For a Seurat object, coordinates are pulled from the FOV whose name matches
#'   \code{slide} (one FOV per slide, as created by \code{TransferPolygonsToSeurat}); if no \code{slide} is given, all
#'   FOVs are combined. Coordinates MUST be in the same frame/units as \code{aoi}.
#' @param cell_id_col Name of the cell identifier column in \code{cell_polygons}. Default "cell_ID".
#' @param aoi_id_col Name of the ROI identifier column in \code{aoi}. Default "aoi_id".
#' @param coords Length-2 character vector naming the x and y columns. Default c("x", "y").
#' @param min_overlap Minimum fraction of a cell's area that must fall inside an ROI to assign it. Default 0.5.
#' @param slide Optional slide identifier; if given, both inputs are filtered to this slide (using \code{slide_col})
#'   and a \code{slide} column is added to the output.
#' @param slide_col Slide column name used for the optional filtering. Default "slide".
#' @param make_valid Logical. Repair invalid polygons with \code{sf::st_make_valid()} before intersecting, which
#'   prevents GEOS "TopologyException / side location conflict" errors from self-intersecting rings. Default TRUE.
#' @return data.frame with one row per cell: \code{cell_id_col}, \code{aoi_id_col} (NA if unassigned),
#'   \code{overlap_prop} (and \code{slide} if provided).
#' @export
AssignCellsToAOIs <- function(aoi, cell_polygons,
                              cell_id_col = "cell_ID",
                              aoi_id_col = "aoi_id",
                              coords = c("x", "y"),
                              min_overlap = 0.5,
                              slide = NULL,
                              slide_col = "slide",
                              make_valid = TRUE) {

  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required for cell-to-AOI assignment.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required for cell-to-AOI assignment.")
  xcol <- coords[1]; ycol <- coords[2]

  # Pull polygon vertices out of a Seurat object if needed. A Seurat object can hold several FOVs
  # (one per slide, as named by TransferPolygonsToSeurat), and GetTissueCoordinates() on the whole
  # object returns only the default FOV. So when a `slide` is given we extract from the FOV whose
  # name matches it; otherwise we combine all FOVs.
  if (inherits(cell_polygons, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Package 'SeuratObject' is required to extract polygons from a Seurat object.")
    }
    imgs <- SeuratObject::Images(cell_polygons)
    if (length(imgs) == 0) stop("The Seurat object has no spatial FOV/Images to extract coordinates from.")
    pick <- if (!is.null(slide) && as.character(slide) %in% imgs) {
      as.character(slide)
    } else if (length(imgs) == 1L) {
      imgs
    } else {
      NULL
    }
    if (!is.null(pick)) {
      cp <- SeuratObject::GetTissueCoordinates(cell_polygons[[pick]])
    } else {
      if (!is.null(slide)) {
        warning("No FOV named '", slide, "' in the object (Images: ", paste(imgs, collapse = ", "),
                "). Combining all FOVs instead.")
      }
      cp <- do.call(rbind, lapply(imgs, function(im) SeuratObject::GetTissueCoordinates(cell_polygons[[im]])))
    }
    names(cp)[names(cp) %in% c("cell", "cells")] <- cell_id_col
    cell_polygons <- cp
  }

  # Coerce to plain data.frames. This is important for data.table inputs, whose `[` semantics
  # differ from base (e.g. DT[, c("x","y")] returns the character vector, not the columns), which
  # otherwise breaks the complete.cases() subsetting below.
  cell_polygons <- as.data.frame(cell_polygons)
  aoi <- as.data.frame(aoi)

  need_cell <- c(cell_id_col, xcol, ycol)
  if (!all(need_cell %in% names(cell_polygons))) {
    stop("cell_polygons must contain columns: ", paste(need_cell, collapse = ", "))
  }
  if (!all(c(aoi_id_col, xcol, ycol) %in% names(aoi))) {
    stop("aoi must contain columns: ", paste(c(aoi_id_col, xcol, ycol), collapse = ", "))
  }

  # Optional per-slide filtering
  if (!is.null(slide)) {
    if (slide_col %in% names(cell_polygons)) {
      cell_polygons <- cell_polygons[as.character(cell_polygons[[slide_col]]) == as.character(slide), , drop = FALSE]
    }
    if (slide_col %in% names(aoi)) {
      aoi <- aoi[as.character(aoi[[slide_col]]) == as.character(slide), , drop = FALSE]
    }
  }

  aoi <- aoi[stats::complete.cases(aoi[, c(xcol, ycol)]), , drop = FALSE]
  cell_polygons <- cell_polygons[stats::complete.cases(cell_polygons[, c(xcol, ycol)]), , drop = FALSE]
  if (nrow(aoi) == 0) stop("No AOI coordinates to work with (after filtering).")
  if (nrow(cell_polygons) == 0) stop("No cell polygon coordinates to work with (after filtering).")

  # ---- Build sf polygons (points -> combined -> polygon), matching the reference pipeline ----
  suppressWarnings({
    sf_cells <- cell_polygons %>%
      sf::st_as_sf(coords = c(xcol, ycol)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cell_id_col))) %>%
      dplyr::summarize(geometry = sf::st_combine(geometry), .groups = "drop") %>%
      sf::st_cast("POLYGON")

    sf_aoi <- aoi %>%
      sf::st_as_sf(coords = c(xcol, ycol)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(aoi_id_col))) %>%
      dplyr::summarize(geometry = sf::st_combine(geometry), .groups = "drop") %>%
      sf::st_cast("POLYGON")

    # Repair invalid polygons before intersecting. Building a polygon from a cell's / ROI's points
    # can produce a self-intersecting ("bowtie") ring when the vertices are not in clean ring order
    # or contain collinear/duplicate points; GEOS then throws
    # "TopologyException: side location conflict ... input geometry is invalid" during st_intersection.
    if (isTRUE(make_valid)) {
      sf_cells <- sf::st_make_valid(sf_cells)
      sf_aoi   <- sf::st_make_valid(sf_aoi)
    }

    sf_cells$orig_area <- as.numeric(sf::st_area(sf_cells))

    inter <- sf::st_intersection(sf_cells, sf_aoi)
    # Intersections along shared boundaries can yield points/lines (or mixed GEOMETRYCOLLECTIONs);
    # keep only the polygonal parts so areas are well-defined.
    if (nrow(inter) > 0) {
      inter <- sf::st_collection_extract(inter, "POLYGON", warn = FALSE)
    }
    inter$intersect_area <- as.numeric(sf::st_area(inter))
  })

  # ---- Majority-overlap assignment (base R, so it is easy to reason about) ----
  d  <- sf::st_drop_geometry(inter)[, c(cell_id_col, aoi_id_col, "intersect_area")]
  oa <- unique(sf::st_drop_geometry(sf_cells)[, c(cell_id_col, "orig_area")])

  res <- oa
  res[[aoi_id_col]] <- NA
  res$overlap_prop <- 0
  if (nrow(d) > 0) {
    # sum intersection area per (cell, aoi) in case an ROI is split into multiple parts
    agg <- aggregate(list(intersect_area = d$intersect_area),
                     by = d[, c(cell_id_col, aoi_id_col)], FUN = sum)
    # keep, per cell, the ROI with the largest overlap
    agg <- agg[order(agg[[cell_id_col]], -agg$intersect_area), ]
    best <- agg[!duplicated(agg[[cell_id_col]]), ]

    m <- merge(oa, best, by = cell_id_col, all.x = TRUE)
    m$intersect_area[is.na(m$intersect_area)] <- 0
    m$overlap_prop <- ifelse(m$orig_area > 0, m$intersect_area / m$orig_area, 0)
    m[[aoi_id_col]] <- ifelse(m$overlap_prop > min_overlap, m[[aoi_id_col]], NA)
    res <- m
  }

  out_cols <- c(cell_id_col, aoi_id_col, "overlap_prop")
  res <- res[, out_cols, drop = FALSE]
  if (!is.null(slide)) res$slide <- slide
  rownames(res) <- NULL
  res
}
