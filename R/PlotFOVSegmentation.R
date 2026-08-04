#!/usr/bin/env Rscript

# Load required libraries if ran on terminal through Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(jpeg)
  library(Polychrome)
})

#' Plot Spatial Transcriptomics Segmentation by FOV
#'
#' @param seurat_obj A Seurat spatial object.
#' @param fovs Numeric or numeric vector specifying the FOVs to plot (e.g., 1, 1:5, c(1,3,5)).
#'             If NULL, iterates through all unique FOVs found dynamically per image.
#' @param fill A string indicating either a metadata column (e.g., "cell_type") or a color (e.g., "red").
#'             It can also be a comma-separated list of up to two columns: the first for fill mapping/filtering,
#'             the second specifically for show filtering. If NULL, defaults to Seurat::Idents(seurat_obj).
#' @param image_name The name of the image slot(s) in the Seurat object. If NULL, automatically uses all available images for each FOV.
#' @param fov_col The name of the metadata column that contains the FOV indices. Defaults to "fov".
#' @param img_col The name of the metadata column linking cells to their specific image. Defaults to "Run_Tissue_name".
#' @param boundary The name of the segmentation boundary to use. Defaults to "segmentation" (for CosMx).
#' @param out_dir The directory to save the PNG images. Defaults to current working directory.
#' @param scale_length Length of the scale bar. Defaults to 0.1.
#' @param scale_label Label for the scale bar. Defaults to "100 µm".
#' @param img_dir Comma-separated list or vector of directories (e.g., up to 4) containing background TMA images. Must match the number and order of images.
#' @param img_interval Interval for image raster mapping. Defaults to 0.5119157.
#' @param alpha Transparency for polygons. Defaults to 0.4 if img_dir is used, otherwise 1.0.
#' @param fov_pos_file Comma-separated list or vector of paths to CosMx fov_positions_file.csv.gz. Must match the number and order of images.
#' @param clusters_show Comma-separated string or vector of category values (e.g. cell-type / cluster names) to
#'               explicitly plot. Matched against the show column (the fill column, or its second element if two are
#'               supplied via \code{fill}). Cells not in this set are dropped. Was previously named \code{cells_show}.
#' @param clusters_fill Comma-separated string or vector of category values to fill. Matched against the fill column.
#'               Cells not in this set are drawn as empty (unfilled) outlines. Was previously named \code{cells_fill}.
#' @param cells_show Vector (or comma-separated string) of individual cell IDs to explicitly plot. Cells whose ID is
#'               not in this set are dropped. IDs are matched against \code{id_col} if supplied, otherwise the cell
#'               names (barcodes / rownames), with auto-detection of the best-matching column as a fallback.
#' @param cells_fill Vector (or comma-separated string) of individual cell IDs to fill. Cells whose ID is not in this
#'               set are drawn as empty outlines (unless they are selected by \code{clusters_fill}). Matched like
#'               \code{cells_show}.
#' @param stitch Logical. If TRUE, all selected FOVs belonging to the same image are combined into a single
#'               stitched montage: every FOV's immunofluorescence tile is placed at its own global position and
#'               all segmentation polygons are overlaid on top, producing one PNG per image instead of one per FOV.
#'               Because cell coordinates and FOV image positions are already in a shared global coordinate frame,
#'               providing \code{fov_pos_file} is strongly recommended for accurate tile alignment when \code{img_dir}
#'               is used. Defaults to FALSE (original per-FOV behavior).
#' @param id_col Optional name of the metadata column that \code{cells_show} / \code{cells_fill} IDs refer to
#'               (e.g. a barcode column such as "updated_cellID"). If omitted, IDs are matched against the cell names
#'               (rownames); if those don't match, the column whose values best cover the supplied IDs is auto-detected.
#'               Note: \code{clusters_show} / \code{clusters_fill} are unaffected by \code{id_col}.
#' @param rotate Numeric. Rotation applied to the whole plot (raster tiles and overlaid polygons together), in
#'               degrees counter-clockwise. Any value is allowed (e.g. 90, 180, 270, or an arbitrary angle like 37).
#'               Defaults to 0.
#' @param flip One of "none", "horizontal", "vertical", or "both". Mirrors the whole plot before rotation.
#'             Combined with \code{rotate}, this covers every orientation/reflection. Defaults to "none".
#' @param return_plots Logical. If TRUE, the ggplot object(s) are collected and returned as a named list
#'             (keys like "<image>_FOV_<fov>" per FOV, or "<image>_stitched" when \code{stitch=TRUE}), so they
#'             can be captured with \code{p <- PlotFOVSegmentation(...)} and further modified or re-saved.
#'             Defaults to FALSE (returns invisible NULL, as before).
#' @param save_png Logical. If TRUE, PNG files are written to \code{out_dir} (original behavior). Set FALSE to skip
#'             writing files, e.g. when you only want the returned objects. Defaults to TRUE.
#' @param show_plot Logical. If TRUE, each plot is printed (rendered in the RStudio Plots pane). Set FALSE to avoid
#'             the slow rendering when you only want to capture or save the plots. Defaults to TRUE.
#' @param extra_polygons Optional data.frame of additional polygons to overlay on every plot (e.g. areas of
#'             interest). Provide the coordinates in their raw frame; the function runs them through the same pipeline
#'             as the cells so they stay aligned. Defaults to NULL (no overlay).
#' @param extra_polygons_group Column in \code{extra_polygons} identifying each polygon (the ggplot group).
#'             Defaults to "aoi_id".
#' @param extra_polygons_coords Length-2 character vector naming the x and y columns in \code{extra_polygons}.
#'             Defaults to c("x", "y").
#' @param extra_polygons_color,extra_polygons_fill,extra_polygons_linewidth,extra_polygons_alpha Aesthetics for the
#'             overlay outline/fill. Defaults: red outline, no fill, linewidth 1, alpha 1.
#' @param extra_polygons_invert_x,extra_polygons_invert_y,extra_polygons_divide_by_1000 Coordinate conversion applied
#'             to \code{extra_polygons} BEFORE the shared flip/rotate, to bring them into the object's frame. These
#'             should match what was used when importing the cell polygons; defaults (invert_y = TRUE,
#'             divide_by_1000 = TRUE, invert_x = FALSE) match \code{TransferPolygonsToSeurat}'s defaults.
#' @param extra_polygons_shift Length-2 numeric c(dx, dy) added to the overlay AFTER all transforms, in final plot
#'             units (mm). Use this to nudge out any small residual offset between the AOI and segmentation sources.
#'             Defaults to c(0, 0).
#' @importFrom dplyr %>%
#' @return Invisible NULL. Saves PNG files to the specified directory.
#' @export
PlotFOVSegmentation <- function(seurat_obj,
                                fovs = NULL,
                                fill = NULL,
                                image_name = NULL,
                                fov_col = "fov",
                                img_col = "Run_Tissue_name",
                                boundary = "segmentation",
                                out_dir = ".",
                                scale_length = 0.1,
                                scale_label = "100 \u00b5m",
                                img_dir = NULL,
                                img_interval = 0.5119157,
                                alpha = NULL,
                                fov_pos_file = NULL,
                                clusters_show = NULL,
                                clusters_fill = NULL,
                                cells_show = NULL,
                                cells_fill = NULL,
                                stitch = FALSE,
                                id_col = NULL,
                                rotate = 0,
                                flip = "none",
                                return_plots = FALSE,
                                save_png = TRUE,
                                show_plot = TRUE,
                                extra_polygons = NULL,
                                extra_polygons_group = "aoi_id",
                                extra_polygons_coords = c("x", "y"),
                                extra_polygons_color = "red",
                                extra_polygons_fill = NA,
                                extra_polygons_linewidth = 1,
                                extra_polygons_alpha = 1,
                                extra_polygons_invert_x = FALSE,
                                extra_polygons_invert_y = TRUE,
                                extra_polygons_divide_by_1000 = TRUE,
                                extra_polygons_shift = c(0, 0)) {

  # Helper to handle parameters that can be a single comma-separated string OR a vector
  parse_list_param <- function(param) {
    if (is.null(param)) return(NULL)
    # If it's a single string, assume it might be comma-separated and split it
    if (length(param) == 1 && is.character(param)) {
      return(trimws(unlist(strsplit(param, ","))))
    }
    # Otherwise, treat it as a vector and format appropriately
    return(trimws(as.character(param)))
  }

  # Validate flip and rotate
  flip <- tolower(as.character(flip)[1])
  if (!flip %in% c("none", "horizontal", "vertical", "both")) {
    stop("`flip` must be one of 'none', 'horizontal', 'vertical', 'both'. Got: '", flip, "'.")
  }
  if (!is.numeric(rotate) || length(rotate) != 1) {
    stop("`rotate` must be a single numeric value (degrees, counter-clockwise).")
  }

  # Helper: apply flip-then-rotate to points about a pivot (cx, cy). Rotation is CCW.
  .gb_transform_pts <- function(x, y, angle_deg, flip, cx, cy) {
    th <- angle_deg * pi / 180
    fx <- if (flip %in% c("horizontal", "both")) -1 else 1
    fy <- if (flip %in% c("vertical", "both"))   -1 else 1
    dx <- (x - cx) * fx
    dy <- (y - cy) * fy
    list(x = cx + dx * cos(th) - dy * sin(th),
         y = cy + dx * sin(th) + dy * cos(th))
  }

  # Helper: add one immunofluorescence tile to plot `p`, applying flip + rotation so it
  # stays aligned with the (also-transformed) polygons. `xmin..ymax` are the tile's raw,
  # pre-transform global-mm bounds; (cx, cy) is the shared pivot.
  .gb_add_tile <- function(p, path, xmin, xmax, ymin, ymax, angle_deg, flip, cx, cy) {
    a <- jpeg::readJPEG(path)
    if (length(dim(a)) == 2L) {
      if (flip %in% c("horizontal", "both")) a <- a[, rev(seq_len(ncol(a))), drop = FALSE]
      if (flip %in% c("vertical", "both"))   a <- a[rev(seq_len(nrow(a))), , drop = FALSE]
    } else {
      if (flip %in% c("horizontal", "both")) a <- a[, rev(seq_len(dim(a)[2])), , drop = FALSE]
      if (flip %in% c("vertical", "both"))   a <- a[rev(seq_len(dim(a)[1])), , , drop = FALSE]
    }
    ras <- grDevices::as.raster(a)
    w <- xmax - xmin
    h <- ymax - ymin
    tcx <- (xmin + xmax) / 2
    tcy <- (ymin + ymax) / 2
    # Transformed center: a flip about the shared pivot moves the tile (e.g. a horizontal
    # flip swaps left/right tiles in a stitched montage), so the tile must be re-placed here,
    # not left at its original bounds. For the no-op case (flip='none', angle 0) nc == center.
    nc <- .gb_transform_pts(tcx, tcy, angle_deg, flip, cx, cy)
    if (angle_deg %% 360 == 0) {
      # Pure flip (or no-op): tile stays axis-aligned; draw pixel-flipped at the repositioned box.
      return(p + ggplot2::annotation_raster(ras,
                                            xmin = nc$x - w / 2, xmax = nc$x + w / 2,
                                            ymin = nc$y - h / 2, ymax = nc$y + h / 2))
    }
    grob <- grid::grobTree(grid::rasterGrob(ras, interpolate = TRUE),
                           vp = grid::viewport(angle = angle_deg))
    p + ggplot2::annotation_custom(grob,
                                   xmin = nc$x - w / 2, xmax = nc$x + w / 2,
                                   ymin = nc$y - h / 2, ymax = nc$y + h / 2)
  }

  # Helper: overlay user-supplied polygons (e.g. areas of interest), running them through the
  # SAME pipeline as the cells so they stay aligned: (1) optional invert x / y, (2) optional
  # divide-by-1000 (raw import-style conversion into the object's coordinate frame), (3) the
  # same flip/rotate about the shared pivot (cx, cy), (4) an optional (dx, dy) shift in final
  # plot units. Drawn on top of the cells; clipped by the plot's coord limits like everything else.
  .gb_add_extra_polys <- function(p, ep, angle_deg, flip, cx, cy) {
    if (is.null(ep)) return(p)
    xcol <- extra_polygons_coords[1]; ycol <- extra_polygons_coords[2]
    if (!all(c(xcol, ycol, extra_polygons_group) %in% names(ep))) {
      warning("extra_polygons must contain the columns '", xcol, "', '", ycol, "' and '",
              extra_polygons_group, "'. Skipping the overlay.")
      return(p)
    }
    epx <- if (isTRUE(extra_polygons_invert_x)) -1 else 1
    epy <- if (isTRUE(extra_polygons_invert_y)) -1 else 1
    epd <- if (isTRUE(extra_polygons_divide_by_1000)) 1000 else 1
    df <- data.frame(
      x   = epx * as.numeric(ep[[xcol]]) / epd,
      y   = epy * as.numeric(ep[[ycol]]) / epd,
      grp = ep[[extra_polygons_group]],
      stringsAsFactors = FALSE
    )
    df <- df[stats::complete.cases(df[, c("x", "y")]), , drop = FALSE]
    if (nrow(df) == 0) return(p)
    if (angle_deg %% 360 != 0 || flip != "none") {
      tp <- .gb_transform_pts(df$x, df$y, angle_deg, flip, cx, cy)
      df$x <- tp$x; df$y <- tp$y
    }
    df$x <- df$x + extra_polygons_shift[1]
    df$y <- df$y + extra_polygons_shift[2]
    p + ggplot2::geom_polygon(
      data = df, ggplot2::aes(x = x, y = y, group = grp),
      fill = extra_polygons_fill, color = extra_polygons_color,
      linewidth = extra_polygons_linewidth, alpha = extra_polygons_alpha,
      inherit.aes = FALSE
    )
  }

  do_transform <- (rotate %% 360 != 0) || (flip != "none")

  # Parse cell subsets if provided
  parsed_clusters_show <- parse_list_param(clusters_show)
  parsed_clusters_fill <- parse_list_param(clusters_fill)
  parsed_cells_show <- parse_list_param(cells_show)
  parsed_cells_fill <- parse_list_param(cells_fill)

  # 1. Determine target metadata and images
  actual_cols <- colnames(seurat_obj@meta.data)

  # Check for exact match first, then fallback to case-insensitive match
  if (!fov_col %in% actual_cols) {
    matched_cols <- grep(paste0("^", fov_col, "$"), actual_cols, ignore.case = TRUE, value = TRUE)
    if (length(matched_cols) > 0) {
      fov_col <- matched_cols[1]
      message(paste("Exact column not found. Using case-insensitive match:", fov_col))
    } else {
      stop(paste("Column matching", fov_col, "(case-insensitive) not found in Seurat object metadata."))
    }
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Prepare global target images
  if (!is.null(image_name)) {
    image_name <- parse_list_param(image_name)
  }
  global_target_images <- if (is.null(image_name)) Seurat::Images(seurat_obj) else image_name

  # Load FOV positions files mappings
  fov_pos_data_list <- list()
  if (!is.null(fov_pos_file)) {
    fov_pos_files <- parse_list_param(fov_pos_file)

    # Handle convenience logic where 1 path is passed for all images
    if (length(fov_pos_files) == 1 && length(global_target_images) > 1) {
      message("Only 1 fov_pos_file provided. Using it for all ", length(global_target_images), " images.")
      fov_pos_files <- rep(fov_pos_files, length(global_target_images))
    } else if (length(fov_pos_files) != length(global_target_images)) {
      stop(paste("The number of fov_pos_file provided (", length(fov_pos_files),
                 ") does not match the number of target images (", length(global_target_images), ")."))
    }

    for (i in seq_along(fov_pos_files)) {
      f_path <- fov_pos_files[i]
      img_name <- global_target_images[i]
      if (file.exists(f_path)) {
        message("Loading FOV positions for image '", img_name, "' from: ", f_path)
        # Check if the file is gzipped to prevent embedded null errors
        if (grepl("\\.gz$", f_path, ignore.case = TRUE)) {
          fov_pos_data_list[[img_name]] <- read.csv(gzfile(f_path))
        } else {
          fov_pos_data_list[[img_name]] <- read.csv(f_path)
        }
      } else {
        warning("FOV positions file not found at ", f_path, ". Falling back to cell minimums for image '", img_name, "'.")
      }
    }
  }

  # Load TMA Background image directories mappings
  img_dir_list <- list()
  if (!is.null(img_dir)) {
    img_dirs <- parse_list_param(img_dir)

    # Handle convenience logic where 1 path is passed for all images
    if (length(img_dirs) == 1 && length(global_target_images) > 1) {
      message("Only 1 img_dir provided. Using it for all ", length(global_target_images), " images.")
      img_dirs <- rep(img_dirs, length(global_target_images))
    } else if (length(img_dirs) != length(global_target_images)) {
      stop(paste("The number of img_dir provided (", length(img_dirs),
                 ") does not match the number of target images (", length(global_target_images), ")."))
    }

    for (i in seq_along(img_dirs)) {
      img_name <- global_target_images[i]
      img_dir_list[[img_name]] <- img_dirs[i]
    }
  }

  # 2. Extract metadata and prepare fill variable
  meta <- seurat_obj@meta.data
  # Safely assign cell names without throwing a duplicate column error
  meta$cell <- rownames(meta)

  # Store Idents explicitly so we can use it for subsetting if needed
  meta$seurat_ident <- as.character(Seurat::Idents(seurat_obj))

  # Parse the 'fill' parameter, which can now be 1 or 2 values
  if (is.null(fill)) {
    fill_col <- "seurat_ident"
    show_col <- "seurat_ident"
  } else {
    fill_parts <- parse_list_param(fill)
    fill_col <- fill_parts[1]
    show_col <- if (length(fill_parts) > 1) fill_parts[2] else fill_col
  }

  # Determine if the parsed strings actually match columns
  is_fill_col <- fill_col %in% colnames(meta)
  is_show_col <- show_col %in% colnames(meta)

  # Determine which columns to actually use for filtering operations
  subset_col_fill <- if (is_fill_col) fill_col else "seurat_ident"
  subset_col_show <- if (is_show_col) show_col else "seurat_ident"

  # clusters_show / clusters_fill are matched against these category columns:
  cluster_col_fill <- subset_col_fill
  cluster_col_show <- subset_col_show

  # --- Resolve the cell-identifier column that cells_show / cells_fill refer to. ---
  # These take literal cell IDs (barcodes / a custom cell-ID column). We (a) honor an
  # explicit id_col, else (b) use the cell names (rownames) if they match, else (c) auto-
  # detect the metadata column whose values best cover the supplied IDs.
  find_id_col <- function(requested, meta_df, intended_col) {
    requested <- unique(as.character(requested))
    if (length(requested) == 0) return(intended_col)
    if (intended_col %in% names(meta_df) &&
        sum(requested %in% as.character(meta_df[[intended_col]])) > 0) {
      return(intended_col)
    }
    cand <- names(meta_df)[vapply(meta_df, function(z) is.character(z) || is.factor(z), logical(1))]
    cand <- unique(c("cell", cand))
    best <- intended_col; best_n <- 0
    for (cc in cand) {
      if (!cc %in% names(meta_df)) next
      n <- sum(requested %in% as.character(meta_df[[cc]]))
      if (n > best_n) { best_n <- n; best <- cc }
    }
    if (best_n == 0) return(intended_col)
    if (best != intended_col) {
      message("cells_show/cells_fill IDs did not match column '", intended_col,
              "'. Auto-detected identifier column '", best, "' (", best_n,
              " matches). Set id_col to override.")
    }
    best
  }

  if (!is.null(id_col) && !(id_col %in% colnames(meta))) {
    warning("id_col '", id_col, "' not found in metadata; ignoring it.")
    id_col <- NULL
  }

  if (!is.null(id_col)) {
    id_col_fill <- id_col
    id_col_show <- id_col
  } else {
    id_col_fill <- if (!is.null(parsed_cells_fill)) find_id_col(parsed_cells_fill, meta, "cell") else "cell"
    id_col_show <- if (!is.null(parsed_cells_show)) find_id_col(parsed_cells_show, meta, "cell") else "cell"
  }

  # Diagnostics on the fill selectors so silent "everything empty" is avoided
  if (!is.null(parsed_cells_fill)) {
    n_cells_match <- sum(unique(as.character(parsed_cells_fill)) %in% as.character(meta[[id_col_fill]]))
    if (n_cells_match == 0) {
      warning("cells_fill matched 0 cells in column '", id_col_fill,
              "'. Check the cells_fill IDs or set id_col.")
    } else {
      message("cells_fill: ", n_cells_match, " cells matched on column '", id_col_fill, "'.")
    }
  }
  if (!is.null(parsed_clusters_fill)) {
    n_clu_match <- sum(as.character(meta[[cluster_col_fill]]) %in% parsed_clusters_fill)
    if (n_clu_match == 0) {
      warning("clusters_fill matched 0 cells in column '", cluster_col_fill,
              "'. Check the clusters_fill values.")
    } else {
      message("clusters_fill: ", n_clu_match, " cells matched on column '", cluster_col_fill, "'.")
    }
  }

  if (!is_fill_col) {
    message("'", fill_col, "' not found in metadata. Assuming it is a specified color. Using Idents for cell_fill filtering options.")
  }

  if (!is.null(fill) && length(fill_parts) > 1 && !is_show_col) {
    message("'", show_col, "' not found in metadata. Using Idents for cell_show filtering options.")
  }

  # --- GLOBALLY MAP COLORS ---
  # Generates a strict mapping dictionary for cell types to guarantee identical colors across all FOVs
  global_color_map <- NULL
  if (is_fill_col) {
    # Extract unique variables for the entire Seurat object to ensure all possible classes are covered
    unique_fills <- sort(unique(as.character(na.omit(meta[[fill_col]]))))
    num_cats <- length(unique_fills)

    if (num_cats > 0) {
      # Fetch or generate Polychrome palette
      poly_pal <- tryCatch({
        # Attempt to natively load the 36-color palette
        unname(as.character(Polychrome::palette36.colors(max(36, num_cats))))
      }, error = function(e) {
        # Fallback: Safely dynamically generate colors if the palette function is missing
        suppressWarnings({
          capture.output({
            set.seed(42) # Seed guarantees the generated colors are reproducible
            p_pal <- Polychrome::createPalette(max(36, num_cats), c("#E41A1C", "#377EB8", "#4DAF4A"))
          })
        })
        unname(as.character(p_pal))
      })

      # Bind the exact colors to the exact class names globally
      global_color_map <- setNames(poly_pal[1:num_cats], unique_fills)
    }
  }

  # Set default alpha based on whether image overlay is used
  poly_alpha <- if (is.null(alpha)) {
    if (!is.null(img_dir)) 0.4 else 1.0
  } else {
    alpha
  }

  # Store user-specified FOVs globally to filter later
  user_fovs <- fovs

  # 3. Iterate over each Image FIRST, then the FOVs corresponding to that image
  # Collector for returned ggplot objects (used when return_plots = TRUE)
  plot_list <- list()

  for (img in global_target_images) {
    if (!img %in% Seurat::Images(seurat_obj)) {
      warning(paste("Image", img, "not found in Seurat object - Skipping."))
      next
    }

    message(paste("Processing Image:", img))

    # Fetch exactly the corresponding directory and position data for this image
    current_img_dir <- img_dir_list[[img]]
    fov_pos_data <- fov_pos_data_list[[img]]

    # Identify which cells actually belong to this image
    if (!is.null(img_col) && img_col %in% colnames(meta)) {
      img_cells <- rownames(meta)[which(meta[[img_col]] == img)]
      if (length(img_cells) == 0) {
        warning(paste("No cells found for image", img, "using metadata column '", img_col, "'. Falling back to Seurat Cells()."))
        img_cells <- Seurat::Cells(seurat_obj[[img]])
      }
    } else {
      img_cells <- Seurat::Cells(seurat_obj[[img]])
    }

    # Apply show filters if specified (each narrows the visible set; intersection)
    if (!is.null(parsed_clusters_show)) {
      keep <- rownames(meta)[as.character(meta[[cluster_col_show]]) %in% parsed_clusters_show]
      img_cells <- intersect(img_cells, keep)
    }
    if (!is.null(parsed_cells_show)) {
      keep <- rownames(meta)[as.character(meta[[id_col_show]]) %in% parsed_cells_show]
      img_cells <- intersect(img_cells, keep)
    }

    if (length(img_cells) == 0) {
      warning(paste("No cells found in image", img, "after applying filters - Skipping."))
      next
    }

    img_meta <- meta[img_cells, , drop = FALSE]

    # Dynamically select FOVs specific to this image if none were provided via arguments
    if (is.null(user_fovs)) {
      current_img_fovs <- unique(na.omit(img_meta[[fov_col]]))
      message("  -> Found ", length(current_img_fovs), " unique FOVs in image '", img, "'.")
    } else {
      current_img_fovs <- user_fovs
    }

    # --- STITCH: containers to accumulate this image's tiles and polygons ---
    stitch_coords_list <- list()
    stitch_rasters <- list()
    if (isTRUE(stitch) && !is.null(current_img_dir) && is.null(fov_pos_data)) {
      warning("stitch=TRUE with a background image but no fov_pos_file for image '", img,
              "'. Tiles will be positioned from per-FOV cell minimums, which can produce seams or ",
              "misalignment between FOVs. Provide fov_pos_file for accurate stitching.")
    }

    for (current_fov in current_img_fovs) {
      message(paste("  -> Rendering FOV:", current_fov))

      # Safely extract exactly the cells for this FOV within the current image
      fov_cells <- rownames(img_meta)[which(as.character(img_meta[[fov_col]]) == as.character(current_fov))]

      if (length(fov_cells) == 0) {
        warning(paste("FOV", current_fov, "has no matching cells for image", img, "- Skipping."))
        next
      }

      # Subset Seurat object using exact cell names (base R subset, compatible with Seurat)
      sub_obj <- tryCatch({
        subset(seurat_obj, cells = fov_cells)
      }, error = function(e) {
        warning(paste("Error subsetting cells for FOV", current_fov, ":", e$message))
        NULL
      })

      # Verify subset has data for the specific image
      if (is.null(sub_obj) || length(Seurat::Cells(sub_obj)) == 0 || !img %in% Seurat::Images(sub_obj)) {
        warning(paste("FOV", current_fov, "subset is empty or missing image", img, "- Skipping."))
        next
      }

      SeuratObject::DefaultBoundary(sub_obj[[img]]) <- boundary

      coords <- tryCatch({
        Seurat::GetTissueCoordinates(sub_obj[[img]])
      }, error = function(e) data.frame())

      if (nrow(coords) == 0) {
        warning(paste("No coordinates retrieved for FOV", current_fov, "image", img, "- Skipping."))
        next
      }

      # Explicitly preserve vertex drawing order to prevent polygon tearing/straight-line artifacts
      coords$vertex_order <- seq_len(nrow(coords))

      coords <- coords %>%
        dplyr::left_join(meta %>% dplyr::select(-c(x,y)), by = "cell") %>%
        dplyr::arrange(vertex_order)

      # 4. Determine plot boundaries
      has_img <- FALSE
      if (!is.null(current_img_dir) && dir.exists(current_img_dir)) {
        fov_str <- sprintf("F%05d", as.numeric(current_fov))
        tma_img_path <- NULL

        # Look strictly within the designated directory for this image
        d <- current_img_dir
        d_base <- basename(d)

        # Attempt to build exact path based on directory basename
        if (grepl("CellComposite", d_base, ignore.case = TRUE)) {
          test_path <- file.path(d, paste0("CellComposite_", fov_str, ".jpg"))
        } else if (grepl("CellOverlay", d_base, ignore.case = TRUE)) {
          test_path <- file.path(d, paste0("CellOverlay_", fov_str, ".jpg"))
        } else {
          test_path <- ""
        }

        # Check if explicitly built path exists
        if (file.exists(test_path)) {
          tma_img_path <- test_path
        } else {
          # Fallback: Find any file ending in FXXXXX.(ext) within this directory
          pattern <- paste0(fov_str, "\\.[a-zA-Z0-9]+$")
          found_files <- list.files(d, pattern = pattern, full.names = TRUE, ignore.case = TRUE)

          if (length(found_files) > 0) {
            tma_img_path <- found_files[1]
          }
        }

        if (!is.null(tma_img_path) && file.exists(tma_img_path)) {
          has_img <- TRUE

          if (!is.null(fov_pos_data) && "fov" %in% tolower(colnames(fov_pos_data))) {
            fov_id_col <- grep("(?i)^fov$", colnames(fov_pos_data), value = TRUE)[1]
            fov_row <- fov_pos_data[fov_pos_data[[fov_id_col]] == as.numeric(current_fov), ]

            if (nrow(fov_row) > 0 && "x_global_mm" %in% colnames(fov_row) && "y_global_mm" %in% colnames(fov_row)) {
              x_min_raster <- as.numeric(fov_row[["x_global_mm"]][1])
              y_max_raster <- as.numeric(fov_row[["y_global_mm"]][1])
              y_min_raster <- y_max_raster - img_interval
            } else {
              warning(paste("Could not find FOV ", current_fov, "in the provided FOV position file for ", img, "- Guessing the coordinates, which may lead to segmentations being unaligned to the FOV image"))
              x_min_raster <- min(coords$x, na.rm = TRUE)
              y_min_raster <- min(coords$y, na.rm = TRUE)
              y_max_raster <- y_min_raster + img_interval
            }
          } else {
            x_min_raster <- min(coords$x, na.rm = TRUE)
            y_min_raster <- min(coords$y, na.rm = TRUE)
            y_max_raster <- y_min_raster + img_interval
          }
          x_max_raster <- x_min_raster + img_interval

          plot_xlim <- c(x_min_raster, x_max_raster)
          plot_ylim <- c(y_min_raster, y_max_raster)
        } else {
          warning(paste("TMA image not found for FOV", current_fov, "in directory", d, "- plotting without overlay."))
          plot_xlim <- c(min(coords$x, na.rm = TRUE), max(coords$x, na.rm = TRUE))
          plot_ylim <- c(min(coords$y, na.rm = TRUE), max(coords$y, na.rm = TRUE))
        }
      } else {
        plot_xlim <- c(min(coords$x, na.rm = TRUE), max(coords$x, na.rm = TRUE))
        plot_ylim <- c(min(coords$y, na.rm = TRUE), max(coords$y, na.rm = TRUE))
      }

      # --- STITCH: accumulate this FOV instead of rendering it on its own ---
      if (isTRUE(stitch)) {
        stitch_coords_list[[as.character(current_fov)]] <- coords
        if (has_img) {
          stitch_rasters[[length(stitch_rasters) + 1L]] <- list(
            path = tma_img_path,
            xmin = plot_xlim[1], xmax = plot_xlim[2],
            ymin = plot_ylim[1], ymax = plot_ylim[2]
          )
        }
        next
      }

      # --- Capture raw (pre-transform) tile bounds and apply rotation/flip ---
      img_xmin <- plot_xlim[1]; img_xmax <- plot_xlim[2]
      img_ymin <- plot_ylim[1]; img_ymax <- plot_ylim[2]
      piv_cx <- mean(plot_xlim); piv_cy <- mean(plot_ylim)
      if (do_transform) {
        tp <- .gb_transform_pts(coords$x, coords$y, rotate, flip, piv_cx, piv_cy)
        coords$x <- tp$x; coords$y <- tp$y
        if (has_img) {
          cc <- .gb_transform_pts(c(img_xmin, img_xmax, img_xmax, img_xmin),
                                  c(img_ymin, img_ymin, img_ymax, img_ymax),
                                  rotate, flip, piv_cx, piv_cy)
          plot_xlim <- range(cc$x); plot_ylim <- range(cc$y)
        } else {
          plot_xlim <- range(coords$x, na.rm = TRUE)
          plot_ylim <- range(coords$y, na.rm = TRUE)
        }
      }

      x_range <- plot_xlim[2] - plot_xlim[1]
      y_range <- plot_ylim[2] - plot_ylim[1]

      # Calculate dynamic positions for the scale bar
      if (has_img) {
        # Place scale bar securely inside the bounding box so it isn't clipped
        x_start <- plot_xlim[1] + (x_range * 0.02)
        x_end <- x_start + scale_length
        y_pos <- plot_ylim[1] + (y_range * 0.03)
        y_text_pos <- y_pos + (y_range * 0.02)
      } else {
        # Place scale bar just outside the standard bounding box
        x_start <- plot_xlim[1]
        x_end <- x_start + scale_length
        y_pos <- plot_ylim[1] - (y_range * 0.02)
        y_text_pos <- y_pos + (y_range * 0.02)
      }

      # 5. Build the ggplot (Initialize without global data)
      p <- ggplot2::ggplot()

      # Add Raster Image if requested
      if (has_img) {
        p <- .gb_add_tile(p, tma_img_path,
                          img_xmin, img_xmax, img_ymin, img_ymax,
                          rotate, flip, piv_cx, piv_cy)
      }

      # Determine which cells to fill: union of clusters_fill and cells_fill selections.
      apply_fill_filter <- !is.null(parsed_clusters_fill) || !is.null(parsed_cells_fill)
      if (apply_fill_filter) {
        fill_mask <- rep(FALSE, nrow(coords))
        if (!is.null(parsed_clusters_fill)) {
          fill_mask <- fill_mask | (as.character(coords[[cluster_col_fill]]) %in% parsed_clusters_fill)
        }
        if (!is.null(parsed_cells_fill)) {
          fill_mask <- fill_mask | (as.character(coords[[id_col_fill]]) %in% parsed_cells_fill)
        }
        coords_fill <- coords[fill_mask, , drop = FALSE]
        coords_empty <- coords[!fill_mask, , drop = FALSE]

        # Add background (empty) cells first
        if (nrow(coords_empty) > 0) {
          p <- p + ggplot2::geom_polygon(data = coords_empty, ggplot2::aes(x = x, y = y, group = cell),
                                         fill = NA, color = "white", linewidth = 0.1, alpha = poly_alpha)
        }

        # Add foreground (filled) cells over the background
        if (nrow(coords_fill) > 0) {
          if (is_fill_col) {
            p <- p + ggplot2::geom_polygon(data = coords_fill, ggplot2::aes(x = x, y = y, group = cell, fill = .data[[fill_col]]),
                                           color = "white", linewidth = 0.1, alpha = poly_alpha)
          } else {
            p <- p + ggplot2::geom_polygon(data = coords_fill, ggplot2::aes(x = x, y = y, group = cell),
                                           fill = fill_col, color = "white", linewidth = 0.1, alpha = poly_alpha)
          }
        }
      } else {
        # Standard plotting all cells normally
        if (is_fill_col) {
          p <- p + ggplot2::geom_polygon(data = coords, ggplot2::aes(x = x, y = y, group = cell, fill = .data[[fill_col]]),
                                         color = "white", linewidth = 0.1, alpha = poly_alpha)
        } else {
          p <- p + ggplot2::geom_polygon(data = coords, ggplot2::aes(x = x, y = y, group = cell),
                                         fill = fill_col, color = "white", linewidth = 0.1, alpha = poly_alpha)
        }
      }

      # Ensure identical colors are mapped identically across all FOV plots
      if (is_fill_col && !is.null(global_color_map)) {
        p <- p + ggplot2::scale_fill_manual(values = global_color_map, name = fill_col)
      }

      # Overlay any user-supplied polygons (areas of interest), transformed to match the cells
      p <- .gb_add_extra_polys(p, extra_polygons, rotate, flip, piv_cx, piv_cy)

      # Add customizations and theme
      p <- p +
        ggplot2::annotate("segment", x = x_start, xend = x_end, y = y_pos, yend = y_pos,
                          color = "white", linewidth = 1.5) +
        ggplot2::annotate("text", x = x_start + (scale_length / 2), y = y_text_pos,
                          label = scale_label, color = "white", size = 4, fontface = "bold") +
        ggplot2::labs(title = paste("FOV:", current_fov, "| Image:", img))

      # Strict geometry clipping
      if (has_img) {
        # Restrict viewport exactly to the raster image bounds
        p <- p + ggplot2::coord_equal(xlim = plot_xlim, ylim = plot_ylim, expand = FALSE)
      } else {
        p <- p + ggplot2::coord_equal()
      }

      p <- p + ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "black", color = "black"),
        plot.background = ggplot2::element_rect(fill = "black", color = "black"),
        plot.title = ggplot2::element_text(color = "white", hjust = 0.5, size = 16),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.background = ggplot2::element_rect(fill = "black"),
        legend.text = ggplot2::element_text(color = "white"),
        legend.title = ggplot2::element_text(color = "white"),
        legend.key = ggplot2::element_rect(fill = "black")
      )
      if (isTRUE(show_plot)) {
        print(p) # renders in the RStudio Plots pane (can be slow); set show_plot=FALSE to skip
      }

      # Collect the plot object if requested
      if (isTRUE(return_plots)) {
        plot_list[[paste0(img, "_FOV_", current_fov)]] <- p
      }

      # 6. Save the plot
      if (isTRUE(save_png)) {
        filename <- file.path(out_dir, paste0(img,"_FOV_", current_fov, img, "_segmentation.png"))
        ggplot2::ggsave(filename = filename, plot = p, width = 8, height = 8, bg = "black", dpi = 300)
        message(paste("Saved:", filename))
      }
    }

    # ==================== STITCHED RENDER (once per image) ====================
    if (isTRUE(stitch)) {
      if (length(stitch_coords_list) == 0) {
        warning(paste("No FOVs with drawable cells to stitch for image", img, "- Skipping stitched output."))
      } else {
        # Combine all per-FOV segmentation polygons (already in the shared global frame)
        coords_all <- dplyr::bind_rows(stitch_coords_list)

        # Union bounding box. Prefer the raster tile extents (true FOV frames) when
        # images are present; otherwise fall back to the span of the polygons.
        if (length(stitch_rasters) > 0) {
          plot_xlim <- c(min(vapply(stitch_rasters, function(r) r$xmin, numeric(1))),
                         max(vapply(stitch_rasters, function(r) r$xmax, numeric(1))))
          plot_ylim <- c(min(vapply(stitch_rasters, function(r) r$ymin, numeric(1))),
                         max(vapply(stitch_rasters, function(r) r$ymax, numeric(1))))
          has_img <- TRUE
        } else {
          plot_xlim <- c(min(coords_all$x, na.rm = TRUE), max(coords_all$x, na.rm = TRUE))
          plot_ylim <- c(min(coords_all$y, na.rm = TRUE), max(coords_all$y, na.rm = TRUE))
          has_img <- FALSE
        }

        # Apply rotation/flip about the montage center (shared pivot for tiles + polygons)
        piv_cx <- mean(plot_xlim); piv_cy <- mean(plot_ylim)
        if (do_transform) {
          tp <- .gb_transform_pts(coords_all$x, coords_all$y, rotate, flip, piv_cx, piv_cy)
          coords_all$x <- tp$x; coords_all$y <- tp$y
          if (length(stitch_rasters) > 0) {
            xs <- numeric(0); ys <- numeric(0)
            for (r in stitch_rasters) {
              cc <- .gb_transform_pts(c(r$xmin, r$xmax, r$xmax, r$xmin),
                                      c(r$ymin, r$ymin, r$ymax, r$ymax),
                                      rotate, flip, piv_cx, piv_cy)
              xs <- c(xs, cc$x); ys <- c(ys, cc$y)
            }
            plot_xlim <- range(xs); plot_ylim <- range(ys)
          } else {
            plot_xlim <- range(coords_all$x, na.rm = TRUE)
            plot_ylim <- range(coords_all$y, na.rm = TRUE)
          }
        }

        x_range <- plot_xlim[2] - plot_xlim[1]
        y_range <- plot_ylim[2] - plot_ylim[1]

        # Scale bar placement (mirrors the per-FOV logic)
        if (has_img) {
          x_start <- plot_xlim[1] + (x_range * 0.02)
          x_end <- x_start + scale_length
          y_pos <- plot_ylim[1] + (y_range * 0.03)
          y_text_pos <- y_pos + (y_range * 0.02)
        } else {
          x_start <- plot_xlim[1]
          x_end <- x_start + scale_length
          y_pos <- plot_ylim[1] - (y_range * 0.02)
          y_text_pos <- y_pos + (y_range * 0.02)
        }

        # Build the stitched ggplot
        p <- ggplot2::ggplot()

        # Lay down every FOV's immunofluorescence tile at its own global position
        for (r in stitch_rasters) {
          p <- .gb_add_tile(p, r$path, r$xmin, r$xmax, r$ymin, r$ymax,
                            rotate, flip, piv_cx, piv_cy)
        }

        # Overlay filled cell polygons: union of clusters_fill and cells_fill selections
        apply_fill_filter <- !is.null(parsed_clusters_fill) || !is.null(parsed_cells_fill)
        if (apply_fill_filter) {
          fill_mask <- rep(FALSE, nrow(coords_all))
          if (!is.null(parsed_clusters_fill)) {
            fill_mask <- fill_mask | (as.character(coords_all[[cluster_col_fill]]) %in% parsed_clusters_fill)
          }
          if (!is.null(parsed_cells_fill)) {
            fill_mask <- fill_mask | (as.character(coords_all[[id_col_fill]]) %in% parsed_cells_fill)
          }
          coords_fill  <- coords_all[fill_mask, , drop = FALSE]
          coords_empty <- coords_all[!fill_mask, , drop = FALSE]

          if (nrow(coords_empty) > 0) {
            p <- p + ggplot2::geom_polygon(data = coords_empty, ggplot2::aes(x = x, y = y, group = cell),
                                           fill = NA, color = "white", linewidth = 0.1, alpha = poly_alpha)
          }
          if (nrow(coords_fill) > 0) {
            if (is_fill_col) {
              p <- p + ggplot2::geom_polygon(data = coords_fill, ggplot2::aes(x = x, y = y, group = cell, fill = .data[[fill_col]]),
                                             color = "white", linewidth = 0.1, alpha = poly_alpha)
            } else {
              p <- p + ggplot2::geom_polygon(data = coords_fill, ggplot2::aes(x = x, y = y, group = cell),
                                             fill = fill_col, color = "white", linewidth = 0.1, alpha = poly_alpha)
            }
          }
        } else {
          if (is_fill_col) {
            p <- p + ggplot2::geom_polygon(data = coords_all, ggplot2::aes(x = x, y = y, group = cell, fill = .data[[fill_col]]),
                                           color = "white", linewidth = 0.1, alpha = poly_alpha)
          } else {
            p <- p + ggplot2::geom_polygon(data = coords_all, ggplot2::aes(x = x, y = y, group = cell),
                                           fill = fill_col, color = "white", linewidth = 0.1, alpha = poly_alpha)
          }
        }

        # Keep colors identical to the per-FOV plots via the same global map
        if (is_fill_col && !is.null(global_color_map)) {
          p <- p + ggplot2::scale_fill_manual(values = global_color_map, name = fill_col)
        }

        # Overlay any user-supplied polygons (areas of interest), transformed to match the cells
        p <- .gb_add_extra_polys(p, extra_polygons, rotate, flip, piv_cx, piv_cy)

        # Build a compact, sorted FOV label for the title/filename
        fov_numeric <- suppressWarnings(as.numeric(names(stitch_coords_list)))
        if (any(is.na(fov_numeric))) {
          fov_label <- paste(names(stitch_coords_list), collapse = ", ")
        } else {
          fov_label <- paste(sort(unique(fov_numeric)), collapse = ", ")
        }

        p <- p +
          ggplot2::annotate("segment", x = x_start, xend = x_end, y = y_pos, yend = y_pos,
                            color = "white", linewidth = 1.5) +
          ggplot2::annotate("text", x = x_start + (scale_length / 2), y = y_text_pos,
                            label = scale_label, color = "white", size = 4, fontface = "bold") +
          ggplot2::labs(title = paste0("Stitched FOVs: ", fov_label, "  |  Image: ", img))

        # Clip strictly to the montage bounds when tiles are present
        if (has_img) {
          p <- p + ggplot2::coord_equal(xlim = plot_xlim, ylim = plot_ylim, expand = FALSE)
        } else {
          p <- p + ggplot2::coord_equal()
        }

        p <- p + ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "black", color = "black"),
          plot.background = ggplot2::element_rect(fill = "black", color = "black"),
          plot.title = ggplot2::element_text(color = "white", hjust = 0.5, size = 16),
          axis.line = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(fill = "black"),
          legend.text = ggplot2::element_text(color = "white"),
          legend.title = ggplot2::element_text(color = "white"),
          legend.key = ggplot2::element_rect(fill = "black")
        )

        if (isTRUE(show_plot)) {
          print(p)
        }

        if (isTRUE(return_plots)) {
          plot_list[[paste0(img, "_stitched")]] <- p
        }

        # Size the canvas to the montage aspect ratio so tiles are not letter-boxed
        base_dim <- 12
        aspect <- if (is.finite(y_range) && is.finite(x_range) && x_range > 0) (y_range / x_range) else 1
        out_w <- max(4, min(base_dim, 40))
        out_h <- max(4, min(base_dim * aspect, 40))

        if (isTRUE(save_png)) {
          safe_label <- gsub("[^0-9A-Za-z]+", "-", fov_label)
          filename <- file.path(out_dir, paste0(img, "_stitched_FOVs_", safe_label, "_segmentation.png"))
          ggplot2::ggsave(filename = filename, plot = p, width = out_w, height = out_h,
                          bg = "black", dpi = 300, limitsize = FALSE)
          message(paste("Saved stitched:", filename))
        }
      }
    }
  }

  message("Pipeline complete.")

  if (isTRUE(return_plots)) {
    # If exactly one plot was produced, return it directly for convenience; else a named list.
    if (length(plot_list) == 1L) return(invisible(plot_list[[1]]))
    return(invisible(plot_list))
  }
  invisible(NULL)
}

# ==============================================================================
# Command Line Execution Block
# ==============================================================================

# sys.nframe() == 0 checks if the script is being executed directly from the terminal
if (sys.nframe() == 0) {
  option_list = list(
    optparse::make_option(c("-i", "--input"), type="character", default=NULL,
                          help="Path to the input Seurat object (.rds file)", metavar="character"),
    optparse::make_option(c("-f", "--fovs"), type="character", default=NULL,
                          help="Comma-separated list of FOVs or intervals to plot (e.g., '1,2,5:10'). If NULL, plots all.", metavar="character"),
    optparse::make_option(c("-c", "--fill"), type="character", default=NULL,
                          help="Metadata column name(s) or static color for cell fill. Can be comma-separated: 'fill_col,show_col'. Defaults to Idents.", metavar="character"),
    optparse::make_option(c("-m", "--image_name"), type="character", default=NULL,
                          help="Name of the image slot. If NULL, uses all available.", metavar="character"),
    optparse::make_option(c("-v", "--fov_col"), type="character", default="fov",
                          help="Metadata column containing FOV indices [default= %default]", metavar="character"),
    optparse::make_option(c("-g", "--img_col"), type="character", default="Run_Tissue_name",
                          help="Metadata column linking cells to their specific image [default= %default]", metavar="character"),
    optparse::make_option(c("-b", "--boundary"), type="character", default="segmentation",
                          help="Segmentation boundary to use (e.g. 'segmentation', 'centroids') [default= %default]", metavar="character"),
    optparse::make_option(c("-o", "--out_dir"), type="character", default=".",
                          help="Output directory [default= %default]", metavar="character"),
    optparse::make_option(c("-s", "--scale_length"), type="numeric", default=0.1,
                          help="Length of the scale bar [default= %default]", metavar="numeric"),
    optparse::make_option(c("-l", "--scale_label"), type="character", default="100 \u00b5m",
                          help="Label for the scale bar [default= %default]", metavar="character"),
    optparse::make_option(c("-d", "--img_dir"), type="character", default=NULL,
                          help="Comma-separated directories containing TMA images for background overlay.", metavar="character"),
    optparse::make_option(c("-t", "--img_interval"), type="numeric", default=0.5119157,
                          help="Interval for the spatial extent of the overlaid image [default= %default]", metavar="numeric"),
    optparse::make_option(c("-a", "--alpha"), type="numeric", default=NULL,
                          help="Alpha transparency for polygons. Defaults to 0.4 if img_dir is provided, else 1.0.", metavar="numeric"),
    optparse::make_option(c("-P", "--fov_pos_file"), type="character", default=NULL,
                          help="Comma-separated paths to fov_positions_file.csv.gz. Must match the number and order of images.", metavar="character"),
    optparse::make_option(c("-S", "--clusters_show"), type="character", default=NULL,
                          help="Comma-separated category/cluster values to explicitly plot (matched against the show column). Others not drawn.", metavar="character"),
    optparse::make_option(c("-F", "--clusters_fill"), type="character", default=NULL,
                          help="Comma-separated category/cluster values to fill (matched against the fill column). Others drawn as empty outlines.", metavar="character"),
    optparse::make_option(c("--cells_show"), type="character", default=NULL,
                          help="Comma-separated cell IDs to explicitly plot. Others not drawn. Matched against id_col / cell names.", metavar="character"),
    optparse::make_option(c("--cells_fill"), type="character", default=NULL,
                          help="Comma-separated cell IDs to fill. Others drawn as empty outlines (unless selected by clusters_fill).", metavar="character"),
    optparse::make_option(c("-q", "--subset"), type="character", default=NULL,
                          help="Logical expression to subset cells (e.g., 'nCount_RNA > 20'). Variables must exist in metadata.", metavar="character"),
    optparse::make_option(c("-w", "--stitch"), action="store_true", default=FALSE,
                          help="Stitch all selected FOVs of each image into a single montage plot (one PNG per image) [default= %default]"),
    optparse::make_option(c("-I", "--id_col"), type="character", default=NULL,
                          help="Metadata column that cells_show/cells_fill IDs refer to (e.g. a cell-ID column). Auto-detected if omitted.", metavar="character"),
    optparse::make_option(c("-r", "--rotate"), type="numeric", default=0,
                          help="Rotate the whole plot by this many degrees counter-clockwise [default= %default]", metavar="numeric"),
    optparse::make_option(c("-x", "--flip"), type="character", default="none",
                          help="Mirror the whole plot: 'none', 'horizontal', 'vertical', or 'both' [default= %default]", metavar="character")
  )

  opt_parser = optparse::OptionParser(option_list=option_list)
  opt = optparse::parse_args(opt_parser)

  if (is.null(opt$input)){
    optparse::print_help(opt_parser)
    stop("Input Seurat object file path must be provided (-i/--input).", call.=FALSE)
  }

  message("Loading Seurat object from: ", opt$input)
  seurat_obj <- readRDS(opt$input)

  # Perform user-defined subsetting if specified
  if (!is.null(opt$subset)) {
    message("Applying subset expression: ", opt$subset)

    # Parse the string into an R expression safely
    expr <- tryCatch({
      parse(text = opt$subset)
    }, error = function(e) {
      stop(paste("Invalid subset expression syntax:", e$message))
    })

    # Extract variable names required by the expression
    req_vars <- all.vars(expr)

    # Verify all required variables exist in the Seurat metadata
    missing_vars <- dplyr::setdiff(req_vars, colnames(seurat_obj@meta.data))
    if (length(missing_vars) > 0) {
      stop(paste("Error: The following variables used in your --subset expression were not found in the Seurat object metadata:\n",
                 paste(missing_vars, collapse = ", ")))
    }

    # Evaluate the logical expression against the metadata context
    keep_cells <- tryCatch({
      eval(expr, envir = seurat_obj@meta.data)
    }, error = function(e) {
      stop(paste("Error evaluating subset expression:", e$message))
    })

    # Handle NA values logically resolving to FALSE
    keep_cells[is.na(keep_cells)] <- FALSE

    if (sum(keep_cells) == 0) {
      stop("Error: Subsetting expression filtered out all cells (0 cells remaining).")
    }

    # Re-subset the Seurat object and keep strictly matching cells
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep_cells])
    message("  -> Cells remaining after subset: ", ncol(seurat_obj))
  }

  # Parse FOVs from string (handles comma-separated and intervals like "1:5, 8, 10:12")
  if (!is.null(opt$fovs)) {
    fov_parts <- unlist(strsplit(opt$fovs, ","))
    fov_list <- lapply(fov_parts, function(x) {
      x <- trimws(x)
      if (grepl(":", x)) {
        bounds <- as.numeric(unlist(strsplit(x, ":")))
        return(seq(bounds[1], bounds[2]))
      } else {
        return(as.numeric(x))
      }
    })
    opt$fovs <- unique(unlist(fov_list))
  }

  PlotFOVSegmentation(
    seurat_obj = seurat_obj,
    fovs = opt$fovs,
    fill = opt$fill,
    image_name = opt$image_name,
    fov_col = opt$fov_col,
    img_col = opt$img_col,
    boundary = opt$boundary,
    out_dir = opt$out_dir,
    scale_length = opt$scale_length,
    scale_label = opt$scale_label,
    img_dir = opt$img_dir,
    img_interval = opt$img_interval,
    alpha = opt$alpha,
    fov_pos_file = opt$fov_pos_file,
    clusters_show = opt$clusters_show,
    clusters_fill = opt$clusters_fill,
    cells_show = opt$cells_show,
    cells_fill = opt$cells_fill,
    stitch = opt$stitch,
    id_col = opt$id_col,
    rotate = opt$rotate,
    flip = opt$flip
  )
}
