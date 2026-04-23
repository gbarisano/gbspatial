#!/usr/bin/env Rscript

# Load required libraries
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
#'             the second specifically for show filtering. If NULL, defaults to Idents(seurat_obj).
#' @param image_name The name of the image slot(s) in the Seurat object. If NULL, automatically uses all available images for each FOV.
#' @param fov_col The name of the metadata column that contains the FOV indices. Defaults to "fov".
#' @param img_col The name of the metadata column linking cells to their specific image. Defaults to "Run_Tissue_name".
#' @param boundary The name of the segmentation boundary to use. Defaults to "segmentation" (for CosMx).
#' @param out_dir The directory to save the PNG images. Defaults to current working directory.
#' @param scale_length Length of the scale bar. Defaults to 0.1.
#' @param scale_label Label for the scale bar. Defaults to "100 µm".
#' @param img_dir Comma-separated list of directories (e.g., up to 4) containing background TMA images. Must match the number and order of images.
#' @param img_interval Interval for image raster mapping. Defaults to 0.5119157.
#' @param alpha Transparency for polygons. Defaults to 0.4 if img_dir is used, otherwise 1.0.
#' @param fov_pos_file Comma-separated list of paths to CosMx fov_positions_file.csv.gz. Must match the number and order of images.
#' @param cells_show Comma-separated string of cell names to explicitly plot. Others will be ignored.
#' @param cells_fill Comma-separated string of cell names to fill. Others will be drawn as empty polygons.
#'
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
                                scale_label = "100 µm",
                                img_dir = NULL,
                                img_interval = 0.5119157,
                                alpha = NULL,
                                fov_pos_file = NULL,
                                cells_show = NULL,
                                cells_fill = NULL) {

  # Parse cell subsets if provided
  parsed_cells_show <- if (!is.null(cells_show)) trimws(unlist(strsplit(cells_show, ","))) else NULL
  parsed_cells_fill <- if (!is.null(cells_fill)) trimws(unlist(strsplit(cells_fill, ","))) else NULL

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
    image_name <- trimws(unlist(strsplit(image_name, ",")))
  }
  global_target_images <- if (is.null(image_name)) Images(seurat_obj) else image_name

  # Load FOV positions files mappings
  fov_pos_data_list <- list()
  if (!is.null(fov_pos_file)) {
    fov_pos_files <- trimws(unlist(strsplit(fov_pos_file, ",")))

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
    img_dirs <- trimws(unlist(strsplit(img_dir, ",")))

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
  meta$seurat_ident <- as.character(Idents(seurat_obj))

  # Parse the 'fill' parameter, which can now be 1 or 2 comma-separated values
  if (is.null(fill)) {
    fill_col <- "seurat_ident"
    show_col <- "seurat_ident"
  } else {
    fill_parts <- trimws(unlist(strsplit(fill, ",")))
    fill_col <- fill_parts[1]
    show_col <- if (length(fill_parts) > 1) fill_parts[2] else fill_col
  }

  # Determine if the parsed strings actually match columns
  is_fill_col <- fill_col %in% colnames(meta)
  is_show_col <- show_col %in% colnames(meta)

  # Determine which columns to actually use for filtering operations
  subset_col_fill <- if (is_fill_col) fill_col else "seurat_ident"
  subset_col_show <- if (is_show_col) show_col else "seurat_ident"

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
  for (img in global_target_images) {
    if (!img %in% Images(seurat_obj)) {
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
        img_cells <- Cells(seurat_obj[[img]])
      }
    } else {
      img_cells <- Cells(seurat_obj[[img]])
    }

    # Apply cells_show filter if specified
    if (!is.null(parsed_cells_show)) {
      show_cells <- rownames(meta)[meta[[subset_col_show]] %in% parsed_cells_show]
      img_cells <- intersect(img_cells, show_cells)
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

    for (current_fov in current_img_fovs) {
      message(paste("  -> Rendering FOV:", current_fov))

      # Safely extract exactly the cells for this FOV within the current image
      fov_cells <- rownames(img_meta)[which(as.character(img_meta[[fov_col]]) == as.character(current_fov))]

      if (length(fov_cells) == 0) {
        warning(paste("FOV", current_fov, "has no matching cells for image", img, "- Skipping."))
        next
      }

      # Subset Seurat object using exact cell names
      sub_obj <- tryCatch({
        subset(seurat_obj, cells = fov_cells)
      }, error = function(e) {
        warning(paste("Error subsetting cells for FOV", current_fov, ":", e$message))
        NULL
      })

      # Verify subset has data for the specific image
      if (is.null(sub_obj) || length(Cells(sub_obj)) == 0 || !img %in% Images(sub_obj)) {
        warning(paste("FOV", current_fov, "subset is empty or missing image", img, "- Skipping."))
        next
      }

      DefaultBoundary(sub_obj[[img]]) <- boundary

      coords <- tryCatch({
        GetTissueCoordinates(sub_obj[[img]])
      }, error = function(e) data.frame())

      if (nrow(coords) == 0) {
        warning(paste("No coordinates retrieved for FOV", current_fov, "image", img, "- Skipping."))
        next
      }

      # Explicitly preserve vertex drawing order to prevent polygon tearing/straight-line artifacts
      coords$vertex_order <- seq_len(nrow(coords))

      coords <- coords %>%
        dplyr::left_join(meta %>% select(-c(x,y)), by = "cell") %>%
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
      p <- ggplot()

      # Add Raster Image if requested
      if (has_img) {
        img_raster <- as.raster(jpeg::readJPEG(tma_img_path))
        p <- p + annotation_raster(img_raster,
                                   xmin = plot_xlim[1],
                                   xmax = plot_xlim[2],
                                   ymin = plot_ylim[1],
                                   ymax = plot_ylim[2])
      }

      # Apply variable fill vs static color, factoring in cells_fill option
      if (!is.null(parsed_cells_fill)) {
        coords_fill <- coords %>% dplyr::filter(.data[[subset_col_fill]] %in% parsed_cells_fill)
        coords_empty <- coords %>% dplyr::filter(!(.data[[subset_col_fill]] %in% parsed_cells_fill))

        # Add background (empty) cells first
        if (nrow(coords_empty) > 0) {
          p <- p + geom_polygon(data = coords_empty, aes(x = x, y = y, group = cell),
                                fill = NA, color = "white", linewidth = 0.1, alpha = poly_alpha)
        }

        # Add foreground (filled) cells over the background
        if (nrow(coords_fill) > 0) {
          if (is_fill_col) {
            p <- p + geom_polygon(data = coords_fill, aes(x = x, y = y, group = cell, fill = .data[[fill_col]]),
                                  color = "white", linewidth = 0.1, alpha = poly_alpha)
          } else {
            p <- p + geom_polygon(data = coords_fill, aes(x = x, y = y, group = cell),
                                  fill = fill_col, color = "white", linewidth = 0.1, alpha = poly_alpha)
          }
        }
      } else {
        # Standard plotting all cells normally
        if (is_fill_col) {
          p <- p + geom_polygon(data = coords, aes(x = x, y = y, group = cell, fill = .data[[fill_col]]),
                                color = "white", linewidth = 0.1, alpha = poly_alpha)
        } else {
          p <- p + geom_polygon(data = coords, aes(x = x, y = y, group = cell),
                                fill = fill_col, color = "white", linewidth = 0.1, alpha = poly_alpha)
        }
      }

      # Ensure identical colors are mapped identically across all FOV plots
      if (is_fill_col && !is.null(global_color_map)) {
        p <- p + scale_fill_manual(values = global_color_map, name = fill_col)
      }

      # Add customizations and theme
      p <- p +
        annotate("segment", x = x_start, xend = x_end, y = y_pos, yend = y_pos,
                 color = "white", linewidth = 1.5) +
        annotate("text", x = x_start + (scale_length / 2), y = y_text_pos,
                 label = scale_label, color = "white", size = 4, fontface = "bold") +
        labs(title = paste("FOV:", current_fov, "| Image:", img))

      # Strict geometry clipping
      if (has_img) {
        # Restrict viewport exactly to the raster image bounds
        p <- p + coord_equal(xlim = plot_xlim, ylim = plot_ylim, expand = FALSE)
      } else {
        p <- p + coord_equal()
      }

      p <- p + theme(
        panel.background = element_rect(fill = "black", color = "black"),
        plot.background = element_rect(fill = "black", color = "black"),
        plot.title = element_text(color = "white", hjust = 0.5, size = 16),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "black"),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "black")
      )
      print(p) #it does not slow-down the job on bash/command line, but it does slow down when launching on the R console in R studio because of the rendering of the image in the "Plots" window

      # 6. Save the plot
      filename <- file.path(out_dir, paste0("FOV_", current_fov, "_", img, "_segmentation.png"))
      ggsave(filename = filename, plot = p, width = 8, height = 8, bg = "black", dpi = 300)
      #filename <- file.path(out_dir, paste0("FOV_", current_fov, "_", img, "_segmentation.pdf"))
      #ggsave(filename = filename, plot = p, width = 8, height = 8, bg = "black")
      #print(nrow(coords))
      message(paste("Saved:", filename))
    }
  }

  message("Pipeline complete.")
}

# ==============================================================================
# Command Line Execution Block
# ==============================================================================

# sys.nframe() == 0 checks if the script is being executed directly from the terminal
if (sys.nframe() == 0) {
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Path to the input Seurat object (.rds file)", metavar="character"),
    make_option(c("-f", "--fovs"), type="character", default=NULL,
                help="Comma-separated list of FOVs or intervals to plot (e.g., '1,2,5:10'). If NULL, plots all.", metavar="character"),
    make_option(c("-c", "--fill"), type="character", default=NULL,
                help="Metadata column name(s) or static color for cell fill. Can be comma-separated: 'fill_col,show_col'. Defaults to Idents.", metavar="character"),
    make_option(c("-m", "--image_name"), type="character", default=NULL,
                help="Name of the image slot. If NULL, uses all available.", metavar="character"),
    make_option(c("-v", "--fov_col"), type="character", default="fov",
                help="Metadata column containing FOV indices [default= %default]", metavar="character"),
    make_option(c("-g", "--img_col"), type="character", default="Run_Tissue_name",
                help="Metadata column linking cells to their specific image [default= %default]", metavar="character"),
    make_option(c("-b", "--boundary"), type="character", default="segmentation",
                help="Segmentation boundary to use (e.g. 'segmentation', 'centroids') [default= %default]", metavar="character"),
    make_option(c("-o", "--out_dir"), type="character", default=".",
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("-s", "--scale_length"), type="numeric", default=0.1,
                help="Length of the scale bar [default= %default]", metavar="numeric"),
    make_option(c("-l", "--scale_label"), type="character", default="100 µm",
                help="Label for the scale bar [default= %default]", metavar="character"),
    make_option(c("-d", "--img_dir"), type="character", default=NULL,
                help="Comma-separated directories containing TMA images for background overlay.", metavar="character"),
    make_option(c("-t", "--img_interval"), type="numeric", default=0.5119157,
                help="Interval for the spatial extent of the overlaid image [default= %default]", metavar="numeric"),
    make_option(c("-a", "--alpha"), type="numeric", default=NULL,
                help="Alpha transparency for polygons. Defaults to 0.4 if img_dir is provided, else 1.0.", metavar="numeric"),
    make_option(c("-P", "--fov_pos_file"), type="character", default=NULL,
                help="Comma-separated paths to fov_positions_file.csv.gz. Must match the number and order of images.", metavar="character"),
    make_option(c("-S", "--cells_show"), type="character", default=NULL,
                help="Comma-separated string of cell names to explicitly plot. Others will not be drawn.", metavar="character"),
    make_option(c("-F", "--cells_fill"), type="character", default=NULL,
                help="Comma-separated string of cell names to fill. Others will be drawn as empty grey polygons.", metavar="character"),
    make_option(c("-q", "--subset"), type="character", default=NULL,
                help="Logical expression to subset cells (e.g., 'nCount_RNA > 20'). Variables must exist in metadata.", metavar="character")
  )

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)

  if (is.null(opt$input)){
    print_help(opt_parser)
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
    missing_vars <- setdiff(req_vars, colnames(seurat_obj@meta.data))
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
    cells_show = opt$cells_show,
    cells_fill = opt$cells_fill
  )
}
