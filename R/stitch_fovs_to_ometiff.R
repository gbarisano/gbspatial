#!/usr/bin/env Rscript

# ==============================================================================
# stitch_fovs_to_ometiff() : CosMx per-FOV morphology images -> stitched OME-TIFF
#
# Part of the {gbspatial} package (https://github.com/gbarisano/gbspatial).
#
# The OME-XML metadata is generated *programmatically* from the actual stitched
# image, so you no longer need a hand-made .xml file. Only the pieces Minerva
# actually needs are written: correct Pixels dimensions, per-Channel Name/Color,
# and the TiffData IFD map. The bulky CosMx JSON metadata block is intentionally
# omitted (Minerva ignores it).
# ==============================================================================


# ---- internal helpers --------------------------------------------------------

# Escape the five XML-special characters for use inside attribute values.
.gb_xml_escape <- function(x) {
  x <- gsub("&",  "&amp;",  x, fixed = TRUE)
  x <- gsub("<",  "&lt;",   x, fixed = TRUE)
  x <- gsub(">",  "&gt;",   x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'",  "&apos;", x, fixed = TRUE)
  x
}

# Convert any R colour (name, "#RRGGBB", or col2rgb-able) to the signed int32
# RGBA value used by the OME "Color" attribute: (R<<24)|(G<<16)|(B<<8)|A .
.gb_col_to_ome <- function(col, alpha = 255L) {
  rgb <- grDevices::col2rgb(col)
  v <- rgb[1, ] * 2^24 + rgb[2, ] * 2^16 + rgb[3, ] * 2^8 + alpha
  v <- ifelse(v >= 2^31, v - 2^32, v)          # wrap to signed 32-bit
  format(v, scientific = FALSE, trim = TRUE)
}

# Built-in channel presets. `cosmx_protein` reproduces the 5-plex protein
# morphology kit (channel/band order B,G,Y,R,U = PanCK,CD68,Membrane,CD45,DAPI).
.gb_channel_preset <- function(preset = c("cosmx_protein")) {
  preset <- match.arg(preset)
  if (preset == "cosmx_protein") {
    data.frame(
      name       = c("PanCK", "CD68", "Membrane", "CD45", "DAPI"),
      color      = c("green", "yellow", "cyan", "magenta", "white"),
      excitation = c(488, 530, 590, 656, 385),
      emission   = c(512, 553, 630, 684, 512),
      stringsAsFactors = FALSE
    )
  }
}

# Normalise the user-supplied `channels` argument into a data.frame with columns
# name / color / excitation / emission of length n. Accepts:
#   NULL          -> cosmx_protein preset when n==5, otherwise generic
#   "cosmx_protein" (character preset name)
#   character vector of channel names (auto colours)
#   data.frame with at least a `name` column (optional color/excitation/emission)
.gb_resolve_channels <- function(channels, n) {
  generic_cols <- c("white", "green", "magenta", "cyan", "yellow",
                    "red", "blue", "orange", "chartreuse", "deeppink")

  if (is.null(channels)) {
    if (n == 5L) {
      message("No 'channels' supplied; assuming the CosMx 5-plex protein ",
              "morphology kit (PanCK, CD68, Membrane, CD45, DAPI). ",
              "Pass 'channels' to override, and check the order matches the ",
              "band order in your TIFFs.")
      return(.gb_channel_preset("cosmx_protein"))
    }
    message("No 'channels' supplied; using generic names Channel1..", n, ".")
    return(data.frame(
      name       = paste0("Channel", seq_len(n)),
      color      = rep(generic_cols, length.out = n),
      excitation = NA_real_, emission = NA_real_, stringsAsFactors = FALSE))
  }

  if (is.character(channels) && length(channels) == 1L &&
      channels %in% c("cosmx_protein")) {
    channels <- .gb_channel_preset(channels)
  }

  if (is.character(channels)) {
    channels <- data.frame(name = channels, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(channels) || !"name" %in% names(channels)) {
    stop("'channels' must be NULL, a preset name, a character vector of names, ",
         "or a data.frame with at least a 'name' column.")
  }
  if (is.null(channels$color)) {
    channels$color <- rep(generic_cols, length.out = nrow(channels))
  }
  if (is.null(channels$excitation)) channels$excitation <- NA_real_
  if (is.null(channels$emission))   channels$emission   <- NA_real_

  if (nrow(channels) != n) {
    stop(sprintf(paste0("'channels' describes %d channel(s) but the images ",
                        "have %d band(s). They must match (and be in the same ",
                        "order as the TIFF bands)."), nrow(channels), n))
  }
  channels
}

# Build a minimal-but-valid OME-XML string for a single-file, single-resolution
# OME-TIFF: one plane (IFD) per channel, uint16, XYCZT.
.gb_build_ome_xml <- function(size_x, size_y, channels, pixel_size_um,
                              image_name = "stitched") {
  n <- nrow(channels)
  col_int <- .gb_col_to_ome(channels$color)

  chan_xml <- vapply(seq_len(n), function(k) {
    ex <- channels$excitation[k]; em <- channels$emission[k]
    ex_attr <- if (!is.na(ex)) sprintf(' ExcitationWavelength="%s" ExcitationWavelengthUnit="nm"', ex) else ""
    em_attr <- if (!is.na(em)) sprintf(' EmissionWavelength="%s" EmissionWavelengthUnit="nm"', em) else ""
    sprintf(paste0('      <Channel ID="Channel:0:%d" Name="%s" ',
                   'SamplesPerPixel="1" Color="%s"%s%s><LightPath/></Channel>'),
            k - 1L, .gb_xml_escape(channels$name[k]), col_int[k], ex_attr, em_attr)
  }, character(1))

  tiff_xml <- vapply(seq_len(n), function(k) {
    sprintf('      <TiffData FirstC="%d" FirstT="0" FirstZ="0" IFD="%d" PlaneCount="1"/>',
            k - 1L, k - 1L)
  }, character(1))

  paste0(
'<?xml version="1.0" encoding="UTF-8"?>\n',
'<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" ',
'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
'xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 ',
'http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">\n',
'  <Image ID="Image:0" Name="', .gb_xml_escape(image_name), '">\n',
'    <Pixels ID="Pixels:0" DimensionOrder="XYCZT" Type="uint16" ',
'SignificantBits="16" BigEndian="false" Interleaved="false" ',
'SizeX="', size_x, '" SizeY="', size_y, '" SizeC="', n, '" SizeZ="1" SizeT="1" ',
'PhysicalSizeX="', pixel_size_um, '" PhysicalSizeXUnit="\u00b5m" ',
'PhysicalSizeY="', pixel_size_um, '" PhysicalSizeYUnit="\u00b5m">\n',
paste(chan_xml, collapse = "\n"), "\n",
paste(tiff_xml, collapse = "\n"), "\n",
'    </Pixels>\n',
'  </Image>\n',
'</OME>\n')
}

# Locate a bftools executable: prefer an explicit dir, then PATH.
.gb_find_bftool <- function(tool, bftools_dir = NULL) {
  if (!is.null(bftools_dir)) {
    cand <- file.path(bftools_dir, tool)
    if (file.exists(cand)) return(cand)
  }
  p <- Sys.which(tool)
  if (nzchar(p)) return(unname(p))
  ""
}

# Read a CosMx *_fov_positions_file and return data.frame(FOV, x_px, y_px).
# Auto-detects px vs mm columns; converts mm->px with pixel_size_um when needed.
.gb_read_fov_positions <- function(pos, pixel_size_um) {
  if (is.character(pos)) {
    if (!file.exists(pos)) stop("FOV positions file not found: ", pos)
    df <- utils::read.csv(pos)
  } else if (is.data.frame(pos)) {
    df <- pos
  } else stop("'fov_positions' must be a file path or a data.frame.")

  nm <- tolower(names(df))
  pick <- function(cands) { i <- which(nm %in% cands); if (length(i)) i[1] else NA_integer_ }

  fov_i <- pick(c("fov", "fov_id"))
  if (is.na(fov_i)) stop("Could not find a FOV column (FOV / fov) in positions file.")

  xpx_i <- pick(c("x_global_px", "x_px")); ypx_i <- pick(c("y_global_px", "y_px"))
  if (!is.na(xpx_i) && !is.na(ypx_i)) {
    x <- df[[xpx_i]]; y <- df[[ypx_i]]
  } else {
    xmm_i <- pick(c("x_global_mm", "x_mm")); ymm_i <- pick(c("y_global_mm", "y_mm"))
    if (is.na(xmm_i) || is.na(ymm_i))
      stop("Positions file needs x/y in px (x_global_px,y_global_px) or mm ",
           "(x_global_mm,y_global_mm).")
    message("Positions file has mm coordinates; converting to px using ",
            "pixel_size_um = ", pixel_size_um, ".")
    x <- df[[xmm_i]] * 1000 / pixel_size_um
    y <- df[[ymm_i]] * 1000 / pixel_size_um
  }
  data.frame(FOV = as.integer(df[[fov_i]]), x_px = as.numeric(x),
             y_px = as.numeric(y), stringsAsFactors = FALSE)
}

# Match image files in `image_dir` to FOV numbers by parsing "...F00001.TIF".
# Falls back to sorted order (with a warning) if no FOV number can be parsed.
.gb_match_images_to_fovs <- function(image_dir, fov_ids) {
  files <- list.files(image_dir, pattern = "\\.tiff?$", ignore.case = TRUE,
                      full.names = TRUE)
  if (!length(files)) stop("No .tif/.tiff files found in: ", image_dir)

  parsed <- suppressWarnings(as.integer(sub(
    ".*[Ff]0*([0-9]+)\\.[Tt][Ii][Ff][Ff]?$", "\\1", basename(files))))

  if (all(is.na(parsed))) {
    warning("Could not parse FOV numbers from filenames in ", image_dir,
            "; falling back to sorted order. Verify the mapping is correct.")
    files <- files[order(basename(files))]
    if (length(files) != length(fov_ids))
      stop("Number of images (", length(files), ") != number of FOVs (",
           length(fov_ids), ") and filenames carry no FOV number to match on.")
    return(stats::setNames(files, fov_ids))
  }

  idx <- match(fov_ids, parsed)
  if (anyNA(idx)) {
    missing <- fov_ids[is.na(idx)]
    warning("No image found for FOV(s): ", paste(missing, collapse = ", "),
            ". They will be skipped.")
  }
  stats::setNames(files[idx], fov_ids)
}


#' Stitch CosMx per-FOV morphology images into a Minerva-ready OME-TIFF
#'
#' @description
#' Takes the per-FOV multi-channel morphology TIFFs produced by a CosMx run and
#' stitches them, using the FOV global coordinates, into a single (optionally
#' downsampled) multi-channel image per slide. The result is written as an
#' OME-TIFF whose OME-XML metadata is generated programmatically so that
#' \href{https://github.com/labsyspharm/minerva-story}{Minerva} reads channels,
#' colours and physical pixel size correctly. No hand-made metadata XML file is
#' required.
#'
#' @details
#' The stitching engine is \pkg{terra}: each FOV is loaded, optionally
#' down-sampled by averaging (\code{aggregate(fun = "mean")}), given a spatial
#' extent from its global pixel coordinates (Y is flipped so tiles assemble on a
#' Cartesian plane), and the tiles are merged; empty background is set to 0.
#'
#' \strong{Metadata.} Only the parts Minerva needs are emitted: the
#' \code{Pixels} dimensions (of the actual stitched, down-sampled image), one
#' \code{Channel} per band with \code{Name}/\code{Color} (and optional
#' excitation/emission), and the \code{TiffData} IFD map (one plane per
#' channel). \code{PhysicalSizeX/Y} is \code{pixel_size_um * downsample_factor}.
#' The large CosMx JSON acquisition block is deliberately dropped.
#'
#' \strong{Backends} (\code{backend}):
#' \itemize{
#'   \item \code{"bftools"} (default, most robust; supports BigTIFF): \pkg{terra}
#'     writes a tiled multi-band TIFF, \code{bfconvert} re-lays it out as a
#'     proper one-plane-per-channel OME-TIFF, then \code{tiffcomment} injects the
#'     generated OME-XML. Requires Bio-Formats command-line tools
#'     (\code{bfconvert}, \code{tiffcomment}).
#'   \item \code{"rbioformats"} (experimental; no command-line tools, but pulls
#'     in the \pkg{RBioFormats} Bioconductor package and Java): writes the
#'     OME-TIFF directly from R. If the installed API does not match, the
#'     function stops with a message suggesting \code{backend = "bftools"}.
#' }
#'
#' @param fov_positions Path to a CosMx \code{*_fov_positions_file.csv(.gz)} (or
#'   an already-read data.frame), or a (named) list/vector of these for multiple
#'   slides. Must contain a FOV column and either px (\code{x_global_px},
#'   \code{y_global_px}) or mm (\code{x_global_mm}, \code{y_global_mm}) columns.
#' @param image_dir Directory containing the per-FOV morphology TIFFs
#'   (e.g. a \code{Morphology2D} folder), or a list/vector of directories, one
#'   per slide, matching \code{fov_positions}.
#' @param out_dir Output directory. Created if it does not exist.
#' @param slide_names Optional character vector of slide names (used for output
#'   file names \code{<slide>.ome.tif}). Defaults to the names of
#'   \code{fov_positions}, or \code{slide1..N}.
#' @param fov_size_px Integer FOV edge length in pixels of the original images.
#'   Default \code{4256}.
#' @param downsample_factor Integer >= 1. Averages blocks of this size; 8 gives
#'   good resolution at a Minerva-friendly size when you have 4 TMAs for example.
#'   Default \code{1}.
#' @param pixel_size_um Physical size of one original pixel, in microns. Default
#'   \code{0.120280945} (CosMx). The written \code{PhysicalSize} is scaled by
#'   \code{downsample_factor}.
#' @param channels Channel definition. \code{NULL} (default) uses the CosMx
#'   5-plex protein preset when there are 5 bands, otherwise generic names; a
#'   preset name (\code{"cosmx_protein"}); a character vector of names; or a
#'   data.frame with columns \code{name} and optionally \code{color},
#'   \code{excitation}, \code{emission}. The order MUST match the TIFF band
#'   order (CosMx protein morphology is B,G,Y,R,U =
#'   PanCK, CD68, Membrane, CD45, DAPI).
#' @param backend One of \code{"bftools"} (default) or \code{"rbioformats"}.
#' @param bftools_dir Directory holding \code{bfconvert}/\code{tiffcomment}
#'   (e.g. \code{"/Applications/bftools"}). If \code{NULL}, they are looked up on
#'   the \code{PATH}. Only used for \code{backend = "bftools"}.
#' @param keep_intermediate Logical; keep the intermediate plain \code{.tif}
#'   (bftools backend). Default \code{FALSE}.
#' @param overwrite Logical; overwrite existing output files. Default \code{TRUE}.
#'
#' @return (Invisibly) a character vector of the OME-TIFF paths written.
#'
#' @examples
#' \dontrun{
#' # Single slide, CosMx 5-plex protein morphology, 8x downsample:
#' stitch_fovs_to_ometiff(
#'   fov_positions = "mytma/flatFiles/tma1/mytma_fov_positions_file.csv.gz",
#'   image_dir     = "mytma/DecodedFiles/tma1/20240404_004321_S2/CellStatsDir/Morphology2D",
#'   out_dir       = "mytma/ometiff/downsampled8"
#' )
#'
#' # Multiple slides with explicit channel names/colours:
#' stitch_fovs_to_ometiff(
#'   fov_positions = list(tma1 = "flatFiles/tma1/fov_positions_file.csv.gz",
#'                        tma2 = "flatFiles/tma2/fov_positions_file.csv.gz"),
#'   image_dir     = list(tma1 = "DecodedFiles/tma1/20240404_004321_S2/CellStatsDir/Morphology2D",
#'                        tma2 = "DecodedFiles/tma2/20250427_010603_S1/CellStatsDir/Morphology2D"),
#'   out_dir       = "mytma/ometiff/downsampled8",
#'   channels      = data.frame(
#'     name  = c("PanCK","CD68","Membrane","CD45","DAPI"),
#'     color = c("green","yellow","cyan","magenta","white"))
#' )
#' }
#'
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @export
stitch_fovs_to_ometiff <- function(fov_positions,
                                   image_dir,
                                   out_dir,
                                   slide_names       = NULL,
                                   fov_size_px       = 4256L,
                                   downsample_factor = 8L,
                                   pixel_size_um     = 0.120280945,
                                   channels          = NULL,
                                   backend           = c("bftools", "rbioformats"),
                                   bftools_dir       = NULL,
                                   keep_intermediate = FALSE,
                                   overwrite         = TRUE) {

  backend <- match.arg(backend)
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required for stitching. Install it with ",
         "install.packages('terra').")

  downsample_factor <- as.integer(downsample_factor)
  if (is.na(downsample_factor) || downsample_factor < 1L)
    stop("'downsample_factor' must be an integer >= 1.")

  # -- normalise inputs into per-slide lists ----------------------------------
  as_list <- function(x) if (is.data.frame(x) || (is.character(x) && length(x) == 1L)) list(x) else as.list(x)
  fov_positions <- as_list(fov_positions)
  image_dir     <- as_list(image_dir)
  n_slides <- length(fov_positions)
  if (length(image_dir) == 1L && n_slides > 1L)
    image_dir <- rep(image_dir, n_slides)
  if (length(image_dir) != n_slides)
    stop("'image_dir' must be length 1 or match the number of slides in ",
         "'fov_positions'.")

  if (!is.null(slide_names)) {
    if (length(slide_names) != n_slides)
      stop("'slide_names' must match the number of slides.")
  } else if (!is.null(names(fov_positions)) && all(nzchar(names(fov_positions)))) {
    slide_names <- names(fov_positions)
  } else {
    slide_names <- paste0("slide", seq_len(n_slides))
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # -- resolve bftools up front (fail fast) -----------------------------------
  if (backend == "bftools") {
    bfconvert   <- .gb_find_bftool("bfconvert",   bftools_dir)
    tiffcomment <- .gb_find_bftool("tiffcomment", bftools_dir)
    if (!nzchar(bfconvert) || !nzchar(tiffcomment))
      stop("Could not find 'bfconvert'/'tiffcomment'. Set 'bftools_dir' (e.g. ",
           "\"/Applications/bftools\") or add bftools to your PATH. ",
           "Alternatively try backend = \"rbioformats\".")
  }

  eff_pixel_um <- pixel_size_um * downsample_factor
  written <- character(0)

  for (s in seq_len(n_slides)) {
    slide <- slide_names[s]
    message("\n=== Slide '", slide, "' (", s, "/", n_slides, ") ===")

    tiles  <- .gb_read_fov_positions(fov_positions[[s]], pixel_size_um)
    imgmap <- .gb_match_images_to_fovs(image_dir[[s]], tiles$FOV)
    keep   <- !is.na(imgmap)
    tiles  <- tiles[keep, , drop = FALSE]
    imgmap <- imgmap[keep]
    n_tiles <- nrow(tiles)
    if (!n_tiles) stop("No FOV images could be matched for slide '", slide, "'.")
    message("Matched ", n_tiles, " FOV images; stitching (downsample x",
            downsample_factor, ")...")

    # -- load + place each tile -----------------------------------------------
    raster_list <- vector("list", n_tiles)
    pb <- utils::txtProgressBar(min = 0, max = n_tiles, style = 3)
    for (i in seq_len(n_tiles)) {
      r <- terra::rast(unname(imgmap[i]))
      if (downsample_factor > 1L)
        r <- terra::aggregate(r, fact = downsample_factor, fun = "mean")
      x_top <- tiles$x_px[i]; y_top <- tiles$y_px[i]
      terra::ext(r) <- c(x_top, x_top + fov_size_px, -(y_top + fov_size_px), -y_top)
      raster_list[[i]] <- r
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)

    n_ch     <- terra::nlyr(raster_list[[1]])
    channels_df <- .gb_resolve_channels(channels, n_ch)

    message("Merging ", n_tiles, " tiles into one image...")
    merged <- terra::merge(terra::sprc(raster_list))
    merged[is.na(merged)] <- 0

    size_x <- terra::ncol(merged); size_y <- terra::nrow(merged)
    ome_xml <- .gb_build_ome_xml(size_x, size_y, channels_df, eff_pixel_um,
                                 image_name = slide)

    ome_tif <- file.path(out_dir, paste0(slide, ".ome.tif"))
    if (file.exists(ome_tif) && !overwrite)
      stop("Output exists and overwrite = FALSE: ", ome_tif)

    if (backend == "bftools") {
      plain_tif <- file.path(out_dir, paste0(slide, ".tif"))
      message("Writing tiled 16-bit TIFF: ", plain_tif)
      terra::writeRaster(
        merged, plain_tif, datatype = "INT2U", overwrite = TRUE,
        wopt = list(gdal = c("COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=YES",
                             "BLOCKXSIZE=512", "BLOCKYSIZE=512")))

      xml_file <- file.path(out_dir, paste0(slide, ".ome.xml"))
      writeLines(ome_xml, xml_file, useBytes = TRUE)

      message("bfconvert -> OME-TIFF...")
      st <- system2(bfconvert,
                    c("-overwrite", shQuote(path.expand(plain_tif)),
                      shQuote(path.expand(ome_tif))),
                    stdout = FALSE, stderr = FALSE)
      if (!identical(st, 0L)) stop("bfconvert failed for slide '", slide, "'.")

      message("tiffcomment -> inject OME-XML...")
      st <- system2(tiffcomment,
                    c("-set", shQuote(path.expand(xml_file)),
                      shQuote(path.expand(ome_tif))),
                    stdout = FALSE, stderr = FALSE)
      if (!identical(st, 0L)) stop("tiffcomment failed for slide '", slide, "'.")

      if (!keep_intermediate) {
        unlink(plain_tif); unlink(xml_file)
      }

    } else { # rbioformats
      if (!requireNamespace("RBioFormats", quietly = TRUE))
        stop("backend = 'rbioformats' needs the 'RBioFormats' package ",
             "(Bioconductor) and Java. Otherwise use backend = 'bftools'.")
      # Reorder to XYCZT and write via Bio-Formats. Channel names/colours are
      # taken from the generated OME-XML where the installed API allows it.
      arr <- terra::as.array(merged)            # rows x cols x channels
      arr <- aperm(arr, c(2, 1, 3))             # -> X, Y, C
      ok <- tryCatch({
        RBioFormats::write.image(arr, file = path.expand(ome_tif),
                                 force = overwrite, pixelType = "uint16")
        TRUE
      }, error = function(e) {
        stop("RBioFormats::write.image failed (", conditionMessage(e),
             "). Your installed API may differ; use backend = 'bftools'.")
      })
      if (ok)
        message("Wrote OME-TIFF via RBioFormats. NOTE: verify channel names/",
                "colours in Minerva; if missing, re-run with backend='bftools'.")
    }

    message("Done: ", ome_tif)
    written <- c(written, ome_tif)
  }

  invisible(written)
}


# ==============================================================================
# Command Line Execution Block  (run directly:  Rscript stitch_fovs_to_ometiff.R ...)
# ==============================================================================
# sys.nframe() == 0 is TRUE only when the file is executed directly, not sourced.
if (sys.nframe() == 0) {

  suppressWarnings(suppressMessages({
    if (!requireNamespace("optparse", quietly = TRUE))
      stop("The CLI needs the 'optparse' package. install.packages('optparse').")
  }))

  option_list <- list(
    optparse::make_option(c("-p", "--fov_positions"), type = "character", default = NULL,
      help = "Comma-separated path(s) to *_fov_positions_file.csv(.gz), one per slide.", metavar = "paths"),
    optparse::make_option(c("-d", "--image_dir"), type = "character", default = NULL,
      help = "Comma-separated Morphology2D director(y/ies), one per slide (or one for all).", metavar = "dirs"),
    optparse::make_option(c("-o", "--out_dir"), type = "character", default = ".",
      help = "Output directory [default= %default]", metavar = "dir"),
    optparse::make_option(c("-n", "--slide_names"), type = "character", default = NULL,
      help = "Comma-separated slide names (output file stems). Optional.", metavar = "names"),
    optparse::make_option(c("-s", "--fov_size_px"), type = "integer", default = 4256L,
      help = "FOV edge length in pixels [default= %default]", metavar = "int"),
    optparse::make_option(c("-f", "--downsample_factor"), type = "integer", default = 8L,
      help = "Downsampling factor (block-mean) [default= %default]", metavar = "int"),
    optparse::make_option(c("-u", "--pixel_size_um"), type = "double", default = 0.120280945,
      help = "Original pixel size in microns [default= %default]", metavar = "num"),
    optparse::make_option(c("-c", "--channels"), type = "character", default = NULL,
      help = paste("Comma-separated channel NAMES in band order, or a preset",
                   "('cosmx_protein'). Colours auto-assigned unless names match the preset."), metavar = "spec"),
    optparse::make_option(c("-b", "--backend"), type = "character", default = "bftools",
      help = "'bftools' or 'rbioformats' [default= %default]", metavar = "str"),
    optparse::make_option(c("-t", "--bftools_dir"), type = "character", default = NULL,
      help = "Directory with bfconvert/tiffcomment (else use PATH).", metavar = "dir"),
    optparse::make_option(c("-k", "--keep_intermediate"), action = "store_true", default = FALSE,
      help = "Keep the intermediate plain .tif and .ome.xml [default= %default]"),
    optparse::make_option(c("--no_overwrite"), action = "store_true", default = FALSE,
      help = "Do not overwrite existing outputs.")
  )

  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

  if (is.null(opt$fov_positions) || is.null(opt$image_dir)) {
    optparse::print_help(optparse::OptionParser(option_list = option_list))
    stop("Both --fov_positions and --image_dir are required.", call. = FALSE)
  }

  split_csv <- function(x) if (is.null(x)) NULL else trimws(unlist(strsplit(x, ",")))

  fov_positions <- split_csv(opt$fov_positions)
  image_dir     <- split_csv(opt$image_dir)
  slide_names   <- split_csv(opt$slide_names)

  # channels: preset name passes through; otherwise treat as a vector of names.
  ch <- split_csv(opt$channels)
  if (!is.null(ch) && length(ch) == 1L && ch %in% c("cosmx_protein")) ch <- ch[[1]]

  # This assumes stitch_fovs_to_ometiff() is defined above (it is, in this file).
  stitch_fovs_to_ometiff(
    fov_positions     = as.list(fov_positions),
    image_dir         = as.list(image_dir),
    out_dir           = opt$out_dir,
    slide_names       = slide_names,
    fov_size_px       = opt$fov_size_px,
    downsample_factor = opt$downsample_factor,
    pixel_size_um     = opt$pixel_size_um,
    channels          = ch,
    backend           = opt$backend,
    bftools_dir       = opt$bftools_dir,
    keep_intermediate = opt$keep_intermediate,
    overwrite         = !opt$no_overwrite
  )
}
