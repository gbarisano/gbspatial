#!/usr/bin/env Rscript

# ==============================================================================
# stitch_fovs_to_ometiff() : CosMx per-FOV morphology images -> stitched OME-TIFF
#
# Part of the {gbspatial} package (https://github.com/gbarisano/gbspatial).
#
# Features
#   * Stitches per-FOV morphology TIFFs into one image per slide (terra).
#   * Optionally combines several slides into a single image on a user-defined
#     grid (rows x cols, input order = fill order) with extra downsampling.
#   * Optionally processes a co-registered H&E image per slide, forced onto the
#     SAME spatial grid as the IF output so the two OME-TIFFs overlay exactly.
#   * OME-XML is generated programmatically (correct dimensions, downsample-
#     scaled physical pixel size, channel names/colours). No external XML file.
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

# Convert any R colour to the signed int32 RGBA used by the OME "Color"
# attribute: (R<<24)|(G<<16)|(B<<8)|A .
.gb_col_to_ome <- function(col, alpha = 255L) {
  rgb <- grDevices::col2rgb(col)
  v <- rgb[1, ] * 2^24 + rgb[2, ] * 2^16 + rgb[3, ] * 2^8 + alpha
  v <- ifelse(v >= 2^31, v - 2^32, v)
  format(v, scientific = FALSE, trim = TRUE)
}

# Built-in channel presets (order = TIFF band order).
.gb_channel_preset <- function(preset = c("cosmx_protein")) {
  preset <- match.arg(preset)
  if (preset == "cosmx_protein") {
    data.frame(
      name       = c("PanCK", "CD68", "Membrane", "CD45", "DAPI"),
      color      = c("green", "yellow", "cyan", "magenta", "white"),
      excitation = c(488, 530, 590, 656, 385),
      emission   = c(512, 553, 630, 684, 512),
      stringsAsFactors = FALSE)
  }
}

# Normalise the `channels` argument into a data.frame of length n.
.gb_resolve_channels <- function(channels, n) {
  generic_cols <- c("white", "green", "magenta", "cyan", "yellow",
                    "red", "blue", "orange", "chartreuse", "deeppink")
  if (is.null(channels)) {
    if (n == 5L) {
      message("No 'channels' supplied; assuming the CosMx 5-plex protein ",
              "morphology kit (PanCK, CD68, Membrane, CD45, DAPI). Pass ",
              "'channels' to override; the order must match the TIFF bands.")
      return(.gb_channel_preset("cosmx_protein"))
    }
    message("No 'channels' supplied; using generic names Channel1..", n, ".")
    return(data.frame(name = paste0("Channel", seq_len(n)),
                      color = rep(generic_cols, length.out = n),
                      excitation = NA_real_, emission = NA_real_,
                      stringsAsFactors = FALSE))
  }
  if (is.character(channels) && length(channels) == 1L &&
      channels %in% c("cosmx_protein")) channels <- .gb_channel_preset(channels)
  if (is.character(channels))
    channels <- data.frame(name = channels, stringsAsFactors = FALSE)
  if (!is.data.frame(channels) || !"name" %in% names(channels))
    stop("'channels' must be NULL, a preset name, a character vector of names, ",
         "or a data.frame with at least a 'name' column.")
  if (is.null(channels$color))
    channels$color <- rep(generic_cols, length.out = nrow(channels))
  if (is.null(channels$excitation)) channels$excitation <- NA_real_
  if (is.null(channels$emission))   channels$emission   <- NA_real_
  if (nrow(channels) != n)
    stop(sprintf("'channels' describes %d channel(s) but the images have %d ",
                 "band(s); they must match and be in the same order.",
                 nrow(channels), n))
  channels
}

# Minimal valid OME-XML for a multi-channel (one plane per channel) uint16 image.
.gb_build_ome_xml <- function(size_x, size_y, channels, pixel_size_um,
                              image_name = "stitched", big_endian = FALSE) {
  n <- nrow(channels); col_int <- .gb_col_to_ome(channels$color)
  chan_xml <- vapply(seq_len(n), function(k) {
    ex <- channels$excitation[k]; em <- channels$emission[k]
    ex_a <- if (!is.na(ex)) sprintf(' ExcitationWavelength="%s" ExcitationWavelengthUnit="nm"', ex) else ""
    em_a <- if (!is.na(em)) sprintf(' EmissionWavelength="%s" EmissionWavelengthUnit="nm"', em) else ""
    sprintf(paste0('      <Channel ID="Channel:0:%d" Name="%s" SamplesPerPixel="1" ',
                   'Color="%s"%s%s><LightPath/></Channel>'),
            k - 1L, .gb_xml_escape(channels$name[k]), col_int[k], ex_a, em_a)
  }, character(1))
  tiff_xml <- vapply(seq_len(n), function(k)
    sprintf('      <TiffData FirstC="%d" FirstT="0" FirstZ="0" IFD="%d" PlaneCount="1"/>',
            k - 1L, k - 1L), character(1))
  paste0(
'<?xml version="1.0" encoding="UTF-8"?>\n',
'<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" ',
'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
'xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 ',
'http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">\n',
'  <Image ID="Image:0" Name="', .gb_xml_escape(image_name), '">\n',
'    <Pixels ID="Pixels:0" DimensionOrder="XYCZT" Type="uint16" SignificantBits="16" ',
'BigEndian="', if (isTRUE(big_endian)) "true" else "false", '" Interleaved="false" ',
'SizeX="', size_x, '" SizeY="', size_y, '" SizeC="', n, '" SizeZ="1" SizeT="1" ',
'PhysicalSizeX="', pixel_size_um, '" PhysicalSizeXUnit="\u00b5m" ',
'PhysicalSizeY="', pixel_size_um, '" PhysicalSizeYUnit="\u00b5m">\n',
paste(chan_xml, collapse = "\n"), "\n", paste(tiff_xml, collapse = "\n"), "\n",
'    </Pixels>\n  </Image>\n</OME>\n')
}

# Minimal valid OME-XML for an interleaved 8-bit RGB image (H&E).
.gb_build_ome_xml_rgb <- function(size_x, size_y, pixel_size_um,
                                  image_name = "he", big_endian = FALSE) {
  paste0(
'<?xml version="1.0" encoding="UTF-8"?>\n',
'<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" ',
'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
'xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 ',
'http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">\n',
'  <Image ID="Image:0" Name="', .gb_xml_escape(image_name), '">\n',
'    <Pixels ID="Pixels:0" DimensionOrder="XYCZT" Type="uint8" SignificantBits="8" ',
'BigEndian="', if (isTRUE(big_endian)) "true" else "false", '" Interleaved="true" ',
'SizeX="', size_x, '" SizeY="', size_y, '" SizeC="3" SizeZ="1" SizeT="1" ',
'PhysicalSizeX="', pixel_size_um, '" PhysicalSizeXUnit="\u00b5m" ',
'PhysicalSizeY="', pixel_size_um, '" PhysicalSizeYUnit="\u00b5m">\n',
'      <Channel ID="Channel:0:0" Name="RGB" SamplesPerPixel="3"><LightPath/></Channel>\n',
'      <TiffData FirstC="0" FirstT="0" FirstZ="0" IFD="0" PlaneCount="1"/>\n',
'    </Pixels>\n  </Image>\n</OME>\n')
}

# Locate a bftools executable.
.gb_find_bftool <- function(tool, bftools_dir = NULL) {
  if (!is.null(bftools_dir)) {
    cand <- file.path(bftools_dir, tool)
    if (file.exists(cand)) return(cand)
  }
  p <- Sys.which(tool); if (nzchar(p)) return(unname(p)); ""
}

# TRUE if a TIFF/BigTIFF file is big-endian ("MM" header), FALSE if little ("II").
.gb_tiff_bigendian <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  identical(readBin(con, "raw", 2L), as.raw(c(0x4d, 0x4d)))
}

# Preflight for backend = "rbioformats": verify RBioFormats, rJava, and a
# working JVM are all present, and fail upfront (before any stitching) if not.
# Critically, the JVM max heap (-Xmx) must be set BEFORE the JVM starts and
# cannot be changed afterwards, so this: (1) checks packages are installed
# WITHOUT loading RBioFormats (which would boot the JVM with a small default
# heap), (2) sets java.parameters -Xmx, then (3) starts the JVM and loads
# RBioFormats.
.gb_check_rbioformats <- function(java_heap = NULL) {
  if (!nzchar(system.file(package = "RBioFormats")))
    stop("backend = 'rbioformats' requires the 'RBioFormats' package, which is ",
         "not installed. Install it from Bioconductor with ",
         "BiocManager::install('RBioFormats'), or use backend = 'bftools'.",
         call. = FALSE)
  if (!nzchar(system.file(package = "rJava")))
    stop("backend = 'rbioformats' needs Java via the 'rJava' package, which is ",
         "not installed. Install a Java runtime and install.packages('rJava'), ",
         "or use backend = 'bftools'.", call. = FALSE)

  # Set the JVM max heap before the JVM boots (only if the user has not already
  # specified -Xmx). This must happen before rJava/RBioFormats start Java.
  heap_mb <- .gb_java_heap_mb(java_heap)
  jp <- getOption("java.parameters")
  if (is.null(jp) || !any(grepl("-Xmx", jp)))
    options(java.parameters = c(jp, paste0("-Xmx", heap_mb, "m")))

  if (!requireNamespace("rJava", quietly = TRUE))
    stop("Failed to load 'rJava'. Use backend = 'bftools'.", call. = FALSE)
  jvm_ok <- tryCatch({ rJava::.jinit(); TRUE },
                     error = function(e) FALSE, warning = function(w) FALSE)
  if (!jvm_ok)
    stop("backend = 'rbioformats' needs a working Java runtime (JVM) but it ",
         "could not be initialised. Install/repair Java (e.g. set JAVA_HOME and ",
         "run 'R CMD javareconf'), or use backend = 'bftools'.", call. = FALSE)
  if (!isTRUE(tryCatch({ loadNamespace("RBioFormats"); TRUE },
                       error = function(e) FALSE)))
    stop("'RBioFormats' is installed but failed to load (usually a Java/",
         "Bio-Formats problem). Use backend = 'bftools'.", call. = FALSE)

  # Report the heap the JVM actually got; if far below what we asked for, the
  # JVM was already running and our setting was ignored -> tell the user.
  actual_mb <- tryCatch({
    rt <- rJava::.jcall("java/lang/Runtime", "Ljava/lang/Runtime;", "getRuntime")
    rJava::.jcall(rt, "J", "maxMemory") / 1024^2
  }, error = function(e) NA_real_)
  if (is.finite(actual_mb)) {
    message(sprintf("  Java max heap ~%.0f MB (requested %d MB).", actual_mb, heap_mb))
    if (actual_mb < 0.8 * heap_mb)
      warning("The JVM was already running, so the larger heap could not be ",
              "applied. If RBioFormats runs out of Java heap, restart R and put ",
              "options(java.parameters = \"-Xmx", heap_mb, "m\") at the very top ",
              "of your script, before loading ANY package -- or use ",
              "backend = 'bftools'.", call. = FALSE)
  }
  invisible(actual_mb)
}

# Total physical RAM in MB (NA if it cannot be determined).
.gb_system_ram_mb <- function() {
  os <- Sys.info()[["sysname"]]
  tryCatch({
    if (identical(os, "Darwin")) {
      as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE)) / 1024^2
    } else if (identical(os, "Linux")) {
      ln <- grep("MemTotal", readLines("/proc/meminfo"), value = TRUE)
      as.numeric(sub(".*?([0-9]+) kB.*", "\\1", ln)) / 1024
    } else if (identical(os, "Windows")) {
      out <- suppressWarnings(system2("wmic",
        c("computersystem", "get", "TotalPhysicalMemory"), stdout = TRUE))
      as.numeric(gsub("\\D", "", paste(out, collapse = ""))) / 1024^2
    } else NA_real_
  }, error = function(e) NA_real_)
}

# Resolve the requested JVM heap to MB. Accepts numeric MB, "24g"/"18000m"
# strings, or NULL (auto: ~70% of RAM, 8 GB fallback).
.gb_java_heap_mb <- function(java_heap) {
  if (!is.null(java_heap)) {
    if (is.numeric(java_heap)) return(max(512L, as.integer(java_heap)))
    val  <- as.numeric(regmatches(java_heap, regexpr("[0-9.]+", java_heap)))
    unit <- tolower(gsub("[0-9. ]", "", java_heap))
    mb <- switch(unit, g = , gb = val * 1024, k = , kb = val / 1024, val) # default MB
    return(max(512L, as.integer(mb)))
  }
  ram <- .gb_system_ram_mb()
  if (is.na(ram)) return(8192L)
  as.integer(max(2048, floor(ram * 0.7)))
}

# Estimate peak RAM (GB) to materialise an in-memory copy of a raster for
# RBioFormats: one double array from as.array() plus an aperm copy plus the
# Java handoff ~= 2.5x the double array.
.gb_rbf_peak_gb <- function(rast) {
  terra::ncell(rast) * terra::nlyr(rast) * 8 / 1024^3 * 2.5
}

# Raise R's vector-memory limit (macOS mem.maxVSize) to cover `need_mb` if the
# current cap is lower. Returns the previous limit (Mb) to restore, or NA if
# nothing was changed / the function is unavailable. Never lowers the cap.
.gb_raise_vsize <- function(need_mb, max_vsize = NULL) {
  if (!exists("mem.maxVSize", mode = "function")) return(NA_real_)
  cur <- tryCatch(mem.maxVSize(), error = function(e) NA_real_)
  if (!is.finite(cur)) return(NA_real_)              # already unlimited
  want <- if (!is.null(max_vsize)) max_vsize else ceiling(need_mb) + 2048
  if (!is.finite(want) || want > cur) {
    old <- tryCatch(mem.maxVSize(want), error = function(e) NA_real_)
    message(sprintf("  Raised R vector-memory limit %.0f -> %s MB (mem.maxVSize).",
                    cur, if (is.finite(want)) format(want) else "Inf"))
    return(cur)
  }
  NA_real_
}

# Inject an OME-XML comment into an existing OME-TIFF using Bio-Formats'
# TiffSaver via rJava -- the same operation the 'tiffcomment' CLI performs, but
# in-process (no external tool). Bio-Formats classes ship inside RBioFormats.
# Returns TRUE on success, FALSE on any failure (so callers can fall back).
.gb_inject_comment_java <- function(ome_tif, xml) {
  tryCatch({
    path <- path.expand(ome_tif)
    in_stream <- rJava::.jnew("loci/common/RandomAccessInputStream", path)
    saver     <- rJava::.jnew("loci/formats/tiff/TiffSaver", path)
    on.exit({
      try(rJava::.jcall(in_stream, "V", "close"), silent = TRUE)
      try(rJava::.jcall(saver,     "V", "close"), silent = TRUE)
    }, add = TRUE)
    jval <- rJava::.jcast(rJava::.jnew("java/lang/String", xml), "java/lang/Object")
    # overwriteComment(RandomAccessInputStream in, Object value): reads the
    # header itself to handle endianness / BigTIFF.
    rJava::.jcall(saver, "V", "overwriteComment",
                  rJava::.jcast(in_stream, "loci/common/RandomAccessInputStream"),
                  jval)
    TRUE
  }, error = function(e) {
    message("  in-process Java metadata injection failed: ", conditionMessage(e))
    FALSE
  })
}

# TRUE if the Bio-Formats classes needed for .gb_inject_comment_java are on the
# (RBioFormats-provided) classpath. Assumes the JVM is already initialised.
.gb_java_injector_available <- function() {
  cls_ok <- function(cls) isTRUE(tryCatch({
    rJava::.jfindClass(gsub("\\.", "/", cls)); TRUE
  }, error = function(e) FALSE))
  cls_ok("loci.formats.tiff.TiffSaver") &&
    cls_ok("loci.common.RandomAccessInputStream")
}

# Read a CosMx *_fov_positions_file -> data.frame(FOV, x_px, y_px).
.gb_read_fov_positions <- function(pos, pixel_size_um) {
  if (is.character(pos)) {
    if (!file.exists(pos)) stop("FOV positions file not found: ", pos)
    df <- utils::read.csv(pos)
  } else if (is.data.frame(pos)) df <- pos
  else stop("'fov_positions' must be a file path or a data.frame.")
  nm <- tolower(names(df))
  pick <- function(c) { i <- which(nm %in% c); if (length(i)) i[1] else NA_integer_ }
  fov_i <- pick(c("fov", "fov_id"))
  if (is.na(fov_i)) stop("Could not find a FOV column in the positions file.")
  xpx <- pick(c("x_global_px", "x_px")); ypx <- pick(c("y_global_px", "y_px"))
  if (!is.na(xpx) && !is.na(ypx)) { x <- df[[xpx]]; y <- df[[ypx]] }
  else {
    xmm <- pick(c("x_global_mm", "x_mm")); ymm <- pick(c("y_global_mm", "y_mm"))
    if (is.na(xmm) || is.na(ymm))
      stop("Positions file needs x/y in px or mm columns.")
    message("Converting mm coordinates to px (pixel_size_um = ", pixel_size_um, ").")
    x <- df[[xmm]] * 1000 / pixel_size_um; y <- df[[ymm]] * 1000 / pixel_size_um
  }
  data.frame(FOV = as.integer(df[[fov_i]]), x_px = as.numeric(x),
             y_px = as.numeric(y), stringsAsFactors = FALSE)
}

# Match image files to FOV numbers by parsing "...F00001.TIF".
.gb_match_images_to_fovs <- function(image_dir, fov_ids) {
  files <- list.files(image_dir, pattern = "\\.tiff?$", ignore.case = TRUE,
                      full.names = TRUE)
  if (!length(files)) stop("No .tif/.tiff files found in: ", image_dir)
  parsed <- suppressWarnings(as.integer(sub(
    ".*[Ff]0*([0-9]+)\\.[Tt][Ii][Ff][Ff]?$", "\\1", basename(files))))
  if (all(is.na(parsed))) {
    warning("Could not parse FOV numbers in ", image_dir,
            "; falling back to sorted order.")
    files <- files[order(basename(files))]
    if (length(files) != length(fov_ids))
      stop("Image count != FOV count and filenames carry no FOV number.")
    return(stats::setNames(files, fov_ids))
  }
  idx <- match(fov_ids, parsed)
  if (anyNA(idx))
    warning("No image for FOV(s): ", paste(fov_ids[is.na(idx)], collapse = ", "),
            "; skipped.")
  stats::setNames(files[idx], fov_ids)
}

# Build one merged (and per-FOV-downsampled) IF SpatRaster for a slide.
.gb_merge_slide <- function(fov_positions, image_dir, pixel_size_um,
                            fov_size_px, downsample_factor, out_dir) {
  tiles  <- .gb_read_fov_positions(fov_positions, pixel_size_um)
  imgmap <- .gb_match_images_to_fovs(image_dir, tiles$FOV)
  keep   <- !is.na(imgmap); tiles <- tiles[keep, , drop = FALSE]
  imgmap <- imgmap[keep]; n_tiles <- nrow(tiles)
  if (!n_tiles) stop("No FOV images could be matched.")
  message("Matched ", n_tiles, " FOV images; loading (downsample x",
          downsample_factor, ")...")
  rl <- vector("list", n_tiles)
  pb <- utils::txtProgressBar(min = 0, max = n_tiles, style = 3)
  for (i in seq_len(n_tiles)) {
    r <- terra::rast(unname(imgmap[i]))
    if (downsample_factor > 1L)
      r <- terra::aggregate(r, fact = downsample_factor, fun = "mean")
    x0 <- tiles$x_px[i]; y0 <- tiles$y_px[i]
    terra::ext(r) <- c(x0, x0 + fov_size_px, -(y0 + fov_size_px), -y0)
    rl[[i]] <- r; utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  tmp <- tempfile(tmpdir = out_dir, fileext = ".tif")
  message("Merging ", n_tiles, " tiles...")
  m <- terra::merge(terra::sprc(rl), filename = tmp, datatype = "INT2U",
                    NAflag = 0, overwrite = TRUE)
  attr(m, "gb_tmp") <- tmp
  m
}

# Grid offsets (dx, dy) for slide placement; row-major fill in input order.
.gb_grid_offsets <- function(rasters, nrow, ncol, gap = 0) {
  n <- length(rasters)
  w <- vapply(rasters, function(r) terra::xmax(r) - terra::xmin(r), numeric(1))
  h <- vapply(rasters, function(r) terra::ymax(r) - terra::ymin(r), numeric(1))
  ri <- ((seq_len(n) - 1L) %/% ncol) + 1L
  ci <- ((seq_len(n) - 1L) %%  ncol) + 1L
  col_w <- vapply(seq_len(ncol), function(c) max(w[ci == c], 0), numeric(1))
  row_h <- vapply(seq_len(nrow), function(r) max(h[ri == r], 0), numeric(1))
  x_start <- c(0, cumsum(col_w + gap))[seq_len(ncol)]
  y_start <- c(0, cumsum(row_h + gap))[seq_len(nrow)]
  dx <- dy <- numeric(n)
  for (k in seq_len(n)) {
    dx[k] <- x_start[ci[k]]  - terra::xmin(rasters[[k]])
    dy[k] <- -y_start[ri[k]] - terra::ymax(rasters[[k]])
  }
  list(dx = dx, dy = dy)
}

# Resolve slide_layout -> c(nrow, ncol).
.gb_resolve_layout <- function(n, layout) {
  if (is.null(layout)) return(c(1L, n))
  layout <- as.integer(layout)
  if (length(layout) == 1L) { nr <- layout; nc <- ceiling(n / nr) }
  else { nr <- layout[1]; nc <- layout[2] }
  if (nr * nc < n)
    stop(sprintf("slide_layout %d x %d cannot hold %d slides.", nr, nc, n))
  c(nr, nc)
}

# ---- polygon helpers ---------------------------------------------------------

# Extract one polygon input (Seurat object or data.frame) into a long table with
# columns cell, x, y (one row per polygon vertex), slide, State.
#   force_slide    : if not NA, assign every row to this slide (positional
#                    per-slide input) and ignore any slide column.
#   fallback_slide : slide to use only when there is no slide column (single
#                    combined input covering one slide).
# Otherwise the slide is read per-cell/row from `slide_col`.
.gb_extract_polys_one <- function(obj, state_col, slide_col, cell_col,
                                  x_col, y_col, force_slide, fallback_slide) {
  assign_slide <- function(col_vals, n) {
    if (!is.na(force_slide)) return(rep(force_slide, n))
    if (!is.null(col_vals))  return(as.character(col_vals))
    if (!is.na(fallback_slide)) return(rep(fallback_slide, n))
    stop("Could not determine each cell's slide: the slide column '", slide_col,
         "' (poly_slide_col) was not found, and there is no per-slide input to ",
         "fall back on. Set 'poly_slide_col' to the correct column, or pass one ",
         "polygon input per slide.", call. = FALSE)
  }
  if (inherits(obj, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE) &&
        !requireNamespace("Seurat", quietly = TRUE))
      stop("A Seurat object was supplied but neither 'SeuratObject' nor ",
           "'Seurat' is installed.", call. = FALSE)
    md <- obj@meta.data
    if (!state_col %in% names(md))
      stop("state_col '", state_col, "' is not a column of the Seurat ",
           "meta.data. Available (first 20): ",
           paste(utils::head(names(md), 20), collapse = ", "), call. = FALSE)
    md$.cell <- rownames(md)
    getimg <- if (requireNamespace("SeuratObject", quietly = TRUE)) SeuratObject::Images else Seurat::Images
    getcoord <- if (requireNamespace("SeuratObject", quietly = TRUE)) SeuratObject::GetTissueCoordinates else Seurat::GetTissueCoordinates
    ims <- tryCatch(getimg(obj), error = function(e) character(0))
    if (!length(ims))
      stop("The Seurat object has no spatial images/FOVs, so segmentation ",
           "polygons cannot be extracted. Provide a polygon table instead ",
           "(see 'polygons').", call. = FALSE)
    coords <- list()
    for (im in ims) {
      gc <- tryCatch(getcoord(obj[[im]], which = "segmentation"),
                     error = function(e) NULL)
      if (!is.null(gc)) coords[[im]] <- as.data.frame(gc)
    }
    if (!length(coords))
      stop("GetTissueCoordinates(..., which = 'segmentation') returned nothing; ",
           "your Seurat/CosMx version may store cell boundaries differently. ",
           "Provide a polygon table instead.", call. = FALSE)
    cd <- do.call(rbind, coords)
    nm <- tolower(names(cd))
    cx <- which(nm == "x")[1]; cy <- which(nm == "y")[1]
    cc <- which(nm %in% c("cell", "cell_id", "cellid"))[1]
    if (anyNA(c(cx, cy, cc)))
      stop("Unexpected columns from GetTissueCoordinates: ",
           paste(names(cd), collapse = ", "),
           " (need x, y, cell).", call. = FALSE)
    out <- data.frame(cell = as.character(cd[[cc]]), x = cd[[cx]], y = cd[[cy]],
                      stringsAsFactors = FALSE)
    col_vals <- if (slide_col %in% names(md))
      md[[slide_col]][match(out$cell, md$.cell)] else NULL
    out$slide <- assign_slide(col_vals, nrow(out))
    out$State <- md[[state_col]][match(out$cell, md$.cell)]
    return(out[!is.na(out$State), , drop = FALSE])
  }
  if (is.data.frame(obj)) {
    need <- c(cell_col, x_col, y_col, state_col)
    miss <- setdiff(need, names(obj))
    if (length(miss))
      stop("Polygon table is missing column(s): ", paste(miss, collapse = ", "),
           ". Available (first 30): ", paste(utils::head(names(obj), 30), collapse = ", "),
           ". Set poly_cell_col / poly_x_col / poly_y_col / state_col.", call. = FALSE)
    col_vals <- if (slide_col %in% names(obj)) obj[[slide_col]] else NULL
    return(data.frame(cell = as.character(obj[[cell_col]]),
                      x = obj[[x_col]], y = obj[[y_col]],
                      slide = assign_slide(col_vals, nrow(obj)),
                      State = obj[[state_col]], stringsAsFactors = FALSE))
  }
  stop("Each 'polygons' input must be a Seurat object or a data.frame (or an ",
       "RDS path to one).", call. = FALSE)
}

# Load + validate all polygon inputs. Returns list(verts = named-by-slide list
# of data.frame(poly_id, x, y in vertex order), state_map = data.frame(poly_id,
# State)). Called during preflight so problems surface before any stitching.
.gb_load_polygons <- function(polygons, slide_names, state_col, slide_col,
                              cell_col, x_col, y_col) {
  if (is.null(state_col) || !nzchar(state_col))
    stop("Polygon generation needs 'state_col': the column to use as the ",
         "Minerva State (e.g. a cluster or cell-type column). Set it, or set ",
         "generate_polygons = FALSE.", call. = FALSE)
  n_slides <- length(slide_names)
  inputs <- polygons
  if (inherits(inputs, "Seurat") || is.data.frame(inputs)) inputs <- list(inputs)
  if (is.character(inputs))
    inputs <- lapply(inputs, function(p) {
      if (!file.exists(p)) stop("Polygon file not found: ", p, call. = FALSE)
      readRDS(p)
    })
  if (!is.list(inputs))
    stop("'polygons' must be a Seurat object, a data.frame, an RDS file path, ",
         "or a list of these.", call. = FALSE)
  n_in <- length(inputs)
  if (!n_in %in% c(1L, n_slides))
    stop("Provide either 1 combined 'polygons' input or exactly ", n_slides,
         " (one per slide); got ", n_in, ".", call. = FALSE)

  long <- vector("list", n_in)
  for (i in seq_len(n_in)) {
    if (n_in == 1L) {
      # One combined input: split by the slide column and keep only the slides
      # being processed. If there is a single slide and no slide column, treat
      # the whole input as that slide.
      force_slide    <- NA_character_
      fallback_slide <- if (n_slides == 1L) slide_names[1] else NA_character_
    } else {
      # One input per slide: assign positionally (input i -> slide_names[i]).
      force_slide    <- slide_names[i]
      fallback_slide <- NA_character_
    }
    long[[i]] <- .gb_extract_polys_one(inputs[[i]], state_col, slide_col,
                                       cell_col, x_col, y_col,
                                       force_slide, fallback_slide)
  }
  df <- do.call(rbind, long)
  df <- df[!is.na(df$slide), , drop = FALSE]
  have <- unique(df$slide)
  miss <- setdiff(slide_names, have)
  if (length(miss) == n_slides)
    stop("None of the requested slides (", paste(slide_names, collapse = ", "),
         ") were found in the polygon slide column '", slide_col, "'. Found: ",
         paste(utils::head(have, 10), collapse = ", "),
         ". Check 'poly_slide_col' and that the names match 'slide_names'.",
         call. = FALSE)
  if (length(miss))
    warning("No polygons found for slide(s): ", paste(miss, collapse = ", "),
            ".", call. = FALSE)
  df <- df[df$slide %in% slide_names, , drop = FALSE]
  if (!nrow(df)) stop("No polygon vertices remain after filtering to the ",
                      "requested slides.", call. = FALSE)

  # Global unique poly_id per (slide, cell), preserving vertex order.
  key <- paste(df$slide, df$cell, sep = "\r")
  df$poly_id <- match(key, unique(key))
  state_map <- df[!duplicated(df$poly_id), c("poly_id", "State")]
  verts <- lapply(slide_names, function(s)
    df[df$slide == s, c("poly_id", "x", "y")])
  names(verts) <- slide_names
  list(verts = verts, state_map = state_map)
}

# Build a SpatVector of polygons from a data.frame(poly_id, x, y). Vertices are
# grouped by poly_id (ascending); cells with < 3 vertices are dropped.
.gb_build_poly_vect <- function(d) {
  d <- d[order(d$poly_id, seq_len(nrow(d))), , drop = FALSE]  # stable, ascending
  tab <- table(d$poly_id)
  ok_ids <- as.numeric(names(tab))[tab >= 3L]
  d <- d[d$poly_id %in% ok_ids, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  v <- terra::vect(as.matrix(d[, c("poly_id", "x", "y")]), type = "polygons")
  terra::values(v) <- data.frame(poly_id = sort(unique(d$poly_id)))
  v
}

# Transform vertex coordinates (per slide), shift, rasterize filled cells +
# borders onto `template`, and write <stem>_polygons.tif and _polygons.csv.
.gb_write_polygons <- function(verts, state_map, template, slide_set, shifts,
                               out_dir, stem, pixel_size_um, fov_size_px,
                               units, flip_y, x_off, y_off, border_value,
                               border_label) {
  factor <- if (units == "mm") 1000 / pixel_size_um else 1
  vects <- list()
  for (s in slide_set) {
    d <- verts[[s]]
    if (is.null(d) || !nrow(d)) next
    x <- round(d$x * factor) + x_off
    yy <- round(d$y * factor)
    y <- (if (isTRUE(flip_y)) -yy else yy) + y_off
    sh <- shifts[[s]]; if (is.null(sh)) sh <- c(0, 0)
    v <- .gb_build_poly_vect(data.frame(poly_id = d$poly_id,
                                        x = x + sh[1], y = y + sh[2]))
    if (!is.null(v)) vects[[length(vects) + 1L]] <- v
  }
  if (!length(vects)) { warning("No polygons to rasterize for '", stem, "'."); return(NULL) }
  master <- if (length(vects) == 1L) vects[[1]] else do.call(rbind, vects)

  tmpl <- terra::rast(template[[1]])          # empty raster, same geometry
  message("  rasterizing ", nrow(master), " cell polygons for '", stem, "' ...")
  tmp_fill <- tempfile(tmpdir = out_dir, fileext = ".tif")
  base_mask <- terra::rasterize(master, tmpl, field = "poly_id", background = 0,
                                filename = tmp_fill, overwrite = TRUE)
  borders <- terra::as.lines(master)
  borders$ub <- border_value
  tmp_upd <- tempfile(tmpdir = out_dir, fileext = ".tif")
  final_mask <- terra::rasterize(borders, base_mask, field = "ub", update = TRUE,
                                 filename = tmp_upd, overwrite = TRUE)
  # Dedicated integer write: rasterize(update=TRUE) inherits the base mask's
  # float type and ignores 'datatype', which yields a float32 TIFF that Minerva
  # cannot convert. writeRaster() with datatype = "INT4U" forces an unsigned
  # 32-bit integer TIFF.
  poly_tif <- file.path(out_dir, paste0(stem, "_polygons.tif"))
  terra::writeRaster(final_mask, poly_tif, datatype = "INT4U", overwrite = TRUE,
                     wopt = list(gdal = c("COMPRESS=LZW", "BIGTIFF=YES")))
  unlink(c(tmp_fill, tmp_upd))

  cs <- state_map[order(state_map$poly_id), , drop = FALSE]
  cs <- data.frame(CellID = cs$poly_id, State = cs$State, stringsAsFactors = FALSE)
  cs <- rbind(cs, data.frame(CellID = border_value, State = border_label))
  poly_csv <- file.path(out_dir, paste0(stem, "_polygons.csv"))
  utils::write.csv(cs, poly_csv, row.names = FALSE, quote = FALSE)
  message("  wrote ", basename(poly_tif), " and ", basename(poly_csv))
  c(poly_tif, poly_csv)
}


#' Stitch CosMx per-FOV morphology images into a Minerva-ready OME-TIFF
#'
#' @description
#' Stitches the per-FOV morphology TIFFs of a CosMx run into one multi-channel
#' image per slide and writes it as an OME-TIFF whose metadata Minerva reads
#' correctly (channel names/colours, physical pixel size). Optionally combines
#' several slides onto a grid, and optionally processes a co-registered H&E
#' image per slide so that its OME-TIFF shares the exact spatial grid and
#' physical size of the IF output.
#'
#' @details
#' Stitching uses \pkg{terra}: each FOV is down-sampled by block-mean, placed by
#' its global pixel coordinates (Y flipped), and the tiles are merged.
#'
#' When \code{stitch_slides = TRUE}, the per-slide mosaics are shifted onto a
#' grid (\code{slide_layout} rows x cols, filled in input order), merged and
#' optionally down-sampled again by \code{combine_downsample}; a single combined
#' OME-TIFF is written.
#'
#' When \code{he_images} is supplied, each H&E image (assumed co-registered to
#' the corresponding IF) is flipped (\code{he_flip}), given the IF slide's
#' coordinate frame, carried through the same grid placement, then
#' \strong{resampled onto the exact IF output grid}. The H&E OME-TIFF therefore
#' has identical \code{SizeX/SizeY} and \code{PhysicalSize} to the IF output
#' (guaranteeing overlay in Minerva) but correct 8-bit RGB channel metadata.
#'
#' Metadata is generated programmatically; no external OME-XML file is needed.
#' Backends: \code{"bftools"} (default; needs \code{bfconvert}/\code{tiffcomment},
#' supports BigTIFF) or \code{"rbioformats"} (experimental, no CLI tools).
#'
#' @param fov_positions Path to a CosMx \code{*_fov_positions_file.csv(.gz)} (or
#'   a data.frame), or a (named) list/vector of these, one per slide.
#' @param image_dir Morphology2D directory (or list/vector, one per slide, or a
#'   single directory reused for all slides).
#' @param out_dir Output directory (created if needed).
#' @param slide_names Optional slide names (also the output stems). If omitted,
#'   taken from the names of \code{fov_positions} when it is a named list,
#'   otherwise derived from the folder containing each positions file (the CosMx
#'   per-slide flatFiles folder, usually equal to \code{Run_Tissue_name}), so
#'   they line up with the polygon slide column; falls back to \code{slide1..N}.
#' @param fov_size_px FOV edge length in pixels. Default \code{4256}.
#' @param downsample_factor Per-FOV block-mean downsampling (>=1). Default \code{8}.
#' @param pixel_size_um Original pixel size in microns. Default \code{0.120280945}.
#'   Written \code{PhysicalSize} = \code{pixel_size_um * downsample_factor}
#'   (times \code{combine_downsample} for the combined image).
#' @param channels Channel definition: \code{NULL} (CosMx 5-plex preset when 5
#'   bands, else generic), a preset name (\code{"cosmx_protein"}), a character
#'   vector of names, or a data.frame with \code{name} and optional
#'   \code{color}/\code{excitation}/\code{emission}. Order must match TIFF bands.
#' @param stitch_slides Logical. If \code{TRUE} and more than one slide is given,
#'   combine them into a single image. Default \code{FALSE}.
#' @param slide_layout Grid for combining: \code{NULL} (single row), a single
#'   integer (number of rows; columns inferred), or \code{c(nrow, ncol)}. Slides
#'   fill the grid row by row in input order.
#' @param combine_downsample Extra block-mean downsampling applied to the
#'   combined image (>=1). Default \code{1}.
#' @param slide_gap Gap between slides on the grid, in pixels of the per-FOV-
#'   downsampled frame. Default \code{0} (abutting).
#' @param combine_name Output stem for the combined image(s). Default
#'   \code{"combined"}.
#' @param he_images Optional co-registered H&E image(s): a path or list/vector of
#'   paths, one per slide (length must equal the number of slides). RGB.
#' @param he_flip How to flip H&E to match the IF Y-orientation: \code{"vertical"}
#'   (default), \code{"horizontal"}, \code{"both"}, or \code{"none"}.
#' @param he_resample_method \pkg{terra} resampling method used to put H&E on the
#'   IF grid. Default \code{"average"} (good for downsampling).
#' @param generate_polygons Logical; if \code{TRUE} (default) and \code{polygons}
#'   is supplied, also write Minerva cell-polygon outputs: a \code{.tif} whose
#'   cells are filled with an integer id (cell borders burned in as
#'   \code{poly_border_value}) and a CSV mapping \code{CellID} to \code{State}
#'   (plus a border row). No OME-TIFF is made for these -- Minerva converts the
#'   plain \code{.tif} itself. If \code{polygons} is \code{NULL}, polygon output
#'   is skipped with a message.
#' @param polygons Cell-boundary input: a \pkg{Seurat} object, a data.frame with
#'   one row per polygon vertex, an \code{.rds} path to either, or a list of
#'   these (one combined input, or one per slide). For a Seurat object, vertices
#'   come from \code{GetTissueCoordinates(obj, which = "segmentation")} and the
#'   slide/State from \code{meta.data}; for a table, from the columns named
#'   below.
#' @param state_col Column holding the Minerva State (e.g. a cluster/cell-type
#'   column) in the Seurat \code{meta.data} or the polygon table. Required when
#'   generating polygons.
#' @param poly_slide_col Slide-identifier column in a combined input, used to
#'   assign each cell to a slide and keep only the slide(s) being processed (so a
#'   multi-slide Seurat object run on one slide contributes only that slide's
#'   cells); its values must match \code{slide_names}. When you instead pass one
#'   input per slide, cells are assigned positionally and this column is ignored.
#'   Default \code{"Run_Tissue_name"}.
#' @param poly_cell_col,poly_x_col,poly_y_col For table input: the cell-id, x and
#'   y column names. Defaults \code{"cell_ID"}, \code{"x_slide_mm"},
#'   \code{"y_slide_mm"}.
#' @param poly_coord_units Units of the polygon x/y: \code{"mm"} (default,
#'   converted with \code{1000 / pixel_size_um}) or \code{"px"}.
#' @param poly_flip_y Logical; negate y to match the raster's flipped frame.
#'   Default \code{TRUE}.
#' @param poly_x_offset_px,poly_y_offset_px Constant offsets (in original px)
#'   added after unit conversion/flip to align polygons to the raster origin.
#'   Defaults \code{0} and \code{-fov_size_px} (\code{poly_y_offset_px = NULL}
#'   resolves to \code{-fov_size_px}).
#' @param poly_border_value,poly_border_label Integer burned into the \code{.tif}
#'   for cell borders and its CSV label. Defaults \code{9999999} and
#'   \code{"Cell Border"}.
#' @param backend \code{"bftools"} (default) or \code{"rbioformats"}. Both build
#'   an OME-TIFF and then inject the generated OME-XML with \code{tiffcomment}
#'   (so channel names/colours/physical size are embedded identically and
#'   Minerva can read the file); \code{tiffcomment} is therefore required for
#'   both. The \code{"bftools"} path builds the container with \code{bfconvert},
#'   streams to disk, and uses roughly constant memory; \code{"rbioformats"}
#'   builds it with \pkg{RBioFormats} but must hold the whole image in RAM, so
#'   prefer \code{"bftools"} for large/combined images.
#' @param metadata_injector How to embed the OME-XML: \code{"auto"} (default)
#'   uses in-process Bio-Formats (\pkg{RBioFormats}, no external tool) when the
#'   backend is \code{"rbioformats"} and the classes are available, otherwise
#'   \code{tiffcomment}; \code{"java"} forces the in-process route (rbioformats
#'   only, no fallback); \code{"tiffcomment"} forces the \code{tiffcomment} CLI.
#'   With \code{backend = "rbioformats"}, \code{metadata_injector = "auto"} and
#'   no \code{he_images}, no bftools tools are required at all.
#' @param bftools_dir Directory with \code{bfconvert}/\code{tiffcomment}
#'   (e.g. \code{"/Applications/bftools"}); otherwise found on \code{PATH}.
#' @param keep_intermediate Keep intermediate plain \code{.tif}/\code{.xml}.
#'   Default \code{FALSE}.
#' @param max_vsize Only used by \code{backend = "rbioformats"} on macOS: the R
#'   vector-memory limit in MB to set via \code{\link[base]{mem.maxVSize}}
#'   (or \code{Inf} to remove it). \code{NULL} (default) auto-raises the cap to
#'   the estimated need when the current cap is lower; the previous value is
#'   restored afterwards. Raising the cap only helps if physical RAM allows it.
#' @param java_heap Only used by \code{backend = "rbioformats"}: the JVM maximum
#'   heap (\code{-Xmx}) as MB (numeric) or a string like \code{"24g"}/
#'   \code{"18000m"}. \code{NULL} (default) uses ~70\% of system RAM. This is
#'   applied only if the JVM has not already started in the session and you have
#'   not set \code{options(java.parameters)} yourself; otherwise set that option
#'   before loading any package. The JVM heap cannot be resized once started.
#' @param overwrite Overwrite existing outputs. Default \code{TRUE}.
#'
#' @return (Invisibly) a character vector of the OME-TIFF paths written.
#'
#' @examples
#' \dontrun{
#' # One slide, IF only:
#' stitch_fovs_to_ometiff("usc1/.../fov_positions_file.csv.gz",
#'                        "~/Morphology2D/usc1", "~/out")
#'
#' # Three slides stitched into a single row, extra 2x downsample, with H&E:
#' stitch_fovs_to_ometiff(
#'   fov_positions = list(usc1 = "usc1/.../fov_positions_file.csv.gz",
#'                        usc2 = "usc2/.../fov_positions_file.csv.gz",
#'                        usc5 = "usc5/.../fov_positions_file.csv.gz"),
#'   image_dir     = list("~/Morphology2D/usc1", "~/Morphology2D/usc2",
#'                        "~/Morphology2D/usc5"),
#'   out_dir       = "~/out",
#'   stitch_slides = TRUE, slide_layout = c(1, 3), combine_downsample = 2,
#'   he_images     = c("~/reg/usc1/he_via_fullIF_smooth.tif",
#'                     "~/reg/usc2/he_via_fullIF_smooth.tif",
#'                     "~/reg/usc5/he_via_fullIF_smooth.tif"))
#' }
#'
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @export
stitch_fovs_to_ometiff <- function(fov_positions,
                                   image_dir,
                                   out_dir,
                                   slide_names        = NULL,
                                   fov_size_px        = 4256L,
                                   downsample_factor  = 8L,
                                   pixel_size_um      = 0.120280945,
                                   channels           = NULL,
                                   stitch_slides      = FALSE,
                                   slide_layout       = NULL,
                                   combine_downsample = 1L,
                                   slide_gap          = 0,
                                   combine_name       = "combined",
                                   he_images          = NULL,
                                   he_flip            = c("vertical", "horizontal", "both", "none"),
                                   he_resample_method = "average",
                                   generate_polygons  = TRUE,
                                   polygons           = NULL,
                                   state_col          = NULL,
                                   poly_slide_col     = "Run_Tissue_name",
                                   poly_cell_col      = "cell_ID",
                                   poly_x_col         = "x_slide_mm",
                                   poly_y_col         = "y_slide_mm",
                                   poly_coord_units   = c("mm", "px"),
                                   poly_flip_y        = TRUE,
                                   poly_x_offset_px   = 0,
                                   poly_y_offset_px   = NULL,
                                   poly_border_value  = 9999999L,
                                   poly_border_label  = "Cell Border",
                                   backend            = c("bftools", "rbioformats"),
                                   metadata_injector  = c("auto", "tiffcomment", "java"),
                                   bftools_dir        = NULL,
                                   max_vsize          = NULL,
                                   java_heap          = NULL,
                                   keep_intermediate  = FALSE,
                                   overwrite          = TRUE) {

  backend <- match.arg(backend)
  metadata_injector <- match.arg(metadata_injector)
  he_flip <- match.arg(he_flip)
  poly_coord_units <- match.arg(poly_coord_units)
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required. install.packages('terra').")
  downsample_factor  <- as.integer(downsample_factor)
  combine_downsample <- as.integer(combine_downsample)
  if (is.na(downsample_factor)  || downsample_factor  < 1L) stop("'downsample_factor' must be >= 1.")
  if (is.na(combine_downsample) || combine_downsample < 1L) stop("'combine_downsample' must be >= 1.")

  as_list <- function(x) if (is.data.frame(x) || (is.character(x) && length(x) == 1L)) list(x) else as.list(x)
  fov_positions <- as_list(fov_positions)
  image_dir     <- as_list(image_dir)
  n_slides <- length(fov_positions)
  if (length(image_dir) == 1L && n_slides > 1L) image_dir <- rep(image_dir, n_slides)
  if (length(image_dir) != n_slides) stop("'image_dir' must be length 1 or match the number of slides.")

  if (!is.null(he_images)) {
    he_images <- if (is.character(he_images) && length(he_images) == 1L) list(he_images) else as.list(he_images)
    if (length(he_images) != n_slides)
      stop("'he_images' (", length(he_images), ") must match the number of ",
           "IF slides (", n_slides, ").")
  }

  if (!is.null(slide_names)) {
    if (length(slide_names) != n_slides) stop("'slide_names' must match the number of slides.")
  } else if (!is.null(names(fov_positions)) && all(nzchar(names(fov_positions)))) {
    slide_names <- names(fov_positions)
  } else {
    # Default: name each slide after the folder containing its positions file
    # (the CosMx per-slide flatFiles folder, usually == Run_Tissue_name), so it
    # matches the slide column of a Seurat object / polygon table. Falls back to
    # slide1..N for non-path inputs or if the folder names are not usable.
    slide_names <- vapply(seq_len(n_slides), function(s) {
      p <- fov_positions[[s]]
      nm <- if (is.character(p)) basename(dirname(p)) else NA_character_
      if (is.na(nm) || !nzchar(nm) || nm %in% c(".", "/")) paste0("slide", s) else nm
    }, character(1))
    if (anyDuplicated(slide_names)) {
      message("Derived slide names are not unique (", paste(slide_names, collapse = ", "),
              "); falling back to slide1..", n_slides,
              ". Set 'slide_names' to control them.")
      slide_names <- paste0("slide", seq_len(n_slides))
    } else {
      message("Using slide name(s) from positions-file folder(s): ",
              paste(slide_names, collapse = ", "),
              " (override with 'slide_names').")
    }
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # ---- preflight (fail upfront, before any stitching) -----------------------
  # The OME-XML is injected either in-process via Bio-Formats (no external tool,
  # only for the rbioformats backend) or with the 'tiffcomment' CLI. bfconvert
  # builds the OME-TIFF container for the bftools backend and for the RGB H&E.
  java_inject <- FALSE
  if (backend == "rbioformats") {
    .gb_check_rbioformats(java_heap)            # inits JVM, loads RBioFormats
    java_inject <- .gb_java_injector_available()
  }
  resolved_injector <- metadata_injector
  if (resolved_injector == "auto")
    resolved_injector <- if (java_inject) "java" else "tiffcomment"
  if (resolved_injector == "java" && !java_inject)
    stop("metadata_injector = 'java' needs the rbioformats backend with a ",
         "working RBioFormats/Java (Bio-Formats TiffSaver classes were not ",
         "found). Use 'auto' or 'tiffcomment', or backend = 'rbioformats'.",
         call. = FALSE)

  need_bfconvert   <- (backend == "bftools") || !is.null(he_images)
  need_tiffcomment <- (resolved_injector == "tiffcomment")
  bfconvert <- tiffcomment <- ""
  if (need_tiffcomment) {
    tiffcomment <- .gb_find_bftool("tiffcomment", bftools_dir)
    if (!nzchar(tiffcomment))
      stop("Could not find 'tiffcomment' (from bftools), needed to embed the ",
           "OME metadata. Set 'bftools_dir' or add bftools to PATH, or (with ",
           "backend = 'rbioformats') use metadata_injector = 'java'/'auto'.",
           call. = FALSE)
  }
  if (need_bfconvert) {
    bfconvert <- .gb_find_bftool("bfconvert", bftools_dir)
    if (!nzchar(bfconvert)) {
      why <- if (backend == "bftools") "backend = 'bftools'" else "the RGB H&E output"
      stop("Could not find 'bfconvert' (from bftools), required for ", why,
           ". Set 'bftools_dir' or add bftools to your PATH.", call. = FALSE)
    }
  }

  # Inject an OME-XML comment into an existing OME-TIFF using the resolved
  # method, with a tiffcomment fallback when 'auto' picked Java but it fails.
  inject_ome_xml <- function(ome_tif, xml, stem) {
    if (resolved_injector == "java") {
      if (.gb_inject_comment_java(ome_tif, xml)) {
        message("  injected OME-XML (in-process, no tiffcomment)")
        return(invisible(TRUE))
      }
      if (metadata_injector == "java")
        stop("In-process Java injection failed for ", basename(ome_tif),
             " and metadata_injector = 'java' forbids fallback.", call. = FALSE)
      tc <- .gb_find_bftool("tiffcomment", bftools_dir)   # auto fallback
      if (!nzchar(tc))
        stop("Java injection failed and no 'tiffcomment' is available to fall ",
             "back on. Install bftools or repair RBioFormats/Java.", call. = FALSE)
      message("  falling back to tiffcomment")
    } else {
      tc <- tiffcomment
    }
    xml_file <- file.path(out_dir, paste0(stem, ".ome.xml"))
    writeLines(xml, xml_file, useBytes = TRUE)
    message("  tiffcomment -> inject OME-XML")
    st <- system2(tc, c("-set", shQuote(path.expand(xml_file)),
                        shQuote(path.expand(ome_tif))), stdout = FALSE, stderr = FALSE)
    if (!identical(st, 0L)) stop("tiffcomment failed for ", basename(ome_tif), ".")
    if (!keep_intermediate) unlink(xml_file)
    invisible(TRUE)
  }

  # ---- polygon preflight (validate/load before any stitching) ---------------
  do_polygons <- FALSE
  poly_data <- NULL
  poly_y_off <- if (is.null(poly_y_offset_px)) -fov_size_px else poly_y_offset_px
  if (isTRUE(generate_polygons)) {
    if (is.null(polygons)) {
      message("generate_polygons = TRUE but no 'polygons' supplied; skipping ",
              "polygon outputs. Provide 'polygons' (a Seurat object / vertex ",
              "table / RDS path, or a list of these) and 'state_col' to enable.")
    } else {
      message("Validating polygon inputs ...")
      poly_data <- .gb_load_polygons(polygons, slide_names, state_col,
                                     poly_slide_col, poly_cell_col,
                                     poly_x_col, poly_y_col)
      do_polygons <- TRUE
    }
  }

  # ---- helper: write an IF SpatRaster to a Minerva OME-TIFF -----------------
  # Both backends now: (1) write an intermediate multi-band <stem>.tif with
  # terra, (2) create the OME-TIFF container (bfconvert, or RBioFormats), then
  # (3) inject the SAME generated OME-XML (in-process via Bio-Formats, or with
  # tiffcomment), so channel names/colours/physical size are identical and
  # Minerva-readable either way.
  write_if_ome <- function(rast, stem, eff_um, channels_df) {
    ome_tif <- file.path(out_dir, paste0(stem, ".ome.tif"))
    if (file.exists(ome_tif) && !overwrite) stop("Exists (overwrite=FALSE): ", ome_tif)

    plain <- file.path(out_dir, paste0(stem, ".tif"))
    message("  writing intermediate ", basename(plain), " ...")
    terra::writeRaster(round(rast), plain, datatype = "INT2U", overwrite = TRUE,
      names = channels_df$name, NAflag = 0,
      wopt = list(gdal = c("COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=YES",
                           "BLOCKXSIZE=512", "BLOCKYSIZE=512",
                           "INTERLEAVE=BAND", "PHOTOMETRIC=MINISBLACK")))

    if (backend == "bftools") {
      message("  bfconvert -> ", basename(ome_tif))
      st <- system2(bfconvert, c("-overwrite", shQuote(path.expand(plain)),
                                 shQuote(path.expand(ome_tif))),
                    stdout = FALSE, stderr = FALSE)
      if (!identical(st, 0L)) stop("bfconvert failed for ", basename(ome_tif), ".")
    } else {                                   # rbioformats
      # RBioFormats has no streaming writer: the whole image must fit in RAM as
      # a double array (plus a transpose copy). Estimate the peak and lift
      # macOS's mem.maxVSize cap when possible; for large mosaics use bftools.
      peak_gb <- .gb_rbf_peak_gb(rast)
      message(sprintf(
        "  RBioFormats loads the full image into RAM: ~%.1f GB peak (%.0f Mpx x %d ch). ",
        peak_gb, terra::ncell(rast) / 1e6, terra::nlyr(rast)),
        "For large mosaics, backend = 'bftools' streams to disk and avoids this.")
      old_vsize <- .gb_raise_vsize(peak_gb * 1024, max_vsize)
      if (is.finite(old_vsize)) on.exit(try(mem.maxVSize(old_vsize), silent = TRUE), add = TRUE)

      arr <- terra::as.array(rast)             # [row, col, band]
      arr <- aperm(arr, c(2, 1, 3))            # -> [col, row, band]
      rbf_fail <- function(e) {
        msg <- conditionMessage(e)
        if (grepl("memory|vsize|heap|OutOfMemory", msg, ignore.case = TRUE))
          stop("RBioFormats ran out of memory (~", sprintf("%.1f", peak_gb),
               " GB needed). Use backend = 'bftools' (streams to disk), ",
               "increase 'max_vsize'/'java_heap', or downsample more. Original: ",
               msg, call. = FALSE)
        stop("RBioFormats::write.image failed (", msg,
             "); use backend = 'bftools'.", call. = FALSE)
      }
      # Prefer little-endian (to match the bftools output); fall back if the
      # installed API does not accept the argument.
      message("  RBioFormats -> ", basename(ome_tif))
      tryCatch(
        RBioFormats::write.image(arr, file = path.expand(ome_tif),
                                 force = overwrite, pixelType = "uint16",
                                 littleEndian = TRUE),
        error = function(e) {
          if (grepl("unused argument|littleEndian", conditionMessage(e)))
            tryCatch(RBioFormats::write.image(arr, file = path.expand(ome_tif),
                                              force = overwrite, pixelType = "uint16"),
                     error = rbf_fail)
          else rbf_fail(e)
        })
      rm(arr); gc(verbose = FALSE)
    }

    # Inject identical, structure-matching OME-XML (endianness detected from the
    # file just written, so the XML always agrees with the pixel byte order).
    big <- tryCatch(.gb_tiff_bigendian(ome_tif), error = function(e) FALSE)
    xml <- .gb_build_ome_xml(terra::ncol(rast), terra::nrow(rast), channels_df,
                             eff_um, image_name = stem, big_endian = big)
    inject_ome_xml(ome_tif, xml, stem)
    if (!keep_intermediate) unlink(plain)
    message("  wrote ", basename(ome_tif))
    ome_tif
  }

  # ---- helper: resample an H&E raster onto an IF grid and write RGB OME ------
  write_he_ome <- function(he_rast, grid_rast, stem, eff_um) {
    ome_tif <- file.path(out_dir, paste0(stem, ".ome.tif"))
    if (file.exists(ome_tif) && !overwrite) stop("Exists (overwrite=FALSE): ", ome_tif)
    message("  resampling H&E onto the IF grid (", he_resample_method, ")...")
    tmp <- tempfile(tmpdir = out_dir, fileext = ".tif")
    he_r <- terra::resample(he_rast, grid_rast, method = he_resample_method,
                            filename = tmp, overwrite = TRUE)
    plain <- file.path(out_dir, paste0(stem, ".tif"))
    terra::writeRaster(round(he_r), plain, datatype = "INT1U", overwrite = TRUE,
      names = c("Red", "Green", "Blue"), NAflag = 255,
      wopt = list(gdal = c("PHOTOMETRIC=RGB", "COMPRESS=LZW",
                           "INTERLEAVE=PIXEL", "BIGTIFF=YES")))
    unlink(tmp)
    message("  bfconvert -> ", basename(ome_tif))
    st <- system2(bfconvert, c("-overwrite", shQuote(path.expand(plain)),
                               shQuote(path.expand(ome_tif))),
                  stdout = FALSE, stderr = FALSE)
    if (!identical(st, 0L)) stop("bfconvert failed for ", basename(ome_tif), ".")
    big <- tryCatch(.gb_tiff_bigendian(ome_tif), error = function(e) FALSE)
    xml <- .gb_build_ome_xml_rgb(terra::ncol(grid_rast), terra::nrow(grid_rast),
                                 eff_um, image_name = stem, big_endian = big)
    inject_ome_xml(ome_tif, xml, stem)
    if (!keep_intermediate) unlink(plain)
    message("  wrote ", basename(ome_tif))
    ome_tif
  }

  # ---- helper: load + orient an H&E image into an IF slide's frame ----------
  load_he <- function(path, if_rast) {
    r <- terra::rast(path)
    if (he_flip != "none") r <- terra::flip(r, direction = he_flip)
    terra::ext(r) <- terra::ext(if_rast)
    terra::crs(r) <- terra::crs(if_rast)
    r
  }

  # ---- Stage 1: per-slide IF (and H&E) rasters ------------------------------
  if_list <- vector("list", n_slides)
  he_list <- if (!is.null(he_images)) vector("list", n_slides) else NULL
  for (s in seq_len(n_slides)) {
    message("\n=== Slide '", slide_names[s], "' (", s, "/", n_slides, ") ===")
    if_list[[s]] <- .gb_merge_slide(fov_positions[[s]], image_dir[[s]],
                                    pixel_size_um, fov_size_px,
                                    downsample_factor, out_dir)
    if (!is.null(he_images)) he_list[[s]] <- load_he(he_images[[s]], if_list[[s]])
  }
  channels_df <- .gb_resolve_channels(channels, terra::nlyr(if_list[[1]]))
  eff_um      <- pixel_size_um * downsample_factor
  written     <- character(0)

  do_combine <- isTRUE(stitch_slides) && n_slides > 1L
  if (isTRUE(stitch_slides) && n_slides == 1L && combine_downsample == 1L)
    message("Only one slide and combine_downsample = 1: nothing to combine.")

  # ---- Stage 2a: combined output --------------------------------------------
  if (do_combine || (isTRUE(stitch_slides) && combine_downsample > 1L)) {
    lay <- .gb_resolve_layout(n_slides, slide_layout)
    message("\nCombining ", n_slides, " slide(s) on a ", lay[1], " x ", lay[2],
            " grid (extra downsample x", combine_downsample, ")...")
    off <- .gb_grid_offsets(if_list, lay[1], lay[2], gap = slide_gap)

    shifted_if <- lapply(seq_len(n_slides), function(k)
      terra::shift(if_list[[k]], dx = off$dx[k], dy = off$dy[k]))
    tmp1 <- tempfile(tmpdir = out_dir, fileext = ".tif")
    merged_if <- terra::merge(terra::sprc(shifted_if), filename = tmp1,
                              datatype = "INT2U", NAflag = 0, overwrite = TRUE)
    if (combine_downsample > 1L) {
      tmp2 <- tempfile(tmpdir = out_dir, fileext = ".tif")
      merged_if <- terra::aggregate(merged_if, fact = combine_downsample,
                                    fun = "mean", filename = tmp2, NAflag = 0)
    }
    eff_um_c <- eff_um * combine_downsample
    written <- c(written, write_if_ome(merged_if, combine_name, eff_um_c, channels_df))

    if (!is.null(he_images)) {
      shifted_he <- lapply(seq_len(n_slides), function(k)
        terra::shift(he_list[[k]], dx = off$dx[k], dy = off$dy[k]))
      tmp3 <- tempfile(tmpdir = out_dir, fileext = ".tif")
      merged_he <- terra::merge(terra::sprc(shifted_he), filename = tmp3,
                                NAflag = 255, overwrite = TRUE)
      written <- c(written, write_he_ome(merged_he, merged_if,
                                         paste0(combine_name, "_he"), eff_um_c))
      unlink(tmp3)
    }

    # Polygons: combined (always) + per-slide (if keep_intermediate). Combined
    # polygons are shifted by the SAME grid offsets as the rasters.
    if (do_polygons) {
      shifts <- stats::setNames(
        lapply(seq_len(n_slides), function(k) c(off$dx[k], off$dy[k])), slide_names)
      written <- c(written, .gb_write_polygons(
        poly_data$verts, poly_data$state_map, merged_if, slide_names, shifts,
        out_dir, combine_name, pixel_size_um, fov_size_px, poly_coord_units,
        poly_flip_y, poly_x_offset_px, poly_y_off, poly_border_value,
        poly_border_label))
      if (isTRUE(keep_intermediate)) {
        for (s in slide_names) {
          k <- match(s, slide_names)
          written <- c(written, .gb_write_polygons(
            poly_data$verts, poly_data$state_map, if_list[[k]], s,
            stats::setNames(list(c(0, 0)), s), out_dir, s, pixel_size_um,
            fov_size_px, poly_coord_units, poly_flip_y, poly_x_offset_px,
            poly_y_off, poly_border_value, poly_border_label))
        }
      }
    }

    unlink(c(tmp1, if (exists("tmp2")) tmp2 else NULL))
    for (r in if_list) unlink(attr(r, "gb_tmp"))
    return(invisible(written))
  }

  # ---- Stage 2b: per-slide output -------------------------------------------
  for (s in seq_len(n_slides)) {
    stem <- slide_names[s]
    written <- c(written, write_if_ome(if_list[[s]], stem, eff_um, channels_df))
    if (!is.null(he_images))
      written <- c(written, write_he_ome(he_list[[s]], if_list[[s]],
                                         paste0(stem, "_he"), eff_um))
    if (do_polygons)
      written <- c(written, .gb_write_polygons(
        poly_data$verts, poly_data$state_map, if_list[[s]], stem,
        stats::setNames(list(c(0, 0)), stem), out_dir, stem, pixel_size_um,
        fov_size_px, poly_coord_units, poly_flip_y, poly_x_offset_px,
        poly_y_off, poly_border_value, poly_border_label))
    unlink(attr(if_list[[s]], "gb_tmp"))
  }
  invisible(written)
}


# ==============================================================================
# Command Line Execution Block  (run directly: Rscript stitch_fovs_to_ometiff.R ...)
# ==============================================================================
if (sys.nframe() == 0) {

  suppressWarnings(suppressMessages(
    if (!requireNamespace("optparse", quietly = TRUE))
      stop("The CLI needs 'optparse'. install.packages('optparse').")))

  option_list <- list(
    optparse::make_option(c("-p", "--fov_positions"), type = "character", default = NULL,
      help = "Comma-separated path(s) to *_fov_positions_file.csv(.gz), one per slide."),
    optparse::make_option(c("-d", "--image_dir"), type = "character", default = NULL,
      help = "Comma-separated Morphology2D dir(s): one per slide, or one for all."),
    optparse::make_option(c("-o", "--out_dir"), type = "character", default = ".",
      help = "Output directory [default= %default]"),
    optparse::make_option(c("-n", "--slide_names"), type = "character", default = NULL,
      help = "Comma-separated slide names (output stems)."),
    optparse::make_option(c("-s", "--fov_size_px"), type = "integer", default = 4256L,
      help = "FOV edge length in pixels [default= %default]"),
    optparse::make_option(c("-f", "--downsample_factor"), type = "integer", default = 8L,
      help = "Per-FOV downsampling factor [default= %default]"),
    optparse::make_option(c("-u", "--pixel_size_um"), type = "double", default = 0.120280945,
      help = "Original pixel size in microns [default= %default]"),
    optparse::make_option(c("-c", "--channels"), type = "character", default = NULL,
      help = "Comma-separated channel NAMES in band order, or preset 'cosmx_protein'."),
    optparse::make_option(c("-S", "--stitch_slides"), action = "store_true", default = FALSE,
      help = "Combine all slides into one image [default= %default]"),
    optparse::make_option(c("-L", "--slide_layout"), type = "character", default = NULL,
      help = "Grid as 'nrow,ncol' (or a single number = nrow). Default: single row."),
    optparse::make_option(c("-D", "--combine_downsample"), type = "integer", default = 1L,
      help = "Extra downsampling of the combined image [default= %default]"),
    optparse::make_option(c("-g", "--slide_gap"), type = "double", default = 0,
      help = "Gap between slides in px [default= %default]"),
    optparse::make_option(c("-N", "--combine_name"), type = "character", default = "combined",
      help = "Output stem for the combined image(s) [default= %default]"),
    optparse::make_option(c("-H", "--he_images"), type = "character", default = NULL,
      help = "Comma-separated co-registered H&E tif(s), one per slide."),
    optparse::make_option(c("--he_flip"), type = "character", default = "vertical",
      help = "Flip H&E: vertical/horizontal/both/none [default= %default]"),
    optparse::make_option(c("--he_resample_method"), type = "character", default = "average",
      help = "terra resampling method for H&E [default= %default]"),
    optparse::make_option(c("-b", "--backend"), type = "character", default = "bftools",
      help = "'bftools' or 'rbioformats' [default= %default]"),
    optparse::make_option(c("--metadata_injector"), type = "character", default = "auto",
      help = "'auto', 'java' (rbioformats, no tiffcomment), or 'tiffcomment' [default= %default]"),
    optparse::make_option(c("-t", "--bftools_dir"), type = "character", default = NULL,
      help = "Directory with bfconvert/tiffcomment (else PATH)."),
    optparse::make_option(c("--max_vsize"), type = "double", default = NULL,
      help = "rbioformats/macOS only: R vector-memory limit in MB (mem.maxVSize)."),
    optparse::make_option(c("--java_heap"), type = "character", default = NULL,
      help = "rbioformats only: JVM max heap, e.g. '24g' or '18000m' (default ~70%% RAM)."),
    optparse::make_option(c("--no_polygons"), action = "store_true", default = FALSE,
      help = "Disable Minerva polygon output."),
    optparse::make_option(c("--polygons"), type = "character", default = NULL,
      help = "Comma-separated .rds path(s) to polygon vertex table(s)/Seurat object(s)."),
    optparse::make_option(c("--state_col"), type = "character", default = NULL,
      help = "Column to use as the Minerva State (required for polygons)."),
    optparse::make_option(c("--poly_slide_col"), type = "character", default = "Run_Tissue_name",
      help = "Slide-id column in the polygon input [default= %default]"),
    optparse::make_option(c("--poly_cell_col"), type = "character", default = "cell_ID",
      help = "Cell-id column in the polygon table [default= %default]"),
    optparse::make_option(c("--poly_x_col"), type = "character", default = "x_slide_mm",
      help = "x column in the polygon table [default= %default]"),
    optparse::make_option(c("--poly_y_col"), type = "character", default = "y_slide_mm",
      help = "y column in the polygon table [default= %default]"),
    optparse::make_option(c("--poly_coord_units"), type = "character", default = "mm",
      help = "'mm' or 'px' for polygon coordinates [default= %default]"),
    optparse::make_option(c("-k", "--keep_intermediate"), action = "store_true", default = FALSE,
      help = "Keep intermediate .tif/.xml (and per-slide polygons) [default= %default]"),
    optparse::make_option(c("--no_overwrite"), action = "store_true", default = FALSE,
      help = "Do not overwrite existing outputs.")
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  if (is.null(opt$fov_positions) || is.null(opt$image_dir)) {
    optparse::print_help(optparse::OptionParser(option_list = option_list))
    stop("Both --fov_positions and --image_dir are required.", call. = FALSE)
  }
  split_csv <- function(x) if (is.null(x)) NULL else trimws(unlist(strsplit(x, ",")))

  ch <- split_csv(opt$channels)
  if (!is.null(ch) && length(ch) == 1L && ch %in% c("cosmx_protein")) ch <- ch[[1]]
  layout <- if (is.null(opt$slide_layout)) NULL else as.integer(split_csv(opt$slide_layout))

  stitch_fovs_to_ometiff(
    fov_positions      = as.list(split_csv(opt$fov_positions)),
    image_dir          = as.list(split_csv(opt$image_dir)),
    out_dir            = opt$out_dir,
    slide_names        = split_csv(opt$slide_names),
    fov_size_px        = opt$fov_size_px,
    downsample_factor  = opt$downsample_factor,
    pixel_size_um      = opt$pixel_size_um,
    channels           = ch,
    stitch_slides      = opt$stitch_slides,
    slide_layout       = layout,
    combine_downsample = opt$combine_downsample,
    slide_gap          = opt$slide_gap,
    combine_name       = opt$combine_name,
    he_images          = split_csv(opt$he_images),
    he_flip            = opt$he_flip,
    he_resample_method = opt$he_resample_method,
    generate_polygons  = !opt$no_polygons,
    polygons           = split_csv(opt$polygons),
    state_col          = opt$state_col,
    poly_slide_col     = opt$poly_slide_col,
    poly_cell_col      = opt$poly_cell_col,
    poly_x_col         = opt$poly_x_col,
    poly_y_col         = opt$poly_y_col,
    poly_coord_units   = opt$poly_coord_units,
    backend            = opt$backend,
    metadata_injector  = opt$metadata_injector,
    bftools_dir        = opt$bftools_dir,
    max_vsize          = opt$max_vsize,
    java_heap          = opt$java_heap,
    keep_intermediate  = opt$keep_intermediate,
    overwrite          = !opt$no_overwrite
  )
}
