#' Arrange tissues in xy space to reduce whitespace
#' 
#' Uses a shelf algorithm: places tallest tissues on the bottom shelf, and so on. 
#'
#' @param xy A matrix of xy coordinates.
#' @param tissue A vector of tissue labels corresponding to the rows of xy.
#' @param tissueorder Optional vector specifying the order tissues are tiled in.
#' @param buffer Space between tissues.
#' @param widthheightratio Desired shape of final xy locations.
#' @return A matrix of updated xy coordinates.
#' @export
condenseTissues <- function(xy, tissue, tissueorder = NULL, buffer = 0.2, widthheightratio = 4/3) {
  
  # get each tissue's dimensions:
  tissdf <- data.frame(tissue = unique(tissue))
  tissdf$width <- sapply(unique(tissue), function(tiss) {
    diff(range(xy[tissue == tiss, 1]))
  })
  tissdf$height <- sapply(unique(tissue), function(tiss) {
    diff(range(xy[tissue == tiss, 2]))
  })
  
  # choose tissue order:
  if (!is.null(tissueorder)) {
    if (length(setdiff(tissdf$tissue, tissueorder)) > 0) {
      stop("values in tissue missing from tissueorder")
    }
    if (length(setdiff(tissueorder, tissdf$tissue)) > 0) {
      stop("values in tissueorder missing from tissue")
    }
    tissdf$order <- match(tissdf$tissue, tissueorder)
  } else {
    tissdf$order <- order(tissdf$height, decreasing = TRUE)
  }
  tissdf <- tissdf[tissdf$order, ]
  
  # choose number of tissues for first shelf:
  tissuesperrow <- round(sqrt(nrow(tissdf)) * widthheightratio * mean(tissdf$height) / mean(tissdf$width))
  targetwidth <- sum(tissdf$width[1:tissuesperrow], na.rm = TRUE) + buffer * (tissuesperrow - 1)
  
  # place tissues:
  tissdf$x <- NA
  tissdf$y <- NA
  tempx <- 0
  tempy <- 0
  tempshelfheight <- 0
  tempshelfwidth <- 0
  for (i in 1:nrow(tissdf)) {
    # place this tissue:
    tissdf$x[i] <- tempx
    tissdf$y[i] <- tempy
    # update the shelf dimensions:
    tempshelfheight <- max(tempshelfheight, tissdf$height[i])
    tempshelfwidth <- tempx + tissdf$width[i]
    # move along the shelf:
    tempx <- tempx + tissdf$width[i] + buffer
    # start a new shelf if it's getting too wide:
    if (i < nrow(tissdf)) {
      if (abs(tempshelfwidth - targetwidth) < abs(tempshelfwidth + buffer + tissdf$width[i+1] - targetwidth)) {
        tempy <- tempy + tempshelfheight + buffer
        tempx <- 0
        tempshelfheight <- 0
        tempshelfwidth <- 0
      }
    }
  }
  # now update xy:
  for (tiss in unique(tissue)) {
    inds <- tissue == tiss
    xy[inds, 1] <- xy[inds, 1] - min(xy[inds, 1]) + tissdf$x[tissdf$tissue == tiss]
    xy[inds, 2] <- xy[inds, 2] - min(xy[inds, 2]) + tissdf$y[tissdf$tissue == tiss]
  }
  return(xy)  
}

#' Data Preparation for CosMx
#'
#' Processes raw CosMx flat files into a structured object containing count matrices,
#' metadata, and condensed spatial coordinates.
#'
#' @param myflatfiledir A character vector of one or more directories containing CosMx slide folders.
#' @param plot_tissues Logical. If TRUE, plots the condensed tissue layouts. Default is FALSE.
#' @return A list containing `"counts"`, `"negcounts"`, `"falsecounts"`, `"metadata"`, and `"xy"`.
#' @importFrom data.table fread rbindlist :=
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom graphics plot text
#' @importFrom stats median
#' @export
dataprep_cosmx <- function(myflatfiledir, plot_tissues = FALSE) {
  
  # Automatically get slide paths from the provided directories
  slide_paths <- character()
  for (dir_path in myflatfiledir) {
    s_names <- dir(dir_path)
    # Ensure we only keep directories (ignoring loose files in the parent dir)
    valid_dirs <- s_names[dir.exists(file.path(dir_path, s_names))]
    slide_paths <- c(slide_paths, file.path(dir_path, valid_dirs))
  }
  
  if (length(slide_paths) == 0) {
    stop("No slide directories found in the provided path(s).")
  }
  
  slidenames <- basename(slide_paths)
  
  # Lists to collect the counts matrices and metadata, one per slide
  countlist <- vector(mode = 'list', length = length(slide_paths)) 
  metadatalist <- vector(mode = 'list', length = length(slide_paths)) 
  
  for(i in seq_along(slide_paths)) {
    
    current_path <- slide_paths[i]
    slidename <- slidenames[i] 
    
    msg <- paste0("Loading slide ", slidename, ", ", i, "/", length(slide_paths), ".")
    message(msg)    
    
    # slide-specific files:
    thisslidesfiles <- dir(current_path)
    
    # load in metadata:
    thisslidesmetadata <- thisslidesfiles[grepl("metadata\\_file", thisslidesfiles)]
    if (length(thisslidesmetadata) == 0) stop(paste("No metadata file found for", slidename))
    
    tempdatatable <- data.table::fread(file.path(current_path, thisslidesmetadata))
    tempdatatable[, slidename := slidename]
    
    # numeric slide ID 
    slide_ID_numeric <- tempdatatable[1,]$slide_ID 
    
    # global cell ID (fixed logic to preserve cell ID mapping before removing it later)
    tempdatatable$global_cell_ID <- paste0("c_", slide_ID_numeric, "_", tempdatatable$fov, "_", tempdatatable$cell_ID)
    
    # load in counts as a data table:
    thisslidescounts <- thisslidesfiles[grepl("exprMat\\_file", thisslidesfiles)]
    if (length(thisslidescounts) == 0) stop(paste("No exprMat file found for", slidename))
    
    countsfile <- file.path(current_path, thisslidescounts)
    nonzero_elements_perchunk <- 5 * 10^7
    
    ### Safely read in the dense (0-filled) counts matrices in chunks.
    lastchunk <- FALSE 
    skiprows <- 0
    chunkid <- 1
    
    required_cols <- data.table::fread(countsfile, select=c("fov", "cell_ID"))
    stopifnot("columns 'fov' and 'cell_ID' are required, but not found in the counts file" = 
                all(c("cell_ID", "fov") %in% colnames(required_cols)))
    number_of_cells <- nrow(required_cols)
    
    number_of_cols <- ncol(data.table::fread(countsfile, nrows = 2))
    number_of_chunks <- ceiling(number_of_cols * number_of_cells / nonzero_elements_perchunk)
    chunk_size <- floor(number_of_cells / number_of_chunks)
    sub_counts_matrix <- vector(mode = 'list', length = number_of_chunks)
    
    pb <- txtProgressBar(min = 0, max = number_of_chunks, initial = 0, char = "=", style = 3)
    cellcount <- 0
    
    while(lastchunk == FALSE) {
      read_header <- ifelse(chunkid == 1, TRUE, FALSE)
      
      countsdatatable <- data.table::fread(countsfile,
                                           nrows = chunk_size,
                                           skip = skiprows + (chunkid > 1),
                                           header = read_header)
      if(chunkid == 1) {
        header <- colnames(countsdatatable)
      } else {
        colnames(countsdatatable) <- header
      }
      
      cellcount <- nrow(countsdatatable) + cellcount     
      if(cellcount == number_of_cells) lastchunk <- TRUE
      
      skiprows <- skiprows + chunk_size
      slide_fov_cell_counts <- paste0("c_", slide_ID_numeric, "_", countsdatatable$fov, "_", countsdatatable$cell_ID)
      
      # Convert to sparseMatrix
      sub_counts_matrix[[chunkid]] <- as(countsdatatable[, -c("fov", "cell_ID"), with = FALSE], "sparseMatrix") 
      rownames(sub_counts_matrix[[chunkid]]) <- slide_fov_cell_counts 
      
      setTxtProgressBar(pb, chunkid)
      chunkid <- chunkid + 1
    }
    
    close(pb)   
    
    countlist[[i]] <- do.call(rbind, sub_counts_matrix) 
    
    # ensure that cell-order in counts matches cell-order in metadata   
    countlist[[i]] <- countlist[[i]][match(tempdatatable$global_cell_ID, rownames(countlist[[i]])), ] 
    metadatalist[[i]] <- tempdatatable 
    
    # track common genes and common metadata columns across slides
    if(i == 1) {
      sharedgenes <- colnames(countlist[[i]]) 
      sharedcolumns <- colnames(tempdatatable)
    } else {
      sharedgenes <- intersect(sharedgenes, colnames(countlist[[i]]))
      sharedcolumns <- intersect(sharedcolumns, colnames(tempdatatable))
    }
  }
  
  # reduce to shared metadata columns and shared genes
  for(i in seq_along(slide_paths)) {
    metadatalist[[i]] <- metadatalist[[i]][, sharedcolumns, with = FALSE]
    countlist[[i]] <- countlist[[i]][, sharedgenes, drop = FALSE]
  }
  
  counts <- do.call(rbind, countlist)
  metadata <- data.table::rbindlist(metadatalist)
  
  # add to metadata: add a global non-slide-specific FOV ID:
  metadata$FOV <- paste0("s", metadata$slide_ID, "f", metadata$fov)
  
  # remove cell_ID metadata column, which only identifies cell within slides, not across slides:
  metadata$cell_ID <- NULL
  
  # add coordinates in mm
  um_per_px <- 0.120280945 # dimension of each pixel in micro-meter
  metadata$x_slide_mm <- um_per_px * metadata$CenterX_global_px / 1e3 
  metadata$y_slide_mm <- um_per_px * metadata$CenterY_global_px / 1e3
  
  # isolate negative control matrices:
  negcounts <- counts[, grepl("Negative", colnames(counts)), drop = FALSE]
  falsecounts <- counts[, grepl("SystemControl", colnames(counts)), drop = FALSE]
  
  # reduce counts matrix to only genes:
  counts <- counts[, !grepl("Negative", colnames(counts)) & !grepl("SystemControl", colnames(counts)), drop = FALSE]
  
  # Then break out cells' xy positions in a distinct data object:
  xy <- as.matrix(metadata[, c("CenterX_global_px", "CenterY_global_px"), with = FALSE])
  
  # Use the generated global_cell_ID for rownames to ensure accuracy 
  rownames(xy) <- metadata$global_cell_ID
  
  # rescale to mm:
  thisinstrument_nanometers_per_pixel = 120.280945   
  xy <- xy * thisinstrument_nanometers_per_pixel / 1000000
  colnames(xy) <- paste0(c("x", "y"), "_mm")
  
  # Condense tissues
  set.seed(1)
  xy <- condenseTissues(xy = xy, 
                        tissue = metadata$Run_Tissue_name, 
                        tissueorder = NULL,  
                        buffer = 1, 
                        widthheightratio = 4/3) 
  
  # Optional: Plot tissues if requested
  if (plot_tissues) {
    sub <- sample(1:nrow(xy), round(nrow(xy) / 20), replace = FALSE)
    plot(xy[sub, ], pch = 16, cex = 0.2, 
         asp = 1, 
         col = as.numeric(as.factor(metadata$Run_Tissue_name[sub])),
         main = "Condensed Tissues Layout")
    
    for (s_name in unique(metadata$Run_Tissue_name)) {
      text(median(xy[metadata$Run_Tissue_name == s_name, 1]), 
           max(xy[metadata$Run_Tissue_name == s_name, 2]), 
           s_name)
    }
  }
  
  # Return final list object
  result_obj <- list(
    counts = counts,
    negcounts = negcounts,
    falsecounts = falsecounts,
    metadata = metadata,
    xy = xy
  )
  
  return(result_obj)
}
