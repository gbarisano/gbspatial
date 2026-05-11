# ==============================================================================
# SPATIAL TRANSCRIPTOMICS QC PIPELINE
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. FOV QC UTILITY FUNCTIONS # Copyright ©2024 Bruker Spatial Biology, Inc. All rights reserved. Subject to additional license terms and conditions provided separately by Bruker Spatial Biology, Inc.
# ------------------------------------------------------------------------------

#' FOV QC 
runFOVQC <- function(counts, xy, fov, tissue = NULL, barcodemap, max_prop_loss = 0.6, max_totalcounts_loss = 0.6) {
  
  if ((max_prop_loss > 1) | (max_prop_loss < 0)) {
    stop("max_prop_loss must fall in range of 0-1.")
  }
  fov <- paste0(tissue, as.character(fov))
  ## create a matrix of barcode bit expression over sub-FOV grids:
  # define grids, get per-square gene expression:
  gridinfo <- makeGrid(xy = xy, fov = fov, squares_per_fov = 49, min_cells_per_square = 10) 
  # convert to per-square barcode bit expression:
  bitcounts <- cellxgene2squarexbit(counts = counts, 
                                    grid = gridinfo$gridid, 
                                    genes = barcodemap$gene, 
                                    barcodes = barcodemap$barcode) 
  
  # normalize:
  bitcounts <- sweep(bitcounts, 1, pmax(rowSums(bitcounts), 1), "/") * mean(pmax(rowSums(bitcounts), 1))

  ## for every grid square, match it to "control" squares from other FOVs, and get its residuals from them:
  # get neighbors:
  comparators <- getNearestNeighborsByFOV(x = bitcounts, 
                                          gridfov = gridinfo$gridfov, 
                                          n_neighbors = 10)
  
  ## Search for FOVs with low outlier total counts:
  totalcounts <- Matrix::rowSums(counts)
  totalcountsgrid <- by(totalcounts, gridinfo$gridid, mean)[rownames(bitcounts)]
  # get expected total counts:
  totalcountshat <- sapply(1:length(totalcountsgrid), function(i) {
    mean(totalcountsgrid[comparators[i, ]], na.rm = TRUE)
  })
  # resids of total counts:
  totalcountsresids <- log2(totalcountsgrid / totalcountshat)
  # flag FOVs:
  flaggedgrids <- totalcountsresids < log2(1 - max_totalcounts_loss)
  flaggedgridsperfov <- by(as.vector(1*flaggedgrids), gridinfo$gridfov[names(flaggedgrids)], mean)
  flaggedfovs_fortotalcounts <- names(which(flaggedgridsperfov > 0.75))
  
  ## Search for failed reporter cycles:
  # get expected:
  yhat <- Matrix::t(sapply(1:nrow(bitcounts), function(i) {
    colMeans(bitcounts[comparators[i, ], ], na.rm = TRUE)
  }))
  # get resids:
  resid <- log2((bitcounts + 1) / (yhat + 1))
  rownames(resid) = rownames(bitcounts)
  
  ## summarize bias per FOV * bit:
  fovstats <- summarizeFOVBias(resid = resid, gridfov = gridinfo$gridfov, max_prop_loss = max_prop_loss)
  
  # count flags per FOV * reportercycle:
  flags_per_fov_x_reportercycle <- c()
  reportercycle <- substr(colnames(fovstats$flag), 1, nchar(colnames(fovstats$flag)) - 1)
  for (rc in unique(reportercycle)) {
    flags_per_fov_x_reportercycle <- cbind(flags_per_fov_x_reportercycle, rowMeans(fovstats$flag[, reportercycle == rc]))
  }
  colnames(flags_per_fov_x_reportercycle) <- unique(reportercycle)
  
  # collate all flagged FOVs:
  flaggedfovs_forbias <- rownames(flags_per_fov_x_reportercycle)[rowSums(flags_per_fov_x_reportercycle >= 0.5) > 0]
  flaggedfovs <- union(flaggedfovs_fortotalcounts, flaggedfovs_forbias)
  
  # report on flagged FOVs:
  if (length(flaggedfovs) > 0) {
    message(paste0("The following FOVs failed QC for one or more barcode positions: ",
                   paste0(flaggedfovs, collapse = ", ")))
  } else {
    message("All FOVs passed QC")
  }
  
  # build a manifest of flagged gene/fov pairs:
  flagged_fov_x_gene <- c()
  for (f in flaggedfovs) {
    flaggedreportercycles <- names(which(flags_per_fov_x_reportercycle[f, ] >= 0.5))
    for (rc in flaggedreportercycles) {
      reporterposition <- as.numeric(substr(rc, 14, nchar(rc) - 1)) * 2
      flaggedgenes <- barcodemap$gene[substr(barcodemap$barcode, reporterposition, reporterposition) != "."]
      # add to growing list:
      tempflags <- cbind(rep(f, length(flaggedgenes)), flaggedgenes)
      colnames(tempflags) <- c("fov", "gene")
      flagged_fov_x_gene <- rbind(flagged_fov_x_gene, tempflags)
    }
  }
  
  return(list(flaggedfovs = flaggedfovs, 
              flaggedfovs_fortotalcounts = flaggedfovs_fortotalcounts, 
              flaggedfovs_forbias = flaggedfovs_forbias,
              flagged_fov_x_gene = flagged_fov_x_gene, 
              flags_per_fov_x_reportercycle = flags_per_fov_x_reportercycle ,
              fovstats = fovstats, resid = resid, totalcountsresids = totalcountsresids, 
              gridinfo = gridinfo, xy = xy, fov = fov))
}

#' Infer colorvals from barcodes:
getcolorvals <- function(barcodes) {
  allvals <- paste0(barcodes, collapse = "")
  uniquevals <- unique(unlist(strsplit(allvals, "")))
  return(setdiff(uniquevals, "."))
}

#' Spatial plots of FOV effects:
FOVEffectsSpatialPlots <- function(res, outdir = NULL, bits = "flagged_reportercycles", plotwidth = NULL, plotheight = NULL) {
  bitnames <- colnames(res$fovstats$p)
  colorvals <- unique(substr(bitnames, nchar(bitnames), nchar(bitnames)))
  
  if (is.null(plotwidth)) {
    plotwidth <- diff(range(res$xy[, 1])) * 1.5
  }
  if (is.null(plotheight)) {
    plotheight <- diff(range(res$xy[, 2])) * 1.5
  }
  bits_to_plot <- match(bits, colnames(res$fovstats$flag))
  if (bits == "flagged_reportercycles") {
    flaggedreportercycles <- colnames(res$flags_per_fov_x_reportercycle)[colSums(res$flags_per_fov_x_reportercycle >= 0.5) > 0]
    names_of_bits_to_plot  <- paste0(rep(flaggedreportercycles, each = 4), rep(colorvals, length(flaggedreportercycles)))
    bits_to_plot <- match(names_of_bits_to_plot, colnames(res$resid))
  }
  if (bits == "flagged_bits") {
    bits_to_plot  <- which(colSums(res$fovstats$flag) > 0)
  }
  if (bits == "all") {
    bits_to_plot <- 1:ncol(res$resid)
  }
  temp <- sapply(bits_to_plot, function(i) {
    if (!is.null(outdir)) {
      png(paste0(outdir, "/", make.names(colnames(res$resid)[i]), ".png"), width = plotwidth, height = plotheight, units = "in", res = 300)
    }
    par(mar = c(0,0,2,0))
    plot(res$xy, cex = 0.2, asp = 1, pch = 16,
         col = colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)[
           pmax(pmin(51 + res$resid[match(res$gridinfo$gridid, rownames(res$resid)), i] * 50, 101), 1)], 
         main = paste0(colnames(res$resid)[i], ": log2(fold-change)\nfrom comparable regions elsewhere"))
    for (f in unique(res$fov)) {
      inds <- res$fov == f
      rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "black")
    }
    for (f in rownames(res$fovstats$flag)[res$fovstats$flag[, i] > 0]) {
      inds <- res$fov == f
      rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), lwd = 2, border = "yellow")
    }
    legend("right", pch = 16,
           col = rev(c("darkblue", "blue", "grey80", "red", "darkred")),
           legend = rev(c("< -1", -0.5, 0, 0.5, "> 1")))
    if (!is.null(outdir)) {
      dev.off()
    }
  })
}

#' Spatial plot of loss in signal strength compared to similar regions:
FOVSignalLossSpatialPlot <- function(res, shownames = TRUE, outdir = NULL, plotwidth = NULL, plotheight = NULL) {
  if (!is.null(outdir)) {
    png(paste0(outdir, "/signal loss.png"), width = plotwidth, height = plotheight, units = "in", res = 300)
  }
  if (is.null(plotwidth)) {
    plotwidth <- diff(range(res$xy[, 1])) * 1.5
  }
  if (is.null(plotheight)) {
    plotheight <- diff(range(res$xy[, 2])) * 1.5
  }
  
  plot(res$xy, cex = 0.2, asp = 1, pch = 16,
       col = colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)[
         pmax(pmin(51 + res$totalcountsresids[match(res$gridinfo$gridid, names(res$totalcountsresids))] * 25, 101), 1)], 
       main = "Log2 fold-change in total counts compared to similar regions")
  for (f in unique(res$fov)) {
    inds <- res$fov == f
    rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "black")
  }
  for (f in res$flaggedfovs_fortotalcounts) {
    inds <- res$fov == f
    rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "yellow", lwd = 2)
    if (shownames) {
      text(median(range(res$xy[inds, 1])), median(range(res$xy[inds, 2])), f, col = "green")
    }
  }
  legend("right", pch = 16,
         col = rev(c("darkblue", "blue", "grey80", "red", "darkred")),
         legend = rev(c("< -2", -1, 0, 1, "> 2")))
  if (!is.null(outdir)) {
    dev.off()
  }
}

#' Map of which FOVs were flagged:
mapFlaggedFOVs <- function(res, shownames = TRUE, outdir = NULL, plotwidth = NULL, plotheight = NULL) {
  if (!is.null(outdir)) {
    png(paste0(outdir, "/flagged FOVs.png"), width = plotwidth, height = plotheight, units = "in", res = 300)
  }
  if (is.null(plotwidth)) {
    plotwidth <- diff(range(res$xy[, 1])) * 1.5
  }
  if (is.null(plotheight)) {
    plotheight <- diff(range(res$xy[, 2])) * 1.5
  }
  
  plot(res$xy, cex = 0.1, asp = 1, pch = 16,
       col = "grey80", 
       main = "Flagged FOVs")
  for (f in unique(res$fov)) {
    inds <- res$fov == f
    rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), col = scales::alpha("dodgerblue2", 0.5))
  }
  for (f in res$flaggedfovs) {
    inds <- res$fov == f
    rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), col = scales::alpha("red", 0.5))
    if (shownames) {
      text(median(range(res$xy[inds, 1])), median(range(res$xy[inds, 2])), f, col = "green")
    }
  }
  if (!is.null(outdir)) {
    dev.off()
  }
}

#' Heatmap of estimated bit bias across FOVs
FOVEffectsHeatmap <- function(res) {
  pheatmap::pheatmap(res$fovstats$bias * (res$fovstats$flag),
                     col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100),
                     breaks = seq(-2,2,length.out = 101),
                     main = "FOV bias: log2(fold-change) from comparable regions in other FOVs")
}

#' Get nearest neighbors, sampling diffusely across other FOVs. 
getNearestNeighborsByFOV <- function(x, gridfov, n_neighbors = 10) {
  gridfov <- gridfov[rownames(x)]
  topneighbors <- FNN::get.knnx(data = x, query = x, k = min(n_neighbors * 50, nrow(x)))$nn.index
  reducedneighbors <- t(sapply(1:nrow(topneighbors), function(i) {
    tempneighbors <- topneighbors[i, ]
    tempneighbors[duplicated(gridfov[topneighbors[i, ]])] <- NA
    tempneighbors[gridfov[topneighbors[i, ]] == gridfov[i]] <- NA
    return(tempneighbors[!is.na(tempneighbors)][1:n_neighbors])
  }))
  return(reducedneighbors)
}

#' Define a grid across FOVs 
makeGrid <- function(xy, fov, squares_per_fov = 49, min_cells_per_square = 25) {
  ncuts <- floor(sqrt(squares_per_fov))
  grid <- rep(NA, nrow(xy))
  gridfov <- c()
  for (fovid in unique(fov)) {   
    inds <- fov == fovid
    grid[inds] <- paste0(fov[inds], "_", cut(xy[inds, 1], ncuts), "_", cut(xy[inds, 2], ncuts))
    gridfov[unique(grid[inds])] <- rep(fovid, length(unique(grid[inds])))
  }
  toosmallgridids <- names(which(table(grid) < min_cells_per_square))
  grid[is.element(grid, toosmallgridids)] <- NA
  return(list(gridid = grid, gridfov = gridfov))
}

#' Convert cell x gene matrix to grid square * bit matrix
cellxgene2squarexbit <- function(counts, grid, genes, barcodes) {
  colorvals <- getcolorvals(barcodes)
  nreportercycles <- nchar(barcodes[1]) / 2
  nbits <- nreportercycles * 4
  bitmat = matrix(0, length(setdiff(unique(grid), NA)), nbits)
  colnames(bitmat) <- paste0("reportercycle", rep(seq_len(nreportercycles), each = 4), colorvals)
  rownames(bitmat) <- setdiff(unique(grid), NA)
  
  inasquare <- !is.na(grid)
  counts <- counts[inasquare, ]
  grid <- grid[inasquare]
  
  gridmap <- Matrix::sparseMatrix(
    i = as.numeric(as.factor(grid)),
    j = seq_len(length(grid)),
    x = 1,
    dims = c(length(levels(as.factor(grid))), length(grid)))
  rownames(gridmap) <- levels(as.factor(grid))
  colnames(gridmap) <- rownames(counts)
  
  gridxgenecounts <- gridmap[, rownames(counts)] %*% counts 
  ncellspergrid <- Matrix::rowSums(gridmap)
  gridxgenecounts <- Matrix::Diagonal(x = 1/ncellspergrid) %*% gridxgenecounts
  rownames(gridxgenecounts) <- rownames(gridmap)
  
  gene2bitmap <- barcode2bitmatrix(barcodes)
  rownames(gene2bitmap) <- genes
  sharedgenes <- intersect(genes, colnames(counts))
  gridxbitcounts <- gridxgenecounts[, sharedgenes]%*% gene2bitmap[sharedgenes, ]
  return(as.matrix(gridxbitcounts))
}

#' Convert the barcode vector to a matrix of bit assignments (genes * bits)
barcode2bitmatrix <- function(barcodes) {
  colorvals <- getcolorvals(barcodes)
  nreportercycles <- nchar(barcodes[1]) / 2
  nbits <- nreportercycles * 4
  bitmap <- matrix(0, length(barcodes), nbits)
  colnames(bitmap) <- paste0("reportercycle", rep(seq_len(nreportercycles), each = 4), rep(colorvals, nreportercycles))
  for (i in seq_len(nreportercycles)) {
    barcodeposition <- i*2
    barcodehere <- substr(barcodes, barcodeposition, barcodeposition)
    for (col in colorvals) {
      bitmap[barcodehere == col, paste0("reportercycle", i, col)] <- 1
    }
  }
  return(bitmap)
}

#' Summarize bias in FOVs
summarizeFOVBias <- function(resid, gridfov, max_prop_loss) {
  gridfov = gridfov[rownames(resid)]
  fovs = unique(gridfov)
  bias <- p <- propagree <- matrix(NA, length(fovs), ncol(resid), dimnames = list(fovs, colnames(resid)))
  
  for (bit in colnames(resid)) {
    mod = summary(lm(resid[, bit] ~ as.factor(gridfov) - 1))$coef
    bias[gsub("as.factor\\(gridfov\\)", "", rownames(mod)), bit] = mod[, "Estimate"]
    p[gsub("as.factor\\(gridfov\\)", "", rownames(mod)), bit] = mod[, "Pr(>|t|)"]
    for (fov in fovs) {
      inds <- gridfov == fov
      propagree[as.character(fov), bit] <- mean(resid[inds, bit] < log2(1 - max_prop_loss) / 3)
    }
  }
  flag <- (bias < log2(1 - max_prop_loss)) * (p < 0.01) * (propagree >= 0.5)
  return(list(flag = flag, bias = bias, p = p, propagree = propagree))
}

# ------------------------------------------------------------------------------
# 2. MAIN PIPELINE WRAPPER
# ------------------------------------------------------------------------------
#' Comprehensive Spatial Transcriptomics QC Pipeline
#'
#' @param counts Matrix or data frame of gene expression counts.
#' @param metadata Data frame containing cell metadata (must contain 'cell_id').
#' @param xy Data frame or matrix of spatial coordinates (required if do_runfovqc or do_regional is TRUE).
#' @param negcounts Matrix or data frame of negative probe counts (required if do_regional is TRUE).
#' @param do_nCount Logical, perform nCount_RNA QC.
#' @param count_threshold Numeric, minimum count threshold.
#' @param do_cellArea Logical, perform cell area QC.
#' @param area_threshold Numeric, maximum allowed cell area.
#' @param do_runfovqc Logical, perform FOV-level QC.
#' @param panel_name Character, panel name for barcodes (e.g., "Hs_6k").
#' @param max_prop_loss Numeric, FOV QC parameter.
#' @param max_totalcounts_loss Numeric, FOV QC parameter.
#' @param do_FOV_boundary Logical, perform FOV boundary QC.
#' @param SplitRatioToLocalThreshold Numeric, threshold for FOV boundary size ratio.
#' @param do_regional Logical, perform regional-level cell QC (SBR).
#' @param bandwidth Numeric, bandwidth for spatial smoothing.
#' @param weight_cutoff Numeric, cutoff for spatial weighting.
#'
#' @return A list containing filtered counts, metadata, xy, negcounts, a detailed flag logic table, and generated plots.
#' @export
run_spatial_qc <- function(counts, 
                           metadata, 
                           xy = NULL, 
                           negcounts = NULL,
                           do_nCount = TRUE, count_threshold = 20,
                           do_cellArea = TRUE, area_threshold = 30000,
                           do_runfovqc = TRUE, panel_name = "Hs_6k", max_prop_loss = 0.6, max_totalcounts_loss = 0.6,
                           do_FOV_boundary = TRUE, SplitRatioToLocalThreshold = 0.5,
                           do_regional = TRUE, bandwidth = 0.010, weight_cutoff = 0.080) {
  
  # Ensure required libraries are available
  required_pkgs <- c("ggplot2", "dplyr", "patchwork", "dbscan", "Matrix", "UpSetR", "purrr", "FNN", "pheatmap", "scales", "grid")
  missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
  if(length(missing_pkgs)) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  
  # Initialize storage
  plots <- list()
  flag_table <- data.frame(
    cell_id = metadata$cell_id,
    flag_nCount = FALSE,
    flag_area = FALSE,
    flag_fovqc = FALSE,
    flag_boundary = FALSE,
    flag_regional = FALSE,
    flag_overall = FALSE
  )
  
  # --- 1. Cell-level QC: nCount_RNA ---
  if (do_nCount) {
    message("Running nCount_RNA QC...")
    flag_table$flag_nCount <- metadata$nCount_RNA < count_threshold
    
    plots$nCount_hist <- ggplot2::ggplot(metadata, ggplot2::aes(x = nCount_RNA)) +
      ggplot2::geom_histogram(bins = 500, fill = "gray50") +
      ggplot2::geom_vline(xintercept = count_threshold, col = "red") +
      ggplot2::coord_cartesian(xlim = c(0, 1000)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "nCount_RNA", y = "Frequency", title = "Cell-level QC: nCount_RNA")
  }
  
  # --- 2. Cell Area QC ---
  if (do_cellArea) {
    message("Running Cell Area QC...")
    flag_table$flag_area <- metadata$Area > area_threshold
    
    plots$area_hist <- ggplot2::ggplot(metadata, ggplot2::aes(x = Area)) +
      ggplot2::geom_histogram(bins = 100, fill = "gray50") +
      ggplot2::geom_vline(xintercept = area_threshold, col = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Cell Area", y = "Frequency", title = "Cell Area QC")
  }
  
  # --- 3. FOV-level QC ---
  if (do_runfovqc) {
    if (is.null(xy)) stop("xy coordinates are required when do_runfovqc = TRUE.")
    message("Running FOV-level QC...")
    
    # Use the internal package data instead of reading from URL
    if(!panel_name %in% names(barcodes_by_panel)) stop("Panel name not found in barcodes list. Ensure 'barcodes_by_panel' is loaded.")
    barcodemap <- barcodes_by_panel[[panel_name]]
    
    fovqcresult <- runFOVQC(counts = counts, xy = xy, fov = metadata$FOV, 
                            barcodemap = barcodemap, 
                            max_prop_loss = max_prop_loss, 
                            max_totalcounts_loss = max_totalcounts_loss)
    
    flag_table$flag_fovqc <- metadata$FOV %in% fovqcresult$flaggedfovs
    
    plots$fov_mapFlagged <- recordPlot(mapFlaggedFOVs(fovqcresult))
    plots$fov_signalLoss <- recordPlot(FOVSignalLossSpatialPlot(fovqcresult))
    if(length(fovqcresult$flaggedfovs_forbias) > 0) {
      plots$fov_effectsHeatmap <- recordPlot(FOVEffectsHeatmap(fovqcresult))
    }
  }
  
  # --- 4. Cells Near FOV Borders QC ---
  if (do_FOV_boundary) {
    message("Running FOV Boundary QC...")
    flag_table$flag_boundary <- (metadata$SplitRatioToLocal > 0) & (metadata$SplitRatioToLocal < SplitRatioToLocalThreshold)
    
    p_detail <- ggplot2::ggplot(
      data = metadata[metadata$SplitRatioToLocal != 0,], 
      ggplot2::aes(x = log2(SplitRatioToLocal))) + 
      ggplot2::geom_histogram(bins = 100, fill = "steelblue") +
      ggplot2::labs(x = expression(log[2]("SplitRatioToLocal")), y = "Number of FOV border cells") + 
      ggplot2::theme_minimal() + 
      ggplot2::geom_vline(xintercept = log2(SplitRatioToLocalThreshold), color = "darkorange", lty = "dashed")
    
    pie_data <- metadata %>%
      dplyr::mutate(Category = dplyr::if_else(SplitRatioToLocal > 0, "Border cell", "Non-border cell")) %>%
      dplyr::count(Category) %>%
      dplyr::mutate(Percentage = (n / sum(n)) * 100, Label = paste0(Category, " (", round(Percentage, 1), "%)"))
    
    p_pie <- ggplot2::ggplot(pie_data, ggplot2::aes(x = 2, y = n, fill = Category)) +
      ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::scale_fill_manual(values = c("Border cell" = "steelblue", "Non-border cell" = "gray90")) +
      ggplot2::geom_text(ggplot2::aes(label = Label), position = ggplot2::position_stack(vjust = 0.5), size = 3, fontface = "bold") +
      ggplot2::theme_void() +
      ggplot2::xlim(0.5, 2.5) + 
      ggplot2::theme(legend.position = "none")
    
    plots$fov_boundary <- p_detail + patchwork::inset_element(p_pie, left = 0.6, bottom = 0.6, right = 1.0, top = 1.0)
  }
  
  # --- 5. Regional-level cell QC ---
  if (do_regional) {
    if (is.null(xy) || is.null(negcounts)) stop("xy and negcounts are required when do_regional = TRUE.")
    message("Running Regional-level Cell QC...")
    
    max_dist <- sqrt(-2 * (bandwidth^2) * log(weight_cutoff))
    nn <- dbscan::frNN(xy, eps = max_dist)
    
    n_neighbors <- sapply(nn$id, length)
    i_idx <- rep(1:nrow(xy), times = n_neighbors) 
    j_idx <- unlist(nn$id)                        
    
    distances <- unlist(nn$dist)
    weights <- exp(-(distances^2) / (2 * bandwidth^2))
    
    connectivities <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = weights, dims = c(nrow(xy), nrow(xy)))
    Matrix::diag(connectivities) <- 1
    
    row_sums <- Matrix::rowSums(connectivities)
    connectivities <- Matrix::Diagonal(x = 1 / row_sums) %*% connectivities
    
    counts_neg_spatially_smooth <- connectivities %*% negcounts
    counts_spatially_smooth <- connectivities %*% counts 
    
    smoothed_mean_neg <- as.vector(Matrix::rowMeans(counts_neg_spatially_smooth))
    smoothed_mean_target <- as.vector(Matrix::rowMeans(counts_spatially_smooth))
    
    epsilon <- 1e-9
    metadata$SBR <- smoothed_mean_target / (smoothed_mean_neg + epsilon)
    
    flag_table$flag_regional <- log2(metadata$SBR) < 0
    
    sbr_breaks <- c(-Inf, -1, 0, 1, 2, 3, Inf)
    sbr_labels <- c("<= -1", "-1 to 0", "0 to 1", "1 to 2", "2 to 3", "> 3")
    sbr_labels <- paste0(sbr_labels, paste0(" (", round(100 * table(cut(log2(metadata$SBR), breaks = sbr_breaks)) / nrow(metadata), 1), "%)"))
    
    manual_colors <- c("<= -1" = "#D73027", "-1 to 0" = "#F46D43", "0 to 1" = "#D9EF8B", "1 to 2" = "#A6D96A", "2 to 3" = "#1A9850", "> 3" = "#006837")
    names(manual_colors) <- sbr_labels
    
    plots$regional_sbr <- ggplot2::ggplot(data = metadata) +
      ggplot2::geom_point(ggplot2::aes(x = x_slide_mm, y = y_slide_mm, color = cut(log2(SBR), breaks = sbr_breaks, labels = sbr_labels)), size = 0.001, alpha = 1) +
      ggplot2::scale_color_manual(values = manual_colors) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 7, alpha = 1))) +
      ggplot2::labs(color = "log2(SBR)") +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw() + 
      ggplot2::facet_wrap(~Run_Tissue_name)
  }
  
  # --- 6. Combine Flags and UpSet Plot ---
  flag_table$flag_overall <- flag_table$flag_nCount | flag_table$flag_area | flag_table$flag_fovqc | flag_table$flag_boundary | flag_table$flag_regional
  
  filter_list <- list(
    `Low Counts` = flag_table$cell_id[flag_table$flag_nCount],
    `High Cell Area` = flag_table$cell_id[flag_table$flag_area],
    `FOV QC` = flag_table$cell_id[flag_table$flag_fovqc],
    `Low SplitRatio` = flag_table$cell_id[flag_table$flag_boundary],
    `Low SBR` = flag_table$cell_id[flag_table$flag_regional]
  )
  filter_list <- purrr::keep(filter_list, ~ length(.) > 0)
  
  if (length(filter_list) > 1) {
    u_plot <- UpSetR::upset(UpSetR::fromList(filter_list), nintersects = 10, order.by = "freq", nsets = length(filter_list), text.scale = 1.5)
    
    pie_data_overall <- data.frame(Category = dplyr::if_else(flag_table$flag_overall, "Flagged", "Kept")) %>%
      dplyr::count(Category) %>%
      dplyr::mutate(
        Percentage = (n / sum(n)) * 100, 
        Label = paste0(Category, "\n(", round(Percentage, 1), "%)")
      )
    
    p_pie_overall <- ggplot2::ggplot(pie_data_overall, ggplot2::aes(x = 2, y = n, fill = Category)) +
      ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::scale_fill_manual(values = c("Flagged" = "#D73027", "Kept" = "gray90")) +
      ggplot2::geom_text(ggplot2::aes(label = Label), position = ggplot2::position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
      ggplot2::theme_void() +
      ggplot2::xlim(0.5, 2.5) + 
      ggplot2::theme(legend.position = "none")
    
    plots$upset_plot <- patchwork::wrap_elements(grid::grid.grabExpr(print(u_plot))) +
      patchwork::inset_element(p_pie_overall, left = 0.65, bottom = 0.65, right = 1.0, top = 1.0)
  }
  
  # --- 7. Subset and Return ---
  keep_idx <- !flag_table$flag_overall
  
  return(list(
    counts = counts[keep_idx, , drop = FALSE],
    metadata = metadata[keep_idx, , drop = FALSE],
    xy = if(!is.null(xy)) xy[keep_idx, , drop = FALSE] else NULL,
    negcounts = if(!is.null(negcounts)) negcounts[keep_idx, , drop = FALSE] else NULL,
    flag_table = flag_table,
    plots = plots
  ))
}
