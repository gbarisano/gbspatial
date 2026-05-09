# Parts of this code have Copyright ©2024 Bruker Spatial Biology, Inc. All rights reserved. Subject to additional license terms and conditions provided separately by Bruker Spatial Biology, Inc.
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
#' @return A list containing filtered counts, metadata, xy, a detailed flag logic table, and generated plots.
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
  required_pkgs <- c("ggplot2", "dplyr", "patchwork", "dbscan", "Matrix", "UpSetR", "purrr")
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
    
    # Source external scripts
    #source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/Main/_code/FOV%20QC/FOV%20QC%20utils.R")
    allbarcodes
    
    if(!panel_name %in% names(allbarcodes)) stop("Panel name not found in barcodes list.")
    barcodemap <- allbarcodes[[panel_name]]
    
    fovqcresult <- runFOVQC(counts = counts, xy = xy, fov = metadata$FOV, 
                            barcodemap = barcodemap, 
                            max_prop_loss = max_prop_loss, 
                            max_totalcounts_loss = max_totalcounts_loss)
    
    flag_table$flag_fovqc <- metadata$FOV %in% fovqcresult$flaggedfovs
    
    # Capture FOV plots (Note: Depending on how runFOVQC is written, these might output to the plotting device directly)
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
    plots$upset_plot <- recordPlot(
      UpSetR::upset(UpSetR::fromList(filter_list), nintersects = 10, order.by = "freq", nsets = length(filter_list), text.scale = 1.5)
    )
  }
  
  # --- 7. Subset and Return ---
  keep_idx <- !flag_table$flag_overall
  
  return(list(
    counts = counts[keep_idx, , drop = FALSE],
    metadata = metadata[keep_idx, , drop = FALSE],
    xy = if(!is.null(xy)) xy[keep_idx, , drop = FALSE] else NULL,
    flag_table = flag_table,
    plots = plots
  ))
}
