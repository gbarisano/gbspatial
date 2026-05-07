#' Assign FOVs to TMA Cores
#'
#' @description 
#' Maps Field of View (FOV) coordinates to a Tissue Microarray (TMA) grid matrix. 
#' Groups adjacent FOVs into spatial cores, aligns them to physical slide dimensions, 
#' and merges them with clinical IDs from a TMA map. Supports batch processing of 
#' lists or vectors of inputs.
#'
#' @param fov_input A file path (character string), dataframe, or list of file paths/dataframes containing FOV coordinate data. Must contain `x_global_px`, `y_global_px`, and `FOV` columns.
#' @param tma_map_input A matrix, dataframe, or list of matrices/dataframes representing the physical TMA layout. If a single matrix is provided alongside multiple `fov_input`s, it will be automatically replicated.
#' @param fov_size Numeric. The size of the FOV in pixels. Default is `4256`.
#' @param tol Numeric. Spatial tolerance in pixels for clustering adjacent FOVs. Default is `150`.
#' @param slidelabels Optional character vector of slide names. Length must match the number of FOV inputs.
#'
#' @return A named list containing two elements:
#' \describe{
#'   \item{mapped_data}{A single combined dataframe containing all merged FOV and core assignments.}
#'   \item{plots}{A nested list of ggplot objects grouped by slide name and TMA row.}
#' }
#' 
#' @importFrom dplyr mutate group_by ungroup bind_rows row_number
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_rect geom_text aes coord_equal theme_minimal labs theme
#' @importFrom stats dist hclust cutree
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' \dontrun{
#' # Process a single slide
#' res <- assign_fovs_to_cores("slide1_fovs.csv", tma_matrix)
#' 
#' # Process a batch of slides with a single replicated TMA map
#' res_batch <- assign_fovs_to_cores(
#'   fov_input = c("s1.csv", "s2.csv", "s3.csv"), 
#'   tma_map_input = tma_matrix,
#'   slidelabels = c("Patient_A", "Patient_B", "Patient_C")
#' )
#' }
assign_fovs_to_cores <- function(fov_input, tma_map_input, fov_size = 4256, tol = 150, slidelabels = NULL) {
  
  # 1. Standardize fov_input
  if (is.data.frame(fov_input) || is.matrix(fov_input)) {
    fov_input <- list(fov_input)
  } else {
    fov_input <- as.list(fov_input)
  }
  
  n_slides <- length(fov_input)
  
  # 2. Standardize tma_map_input
  if (is.data.frame(tma_map_input) || is.matrix(tma_map_input)) {
    tma_map_input <- list(tma_map_input)
  } else {
    tma_map_input <- as.list(tma_map_input)
  }
  
  # 3. Handle single map for multiple slides
  if (length(tma_map_input) == 1 && n_slides > 1) {
    tma_map_input <- rep(tma_map_input, n_slides)
  }
  
  if (n_slides != length(tma_map_input)) {
    stop("The number of TMA map inputs must be 1, or exactly match the number of FOV inputs.")
  }
  
  # 4. Determine slidenames based on precedence rules
  if (!is.null(slidelabels)) {
    if (length(slidelabels) != n_slides) {
      stop("Length of 'slidelabels' must match the number of FOV inputs.")
    }
    slide_names <- slidelabels
  } else if (!is.null(names(fov_input)) && all(names(fov_input) != "")) {
    slide_names <- names(fov_input)
  } else {
    slide_names <- paste0("s", seq_len(n_slides))
  }
  
  all_mapped_data <- list()
  all_plots <- list()
  
  # 5. Iterate through each slide/input pair
  for (i in seq_len(n_slides)) {
    
    current_slide <- slide_names[i]
    curr_fov_input <- fov_input[[i]]
    curr_map_input <- tma_map_input[[i]]
    
    # Read FOV data safely
    if (is.character(curr_fov_input)) {
      fov <- utils::read.csv(curr_fov_input)
    } else {
      fov <- as.data.frame(curr_fov_input)
    }
    
    fov$slidename <- current_slide
    
    tma_map_df <- as.data.frame(curr_map_input)
    n_rows <- nrow(tma_map_df)
    n_cols <- ncol(tma_map_df)
    
    # Group adjacent FOVs into Cores
    dist_mat <- stats::dist(fov[, c("x_global_px", "y_global_px")], method = "maximum")
    hc <- stats::hclust(dist_mat, method = "single")
    fov$core_cluster <- stats::cutree(hc, h = fov_size + tol)
    
    # Calculate Centers for each Core Cluster
    # Using the .data pronoun is best practice in packages to avoid R CMD check warnings
    fov <- fov |>
      dplyr::group_by(.data$core_cluster) |>
      dplyr::mutate(
        xcenter = (min(.data$x_global_px) + max(.data$x_global_px) + fov_size) / 2,
        ycenter = (min(.data$y_global_px) + max(.data$y_global_px) + fov_size) / 2,
        fov_core = paste(.data$FOV, collapse = "+") 
      ) |>
      dplyr::ungroup()
    
    # Grid Assignment
    xmin_global <- min(fov$x_global_px)
    xmax_global <- max(fov$x_global_px) + fov_size
    ymin_global <- min(fov$y_global_px)
    ymax_global <- max(fov$y_global_px) + fov_size
    
    xinter <- (xmax_global - xmin_global) / n_cols
    yinter <- (ymax_global - ymin_global) / n_rows
    
    fov <- fov |>
      dplyr::mutate(
        core_col = pmin(floor((.data$xcenter - xmin_global) / xinter) + 1, n_cols),
        core_row = pmin(floor((.data$ycenter - ymin_global) / yinter) + 1, n_rows),
        core_str = paste0("C", .data$core_col, "R", .data$core_row)
      )
    
    # Parse TMA matrix safely
    colnames(tma_map_df) <- 1:n_cols
    tma_map_long <- tma_map_df |>
      dplyr::mutate(core_row = rev(dplyr::row_number())) |> 
      tidyr::pivot_longer(cols = -.data$core_row, names_to = "core_col", values_to = "id") |>
      dplyr::mutate(
        core_col = as.numeric(.data$core_col),
        core_str = paste0("C", .data$core_col, "R", .data$core_row)
      )
    
    # Merge FOV data with Clinical IDs
    res <- merge(fov, tma_map_long, by = c("core_col", "core_row", "core_str"), all.x = TRUE)
    
    # Update FOV to be a unique ID
    res$original_FOV <- res$FOV
    res$FOV <- paste(res$slidename, res$original_FOV, sep = "_")
    
    # Generate Plots per TMA Row
    plot_list <- list()
    unique_rows <- sort(unique(res$core_row))
    
    for (c_row in unique_rows) {
      p <- ggplot2::ggplot(res[res$core_row == c_row, ]) +
        ggplot2::geom_rect(
          ggplot2::aes(xmin = .data$x_global_px, xmax = .data$x_global_px + fov_size, 
                       ymin = .data$y_global_px - fov_size, ymax = .data$y_global_px, 
                       fill = factor(.data$core_str)), 
          color = "black", alpha = 0.6
        ) +
        ggplot2::geom_text(
          ggplot2::aes(x = .data$x_global_px + (fov_size / 2), 
                       y = .data$y_global_px - (fov_size / 2), 
                       label = .data$original_FOV), 
          color = "black", size = 3.5, fontface = "bold"
        ) +
        ggplot2::coord_equal() +  
        ggplot2::theme_minimal() + 
        ggplot2::labs(
          title = paste("Slide:", current_slide, "| TMA Row:", c_row), 
          fill = "Core Assignment"
        ) +
        ggplot2::theme(legend.position = "right")
      
      plot_list[[paste0("Row_", c_row)]] <- p
    }
    
    all_mapped_data[[current_slide]] <- res
    all_plots[[current_slide]] <- plot_list
  }
  
  final_mapped_data <- dplyr::bind_rows(all_mapped_data)
  
  return(list(mapped_data = final_mapped_data, plots = all_plots))
}