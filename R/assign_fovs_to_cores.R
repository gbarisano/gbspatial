#' Assign FOVs to TMA Cores with Strict Geographic Lock and Conditional Cohesion
#'
#' @description 
#' Maps FOVs to a TMA grid. Drops anchors that drift too far or capture insufficient cells. 
#' Enforces a Strict Geographic Lock: FOVs fully contained within a single grid cell are 
#' assigned purely based on geometry. FOVs that straddle grid lines receive a provisional 
#' cell-majority vote followed by a Spatial Cohesion Pass to group them with contiguous tissue.
#'
#' @param fov_input Data source for FOVs. Must contain `x_global_px`, `y_global_px`, and `FOV`.
#' @param cell_input Data source for Cells. Must contain `CenterX_global_px`, `CenterY_global_px`, and `FOV`.
#' @param tma_map_input A matrix/dataframe representing the physical TMA layout.
#' @param fov_size Numeric. The size of the FOV in pixels. Default is `4256`.
#' @param core_drift_tolerance Numeric (0 to 1). Max fraction of grid spacing an anchor can drift. Default `0.40`.
#' @param min_cells Numeric. Minimum number of cells required to validate a core anchor. Default `50`.
#' @param slidelabels Optional character vector of slide names.
#'
#' @return A named list containing `mapped_data` and `plots`.
#' 
#' @importFrom dplyr mutate group_by ungroup bind_rows row_number select summarise left_join coalesce slice_sample filter slice_max
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_rect geom_text geom_point geom_hline geom_vline aes coord_equal theme_minimal labs theme scale_x_continuous scale_y_continuous expansion element_text element_blank
#' @importFrom stats dist
#' @importFrom utils read.csv
#' @export
assign_fovs_to_cores <- function(fov_input, cell_input, tma_map_input, fov_size = 4256, 
                                 core_drift_tolerance = 0.40, min_cells = 50, slidelabels = NULL) {
  
  # 1. Standardize Inputs
  if (is.data.frame(fov_input) || is.matrix(fov_input)) fov_input <- list(fov_input)
  else fov_input <- as.list(fov_input)
  
  if (is.data.frame(cell_input) || is.matrix(cell_input)) cell_input <- list(cell_input)
  else cell_input <- as.list(cell_input)
  
  n_slides <- length(fov_input)
  if (length(cell_input) != n_slides) stop("Number of FOV inputs must match cell inputs.")
  
  if (is.data.frame(tma_map_input) || is.matrix(tma_map_input)) tma_map_input <- list(tma_map_input)
  else tma_map_input <- as.list(tma_map_input)
  
  if (length(tma_map_input) == 1 && n_slides > 1) tma_map_input <- rep(tma_map_input, n_slides)
  if (n_slides != length(tma_map_input)) stop("TMA map inputs must be 1 or match the number of FOV/cell inputs.")
  
  if (!is.null(slidelabels)) {
    if (length(slidelabels) != n_slides) stop("Length of 'slidelabels' must match inputs.")
    slide_names <- slidelabels
  } else if (!is.null(names(fov_input)) && all(names(fov_input) != "")) {
    slide_names <- names(fov_input)
  } else {
    slide_names <- paste0("s", seq_len(n_slides))
  }
  
  all_mapped_data <- list()
  all_plots <- list()
  
  for (i in seq_len(n_slides)) {
    current_slide <- slide_names[i]
    
    if (is.character(fov_input[[i]])) fov <- utils::read.csv(fov_input[[i]]) else fov <- as.data.frame(fov_input[[i]])
    if (is.character(cell_input[[i]])) cells <- utils::read.csv(cell_input[[i]]) else cells <- as.data.frame(cell_input[[i]])
    
    # Clean whitespace aggressively right out of the gate
    fov$original_FOV <- trimws(as.character(fov$FOV))
    cells$original_FOV <- trimws(as.character(cells$fov))
    
    fov$slidename <- trimws(as.character(current_slide))
    cells$slidename <- trimws(as.character(current_slide))
    
    # Create composite ID
    fov$FOV <- paste(fov$slidename, fov$original_FOV, sep = "_")
    cells$FOV <- paste(cells$slidename, cells$original_FOV, sep = "_")
    
    tma_map_df <- as.data.frame(tma_map_input[[i]])
    n_rows <- nrow(tma_map_df)
    n_cols <- ncol(tma_map_df)
    
    # 2. Establish Grid from True FOV Coordinates
    xmin_global <- min(fov$x_global_px, na.rm = TRUE)
    xmax_global <- max(fov$x_global_px, na.rm = TRUE) + fov_size
    ymax_global <- max(fov$y_global_px, na.rm = TRUE)
    ymin_global <- min(fov$y_global_px, na.rm = TRUE) - fov_size 
    
    xinter <- (xmax_global - xmin_global) / n_cols
    yinter <- (ymax_global - ymin_global) / n_rows
    
    v_lines <- seq(xmin_global, xmax_global, by = xinter)
    h_lines <- seq(ymin_global, ymax_global, by = yinter)
    
    # The true unadjusted grid anchors
    grid_anchors <- expand.grid(core_col = 1:n_cols, core_row = 1:n_rows) |>
      dplyr::mutate(
        anchor_x = xmin_global + (.data$core_col - 0.5) * xinter,
        anchor_y = ymax_global - (.data$core_row - 0.5) * yinter,
        core_str = paste0("C", .data$core_col, "R", .data$core_row)
      )
    
    # 3. Fast Initial Cell Assignment
    dist_mat <- matrix(NA, nrow = nrow(cells), ncol = nrow(grid_anchors))
    for (k in 1:nrow(grid_anchors)) {
      dist_mat[, k] <- (cells$CenterX_global_px - grid_anchors$anchor_x[k])^2 + 
                       (cells$CenterY_global_px - grid_anchors$anchor_y[k])^2
    }
    cells$initial_idx <- max.col(-dist_mat, ties.method = "first")
    cells$core_str <- grid_anchors$core_str[cells$initial_idx]
    
    # 4. Refine Anchors AND Validate
    refined_anchors <- cells |>
      dplyr::group_by(.data$core_str) |>
      dplyr::summarise(
        refined_x = mean(.data$CenterX_global_px, na.rm = TRUE),
        refined_y = mean(.data$CenterY_global_px, na.rm = TRUE),
        cell_count = dplyr::n(),
        .groups = "drop"
      )
    
    final_anchors <- grid_anchors |>
      dplyr::left_join(refined_anchors, by = "core_str") |>
      dplyr::mutate(
        drift_x = abs(.data$refined_x - .data$anchor_x),
        drift_y = abs(.data$refined_y - .data$anchor_y),
        is_valid = !is.na(.data$refined_x) & 
                   .data$cell_count >= min_cells &
                   .data$drift_x <= (xinter * core_drift_tolerance) &
                   .data$drift_y <= (yinter * core_drift_tolerance)
      ) |>
      dplyr::filter(.data$is_valid) |>
      dplyr::mutate(
        final_x = .data$refined_x,
        final_y = .data$refined_y
      )
    
    if (nrow(final_anchors) == 0) stop(paste("No valid cores found for slide", current_slide))
    
    # 5. Provisional Cell-Majority Voting for Straddling FOVs
    cell_dist_mat2 <- matrix(NA, nrow = nrow(cells), ncol = nrow(final_anchors))
    for (k in 1:nrow(final_anchors)) {
      cell_dist_mat2[, k] <- (cells$CenterX_global_px - final_anchors$final_x[k])^2 + 
                             (cells$CenterY_global_px - final_anchors$final_y[k])^2
    }
    cells$closest_valid_idx <- max.col(-cell_dist_mat2, ties.method = "first")
    cells$valid_core_str <- final_anchors$core_str[cells$closest_valid_idx]
    
    fov_votes <- cells |>
      dplyr::group_by(.data$FOV, .data$valid_core_str) |>
      dplyr::summarise(vote_count = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(.data$FOV) |>
      dplyr::slice_max(order_by = .data$vote_count, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
    
    fov <- fov |>
      dplyr::mutate(
        fov_xmin = .data$x_global_px,
        fov_xmax = .data$x_global_px + fov_size,
        fov_ymax = .data$y_global_px,
        fov_ymin = .data$y_global_px - fov_size,
        fov_xcenter = .data$x_global_px + (fov_size / 2),
        fov_ycenter = .data$y_global_px - (fov_size / 2)
      )
    
    # 6. Flag Grid Straddlers & Calculate Strict Geographic Assignment
    crosses_v <- sapply(1:nrow(fov), function(j) {
      any(v_lines > fov$fov_xmin[j] + 1e-3 & v_lines < fov$fov_xmax[j] - 1e-3)
    })
    crosses_h <- sapply(1:nrow(fov), function(j) {
      any(h_lines > fov$fov_ymin[j] + 1e-3 & h_lines < fov$fov_ymax[j] - 1e-3)
    })
    fov$straddles_grid <- crosses_v | crosses_h
    
    # Mathematical grid assignment (for FOVs that don't straddle)
    geo_dist_mat <- matrix(NA, nrow = nrow(fov), ncol = nrow(grid_anchors))
    for (k in 1:nrow(grid_anchors)) {
      geo_dist_mat[, k] <- (fov$fov_xcenter - grid_anchors$anchor_x[k])^2 + 
                           (fov$fov_ycenter - grid_anchors$anchor_y[k])^2
    }
    fov$geometric_core_str <- grid_anchors$core_str[max.col(-geo_dist_mat, ties.method = "first")]
    
    # Fallback to nearest valid anchor for provisional voting
    fov_dist_mat <- matrix(NA, nrow = nrow(fov), ncol = nrow(final_anchors))
    for (k in 1:nrow(final_anchors)) {
      fov_dist_mat[, k] <- (fov$fov_xcenter - final_anchors$final_x[k])^2 + 
                           (fov$fov_ycenter - final_anchors$final_y[k])^2
    }
    fov$fallback_idx <- max.col(-fov_dist_mat, ties.method = "first")
    fov$fallback_core_str <- final_anchors$core_str[fov$fallback_idx]
    
    fov <- fov |>
      dplyr::left_join(fov_votes |> dplyr::select(.data$FOV, .data$valid_core_str), by = "FOV") |>
      dplyr::mutate(
        prov_core_str = dplyr::coalesce(.data$valid_core_str, .data$fallback_core_str)
      )
      
    # 7. Strict Lock & Conditional Spatial Cohesion Pass
    fov_centers <- cells |>
      dplyr::group_by(.data$FOV) |>
      dplyr::summarise(
        c_x = mean(.data$CenterX_global_px, na.rm = TRUE),
        c_y = mean(.data$CenterY_global_px, na.rm = TRUE),
        .groups = "drop"
      )
      
    fov <- fov |>
      dplyr::left_join(fov_centers, by = "FOV") |>
      dplyr::mutate(
        c_x = dplyr::coalesce(.data$c_x, .data$fov_xcenter),
        c_y = dplyr::coalesce(.data$c_y, .data$fov_ycenter)
      )
      
    fov_coords <- as.matrix(fov[, c("c_x", "c_y")])
    dist_mat_fovs <- as.matrix(stats::dist(fov_coords))
    neighborhood_radius <- fov_size * 1.6 
    
    smoothed_core_str <- character(nrow(fov))
    for (k in 1:nrow(fov)) {
      if (fov$straddles_grid[k]) {
        # Straddlers: Apply neighborhood smoothing based on provisional core
        neighbor_idx <- which(dist_mat_fovs[k, ] <= neighborhood_radius)
        neighbor_cores <- fov$prov_core_str[neighbor_idx]
        neighbor_cores <- neighbor_cores[!is.na(neighbor_cores)]
        
        if (length(neighbor_cores) > 0) {
          tally <- table(neighbor_cores)
          smoothed_core_str[k] <- names(tally)[which.max(tally)]
        } else {
          smoothed_core_str[k] <- fov$prov_core_str[k]
        }
      } else {
        # Strict Geographic Lock: Fully contained FOVs are locked to their mathematical grid
        smoothed_core_str[k] <- fov$geometric_core_str[k]
      }
    }
    
    fov$assigned_core_str <- smoothed_core_str

    # 8. Resolve Final Assignments
    fov <- fov |>
      dplyr::left_join(
        final_anchors |> dplyr::select(.data$core_str, .data$core_col, .data$core_row, .data$final_x, .data$final_y), 
        by = c("assigned_core_str" = "core_str")
      )
    
    max_fov_dist <- sqrt(xinter^2 + yinter^2) * 0.90 
    
    fov <- fov |>
      dplyr::mutate(
        dist_to_anchor = sqrt((.data$fov_xcenter - .data$final_x)^2 + (.data$fov_ycenter - .data$final_y)^2),
        core_col = ifelse(.data$dist_to_anchor <= max_fov_dist, .data$core_col, NA),
        core_row = ifelse(.data$dist_to_anchor <= max_fov_dist, .data$core_row, NA),
        core_str = ifelse(.data$dist_to_anchor <= max_fov_dist, .data$assigned_core_str, NA)
      )
    
    # 9. Merge with Clinical IDs
    colnames(tma_map_df) <- 1:n_cols
    tma_map_long <- tma_map_df |>
      dplyr::mutate(core_row = rev(dplyr::row_number())) |> 
      tidyr::pivot_longer(cols = -.data$core_row, names_to = "core_col", values_to = "id") |>
      dplyr::mutate(
        core_col = as.numeric(.data$core_col),
        core_str = paste0("C", .data$core_col, "R", .data$core_row)
      )
    
    res <- merge(fov, tma_map_long, by = c("core_col", "core_row", "core_str"), all.x = TRUE)

    
    # 10. Assign Cells
    cells_merged <- merge(cells, res[, c("original_FOV", "core_str")], by.x = "FOV", by.y = "original_FOV", all.x = TRUE)
    
    plot_list <- list()
    
    # 11. Prepare axes formatting for labels
    x_breaks <- xmin_global + (1:n_cols - 0.5) * xinter
    y_breaks <- ymax_global - (1:n_rows - 0.5) * yinter
    x_labels <- paste0("C", 1:n_cols)
    y_labels <- paste0("R", 1:n_rows)
    
    sampled_cells_whole <- cells_merged |>
      dplyr::group_by(.data$FOV) |>
      dplyr::slice_sample(n = 200) |> 
      dplyr::ungroup()
      
    # 12. Whole-Slide Plot (Colored by Column)
    p_col <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = h_lines, color = "grey40", linetype = "dashed", alpha = 0.6) +
      ggplot2::geom_vline(xintercept = v_lines, color = "grey40", linetype = "dashed", alpha = 0.6) +
      ggplot2::geom_point(
        data = sampled_cells_whole,
        ggplot2::aes(x = .data$CenterX_global_px, y = .data$CenterY_global_px), 
        color = "grey60", size = 0.2, alpha = 0.4, na.rm = TRUE
      ) +
      ggplot2::geom_rect(
        data = res,
        ggplot2::aes(xmin = .data$fov_xmin, xmax = .data$fov_xmax, 
                     ymin = .data$fov_ymin, ymax = .data$fov_ymax, 
                     fill = factor(.data$core_col)), 
        color = "black", alpha = 0.3
      ) +
      ggplot2::geom_text(
        data = res,
        ggplot2::aes(x = .data$fov_xcenter, y = .data$fov_ycenter, 
                     label = ifelse(is.na(.data$core_str), 
                                    paste0(.data$original_FOV, "\n[NA]"), 
                                    paste0(.data$original_FOV, "\n[", .data$core_str, "]"))), 
        color = "black", size = 2.2, fontface = "bold"
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = ggplot2::expansion(mult = 0.05)) +
      ggplot2::scale_y_continuous(breaks = y_breaks, labels = y_labels, expand = ggplot2::expansion(mult = 0.05)) +
      ggplot2::coord_equal() +  
      ggplot2::theme_minimal() + 
      ggplot2::labs(title = paste("Whole Slide Overview (By Column):", current_slide),
                    subtitle = "Cells: Grey | Colors: Column Assignments") +
      ggplot2::theme(legend.position = "none",
                     axis.text = ggplot2::element_text(face = "bold", size = 10),
                     axis.title = ggplot2::element_blank()) 
    
    plot_list[["Whole_Slide_by_Column"]] <- p_col
    
    # 13. Whole-Slide Plot (Colored by Row)
    p_row <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = h_lines, color = "grey40", linetype = "dashed", alpha = 0.6) +
      ggplot2::geom_vline(xintercept = v_lines, color = "grey40", linetype = "dashed", alpha = 0.6) +
      ggplot2::geom_point(
        data = sampled_cells_whole,
        ggplot2::aes(x = .data$CenterX_global_px, y = .data$CenterY_global_px), 
        color = "grey60", size = 0.2, alpha = 0.4, na.rm = TRUE
      ) +
      ggplot2::geom_rect(
        data = res,
        ggplot2::aes(xmin = .data$fov_xmin, xmax = .data$fov_xmax, 
                     ymin = .data$fov_ymin, ymax = .data$fov_ymax, 
                     fill = factor(.data$core_row)), 
        color = "black", alpha = 0.3
      ) +
      ggplot2::geom_text(
        data = res,
        ggplot2::aes(x = .data$fov_xcenter, y = .data$fov_ycenter, 
                     label = ifelse(is.na(.data$core_str), 
                                    paste0(.data$original_FOV, "\n[NA]"), 
                                    paste0(.data$original_FOV, "\n[", .data$core_str, "]"))), 
        color = "black", size = 2.2, fontface = "bold"
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = ggplot2::expansion(mult = 0.05)) +
      ggplot2::scale_y_continuous(breaks = y_breaks, labels = y_labels, expand = ggplot2::expansion(mult = 0.05)) +
      ggplot2::coord_equal() +  
      ggplot2::theme_minimal() + 
      ggplot2::labs(title = paste("Whole Slide Overview (By Row):", current_slide),
                    subtitle = "Cells: Grey | Colors: Row Assignments") +
      ggplot2::theme(legend.position = "none",
                     axis.text = ggplot2::element_text(face = "bold", size = 10),
                     axis.title = ggplot2::element_blank()) 
    
    plot_list[["Whole_Slide_by_Row"]] <- p_row
    
    all_mapped_data[[current_slide]] <- res
    all_plots[[current_slide]] <- plot_list
  }
  
  return(list(mapped_data = dplyr::bind_rows(all_mapped_data), plots = all_plots))
}
