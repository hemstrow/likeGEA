
plot_migration <- function(M, data, coordinates,
                           threshold_value = 0.4,
                           topn = 2,
                           bin_size = 1){

  DistM = data$DistM

  ## N_ind_per_patch comes from the migration data (already in coordinates order)
  N_ind_per_patch = data$N_table

  ## ---- Build migration network -----------------------------------------------

  diag(M) <- 0
  n = nrow(M)

  migration_lines <- data.frame()
  for (i in 1:n) {
    row_values <- M[i, ]
    row_values[i] <- NA

    topn_indices <- order(row_values, decreasing = TRUE)[1:topn]
    threshold_indices <- which(M[i, ] >= threshold_value)
    union_indices <- union(topn_indices, threshold_indices)

    for (j in union_indices) {
      migration_lines <- rbind(migration_lines, data.frame(
        patch_from = coordinates$patch[i],
        patch_to   = coordinates$patch[j],
        lon_from   = coordinates$lon[i],
        lat_from   = coordinates$lat[i],
        lon_to     = coordinates$lon[j],
        lat_to     = coordinates$lat[j],
        migration_rate = M[i, j],
        mig_rate = M[i, j]
      ))
    }
  }

  migration_lines <- migration_lines |>
    dplyr::mutate(p1 = pmin(patch_from, patch_to),
                  p2 = pmax(patch_from, patch_to)) |>
    dplyr::distinct(p1, p2, .keep_all = TRUE)

  migration_lines <- migration_lines |> dplyr::arrange(migration_rate)

  ## ---- Reduce dense connections ----------------------------------------------

  migration_lines <- migration_lines |>
    dplyr::mutate(mid_lon = (lon_from + lon_to) / 2,
           mid_lat = (lat_from + lat_to) / 2)

  migration_lines <- migration_lines |>
    dplyr::mutate(mid_lon_bin = round(mid_lon / bin_size) * bin_size,
           mid_lat_bin = round(mid_lat / bin_size) * bin_size)

  migration_lines <- migration_lines |>
    dplyr::group_by(mid_lon_bin, mid_lat_bin) |>
    dplyr::filter(migration_rate == max(migration_rate)) |>
    dplyr::ungroup()

  ## ---- Map plot --------------------------------------------------------------

  ## Bounding box: pad by the larger of 20% of the coordinate range or 0.5 degrees,
  ## so that small study areas still show enough geographic context.
  lon_range <- range(coordinates$lon)
  lat_range <- range(coordinates$lat)
  pad_lon   <- max(diff(lon_range) * 0.20, 0.5)
  pad_lat   <- max(diff(lat_range) * 0.20, 0.5)
  xlim      <- lon_range + c(-pad_lon, pad_lon)
  ylim      <- lat_range + c(-pad_lat, pad_lat)

  ## Political boundaries (vector data — loads instantly from rnaturalearth)
  world  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  sf::sf_use_s2(FALSE)   # disable spherical geometry to avoid invalid-edge errors on ne_states
  states <- rnaturalearth::ne_states(returnclass = "sf") |>
    sf::st_crop(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
  sf::sf_use_s2(TRUE)    # re-enable for everything else

  ## N_ind_per_patch already loaded from data$N_table above

  map <- ggplot2::ggplot() +
    ## Country fill + outlines
    ggplot2::geom_sf(data = world, fill = "grey92", color = "grey40", linewidth = 0.4) +
    ## State / province borders for within-country geographic context
    ggplot2::geom_sf(data = states, fill = NA, color = "grey60", linewidth = 0.25, linetype = "dashed") +
    ## Migration network
    ggplot2::geom_curve(
      data = migration_lines,
      ggplot2::aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to,
                   linewidth = mig_rate, color = migration_rate),
      curvature = 0.2, alpha = 0.8,
      show.legend = c(linewidth = FALSE)
    ) +
    ggplot2::scale_linewidth_continuous(trans = "log10", range = c(0.1, 2)) +
    ggplot2::scale_color_gradient(name = "Rel. Migration Rate", low = "yellow", high = "red", trans = "log10") +
    ## Sampling sites
    ggplot2::geom_point(
      data = coordinates,
      ggplot2::aes(x = lon, y = lat, size = N_ind_per_patch),
      fill = "grey50", color = "black", shape = 21, stroke = 1
    ) +
    ggplot2::scale_size_continuous(range = c(0.5, 4)) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(
      x = "Longitude",
      y = "Latitude",
      size = "Individuals \n Sampled"
    ) +
    ggplot2::guides(
      size  = ggplot2::guide_legend(title.position = "top", title.hjust = 0.5),
      color = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5),
      linewidth = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )

  ## ---- Isolation by distance plot --------------------------------------------

  dist_vec <- as.vector(DistM)
  mig_vec <- as.vector(M)
  valid_idx <- which(dist_vec != 0)
  dist_vec <- dist_vec[valid_idx]
  mig_vec <- mig_vec[valid_idx]

  Dist_mig <- ggplot2::ggplot(data.frame(Distance = dist_vec*1e-3, Migration = mig_vec),
                              ggplot2::aes(x = Distance, y = log(Migration))) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Distance (km)", y = "log(Rel. Migration Rate)") +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black")
    )

  ## ---- Combine and save ------------------------------------------------------

  combined_plot <- cowplot::plot_grid(
    Dist_mig, map,
    nrow = 1,
    rel_widths = c(2.8, 3)
  )

  return(combined_plot)
}
