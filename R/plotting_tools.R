#' Plot migration stan model results.
#'
#' Consider running via \code{\link{likeGEA}} for the full pipeline. Generally
#' for internal use or when partial re-runs are needed.
#'
#' @export
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
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
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

#' Plot migration core stan GEA results.
#'
#' Consider running via \code{\link{likeGEA}} for the full pipeline. Generally
#' for internal use or when partial re-runs are needed.
#'
#' @export
plot_gea <- function(stan_data,
                     migration_results,
                     beta_draw,
                     beta_summary,
                     rho_kappa_mu,
                     env_names,
                     CI = 0.95){

  # DEV NOTES: Plots should probably all be swapped over to ggplots and the raw data
  # returned to allow for much easier customization of plots!

  ## ---- Load -----------------------------------------------------------------
  data         <- stan_data; rm(stan_data)
  Migration    <- migration_results; rm(migration_results)

  ## ---- Identify adaptive loci ----------------------------------------------

  result         <- identify_adaptive_loci(beta_draw, beta_summary, data$N_locus,
                                           data$N_env, env_names, CI)
  imp_loci_table <- result$imp_loci_table
  beta_Mlogp     <- result$beta_Mlogp

  if (is.null(imp_loci_table)) {
    inferred <- data.table::data.table(Env_var = character(),
                                       Locus   = integer())
  } else {
    inferred <- imp_loci_table[, .(Env_var, Locus)]
  }

  dir.create("files", showWarnings = FALSE, recursive = TRUE)
  dir.create("plot",  showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(inferred = inferred, CI = CI),
          "files/5_adaptive_loci.RDS")
  cat(sprintf("Inferred %d adaptive (locus, env) pairs at CI = %.2f.\n",
              nrow(inferred), CI))

  ## ---- Env color palette (Okabe-Ito CVD-safe) ------------------------------
  okabe_ito <- grDevices::palette.colors(palette = "Okabe-Ito")
  cvd_safe  <- unname(okabe_ito[!names(okabe_ito) %in% c("black", "gray")])
  if (data$N_env > length(cvd_safe))
    cvd_safe <- grDevices::palette.colors(data$N_env, palette = "Polychrome 36",
                                          recycle = FALSE)
  env_colors <- setNames(cvd_safe[seq_len(data$N_env)], env_names)

  ## ---- Manhattan plot ------------------------------------------------------

  # DEV NOTE: should convert this over to ggplot for consistancy; can pull scripting
  # from snpR for a nice plot.
  ymin <- min(beta_Mlogp, na.rm = TRUE)
  ymax <- max(c(beta_Mlogp, -log(1 - CI) + 1), na.rm = TRUE) + 1

  pdf("plot/manhattan.pdf", width = 10, height = 5)
  op <- par(mai = c(0.9, 0.9, 0.5, 0.3))

  plot(NA, xlim = c(1, data$N_locus), ylim = c(ymin, ymax),
       xaxs = "i", yaxs = "i",
       xlab = "Locus index (candidate set)",
       ylab = expression(-log(p[th])))

  ## background: every (locus, env) point as small grey dot
  for (e in seq_len(data$N_env)) {
    points(seq_len(data$N_locus), beta_Mlogp[, e],
           cex = 0.5, col = "grey70", pch = 16)
  }

  ## inferred adaptive loci as filled colored dots
  if (nrow(inferred) > 0L) {
    inferred[, Env := match(Env_var, env_names)]
    inf_y   <- beta_Mlogp[cbind(inferred$Locus, inferred$Env)]
    inf_col <- env_colors[inferred$Env_var]
    points(inferred$Locus, inf_y, pch = 21, bg = inf_col,
           col = "black", cex = 1.3)
  }

  ## threshold line
  abline(h = -log(1 - CI), lty = 1, lwd = 2, col = "grey30")

  ## env color legend
  legend("topright",
         legend = gsub("_", " ", env_names),
         pt.bg  = env_colors[env_names],
         col    = "black", pch = 21,
         cex    = 0.85, box.lty = 0, bg = NA,
         title  = "env variable")

  par(op)
  dev.off()
  cat("Saved: plot/manhattan.pdf\n")

  ## ============================================================================
  ## MIGRATION-SELECTION BALANCE for inferred adaptive loci
  ## ============================================================================

  if (nrow(inferred) > 0L) {
    cat("Compiling and exposing migration_off / migration_on ...\n")
    .calc_sel_code <- "
  functions{
    real log_kummer_1F1_pos(real a, real b, real z){
      real log_term = 0;
      real log_sum  = 0;
      for(n in 1:500){
        log_term += log((a + n - 1) * z / ((b + n - 1) * n));
        log_sum   = log_sum_exp(log_sum, log_term);
        if(log_term < log_sum - 30) break;
      }
      return log_sum;
    }
    real log_1F1_full(real a, real b, real z){
      if(z >= 0) return log_kummer_1F1_pos(a, b, z);
      return z + log_kummer_1F1_pos(b - a, b, -z);
    }
    vector pbar0_vec(real mu, real kappa, vector s){
      int N = num_elements(s);
      vector[N] out;
      real a = mu * kappa;
      real b = kappa;
      for(k in 1:N) out[k] = mu * exp(log_1F1_full(a + 1, b + 1, s[k])
                                    - log_1F1_full(a,     b,     s[k]));
      return out;
    }
    matrix migration_off(int N_patch, int N_iter,
                         vector kappa, vector rho, vector mu,
                         matrix migration, matrix s_coef){
      matrix [N_patch, N_iter] p_off;
      for(i in 1:N_iter) p_off[,i] = pbar0_vec(mu[i], kappa[i], s_coef[,i]);
      return p_off;
    }
    matrix migration_on(int N_patch, int N_iter,
                        vector kappa, vector rho, vector mu,
                        matrix migration, matrix s_coef){
      matrix [N_patch, N_iter] p_on;
      vector[N_patch] diag_M = diagonal(migration);
      vector[N_patch] S      = -diag_M;
      for(i in 1:N_iter){
        vector[N_patch] kappa_star      = kappa[i] + rho[i] * S;
        vector[N_patch] pbar0           = pbar0_vec(mu[i], kappa[i], s_coef[,i]);
        vector[N_patch] M_offdiag_pbar0 = migration * pbar0 - diag_M .* pbar0;
        vector[N_patch] mu_star = (mu[i] * kappa[i]
                                    + rho[i] * M_offdiag_pbar0) ./ kappa_star;
        for(k in 1:N_patch){
          real ak = mu_star[k] * kappa_star[k];
          real bk = kappa_star[k];
          p_on[k,i] = mu_star[k] * exp(log_1F1_full(ak + 1, bk + 1, s_coef[k,i])
                                     - log_1F1_full(ak,     bk,     s_coef[k,i]));
        }
      }
      return p_on;
    }
  }
  "
    rstan::expose_stan_functions(rstan::stanc(model_code = .calc_sel_code))

    N_locus_all <- data$N_locus
    N_iter      <- nrow(beta_draw)
    .rho   <- rho_kappa_mu[, 1:N_locus_all]
    .kappa <- rho_kappa_mu[, (N_locus_all + 1):(2 * N_locus_all)]
    .mu    <- rho_kappa_mu[, (2 * N_locus_all + 1):(3 * N_locus_all)]

    bal_table     <- as.data.frame(inferred[order(Locus)])
    bal_table$Env <- match(bal_table$Env_var, env_names)

    p_mig_off <- array(0, c(data$N_patch, N_iter, nrow(bal_table)))
    p_mig_on  <- array(0, c(data$N_patch, N_iter, nrow(bal_table)))

    cat("Computing per-deme means for each adaptive locus ...\n")
    for (i in seq_len(nrow(bal_table))) {
      l <- bal_table$Locus[i]
      e <- bal_table$Env[i]
      beta_col <- (l - 1) * data$N_env + e
      selection <- matrix(beta_draw[, beta_col],
                          nrow = data$N_patch, ncol = N_iter, byrow = TRUE) *
        matrix(as.vector(unlist(data$Env[, e])),
               nrow = data$N_patch, ncol = N_iter)
      p_mig_off[, , i] <- migration_off(data$N_patch, N_iter,
                                        .kappa[, l], .rho[, l], .mu[, l],
                                        Migration, selection)
      p_mig_on [, , i] <- migration_on (data$N_patch, N_iter,
                                        .kappa[, l], .rho[, l], .mu[, l],
                                        Migration, selection)
      if (i %% 5L == 0L || i == nrow(bal_table))
        cat(sprintf("  locus %2d / %2d done\n", i, nrow(bal_table)))
    }

    L               <- 4L
    panels_per_page <- L * L
    total_panels    <- nrow(bal_table)
    n_pages         <- ceiling(total_panels / panels_per_page)

    pdf("plot/migration_sel_balance.pdf", width = 8.5, height = 9.4)

    col_off    <- "#0072B2"
    col_on     <- "#E69F00"
    col_arrow1 <- "grey25"

    for (page in seq_len(n_pages)) {
      if (page == 1L) {
        layout(rbind(rep(1, L), matrix(2:(L * L + 1), L, L, byrow = TRUE)),
               heights = c(0.55, rep(1, L)))
        par(mar = c(0, 0, 0, 0), xpd = NA)
        plot.new()
        text(0.5, 0.80, "Migration-selection balance for adaptive loci",
             font = 2, cex = 1.55)
        legend("bottom", inset = c(0, 0.04),
               legend = c("No migration", "Mean-field  (1st correction)"),
               pch    = 21,
               pt.bg  = c(col_off, col_on),
               col    = "white",
               pt.cex = c(1.7, 1.7),
               horiz  = TRUE, cex = 1.10, bty = "n")
      } else {
        layout(matrix(seq_len(L * L), L, L, byrow = TRUE))
      }

      par(mar      = c(4.4, 4.7, 2.4, 0.7),
          mgp      = c(2.7, 0.7, 0),
          tcl      = -0.35,
          cex.axis = 1.05, cex.lab = 1.20,
          family   = "sans")

      i_start  <- (page - 1L) * panels_per_page + 1L
      i_end    <- min(page * panels_per_page, total_panels)
      page_idx <- seq(i_start, i_end)

      for (i in page_idx) {
        env_vec <- as.vector(unlist(data$Env[, bal_table$Env[i]]))
        y_off   <- rowMeans(p_mig_off[, , i])
        y_on    <- rowMeans(p_mig_on [, , i])

        yrng <- range(c(y_off, y_on))
        pad  <- max(0.025, 0.18 * diff(yrng))
        ylim <- c(max(0, yrng[1] - pad), min(1, yrng[2] + pad))

        plot(env_vec, y_on, type = "n",
             xlab = gsub("_", " ", colnames(data$Env)[bal_table$Env[i]]),
             ylab = expression(paste("Per-deme mean  ", mu[k])),
             ylim = ylim, las = 1)

        points(env_vec, y_off, pch = 21, bg = col_off, col = "white",
               cex = 1.5, lwd = 0.7)
        points(env_vec, y_on,  pch = 21, bg = col_on,  col = "white",
               cex = 1.5, lwd = 0.7)
        arrows(x0 = env_vec, y0 = y_off, x1 = env_vec, y1 = y_on,
               length = 0.07, angle = 22, code = 2, lwd = 1.3, col = col_arrow1)

        title(main = paste0("Locus ", bal_table$Locus[i]),
              line = 0.6, font.main = 2, cex.main = 1.20)
      }

      n_unused <- panels_per_page - length(page_idx)
      if (n_unused > 0L) for (k in seq_len(n_unused)) plot.new()
    }

    dev.off()
    cat("Saved: plot/migration_sel_balance.pdf\n")
  } else {
    cat("No inferred adaptive loci -- skipping migration-sel balance plot.\n")
  }

  ## ============================================================================
  ## GE-ASSOCIATION via posterior_check.stan
  ## ============================================================================

  if (nrow(inferred) == 0L) {
    cat("No inferred adaptive loci -- skipping posterior_check.stan.\n")
  } else {
    inferred      <- inferred[order(Locus)]
    imp_env_index <- match(inferred$Env_var, env_names)

    imp_loci_data <- list(
      N_patch  = data$N_patch,
      N_locus  = nrow(inferred),
      G_eff    = data$G_call[inferred$Locus, , drop = FALSE],
      TwoN_eff = data$TwoN  [inferred$Locus, , drop = FALSE],
      Env      = matrix(NA_real_, nrow = nrow(inferred), ncol = data$N_patch)
    )
    for (i in seq_len(nrow(inferred))) {
      imp_loci_data$Env[i, ] <- data$Env[, imp_env_index[i]]
    }

    STAN_MODEL <- .source_model("posterior_check.stan")
    if(!dir.exists("./stan")){dir.create("stan")}
    if(!file.exists("stan/posterior_check.stan")){
      file.copy(STAN_MODEL, "./stan/posterior_check.stan")
    }


    mod_pc <- cmdstan_model(
      "stan/posterior_check.stan",
      cpp_options = list(stan_threads = TRUE)
    )

    # DEV NOTE: is it worth allowing these to change?
    n_cores_pc <- parallel::detectCores()
    chain_pc   <- 2
    threads_pc <- max(1L, ceiling(n_cores_pc / chain_pc))

    fit_pc <- mod_pc$sample(
      data              = imp_loci_data,
      chains            = chain_pc,
      seed              = 123,
      parallel_chains   = chain_pc,
      threads_per_chain = threads_pc,
      refresh           = 1000,
      iter_warmup       = 350,
      iter_sampling     = 800,
      show_messages     = TRUE
    )

    p_mean <- matrix(fit_pc$summary(variables = "p")$mean,
                     nrow = imp_loci_data$N_locus, ncol = imp_loci_data$N_patch)
    ymin_l <- LaplacesDemon::logit(min(apply(p_mean, 1, min)))
    ymax_l <- LaplacesDemon::logit(max(apply(p_mean, 1, max)))

    L_grid <- 4L
    pdf("plot/GE_association.pdf", width = 8, height = 8)
    op <- par(mfrow = c(L_grid, L_grid), mai = c(0.55, 0.7, 0.3, 0.1))
    for (i in seq_len(nrow(inferred))) {
      plot(imp_loci_data$Env[i, ], LaplacesDemon::logit(p_mean[i, ]),
           pch = 16, col = env_colors[inferred$Env_var[i]],
           xlab = inferred$Env_var[i],
           ylab = sprintf("logit(p)  Locus %d", inferred$Locus[i]),
           ylim = c(ymin_l, ymax_l))
      fit_lm <- lm(LaplacesDemon::logit(p_mean[i, ]) ~ imp_loci_data$Env[i, ])
      abline(fit_lm, col = adjustcolor(env_colors[inferred$Env_var[i]], 0.6),
             lwd = 1.5)
    }
    par(op)
    dev.off()
    cat("Saved: plot/GE_association.pdf\n")
  }

}
