#' Likelihood-based genome-environment association test with migration pipeline.
#'
#' Runs the core likelihood-based GEA-with-migration pipeline, producing
#' output files and plots in the working directory (and returned as an object).
#' For large data (e.g. many SNPs), consider running this on small \code{Phred}
#' batches and then combining.
#'
#' @param Phred Phred-scaled genotype liklelihoods, as generated from a
#'   \code{.beagle} or \code{.beagle.gz} file by \code{\link{beagle2pl}}. See
#'   \code{\link{example_gl}}\code{$likelihoods} for an example.
#'
#' @author William Hemstrom, Nikunj Goel
#'
#' @export
#' @references Insert Goel et. al. ref here.
likeGEA <- function(Phred, patch_info, patches_xy, env,
                    N_sim_mig = 100,
                    migration_loci = 1000,
                    k_sd = 4,
                    p0 = 20,
                    seed = 42,
                    verbose = TRUE,
                    migration_model = "migration_fast",
                    max_attempts_mig = 10,
                    target_div_pct_mig = 1.0,
                    target_rhat_pct_mig = 1.0,
                    chains_mig = 4,
                    iter_sampling_mig = 500,
                    iter_warmup_mig = 250,
                    max_treedepth_mig = 12,
                    refresh_mig = 200,
                    mig_plot_threshold = 0.4,
                    mig_plot_topn = 2,
                    mig_plot_binsize = 1,
                    chains_gea = 2,
                    refresh_gea = 200,
                    iter_warmup_gea = 250,
                    iter_sampling_gea = 500,
                    plotting_CI = 0.95){

  #=====================1_prepare_data==============
  mig_data <- prepare_mig_data(Phred = Phred,
                             individuals = patch_info,
                             coordinates = patches_xy, env_data = env,
                             migration_loci = migration_loci,
                             N_sim = N_sim_mig,
                             seed = seed,
                             verbose = verbose)

  #=====================2_estimate_migration========
  Migration <- estimate_migration(data = mig_data$migration_data,
                                  model = migration_model,
                                  max_attempts = max_attempts_mig,
                                  target_div_pct = target_div_pct_mig,
                                  target_rhat_pct = target_rhat_pct_mig,
                                  chains = chains_mig,
                                  iter_sampling = iter_sampling_mig,
                                  iter_warmup = iter_warmup_mig,
                                  max_treedepth = max_treedepth_mig,
                                  refresh = refresh_mig,
                                  par = parallel::detectCores() - 2)

  #=====================3_plot_migration========
  mig_plot <- plot_migration(M = Migration,
                             data = mig_data$migration_data,
                             coordinates = patches_xy,
                             threshold_value = mig_plot_threshold,
                             topn = mig_plot_topn,
                             bin_size = mig_plot_binsize)

  #=====================4_run_stan_model========
  GEA_res <- run_stan_GEA_model(stan_data = mig_data$stan_data,
                                migration_results = Migration,
                                k_sd = k_sd,
                                p0 = p0,
                                chains = chains_gea,
                                refresh = refresh_gea,
                                iter_warmup = iter_warmup_gea,
                                iter_sampling = iter_sampling_gea,
                                seed = seed,
                                par = par)

  #=====================5_plot_results========
  plot_gea(stan_data = mig_data$stan_data,
           migration_results = Migration,
           beta_draw = GEA_res$beta_draw,
           beta_summary = GEA_res$beta_summary,
           rho_kappa_mu =  GEA_res$rho_kappa_mu,
           env_names = colnames(env),
           CI = plotting_CI)

  #=====================6_estimate_neff========
  neff <- estimate_neff(mig_data$stan_data)

  return(list(GEA_results = GEA_res,
              migration_results = Migration,
              Neff_results = neff,
              migration_plot = mig_plot,
              stan_and_migration_data = mig_data))
}

#' Prepare migration data for stan.
#'
#' Consider running via \code{\link{likeGEA}} for the full pipeline. Generally
#' for internal use.
#'
#' @export
prepare_mig_data <- function(Phred, individuals, coordinates, env_data, migration_loci,
                            N_sim, seed = 42, verbose = TRUE) {


  if(verbose){
    cat("Loaded:\n")
    cat("  Phred:       ", nrow(Phred), "loci x", ncol(Phred), "individuals\n")
    cat("  individuals: ", nrow(individuals), "rows\n")
    cat("  coordinates: ", nrow(coordinates), "patches\n")
    cat("  env_data:    ", nrow(env_data), "patches x", ncol(env_data) - 1, "env\n\n")
    cat("Running input validation ...\n")
  }


  validate_inputs(Phred, individuals, coordinates, env_data)
  if(verbose){cat("  [OK] all sanity checks passed\n\n")}


  #========FAST PHRED VALIDATION (treat malformed cells as missing)========


  if(verbose){cat("Validating Phred cell format ...\n")}
  Phred_vec  <- as.vector(Phred)
  valid_mask <- stringi::stri_detect_regex(Phred_vec, "^[0-9]+,[0-9]+,[0-9]+$")
  n_invalid  <- sum(!valid_mask)
  if (n_invalid > 0L) {
    Phred_vec[!valid_mask] <- "."
    Phred[] <- Phred_vec
  }
  n_total   <- length(Phred_vec)
  n_missing <- sum(Phred_vec == ".")
  if(verbose){cat(sprintf("  %.2f%% of cells are missing (%d / %d)\n\n",
                          100 * n_missing / n_total, n_missing, n_total))}


  #========SHARED SETUP========


  patches_in_individuals <- rle(as.character(individuals$patch))$values
  N_ind_table <- as.data.frame(table(factor(individuals$patch,
                                            levels = patches_in_individuals)))
  colnames(N_ind_table) <- c("patch", "individuals")
  N_patch <- nrow(N_ind_table)
  N_table <- as.integer(N_ind_table$individuals)
  N_env   <- ncol(env_data) - 1L
  Env     <- as.matrix(env_data[, 2:(N_env + 1L)])

  if(verbose){cat(sprintf("N_patch = %d, N_env = %d, sum(N_table) = %d (= N_ind)\n",
                          N_patch, N_env, sum(N_table)))}


  #========LOCI SUBSAMPLE FOR MIGRATION ESTIMATION========

  ## The migration model is fit on a subsampled set of loci (default 1000) to
  ## estimate pairwise diversity (pi_off, pi_diag).  The GEA model is then fit
  ## on the FULL locus set -- the migration-subsample loci are NOT excluded
  ## from the GEA scan, so every locus has a chance of being called adaptive.

  N_locus <- nrow(Phred)
  if (migration_loci > 0 && migration_loci < 1) {
    N_migration <- round(migration_loci * N_locus)
  } else if (migration_loci >= 1 && migration_loci == round(migration_loci)) {
    N_migration <- as.integer(migration_loci)
  } else {
    stop("migration_loci must be a proportion in (0, 1) or a positive integer.")
  }
  if (N_migration < 1000L) {
    N_migration <- min(1000L, N_locus)
    if(verbose){cat(sprintf("  Note: migration set adjusted to minimum of %d loci.\n",
                            N_migration))}
  }
  if (N_migration > N_locus)
    stop("migration_loci exceeds the total number of loci.")

  if(verbose){cat(sprintf("Total loci: %d\n", N_locus))}
  if (N_locus > 10000L)
    message("  Tip: 1k-10k migration loci is the recommended range.\n")
  if(verbose){cat(sprintf("  Migration set: %d loci (randomly sampled subset)\n", N_migration))}
  if(verbose){cat(sprintf("  GEA scan set : %d loci (full set)\n\n", N_locus))}

  ## Random uniform-without-replacement subsample of N_migration loci from
  ## across the entire genome.  Seeded for reproducibility.
  set.seed(seed)
  migration_idx <- sort(sample(N_locus, N_migration))

  ## Convert Phred -> log-likelihood array (full set; migration block subsets).
  if(verbose){cat("Parsing Phred -> log-likelihoods ...\n")}
  gl_all <- phred_to_loglik(Phred)
  if(verbose){cat(sprintf("  fraction NA: %.3f\n\n", mean(is.na(gl_all[, , 1]))))}


  #========PART A: STAN DATA  (G_eff, 2*N_eff for ALL loci)========


  if(verbose){cat("=== Part A: (G_eff, 2*N_eff) for all loci (GEA scan set) ===\n")}
  tic("G_eff / 2*N_eff (gl moment match)")
  eff <- compute_N_eff_G_eff_gl(gl_all, N_table = N_table)
  toc(quiet = FALSE, func.toc = toc_hours)

  G_eff    <- eff$G_eff
  TwoN_eff <- eff$TwoN_eff
  G_eff[is.na(G_eff)] <- 0   # zero-coverage cells -> zero-evidence in structured.stan

  if(verbose){cat("Computing BetaBin_theta + Group (gl_p_est.stan input for 6_estimate_Neff) ...\n")}
  tic("BetaBin_theta + Group")
  bb <- compute_BetaBin_theta_gl(gl_all, N_table)
  toc(quiet = FALSE, func.toc = toc_hours)

  stan_data <- list(
    N_patch       = N_patch,
    N_locus       = N_locus,
    N_env         = N_env,
    Env           = Env,
    G_call        = G_eff,
    TwoN          = TwoN_eff,
    BetaBin_theta = bb$BetaBin_theta,
    Group         = bb$Group,
    N_table       = N_table
  )

  dir.create("files", showWarnings = FALSE, recursive = TRUE)
  saveRDS(stan_data, "files/1_stan_data.RDS")
  cat("Saved: files/1_stan_data.RDS\n\n") # always cat this, even if verbose = FALSE.

  #========PART B: MIGRATION DATA  (pi_off, pi_diag, DistM from migration loci)========

  if(verbose){cat("=== Part B: (pi_off, pi_diag) from migration loci ===\n")}
  if(verbose){cat(sprintf("Averaging genotype resamples (N_sim = %d) ...\n", N_sim))}

  set.seed(seed + 1)
  gl_migration <- gl_all[migration_idx, , , drop = FALSE]
  pi_off_sum   <- matrix(0, N_patch, N_patch)
  pi_diag_sum  <- numeric(N_patch)
  tic("N_sim averaging")
  for (s in seq_len(N_sim)) {
    geno_s      <- genotype_from_gl(gl_migration)
    pi_s        <- compute_pi_genotype(geno_s, N_table)
    pi_off_sum  <- pi_off_sum  + pi_s$pi_off
    pi_diag_sum <- pi_diag_sum + pi_s$pi_diag
    if(verbose){if (s %% 10L == 0L) cat(sprintf("  %d / %d\n", s, N_sim))}
  }
  toc(quiet = FALSE, func.toc = toc_hours)
  pi_off  <- pi_off_sum  / N_sim
  pi_diag <- pi_diag_sum / N_sim
  if(verbose){cat(sprintf("  pi_diag range          : [%.4f, %.4f]\n",
                          min(pi_diag, na.rm = TRUE), max(pi_diag, na.rm = TRUE)))}
  if(verbose){cat(sprintf("  pi_off (upper.tri) rng : [%.4f, %.4f]\n",
                          min(pi_off[upper.tri(pi_off)], na.rm = TRUE),
                          max(pi_off[upper.tri(pi_off)], na.rm = TRUE)))}

  ## Geographic distance matrix from patch coordinates (geodesic; real data uses
  ## lat/lon, so this is the right measure -- not Euclidean).
  DistM <- geodist::geodist(as.data.frame(coordinates), measure = "geodesic")

  migration_data <- list(
    N_patch         = N_patch,
    N_locus_neutral = N_migration,
    N_table         = N_table,
    DistM           = DistM,
    pi_off          = pi_off,
    pi_diag         = pi_diag
  )
  saveRDS(migration_data, "files/1_migration_data.RDS")
  cat("Saved: files/1_migration_data.RDS\n\n")

  toc(quiet = FALSE, func.toc = toc_hours)

  return(list(migration_data = migration_data, stan_data = stan_data))
}

#' Estimate migration rates.
#'
#' Consider running via \code{\link{likeGEA}} for the full pipeline. Generally
#' for internal use or when partial re-runs are needed.
#'
#' @param model character, default "migration_fast". Choice of model for
#'   estimating migration. Options are "migration_fast" or "migration_M3b".
#'   See details.
#' @export
estimate_migration <- function(data,
                               model = "migration_fast",
                               max_attempts = 10,
                               target_div_pct = 1.0,
                               target_rhat_pct = 1.0,
                               chains = 4,
                               iter_sampling = 500,
                               iter_warmup = 250,
                               max_treedepth = 12,
                               refresh = 200,
                               par = parallel::detectCores() - 2){
  tic("migration estimation")

  ## ---- Stan model -- swap by uncommenting (both consume the same data) ----
  ##   migration_fast.stan -- Petkova-centered combination (default)
  ##   migration_M3b.stan  -- joint (pi_off, pi_diag) fit with diagonal anchor
  model <- paste0(model, ".stan")
  STAN_MODEL <- .source_model(model)
  if(!dir.exists("./stan")){dir.create("stan")}
  file.copy(STAN_MODEL, "./stan")
  STAN_MODEL <- file.path("stan", model)
  # STAN_MODEL <- "stan/migration_M3b.stan"

  ## ---- Load data ------------------------------------------------------------

  N_patch <- data$N_patch
  N_pair  <- N_patch * (N_patch - 1L) %/% 2L

  stopifnot(is.matrix(data$pi_off),
            all(dim(data$pi_off) == c(N_patch, N_patch)),
            length(data$pi_diag) == N_patch,
            all(!is.na(data$pi_off[upper.tri(data$pi_off)])),
            all(!is.na(data$pi_diag)))

  ## Standardise DistM by SD of off-diagonal entries.
  DistM_std  <- data$DistM / sd(data$DistM[upper.tri(data$DistM)])
  covariates <- matrix(convert_T_to_t(DistM_std), ncol = 1)

  ## ---- Compile Stan model ---------------------------------------------------
  hpp <- .source_model("coalescent_time.hpp")
  if(!file.exists(file.path("stan", "coalescent_time.hpp"))){
    file.copy(hpp, "./stan")
  } else {
    cat("Using existing stan/coalescent_time.hpp\n")
  }
  hpp <- file.path("stan", hpp)
  hpp_file <- normalizePath(hpp)

  mod <- cmdstan_model(
    STAN_MODEL,
    cpp_options   = list(USER_HEADER = hpp_file, stan_threads = TRUE),
    stanc_options = list("allow-undefined" = TRUE, "O1")
  )

  ## ---- Sample ---------------------------------------------------------------

  n_cores           <- par
  threads_per_chain <- max(1L, floor(n_cores / chains))

  dat <- list(
    pi_off     = convert_T_to_t(data$pi_off),    # length N_pair
    pi_diag    = data$pi_diag,                   # length N_patch
    covariates = covariates,
    N_cov      = 1L,
    N_patch    = N_patch,
    gamma_sd   = 1
  )

  ## ---- Retry loop: resample with fresh seed until divergences < target ----
  total_iter    <- chains * iter_sampling

  best_fit       <- NULL
  best_score     <- Inf
  best_div_pct   <- NA_real_
  best_rhat_pct  <- NA_real_
  best_seed      <- NA_integer_

  model_label <- sub("\\.stan$", "", basename(STAN_MODEL))
  for (attempt in seq_len(max_attempts)) {
    seed <- as.integer(Sys.time()) + attempt
    cat(sprintf("\n=== %s attempt %d/%d  (seed=%d) ===\n",
                model_label, attempt, max_attempts, seed))

    fit <- mod$sample(
      data              = dat,
      chains            = chains,
      parallel_chains   = chains,
      seed              = seed,
      iter_warmup       = iter_warmup,
      iter_sampling     = iter_sampling,
      threads_per_chain = threads_per_chain,
      max_treedepth     = max_treedepth,
      refresh           = refresh,
      show_messages     = TRUE,
      show_exceptions   = FALSE
    )

    n_div    <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
    div_pct  <- 100 * n_div / total_iter

    rhat_vec   <- fit$summary()$rhat
    rhat_vec   <- rhat_vec[is.finite(rhat_vec)]
    rhat_pct   <- 100 * mean(rhat_vec > 1.05)
    n_bad_rhat <- sum(rhat_vec > 1.05)

    cat(sprintf("  -> div: %d (%.2f%% of %d), target < %.1f%%\n",
                n_div, div_pct, total_iter, target_div_pct))
    cat(sprintf("  -> Rhat>1.05: %d / %d params (%.2f%%), target < %.1f%%\n",
                n_bad_rhat, length(rhat_vec), rhat_pct, target_rhat_pct))

    score <- max(div_pct, rhat_pct)
    if (score < best_score) {
      best_fit      <- fit
      best_score    <- score
      best_div_pct  <- div_pct
      best_rhat_pct <- rhat_pct
      best_seed     <- seed
    }

    if (div_pct < target_div_pct && rhat_pct < target_rhat_pct) {
      cat(sprintf("  ACCEPTED after %d attempt(s); seed=%d, div=%.2f%%, Rhat>1.05=%.2f%%.\n",
                  attempt, seed, div_pct, rhat_pct))
      break
    }
  }

  fit <- best_fit
  if (best_div_pct >= target_div_pct || best_rhat_pct >= target_rhat_pct) {
    warning(sprintf("\nWARNING: %d attempts exhausted without hitting targets (<%.1f%% div, <%.1f%% Rhat>1.05).\n",
                    max_attempts, target_div_pct, target_rhat_pct))
    cat(sprintf("         Using best run: seed=%d, div=%.2f%%, Rhat>1.05=%.2f%%.\n",
                best_seed, best_div_pct, best_rhat_pct))
  }

  ## ---- Extract + save migration matrix --------------------------------------

  M_summary <- fit$summary(variables = "Migration")
  Migration <- matrix(M_summary$mean, nrow = N_patch, ncol = N_patch)
  saveRDS(Migration, "files/2_migration.RDS")
  cat("Saved: files/2_migration.RDS\n")

  ## ---- Diagnostics ----------------------------------------------------------
  all_vars    <- fit$metadata()$variables
  scalar_pars <- setdiff(all_vars[!grepl("\\[", all_vars) |
                                    grepl("^gamma\\[", all_vars)],
                         c("lp__", "Migration"))
  print(fit$summary(variables = scalar_pars))
  cat(sprintf("Migration[i,j]:  max Rhat = %.3f,  min ESS_bulk = %.0f\n",
              max(M_summary$rhat,     na.rm = TRUE),
              min(M_summary$ess_bulk, na.rm = TRUE)))

  toc(quiet = FALSE, func.toc = toc_hours)

  return(Migration)
}

#' Run the core STAN model to estimate genotype-environment association with
#' likelihoods.
#'
#' Consider running via \code{\link{likeGEA}} for the full pipeline. Generally
#' for internal use or when partial re-runs are needed.
#'
#' @export
run_stan_GEA_model <- function(stan_data,
                               migration_results,
                               k_sd = 4,
                               p0 = 20,
                               refresh = 100,
                               chains = 2,
                               iter_warmup = 250,
                               iter_sampling = 500,
                               seed = 42,
                               par = parallel::detectCores() - 2){
  tic("estimating Stan time")

  ##   migration_fast.stan -- Petkova-centered combination (default)
  ##   migration_M3b.stan  -- joint (pi_off, pi_diag) fit with diagonal anchor
  STAN_MODEL <- .source_model("stan/structured.stan")
  if(!dir.exists("./stan")){dir.create("stan")}

  if(!file.exists(file.path("stan", "structured.stan"))){
    file.copy(STAN_MODEL, "./stan")
  }

  STAN_MODEL <- file.path("stan", model)
  ## ---- Load data ------------------------------------------------------------

  data      <- stan_data; rm(stan_data)
  Migration <- migration_results; rm(migration_results)

  ## ---- Append hyperparameters -----------------------------------------------

  data$k_sd <- k_sd
  data$p0   <- p0

  ## structured.stan declares G_call and TwoN as real (lbeta + 1F1 arithmetic).
  data$G_call <- as.double(data$G_call)
  data$TwoN   <- as.double(data$TwoN)
  dim(data$G_call) <- c(data$N_locus, data$N_patch)
  dim(data$TwoN)   <- c(data$N_locus, data$N_patch)

  ## ---- Write data include file ----------------------------------------------

  write_common_data_stan(
    as.matrix(data$Env),
    as.matrix(Migration),
    data$k_sd,
    "stan/structured_data.stan"
  )

  ## ---- Compile Stan model ---------------------------------------------------
  hpp <- .source_model("1F1_fast.hpp")
  nhpp <- file.path("stan", "1F1_fast.hpp")
  if(!file.exists(nhpp)){
    file.copy(hpp, nhpp)
  } else {
    cat("Using existing stan/1F1_fast.hpp\n")
  }
  hpp <- nhpp
  hpp_file <- normalizePath(hpp)

  mod <- cmdstan_model(
    STAN_MODEL,
    cpp_options     = list(USER_HEADER = hpp_file, stan_threads = TRUE),
    stanc_options   = list("allow-undefined" = TRUE),
    include_paths   = normalizePath("stan"),
    force_recompile = TRUE
  )

  ## ---- Sample ---------------------------------------------------------------

  n_cores           <- par
  threads_per_chain <- floor(par/chains)
  threads_per_chain <- max(c(threads_per_chain, 1))


  # DEV NOTE: should modify the function to take these as args at some point
  init_fun <- function() {
    list(
      mu          = runif(data$N_locus, 0.4, 0.6),
      kappa       = exp(rnorm(data$N_locus, log(1), 0.3)),
      rho_tilda   = exp(rnorm(data$N_locus, log(1), 0.2)),
      tau_unif    = runif(1, base::pi / 10, base::pi / 4),
      c_square    = abs(rnorm(1, 1.3, 0.3)),
      z_beta      = matrix(rnorm(data$N_env * data$N_locus, 0, 0.1),
                           nrow = data$N_env, ncol = data$N_locus),
      lambda_unif = matrix(runif(data$N_env * data$N_locus,
                                 base::pi / 10, base::pi / 5),
                           nrow = data$N_env, ncol = data$N_locus)
    )
  }

  fit <- mod$sample(
    data              = data,
    chains            = chains,
    seed              = seed,
    parallel_chains   = chains,
    init              = init_fun,
    show_exceptions   = TRUE,
    show_messages     = TRUE,
    refresh           = refresh,
    iter_warmup       = iter_warmup,
    threads_per_chain = threads_per_chain,
    iter_sampling     = iter_sampling
  )

  ## ---- Save posterior draws -------------------------------------------------

  dir.create("files", showWarnings = FALSE, recursive = TRUE)

  rho_kappa_mu <- fit$draws(variables = c("rho", "kappa", "mu"),
                            format = "draws_matrix")
  saveRDS(rho_kappa_mu, "files/4_rho_kappa_mu_draw.RDS")

  beta_draw <- fit$draws(variables = c("Beta"), format = "draws_matrix")
  saveRDS(beta_draw, "files/4_beta_draw.RDS")

  tau_draw <- fit$draws(variables = c("Tau"), format = "draws_matrix")
  saveRDS(tau_draw, "files/4_tau_draw.RDS")

  beta_summary <- summarise_draws_parallel(beta_draw)
  saveRDS(beta_summary, "files/4_beta_summary.RDS")

  mlog_p <- score_stan(beta_draw, data$N_locus, data$N_env)
  saveRDS(mlog_p, "files/4_mlog_p.RDS")

  cat("Saved files/4_{rho_kappa_mu_draw, beta_draw, tau_draw,",
      "beta_summary, mlog_p}.RDS\n")

  toc(quiet = FALSE, func.toc = toc_hours)

  return(list(rho_kappa_mu = rho_kappa_mu,
              beta_draw = beta_draw,
              tau = tau_draw,
              beta_summary = beta_summary,
              mlog_p = mlog_p))

}

#' Estimate effective sample sizes.
#'
#' Consider running via \code{\link{likeGEA}} for the full pipeline. Generally
#' for internal use or when partial re-runs are needed.
#'
#' @export
estimate_neff <- function(stan_data){
  tic("estimating Neff time")

  ## ---- Load -----------------------------------------------------------------
  N_patch      <- stan_data$N_patch
  N_locus_tot  <- stan_data$N_locus
  G_eff_ana    <- stan_data$G_call         # analytical
  TwoN_eff_ana <- stan_data$TwoN           # analytical

  ## Subsample loci for speed (gl_p_est sampling is the cost driver).
  N_loci_use <- min(1000L, N_locus_tot)
  if (N_loci_use < N_locus_tot) {
    set.seed(789)
    loci_used <- sort(sample(N_locus_tot, N_loci_use))
    cat(sprintf("Subsampled %d of %d loci for N_eff estimation.\n",
                N_loci_use, N_locus_tot))
  } else {
    loci_used <- seq_len(N_locus_tot)
  }

  dat <- list(
    N_patch       = N_patch,
    N_locus       = length(loci_used),
    N_table       = as.integer(stan_data$N_table),
    BetaBin_theta = stan_data$BetaBin_theta[loci_used, , drop = FALSE],
    Group         = stan_data$Group        [loci_used, , drop = FALSE]
  )

  ## ---- Compile + sample gl_p_est.stan ---------------------------------------

  STAN_MODEL <- .source_model("gl_p_est.stan")
  if(!dir.exists("./stan")){dir.create("stan")}
  if(!file.exists("stan/gl_p_est.stan")){
    file.copy(STAN_MODEL, "./stan/gl_p_est.stan")
  }

  mod <- cmdstan_model(
    "stan/gl_p_est.stan",
    cpp_options   = list(stan_threads = TRUE),
    stanc_options = list("O1")
  )

  n_cores           <- parallel::detectCores() - 2
  chain             <- 2
  threads_per_chain <- max(1L, floor(n_cores / chain))

  # DEV NOTE: Is it worth it to allow all of these to change?
  tic("NUTS posterior on p")
  fit <- mod$sample(
    data              = dat,
    chains            = chain,
    parallel_chains   = chain,
    seed              = 546,
    threads_per_chain = threads_per_chain,
    refresh           = 100,
    iter_warmup       = 200,
    iter_sampling     = 500,
    show_messages     = TRUE
  )
  toc(quiet = FALSE, func.toc = toc_hours)

  ## ---- Extract per-(locus, patch) posterior moments -------------------------

  p_summary <- summarise_draws_parallel(fit$draws("p", format = "draws_matrix"))
  idx_mat   <- do.call(rbind, regmatches(
    p_summary$variable,
    regexec("p\\[([0-9]+),([0-9]+)\\]", p_summary$variable)))
  i_idx     <- as.integer(idx_mat[, 2])
  k_idx     <- as.integer(idx_mat[, 3])
  stopifnot(length(i_idx) == length(loci_used) * N_patch,
            all(!is.na(i_idx)), all(!is.na(k_idx)))

  mu_post  <- matrix(NA_real_, length(loci_used), N_patch)
  sd_post  <- matrix(NA_real_, length(loci_used), N_patch)
  mu_post[cbind(i_idx, k_idx)] <- p_summary$mean
  sd_post[cbind(i_idx, k_idx)] <- p_summary$sd
  var_post <- sd_post^2

  ## ---- Beta(G+1, 2N-G+1) moment-match  (same convention as analytical) ------

  TwoN_eff_stan <- pmax(mu_post * (1 - mu_post) / pmax(var_post, 1e-20) - 3, 1e-6)
  G_eff_stan    <- (2 + TwoN_eff_stan) * mu_post - 1
  G_eff_stan    <- pmax(pmin(G_eff_stan, TwoN_eff_stan), 0)

  ## ---- Save -----------------------------------------------------------------



  ## Mean per-individual N_eff diagnostic (Stan estimate, like the original).
  Neff_mean <- mean(rowSums(TwoN_eff_stan) / 2) / sum(stan_data$N_table)
  cat(sprintf("Mean N_eff / N (Stan, across subsampled loci): %.4f\n", Neff_mean))

  ret <- list(
    loci_used     = loci_used,
    G_eff_stan    = G_eff_stan,
    TwoN_eff_stan = TwoN_eff_stan,
    G_eff_ana     = G_eff_ana   [loci_used, , drop = FALSE],
    TwoN_eff_ana  = TwoN_eff_ana[loci_used, , drop = FALSE],
    Neff_mean     = Neff_mean
  )

  dir.create("files", showWarnings = FALSE, recursive = TRUE)
  saveRDS(ret, "files/6_Neff.RDS")
  cat("Saved: files/6_Neff.RDS\n")

  toc(quiet = FALSE, func.toc = toc_hours)

  return(ret)
}
