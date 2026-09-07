###############################################################################
## 0_functions.R
##
## Self-contained helpers + package manager for the analyze_your_data
## pipeline.  Sourced by every numbered script:
##   source("0_functions.R")
## (Working directory is the analyze_your_data/ folder.)
##
## This file is intentionally self-contained — it does NOT source anything
## outside the folder.  The R helpers below are direct ports of the
## versions in code/synthetic_simulations/0_functions.R (snapshotted at
## refactor time); update them in lockstep if the upstream versions evolve.
##
## Helpers provided:
##   validate_inputs           -- sanity-check user-supplied data
##   phred_to_loglik           -- Phred string matrix -> log-likelihood array
##   genotype_from_gl          -- NA-safe genotype resample from GL
##   compute_N_eff_G_eff_gl    -- log-seq Gaussian moment match (parallel)
##   compute_BetaBin_theta_gl  -- polynomial log-weights (gl_p_est.stan input)
##   compute_pi_genotype       -- pi_off / pi_diag from genotypes
##   convert_T_to_t            -- square distance matrix -> upper-tri vector
##   write_common_data_stan    -- eigendecomp + Env -> structured_data.stan
##   summarise_draws_parallel  -- PSOCK-parallel posterior::summarise_draws
##   identify_adaptive_loci    -- CI-threshold -> adaptive loci table
##   score_stan                -- [N_locus x N_env] mlog_p matrix
##   toc_hours                 -- adaptive-unit elapsed-time formatter
###############################################################################

validate_inputs <- function(phred, individuals, coordinates, env_data) {
  ## Sanity-checks for the four user-supplied input files.
  N_locus <- nrow(phred)
  N_ind   <- ncol(phred)
  if (nrow(individuals) != N_ind)
    stop(sprintf("individuals has %d rows but phred has %d columns",
                 nrow(individuals), N_ind))
  if (!all(c("ind_id", "patch") %in% names(individuals)))
    stop("individuals must have columns: ind_id, patch")
  if (!all(c("patch", "lat", "lon") %in% names(coordinates)))
    stop("coordinates must have columns: patch, lat, lon")
  if (!"patch" %in% names(env_data))
    stop("env_data must have a patch column")

  ## Block-ordering: all rows for a given patch must be consecutive.
  patch_runs <- rle(as.character(individuals$patch))
  if (any(duplicated(patch_runs$values)))
    stop("individuals is not block-ordered by patch ",
         "(patches appear in multiple non-consecutive runs)")

  ## Patches in coordinates and env_data must match individuals' patches.
  ind_patches <- unique(individuals$patch)
  if (!setequal(ind_patches, coordinates$patch))
    stop("patches in individuals do not match patches in coordinates")
  if (!setequal(ind_patches, env_data$patch))
    stop("patches in individuals do not match patches in env_data")

  ## Phred cells: "Q1,Q2,Q3" or ".".
  bad <- !grepl("^([0-9]+,[0-9]+,[0-9]+|\\.)$", phred[1:min(100L, length(phred))])
  if (any(bad))
    stop("phred matrix contains malformed cells (expected 'Q1,Q2,Q3' or '.')")

  invisible(TRUE)
}


## ============================================================================
## PHRED -> LOG-LIKELIHOOD CODEC
## ============================================================================

phred_to_loglik <- function(phred_mat) {
  ## phred_mat : CHARACTER matrix [N_locus x N_ind]; cells "Q_aa,Q_Aa,Q_AA" or "."
  ## Returns   : [N_locus x N_ind x 3] log-likelihood array, row-normalised so
  ##             the three entries sum to 1 in linear space.  Cells that were
  ##             "." get NA on all three entries.
  N_locus <- nrow(phred_mat)
  N_ind   <- ncol(phred_mat)
  v <- as.vector(phred_mat)
  na_mask <- (v == ".")
  v[na_mask] <- "NA,NA,NA"

  triples <- stringi::stri_split_fixed(v, ",", simplify = TRUE)
  storage.mode(triples) <- "double"

  log_raw <- -0.1 * triples * log(10)
  m <- matrixStats::rowMaxs(log_raw)
  log_norm <- m + log(rowSums(exp(log_raw - m)))
  ll <- log_raw - log_norm
  ll[na_mask, ] <- NA_real_

  array(c(ll[, 1], ll[, 2], ll[, 3]), dim = c(N_locus, N_ind, 3))
}


## ============================================================================
## GL GENOTYPE RESAMPLING (NA-safe)
## ============================================================================

genotype_from_gl <- function(gl_log) {
  ## gl_log : [N_locus x N_ind x 3] log-likelihood array (NA at empty depth)
  ## Returns: [N_locus x N_ind] integer matrix (0/1/2), NA where input was NA.
  A <- exp(gl_log)
  d <- dim(A)
  m <- matrix(A, ncol = d[3])
  has_data <- rowSums(is.na(m)) == 0
  res <- rep(NA_integer_, nrow(m))
  if (any(has_data)) {
    r <- runif(sum(has_data), 0, rowSums(m[has_data, , drop = FALSE]))
    sub <- m[has_data, , drop = FALSE]
    res[has_data] <- ifelse(r < sub[, 1], 0L,
                     ifelse(r < sub[, 1] + sub[, 2], 1L, 2L))
  }
  matrix(res, d[1], d[2])
}


## ============================================================================
## EFFECTIVE PAIR (G_eff, 2*N_eff) -- log-seq Gaussian moment match
## ============================================================================

.log_seq_moments_one <- function(L) {
  N <- nrow(L)
  if (N == 0L) return(c(mu = NA_real_, sigma2 = NA_real_))
  L_safe <- pmax(L, 1e-300)
  b0 <- log(L_safe[, 1])
  b1 <- log(2 * L_safe[, 2])
  b2 <- log(L_safe[, 3])
  log_A <- c(b0[1], b1[1], b2[1])
  if (N >= 2) for (i in 2:N) {
    cur <- length(log_A); new <- cur + 2L
    x0 <- rep(-Inf, new); x0[1:cur]       <- log_A + b0[i]
    x1 <- rep(-Inf, new); x1[2:(cur + 1)] <- log_A + b1[i]
    x2 <- rep(-Inf, new); x2[3:(cur + 2)] <- log_A + b2[i]
    m <- pmax(x0, x1, x2)
    log_A <- m + log(exp(x0 - m) + exp(x1 - m) + exp(x2 - m))
  }
  G_idx    <- 0:(2 * N)
  log_pi   <- log_A - lchoose(2 * N, G_idx)
  pi_G     <- exp(log_pi - max(log_pi)); pi_G <- pi_G / sum(pi_G)
  EG       <- sum(pi_G * G_idx)
  EGGp1Gp2 <- sum(pi_G * (G_idx + 1) * (G_idx + 2))
  mu_      <- (EG + 1) / (2 * N + 2)
  sigma2   <- EGGp1Gp2 / ((2 * N + 2) * (2 * N + 3)) - mu_^2
  c(mu = mu_, sigma2 = sigma2)
}

compute_N_eff_G_eff_gl <- function(gl, N_table, n_cores = NULL) {
  ## gl      : [N_locus x N_ind x 3] log-likelihood array (NA at empty depth)
  ## N_table : integer vector of individuals per patch (block-ordered)
  ## Returns : list(G_eff, TwoN_eff, N_sample), each [N_locus x N_patch].
  if (is.null(n_cores)) n_cores <- max(1, min(parallel::detectCores(), 12))
  N_locus <- dim(gl)[1]
  N_patch <- length(N_table)
  patch_start <- cumsum(c(1L, head(N_table, -1)))
  patch_end   <- cumsum(N_table)
  patch_idx <- lapply(seq_len(N_patch), function(k) patch_start[k]:patch_end[k])

  per_locus <- function(l) {
    G  <- numeric(N_patch)
    Tn <- numeric(N_patch)
    Ns <- integer(N_patch)
    for (k in seq_len(N_patch)) {
      L <- gl[l, patch_idx[[k]], , drop = FALSE]
      L <- matrix(L, ncol = 3)
      keep <- !rowSums(is.na(L)) > 0
      L <- L[keep, , drop = FALSE]
      Ns[k] <- as.integer(N_table[k])
      if (nrow(L) == 0L) {
        G[k] <- NA_real_; Tn[k] <- 0
        next
      }
      L_lin <- exp(L)
      L_lin <- L_lin / rowSums(L_lin)
      mom <- .log_seq_moments_one(L_lin)
      mu_ <- mom["mu"]; s2 <- mom["sigma2"]
      twoN <- max(mu_ * (1 - mu_) / max(s2, 1e-20) - 3, 1e-6)
      g_eff <- (2 + twoN) * mu_ - 1
      g_eff <- max(min(g_eff, twoN), 0)
      G[k] <- g_eff
      Tn[k] <- twoN
    }
    list(G = G, Tn = Tn, Ns = Ns)
  }

  chunks <- parallel::splitIndices(N_locus, n_cores)
  parts <- parallel::mclapply(chunks, function(idx) {
    out <- vector("list", length(idx))
    for (k in seq_along(idx)) out[[k]] <- per_locus(idx[k])
    out
  }, mc.cores = n_cores)

  G_eff    <- matrix(NA_real_, N_locus, N_patch)
  TwoN_eff <- matrix(NA_real_, N_locus, N_patch)
  N_sample <- matrix(0L,       N_locus, N_patch)
  pos <- 1L
  for (j in seq_along(parts)) for (item in parts[[j]]) {
    G_eff   [pos, ] <- item$G
    TwoN_eff[pos, ] <- item$Tn
    N_sample[pos, ] <- item$Ns
    pos <- pos + 1L
  }

  list(G_eff = G_eff, TwoN_eff = TwoN_eff, N_sample = N_sample)
}


## ============================================================================
## GL-ONLY BetaBin_theta POLYNOMIAL BUILDER  (input to gl_p_est.stan)
## ============================================================================

.compute_log_coeffs <- function(loga, logb, logc) {
  stopifnot(length(loga) == length(logb), length(logb) == length(logc))
  lse3 <- function(x, y, z) {
    m <- pmax(x, y, z)
    m + log(exp(x - m) + exp(y - m) + exp(z - m))
  }
  logcfs <- 0
  for (i in seq_along(loga)) {
    s0 <- c(logcfs + loga[i], -Inf, -Inf)
    s1 <- c(-Inf, logcfs + logb[i], -Inf)
    s2 <- c(-Inf, -Inf, logcfs + logc[i])
    logcfs <- mapply(lse3, s0, s1, s2)
  }
  logcfs
}

compute_BetaBin_theta_gl <- function(gl, N_table, n_cores = NULL) {
  ## gl      : [N_locus x N_ind x 3] log-likelihood array (NA at empty depth)
  ## N_table : integer vector of individuals per patch
  ## Returns : list(BetaBin_theta = [N_locus x (2*sum(N_table)+N_patch)],
  ##               Group         = [N_locus x N_patch] non-NA count)
  if (is.null(n_cores)) n_cores <- max(1, min(parallel::detectCores(), 12))
  N_locus <- dim(gl)[1]
  N_patch <- length(N_table)
  patch_start <- cumsum(c(1L, head(N_table, -1)))
  patch_end   <- cumsum(N_table)
  patch_idx   <- lapply(seq_len(N_patch), function(k) patch_start[k]:patch_end[k])

  Group <- matrix(0L, N_locus, N_patch)
  for (k in seq_len(N_patch)) {
    sub <- gl[, patch_idx[[k]], 1]
    if (length(patch_idx[[k]]) == 1L) sub <- matrix(sub, ncol = 1)
    Group[, k] <- rowSums(!is.na(sub))
  }

  tot_cols   <- sum(2L * N_table + 1L)
  lens       <- 2L * Group + 1L
  pos_offset <- t(apply(lens, 1L, function(x) cumsum(c(1L, x[-length(x)]))))

  per_locus <- function(l) {
    row_l <- numeric(tot_cols)
    for (j in seq_len(N_patch)) {
      if (Group[l, j] > 0L) {
        start_col <- pos_offset[l, j]
        end_col   <- start_col + 2L * Group[l, j]
        g1 <- na.omit(gl[l, patch_idx[[j]], 1])
        g2 <- na.omit(gl[l, patch_idx[[j]], 2]) + log(2)
        g3 <- na.omit(gl[l, patch_idx[[j]], 3])
        row_l[start_col:end_col] <- .compute_log_coeffs(g1, g2, g3)
      }
    }
    row_l
  }

  chunks <- parallel::splitIndices(N_locus, n_cores)
  parts <- parallel::mclapply(chunks, function(idx) {
    out <- vector("list", length(idx))
    for (k in seq_along(idx)) out[[k]] <- per_locus(idx[k])
    out
  }, mc.cores = n_cores)

  BetaBin_theta <- matrix(0, N_locus, tot_cols)
  pos <- 1L
  for (j in seq_along(parts)) for (row_l in parts[[j]]) {
    BetaBin_theta[pos, ] <- row_l
    pos <- pos + 1L
  }

  list(BetaBin_theta = BetaBin_theta, Group = Group)
}


## ============================================================================
## PAIRWISE DIVERSITIES (pi_off, pi_diag)
## ============================================================================

compute_pi_genotype <- function(genotype, N_table) {
  ## genotype : [N_locus x N_ind] integer 0/1/2/NA, block-ordered by patch
  ## N_table  : integer vector of individuals per patch
  ## Returns  : list(pi_off  = [N_patch x N_patch] symmetric (diag 0),
  ##                 pi_diag = length N_patch)
  N_locus <- nrow(genotype)
  N_patch <- length(N_table)
  patch_start <- cumsum(c(1L, head(N_table, -1)))
  patch_end   <- cumsum(N_table)

  p_hat <- matrix(NA_real_, N_locus, N_patch)
  for (k in seq_len(N_patch)) {
    gs <- genotype[, patch_start[k]:patch_end[k], drop = FALSE]
    p_hat[, k] <- rowMeans(gs, na.rm = TRUE) / 2
  }

  pi_diag <- vapply(seq_len(N_patch), function(k)
    mean(2 * p_hat[, k] * (1 - p_hat[, k]), na.rm = TRUE),
    numeric(1))

  pi_off <- matrix(0, N_patch, N_patch)
  for (i in seq_len(N_patch - 1L)) for (j in (i + 1L):N_patch) {
    pi_off[i, j] <- mean(p_hat[, i] * (1 - p_hat[, j]) +
                         (1 - p_hat[, i]) * p_hat[, j], na.rm = TRUE)
    pi_off[j, i] <- pi_off[i, j]
  }

  list(pi_off = pi_off, pi_diag = pi_diag)
}


## ============================================================================
## MIGRATION-FIT HELPERS
## ============================================================================

convert_T_to_t <- function(A) {
  ## A : N_patch x N_patch symmetric matrix.
  ## Returns the strict-upper-triangle in row-major order as a length-
  ## (N_patch choose 2) vector.
  N_patch <- nrow(A)
  n <- choose(N_patch, 2)
  t <- rep(0, n)
  pos <- 1
  for (i in 1:(N_patch - 1)) {
    for (j in (i + 1):N_patch) {
      t[pos] <- A[i, j]
      pos <- pos + 1
    }
  }
  t
}


## ============================================================================
## STRUCTURED-FIT HELPERS
## ============================================================================

write_common_data_stan <- function(env_matrix, mig_matrix, k_sd, filename) {
  ## Writes the data-include file used by structured.stan
  ## (`structured_data.stan`).  Contains the eigendecomp of mig_matrix
  ## together with env_matrix and k_sd.
  N_patch <- nrow(env_matrix)
  N_env   <- ncol(env_matrix)
  eig <- eigen(mig_matrix, symmetric = TRUE)

  eigenval_str <- glue::glue_collapse(format(eig$values, digits = 15), sep = ", ")
  row_strings_Q <- apply(eig$vectors, 1, function(row)
    glue::glue("[{glue::glue_collapse(format(row, digits = 15), sep = ', ')}]"))
  Q_str <- glue::glue_collapse(row_strings_Q, sep = ',\n  ')
  row_strings_Env <- apply(env_matrix, 1, function(row)
    glue::glue("[{glue::glue_collapse(format(row, digits = 15), sep = ', ')}]"))
  Env_str <- glue::glue_collapse(row_strings_Env, sep = ',\n  ')

  out_text <- glue::glue(
    "real k_sd = {k_sd};\n",
    "int N_env = {N_env};\n",
    "int N_patch = {N_patch};\n",
    "vector[N_patch] d = [{eigenval_str}]';\n",
    "matrix[N_patch, N_patch] Q = [{Q_str}];\n",
    "matrix[N_patch, N_env] Env = [{Env_str}];\n"
  )
  writeLines(out_text, filename)
  out_text
}

summarise_draws_parallel <- function(draws, n_cores = NULL,
                                       chunk_min = 50L) {
  ## SAFE PARALLEL wrapper over posterior::summarise_draws via PSOCK
  ## (avoids macOS mclapply fork hangs after a cmdstanr fit).
  if (!inherits(draws, "draws_matrix"))
    draws <- posterior::as_draws_matrix(draws)

  sequential_summary <- function(d) {
    res <- posterior::summarise_draws(
      d,
      mean, median, sd, mad,
      ~posterior::quantile2(.x, probs = c(0.05, 0.95)),
      posterior::rhat, posterior::ess_bulk, posterior::ess_tail
    )
    data.table::setDT(res)
    res
  }

  n_params <- ncol(draws)
  if (is.null(n_cores)) n_cores <- min(parallel::detectCores(), 8L)
  if (n_cores <= 1L || n_params < chunk_min * 2L)
    return(sequential_summary(draws))

  tryCatch({
    chunks       <- parallel::splitIndices(n_params, n_cores)
    draws_chunks <- lapply(chunks, function(idx)
      posterior::as_draws_matrix(draws[, idx, drop = FALSE]))

    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(posterior))
    })

    parts <- parallel::parLapply(cl, draws_chunks, function(d) {
      posterior::summarise_draws(
        d,
        mean, median, sd, mad,
        ~posterior::quantile2(.x, probs = c(0.05, 0.95)),
        posterior::rhat, posterior::ess_bulk, posterior::ess_tail
      )
    })

    result <- do.call(rbind, parts)
    data.table::setDT(result)
    result
  }, error = function(e) {
    message("summarise_draws_parallel: PSOCK failed (",
            conditionMessage(e), "); falling back to sequential")
    sequential_summary(draws)
  })
}

.split_string <- function(x) strsplit(x, ",|\\]")[[1]][2]

identify_adaptive_loci <- function(beta_draw, beta_summary, N_locus, N_env,
                                   env_names, CI = 0.95) {
  ## Returns list(imp_loci_table, beta_Mlogp, Q95, Q5).
  n_draws <- nrow(beta_draw)
  Q95     <- pmax(colMeans(beta_draw > 0), 1 / n_draws)
  Q5      <- pmax(colMeans(beta_draw < 0), 1 / n_draws)
  mlog_p  <- -log(1 - pmax(Q95, Q5))

  beta_95 <- matrix(Q95, nrow = N_locus, ncol = N_env, byrow = TRUE)
  beta_5  <- matrix(Q5,  nrow = N_locus, ncol = N_env, byrow = TRUE)
  beta_Mlogp <- -log(1 - pmax(beta_95, beta_5))
  beta_Mlogp[is.infinite(beta_Mlogp)] <- log(n_draws)

  imp_loci_env <- which(mlog_p > -log(1 - CI))
  if (length(imp_loci_env) == 0L) {
    cat("No adaptive loci found at CI =", CI, "\n")
    return(list(imp_loci_table = NULL, beta_Mlogp = beta_Mlogp,
                Q95 = Q95, Q5 = Q5))
  }

  imp_loci_summary <- beta_summary[imp_loci_env, ]
  imp_loci_table <- data.table::data.table(
    Env_var = env_names[as.integer(sub(".*\\[(.*?),.*", "\\1",
                                       imp_loci_summary$variable))],
    Locus   = as.integer(unlist(lapply(imp_loci_summary$variable,
                                       .split_string)))
  )

  list(imp_loci_table = imp_loci_table, beta_Mlogp = beta_Mlogp,
       Q95 = Q95, Q5 = Q5)
}

score_stan <- function(beta_draw, N_locus, N_env) {
  ## [N_locus x N_env] mlog_p = -log(1 - max_tail) matrix from beta draws.
  if (!inherits(beta_draw, "draws_matrix"))
    beta_draw <- posterior::as_draws_matrix(beta_draw)

  pn      <- colnames(beta_draw)
  is_beta <- grepl("^Beta\\[", pn)
  bd      <- beta_draw[, is_beta, drop = FALSE]
  pn_beta <- pn[is_beta]
  env_idx <- as.integer(sub(".*\\[([0-9]+),.*", "\\1", pn_beta))
  loc_idx <- as.integer(sub(".*,([0-9]+)\\].*", "\\1", pn_beta))

  n_draws <- nrow(bd)
  Q95     <- pmax(colMeans(bd > 0), 1 / n_draws)
  Q5      <- pmax(colMeans(bd < 0), 1 / n_draws)
  mlog_p  <- -log(1 - pmax(Q95, Q5))
  mlog_p[is.infinite(mlog_p)] <- log(n_draws)

  score <- matrix(NA_real_, nrow = N_locus, ncol = N_env)
  score[cbind(loc_idx, env_idx)] <- mlog_p
  score
}


## ============================================================================
## TIMING UTILITY
## ============================================================================

toc_hours <- function(tic, toc, msg) {
  ## Adaptive-unit elapsed-time formatter for tictoc::toc(func.toc = ...).
  secs  <- as.numeric(toc - tic)
  label <- if (is.null(msg) || msg == "") "" else paste0(msg, ": ")
  if (secs < 60) {
    paste0(label, signif(secs, 3), " sec elapsed")
  } else if (secs < 3600) {
    paste0(label, signif(secs / 60, 3), " min elapsed")
  } else {
    paste0(label, signif(secs / 3600, 3), " hr elapsed")
  }
}
