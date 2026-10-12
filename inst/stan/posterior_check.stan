// gl_posterior_check.stan -- N_eff / G_eff based posterior-predictive check.
//
// Simplified from the full gl-weighted Beta-Binomial version to use the
// moment-matched (G_eff, 2*N_eff) sufficient statistics that
// structured.stan uses for inference.  Posterior check is
// now internally consistent with the main fit: both use the same data
// representation, so this script tests "does the model's posterior on p
// fit the moment-matched data?" without re-introducing the moment-match
// approximation issue.
//
// Per-locus likelihood at patch k (continuous-G binomial, up to a constant
// in p):
//   G_eff[i, k] | p[i, k] ~ Binomial(2*N_eff[i, k], p[i, k])
//
// Logit-Normal regression on the locus's env covariate (Jacobian-
// corrected so the <lower=0,upper=1> p transform cancels):
//   logit(p[i, k]) ~ Normal(alpha[i] + beta[i] * Env[i, k], sigma)
//
// Each map_rect shard receives one locus:
//   p_alpha_beta = c(p[i, ], alpha[i], beta[i])         (vector length N_patch+2)
//   real_data    = c(G_eff[i, ], TwoN_eff[i, ], Env[i, ])    (length 3*N_patch)
//   empty int data

functions{
  vector low_coverage_whole_genome_obs_pdf(vector sigma_shared, vector p_alpha_beta, array[] real real_data, array[] int empty_int_data){
    int N_patch = num_elements(p_alpha_beta) - 2;
    real sigma  = sigma_shared[1];

    vector [N_patch] p_slice = segment(p_alpha_beta, 1, N_patch);
    real alpha_i = p_alpha_beta[N_patch + 1];
    real beta_i  = p_alpha_beta[N_patch + 2];

    vector [N_patch] G_eff    = to_vector(real_data[1:N_patch]);
    vector [N_patch] TwoN_eff = to_vector(real_data[(N_patch + 1):(2*N_patch)]);
    vector [N_patch] Env_vec  = to_vector(real_data[(2*N_patch + 1):(3*N_patch)]);

    // Pre-compute logit(p) and log(1-p) once; each is reused twice below.
    vector [N_patch] logit_p = logit(p_slice);
    vector [N_patch] log1m_p = log1m(p_slice);

    // Continuous-G binomial log-density in exponential-family form (up to a constant in p). The identity
    //   G*log(p) + (2N - G)*log(1-p)  =  G*logit(p) + 2N*log(1-p)
    // makes logit(p) the natural parameter and log(1-p) the log-partition per trial. Patches with TwoN_eff = 0 contribute zero.
    real lq = sum(G_eff .* logit_p + TwoN_eff .* log1m_p);

    // Logit-Normal regression on env (folded into the per-locus shard).
    lq += normal_lpdf(logit_p | alpha_i + beta_i * Env_vec, sigma);

    // Jacobian correction: cancel Stan's <lower=0,upper=1> p transform so the prior on p is proper LogitNormal on the absolute scale.
    //   -log(p) - log(1-p) = -[logit(p) + log(1-p)] - log(1-p)
    //                      = -logit(p) - 2*log(1-p)
    lq += -sum(logit_p) - 2 * sum(log1m_p);

    return [lq]';
  }
}

data {
  int <lower = 1> N_patch;
  int <lower = 1> N_locus;
  array [N_locus, N_patch] real <lower = 0> G_eff;
  array [N_locus, N_patch] real <lower = 0> TwoN_eff;
  array [N_locus] vector [N_patch] Env;
}

transformed data {
  // Pack (G_eff, TwoN_eff, Env) per locus into a single 3*N_patch real array so map_rect can read all three in one shard.
  array [N_locus, 3*N_patch] real real_data;
  array [N_locus, 0]         int  empty_int_data;
  for (i in 1:N_locus) {
    for (j in 1:N_patch) {
      real_data[i, j]             = G_eff[i, j];
      real_data[i, N_patch + j]   = TwoN_eff[i, j];
      real_data[i, 2*N_patch + j] = Env[i, j];
    }
  }
}

parameters {
  array [N_locus] vector <lower = 0, upper = 1> [N_patch] p;
  vector [N_locus] beta;
  vector [N_locus] alpha;
  real <lower = 0> sigma;
}

model {
  array [N_locus] vector [N_patch + 2] p_alpha_beta;
  for (i in 1:N_locus) p_alpha_beta[i] = append_row(p[i], [alpha[i], beta[i]]');

  target += sum(map_rect(low_coverage_whole_genome_obs_pdf, [sigma]', p_alpha_beta, real_data, empty_int_data));

  beta  ~ normal(0, 1);
  alpha ~ normal(0, 1);
  sigma ~ exponential(1);
}
