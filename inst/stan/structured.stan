functions{
  matrix sparce_beta(int N_env, int N_locus, real c_square, real tau, matrix lambda_unif, matrix z_beta){
    matrix [N_env, N_locus] lambda = tan(lambda_unif);
    matrix [N_env, N_locus] lambda_tilde = sqrt(c_square*square(lambda) .* inv(c_square + square (tau* lambda)));
    return (z_beta.*lambda_tilde)*tau ;
  }

  vector compute_rho(int N_env, int N_locus, real c_square, real tau, matrix lambda_unif, vector rho_tilda){
    matrix [N_env, N_locus] lambda = tan(lambda_unif);
    matrix [N_env, N_locus] lambda_tilde = sqrt(c_square*square(lambda) .* inv(c_square + square (tau* lambda)));
    vector [N_locus] rho;
    for (l in 1:N_locus) rho[l] = rho_tilda[l] * max(lambda_tilde[, l]) * tau * inv_sqrt(c_square);
    return rho;
  }

  vector ObsEvoCountPdf(vector tau_csq, vector zbeta_lambda_mu_kappa_rhotilda, data array[] real Gcall_TwoN, array[] int empty_int_data){
    // Unpacking data and parameters
    real tau      = tau_csq[1];
    real c_square = tau_csq[2];
    #include structured_data.stan
    vector [N_env] z_beta = segment(zbeta_lambda_mu_kappa_rhotilda, 1, N_env);
    vector [N_env] lambda = tan(segment(zbeta_lambda_mu_kappa_rhotilda, N_env + 1, N_env));
    real mu = zbeta_lambda_mu_kappa_rhotilda[2*N_env + 1];
    real kappa = zbeta_lambda_mu_kappa_rhotilda[2*N_env + 2];
    real rho_tilda = zbeta_lambda_mu_kappa_rhotilda[2*N_env + 3];
    vector [N_patch] TwoN   = to_vector(segment(Gcall_TwoN, N_patch+1, N_patch));
    vector [N_patch] G_call = to_vector(segment(Gcall_TwoN, 1, N_patch));

    // regularize sensitivity coefficients and rho
    vector [N_env] lambda_tilda = sqrt(c_square*square(lambda) .* inv(c_square + square (tau* lambda)));
    vector [N_env] beta_sigma = lambda_tilda*tau;
    vector [N_env] beta = z_beta.*beta_sigma;
    real rho  = rho_tilda * max(lambda_tilda) * tau * inv_sqrt(c_square);

    // Selection coefficient per deme  s_k = (E beta)_k  (identity Q Q' = I)
    vector [N_patch] s_coef = Env * beta;

    // ---- Step 1: Off-diagonal row sums  S_k = -sum_{j} Q_{kj}^2 d_j ----------
    //            = sum_{k' != k} m_tilde_{k,k'}   (positive)
    vector [N_patch] S          = -(Q .* Q) * d;
    vector [N_patch] kappa_star = kappa + rho * S;

    // ---- Step 2: Zero-order 1F1 means  pbar0_k = (a/b) 1F1(a+1,b+1,s)/1F1(a,b,s)
    //             with a = mu*kappa, b = kappa.  FPF_1F1 returns (a/b)·R_1.
    real a = mu * kappa;
    real b = kappa;
    vector [N_patch] pbar0 = FPF_1F1(rep_vector(a, N_patch), rep_vector(b, N_patch), s_coef);
    // vector[2*N_patch] b1_a = append_row(rep_vector(a,     N_patch), rep_vector(a + 1, N_patch));
    // vector[2*N_patch] b1_b = append_row(rep_vector(b,     N_patch), rep_vector(b + 1, N_patch));
    // vector[2*N_patch] b1_z = append_row(s_coef, s_coef);
    // vector[2*N_patch] log_1F1_b1 = Log_1F1(b1_a, b1_b, b1_z);
    // vector [N_patch] pbar0 = (a / b) * exp(segment(log_1F1_b1, N_patch+1, N_patch) - segment(log_1F1_b1, 1, N_patch));

    // ---- Step 3: 1F1-based mean-field shape (mu_star, kappa_star) -------------
    //   mu_star_k = (mu*kappa + rho * [M_offdiag pbar0]_k) / kappa_star_k
    // At rho = 0 this collapses to mu, as required for the unstructured limit.
    vector [N_patch] M_pbar0 = Q * (d .* (Q' * pbar0)) + S .* pbar0;
    vector [N_patch] mu_star = (mu * kappa + rho * M_pbar0) ./ kappa_star;
    vector [N_patch] a_star  = mu_star .* kappa_star;
    vector [N_patch] b_star  = kappa_star;

    // kappa_star - a_star = kappa_star * (1 - mu_star) = (1 - mu) * kappa + rho * (S - M_pbar0)
    vector [N_patch] b_minus_a = (1 - mu) * kappa + rho * (S - M_pbar0);

    // ---- Step 6: 1F1-BetaBinomial pseudo-likelihood at (a_star, b_star) -------
    //   log P = lbeta(G + a_star, 2N - G + b_star - a_star) - lbeta(a_star, b_star - a_star)
    //           + log 1F1(G + a_star, 2N + b_star, s) - log 1F1(a_star, b_star, s)
    vector[2*N_patch] b3_a = append_row(a_star, a_star + G_call);
    vector[2*N_patch] b3_b = append_row(kappa_star, TwoN + kappa_star);
    vector[2*N_patch] b3_z = append_row(s_coef, s_coef);
    vector[2*N_patch] log_1F1_b3 = Log_1F1(b3_a, b3_b, b3_z);
    real lq = sum(lbeta(G_call + a_star, TwoN - G_call + b_minus_a) - lbeta(a_star, b_minus_a));
    lq -= sum(segment(log_1F1_b3, 1, N_patch));   // subtract denominator 1F1(a_star, b_star, s)
    lq += sum(segment(log_1F1_b3, N_patch+1, N_patch));   // add numerator    1F1(G + a_star, 2N + b_star, s)
    lq += normal_lpdf(kappa | 0, k_sd) + std_normal_lpdf(z_beta) + normal_lpdf(rho_tilda | 0, 1);
    return [lq]';
  }

  // External C++ implementations via Boost with analytical gradients (1F1_fast.hpp).
  vector Log_1F1(vector a, vector b, vector z); //   Log_1F1 : log 1F1(a,b,z)
  vector FPF_1F1(vector a, vector b, vector z); //   FPF_1F1 : 1F1'(a,b,z) / 1F1(a,b,z) = (a/b)·R_1
}

data {
  real p0;
  real <lower = 0> k_sd;           // kappa ~ normal(0, k_sd)
  int < lower = 1 > N_patch, N_locus, N_env ;
  array [N_locus, N_patch] real <lower = 0> G_call, TwoN;
}

transformed data{
  real f = p0*inv(N_locus);
  real scale_global = f*inv(sqrt(N_patch)*N_env);

  array [N_locus, 2*N_patch] real <lower = 0> Gcall_TwoN;
  array [N_locus, 0] int empty_int_data;

  Gcall_TwoN[1:N_locus, 1:N_patch] = G_call;
  Gcall_TwoN[1:N_locus, (N_patch+1):(2*N_patch)] = TwoN;
}

parameters {
  vector <lower = 0, upper = 1> [N_locus] mu;
  vector <lower = 0> [N_locus] kappa;
  vector <lower = 0> [N_locus] rho_tilda;
  real <lower = 0, upper = pi()/2> tau_unif;
  real <lower = 0> c_square;
  matrix [N_env, N_locus] z_beta;
  matrix <lower = 0, upper = pi()/2> [N_env, N_locus] lambda_unif ;
}

model {
  array [N_locus] vector [2*N_env + 3] zbeta_lambda_mu_kappa_rhotilda;
  for (i in 1:N_locus) zbeta_lambda_mu_kappa_rhotilda[i] = append_row(append_row(z_beta[, i], lambda_unif[,i]), [mu[i], kappa[i], rho_tilda[i]]');

  real tau = scale_global * tan(tau_unif);
  target += sum(map_rect(ObsEvoCountPdf, [tau, c_square]', zbeta_lambda_mu_kappa_rhotilda, Gcall_TwoN, empty_int_data));
  c_square ~ inv_gamma(2, 4);
}

generated quantities{
  real Tau = scale_global * tan(tau_unif);
  matrix [N_env, N_locus] Beta = sparce_beta(N_env, N_locus, c_square, Tau, lambda_unif, z_beta);
  vector [N_locus] rho = compute_rho(N_env, N_locus, c_square, Tau, lambda_unif, rho_tilda);
}
