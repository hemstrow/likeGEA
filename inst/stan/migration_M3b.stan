// migration_M3b.stan — M3b migration estimator (raw pi diagonal anchor)
//
// Derived from migration_fast.stan (same Lyapunov-based O(N^4) coalescent_time
// algorithm via the external coalescent_time.hpp header). The likelihood is
// rewritten per derivation_M3b_migration.pdf, §6:
//
//   (1) Two data blocks instead of the Petkova centered combination:
//         pi_off  : length N_pair, off-diagonal pairwise diversity (i<j)
//         pi_diag : length N_patch, within-deme heterozygosity
//
//   (2) Single shared scale c = 2 mu fitted jointly:
//         pi_off  ~ Normal(c * T_off,  sigma_off)
//         pi_diag ~ Normal(c * T_diag, sigma_diag)
//       The diagonal anchor breaks the (m, alpha) -> (k m, k alpha) scale
//       ridge that the Petkova combination is invariant under.
//
//   (3) Split residual SDs (sigma_off, sigma_diag) absorb the differing
//       per-locus sampling geometry of the diagonal vs off-diagonal estimators
//       (derivation §6.4).
//
//   (4) Distance-decay prior on log m_ij is unchanged from migration_fast:
//         log m_ij = gamma0 + covariates * gamma + z, z ~ Exp(1) non-centred.
//
// Compile with:
//   cpp_options   = list(USER_HEADER = "<abs path>/coalescent_time.hpp",
//                        stan_threads = TRUE)
//   stanc_options = list("allow-undefined")

functions {
  matrix construct_M (vector m_ij, int N_patch){
    matrix[N_patch, N_patch] M = rep_matrix(0.0, N_patch, N_patch);
    int k = 1;
    for(i in 1:(N_patch-1)){
      for(j in (i+1):N_patch){
        M[i,j] = m_ij[k];
        M[j,i] = M[i,j];
        k += 1;
      }
    }
    vector[N_patch] diag_M = M*rep_vector(1,N_patch);
    for(i in 1:N_patch) M[i,i] = -diag_M[i];
    return M;
  }

  // Upper-triangular pair indices (i, j) with i < j, in column-major order.
  // Returns the two length-N_pair int arrays as a tuple.
  tuple(array[] int, array[] int) build_pair_indices(int N_patch){
    int N_pair = (N_patch * (N_patch - 1)) %/% 2;
    array[N_pair] int idx_i;
    array[N_pair] int idx_j;
    int pos = 1;
    for(i in 1:(N_patch-1)){
      for(j in (i+1):N_patch){
        idx_i[pos] = i;
        idx_j[pos] = j;
        pos += 1;
      }
    }
    return (idx_i, idx_j);
  }

  matrix coalescent_time(matrix M);
}

data {
  int N_patch, N_cov;
  vector[(N_patch * (N_patch - 1)) %/% 2] pi_off;
  vector[N_patch] pi_diag;
  matrix[(N_patch * (N_patch - 1)) %/% 2, N_cov] covariates;
  real <lower = 0> gamma_sd;
}

transformed data{
  int N_pair = (N_patch * (N_patch - 1)) %/% 2;
  array[N_pair] int idx_i;
  array[N_pair] int idx_j;
  (idx_i, idx_j) = build_pair_indices(N_patch);
}

parameters{
  vector <lower = 0> [N_pair] z;
  real <lower = 0> c;
  real <lower = 0> sigma_off;
  real <lower = 0> sigma_diag;
  real gamma0;
  vector [N_cov] gamma;
}

transformed parameters{
  vector[N_pair] mij_tilda = gamma0 + covariates*gamma;
  vector <lower = 0> [N_pair] m_ij = z .* exp(-mij_tilda);
}

model{
  matrix[N_patch, N_patch] M = construct_M(m_ij, N_patch);
  matrix[N_patch, N_patch] T = coalescent_time(M);

  vector[N_pair] expected_pi_off;
  for(p in 1:N_pair) expected_pi_off[p] = c * T[idx_i[p], idx_j[p]];
  vector[N_patch] expected_pi_diag = c * diagonal(T);

  pi_off  ~ normal(expected_pi_off,  sigma_off);
  pi_diag ~ normal(expected_pi_diag, sigma_diag);

  z                   ~ exponential(1);
  append_row([gamma0]', gamma)     ~ normal(0, gamma_sd);
  c                   ~ exponential(1);
  sigma_off           ~ exponential(1);
  sigma_diag          ~ exponential(1);
}

generated quantities{
  matrix[N_patch, N_patch] Migration = construct_M(m_ij, N_patch);
}
