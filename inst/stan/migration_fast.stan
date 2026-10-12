// migration_fast.stan — faster drop-in for migration_eig.stan
//
// Uses the same Lyapunov-based O(N^4) coalescent_time algorithm as
// migration_eig.stan (see that file and migration_eig_derivation.pdf).
//
//   LYAPUNOV EQUATION
//   The two-lineage structured-coalescent system Q_rate t = -1, re-indexed
//   into the symmetric N x N matrix T, yields
//     M T + T M = 2(D - J),   D = diag(T_11,...,T_NN),  J = 11'.
//
//   EIGENBASIS DECOUPLING
//   With M = V Lambda V' (symmetric eigendecomposition):
//     T_tilde[a,b] = 2 D_tilde[a,b] / (lam[a] + lam[b]),  (a,b) != (z,z)
//   where z is the zero-eigenvalue index (M 1 = 0 => lam[z] = 0).
//
//   ZERO-MODE SINGULARITY (trace constraint sum_k T_kk = N^2)
//     A[i,k] = 2 * sum_{(a,b)!=(z,z)} V[i,a]*V[i,b]*V[k,a]*V[k,b] / (lam[a]+lam[b])
//     solve (I - A) v = 1,  T_tilde[z,z] = N^3 / (1'v),  d[k] = (T_tilde[z,z]/N)*v[k]
//
//   COMPLEXITY (vs migration.stan)
//     N =   8:  ~4x     N =  32:  ~47x     N = 128:  ~700x
//     Break-even around N = 5-6.
//
// Additional improvements over migration_eig.stan:
//
//   (1) coalescent_time() is in external C++ (coalescent_time.hpp) with a
//       custom analytical adjoint. Solves the adjoint Lyapunov equation
//         M Lambda + Lambda M - 2 diag(diag Lambda) = T_bar
//         trace(Lambda) = -(1/2) 1' T_bar 1
//       by reusing the cached eigendecomposition and (I-A) LU from the
//       forward pass. Eliminates O(N^3)-per-op autodiff through eigenvectors
//       entirely. Compile with:
//         cpp_options   = list(USER_HEADER = "<abs path>/coalescent_time.hpp")
//         stanc_options = list("allow-undefined")
//
//   (2) expected_D built from a single matrix expression (Edist = 4T - diag
//       broadcast) via precomputed index arrays (idx_i, idx_j in transformed
//       data), removing the nested loop from the autodiff tape.
//
//   (3) Non-centred parameterisation of m_ij (same as migration_eig.stan):
//         z ~ Exp(1),  m_ij = z .* exp(-mij_tilda)
//       Same model as m_ij ~ Exp(exp(mij_tilda)) but removes the funnel
//       between m_ij and gamma, typically cutting leapfrog steps by 2-10x.

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
  // Canonical pairwise-diversity parents (shared with migration_M3b.stan).
  vector[(N_patch * (N_patch - 1)) %/% 2] pi_off;
  vector[N_patch] pi_diag;
  matrix[(N_patch * (N_patch - 1)) %/% 2, N_cov] covariates;
  real <lower=0> gamma_sd;
}

transformed data{
  int N_pair = (N_patch * (N_patch - 1)) %/% 2;
  array[N_pair] int idx_i;
  array[N_pair] int idx_j;
  (idx_i, idx_j) = build_pair_indices(N_patch);

  // Petkova-centered combination derived from the pi parents: Dij = 4 * pi_off - pi_diag_i - pi_diag_j
  vector[N_pair] Dij;
  for (p in 1:N_pair) Dij[p] = 4 * pi_off[p] - pi_diag[idx_i[p]] - pi_diag[idx_j[p]];
}

parameters{
  vector <lower = 0> [N_pair] z;
  real <lower = 0> alpha, sigma;
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

  vector[N_patch] Tdiag = diagonal(T);
  matrix[N_patch, N_patch] Edist = 4*T - rep_matrix(Tdiag, N_patch) - rep_matrix(Tdiag', N_patch);
  vector[N_pair] expected_D;
  for(p in 1:N_pair) expected_D[p] = alpha * Edist[idx_i[p], idx_j[p]];

  Dij ~ normal(expected_D, sigma);
  z ~ exponential(1);
  append_row([gamma0]', gamma) ~ normal(0, gamma_sd);
  [alpha, sigma]' ~ exponential(1);
}

generated quantities{
  matrix[N_patch, N_patch] Migration = construct_M(m_ij, N_patch);
}
