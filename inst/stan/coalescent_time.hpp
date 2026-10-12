// coalescent_time.hpp
//
// External C++ implementation of the `coalescent_time(matrix M)` Stan
// function declared (without a body) in migration_fast.stan.
// See migration_eig.stan and migration_eig_derivation.pdf for the full
// algorithm derivation. Summary:
//
//   FORWARD PASS (on plain doubles, off the autodiff tape)
//   The two-lineage coalescence-time matrix T solves the Lyapunov equation
//     M T + T M = 2(D - J),   D = diag(T_11,...,T_NN),  J = 11'.
//   In the eigenbasis M = V Lambda V' this decouples to
//     T_tilde[a,b] = 2 D_tilde[a,b] / (lam[a] + lam[b]),  (a,b) != (z,z).
//   The zero-mode T_tilde[z,z] is pinned by the trace constraint
//   sum_k T_kk = N^2 via the closure solve:
//     build A (O(N^4), done as a single DSYRK call  A = -W W'),
//     factor (I - A) once with LDLT (cached for backward pass),
//     solve (I - A) v = 1,  T_tilde[z,z] = N^3 / (1' v),
//     d[k] = (T_tilde[z,z] / N) * v[k].
//   Then T_tilde[a,b] = 2 D_tilde[a,b] / (lam[a]+lam[b]) for off-modes,
//   and T = V T_tilde V'.
//
//   BACKWARD PASS (custom analytical adjoint, reuses cached forward state)
//   The adjoint Lyapunov equation has the same operator as the forward:
//     M Lambda + Lambda M - 2 diag(diag Lambda) = T_bar,
//     trace(Lambda) = -(1/2) 1' T_bar 1,
//   with solution  M_bar = -(Lambda T + T Lambda).
//   Solved by reusing the cached eigendecomposition, inverse-denominator
//   matrix 1/(lam[a]+lam[b]), and (I-A) LDLT from the forward pass.
//   No new eigendecomp or LDLT per HMC gradient step.
//
//   Stan-visible templated overloads follow the standard stan::math idiom:
//   store inputs/outputs on the autodiff arena via make_chainable_ptr,
//   register the gradient via reverse_pass_callback, accumulate adjoints
//   via whole-matrix Eigen expressions instead of N^2 scalar loops.
//
// Build
// -----
// Supplied to cmdstanr / cmdstan via
//   cpp_options   = list(USER_HEADER = "<abs path to this file>"),
//   stanc_options = list("allow-undefined")

#ifndef COALESCENT_TIME_HPP
#define COALESCENT_TIME_HPP

#include <stan/math.hpp>
#include <stan/math/rev.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <ostream>
#include <utility>

// Helpers live in a generic, filename-agnostic sub-namespace. The public
// Stan-callable overloads of coalescent_time(...) are defined at the top
// level below, so that cmdstanr's USER_HEADER include drops them into
// whichever <stan_filename>_model_namespace is generated. This decouples
// the .hpp from the .stan filename.
namespace coalescent_internal {

// ----------------------------------------------------------------------
// Cached forward state. Everything the backward pass needs from the
// forward pass lives here. Stashed on the autodiff arena via
// stan::math::make_chainable_ptr so it survives until reverse pass.
// ----------------------------------------------------------------------
struct CoalCache {
  int N;
  int z;                                                // index of zero eigenvalue
  Eigen::VectorXd lam;                                  // eigenvalues, ascending
  Eigen::MatrixXd V;                                    // eigenvectors (columns)
  Eigen::MatrixXd inv_denom;                            // 1/(λα+λβ) with the (z,z) entry zeroed
  Eigen::LDLT<Eigen::MatrixXd> IA_ldlt;                 // LDLT factorisation of I − A (PD)
  Eigen::VectorXd v_sol;                                // (I − A)^{-1} 1
  double S;                                             // 1ᵀ v_sol
  Eigen::MatrixXd T;                                    // forward output (cached for the M̄ formula)
};

// ----------------------------------------------------------------------
// Forward pass on plain doubles. If `cache` is non-null, the eigendecomp,
// inverse-denominator matrix, LDLT factorisation, v_sol, S and T are
// stored in it for reuse in the backward pass.
// ----------------------------------------------------------------------
inline Eigen::MatrixXd coalescent_time_forward(const Eigen::MatrixXd& M,
                                               CoalCache* cache) {
  const int N = static_cast<int>(M.rows());

  // (1) Symmetric eigendecomposition. SelfAdjointEigenSolver returns
  // eigenvalues in ascending order; since M is negative semidefinite
  // with exactly one zero mode, that zero lives at index N-1. No search.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  Eigen::VectorXd lam = es.eigenvalues();
  Eigen::MatrixXd V   = es.eigenvectors();
  const int z = N - 1;

  // (2) Inverse denominator 1/(λα+λβ), with the (z,z) entry zeroed so
  // that downstream cwiseProduct calls do not need a branch. Vectorised
  // broadcast + elementwise inverse; reused twice in the backward pass.
  Eigen::MatrixXd inv_denom
      = (lam.replicate(1, N) + lam.transpose().replicate(N, 1))
            .array()
            .inverse()
            .matrix();
  inv_denom(z, z) = 0.0;

  // (3) Build A via ONE symmetric rank-k update (SYRK) instead of a
  // loop of N(N+1)/2 small rank-1 updates.
  //
  // The closure matrix is  A = Σ_{(a,b)≠(z,z)} c_{ab} · u_{ab} u_{ab}ᵀ
  // with u_{aa} = V[:,a]⊙V[:,a], c_{aa}=1/λ_a
  // and  u_{ab} = V[:,a]⊙V[:,b], c_{ab}=4/(λ_a+λ_b) for a<b.
  //
  // Every c_{ab} is strictly negative (M neg semidef ⇒ λ_a < 0 for
  // a≠z). Absorbing sign and magnitude into a single column scaling,
  //       W = U · diag(√|c|)   ⇒   A = −W Wᵀ,
  // which is a single DSYRK call:
  //       A.selfadjointView<Lower>().rankUpdate(W, −1.0)
  // This touches N(N+1)/2 output entries and is one BLAS call, versus
  // O(N²) tiny rank-1 updates in the previous version.
  const int P = N * (N + 1) / 2 - 1;           // exclude (z,z)
  Eigen::MatrixXd W(N, P);
  {
    int col = 0;
    for (int a = 0; a < N; ++a) {
      for (int b = a; b < N; ++b) {
        if (a == z && b == z) continue;
        const double c_ab = (a == b)
                              ? 1.0 / lam[a]
                              : 4.0 / (lam[a] + lam[b]);
        // c_ab < 0 ⇒ scale by sqrt(-c_ab)
        W.col(col) =
            V.col(a).cwiseProduct(V.col(b)) * std::sqrt(-c_ab);
        ++col;
      }
    }
  }

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);
  A.selfadjointView<Eigen::Lower>().rankUpdate(W, -1.0);   // A += -W Wᵀ
  A.triangularView<Eigen::StrictlyUpper>() = A.transpose();

  // (4) Closure solve: factor (I − A) once and cache the LDLT for reuse
  // in the backward pass. A = -W Wᵀ is negative semi-definite (W Wᵀ is
  // PSD), so (I − A) = I + W Wᵀ is strictly positive definite with
  // smallest eigenvalue ≥ 1. LDLT is ~2× faster than PartialPivLU for
  // both factorisation and solve, and its pivoting keeps it robust.
  Eigen::MatrixXd IA = -A;
  IA.diagonal().array() += 1.0;
  Eigen::LDLT<Eigen::MatrixXd> IA_ldlt(IA);
  Eigen::VectorXd v_sol   = IA_ldlt.solve(Eigen::VectorXd::Ones(N));
  const double S          = v_sol.sum();
  const double T_tilde_zz = static_cast<double>(N) * N * N / S;

  // Diagonal of T in the original basis.
  Eigen::VectorXd d = (T_tilde_zz / static_cast<double>(N)) * v_sol;

  // (5) Build T̃ = 2 D̃ / (λa+λb), then rotate back T = V T̃ Vᵀ.
  //   D̃ = Vᵀ diag(d) V    — one GEMM after a row-scale.
  //   T̃ symmetric ⇒ V·T̃ uses SSYMM (saves ~⅓ flops vs GEMM).
  Eigen::MatrixXd D_tilde;
  D_tilde.noalias() = V.transpose() * d.asDiagonal() * V;
  Eigen::MatrixXd T_tilde = 2.0 * D_tilde.cwiseProduct(inv_denom);
  T_tilde(z, z) = T_tilde_zz;

  Eigen::MatrixXd VTt;
  VTt.noalias() = V * T_tilde.selfadjointView<Eigen::Lower>();  // SSYMM
  Eigen::MatrixXd T;
  T.noalias() = VTt * V.transpose();                             // GEMM
  // Rounding-level symmetrisation (T is symmetric by construction).
  T = 0.5 * (T + T.transpose());

  if (cache) {
    cache->N         = N;
    cache->z         = z;
    cache->lam       = std::move(lam);
    cache->V         = std::move(V);
    cache->inv_denom = std::move(inv_denom);
    cache->IA_ldlt   = std::move(IA_ldlt);
    cache->v_sol     = std::move(v_sol);
    cache->S         = S;
    cache->T         = T;
  }
  return T;
}

// ----------------------------------------------------------------------
// Backward pass: solve the adjoint Lyapunov equation (same operator as
// the forward Lyapunov equation, so we reuse `c.lam`, `c.V`,
// `c.inv_denom`, `c.IA_ldlt`, `c.v_sol`, `c.S`) and return M̄.
// ----------------------------------------------------------------------
inline Eigen::MatrixXd coalescent_time_backward(const Eigen::MatrixXd& Tbar_in,
                                                const CoalCache& c) {
  const int N                      = c.N;
  const int z                      = c.z;
  const Eigen::MatrixXd& V         = c.V;
  const Eigen::MatrixXd& inv_denom = c.inv_denom;

  // Symmetric part of incoming adjoint — only the symmetric part of
  // T̄ contributes to the gradient since T is symmetric.
  Eigen::MatrixXd Tbar = 0.5 * (Tbar_in + Tbar_in.transpose());

  // Rotate to eigenbasis: T̃̄ = Vᵀ T̄ V. T̄ is symmetric, so the first
  // product uses SSYMM (~⅔ the flops of a general GEMM).
  Eigen::MatrixXd TbarV;
  TbarV.noalias() = Tbar.selfadjointView<Eigen::Lower>() * V;    // SSYMM
  Eigen::MatrixXd Tbar_tilde;
  Tbar_tilde.noalias() = V.transpose() * TbarV;                  // GEMM

  // H = T̃̄ / (λa+λb), (z,z) zeroed via cached inv_denom.
  Eigen::MatrixXd H = Tbar_tilde.cwiseProduct(inv_denom);

  // g_i = (V H Vᵀ)_{ii}. Computed as row-sum of (V H) ⊙ V, replacing
  // the second matmul with an O(N²) elementwise op.
  Eigen::MatrixXd VH;
  VH.noalias() = V * H;
  Eigen::VectorXd g = VH.cwiseProduct(V).rowwise().sum();

  // (I − A) w = g  (reuses the cached forward LDLT; one back-solve).
  Eigen::VectorXd w = c.IA_ldlt.solve(g);

  // Trace constraint:  trace(Λ) = -(1/2) 1ᵀ T̄ 1.
  const double cval = -0.5 * Tbar.sum();
  const double Lambda_tilde_zz =
      static_cast<double>(N) * (cval - w.sum()) / c.S;

  // Diagonal of Λ in the original basis.
  Eigen::VectorXd e = (Lambda_tilde_zz / static_cast<double>(N)) * c.v_sol + w;

  // Λ̃ in eigenbasis. Ẽ = Vᵀ diag(e) V — single GEMM after row-scale.
  Eigen::MatrixXd E_tilde;
  E_tilde.noalias() = V.transpose() * e.asDiagonal() * V;
  Eigen::MatrixXd Lambda_tilde =
      (Tbar_tilde + 2.0 * E_tilde).cwiseProduct(inv_denom);
  Lambda_tilde(z, z) = Lambda_tilde_zz;

  // Rotate back Λ = V Λ̃ Vᵀ. Λ̃ symmetric ⇒ V·Λ̃ uses SSYMM.
  Eigen::MatrixXd VLt;
  VLt.noalias() = V * Lambda_tilde.selfadjointView<Eigen::Lower>(); // SSYMM
  Eigen::MatrixXd Lambda;
  Lambda.noalias() = VLt * V.transpose();                            // GEMM
  Lambda = 0.5 * (Lambda + Lambda.transpose());

  // M̄ = -(Λ T + T Λ). Λ and T are both symmetric, so T Λ = (Λ T)ᵀ
  // and only ONE matmul is needed. Since Λ is symmetric, use SSYMM.
  Eigen::MatrixXd LT;
  LT.noalias() = Lambda.selfadjointView<Eigen::Lower>() * c.T;       // SSYMM
  Eigen::MatrixXd Mbar = -(LT + LT.transpose());
  return Mbar;
}

}  // namespace coalescent_internal

// ----------------------------------------------------------------------
// Stan-visible overloads. Signature follows the Stan external-function
// ABI:  T  coalescent_time(const M& M, std::ostream* pstream__);
//
// These sit at top level so cmdstanr's USER_HEADER include lands them
// inside whichever model namespace the generated .cpp opens — i.e.,
// they work regardless of the consuming .stan filename.
// ----------------------------------------------------------------------

// Plain double overload (Stan calls this from generated_quantities,
// transformed_data, double-fast paths, etc.).
template <typename T0__,
          stan::require_eigen_matrix_dynamic_t<T0__>* = nullptr,
          stan::require_not_st_var<T0__>* = nullptr>
inline Eigen::Matrix<double, -1, -1>
coalescent_time(const T0__& M, std::ostream* /*pstream__*/) {
  Eigen::MatrixXd Md = M.template cast<double>();
  return coalescent_internal::coalescent_time_forward(Md, nullptr);
}

// Autodiff (var) overload. This is what is called from the model block.
//
// We deliberately return Eigen::Matrix<var, -1, -1> (one vari per
// entry) rather than var_value<Eigen::MatrixXd>: Stan's generated model
// code declares local `matrix[N,N]` variables as
// Eigen::Matrix<local_scalar_t__, -1, -1>, and the whole-matrix
// `stan::model::assign(Matrix<var>, var_value<MatrixXd>, name)`
// overload does not exist in cmdstan 2.38, so returning a var_value
// triggers an overload-resolution failure on the LHS=RHS type mismatch.
// The N² constant from individual vari nodes is negligible vs. the
// O(N³) work this op does.
template <typename T0__,
          stan::require_eigen_matrix_dynamic_t<T0__>* = nullptr,
          stan::require_st_var<T0__>* = nullptr>
inline Eigen::Matrix<stan::math::var, -1, -1>
coalescent_time(const T0__& M_in, std::ostream* /*pstream__*/) {
  using stan::math::var;
  using stan::arena_t;                        // note: stan::, not stan::math::
  using stan::math::make_chainable_ptr;
  using stan::math::reverse_pass_callback;

  // Pin the (var-valued) input on the autodiff arena so we can read /
  // accumulate adjoints from the reverse-pass callback.
  arena_t<stan::math::promote_scalar_t<var, T0__>> arena_M = M_in;
  const int N = static_cast<int>(arena_M.rows());

  // Forward pass on doubles. Stan's Eigen plugin gives arena_M.val() as
  // a CwiseUnaryOp view of the underlying doubles, so this is just an
  // expression assignment — no manual copy.
  auto cache_ptr = make_chainable_ptr(coalescent_internal::CoalCache{});
  Eigen::MatrixXd Tval = coalescent_internal::coalescent_time_forward(arena_M.val(), cache_ptr);

  // Allocate the output as fresh vars on the autodiff arena via a
  // single Eigen expression: each var(double) constructor pushes a
  // vari onto the chainable stack with .val_ = Tval(i,j) and
  // .adj_ = 0. The arena_matrix `=` operator copies the resulting
  // Matrix<var,-1,-1> into arena-allocated storage.
  arena_t<Eigen::Matrix<var, -1, -1>> arena_T
      = Tval.unaryExpr([](double x) { return var(x); });

  reverse_pass_callback([arena_M, arena_T, cache_ptr]() mutable {
    // Pull the upstream adjoint as a single Eigen expression. The
    // .adj() view is provided by Stan's Eigen plugin and works for any
    // Matrix<var> derived expression including arena_matrix.
    Eigen::MatrixXd Tbar = arena_T.adj();

    Eigen::MatrixXd Mbar = coalescent_internal::coalescent_time_backward(Tbar, *cache_ptr);

    // One whole-matrix accumulation into the input adjoints; same idiom
    // as stan/math/rev/fun/trace_quad_form.hpp et al.
    arena_M.adj() += Mbar;
  });

  return Eigen::Matrix<var, -1, -1>(arena_T);
}

#endif  // COALESCENT_TIME_HPP
