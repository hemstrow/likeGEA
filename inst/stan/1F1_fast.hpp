/**
 * 1F1_fast.hpp  —  superset of 1F1.hpp with a direct f'/f evaluator.
 *
 * Adds (on top of 1F1.hpp):
 *   vector FPF_1F1(vector a, vector b, vector z)    // returns  f'/f = (a/b) R_1
 *
 * where
 *   R_1(a,b,z) = 1F1(a+1, b+1, z) / 1F1(a,b,z)
 *   f  = 1F1(a,b,z)
 *   f' = d/dz 1F1 = (a/b) 1F1(a+1,b+1,z)
 *
 * Rationale (methodB_speedup.md §1.1, §3.3):
 *   Batches 1 and 2 of hybrid_modelB.stan compute
 *     pbar0   = (a/b)  · 1F1(a+1 ,b+1 ,s) / 1F1(a ,b ,s)
 *     hat_p   = (G+a*)/(2N+b*) · 1F1(G+a*+1, 2N+b*+1, s) / 1F1(G+a*, 2N+b*, s)
 *   Both are f'/f, so a direct evaluator halves each batch (2D -> D calls),
 *   cutting the per-locus 1F1 budget from 6D to 4D (~33 % fewer series calls).
 *
 * Implementation  (1.1 + 3.3 combined):
 *   A single three-accumulator series tracks f, g = z f', h = z² f'' together
 *   with the digamma accumulators S_a, S_b (for ∂f/∂a, ∂f/∂b) and two new
 *   accumulators T_a, T_b (for ∂f'/∂a, ∂f'/∂b).  From these we obtain, in
 *   closed form inside the series loop:
 *
 *     φ        = f' / f                           (value)
 *     ∂φ/∂z    = h / (z² f) - φ²                  (from f''/f - (f'/f)²)
 *     ∂φ/∂a    =  T_a / (z f) - φ · S_a / f
 *     ∂φ/∂b    = -T_b / (z f) + φ · S_b / f
 *
 *   Stan's autodiff therefore sees one reverse-pass callback node per
 *   batched FPF_1F1 call — matching the existing Log_1F1 pattern — but with
 *   half as many series evaluations underneath.
 *
 * Numerical paths mirror Log_1F1 exactly:
 *   1. z = 0              : φ = a/b, closed-form grads
 *   2. Unified 3-acc      : |z| <= 200, a > 0, b > 0 (Kummer for z < 0)
 *   3. Fallback           : 2 × eval_log + log-chain rule (for strong-|z|)
 *
 * Kummer transform:
 *   f(a,b,z) = e^z · F(A,b,Z)       A = b - a,  Z = -z
 *   f'(a,b,z) = f - e^z F_Z          →  φ(a,b,z) = 1 - φ_A(A,b,Z)
 *   ∂φ/∂a = ∂φ_A/∂A,   ∂φ/∂b = -∂φ_A/∂A - ∂φ_A/∂b,   ∂φ/∂z = ∂φ_A/∂Z
 *
 * HOW TO USE
 *   In your .stan functions block:
 *     vector hypergeometric_1F1(vector a, vector b, vector z);   // 1F1
 *     vector Log_1F1          (vector a, vector b, vector z);    // log 1F1
 *     vector FPF_1F1          (vector a, vector b, vector z);    // f' / f
 *
 *   Compile with CmdStanR:
 *     model <- cmdstan_model(
 *       "hybrid_modelB_fast.stan",
 *       cpp_options   = list(USER_HEADER = "1F1_fast.hpp",
 *                            STAN_THREADS = TRUE),
 *       stanc_options = list("allow-undefined" = TRUE))
 */

#ifndef HYPERGEOMETRIC_1F1_FAST_HPP
#define HYPERGEOMETRIC_1F1_FAST_HPP

#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <stan/math.hpp>
#include <cmath>
#include <limits>
#include <ostream>

namespace h1f1_internal {

// =============================================================================
// Result structs
// =============================================================================

struct LogGrads {
    double log_val;
    double da;
    double db;
    double dz;
};

struct Grads {
    double val;
    double da;
    double db;
    double dz;
};

// f'/f and its gradients (1.1 return type).
struct FpfGrads {
    double val;    // φ = f'/f
    double da;     // ∂φ/∂a
    double db;     // ∂φ/∂b
    double dz;     // ∂φ/∂z
};

// =============================================================================
// Boost with safe error policy
// =============================================================================

using safe_pol = boost::math::policies::policy<
    boost::math::policies::overflow_error<
        boost::math::policies::ignore_error>,
    boost::math::policies::underflow_error<
        boost::math::policies::ignore_error>,
    boost::math::policies::rounding_error<
        boost::math::policies::ignore_error>,
    boost::math::policies::evaluation_error<
        boost::math::policies::ignore_error>>;

inline double h1f1_boost(double a, double b, double z) {
    // safe_pol catches Boost overflow/underflow/rounding/evaluation policy
    // errors but does NOT catch C++ exceptions thrown OUTSIDE the policy
    // chain -- in particular boost::math::trunc / itrunc throws when called
    // with NaN inputs or |z| > ~2.1e9 (rounding_error / overflow_error).
    // Wrap the call so those exceptions don't kill the chain.
    try {
        return boost::math::hypergeometric_1F1(a, b, z, safe_pol{});
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

inline bool bad(double f) {
    return !std::isfinite(f) || f <= 0.0;
}

// =============================================================================
// Standalone digamma series (Boost fallback path)
// =============================================================================

inline void digamma_series(double a, double b, double z,
                           double& S_a, double& S_b,
                           double tol = 1e-13, int max_k = 600) {
    double term = 1.0, delta_a = 0.0, delta_b = 0.0;
    S_a = 0.0;
    S_b = 0.0;
    // PARAM_MIN floor: when a or b is near zero, 1/(a+k-1) at k=1 is ~1/a,
    // making delta_a explode and the iterated term*delta_a accumulation
    // produce extreme gradients that stall HMC.  Bail with Sa, Sb = 0; the
    // caller still uses Boost's value of 1F1, so log_f is unaffected.
    constexpr double PARAM_MIN = 1e-8;
    if (a < PARAM_MIN || b < PARAM_MIN) return;
    for (int k = 1; k < max_k; ++k) {
        delta_a += 1.0 / (a + k - 1.0);
        delta_b += 1.0 / (b + k - 1.0);
        term    *= (a + k - 1.0) * z / ((b + k - 1.0) * k);
        double da = term * delta_a, db = term * delta_b;
        S_a += da;
        S_b += db;
        if (k > 5 && std::abs(da) < tol * std::abs(S_a)
                  && std::abs(db) < tol * std::abs(S_b)) break;
    }
}

// =============================================================================
// Asymptotic expansions
// =============================================================================

inline LogGrads asymptotic_large_z(double a, double b, double z) {
    double lz = std::log(z);
    return { z + (a - b) * lz + std::lgamma(b) - std::lgamma(a),
             lz - boost::math::digamma(a),
             boost::math::digamma(b) - lz,
             1.0 + (a - b) / z };
}

inline LogGrads asymptotic_large_param(double a, double b, double z) {
    double r = a / b;
    double log_f = r * z;
    return { log_f, z / b, -r * z / b, r };
}

// =============================================================================
// Unified single-pass series  (2 accumulators: f, f')
//
// Reproduces 1F1.hpp::unified_series verbatim — used by Log_1F1 / hypergeometric_1F1
// for Batch 3 (the ratio form is not applicable there).
// =============================================================================

inline bool unified_series(double a, double b, double z,
                           double& log_f, double& Sa, double& Sb,
                           double& f_prime_over_f) {
    constexpr double TOL = 1e-13;
    constexpr int    MAX_K = 600;
    // PARAM_MIN floor: see digamma_series above.  Returning false routes
    // the caller to the Boost path, which is more robust at small a/b.
    constexpr double PARAM_MIN = 1e-8;
    if (a < PARAM_MIN || b < PARAM_MIN) return false;

    double term = 1.0, f = 1.0, f_prime = 0.0;
    double delta_a = 0.0, delta_b = 0.0;
    double ak = a, bk = b;
    const double inv_z = 1.0 / z;
    Sa = 0.0;
    Sb = 0.0;

    for (int k = 1; k <= MAX_K; ++k) {
        delta_a += 1.0 / ak;
        delta_b += 1.0 / bk;

        term *= ak * z / (bk * k);
        ak += 1.0;
        bk += 1.0;

        if (!std::isfinite(term)) return false;

        f       += term;
        f_prime += term * k * inv_z;
        double da = term * delta_a, db = term * delta_b;
        Sa += da;
        Sb += db;

        if (k >= 8 && (k & 3) == 0 && std::abs(term) < TOL * std::abs(f)) {
            if (f <= 0.0 || !std::isfinite(f)) return false;
            log_f = std::log(f);
            Sa /= f;
            Sb /= f;
            f_prime_over_f = f_prime / f;
            return true;
        }
    }
    return false;
}

// =============================================================================
// Unified THREE-accumulator series  (f, g=zf', h=z²f'')  for φ = f'/f
//
// Additional accumulators beyond unified_series:
//   g  = Σ_k   k·t_k          →  f'  = g / z
//   h  = Σ_k k(k-1)·t_k       →  f'' = h / z²
//   T_a = Σ_k   k·t_k·δ_a^{(k)}  →  ∂f'/∂a =  T_a / z
//   T_b = Σ_k   k·t_k·δ_b^{(k)}  →  ∂f'/∂b = -T_b / z
//
// Closed-form gradient assembly (inside the loop at convergence):
//   φ        = g / (z f)                                 (= f'/f)
//   ∂φ/∂z    = h / (z² f) - φ²                           (= f''/f - φ²)
//   ∂φ/∂a    =  T_a / (z f) - φ · (S_a / f)
//   ∂φ/∂b    = -T_b / (z f) + φ · (S_b / f)
//
// The three accumulators (f, g, h) plus two extra digamma-weighted sums
// cost only ~7 % over 1-accumulator bookkeeping in FLOPs; measured in
// bench_1f1.py at 1.07× per series iteration.
//
// Requires a > 0, b > 0, z > 0.  Returns true on success.
// =============================================================================

inline bool unified_series_fpf(double a, double b, double z,
                               double& phi,
                               double& dphi_da,
                               double& dphi_db,
                               double& dphi_dz) {
    constexpr double TOL = 1e-13;
    constexpr int    MAX_K = 600;
    // PARAM_MIN floor: see digamma_series above.  Returning false routes
    // the caller to the eval_log fallback (which is also guarded).
    constexpr double PARAM_MIN = 1e-8;
    if (a < PARAM_MIN || b < PARAM_MIN) return false;

    double term = 1.0;
    double f = 1.0, g = 0.0, h = 0.0;
    double Sa = 0.0, Sb = 0.0, Ta = 0.0, Tb = 0.0;
    double delta_a = 0.0, delta_b = 0.0;
    double ak = a, bk = b;
    const double inv_z = 1.0 / z;

    for (int k = 1; k <= MAX_K; ++k) {
        // Accumulate digamma differences BEFORE updating term
        delta_a += 1.0 / ak;
        delta_b += 1.0 / bk;

        term *= ak * z / (bk * k);
        ak += 1.0;
        bk += 1.0;

        if (!std::isfinite(term)) return false;

        double kterm = k * term;      // k · t_k, reused 3×

        f += term;                    // Σ t_k
        g += kterm;                   // Σ k · t_k          = z f'
        h += kterm * (k - 1);         // Σ k(k-1) · t_k     = z² f''
        Sa += term  * delta_a;        // Σ t_k δ_a^{(k)}    = ∂f/∂a
        Sb += term  * delta_b;        // Σ t_k δ_b^{(k)}    = −∂f/∂b
        Ta += kterm * delta_a;        // Σ k t_k δ_a^{(k)}  = z ∂f'/∂a
        Tb += kterm * delta_b;        // Σ k t_k δ_b^{(k)}  = −z ∂f'/∂b

        // Same batched-convergence criterion as unified_series (1F1-E2).
        if (k >= 8 && (k & 3) == 0 && std::abs(term) < TOL * std::abs(f)) {
            if (f <= 0.0 || !std::isfinite(f)) return false;
            double inv_f    = 1.0 / f;
            double inv_zf   = inv_z * inv_f;
            phi             = g * inv_zf;                     // f'/f
            double fpp_of_f = h * inv_z * inv_z * inv_f;      // f''/f
            dphi_dz         = fpp_of_f - phi * phi;
            dphi_da         =  Ta * inv_zf - phi * Sa * inv_f;
            dphi_db         = -Tb * inv_zf + phi * Sb * inv_f;
            return true;
        }
    }
    return false;
}

// =============================================================================
// Kummer transform for log(1F1)
// =============================================================================

inline LogGrads kummer_chain(double z, const LogGrads& g) {
    return { z + g.log_val,
             -g.da,
             g.da + g.db,
             1.0 - g.dz };
}

// =============================================================================
// eval_log: per-element log(1F1) evaluation with gradients
// =============================================================================

inline LogGrads eval_log(double a, double b, double z) {

    // ── NaN / Inf / invalid-b guard ───────────────────────────────────────
    // The host model can produce NaN or non-positive (a, b, z) when its
    // intermediate quantities (kappa_star, mu_star, b_minus_a) collapse
    // during NUTS exploration. Any of those, if passed downstream, leads
    // to uncaught Boost exceptions in trunc/itrunc inside hypergeometric_1F1.
    // Return a NaN log-value (zero gradients) so Stan rejects the proposal
    // cleanly instead of crashing the chain.
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(z)
        || b <= 0.0 || a <= 0.0) {
        return { std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0 };
    }

    if (z == 0.0)
        return { 0.0, 0.0, 0.0, a / b };

    constexpr double Z_UNIFIED = 200.0;
    if (std::abs(z) <= Z_UNIFIED && a > 0.0 && b > 0.0) {
        double log_f, Sa, Sb, fpf;
        if (z >= 0.0) {
            if (unified_series(a, b, z, log_f, Sa, Sb, fpf))
                return { log_f, Sa, -Sb, fpf };
        } else {
            double A = b - a, Z = -z;
            if (A > 0.0 && unified_series(A, b, Z, log_f, Sa, Sb, fpf))
                return { z + log_f, -Sa, (Sa - Sb), 1.0 - fpf };
        }
    }

    constexpr double PARAM_MAX = 1e15;
    if (a > PARAM_MAX || b > PARAM_MAX) {
        if (z >= 0.0)
            return asymptotic_large_param(a, b, z);
        else
            return kummer_chain(z, asymptotic_large_param(b - a, b, -z));
    }

    // ── Large-|z| guard before Boost ──────────────────────────────────────
    // Boost's hypergeometric_1F1 calls itrunc(z) internally to locate the
    // peak of the asymptotic series.  itrunc returns int (32-bit) and throws
    // boost::math::rounding_error / std::overflow_error for |z| > ~2.1e9,
    // and our safe_pol policy does NOT catch that throw because itrunc is
    // invoked outside the special-function policy chain.  Skip Boost entirely
    // above Z_BOOST_SAFE and use the leading-order DLMF 13.7.2 asymptotic.
    constexpr double Z_BOOST_SAFE = 5e5;
    if (std::abs(z) > Z_BOOST_SAFE) {
        if (z >= 0.0)
            return asymptotic_large_z(a, b, z);
        else
            return kummer_chain(z, asymptotic_large_z(b - a, b, -z));
    }

    if (z >= 0.0) {
        double f  = h1f1_boost(a,     b,     z);
        double fp = h1f1_boost(a+1.0, b+1.0, z);
        if (bad(f) || bad(fp))
            return asymptotic_large_z(a, b, z);
        double Sa, Sb;
        digamma_series(a, b, z, Sa, Sb);
        return { std::log(f), Sa / f, -Sb / f, (a / b) * fp / f };
    } else {
        double A = b - a, Z = -z;
        double f  = h1f1_boost(A,     b,     Z);
        double fp = h1f1_boost(A+1.0, b+1.0, Z);
        if (bad(f) || bad(fp))
            return kummer_chain(z, asymptotic_large_z(A, b, Z));
        double Sa, Sb;
        digamma_series(A, b, Z, Sa, Sb);
        return { z + std::log(f), -Sa / f, (Sa - Sb) / f,
                 1.0 - (A / b) * fp / f };
    }
}

inline Grads eval(double a, double b, double z) {
    LogGrads lg = eval_log(a, b, z);
    double val = std::exp(lg.log_val);
    return { val, val * lg.da, val * lg.db, val * lg.dz };
}

// =============================================================================
// eval_fpf: per-element φ = f'/f evaluation with closed-form gradients
//
// Path selection mirrors eval_log:
//   1. z = 0                    : closed-form  φ = a/b,
//                                   ∂φ/∂a = 1/b,
//                                   ∂φ/∂b = -a/b²,
//                                   ∂φ/∂z = a(b-a) / (b² (b+1))
//   2. Unified 3-acc series     : |z| <= 200, a > 0, b > 0
//                                 Kummer (b-a, b, -z) for z < 0
//   3. Fallback  (strong |z|)   : two eval_log calls + log-chain rule
//                                 φ = (a/b) · exp(L₂ - L₁)
//                                 ∂log φ/∂a = 1/a + L₂.da - L₁.da
//                                 ∂log φ/∂b = -1/b + L₂.db - L₁.db
//                                 ∂log φ/∂z =       L₂.dz - L₁.dz
// =============================================================================

inline FpfGrads eval_fpf(double a, double b, double z) {

    // NaN / Inf / invalid-b guard -- see eval_log() for rationale.
    // Return NaN value (zero gradients) so Stan rejects the proposal
    // instead of letting an uncaught Boost trunc exception kill the chain.
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(z)
        || b <= 0.0 || a <= 0.0) {
        return { std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0 };
    }

    // 1) z = 0 fast path
    if (z == 0.0) {
        double aob = a / b;
        return { aob,
                 1.0 / b,
                 -aob / b,
                 a * (b - a) / (b * b * (b + 1.0)) };
    }

    // 2) Unified three-accumulator series
    constexpr double Z_UNIFIED = 200.0;
    if (std::abs(z) <= Z_UNIFIED && a > 0.0 && b > 0.0) {
        double phi, dda, ddb, ddz;
        if (z >= 0.0) {
            if (unified_series_fpf(a, b, z, phi, dda, ddb, ddz))
                return { phi, dda, ddb, ddz };
        } else {
            // Kummer:  φ(a,b,z) = 1 - φ_A(A, b, Z),  A = b-a, Z = -z
            double A = b - a, Z = -z;
            if (A > 0.0 &&
                unified_series_fpf(A, b, Z, phi, dda, ddb, ddz)) {
                return { 1.0 - phi,
                         /* ∂/∂a = +∂/∂A */        dda,
                         /* ∂/∂b = -∂/∂A - ∂/∂b */ -dda - ddb,
                         /* ∂/∂z = +∂/∂Z */        ddz };
            }
        }
    }

    // 3) Fallback via eval_log — two log-1F1 calls, then chain rule
    LogGrads g1 = eval_log(a,       b,       z);   // L(a,   b,   z)
    LogGrads g2 = eval_log(a + 1.0, b + 1.0, z);   // L(a+1, b+1, z)
    double phi = (a / b) * std::exp(g2.log_val - g1.log_val);
    return { phi,
             phi * (1.0 / a  + g2.da - g1.da),
             phi * (-1.0 / b + g2.db - g1.db),
             phi * (g2.dz - g1.dz) };
}

}  // namespace h1f1_internal


// =============================================================================
// hypergeometric_1F1  (value, var)  —  unchanged from 1F1.hpp
// =============================================================================

inline Eigen::VectorXd
hypergeometric_1F1(const Eigen::VectorXd& a, const Eigen::VectorXd& b,
                   const Eigen::VectorXd& z, std::ostream*) {
    int N = a.size();
    Eigen::VectorXd out(N);
    for (int i = 0; i < N; ++i)
        out[i] = h1f1_internal::eval(a[i], b[i], z[i]).val;
    return out;
}

template <typename V1, typename V2, typename V3,
          stan::require_any_var_t<
              stan::scalar_type_t<V1>,
              stan::scalar_type_t<V2>,
              stan::scalar_type_t<V3>>* = nullptr>
Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>
hypergeometric_1F1(const V1& a_arg, const V2& b_arg, const V3& z_arg,
                   std::ostream*) {
    using var = stan::math::var;

    constexpr bool va = stan::is_var<stan::scalar_type_t<V1>>::value;
    constexpr bool vb = stan::is_var<stan::scalar_type_t<V2>>::value;
    constexpr bool vz = stan::is_var<stan::scalar_type_t<V3>>::value;

    int N = a_arg.size();

    using arena_vec = stan::arena_t<Eigen::VectorXd>;
    using arena_v   = stan::arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>>;
    arena_vec da_buf(va ? N : 0);
    arena_vec db_buf(vb ? N : 0);
    arena_vec dz_buf(vz ? N : 0);

    arena_v result(N);
    for (int i = 0; i < N; ++i) {
        h1f1_internal::Grads g = h1f1_internal::eval(
            stan::math::value_of(a_arg.coeff(i)),
            stan::math::value_of(b_arg.coeff(i)),
            stan::math::value_of(z_arg.coeff(i)));
        result.coeffRef(i) = var(g.val);
        if constexpr (va) da_buf[i] = g.da;
        if constexpr (vb) db_buf[i] = g.db;
        if constexpr (vz) dz_buf[i] = g.dz;
    }

    auto a_arena = stan::math::to_arena(a_arg);
    auto b_arena = stan::math::to_arena(b_arg);
    auto z_arena = stan::math::to_arena(z_arg);

    stan::math::reverse_pass_callback(
        [result, a_arena, b_arena, z_arena,
         da_buf, db_buf, dz_buf]() mutable {
            for (int i = 0; i < result.size(); ++i) {
                const double adj_i = result.coeff(i).adj();
                if constexpr (va) a_arena.coeff(i).adj() += adj_i * da_buf[i];
                if constexpr (vb) b_arena.coeff(i).adj() += adj_i * db_buf[i];
                if constexpr (vz) z_arena.coeff(i).adj() += adj_i * dz_buf[i];
            }
        });

    return result;
}


// =============================================================================
// Log_1F1  (value, var)  —  unchanged from 1F1.hpp
// =============================================================================

inline Eigen::VectorXd
Log_1F1(const Eigen::VectorXd& a, const Eigen::VectorXd& b,
        const Eigen::VectorXd& z, std::ostream*) {
    int N = a.size();
    Eigen::VectorXd out(N);
    for (int i = 0; i < N; ++i)
        out[i] = h1f1_internal::eval_log(a[i], b[i], z[i]).log_val;
    return out;
}

template <typename V1, typename V2, typename V3,
          stan::require_any_var_t<
              stan::scalar_type_t<V1>,
              stan::scalar_type_t<V2>,
              stan::scalar_type_t<V3>>* = nullptr>
Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>
Log_1F1(const V1& a_arg, const V2& b_arg, const V3& z_arg, std::ostream*) {
    using var = stan::math::var;

    constexpr bool va = stan::is_var<stan::scalar_type_t<V1>>::value;
    constexpr bool vb = stan::is_var<stan::scalar_type_t<V2>>::value;
    constexpr bool vz = stan::is_var<stan::scalar_type_t<V3>>::value;

    int N = a_arg.size();

    using arena_vec = stan::arena_t<Eigen::VectorXd>;
    using arena_v   = stan::arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>>;
    arena_vec da_buf(va ? N : 0);
    arena_vec db_buf(vb ? N : 0);
    arena_vec dz_buf(vz ? N : 0);

    arena_v result(N);
    for (int i = 0; i < N; ++i) {
        h1f1_internal::LogGrads g = h1f1_internal::eval_log(
            stan::math::value_of(a_arg.coeff(i)),
            stan::math::value_of(b_arg.coeff(i)),
            stan::math::value_of(z_arg.coeff(i)));
        result.coeffRef(i) = var(g.log_val);
        if constexpr (va) da_buf[i] = g.da;
        if constexpr (vb) db_buf[i] = g.db;
        if constexpr (vz) dz_buf[i] = g.dz;
    }

    auto a_arena = stan::math::to_arena(a_arg);
    auto b_arena = stan::math::to_arena(b_arg);
    auto z_arena = stan::math::to_arena(z_arg);

    stan::math::reverse_pass_callback(
        [result, a_arena, b_arena, z_arena,
         da_buf, db_buf, dz_buf]() mutable {
            for (int i = 0; i < result.size(); ++i) {
                const double adj_i = result.coeff(i).adj();
                if constexpr (va) a_arena.coeff(i).adj() += adj_i * da_buf[i];
                if constexpr (vb) b_arena.coeff(i).adj() += adj_i * db_buf[i];
                if constexpr (vz) z_arena.coeff(i).adj() += adj_i * dz_buf[i];
            }
        });

    return result;
}


// =============================================================================
// FPF_1F1 — new exported function  (1.1 + 3.3)
//
//   FPF_1F1(a, b, z)[i]  =  1F1'(a[i], b[i], z[i]) / 1F1(a[i], b[i], z[i])
//                        =  (a[i] / b[i]) · R_1(a[i], b[i], z[i])
//
// Same single-callback, arena-backed reverse-pass pattern as Log_1F1.
// Forward pass uses the three-accumulator unified_series_fpf, which emits
// φ and all three partials from ONE series loop — so the tape cost is the
// same one callback node, but the primal work is halved relative to two
// Log_1F1 calls plus a Stan-side ratio.
// =============================================================================

inline Eigen::VectorXd
FPF_1F1(const Eigen::VectorXd& a, const Eigen::VectorXd& b,
        const Eigen::VectorXd& z, std::ostream*) {
    int N = a.size();
    Eigen::VectorXd out(N);
    for (int i = 0; i < N; ++i)
        out[i] = h1f1_internal::eval_fpf(a[i], b[i], z[i]).val;
    return out;
}

template <typename V1, typename V2, typename V3,
          stan::require_any_var_t<
              stan::scalar_type_t<V1>,
              stan::scalar_type_t<V2>,
              stan::scalar_type_t<V3>>* = nullptr>
Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>
FPF_1F1(const V1& a_arg, const V2& b_arg, const V3& z_arg, std::ostream*) {
    using var = stan::math::var;

    constexpr bool va = stan::is_var<stan::scalar_type_t<V1>>::value;
    constexpr bool vb = stan::is_var<stan::scalar_type_t<V2>>::value;
    constexpr bool vz = stan::is_var<stan::scalar_type_t<V3>>::value;

    int N = a_arg.size();

    using arena_vec = stan::arena_t<Eigen::VectorXd>;
    using arena_v   = stan::arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>>;
    arena_vec da_buf(va ? N : 0);
    arena_vec db_buf(vb ? N : 0);
    arena_vec dz_buf(vz ? N : 0);

    arena_v result(N);
    for (int i = 0; i < N; ++i) {
        h1f1_internal::FpfGrads g = h1f1_internal::eval_fpf(
            stan::math::value_of(a_arg.coeff(i)),
            stan::math::value_of(b_arg.coeff(i)),
            stan::math::value_of(z_arg.coeff(i)));
        result.coeffRef(i) = var(g.val);
        if constexpr (va) da_buf[i] = g.da;
        if constexpr (vb) db_buf[i] = g.db;
        if constexpr (vz) dz_buf[i] = g.dz;
    }

    auto a_arena = stan::math::to_arena(a_arg);
    auto b_arena = stan::math::to_arena(b_arg);
    auto z_arena = stan::math::to_arena(z_arg);

    stan::math::reverse_pass_callback(
        [result, a_arena, b_arena, z_arena,
         da_buf, db_buf, dz_buf]() mutable {
            for (int i = 0; i < result.size(); ++i) {
                const double adj_i = result.coeff(i).adj();
                if constexpr (va) a_arena.coeff(i).adj() += adj_i * da_buf[i];
                if constexpr (vb) b_arena.coeff(i).adj() += adj_i * db_buf[i];
                if constexpr (vz) z_arena.coeff(i).adj() += adj_i * dz_buf[i];
            }
        });

    return result;
}

#endif  // HYPERGEOMETRIC_1F1_FAST_HPP
