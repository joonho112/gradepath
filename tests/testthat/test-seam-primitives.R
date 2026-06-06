# =============================================================================
# test-seam-primitives.R  --  parity of the ebrecipe seam
# -----------------------------------------------------------------------------
# For EACH reused primitive, call the gradepath `.gp_eb_*` wrapper AND the raw
# `getFromNamespace(".eb_*", "ebrecipe")` internal on the SAME frozen small
# input, and assert they are bit-for-bit identical (`expect_identical`). The
# wrapper is a thin pass-through, so identity must hold exactly -- any drift
# means the wrapper is transforming, which Ch7 forbids.
#
# Also asserts `.gp_eb_input` returns an `eb_estimates` on krw_firms race cols.
#
# All frozen inputs match the REAL installed-ebrecipe (0.5.0) signatures
# (verified via args() + bodies): basis is .eb_spline_basis(support, n_knots);
# the softmax density / penalized-loglik / constraint solves take Q + alpha; the
# softmax density returns list(g, log_g); the pushforward is single-support
# (support, g, s, psi_1, psi_2) and REQUIRES characteristic/standardization_model;
# the W primitive takes an eb_estimates + an eb_prior.
# =============================================================================

skip_if_not_installed("ebrecipe")

# ---- frozen inputs ---------------------------------------------------------

eb <- function(nm) getFromNamespace(nm, "ebrecipe")

support <- seq(-2, 2, length.out = 9L)

# natural-spline (ns) basis on the support: REAL sig is (support, n_knots)
basis <- eb(".eb_spline_basis")(support, n_knots = 5L)

# spline coefficients: one per basis column
set.seed(2718L)
alpha <- stats::rnorm(ncol(basis))

# firm-level estimates (small frozen vectors)
theta_hat <- c(-0.8, -0.2, 0.1, 0.5, 1.1)
s         <- c(0.30, 0.45, 0.25, 0.60, 0.40)

# normal-mixture (log-)likelihood matrix on the support
log_P <- eb(".eb_normal_mixture_matrix")(theta_hat, s, support, log = TRUE)

# free coefficients for the mean-constraint solves: length = ncol(basis) - 1
set.seed(1414L)
alpha_free  <- stats::rnorm(ncol(basis) - 1L)
target_mean <- 0

# a small matrix for row-wise log-sum-exp
lse_mat <- matrix(
  c(0, log(2), log(3), 0, -1, 2),
  nrow = 2L, byrow = TRUE
)

# a normalized density on the support for the one-level pushforward cross-check.
# NOTE: .eb_softmax_density returns list(g, log_g); the pushforward kernel takes
# the bare normalized mass vector `g`.
g_dens <- eb(".eb_softmax_density")(basis, alpha)$g

# a real eb_estimates + eb_prior for the W primitive
est_in   <- ebrecipe::eb_input(theta_hat = theta_hat, s = s)
prior_in <- ebrecipe::eb_deconvolve(est_in)


# ---- parity tests, one per reused primitive --------------------------------

test_that(".gp_eb_spline_basis matches the ebrecipe internal", {
  raw <- eb(".eb_spline_basis")(support, n_knots = 5L)
  expect_identical(.gp_eb_spline_basis(support, n_knots = 5L), raw)
})

test_that(".gp_eb_softmax_density matches the ebrecipe internal", {
  raw <- eb(".eb_softmax_density")(basis, alpha)
  expect_identical(.gp_eb_softmax_density(basis, alpha), raw)
})

test_that(".gp_eb_normal_mixture_matrix matches the ebrecipe internal", {
  raw <- eb(".eb_normal_mixture_matrix")(theta_hat, s, support, log = TRUE)
  expect_identical(
    .gp_eb_normal_mixture_matrix(theta_hat, s, support, log = TRUE),
    raw
  )
})

test_that(".gp_eb_penalized_loglik matches the ebrecipe internal", {
  raw <- eb(".eb_penalized_loglik")(alpha, basis, log_P, penalty = 1.5)
  expect_identical(
    .gp_eb_penalized_loglik(alpha, basis, log_P, penalty = 1.5),
    raw
  )
})

test_that(".gp_eb_solve_alpha_T matches the ebrecipe internal", {
  raw <- eb(".eb_solve_alpha_T")(alpha_free, basis, support, target_mean)
  expect_identical(
    .gp_eb_solve_alpha_T(alpha_free, basis, support, target_mean),
    raw
  )
})

test_that(".gp_eb_full_alpha matches the ebrecipe internal", {
  raw <- eb(".eb_full_alpha")(alpha_free, basis, support, target_mean)
  expect_identical(
    .gp_eb_full_alpha(alpha_free, basis, support, target_mean),
    raw
  )
})

test_that(".gp_eb_row_log_sum_exp matches the ebrecipe internal", {
  raw <- eb(".eb_row_log_sum_exp")(lse_mat)
  expect_identical(.gp_eb_row_log_sum_exp(lse_mat), raw)
})

test_that(".gp_eb_pushforward_theta matches the ebrecipe internal (one-level cross-check, invariant 8)", {
  # .eb_pushforward_theta REQUIRES characteristic or standardization_model to
  # pick the scale (white/multiplicative vs male/additive); it errors with
  # neither. Use characteristic = "white" (the race / multiplicative path).
  raw <- eb(".eb_pushforward_theta")(support, g_dens, s = 0.4,
                                     psi_1 = 0, psi_2 = 0,
                                     characteristic = "white")
  expect_identical(
    .gp_eb_pushforward_theta(support, g_dens, s = 0.4, psi_1 = 0, psi_2 = 0,
                             characteristic = "white"),
    raw
  )
})

test_that(".gp_eb_posterior_weights matches the ebrecipe internal (the W primitive)", {
  raw <- eb(".eb_posterior_weights")(est_in, prior_in)
  expect_identical(.gp_eb_posterior_weights(est_in, prior_in), raw)
})


# ---- stage-1 entry: .gp_eb_input returns an eb_estimates -------------------

test_that(".gp_eb_input returns an eb_estimates on krw_firms race columns", {
  skip_if_not(exists("krw_firms"))
  est <- .gp_eb_input(
    theta_hat = krw_firms$theta_hat_race,
    s         = krw_firms$se_race,
    unit_id   = krw_firms$firm_id
  )
  expect_s3_class(est, "eb_estimates")
  # faithful pass-through: identical to a direct ebrecipe::eb_input call
  expect_identical(
    est,
    ebrecipe::eb_input(
      theta_hat = krw_firms$theta_hat_race,
      s         = krw_firms$se_race,
      unit_id   = krw_firms$firm_id
    )
  )
})
