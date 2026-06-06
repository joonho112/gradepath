# Tests for the native beta-GMM estimation core (W-seam).
#
# STRUCTURAL correctness on the bundled example (ebrecipe::krw_firms).
# These tests assert the SHAPE, ALGEBRA, and INVARIANTS of the estimator (moment
# formulas, Omega PSD, the single-N sandwich, df == 1, coupling carriers, error
# classes, determinism) using the bundled krw_firms example as fixture data.
# NOTE: krw_firms is NOT KRW's actual GMM input (see estimation-input.R; sorted
# correlation ~0.94, different scale), so the beta VALUE here is not the published
# 0.510 and is deliberately NOT asserted in this file. KRW Table-3 parity (race
# beta = 0.510, gender = 1.255) is gated in test-estimation-parity.R, which feeds
# KRW's real GMM input via gp_krw_gmm_input().
#
# NOTE on inputs: an `eb_estimates` (the `ebrecipe::eb_input` container) carries
# `theta_hat` and `s` but NO `characteristic` and NO group column, so:
#   * make_*_input() call .gp_eb_input(theta_hat=, s=) (the ebrecipe arg is `s`,
#     not `se`), with NO `group=` / `characteristic=` args;
#   * `characteristic` is passed explicitly to gp_estimation_core();
#   * the moment-conditions data matrix is built with .gp_estimation_data_matrix().

if (!exists("make_race_input", mode = "function")) {
  make_race_input <- function() {
    .gp_eb_input(
      theta_hat = gradepath::krw_firms$theta_hat_race,
      s = gradepath::krw_firms$se_race,
      unit_id = gradepath::krw_firms$firm_id
    )
  }
}
if (!exists("make_gender_input", mode = "function")) {
  make_gender_input <- function() {
    .gp_eb_input(
      theta_hat = gradepath::krw_firms$theta_hat_gender,
      s = gradepath::krw_firms$se_gender,
      unit_id = gradepath::krw_firms$firm_id
    )
  }
}

is_psd <- function(M, tol = 1e-8) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  all(ev >= -tol)
}

# --- moment conditions ------------------------------------------------------

test_that("race moment vector has length 4 and g_i is N x 4", {
  inp <- make_race_input()
  data <- .gp_estimation_data_matrix(inp)
  delta <- gp_gmm_start("race")
  mc <- gp_moment_conditions(delta, data, "race")
  expect_length(mc$g, 4L)
  expect_equal(dim(mc$g_i), c(nrow(data), 4L))
  expect_true(all(is.finite(mc$g)))
  expect_true(all(is.finite(mc$g_i)))
})

test_that("gender moment vector has length 4 and g_i is N x 4", {
  inp <- make_gender_input()
  data <- .gp_estimation_data_matrix(inp)
  delta <- gp_gmm_start("gender")
  mc <- gp_moment_conditions(delta, data, "gender")
  expect_length(mc$g, 4L)
  expect_equal(dim(mc$g_i), c(nrow(data), 4L))
  expect_true(all(is.finite(mc$g)))
})

test_that("moment conditions match the hand-coded Matlab formulas (race)", {
  ## Independent re-derivation of moment_conditions.m for race.
  inp <- make_race_input()
  data <- .gp_estimation_data_matrix(inp)
  theta_hat <- data[, 2]; s <- data[, 3]
  delta <- c(-1.0, -1.4, 0.5)
  mu <- exp(delta[1]); sigma_xi <- exp(delta[2]); beta <- delta[3]
  sigma <- sqrt(sigma_xi^2)            # group_fx == 0
  v_hat <- theta_hat / s^beta
  s_v <- s^(1 - beta)
  r <- (v_hat - mu) / sqrt(sigma^2 + s_v^2)
  g_ref <- c(mean(r), mean(r * s), mean(r^2 - 1), mean((r^2 - 1) * s))
  mc <- gp_moment_conditions(delta, data, "race")
  expect_equal(mc$g, g_ref, tolerance = 1e-12)
  expect_equal(mc$v_hat, v_hat, tolerance = 1e-12)
  expect_equal(mc$s_v, s_v, tolerance = 1e-12)
})

test_that("moment conditions match the hand-coded Matlab formulas (gender)", {
  inp <- make_gender_input()
  data <- .gp_estimation_data_matrix(inp)
  theta_hat <- data[, 2]; s <- data[, 3]
  delta <- c(0.1, -1.5, 0.5)
  mu <- delta[1]; sigma_xi <- exp(delta[2]); beta <- delta[3]
  sigma <- sqrt(sigma_xi^2)
  v_hat <- (theta_hat - mu) / s^beta
  s_v <- s^(1 - beta)
  r <- v_hat / sqrt(sigma^2 + s_v^2)
  g_ref <- c(mean(r), mean(r * s), mean(r^2 - 1), mean((r^2 - 1) * s))
  mc <- gp_moment_conditions(delta, data, "gender")
  expect_equal(mc$g, g_ref, tolerance = 1e-12)
})

test_that("race and gender residual formulas are genuinely asymmetric", {
  ## Multiplicative (race) vs additive (gender) must NOT coincide: same delta,
  ## same data, different r. Guards against accidental 'symmetrization'.
  inp <- make_race_input()
  data <- .gp_estimation_data_matrix(inp)
  delta <- c(-1.0, -1.4, 0.5)
  r_race <- gp_moment_conditions(delta, data, "race")$r
  r_gen  <- gp_moment_conditions(delta, data, "gender")$r
  expect_false(isTRUE(all.equal(r_race, r_gen)))
})

test_that("group_fx == 1 is rejected as out of scope", {
  inp <- make_race_input()
  data <- .gp_estimation_data_matrix(inp)
  expect_error(
    gp_moment_conditions(gp_gmm_start("race"), data, "race", group_fx = 1L),
    class = "gradepath_scope_error"
  )
})

test_that("unknown characteristic in moment conditions errors", {
  inp <- make_race_input()
  data <- .gp_estimation_data_matrix(inp)
  expect_error(
    gp_moment_conditions(c(-1, -1.5, 0.5), data, "ethnicity"),
    class = "gradepath_validation_error"
  )
})

# --- Omega ------------------------------------------------------------------

test_that("Omega is 4x4, symmetric, and PSD (race)", {
  inp <- make_race_input()
  data <- .gp_estimation_data_matrix(inp)
  mc <- gp_moment_conditions(gp_gmm_start("race"), data, "race", get_cov = TRUE)
  expect_equal(dim(mc$Omega), c(4L, 4L))
  expect_equal(mc$Omega, t(mc$Omega), tolerance = 1e-12)
  expect_true(is_psd(mc$Omega))
})

test_that("Omega is 4x4, symmetric, and PSD (gender)", {
  inp <- make_gender_input()
  data <- .gp_estimation_data_matrix(inp)
  mc <- gp_moment_conditions(gp_gmm_start("gender"), data, "gender", get_cov = TRUE)
  expect_equal(dim(mc$Omega), c(4L, 4L))
  expect_equal(mc$Omega, t(mc$Omega), tolerance = 1e-12)
  expect_true(is_psd(mc$Omega))
})

test_that("gp_moment_cov equals the get_cov Omega", {
  inp <- make_race_input()
  data <- .gp_estimation_data_matrix(inp)
  mc <- gp_moment_conditions(gp_gmm_start("race"), data, "race", get_cov = TRUE)
  expect_equal(gp_moment_cov(mc$g_i), mc$Omega, tolerance = 1e-12)
})

# --- start values -----------------------------------------------------------

test_that("race start is the exact Matlab one-step start", {
  expect_equal(gp_gmm_start("race"), c(-1.1789, -1.5761, 0.5099))
})

test_that("gender start is length 3 and finite", {
  st <- gp_gmm_start("gender")
  expect_length(st, 3L)
  expect_true(all(is.finite(st)))
})

# --- single-weighting minimizer ---------------------------------------------

test_that("gp_gmm_min runs under identity weight and returns delta of length 3", {
  inp <- make_race_input()
  data <- .gp_estimation_data_matrix(inp)
  fit <- gp_gmm_min(data, "race", diag(4), gp_gmm_start("race"))
  expect_length(fit$delta, 3L)
  expect_true(is.finite(fit$objective))
  expect_true(fit$objective >= 0)
})

# --- two-step driver: RACE --------------------------------------------------

test_that("two-step race fit runs, returns a finite beta and correct shapes", {
  fit <- gp_estimation_core(make_race_input(), characteristic = "race")
  expect_length(fit$delta, 3L)
  expect_true(is.finite(fit$beta))
  expect_true(is.finite(fit$mu) && fit$mu > 0)        # race mu = exp(.) > 0
  expect_true(is.finite(fit$sigma_xi) && fit$sigma_xi > 0)
  expect_equal(fit$characteristic, "race")
})

test_that("race fit is deterministic across repeated runs", {
  f1 <- gp_estimation_core(make_race_input(), characteristic = "race")
  f2 <- gp_estimation_core(make_race_input(), characteristic = "race")
  expect_equal(f1$beta, f2$beta, tolerance = 1e-10)
  expect_equal(f1$delta, f2$delta, tolerance = 1e-10)
})

test_that("race J-stat is finite, non-negative, df == 1", {
  fit <- gp_estimation_core(make_race_input(), characteristic = "race")
  expect_true(is.finite(fit$J_stat))
  expect_gte(fit$J_stat, 0)
  expect_equal(fit$df, 1L)
})

test_that("race sandwich C is symmetric PSD and SEs are positive finite", {
  fit <- gp_estimation_core(make_race_input(), characteristic = "race")
  expect_equal(dim(fit$C), c(3L, 3L))
  expect_equal(fit$C, t(fit$C), tolerance = 1e-10)
  expect_true(is_psd(fit$C))
  expect_length(fit$SE, 3L)
  expect_true(all(is.finite(fit$SE)))
  expect_true(all(fit$SE > 0))
})

test_that("race coupling: m_hat length 2 = [exp(mu_raw), sigma_xi]; V_m sym finite", {
  fit <- gp_estimation_core(make_race_input(), characteristic = "race")
  expect_length(fit$m_hat, 2L)
  expect_true(all(is.finite(fit$m_hat)))
  expect_equal(fit$m_hat[1], exp(fit$delta[1]), tolerance = 1e-10)  # exp(mu_raw)
  expect_equal(fit$m_hat[2], exp(fit$delta[2]), tolerance = 1e-10)  # sigma_xi
  expect_equal(dim(fit$V_m), c(2L, 2L))
  expect_equal(fit$V_m, t(fit$V_m), tolerance = 1e-10)
  expect_true(all(is.finite(fit$V_m)))
})

test_that("race v_hat / s_v carriers are length N and finite", {
  inp <- make_race_input()
  fit <- gp_estimation_core(inp, characteristic = "race")
  expect_length(fit$v_hat, length(inp$theta_hat))
  expect_length(fit$s_v, length(inp$theta_hat))
  expect_true(all(is.finite(fit$v_hat)))
  expect_true(all(is.finite(fit$s_v)))
})

test_that("race fit carries a provenance stamp with n_moments == 4", {
  fit <- gp_estimation_core(make_race_input(), characteristic = "race")
  expect_equal(fit$provenance$step, "w-seam:beta-GMM")
  expect_equal(fit$provenance$n_moments, 4L)
  expect_equal(fit$provenance$characteristic, "race")
})

test_that("gp_estimation_core requires an explicit characteristic", {
  expect_error(
    gp_estimation_core(make_race_input()),
    class = "gradepath_validation_error"
  )
})

# --- two-step driver: GENDER ------------------------------------------------

test_that("two-step gender fit runs and returns a finite beta", {
  fit <- gp_estimation_core(make_gender_input(), characteristic = "gender")
  expect_true(is.finite(fit$beta))
  expect_equal(fit$characteristic, "gender")
})

test_that("gender fit is deterministic across repeated runs (seeded multistart)", {
  f1 <- gp_estimation_core(make_gender_input(), characteristic = "gender")
  f2 <- gp_estimation_core(make_gender_input(), characteristic = "gender")
  expect_equal(f1$beta, f2$beta, tolerance = 1e-10)
})

test_that("gender J-stat finite non-negative, df == 1; SEs positive", {
  fit <- gp_estimation_core(make_gender_input(), characteristic = "gender")
  expect_true(is.finite(fit$J_stat))
  expect_gte(fit$J_stat, 0)
  expect_equal(fit$df, 1L)
  expect_true(all(fit$SE > 0))
})

test_that("gender coupling: m_hat length 1 = [sigma_xi]; V_m sym finite", {
  fit <- gp_estimation_core(make_gender_input(), characteristic = "gender")
  expect_length(fit$m_hat, 1L)
  expect_true(is.finite(fit$m_hat))
  expect_equal(fit$m_hat[1], exp(fit$delta[2]), tolerance = 1e-10)  # sigma_xi
  expect_equal(dim(fit$V_m), c(1L, 1L))
  expect_equal(fit$V_m, t(fit$V_m), tolerance = 1e-10)
  expect_true(is.finite(as.numeric(fit$V_m)))
})

# --- jacobian backend -------------------------------------------------------

test_that("gp_jacobian computes a correct Jacobian (vs analytic linear map)", {
  fn <- function(x) c(2 * x[1] + x[2], x[1] - 3 * x[3], x[2] * x[3])
  J <- gp_jacobian(fn, c(1, 2, 3))
  J_ref <- rbind(c(2, 1, 0), c(1, 0, -3), c(0, 3, 2))  # at (1,2,3)
  expect_equal(unname(J), J_ref, tolerance = 1e-6)
})
