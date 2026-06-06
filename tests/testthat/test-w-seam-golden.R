# =============================================================================
# test-w-seam-golden.R  --  GP-W-EXACT: the posterior-weight golden master
# -----------------------------------------------------------------------------
# Chapter 8 sec-ch08-golden / GP-DEC-08-B. Pins gp_posterior_weights()'s native
# recompute of W bit-for-bit against ebrecipe's internal .eb_posterior_weights on
# the KRW race/multiplicative one-level W (97 x 1000):
#
#     GP-W-EXACT : max(abs(W_gp - W_eb)) == 0   (tolerance class: exact)
#
# The ebrecipe internal is touched ONLY here, through the seam wrapper
# .gp_eb_posterior_weights() (getFromNamespace, never `:::`), under a
# skip_if_not_installed guard -- never in R/ (GP-DEC-08-B).
#
# Uses the REAL gradepath API: gp_estimation_core(input, characteristic=) returns
# a flat list (beta, mu, v_hat, s_v, characteristic, ...); gp_krw_gmm_input()
# supplies the KRW GMM input whose $s is the pre-standardization original_s.
# =============================================================================

test_that("GP-W-EXACT: recomputed W is bit-exact vs the ebrecipe internal (race)", {
  skip_if_not_installed("ebrecipe")

  inp <- gp_krw_gmm_input("race")
  fit <- gp_estimation_core(inp, characteristic = "race")
  # r-scale estimates: v_hat / s_v are the standardized carriers; original_s is
  # the pre-standardization SE = the GMM input s.
  r_estimates <- list(theta_hat = fit$v_hat, s = fit$s_v, original_s = inp$s)

  estimates <- ebrecipe::eb_input(theta_hat = fit$v_hat, s = fit$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  expect_identical(prior$scale, "r")           # W3: deconvolution is on the r-scale

  J <- length(fit$v_hat)
  M <- length(prior$support)
  expect_identical(J, 97L)

  parts <- gp_posterior_weights(prior, r_estimates, fit)
  W_gp  <- parts$W
  W_eb  <- .gp_eb_posterior_weights(estimates, prior)   # seam wrapper, test-only

  # ---- GP-W-EXACT: max|diff| == 0, exact ------------------------------------
  expect_equal(max(abs(W_gp - W_eb)), 0, tolerance = 0)
  expect_true(identical(W_gp, W_eb))           # byte-identical doubles

  # ---- W1 orientation: J x M, units = ROWS ----------------------------------
  expect_identical(dim(W_gp), dim(W_eb))
  expect_identical(dim(W_gp), c(J, M))         # 97 x 1000 (inv #5 / W1)
  expect_identical(nrow(W_gp), 97L)
  expect_false(identical(dim(W_gp), c(M, J)))  # not the v1 transpose

  # ---- W2 row-stochastic; W5 non-negative finite ----------------------------
  expect_true(all(abs(rowSums(W_gp) - 1) < 1e-12))
  expect_false(isTRUE(all(abs(colSums(W_gp) - 1) < 1e-12)))
  expect_true(all(W_gp >= 0))
  expect_true(all(is.finite(W_gp)))

  # ---- support is the r-scale grid ------------------------------------------
  expect_identical(parts$support, prior$support)

  # ---- reporting support: M x J, race multiplicative back-transform (W4) -----
  expect_true(!is.null(parts$reporting_support))
  expect_identical(dim(parts$reporting_support), c(M, J))
  expect_true(all(is.finite(parts$reporting_support)))
  expected_rs <- outer(prior$support, inp$s^fit$beta, `*`)
  expect_equal(parts$reporting_support, expected_rs, tolerance = 0)
})

test_that("gp_posterior_weights matches the oracle on the native gp_prior's public fields", {
  skip_if_not_installed("ebrecipe")

  inp <- gp_krw_gmm_input("race")
  # gp_deconvolve() requires the CLASSED coupling object (gp_estimation_fit),
  # which gp_w_seam() returns; the raw gp_estimation_core() output is a plain list.
  cpl <- gp_w_seam(inp, "race")
  expect_s3_class(cpl, "gp_estimation_fit")
  r_estimates <- list(theta_hat = cpl$v_hat, s = cpl$s_v, original_s = inp$s)
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)

  gp_prior <- gp_deconvolve(cpl)               # gradepath native prior
  expect_identical(gp_prior$scale, "r")

  parts <- gp_posterior_weights(gp_prior, r_estimates, cpl)

  # Build a valid eb_prior, then overwrite support/density with the native
  # prior's so the oracle reads identical public fields (apples-to-apples).
  eb_prior <- ebrecipe::eb_deconvolve(estimates)
  eb_prior$support <- gp_prior$support
  eb_prior$density <- gp_prior$density
  W_eb <- .gp_eb_posterior_weights(estimates, eb_prior)

  expect_identical(dim(parts$W), c(97L, length(gp_prior$support)))
  expect_equal(max(abs(parts$W - W_eb)), 0, tolerance = 0)
  expect_true(all(abs(rowSums(parts$W) - 1) < 1e-12))
})

test_that("gp_posterior_weights enforces the W contract (orientation + scale guards)", {
  support <- seq(-3, 3, length.out = 40)
  density <- dnorm(support, 0, 1)
  prior   <- list(support = support, density = density, scale = "r")
  r_est   <- list(theta_hat = c(-1, 0, 0.5, 1.2), s = c(0.5, 0.4, 0.6, 0.5))

  parts <- gp_posterior_weights(prior, r_est, precision_fit = NULL)
  expect_identical(dim(parts$W), c(4L, 40L))            # W1: J x M
  expect_true(all(abs(rowSums(parts$W) - 1) < 1e-12))   # W2
  expect_true(all(parts$W >= 0) && all(is.finite(parts$W)))  # W5
  expect_null(parts$reporting_support)                  # no metadata -> omitted

  bad_prior <- prior; bad_prior$scale <- "theta"        # W3: non-'r' rejected
  expect_error(gp_posterior_weights(bad_prior, r_est, NULL))

  bad_est <- list(theta_hat = c(0, 1), s = c(0.5))       # length mismatch rejected
  expect_error(gp_posterior_weights(prior, bad_est, NULL))
})
