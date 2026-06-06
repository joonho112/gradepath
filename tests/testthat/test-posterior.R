# =============================================================================
# test-posterior.R  --  one-level posterior summaries (gp_posterior)
# -----------------------------------------------------------------------------
# gp_posterior_onelevel(): per-unit r-scale posterior mean/sd + central interval
# from the W-seam output, eb_posterior-shaped, cross-checked bit-for-bit against
# ebrecipe::eb_shrink (nonparametric one-level oracle, frozen invariant #8).
# =============================================================================

.gpp_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}
.krw   <- function(...) .gpp_get("gp_krw_gmm_input")(...)
.core  <- function(...) .gpp_get("gp_estimation_core")(...)
.wseam <- function(...) .gpp_get("gp_w_seam")(...)
.pw    <- function(...) .gpp_get("gp_posterior_weights")(...)
.post  <- function(...) .gpp_get("gp_posterior_onelevel")(...)
.valid <- function(...) .gpp_get("validate_gp_posterior")(...)

.build_race <- function() {
  inp <- .krw("race")
  fit <- .core(inp, characteristic = "race")
  cpl <- .wseam(inp, "race")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  r_est <- list(theta_hat = cpl$v_hat, s = cpl$s_v,
                original_s = inp$s, id = inp$unit_id)
  parts <- .pw(prior, r_est, fit)
  list(inp = inp, fit = fit, prior = prior, estimates = estimates,
       r_est = r_est, parts = parts)
}

test_that("gp_posterior_onelevel returns a valid r-scale gp_posterior", {
  skip_if_not_installed("ebrecipe")
  b <- .build_race()
  post <- .post(b$parts, b$prior, b$r_est, b$fit)
  expect_s3_class(post, "gp_posterior")
  expect_silent(.valid(post))
  expect_identical(post$scale, "r")
  J <- length(b$r_est$theta_hat)
  for (nm in c("estimate", "se", "posterior_mean", "posterior_sd",
               "lower", "upper", "id", "label")) {
    expect_length(post[[nm]], J)
  }
  expect_true(all(is.finite(post$posterior_mean)))
  expect_true(all(is.finite(post$posterior_sd)) && all(post$posterior_sd >= 0))
})

test_that("posterior mean equals W %*% support (consumes W, no refit)", {
  skip_if_not_installed("ebrecipe")
  b <- .build_race()
  post <- .post(b$parts, b$prior, b$r_est, b$fit)
  pm_recompute <- as.numeric(b$parts$W %*% b$prior$support)
  expect_equal(post$posterior_mean, pm_recompute, tolerance = 0)
  # sd from the per-unit second moment
  e2 <- as.numeric(b$parts$W %*% (b$prior$support^2))
  psd_recompute <- sqrt(pmax(e2 - pm_recompute^2, 0))
  expect_equal(post$posterior_sd, psd_recompute, tolerance = 0)
})

test_that("central interval brackets the mean and CDF is monotone", {
  skip_if_not_installed("ebrecipe")
  b <- .build_race()
  post <- .post(b$parts, b$prior, b$r_est, b$fit, interval_level = 0.90)
  expect_true(all(post$lower <= post$posterior_mean + 1e-8))
  expect_true(all(post$posterior_mean <= post$upper + 1e-8))
  expect_true(all(post$lower <= post$upper))
  # per-unit CDF is non-decreasing and ends at 1
  for (j in seq_len(nrow(b$parts$W))) {
    cdf <- cumsum(b$parts$W[j, ])
    expect_true(all(diff(cdf) >= -1e-12))
    expect_equal(cdf[length(cdf)], 1, tolerance = 1e-10)
  }
})

test_that("invariant #8: pm_r matches ebrecipe::eb_shrink nonparametric (one-level)", {
  skip_if_not_installed("ebrecipe")
  b <- .build_race()
  post <- .post(b$parts, b$prior, b$r_est, b$fit)
  shr <- ebrecipe::eb_shrink(b$estimates, b$prior,
                             method = "nonparametric", unstandardize = FALSE)
  eb_pm <- shr$posterior$.posterior_mean
  expect_equal(post$posterior_mean, eb_pm, tolerance = 1e-8)
})

test_that("reporting-scale summaries ride in metadata (theta axis, W4)", {
  skip_if_not_installed("ebrecipe")
  b <- .build_race()
  post <- .post(b$parts, b$prior, b$r_est, b$fit)
  rep <- post$metadata$reporting
  expect_false(is.null(rep))
  J <- length(b$r_est$theta_hat)
  expect_length(rep$posterior_mean, J)
  expect_identical(rep$scale, "theta")
  # theta-scale mean = sum_m W[j,m] * reporting_support[m, j]
  rs <- b$parts$reporting_support
  pm_t <- vapply(seq_len(J), function(j) sum(b$parts$W[j, ] * rs[, j]), numeric(1))
  expect_equal(rep$posterior_mean, pm_t, tolerance = 0)
  # the r-scale object's own posterior_mean is NOT the theta-scale one
  expect_false(isTRUE(all.equal(post$posterior_mean, rep$posterior_mean)))
})

test_that("reporting metadata is absent when no reporting_support is supplied", {
  skip_if_not_installed("ebrecipe")
  b <- .build_race()
  # feed a bare W matrix + prior (no reporting_support)
  post <- .post(b$parts$W, b$prior, b$r_est, b$fit)
  expect_true(is.null(post$metadata$reporting))
  expect_identical(post$scale, "r")
  expect_length(post$posterior_mean, length(b$r_est$theta_hat))
})

test_that("gender posterior is valid and matches eb_shrink (one-level)", {
  skip_if_not_installed("ebrecipe")
  inp <- .krw("gender")
  fit <- .core(inp, characteristic = "gender")
  cpl <- .wseam(inp, "gender")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  r_est <- list(theta_hat = cpl$v_hat, s = cpl$s_v,
                original_s = inp$s, id = inp$unit_id)
  parts <- .pw(prior, r_est, fit)
  post <- .post(parts, prior, r_est, fit)
  expect_silent(.valid(post))
  shr <- ebrecipe::eb_shrink(estimates, prior,
                             method = "nonparametric", unstandardize = FALSE)
  expect_equal(post$posterior_mean, shr$posterior$.posterior_mean,
               tolerance = 1e-8)
})
