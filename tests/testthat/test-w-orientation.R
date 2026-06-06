# =============================================================================
# test-w-orientation.R  --  W orientation + reporting-metadata guards
# -----------------------------------------------------------------------------
# Hazard 1: W is J x M, units = ROWS; reject the
# v1 M x J transpose. Hazard 2: the per-unit reporting
# back-transform reads KRW beta-GMM metadata from the FIT [precision_fit], plus
# r_estimates$original_s, NEVER prior$spline_info). The
# M x J per-unit theta axis: race multiplicative original_s^beta * support;
# gender additive mu + original_s^beta * support).
#
# These are REGRESSION guards: the production behaviour is already implemented in
# R/posterior-weights.R (gp_posterior_weights + .gp_reporting_support +
# .gp_assert_w_row_stochastic). These regression guards pin each guarantee so a
# future transpose, scale slip, metadata-source swap, or back-transform drift fails
# loudly in CI.
# =============================================================================

# Resolve internals namespace-first (works under devtools::load_all + installed).
.gpo_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}
.assert_w   <- function(...) .gpo_get(".gp_assert_w_row_stochastic")(...)
.krw        <- function(...) .gpo_get("gp_krw_gmm_input")(...)
.core       <- function(...) .gpo_get("gp_estimation_core")(...)
.wseam      <- function(...) .gpo_get("gp_w_seam")(...)
.pw         <- function(...) .gpo_get("gp_posterior_weights")(...)

# ---------------------------------------------------------------------------
# Hazard 1 -- orientation guard: W is J x M (units = rows); reject the transpose
# ---------------------------------------------------------------------------

test_that(".gp_assert_w_row_stochastic accepts a correct J x M row-stochastic W", {
  J <- 5L; M <- 8L
  W <- matrix(1 / M, nrow = J, ncol = M)        # rows sum to 1
  expect_silent(.assert_w(W, n_units = J, n_support = M))
})

test_that(".gp_assert_w_row_stochastic rejects the v1 M x J transpose", {
  J <- 5L; M <- 8L
  Wt <- matrix(1 / M, nrow = M, ncol = J)       # transposed: M x J
  # Asked for J x M but handed M x J -> located dimension error.
  err <- tryCatch(.assert_w(Wt, n_units = J, n_support = M),
                  error = function(e) conditionMessage(e))
  expect_true(is.character(err))
  expect_match(err, "J x M|units|expected", ignore.case = TRUE)
})

test_that(".gp_assert_w_row_stochastic rejects a non-row-stochastic matrix", {
  J <- 4L; M <- 6L
  W <- matrix(1, nrow = J, ncol = M)            # rows sum to M, not 1
  expect_error(.assert_w(W, n_units = J, n_support = M))
})

test_that("a square W cannot hide a transpose (rows vs cols both checked)", {
  # When J == M a transpose is shape-invisible; the row-sum check must still
  # catch a column-stochastic matrix fed as if row-stochastic.
  n <- 7L
  col_stoch <- matrix(0, n, n)
  col_stoch[1, ] <- 1                            # each COLUMN sums to 1, rows do not
  expect_error(.assert_w(col_stoch, n_units = n, n_support = n))
})

# ---------------------------------------------------------------------------
# Hazard 2 + reporting -- the per-unit reporting support (M x J), W4-safe
# ---------------------------------------------------------------------------

test_that("race reporting_support is M x J and the exact multiplicative back-transform", {
  skip_if_not_installed("ebrecipe")
  inp <- .krw("race")
  fit <- .core(inp, characteristic = "race")
  cpl <- .wseam(inp, "race")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  r_est <- list(theta_hat = cpl$v_hat, s = cpl$s_v, original_s = inp$s)

  parts <- .pw(prior, r_est, fit)
  rs <- parts$reporting_support
  M <- length(prior$support); J <- length(cpl$v_hat)
  expect_identical(dim(rs), c(M, J))             # M x J: column j = unit j's axis
  expect_true(all(is.finite(rs)))
  # EXACT race back-transform: reporting_support[m, j] = original_s[j]^beta * support[m]
  expect_equal(rs, outer(prior$support, inp$s^fit$beta, `*`), tolerance = 0)
})

test_that("gender reporting_support is the exact additive back-transform (mu shift)", {
  skip_if_not_installed("ebrecipe")
  inp <- .krw("gender")
  fit <- .core(inp, characteristic = "gender")
  cpl <- .wseam(inp, "gender")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  r_est <- list(theta_hat = cpl$v_hat, s = cpl$s_v, original_s = inp$s)

  rs <- .pw(prior, r_est, fit)$reporting_support
  M <- length(prior$support); J <- length(cpl$v_hat)
  expect_identical(dim(rs), c(M, J))
  # additive: mu + original_s[j]^beta * support[m]
  expect_equal(rs, fit$mu + outer(prior$support, inp$s^fit$beta, `*`), tolerance = 0)
})

test_that("classed gp_estimation_fit carries enough metadata for gender reporting_support", {
  skip_if_not_installed("ebrecipe")
  inp <- .krw("gender")
  fit <- .core(inp, characteristic = "gender")
  cpl <- .wseam(inp, "gender")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  r_est <- list(theta_hat = cpl$v_hat, s = cpl$s_v, original_s = inp$s)

  expect_s3_class(cpl, "gp_estimation_fit")
  expect_equal(cpl$provenance$mu, fit$mu, tolerance = 0)
  rs <- .pw(prior, r_est, cpl)$reporting_support
  expect_identical(dim(rs), c(length(prior$support), length(cpl$v_hat)))
  expect_equal(rs, fit$mu + outer(prior$support, inp$s^fit$beta, `*`), tolerance = 0)
})

test_that("reporting support uses original_s (the GMM input s), not a recompute", {
  skip_if_not_installed("ebrecipe")
  inp <- .krw("race")
  fit <- .core(inp, characteristic = "race")
  cpl <- .wseam(inp, "race")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)

  # Feed a DISTINCT original_s; reporting_support must move with it (proving it is
  # the source, not r_estimates$s = s_v or any recomputed scale).
  bumped <- inp$s * 1.5
  r_est <- list(theta_hat = cpl$v_hat, s = cpl$s_v, original_s = bumped)
  rs <- .pw(prior, r_est, fit)$reporting_support
  expect_equal(rs, outer(prior$support, bumped^fit$beta, `*`), tolerance = 0)
  # and it is NOT what s_v (the likelihood scale) would give
  rs_wrong <- outer(prior$support, cpl$s_v^fit$beta, `*`)
  expect_false(isTRUE(all.equal(rs, rs_wrong)))
})

test_that("psitrap: prior$spline_info does NOT influence reporting_support", {
  skip_if_not_installed("ebrecipe")
  inp <- .krw("race")
  fit <- .core(inp, characteristic = "race")
  cpl <- .wseam(inp, "race")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  r_est <- list(theta_hat = cpl$v_hat, s = cpl$s_v, original_s = inp$s)

  base_rs <- .pw(prior, r_est, fit)$reporting_support
  prior_bogus <- prior
  prior_bogus$spline_info <- list(psi_1 = 99, psi_2 = 99)   # poisoned metadata
  bogus_rs <- .pw(prior_bogus, r_est, fit)$reporting_support
  expect_equal(bogus_rs, base_rs, tolerance = 0)            # unaffected (W4 / Hazard 2)
})

# ---------------------------------------------------------------------------
# NULL-metadata omission -- reporting support only when the metadata is present
# ---------------------------------------------------------------------------

test_that("reporting_support is omitted when precision_fit is NULL", {
  skip_if_not_installed("ebrecipe")
  inp <- .krw("race")
  cpl <- .wseam(inp, "race")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  parts <- .pw(prior, list(theta_hat = cpl$v_hat, s = cpl$s_v), precision_fit = NULL)
  expect_null(parts$reporting_support)
  # the r-scale W is still fully usable on its own
  expect_identical(dim(parts$W), c(length(cpl$v_hat), length(prior$support)))
})

test_that("reporting_support is omitted when original_s is absent", {
  skip_if_not_installed("ebrecipe")
  inp <- .krw("race")
  fit <- .core(inp, characteristic = "race")
  cpl <- .wseam(inp, "race")
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  parts <- .pw(prior, list(theta_hat = cpl$v_hat, s = cpl$s_v), fit)  # no original_s
  expect_null(parts$reporting_support)
})
