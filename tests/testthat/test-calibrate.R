# =============================================================================
# test-calibrate.R -- gp_calibrate() / the gp_calibration object
#                     (GP-DEC-18-A).
# -----------------------------------------------------------------------------
# The whole file is gated on an open solver (highs): the seed fit and every
# refit run the full one-level grade IP. The seed fit (the deconvolution penalty
# grid) is the dominant cost, so it is built ONCE and memoized; the deterministic
# machinery tests then mock the heavy per-draw refit (`.gp_cal_refit`) to return
# that real fit, exercising the full oracle / coverage / regret / aggregation
# machinery on a genuine `gp_fit` WITHOUT paying for (or risking the small-N
# fragility of) a fresh per-draw beta-GMM. A separate, tolerant smoke test runs
# the REAL unmocked refit path (n_sim = 2) and accepts either a valid object or
# the informative abort -- the KRW beta-GMM legitimately fails on a fraction of
# small synthetic draws, which is the behaviour the harness is built to absorb.
#
# n_sim is kept TINY (2-3); we never run the 200-draw default in the suite. The
# assertions check well-formedness and finite / in-range scalars -- NOT that the
# two boolean verdicts pass (on a tiny fit they usually do not, which is honest).
# =============================================================================

# A small, well-behaved race fit the real one-level beta-GMM accepts (a spread of
# theta_hat with heteroskedastic, non-collinear s). Built once and memoized.
gp_cal_test_fit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }
    J <- 8L
    data <- data.frame(
      theta_hat = c(-0.50, -0.30, -0.12, -0.04, 0.06, 0.18, 0.34, 0.52),
      s = c(0.21, 0.25, 0.20, 0.28, 0.19, 0.26, 0.22, 0.24),
      unit_id = paste0("u", seq_len(J)),
      label = paste0("U", seq_len(J)),
      stringsAsFactors = FALSE
    )
    cached <<- krw_report_card(
      data,
      demographic = "race",
      control = gp_control(backend = "highs", precision_rule = "krw_gmm")
    )
    cached
  }
})

# A deterministic refit mock: every "refit" returns the same precomputed real
# fit, so the downstream oracle / coverage / regret / aggregation run for real on
# a genuine gp_fit while the heavy, fragile per-draw beta-GMM is bypassed. The
# signature matches .gp_cal_refit(syn, dgp, ci_level).
gp_cal_const_refit <- function(fit) {
  force(fit)
  function(syn, dgp, ci_level) fit
}


test_that("gp_calibrate() returns a well-formed compact gp_calibration", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()

  testthat::local_mocked_bindings(.gp_cal_refit = gp_cal_const_refit(fit))
  cal <- gp_calibrate(fit, n_sim = 2L, seed = 1L)

  expect_s3_class(cal, "gp_calibration")
  for (nm in c("n_sim", "seed", "characteristic", "n_ok", "n_failed",
               "dr_mean", "dr_target", "coverage", "ci_level", "regret_mean",
               "dr_ok", "coverage_ok")) {
    expect_true(nm %in% names(cal), info = nm)
  }
  # the two boolean verdicts
  expect_type(cal$dr_ok, "logical")
  expect_type(cal$coverage_ok, "logical")
  expect_length(cal$dr_ok, 1L)
  expect_length(cal$coverage_ok, 1L)
  expect_false(is.na(cal$dr_ok))
  expect_false(is.na(cal$coverage_ok))

  # provenance scalars + the success/skip split
  expect_identical(cal$n_sim, 2L)
  expect_identical(cal$seed, 1L)
  expect_identical(cal$dr_target, 0.05)
  expect_identical(cal$ci_level, 0.90)
  expect_identical(cal$characteristic, "race")
  expect_identical(cal$n_ok + cal$n_failed, cal$n_sim)
  expect_gte(cal$n_ok, 1L)

  # dr_mean / coverage finite in [0, 1]; regret finite and >= 0 by construction.
  expect_true(is.finite(cal$dr_mean) && cal$dr_mean >= 0 && cal$dr_mean <= 1)
  expect_true(is.finite(cal$coverage) && cal$coverage >= 0 && cal$coverage <= 1)
  expect_true(is.finite(cal$regret_mean))
  expect_gte(cal$regret_mean, -1e-6)

  # COMPACT: no matrices, no per-draw archive (only scalars + verdicts + small
  # provenance list).
  expect_false(any(vapply(cal[setdiff(names(cal), "provenance")],
                          is.matrix, logical(1))))
})


test_that("gp_calibrate() is reproducible for a fixed seed", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()
  testthat::local_mocked_bindings(.gp_cal_refit = gp_cal_const_refit(fit))

  a <- gp_calibrate(fit, n_sim = 2L, seed = 1L)
  b <- gp_calibrate(fit, n_sim = 2L, seed = 1L)
  expect_identical(a, b)

  # The caller's global RNG state is left untouched.
  set.seed(99)
  before <- get(".Random.seed", envir = .GlobalEnv)
  invisible(gp_calibrate(fit, n_sim = 2L, seed = 1L))
  expect_identical(before, get(".Random.seed", envir = .GlobalEnv))

  # A different seed draws different synthetic data; the coverage (which depends
  # on each draw's true theta) generically moves.
  c2 <- gp_calibrate(fit, n_sim = 2L, seed = 7L)
  expect_s3_class(c2, "gp_calibration")
})


test_that("print.gp_calibration is a compact result-first one-liner returning invisibly", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()
  testthat::local_mocked_bindings(.gp_cal_refit = gp_cal_const_refit(fit))
  cal <- gp_calibrate(fit, n_sim = 2L, seed = 1L)

  out <- capture.output(ret <- withVisible(print(cal)))
  expect_false(ret$visible)
  expect_identical(ret$value, cal)

  expect_lte(length(out), 2L)
  joined <- paste(out, collapse = " ")
  expect_match(joined, "^<gp_calibration>")
  expect_match(joined, "DR mean")
  expect_match(joined, "coverage")
  expect_match(joined, "regret")
  # No claim / tier language; ASCII only.
  expect_false(grepl("claim", joined, ignore.case = TRUE))
  expect_false(grepl("tier", joined, ignore.case = TRUE))
  expect_false(any(grepl("[^\x01-\x7F]", out)))
})


test_that("summary.gp_calibration returns a one-row typed data frame", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()
  testthat::local_mocked_bindings(.gp_cal_refit = gp_cal_const_refit(fit))
  cal <- gp_calibrate(fit, n_sim = 2L, seed = 1L)

  s <- summary(cal)
  expect_true(is.data.frame(s))
  expect_identical(nrow(s), 1L)
  expect_true(all(c("dr_mean", "coverage", "regret_mean",
                    "dr_ok", "coverage_ok") %in% names(s)))
})


test_that("the DGP / simulate helpers expose the truth ingredients", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()

  dgp <- .gp_cal_dgp_from_fit(fit)
  expect_true(all(c("support", "density", "beta", "s",
                    "characteristic", "prior") %in% names(dgp)))
  expect_length(dgp$density, length(dgp$support))
  expect_equal(sum(dgp$density), 1, tolerance = 1e-8)
  expect_true(is.finite(dgp$beta))
  expect_true(all(dgp$s > 0))
  expect_identical(dgp$characteristic, "race")
  expect_identical(dgp$J, length(.gp_estimates_se(fit$estimates)))

  # A simulated dataset has J synthetic theta_hat plus the true theta / v of len J.
  syn <- .gp_cal_simulate_dataset(dgp)
  expect_identical(nrow(syn$data), dgp$J)
  expect_length(syn$theta_true, dgp$J)
  expect_length(syn$v_true, dgp$J)
  expect_true(all(is.finite(syn$data$theta_hat)))
})


test_that("oracle pairwise uses true-DGP precision, not refit beta", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()
  dgp <- .gp_cal_dgp_from_fit(fit)
  set.seed(123)
  syn <- .gp_cal_simulate_dataset(dgp)

  base <- .gp_cal_oracle_pairwise(fit, syn, dgp)
  expect_s3_class(base, "gp_pairwise")

  shifted <- fit
  beta0 <- shifted$precision_fit$parameters$beta %gp_or% shifted$precision_fit$beta
  if (!is.null(shifted$precision_fit$parameters$beta)) {
    shifted$precision_fit$parameters$beta <- beta0 + 0.5
  }
  if (!is.null(shifted$precision_fit$beta)) {
    shifted$precision_fit$beta <- beta0 + 0.5
  }

  perturbed <- .gp_cal_oracle_pairwise(shifted, syn, dgp)
  expect_s3_class(perturbed, "gp_pairwise")
  expect_equal(perturbed$matrix, base$matrix, tolerance = 1e-12)
})


test_that("oracle pairwise matches an explicit true-DGP construction", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()
  dgp <- .gp_cal_dgp_from_fit(fit)
  set.seed(456)
  syn <- .gp_cal_simulate_dataset(dgp)

  shifted <- fit
  beta0 <- shifted$precision_fit$parameters$beta %gp_or% shifted$precision_fit$beta
  if (!is.null(shifted$precision_fit$parameters$beta)) {
    shifted$precision_fit$parameters$beta <- beta0 + 0.5
  }
  if (!is.null(shifted$precision_fit$beta)) {
    shifted$precision_fit$beta <- beta0 + 0.5
  }

  theta_hat <- syn$data$theta_hat
  s <- syn$data$s
  v_hat <- if (identical(dgp$characteristic, "gender")) {
    (theta_hat - dgp$mu) / (s ^ dgp$beta)
  } else {
    theta_hat / (s ^ dgp$beta)
  }
  manual_r_estimates <- list(
    theta_hat = v_hat,
    s = s ^ (1 - dgp$beta),
    original_s = s,
    id = as.character(syn$data$unit_id),
    label = as.character(syn$data$label)
  )
  truth_precision <- list(
    beta = dgp$beta,
    mu = dgp$mu,
    characteristic = dgp$characteristic,
    model_form = if (identical(dgp$characteristic, "gender")) "additive" else "multiplicative"
  )

  expected <- gp_pairwise(
    gp_posterior_weights(
      prior = dgp$prior,
      r_estimates = manual_r_estimates,
      precision_fit = truth_precision
    ),
    ids = as.character(syn$data$unit_id),
    control = dgp$control
  )
  got <- .gp_cal_oracle_pairwise(shifted, syn, dgp)

  expect_s3_class(got, "gp_pairwise")
  expect_equal(got$matrix, expected$matrix, tolerance = 1e-12)
})


test_that("a refit failure is skipped, not fatal, and counted in n_failed", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()

  # Fail the FIRST refit only; the rest return the real fit. The failed draw must
  # be skipped and counted, and the run must still summarize over the survivors.
  const <- gp_cal_const_refit(fit)
  calls <- 0L
  testthat::local_mocked_bindings(
    .gp_cal_refit = function(syn, dgp, ci_level) {
      calls <<- calls + 1L
      if (calls == 1L) NULL else const(syn, dgp, ci_level)
    }
  )
  cal <- gp_calibrate(fit, n_sim = 3L, seed = 1L, min_ok = 1L)
  expect_s3_class(cal, "gp_calibration")
  expect_gte(cal$n_failed, 1L)
  expect_identical(cal$n_ok + cal$n_failed, 3L)
})


test_that("gp_calibrate() aborts informatively when every draw fails", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()

  testthat::local_mocked_bindings(.gp_cal_refit = function(syn, dgp, ci_level) NULL)
  expect_error(
    gp_calibrate(fit, n_sim = 2L, seed = 1L, min_ok = 1L),
    regexp = "succeeded",
    class = "gp_calibration_error"
  )
})


test_that("gp_calibrate() validates its scalar arguments and the fit", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()
  expect_error(gp_calibrate(list()), class = "gradepath_error")
  expect_error(gp_calibrate(fit, n_sim = 0L), class = "gp_calibration_error")
  expect_error(gp_calibrate(fit, n_sim = 2L, dr_target = 1.5), class = "gp_calibration_error")
  expect_error(gp_calibrate(fit, n_sim = 2L, ci_level = 1), class = "gp_calibration_error")
})


test_that("the REAL unmocked refit path runs (tolerant of beta-GMM fragility)", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fit <- gp_cal_test_fit()

  # No mock: exercise the genuine simulate -> refit -> oracle -> metrics path.
  # The small-N beta-GMM may fail on every draw; accept either a valid result or
  # the informative abort. This proves the real wiring runs end-to-end.
  res <- tryCatch(
    gp_calibrate(fit, n_sim = 2L, seed = 1L, min_ok = 1L),
    gp_calibration_error = function(e) e
  )
  if (inherits(res, "gp_calibration")) {
    expect_true(is.finite(res$dr_mean))
    expect_true(is.finite(res$coverage))
    expect_gte(res$regret_mean, -1e-6)
    expect_identical(res$n_ok + res$n_failed, res$n_sim)
  } else {
    expect_s3_class(res, "gp_calibration_error")
  }
})
