# =============================================================================
# test-deconvolution.R  --  native one-level deconvolution tests
# -----------------------------------------------------------------------------
# Structural assertions on the gp_prior produced by gp_deconvolve() PLUS the
# frozen-invariant-#8 ONE-LEVEL cross-check against ebrecipe::eb_deconvolve on
# the SAME one-level carriers. Registry rows that ACTUALLY exist (verified
# against gp_registry$id) are gated; non-existent rows are NOT invented.
#
# These tests assume gp_deconvolve() (R/deconvolution.R) and the
# gp_w_seam()/gp_estimation_fit are loaded into the package namespace.
# =============================================================================

## ---- shared fixtures --------------------------------------------------------

# Build a KRW estimation input via the package's own loader. `gp_krw_gmm_input()`
# reads the shipped Matlab CSV (inst/extdata/krw-gmm-input/) and returns the
# `eb_estimates` container the W-seam core consumes. (NB: the ebrecipe input
# wrapper takes `theta_hat`/`s`, NOT `estimate`/`se` -- see seam-ebrecipe.R.)
.gp_decon_test_input <- function(characteristic) {
  gp_krw_gmm_input(characteristic)
}

# A fast control with a SMALL penalty grid (full 200/40-node grids are accurate
# but slow; the structural contract is grid-size-independent). The small grid
# brackets the Matlab optimum (race c_xi ~ 0.067; gender c_xi ~ 0.0037).
.gp_decon_test_control <- function(characteristic) {
  ctl <- tryCatch(gp_control(), error = function(e) NULL)
  if (is.null(ctl)) return(NULL)
  if (identical(characteristic, "race")) {
    ctl$deconv_penalty_grid_race <- seq(0.02, 0.12, by = 0.02)
  } else {
    ctl$deconv_penalty_grid_gender <- seq(0.001, 0.006, by = 0.001)
  }
  ctl
}

.gp_decon_test_prior <- function(characteristic) {
  fit <- gp_w_seam(.gp_decon_test_input(characteristic),
                   characteristic = characteristic)
  gp_deconvolve(fit, control = .gp_decon_test_control(characteristic))
}

.gp_decon_matlab_xi_caps <- function(fit, characteristic) {
  v <- fit$v_hat
  if (identical(characteristic, "race")) {
    mu <- fit$m_hat[1L]
    sig <- fit$m_hat[2L]
    return(c(lo = 0, hi = min(max(max(v), mu + 5 * sig), mu + 7 * sig)))
  }
  sig <- fit$m_hat[1L]
  c(lo = max(min(min(v), -5 * sig), -7 * sig),
    hi = min(max(max(v),  5 * sig),  7 * sig))
}

# Registry lookup by id (gp_registry is keyed by `id`, value `paper_value`,
# absolute tol `tolerance`).
.gp_decon_reg <- function(id) {
  e <- new.env()
  utils::data("gp_registry", package = "gradepath", envir = e)
  reg <- get("gp_registry", envir = e)
  row <- reg[reg$id == id, , drop = FALSE]
  if (nrow(row) != 1L) return(NULL)
  list(paper_value = as.numeric(row$paper_value),
       tolerance = as.numeric(row$tolerance))
}

## ---- structural -------------------------------------------------------------

test_that("gp_deconvolve returns a valid gp_prior (race)", {
  skip_if_not_installed("ebrecipe")
  prior <- .gp_decon_test_prior("race")

  expect_s3_class(prior, "gp_prior")
  expect_silent(validate_gp_prior(prior))
  expect_setequal(names(prior),
                  c("support", "density", "mean", "scale",
                    "diagnostics", "metadata"))
  expect_identical(prior$scale, "r")
})

test_that("density is nonnegative and sums to 1 on the grid (race)", {
  skip_if_not_installed("ebrecipe")
  prior <- .gp_decon_test_prior("race")
  expect_true(all(prior$density >= 0))
  expect_equal(sum(prior$density), 1, tolerance = 1e-8)
  expect_true(all(is.finite(prior$density)))
})

test_that("support spans the coupling caps with 1000 one-level points (race)", {
  skip_if_not_installed("ebrecipe")
  fit <- gp_w_seam(.gp_decon_test_input("race"), characteristic = "race")
  prior <- gp_deconvolve(fit, control = .gp_decon_test_control("race"))
  expect_length(prior$support, 1000L)
  expect_equal(min(prior$support), unname(fit$caps[["lo"]]), tolerance = 1e-10)
  expect_equal(max(prior$support), unname(fit$caps[["hi"]]), tolerance = 1e-10)
  ## race one-level support is one-sided, floored at 0.
  expect_equal(min(prior$support), 0, tolerance = 1e-10)
})

test_that("prior-implied moments match the coupling m_hat (race)", {
  skip_if_not_installed("ebrecipe")
  fit <- gp_w_seam(.gp_decon_test_input("race"), characteristic = "race")
  prior <- gp_deconvolve(fit, control = .gp_decon_test_control("race"))
  mom <- prior$diagnostics$model_moments
  ## race m_hat = c(mu, sigma_xi); the prior matches (mean, sd). The penalty is
  ## CHOSEN to minimise the V_m-weighted distance, so the match is close but not
  ## exact (regularised). Loose tolerance on the moments themselves.
  expect_equal(unname(mom[["mean"]]), fit$m_hat[1L], tolerance = 0.15)
  expect_equal(unname(mom[["sd"]]),   fit$m_hat[2L], tolerance = 0.15)
  ## prior$mean is the grid mean of the density and must equal mom["mean"].
  expect_equal(prior$mean, unname(mom[["mean"]]), tolerance = 1e-8)
})

test_that("selected penalty + J are finite and on the grid (race)", {
  skip_if_not_installed("ebrecipe")
  prior <- .gp_decon_test_prior("race")
  d <- prior$diagnostics
  expect_true(is.finite(d$penalty))
  expect_true(d$penalty %in% d$penalty_grid)
  expect_true(is.finite(d$J) && d$J >= 0)
  ## the selected penalty minimises J over the grid.
  expect_equal(d$J, min(d$J_grid[is.finite(d$J_grid)]), tolerance = 1e-10)
  expect_identical(d$penalty, d$penalty_grid[d$penalty_index])
})

test_that("optimizer convergence diagnostics are exposed for the selected penalty and grid", {
  skip_if_not_installed("ebrecipe")
  prior <- .gp_decon_test_prior("race")
  d <- prior$diagnostics
  expect_true("convergence_code" %in% names(d))
  expect_true("convergence_message" %in% names(d))
  expect_true("convergence_code_grid" %in% names(d))
  expect_true("convergence_message_grid" %in% names(d))
  expect_length(d$convergence_code_grid, length(d$penalty_grid))
  expect_length(d$convergence_message_grid, length(d$penalty_grid))
  expect_identical(d$converged, d$converged_grid[d$penalty_index])
  expect_identical(d$convergence_code, d$convergence_code_grid[d$penalty_index])
})

test_that("gp_deconvolve returns a valid gp_prior (gender)", {
  skip_if_not_installed("ebrecipe")
  prior <- .gp_decon_test_prior("gender")
  expect_s3_class(prior, "gp_prior")
  expect_silent(validate_gp_prior(prior))
  expect_equal(sum(prior$density), 1, tolerance = 1e-8)
  expect_true(all(prior$density >= 0))
  expect_length(prior$support, 1000L)
})

test_that("gender prior support uses standardized xi caps and mean ~ 0", {
  skip_if_not_installed("ebrecipe")
  fit <- gp_w_seam(.gp_decon_test_input("gender"), characteristic = "gender")
  prior <- gp_deconvolve(fit, control = .gp_decon_test_control("gender"))
  ## Matlab sets xi_hat = v_hat and, in the one-level gender branch, centers the
  ## xi support at zero. Raw theta_hat and the additive location mu are not the
  ## cap scale here.
  ref <- .gp_decon_matlab_xi_caps(fit, "gender")
  expect_equal(min(prior$support), unname(ref[["lo"]]), tolerance = 1e-10)
  expect_equal(max(prior$support), unname(ref[["hi"]]), tolerance = 1e-10)
  expect_equal(ref[["lo"]], -ref[["hi"]], tolerance = 1e-10)
  ## the mean-constraint fixes the xi mixing-density mean at 0 (likelihood.m
  ## fsolve mean_xi == 0); this holds regardless of the raw additive mu.
  expect_equal(prior$mean, 0, tolerance = 1e-3)
})

test_that("gender prior sd matches the coupling sigma_xi (m_hat length 1)", {
  skip_if_not_installed("ebrecipe")
  fit <- gp_w_seam(.gp_decon_test_input("gender"), characteristic = "gender")
  prior <- gp_deconvolve(fit, control = .gp_decon_test_control("gender"))
  expect_length(fit$m_hat, 1L)
  expect_equal(unname(prior$diagnostics$model_moments[["sd"]]),
               fit$m_hat[1L], tolerance = 0.1)
})

test_that("two-level fit is refused (invariant #8: one-level only)", {
  skip_if_not_installed("ebrecipe")
  fit <- gp_w_seam(.gp_decon_test_input("race"), characteristic = "race")
  ## fabricate a two-level-shaped m_hat (length 3) and expect a refusal.
  bad <- fit
  bad$m_hat <- c(fit$m_hat, 0.1)
  bad$V_m <- diag(3)
  expect_error(gp_deconvolve(bad), class = "gradepath_validation_error")
})

## ---- PARITY / cross-check (frozen invariant #8) -----------------------------

test_that("native one-level prior cross-checks vs ebrecipe::eb_deconvolve", {
  skip_if_not_installed("ebrecipe")
  fit <- gp_w_seam(.gp_decon_test_input("race"), characteristic = "race")
  prior <- gp_deconvolve(fit, control = .gp_decon_test_control("race"))

  ## ONE-LEVEL oracle: eb_deconvolve on the SAME carriers (v_hat, s_v). The
  ## ebrecipe input wrapper takes theta_hat / s (NOT estimate / se); reach it
  ## through the .gp_eb_input seam wrapper.
  eb_in <- get(".gp_eb_input", envir = asNamespace("gradepath"))
  oracle <- ebrecipe::eb_deconvolve(eb_in(theta_hat = fit$v_hat, s = fit$s_v))

  expect_s3_class(oracle, "eb_prior")
  expect_true(all(oracle$density >= 0))
  ## NB: the eb_prior `density` is a per-point PDF, NOT probability masses -- it
  ## does NOT sum to 1 (verified live: sum ~ 478), and its grid spacing/normal-
  ## isation is an ebrecipe-internal detail (a naive sum(density)*diff(supp)[1]
  ## gives ~1.34, not 1). The NATIVE gp_prior density IS probability masses
  ## summing to 1 (Matlab g/sum(g) convention) -- a deliberate, correct contract
  ## difference. We therefore renormalise the oracle to masses before comparing.

  ## The two priors use DIFFERENT support grids (native = GMM caps;
  ## oracle = data-driven). The canonical inv#8 cross-check is the POINTWISE max
  ## abs diff after putting BOTH on the native grid as probability masses (renorm
  ## to sum 1) and comparing. Verified live: race ~0.001, gender ~0.005.
  oi <- stats::approx(oracle$support, oracle$density,
                      xout = prior$support, rule = 2, ties = mean)$y
  oi[!is.finite(oi)] <- 0
  if (sum(oi) > 0) oi <- oi / sum(oi)
  max_abs_diff <- max(abs(prior$density - oi))
  expect_true(is.finite(max_abs_diff))
  expect_lt(max_abs_diff, 0.05)           # one-level priors agree closely

  ## Secondary: the first two MOMENTS are close (loose -- native is GMM-coupled,
  ## oracle is data-driven). Verified live: native (mean 0.323, sd 0.221) vs
  ## oracle (mean 0.316, sd 0.215) for race.
  od <- oracle$density / sum(oracle$density)
  omean <- sum(oracle$support * od)
  osd   <- sqrt(sum((oracle$support - omean)^2 * od))
  expect_equal(prior$mean, omean, tolerance = 0.1)
  expect_equal(prior$diagnostics$model_moments[["sd"]], osd, tolerance = 0.1)
})

## ---- registry gating (ONLY rows that ACTUALLY exist) ------------------------
# VERIFIED against gp_registry$id (53 rows; keyed by `id`, value `paper_value`,
# abs tol `tolerance`). There is NO `race_sigma_xi` / `gender_sigma_xi` /
# `*_E_theta` / `*_SD_theta` row -- earlier fabricated reads invented those, so
# we do NOT gate them. The dispersion-related ids that DO exist are all
# TWO-LEVEL ("wi" = with-industry), hence OUT of one-level scope (#8):
#   t3_race_wi_sigmaxi (0.113), t3_race_wi_sigmaeta (0.528),
#   t3b_race_withinshare (0.366), t3b_gender_withinshare (0.562).
# The one-level ("ni" = no-industry) row is `t3_race_ni_beta` (0.510, tol
# 0.001) -- the precision exponent the deconvolution consumes verbatim
# (invariant #3, carried on prior$metadata$beta). We gate THAT, and add a
# documented-absence test so no future fabricated read silently re-introduces a
# bogus one-level sigma_xi row.

test_that("registry: race no-industry beta (t3_race_ni_beta) matches paper", {
  skip_if_not_installed("ebrecipe")
  reg <- .gp_decon_reg("t3_race_ni_beta")
  skip_if(is.null(reg), "gp_registry row 't3_race_ni_beta' not present.")
  prior <- .gp_decon_test_prior("race")
  ## the deconvolution consumes fit$beta verbatim (invariant #3) and stamps it
  ## on the prior metadata. Gate it against the paper's one-level beta = 0.510.
  expect_equal(prior$metadata$beta, reg$paper_value, tolerance = reg$tolerance)
})

test_that("registry: no one-level sigma_xi row exists (documented absence)", {
  ## DEFENSIVE. If a future registry adds a dedicated one-level prior-dispersion
  ## row, this fails loudly and the gate above should be extended. Today there
  ## is no such row; the two-level dispersion rows DO exist (out of #8 scope).
  expect_null(.gp_decon_reg("race_sigma_xi"))
  expect_null(.gp_decon_reg("gender_sigma_xi"))
  expect_false(is.null(.gp_decon_reg("t3_race_wi_sigmaxi")))
})

## ---- CCR-04: opt-in seeded per-penalty multistart (n_starts) ----------------
# The penalty loop fits the spline density once per penalty node. The historical
# default is a SINGLE deterministic zeros-start (a determinism substitution for
# KRW's per-penalty `randn` restart, estimate_lsqnonlin.m:356). `n_starts > 1`
# is an OPT-IN seeded multistart: start #1 is STILL the zeros vector, plus
# `n_starts - 1` seeded `rnorm` starts per node, keeping the lowest-objective
# fit. Because start #1 is always zeros, the per-node objective can only improve
# or tie -- the default is never degraded.
#
# EMPIRICAL FINDING (small test grids; race seq(0.02,0.12,0.02), gender
# seq(0.001,0.006,0.001), default seed 1L):
#   * RACE: multistart lowers the per-node penalised objective at every node
#     (max ~1.4e-4), but the selected penalty is unchanged (0.06) and the
#     coupling-J is NOT improved (J actually rises ~5.4e-4: the per-node
#     log-likelihood objective and the coupling-J selection metric are DIFFERENT
#     objectives). => race is a STABILITY case (no improvement; selection robust).
#   * GENDER: with the zeros-start, 5 of 6 nodes fail to a non-finite density and
#     only c_xi = 0.006 yields a finite fit, so n_starts=1 selects 0.006 almost
#     by default. The multistart RESCUES the failed nodes (all 6 finite) and
#     finds a STRICTLY LOWER coupling-J: J drops 0.001267 -> 0.001108 and the
#     selected penalty moves 0.006 -> 0.004. => gender is a documented WIN.
# The robust, characteristic-INDEPENDENT invariant we assert is ONLY: the
# per-node PENALISED OBJECTIVE is never worse under multistart (guaranteed --
# start #1 is the same zeros vector). The selected coupling-J is a SEPARATE
# metric: it is NOT guaranteed to improve and on race it rises slightly (the
# stability case); on gender it strictly improves (the documented win). We do
# NOT assert a universal "selected J never worse" claim -- that is false (race).

test_that("n_starts validates as a positive integer", {
  skip_if_not_installed("ebrecipe")
  fit <- gp_w_seam(.gp_decon_test_input("race"), characteristic = "race")
  expect_error(gp_deconvolve(fit, n_starts = 0L),
               class = "gradepath_validation_error")
  expect_error(gp_deconvolve(fit, n_starts = -2L),
               class = "gradepath_validation_error")
  expect_error(gp_deconvolve(fit, n_starts = 2.5),
               class = "gradepath_validation_error")
  expect_error(gp_deconvolve(fit, n_starts = c(2L, 3L)),
               class = "gradepath_validation_error")
  expect_error(gp_deconvolve(fit, n_starts = NA_integer_),
               class = "gradepath_validation_error")
})

# Slow-test gate: the per-node multistart over the full penalty grid (n_starts > 1)
# is compute-heavy (race ~30s+ per call; the whole empirical suite OOM/time-limits
# in constrained CI). These EMPIRICAL multistart tests are SKIPPED by default so
# the suite stays light, and are run with GRADEPATH_RUN_SLOW_TESTS=1. The always-on tests
# (validation + the bit-identical n_starts = 1L default below) guarantee the
# parity-critical contract regardless.
.gp_decon_skip_slow <- function() {
  if (!nzchar(Sys.getenv("GRADEPATH_RUN_SLOW_TESTS"))) {
    testthat::skip("slow multistart test (set GRADEPATH_RUN_SLOW_TESTS=1 to run)")
  }
}

test_that("n_starts = 1L reproduces the default (single zeros-start) bit-for-bit", {
  skip_if_not_installed("ebrecipe")
  for (ch in c("race", "gender")) {
    ctl <- .gp_decon_test_control(ch)
    fit <- gp_w_seam(.gp_decon_test_input(ch), characteristic = ch)
    p_default <- gp_deconvolve(fit, control = ctl)              # no n_starts arg
    p_one     <- gp_deconvolve(fit, control = ctl, n_starts = 1L)
    ## bit-identical density, mean, penalty, J, alpha, full J grid.
    expect_identical(p_one$density, p_default$density)
    expect_identical(p_one$mean, p_default$mean)
    expect_identical(p_one$diagnostics$penalty, p_default$diagnostics$penalty)
    expect_identical(p_one$diagnostics$J, p_default$diagnostics$J)
    expect_identical(p_one$diagnostics$J_grid, p_default$diagnostics$J_grid)
    expect_identical(p_one$diagnostics$objective_grid,
                     p_default$diagnostics$objective_grid)
    expect_identical(p_one$metadata$alpha, p_default$metadata$alpha)
  }
  ## and the historical published-scale values are preserved by the default.
  fit_r <- gp_w_seam(.gp_decon_test_input("race"), characteristic = "race")
  pr <- gp_deconvolve(fit_r, control = .gp_decon_test_control("race"))
  expect_equal(pr$mean, 0.323, tolerance = 1e-2)   # race prior mean ~0.323
  expect_equal(sum(pr$density), 1, tolerance = 1e-10)
  fit_g <- gp_w_seam(.gp_decon_test_input("gender"), characteristic = "gender")
  pg <- gp_deconvolve(fit_g, control = .gp_decon_test_control("gender"))
  expect_equal(pg$mean, 0, tolerance = 1e-3)        # gender prior mean ~0
  expect_equal(sum(pg$density), 1, tolerance = 1e-10)
})

test_that("multistart is deterministic given the seed (n_starts > 1)", {
  skip_if_not_installed("ebrecipe")
  .gp_decon_skip_slow()
  for (ch in c("race", "gender")) {
    ctl <- .gp_decon_test_control(ch)
    fit <- gp_w_seam(.gp_decon_test_input(ch), characteristic = ch)
    ## default seed (control$seed unset -> fixed default 1L): two calls identical.
    a <- gp_deconvolve(fit, control = ctl, n_starts = 4L)
    b <- gp_deconvolve(fit, control = ctl, n_starts = 4L)
    expect_identical(a$density, b$density)
    expect_identical(a$diagnostics$J, b$diagnostics$J)
    expect_identical(a$diagnostics$objective_grid, b$diagnostics$objective_grid)
    expect_identical(a$metadata$alpha, b$metadata$alpha)
    ## explicit control$seed: reproducible, and a DIFFERENT seed changes the draws
    ## (so the multistart genuinely uses the seed, not a constant).
    ctl1 <- ctl; ctl1$seed <- 42L
    ctl2 <- ctl; ctl2$seed <- 99L
    s1  <- gp_deconvolve(fit, control = ctl1, n_starts = 4L)
    s1b <- gp_deconvolve(fit, control = ctl1, n_starts = 4L)
    s2  <- gp_deconvolve(fit, control = ctl2, n_starts = 4L)
    expect_identical(s1$diagnostics$objective_grid, s1b$diagnostics$objective_grid)
    expect_false(identical(s1$diagnostics$objective_grid,
                           s2$diagnostics$objective_grid))
  }
})

test_that("multistart never raises the per-node penalised objective (the robust invariant)", {
  skip_if_not_installed("ebrecipe")
  .gp_decon_skip_slow()
  ## The ROCK-SOLID invariant: because start #1 is ALWAYS the same zeros vector,
  ## the per-node penalised log-likelihood objective under multistart can only
  ## tie or IMPROVE -- never worsen. (Note this is the OPTIMISER objective, NOT
  ## the coupling-J selection metric: see the race-stability test below for why
  ## a lower objective does NOT imply a lower selected J.)
  for (ch in c("race", "gender")) {
    ctl <- .gp_decon_test_control(ch)
    fit <- gp_w_seam(.gp_decon_test_input(ch), characteristic = ch)
    p1 <- gp_deconvolve(fit, control = ctl, n_starts = 1L)
    p8 <- gp_deconvolve(fit, control = ctl, n_starts = 4L)
    o1 <- p1$diagnostics$objective_grid
    o8 <- p8$diagnostics$objective_grid
    fin <- is.finite(o1) & is.finite(o8)
    expect_true(any(fin))
    expect_true(all(o8[fin] <= o1[fin] + 1e-9))
    ## multistart never has FEWER finite nodes than the single zeros-start
    ## (it can only rescue nodes that the zeros-start failed on).
    expect_true(sum(is.finite(o8)) >= sum(is.finite(o1)))
  }
})

test_that("race multistart improves per-node objectives; selected-J direction is draw-dependent", {
  skip_if_not_installed("ebrecipe")
  .gp_decon_skip_slow()
  ## EMPIRICAL REALITY (race, small grid): all 6 nodes already fit finitely from
  ## the zeros-start, and the seeded multistart finds a strictly lower per-node
  ## penalised objective at >= 1 node. The DIRECTION of the *selected* coupling-J
  ## change is NOT a robust invariant: it depends on the RNG draw structure. Two
  ## independent implementations with different per-node draw schemes but the same
  ## seed observed OPPOSITE signs
  ## for race (one J ticked up ~+5.4e-4, the other down ~-1.6e-4) -- because the
  ## optimiser minimises the penalised NEG-LOG-LIKELIHOOD while the penalty
  ## SELECTION minimises coupling-J (a different metric). We therefore assert only
  ## the robust, draw-independent facts and NOT a J direction. This draw-
  ## dependence is exactly why the default stays single-start (bit-identical) and
  ## the multistart is strictly opt-in.
  ctl <- .gp_decon_test_control("race")
  fit <- gp_w_seam(.gp_decon_test_input("race"), characteristic = "race")
  p1 <- gp_deconvolve(fit, control = ctl, n_starts = 1L)
  p8 <- gp_deconvolve(fit, control = ctl, n_starts = 4L)
  o1 <- p1$diagnostics$objective_grid; o8 <- p8$diagnostics$objective_grid
  fin <- is.finite(o1) & is.finite(o8)
  expect_true(all(o8[fin] <= o1[fin] + 1e-9))   # never worse (start #1 is the zeros start)
  expect_true(any(o1[fin] - o8[fin] > 1e-8))    # improves somewhere
  ## NB: not even the selected *penalty* is a robust invariant. On the FULL
  ## default grid race n>1 selects a DIFFERENT penalty (0.066 -> 0.067) with
  ## higher J; on a fast grid with different draws the J directions were opposite.
  ## So the only draw/grid-independent facts are the
  ## per-node-objective ones above plus a finite selected J. We assert nothing
  ## about the selected J value, its direction, or the selected penalty.
  expect_true(is.finite(p8$diagnostics$J))      # selected J finite (value/direction/penalty NOT asserted)
})

test_that("gender multistart strictly improves the selected coupling-J (documented win)", {
  skip_if_not_installed("ebrecipe")
  .gp_decon_skip_slow()
  ## On the gender small test grid the zeros-start leaves 5/6 nodes non-finite,
  ## so n_starts=1 effectively selects c_xi = 0.006 by default. The seeded
  ## multistart rescues the failed nodes and finds a strictly lower coupling-J
  ## at a different penalty (J ~0.001267 -> ~0.001108; penalty 0.006 -> 0.004).
  ctl <- .gp_decon_test_control("gender")
  fit <- gp_w_seam(.gp_decon_test_input("gender"), characteristic = "gender")
  p1 <- gp_deconvolve(fit, control = ctl, n_starts = 1L)
  p8 <- gp_deconvolve(fit, control = ctl, n_starts = 4L)
  ## n_starts=1 finds far fewer finite J nodes than the multistart.
  expect_lt(sum(is.finite(p1$diagnostics$J_grid)),
            sum(is.finite(p8$diagnostics$J_grid)))
  ## strict improvement in the selected coupling-J.
  expect_lt(p8$diagnostics$J, p1$diagnostics$J)
  ## the gender mean-constraint (mean_xi == 0) is preserved under multistart.
  expect_equal(p8$mean, 0, tolerance = 1e-3)
  expect_equal(sum(p8$density), 1, tolerance = 1e-10)
})

test_that("multistart leaves the global RNG state untouched", {
  skip_if_not_installed("ebrecipe")
  .gp_decon_skip_slow()
  ctl <- .gp_decon_test_control("race"); ctl$seed <- 7L
  fit <- gp_w_seam(.gp_decon_test_input("race"), characteristic = "race")
  set.seed(12345L)
  before <- .Random.seed
  invisible(gp_deconvolve(fit, control = ctl, n_starts = 4L))
  expect_identical(.Random.seed, before)
})
