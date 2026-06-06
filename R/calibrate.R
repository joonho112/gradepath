# =============================================================================
# calibrate.R -- gp_calibrate(), a seeded Monte-Carlo calibration
#                harness for the one-level KRW grading method (GP-DEC-18-A).
# -----------------------------------------------------------------------------
# gp_calibrate() treats a FITTED prior as ground truth, simulates synthetic
# datasets from it, RE-RUNS the full one-level pipeline (krw_report_card) on each,
# and checks -- QUALITATIVELY -- that the grading method's own frequentist
# guarantees hold: that the empirical discordance rate (DR) tracks the target,
# that credible intervals cover the truth near the nominal level, and that the
# empirical grading's regret against an oracle that KNOWS the true prior is small
# and non-negative.
#
# It is NOT a replication of KRW Table E1, it is NOT seeded to match Matlab (the
# 15238 default is a cosmetic homage to the reference RNG seed -- in-package
# reproducibility only), and it is NOT part of any default scorecard. It is a
# standalone, exported diagnostic verb that never runs on package load or on
# vignette render.
#
# The return is a compact `gp_calibration` object: a handful of scalars plus two
# boolean verdicts. It carries NO J x J matrices and NO per-draw archive.
#
# -- The two genuine pitfalls, and how this file resolves them -----------------
#
#   (1) ORACLE grading uses the KNOWN true prior WITHOUT re-deconvolving.
#       We build the oracle pairwise `Pi_truth` directly from the truth's prior
#       density via gp_posterior_weights(prior = truth$prior, ...) -> gp_pairwise,
#       reusing the refit's standardized r-estimates and precision fit. We then
#       grade `Pi_truth` and score BOTH the empirical and the oracle grading on
#       this same truth pairwise. With risk R = (1 - lambda) * DR - lambda * tau
#       at the selected lambda, the oracle grading is the risk-minimizer ON
#       `Pi_truth`, so regret = R(empirical) - R(oracle) >= 0 by construction.
#
#   (2) COVERAGE scale. The true discrimination parameter is theta_i = s_i^beta *
#       v_i (race) or mu + s_i^beta * v_i (gender), generated on the INPUT fit's
#       r-scale (beta_truth). A refit re-estimates its OWN beta, so its r/v scale
#       is generally a DIFFERENT standardization -- comparing the true v against
#       the refit's r-scale interval would be a scale mismatch whenever
#       beta_refit != beta_truth. We therefore assess coverage on the common,
#       beta-invariant THETA scale: does the refit's theta-scale REPORTING
#       interval (already back-transformed per unit by the refit's s^beta) cover
#       the true theta? The refit is run at interval_level = ci_level so the
#       reported interval is the requested level.
#
# -- Robustness ----------------------------------------------------------------
#   The KRW beta-GMM is genuinely fragile on small synthetic draws: a fraction of
#   refits fail (singular step-1 moment covariance) or land in a degenerate
#   large-|beta| basin (s^beta -> 0) that corrupts the r-scale standardization,
#   the posterior, and the grading. gp_calibrate is robust: each draw's refit is
#   wrapped in tryCatch, failures AND degenerate-beta refits are SKIPPED and
#   counted in `n_failed`, and the reported means are over the successful draws
#   only. If fewer than `min_ok` draws succeed it errors informatively rather than
#   averaging noise.
#
# Source is ASCII-only; errors go through .gradepath_abort(); helpers are .gp_*.
# =============================================================================


# ---- small internal helpers -------------------------------------------------

# Largest plausible |beta| for a NON-degenerate one-level fit. The KRW published
# race/gender betas are 0.51 / 1.26; the small-N beta-GMM occasionally converges
# into a spurious large-|beta| basin (s^beta -> 0) that makes the r-scale -- and
# therefore the posterior, the standardized v_hat, and the grading -- meaningless.
# Such a refit is treated as a (soft) non-convergence and skipped. Generous bound.
.gp_cal_beta_sane <- 20

.gp_cal_abort <- function(message, ...) {
  .gradepath_abort(
    message,
    ...,
    class = c("gp_calibration_error", "gradepath_error")
  )
}

# A finite scalar in [lo, hi] (inclusive) or (lo, hi) (exclusive)?
.gp_cal_validate_scalar <- function(x, name, lo, hi, inclusive = TRUE) {
  ok_num <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!ok_num) {
    .gp_cal_abort("`%s` must be a single finite numeric value.", name)
  }
  in_range <- if (inclusive) (x >= lo && x <= hi) else (x > lo && x < hi)
  if (!in_range) {
    .gp_cal_abort("`%s` must lie in [%s, %s].", name, format(lo), format(hi))
  }
  as.numeric(x)
}

# Coerce a positive-integer count argument (n_sim, seed, min_ok).
.gp_cal_validate_count <- function(x, name, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      abs(x - round(x)) > 1e-8) {
    .gp_cal_abort("`%s` must be a single integer.", name)
  }
  x <- as.integer(round(x))
  if (x < min) {
    .gp_cal_abort("`%s` must be at least %d.", name, as.integer(min))
  }
  x
}

# Budget heads-up: gp_calibrate() re-solves a grade integer program once per draw,
# so its cost is ~ n_sim x one refit+grade. Returns a one-line warning string when
# a run is likely slow (large n_sim and/or a large fit), else NULL. Heuristic
# threshold; the message is informational (suppressible) and never changes results.
.gp_cal_budget_message <- function(n_sim, n_units) {
  if (!is.finite(n_units) || n_sim * n_units <= 600L) {
    return(NULL)
  }
  sprintf(paste0(
    "gp_calibrate(): about to run %d refit+grade solves on %d units. The grade ",
    "step is an integer program (often minutes per solve on a large fit), so this ",
    "can take a long time. For a quick check use a smaller n_sim or a smaller, ",
    "well-identified fit; gp_preview() gauges one solve's cost."),
    n_sim, n_units)
}

# Risk R = (1 - lambda) * DR - lambda * tau, scored for a grading on a pairwise
# matrix Pi at the selected lambda. The canonical reported risk the grade IP
# minimizes (see .gradepath_grade_objective_value / frontier.R).
.gp_cal_risk <- function(grades, pi_matrix, lambda) {
  m <- .gradepath_frontier_metrics(grades, pi_matrix)
  (1 - lambda) * m$discordance_rate - lambda * m$tau_bar
}


# ---- DGP from a fitted prior ------------------------------------------------

#' Extract the data-generating truth from a fitted one-level prior
#'
#' The fitted prior `g_hat`, the precision exponent `beta` (and additive `mu` for
#' gender), and the FIXED observed standard errors `s_i` become the "truth": the
#' simulation draws latent v_i ~ g_hat and forms theta_i = s_i^beta * v_i (race)
#' or mu + s_i^beta * v_i (gender).
#'
#' @param fit A `gp_fit` from `krw_report_card()` / `gradepath()`.
#' @return A list `gp_dgp` carrying `support`, `density`, `beta`, `mu`, `s`,
#'   `characteristic`, `prior` (the original `gp_prior` for the oracle), `ids`,
#'   `precision_fit`, and `control` (cloned from the fit).
#' @keywords internal
#' @noRd
.gp_cal_dgp_from_fit <- function(fit) {
  fit <- validate_gp_fit(fit)

  prior <- fit$prior
  support <- as.numeric(prior$support)
  density <- as.numeric(prior$density)
  if (length(support) < 2L || length(support) != length(density) ||
      any(!is.finite(support)) || any(!is.finite(density)) ||
      any(density < 0) || sum(density) <= 0) {
    .gp_cal_abort("`fit$prior` must carry a valid support/density on the r scale.")
  }
  if (!identical(prior$scale, "r")) {
    .gp_cal_abort("`fit$prior$scale` must be 'r' for the calibration DGP.")
  }

  pf <- fit$precision_fit
  pars <- if (!is.null(pf$parameters)) pf$parameters else pf
  beta <- pars$beta %gp_or% pf$beta
  if (is.null(beta) || !is.finite(beta)) {
    .gp_cal_abort("`fit$precision_fit` does not carry a finite `beta`.")
  }
  beta <- as.numeric(beta)

  characteristic <- prior$metadata$characteristic %gp_or%
    pars$characteristic %gp_or%
    pf$characteristic %gp_or%
    fit$provenance$demographic %gp_or% "race"
  characteristic <- match.arg(characteristic, c("race", "gender"))

  # Additive location for gender only; race is purely multiplicative.
  mu <- if (identical(characteristic, "gender")) {
    as.numeric(pars$mu %gp_or% pf$mu %gp_or% pf$provenance$mu %gp_or% 0)
  } else {
    0
  }

  s <- as.numeric(.gp_estimates_se(fit$estimates))
  if (length(s) < 2L || any(!is.finite(s)) || any(s <= 0)) {
    .gp_cal_abort("`fit$estimates` must carry finite positive standard errors `s`.")
  }

  ids <- as.character(.gp_estimates_id(fit$estimates))
  labels <- fit$estimates$covariates$label %gp_or% ids

  list(
    support = support,
    density = density / sum(density),
    beta = beta,
    mu = mu,
    s = s,
    J = length(s),
    characteristic = characteristic,
    prior = prior,
    precision_fit = pf,
    ids = ids,
    labels = as.character(labels),
    control = validate_gp_control(fit$control)
  )
}

# Map standardized v -> reporting-scale theta under the DGP characteristic.
# Race: theta = s^beta * v ; Gender: theta = mu + s^beta * v.
.gp_cal_theta_from_v <- function(v, dgp) {
  base <- (dgp$s ^ dgp$beta) * v
  if (identical(dgp$characteristic, "gender")) dgp$mu + base else base
}

#' Draw one synthetic dataset from the DGP truth
#'
#' Samples latent v_i ~ g_hat, forms theta_i = s_i^beta * v_i (+ mu for gender),
#' and adds reporting-scale sampling noise theta_hat_i = theta_i + s_i * eps_i.
#' Returns the synthetic `data.frame` (the `krw_report_card` input) plus the true
#' latent v and true theta (for the oracle / coverage scoring).
#'
#' @keywords internal
#' @noRd
.gp_cal_simulate_dataset <- function(dgp) {
  J <- dgp$J
  v_true <- sample(dgp$support, size = J, replace = TRUE, prob = dgp$density)
  theta_true <- .gp_cal_theta_from_v(v_true, dgp)
  theta_hat <- theta_true + dgp$s * stats::rnorm(J)
  data <- data.frame(
    theta_hat = theta_hat,
    s = dgp$s,
    unit_id = dgp$ids,
    label = dgp$labels,
    stringsAsFactors = FALSE
  )
  list(data = data, v_true = v_true, theta_true = theta_true)
}

#' Re-estimate the one-level pipeline on a synthetic dataset
#'
#' Runs the full native KRW one-level chain (`krw_report_card`) on the synthetic
#' (theta_hat, s) at `interval_level = ci_level` (so the reporting interval is the
#' requested level) and the truth's solver backend / lambda grid. Returns a
#' `gp_fit` or NULL on failure (the caller treats NULL as a skipped draw).
#'
#' @keywords internal
#' @noRd
.gp_cal_refit <- function(syn, dgp, ci_level) {
  ctrl <- gp_control(
    lambda_grid    = dgp$control$lambda_grid,
    backend        = dgp$control$backend,
    precision_rule = "krw_gmm",
    interval_level = ci_level,
    solver_options = dgp$control$solver_options
  )
  # Calibration deliberately stresses the pipeline on many synthetic draws, some
  # of which do not fully converge; those warnings are expected noise (the
  # non-convergent draws are skipped downstream), so they are suppressed here to
  # keep a clean Monte-Carlo run. Genuine errors still abort the refit -> NULL.
  tryCatch(
    suppressWarnings(
      krw_report_card(syn$data, demographic = dgp$characteristic, control = ctrl)
    ),
    error = function(e) NULL
  )
}

# Reconstruct the TRUE-DGP standardized r-estimates (v_hat, s_v) from the
# synthetic (theta_hat, s). The Monte-Carlo oracle knows the DGP precision law
# that generated the draw, so it must use dgp$beta/dgp$mu here, not the synthetic
# refit's re-estimated beta/mu.
.gp_cal_truth_r_estimates <- function(syn, dgp) {
  beta <- as.numeric(dgp$beta)
  characteristic <- dgp$characteristic
  theta_hat <- syn$data$theta_hat
  s <- syn$data$s
  v_hat <- if (identical(characteristic, "gender")) {
    (theta_hat - dgp$mu) / (s ^ beta)
  } else {
    theta_hat / (s ^ beta)
  }
  list(
    theta_hat = v_hat,
    s = s ^ (1 - beta),
    original_s = s,
    id = as.character(syn$data$unit_id),
    label = as.character(syn$data$label)
  )
}

#' Oracle pairwise from the TRUE prior (no re-deconvolution)
#'
#' Builds the pairwise outranking matrix the grading SHOULD see if the true prior
#' were known: recompute the posterior weights from `dgp$prior` (the truth) against
#' the true-DGP standardized r-estimates, then form pi_ij. Never re-runs the
#' deconvolution -- it injects the known truth density and precision metadata
#' directly.
#'
#' @return A `gp_pairwise` object, or NULL if the oracle wiring fails.
#' @keywords internal
#' @noRd
.gp_cal_oracle_pairwise <- function(refit, syn, dgp) {
  tryCatch({
    r_estimates <- .gp_cal_truth_r_estimates(syn, dgp)
    truth_precision <- list(
      beta = dgp$beta,
      mu = dgp$mu,
      characteristic = dgp$characteristic,
      model_form = if (identical(dgp$characteristic, "gender")) "additive" else "multiplicative"
    )
    weights <- gp_posterior_weights(
      prior = dgp$prior,
      r_estimates = r_estimates,
      precision_fit = truth_precision
    )
    gp_pairwise(
      weights,
      ids = as.character(syn$data$unit_id),
      control = dgp$control
    )
  }, error = function(e) NULL)
}

# Fraction of units whose refit THETA-scale reporting interval covers the TRUE
# theta. The theta scale is beta-invariant, so this stays coherent even when the
# refit's estimated beta differs from the truth's. Returns NA when the reporting
# interval is unavailable or misaligned.
.gp_cal_coverage <- function(refit, theta_true) {
  reporting <- refit$posterior$metadata$reporting
  if (!is.list(reporting) || is.null(reporting$lower) || is.null(reporting$upper)) {
    return(NA_real_)
  }
  lo <- as.numeric(reporting$lower)
  hi <- as.numeric(reporting$upper)
  target <- as.numeric(theta_true)
  if (length(lo) != length(target) || length(hi) != length(target) ||
      any(!is.finite(lo)) || any(!is.finite(hi)) || any(!is.finite(target))) {
    return(NA_real_)
  }
  mean(target >= lo & target <= hi)
}

# Per-draw metrics on a SUCCESSFUL, non-degenerate refit. Returns NULL when the
# refit is degenerate (extreme beta) or the oracle wiring fails -- the caller then
# skips the draw. Otherwise a length-3 named list: dr_emp, covered, regret.
.gp_cal_draw_metrics <- function(refit, syn, dgp, lambda) {
  beta_refit <- as.numeric(
    refit$precision_fit$parameters$beta %gp_or% refit$precision_fit$beta
  )
  if (!is.finite(beta_refit) || abs(beta_refit) > .gp_cal_beta_sane) {
    return(NULL)                       # degenerate r-scale: treat as non-convergence
  }

  emp_grades <- refit$selected_grade$assignment$grade
  pi_emp <- refit$pairwise$matrix

  # (1) ORACLE: regret of the empirical vs oracle grading, both scored on Pi_truth.
  pi_truth_obj <- .gp_cal_oracle_pairwise(refit, syn, dgp)
  if (is.null(pi_truth_obj)) return(NULL)
  pi_truth <- pi_truth_obj$matrix
  oracle_fit <- tryCatch(
    gp_grade(pi_truth_obj, lambda = lambda, control = dgp$control),
    error = function(e) NULL
  )
  if (is.null(oracle_fit)) return(NULL)
  oracle_grades <- oracle_fit$assignment$grade
  regret <- .gp_cal_risk(emp_grades, pi_truth, lambda) -
    .gp_cal_risk(oracle_grades, pi_truth, lambda)

  # Empirical DR: the refit's own grading scored on its own empirical pairwise.
  dr_emp <- .gradepath_frontier_metrics(emp_grades, pi_emp)$discordance_rate

  # (2) COVERAGE on the beta-invariant theta scale.
  covered <- .gp_cal_coverage(refit, syn$theta_true)

  if (!is.finite(dr_emp) || !is.finite(regret) || !is.finite(covered)) {
    return(NULL)
  }
  list(dr_emp = dr_emp, covered = covered, regret = regret)
}


# ---- the gp_calibration object ----------------------------------------------

.gp_calibration_fields <- c(
  "n_sim", "seed", "characteristic",
  "n_ok", "n_failed",
  "dr_mean", "dr_target",
  "coverage", "ci_level",
  "regret_mean",
  "dr_ok", "coverage_ok",
  "provenance"
)

#' Construct a `gp_calibration` result object
#'
#' A compact carrier: scalar summaries plus two boolean verdicts and light
#' provenance. No matrices, no per-draw archive.
#'
#' @keywords internal
#' @noRd
new_gp_calibration <- function(n_sim, seed, characteristic,
                               n_ok, n_failed,
                               dr_mean, dr_target,
                               coverage, ci_level,
                               regret_mean,
                               dr_ok, coverage_ok,
                               provenance = list()) {
  structure(
    list(
      n_sim = as.integer(n_sim),
      seed = as.integer(seed),
      characteristic = as.character(characteristic),
      n_ok = as.integer(n_ok),
      n_failed = as.integer(n_failed),
      dr_mean = as.numeric(dr_mean),
      dr_target = as.numeric(dr_target),
      coverage = as.numeric(coverage),
      ci_level = as.numeric(ci_level),
      regret_mean = as.numeric(regret_mean),
      dr_ok = isTRUE(dr_ok),
      coverage_ok = isTRUE(coverage_ok),
      provenance = provenance
    ),
    class = c("gp_calibration", "list")
  )
}

#' Validate a `gp_calibration` object
#' @keywords internal
#' @noRd
validate_gp_calibration <- function(x) {
  if (!inherits(x, "gp_calibration")) {
    .gp_cal_abort("Expected a `gp_calibration` object.")
  }
  if (any(!.gp_calibration_fields %in% names(x))) {
    missing <- setdiff(.gp_calibration_fields, names(x))
    .gp_cal_abort(
      "`gp_calibration` is missing required field(s): %s.",
      paste(missing, collapse = ", ")
    )
  }
  for (nm in c("n_sim", "seed", "n_ok", "n_failed")) {
    v <- x[[nm]]
    if (!is.integer(v) || length(v) != 1L || is.na(v)) {
      .gp_cal_abort("`gp_calibration$%s` must be a length-1 integer.", nm)
    }
  }
  for (nm in c("dr_mean", "dr_target", "coverage", "ci_level", "regret_mean")) {
    v <- x[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      .gp_cal_abort("`gp_calibration$%s` must be a length-1 finite numeric.", nm)
    }
  }
  for (nm in c("dr_ok", "coverage_ok")) {
    v <- x[[nm]]
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      .gp_cal_abort("`gp_calibration$%s` must be a length-1 logical.", nm)
    }
  }
  if (x$n_ok + x$n_failed != x$n_sim) {
    .gp_cal_abort("`gp_calibration`: n_ok + n_failed must equal n_sim.")
  }
  if (x$dr_mean < 0 || x$dr_mean > 1 || x$coverage < 0 || x$coverage > 1) {
    .gp_cal_abort("`gp_calibration`: dr_mean and coverage must lie in [0, 1].")
  }
  x
}

#' Print a `gp_calibration` -- a compact, result-first one-liner
#'
#' ASCII-only, returns invisibly. House-style result vocabulary only (GP-DEC-14-A).
#'
#' @param x A `gp_calibration`.
#' @param ... Ignored.
#' @method print gp_calibration
#' @export
print.gp_calibration <- function(x, ...) {
  x <- validate_gp_calibration(x)
  flag <- function(ok) if (isTRUE(ok)) "[OK]" else "[--]"
  pct <- round(100 * x$ci_level)
  skipped <- if (x$n_failed > 0L) sprintf(" (%d skipped)", x$n_failed) else ""
  cat(sprintf(
    paste0(
      "<gp_calibration> n_sim=%d seed=%d %s%s | ",
      "DR mean=%.3f (target %.3f) %s | ",
      "coverage=%.3f (%d%%) %s | regret mean=%.4f\n"
    ),
    x$n_sim, x$seed, x$characteristic, skipped,
    x$dr_mean, x$dr_target, flag(x$dr_ok),
    x$coverage, pct, flag(x$coverage_ok),
    x$regret_mean
  ))
  invisible(x)
}

#' Summarize a `gp_calibration` (typed one-row data frame)
#'
#' @param object A `gp_calibration`.
#' @param ... Ignored.
#' @method summary gp_calibration
#' @export
summary.gp_calibration <- function(object, ...) {
  object <- validate_gp_calibration(object)
  data.frame(
    n_sim = object$n_sim,
    n_ok = object$n_ok,
    n_failed = object$n_failed,
    characteristic = object$characteristic,
    dr_mean = object$dr_mean,
    dr_target = object$dr_target,
    dr_ok = object$dr_ok,
    coverage = object$coverage,
    ci_level = object$ci_level,
    coverage_ok = object$coverage_ok,
    regret_mean = object$regret_mean,
    stringsAsFactors = FALSE
  )
}


# ---- the harness ------------------------------------------------------------

#' Seeded Monte-Carlo calibration of the one-level KRW grading method
#'
#' @description
#' `gp_calibrate()` is a standalone, seeded, qualitative Monte-Carlo calibration
#' harness. It treats a fitted prior as ground
#' truth, simulates `n_sim` synthetic datasets from it, re-runs the full one-level
#' pipeline ([krw_report_card()]) on each, and checks that the grading method's own
#' frequentist guarantees hold qualitatively:
#'
#' \itemize{
#'   \item the empirical discordance rate (DR) tracks `dr_target`;
#'   \item posterior credible intervals cover the truth near `ci_level`;
#'   \item the empirical grading's regret against an oracle that KNOWS the true
#'     prior is small and non-negative.
#' }
#'
#' It is NOT a replication of KRW Table E1 (there is no parseable E1 artifact and
#' the Matlab Monte-Carlo stream cannot be reproduced; see the calibration scope
#' in the `m5-two-level-and-calibration` vignette), is NOT seeded to match Matlab (the `15238`
#' default is a cosmetic homage; the seed only fixes in-package reproducibility),
#' and is NOT a headline result or part of any default scorecard. It never runs on
#' package load or on vignette render. The result is a compact `gp_calibration`
#' object: scalar summaries plus the two boolean verdicts `dr_ok` / `coverage_ok`.
#'
#' @details
#' The harness is robust to refit non-convergence. Each draw's refit is wrapped in
#' `tryCatch`; a failed refit (e.g. a singular step-1 GMM moment covariance) or a
#' degenerate refit (an implausibly large `|beta|` that corrupts the r-scale) is
#' skipped and counted in `n_failed`. The reported means are over the successful
#' draws only. If fewer than `min_ok` draws succeed, `gp_calibrate()` errors rather
#' than averaging noise.
#'
#' Two scale subtleties are handled internally: the ORACLE grading is built from
#' the known true prior WITHOUT re-deconvolving (its pairwise comes from
#' `gp_posterior_weights(prior = truth, ...)`), and both gradings are scored on
#' that same truth pairwise, so regret is non-negative by construction. COVERAGE is
#' assessed on the beta-invariant theta scale (the refit's theta-scale reporting
#' interval against the true `theta_i = s_i^beta * v_i`), which stays coherent even
#' when a refit's estimated beta differs from the truth's.
#'
#' Reproducibility: the whole draw loop runs under a single state-preserving seeded
#' context, so the same `seed` yields identical results and the caller's global RNG
#' stream is left untouched.
#'
#' Cost: each draw re-solves a grade integer program once, so a run costs about
#' `n_sim` refit-and-grade solves. The function prints a one-line budget heads-up
#' (an ordinary, suppressible `message()`) before the loop whenever
#' `n_sim * <number of units>` exceeds 600; the message is informational and never
#' changes the result. For a quick check use a small `n_sim`.
#'
#' On tiny or weakly identified fits the verdicts are frequently `FALSE` -- the
#' small-N KRW beta-GMM is noisy and its r-scale can be ill-identified -- so
#' meaningful calibration needs a well-identified fit (KRW's real multi-firm input).
#' The harness reports honestly either way; it never tunes the data to force a pass.
#'
#' @param fit A `gp_fit` from [krw_report_card()] / [gradepath()] whose fitted
#'   prior, precision exponent `beta` (and additive `mu` for gender), and observed
#'   standard errors `s_i` define the data-generating truth.
#' @param n_sim Integer; number of Monte-Carlo synthetic datasets to draw and
#'   re-fit. Default `200L`. Use a small value (e.g. `5L`) for a fast check; each
#'   draw triggers one grade-IP solve.
#' @param seed Integer; RNG seed for in-package reproducibility. Default `15238L`
#'   (a cosmetic homage to the reference seed -- it does NOT reproduce the Matlab
#'   stream). The same `seed` yields identical output; the caller's RNG stream is
#'   restored afterward.
#' @param dr_target Numeric in `[0, 1]`; the target discordance rate the mean
#'   empirical DR is compared against for the `dr_ok` verdict (passes when within
#'   `0.01`). Default `0.05`.
#' @param ci_level Numeric in the open interval `(0, 1)`; the nominal
#'   credible-interval level, also used as each refit's reporting-interval level.
#'   The `coverage_ok` verdict passes when mean coverage is within `0.02` of this.
#'   Default `0.90`.
#' @param min_ok Integer `>= 1`; the minimum number of successful (non-skipped)
#'   draws required before the harness reports. If fewer draws succeed it errors
#'   rather than averaging noise. Default `1L`.
#'
#' @return A validated `gp_calibration` object (a list of class
#'   `c("gp_calibration", "list")`) summarizing the run: \describe{
#'   \item{`n_sim`, `seed`}{Integer; the requested number of draws and the RNG seed
#'     used.}
#'   \item{`characteristic`}{`"race"` or `"gender"`; the truth's demographic.}
#'   \item{`n_ok`, `n_failed`}{Integer; successful and skipped draw counts
#'     (`n_ok + n_failed == n_sim`).}
#'   \item{`dr_mean`, `dr_target`}{Numeric in `[0, 1]`; the mean empirical
#'     discordance rate over successful draws and the target it is compared with.}
#'   \item{`coverage`, `ci_level`}{Numeric; mean theta-scale interval coverage over
#'     successful draws and the nominal level it is compared with.}
#'   \item{`regret_mean`}{Numeric; mean regret of the empirical grading versus the
#'     true-prior oracle (non-negative by construction).}
#'   \item{`dr_ok`, `coverage_ok`}{Logical; the two qualitative verdicts (DR within
#'     `0.01` of target; coverage within `0.02` of `ci_level`).}
#'   \item{`provenance`}{Named list: producer, selected `lambda`, `beta_truth`,
#'     `n_units`, oracle and coverage-scale tags, and `matlab_rng_parity = FALSE`.}
#' }
#'
#' @note Not a headline result: the published Table E1 Monte-Carlo is deferred --
#'   gp_calibrate is a qualitative in-package diagnostic, not an E1 reproduction.
#'   See the `m5-two-level-and-calibration` vignette.
#'
#' @examples
#' # A pre-solved tiny one-level fit bundled with the package.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#'
#' # A live run re-solves the grade IP once per draw, so it needs a backend and is
#' # slow; keep n_sim small. (The harness also prints a budget heads-up before the
#' # loop when n_sim * n_units > 600.)
#' \donttest{
#' cal <- gp_calibrate(fit, n_sim = 5, seed = 15238)
#' cal                       # compact, result-first one-liner
#' summary(cal)              # one-row typed data frame
#' }
#'
#' @seealso [krw_report_card()], [gp_grade()], [gp_report_card()]
#' @family gradepath-calibrate
#' @export
gp_calibrate <- function(fit, n_sim = 200L, seed = 15238L,
                         dr_target = 0.05, ci_level = 0.90, min_ok = 1L) {
  dgp <- .gp_cal_dgp_from_fit(fit)
  n_sim <- .gp_cal_validate_count(n_sim, "n_sim", min = 1L)
  seed <- .gp_cal_validate_count(seed, "seed", min = 0L)
  dr_target <- .gp_cal_validate_scalar(dr_target, "dr_target", 0, 1)
  ci_level <- .gp_cal_validate_scalar(ci_level, "ci_level", 0, 1, inclusive = FALSE)
  min_ok <- .gp_cal_validate_count(min_ok, "min_ok", min = 1L)

  lambda <- as.numeric(fit$selected_grade$lambda %gp_or% 0.25)

  # Budget heads-up before a potentially long run (one grade-IP solve per draw).
  # Informational only -- suppressible with suppressMessages(), never changes the
  # result.
  budget_msg <- .gp_cal_budget_message(n_sim, length(fit$ids))
  if (!is.null(budget_msg)) {
    message(budget_msg)
  }

  # Seed once; preserve the caller's RNG stream via the package's TLS helper.
  results <- .gp_tls_with_seed(seed, {
    out <- vector("list", n_sim)
    for (i in seq_len(n_sim)) {
      syn <- .gp_cal_simulate_dataset(dgp)
      refit <- .gp_cal_refit(syn, dgp, ci_level)
      out[[i]] <- if (is.null(refit)) {
        NULL
      } else {
        .gp_cal_draw_metrics(refit, syn, dgp, lambda)
      }
    }
    out
  })

  ok <- !vapply(results, is.null, logical(1))
  n_ok <- sum(ok)
  n_failed <- n_sim - n_ok
  if (n_ok < min_ok) {
    .gp_cal_abort(
      paste0(
        "Only %d of %d calibration refits succeeded (need at least %d). ",
        "The synthetic one-level beta-GMM did not converge often enough; ",
        "try a larger / better-identified fit or more draws."
      ),
      n_ok, n_sim, min_ok
    )
  }

  good <- results[ok]
  dr_mean <- mean(vapply(good, function(r) r$dr_emp, numeric(1)))
  coverage <- mean(vapply(good, function(r) r$covered, numeric(1)))
  regret_mean <- mean(vapply(good, function(r) r$regret, numeric(1)))

  validate_gp_calibration(new_gp_calibration(
    n_sim = n_sim,
    seed = seed,
    characteristic = dgp$characteristic,
    n_ok = n_ok,
    n_failed = n_failed,
    dr_mean = dr_mean,
    dr_target = dr_target,
    coverage = coverage,
    ci_level = ci_level,
    regret_mean = regret_mean,
    dr_ok = abs(dr_mean - dr_target) <= 0.01,
    coverage_ok = abs(coverage - ci_level) <= 0.02,
    provenance = list(
      producer = "gp_calibrate",
      lambda = lambda,
      beta_truth = dgp$beta,
      n_units = dgp$J,
      min_ok = min_ok,
      oracle = "true_prior_pairwise_no_redeconvolution",
      coverage_scale = "theta_reporting_scale",
      matlab_rng_parity = FALSE
    )
  ))
}
