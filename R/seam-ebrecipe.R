# =============================================================================
# seam-ebrecipe.R  --  THE ebrecipe boundary
# -----------------------------------------------------------------------------
# This file is the SINGLE controlled point of contact between gradepath and
# ebrecipe's INTERNALS (GP-DEC-07-A). gradepath reaches
# ebrecipe internals through exactly one file of `.gp_eb_*` wrappers and nowhere
# else; every other gradepath function calls these wrappers, never an `ebrecipe`
# internal directly.
#
# Each reused internal is resolved LAZILY (at call time) via
#   getFromNamespace("<name>", "ebrecipe")
# and NEVER via the triple-colon operator. Why getFromNamespace and not the
# triple-colon form:
#   * triple-colon access to an unexported object triggers an `R CMD check` NOTE
#     ("Unexported objects imported by ':::'"); getFromNamespace is check-clean.
#   * ebrecipe is a hard `Imports` whose internal numerical building blocks
#     gradepath's stage-2 core reuses. Pinning the contact to this one file makes
#     any upstream drift a single-file, checksum-guarded failure instead of
#     silent numerical drift scattered
#     across the package.
#   * Lazy resolution avoids load-order coupling between the two namespaces.
#
# Each wrapper is a THIN, FAITHFUL pass-through: identical arguments, identical
# return value, NO transformation. The wrapper exists only to name and localise
# the dependency -- which is what lets the parity test (test-seam-primitives.R)
# assert bit-for-bit identity against the raw internal.
#
# Naming: reused ebrecipe internals are wrapped as `.gp_eb_<primitive>`
# (drop the leading `.eb_`, prefix `.gp_eb_`). These are ebrecipe's functions,
# not gradepath's. gradepath's OWN native helpers are `.gp_*`/`gp_*` (no `eb_`)
# and live elsewhere -- e.g. the NATIVE two-level `gp_pushforward_theta`
# is distinct from the reused one-level `.gp_eb_pushforward_theta`
# cross-check below (frozen invariant 8).
#
# NOTE on signatures: every signature below is the REAL installed-ebrecipe
# (0.5.0) signature, verified against `args(getFromNamespace(...))` and the
# function bodies. Several differ from earlier drafts; the wrappers
# track the real internals.
# =============================================================================


# ---- Spline basis ----------------------------------------------------------

#' Reused ebrecipe spline-basis primitive
#'
# Wraps ebrecipe internal `.eb_spline_basis(support, n_knots = 5L)`: the
# natural-spline (ns) basis matrix Q on `support`. Reused so gradepath builds the
# SAME basis ebrecipe fits in stage 1 (crosswalk: ns, df = 5, center +
# unit-norm).
# NOTE: the basis primitive is `.eb_spline_basis` ONLY -- there is NO
# `.eb_spline_basis_matrix` (an earlier draft named a non-existent
# `_matrix` suffix; do not wrap it).
#' @keywords internal
#' @noRd
.gp_eb_spline_basis <- function(support, n_knots = 5L) {
  fn <- utils::getFromNamespace(".eb_spline_basis", "ebrecipe")
  fn(support, n_knots = n_knots)
}


# ---- Softmax density -------------------------------------------------------

#' Reused ebrecipe softmax-density primitive
#'
# Wraps ebrecipe internal `.eb_softmax_density(Q, alpha)`: the softmax
# exp-spline density map; returns a list(g, log_g) where g is the normalized
# grid mass (crosswalk: the softmax map of likelihood.m). The fitted g used
# throughout stage 2.
#' @keywords internal
#' @noRd
.gp_eb_softmax_density <- function(Q, alpha) {
  fn <- utils::getFromNamespace(".eb_softmax_density", "ebrecipe")
  fn(Q, alpha)
}


# ---- Normal-mixture (likelihood) matrix ------------------------------------

#' Reused ebrecipe normal-mixture-matrix primitive
#'
# Wraps ebrecipe internal
#   `.eb_normal_mixture_matrix(theta_hat, s, support, log = TRUE, warn_threshold = 1e7)`:
# the J x M per-firm normal (log-)likelihood matrix in the posterior
# (crosswalk: the one-level building block of likelihood.m / get_posteriors.m).
#' @keywords internal
#' @noRd
.gp_eb_normal_mixture_matrix <- function(theta_hat, s, support,
                                         log = TRUE, warn_threshold = 1e7) {
  fn <- utils::getFromNamespace(".eb_normal_mixture_matrix", "ebrecipe")
  fn(theta_hat, s, support, log = log, warn_threshold = warn_threshold)
}


# ---- Penalized log-likelihood objective ------------------------------------

#' Reused ebrecipe penalized-log-likelihood primitive
#'
# Wraps ebrecipe internal `.eb_penalized_loglik(alpha, Q, log_P, penalty = 0)`:
# the penalized log-likelihood objective in the spline coefficients `alpha`
# (log_P is the normal-mixture LOG-likelihood matrix). gradepath's deconvolution
# reuses the SAME objective core so stage-2 estimation sees
# identical arithmetic.
#' @keywords internal
#' @noRd
.gp_eb_penalized_loglik <- function(alpha, Q, log_P, penalty = 0) {
  fn <- utils::getFromNamespace(".eb_penalized_loglik", "ebrecipe")
  fn(alpha, Q, log_P, penalty = penalty)
}


# ---- Mean-constraint root-solve helpers ------------------------------------

#' Reused ebrecipe mean-constraint solve (last coefficient)
#'
# Wraps ebrecipe internal
#   `.eb_solve_alpha_T(alpha_free, Q, support, target_mean, interval = c(-10, 10),
#                      max_expansions = 3L, tol = 1e-11)`:
# the scalar mean-constraint root-solve (KRW's fsolve on the last spline
# coefficient) pinning E[xi] = target_mean. Reused so gradepath uses the
# IDENTICAL constraint solve (crosswalk).
#' @keywords internal
#' @noRd
.gp_eb_solve_alpha_T <- function(alpha_free, Q, support, target_mean,
                                 interval = c(-10, 10),
                                 max_expansions = 3L, tol = 1e-11) {
  fn <- utils::getFromNamespace(".eb_solve_alpha_T", "ebrecipe")
  fn(alpha_free, Q, support, target_mean,
     interval = interval, max_expansions = max_expansions, tol = tol)
}

#' Reused ebrecipe free->full alpha expansion under the mean constraint
#'
# Wraps ebrecipe internal `.eb_full_alpha(alpha_free, Q, support, target_mean, ...)`:
# expands the free coefficient vector to the full constrained `alpha` (appending
# the mean-pinning last coefficient from `.eb_solve_alpha_T`). Same constraint
# machinery reused by gradepath (crosswalk).
#' @keywords internal
#' @noRd
.gp_eb_full_alpha <- function(alpha_free, Q, support, target_mean, ...) {
  fn <- utils::getFromNamespace(".eb_full_alpha", "ebrecipe")
  fn(alpha_free, Q, support, target_mean, ...)
}


# ---- Row-wise log-sum-exp --------------------------------------------------

#' Reused ebrecipe row-wise log-sum-exp
#'
# Wraps ebrecipe internal `.eb_row_log_sum_exp(x)`: numerically-stable row-wise
# log-sum-exp over the matrix `x`, reused for stable normalization in stage 2
# (crosswalk: the row-logsumexp building block).
#' @keywords internal
#' @noRd
.gp_eb_row_log_sum_exp <- function(x) {
  fn <- utils::getFromNamespace(".eb_row_log_sum_exp", "ebrecipe")
  fn(x)
}


# ---- One-level pushforward (CROSS-CHECK ONLY -- invariant 8) ----------------

#' Reused ebrecipe one-level theta pushforward (cross-check ONLY)
#'
# Wraps ebrecipe internal
#   `.eb_pushforward_theta(support, g, s, psi_1, psi_2, characteristic = NULL,
#                          standardization_model = NULL)`:
# the ONE-LEVEL (single-support) r -> theta pushforward kernel. It REQUIRES one
# of `characteristic` (e.g. "white"=multiplicative / "male"=additive) or
# `standardization_model` to pick the scale, and errors if both are NULL.
# FROZEN INVARIANT 8: this reused one-level version is used in gradepath ONLY as
# a `group_fx==0` parity / cross-check anchor (it coincides with KRW's
# get_g_theta.m only when the psi/eta-support degenerates to a point). gradepath's
# OWN estimand requires a NATIVE TWO-LEVEL pushforward `gp_pushforward_theta`,
# built natively (named `gp_*`, NOT `.gp_eb_*`). Do NOT conflate them and do
# NOT route the two-level path through this one-level kernel
# (crosswalk; frozen invariant 8).
#' @keywords internal
#' @noRd
.gp_eb_pushforward_theta <- function(support, g, s, psi_1, psi_2,
                                     characteristic = NULL,
                                     standardization_model = NULL) {
  fn <- utils::getFromNamespace(".eb_pushforward_theta", "ebrecipe")
  fn(support, g, s, psi_1, psi_2,
     characteristic = characteristic,
     standardization_model = standardization_model)
}


# ---- Posterior weights (the W primitive) -----------------------------------

#' Reused ebrecipe posterior-weights primitive (the W primitive)
#'
# Wraps ebrecipe internal `.eb_posterior_weights(estimates, prior)`: the per-firm
# posterior weight matrix W from an `eb_estimates` and an `eb_prior`. Reused here
# because the GP-W-EXACT check compares gradepath's native W recompute against
# THIS exact primitive (max|W_gp - this| == 0) (crosswalk; invariant 5).
#' @keywords internal
#' @noRd
.gp_eb_posterior_weights <- function(estimates, prior) {
  fn <- utils::getFromNamespace(".eb_posterior_weights", "ebrecipe")
  fn(estimates, prior)
}


# ---- Stage-1 entry: the public eb_input boundary call ----------------------

#' gradepath stage-1 entry over the exported `ebrecipe::eb_input`
#'
# The single stage-1 call gradepath makes. Unlike the primitives above (ebrecipe
# INTERNALS reached via getFromNamespace), `eb_input` is EXPORTED, so this is a
# public boundary call -- but it is kept in this same seam file so the entire
# ebrecipe contact surface lives in one place ("Stage 1 (CALL)").
#
# Takes firm-level estimates theta_hat / s (and optional unit ids, sample sizes,
# covariates, description) and returns the stage-1 `eb_estimates` container that
# everything downstream consumes. The argument names mirror
#   ebrecipe::eb_input(theta_hat, s, unit_id = NULL, n = NULL,
#                      covariates = NULL, description = NULL)
# exactly (verified against the installed signature -- the id arg is `unit_id`,
# not `ids`); `...` is forwarded faithfully.
#' @keywords internal
#' @noRd
.gp_eb_input <- function(theta_hat, s, unit_id = NULL, ...) {
  ebrecipe::eb_input(theta_hat = theta_hat, s = s, unit_id = unit_id, ...)
}

#' Reused ebrecipe estimates validator
#'
#' This is the only production contact point for ebrecipe's internal
#' `validate_eb_estimates()` validator. Composite gradepath validators call this
#' seam wrapper rather than performing namespace introspection outside the seam.
#' @keywords internal
#' @noRd
.gp_eb_validate_estimates <- function(estimates) {
  fn <- utils::getFromNamespace("validate_eb_estimates", "ebrecipe")
  fn(estimates)
}


# ===========================================================================
# Boundary shape-assertion (the middle layer of the GP-DEC-07-A guard)
# ===========================================================================
#
# `gp_assert_ebrecipe_boundary()` is the seam's loud-failure guard. It runs at
# the START of the native core (before any
# deconvolution), so a *schema* drift in the imported ebrecipe -- a
# renamed/removed eb_estimates/eb_prior field, a changed posterior-weight
# orientation, a non-'r' prior scale -- is caught AT THE BOUNDARY with a named
# error, not three stages downstream as a silently-wrong number. It is
# DEFENSIVE, not corrective: it never repairs a drifted object, it stops the run
# so the maintainer re-pins or vendors (the GP-DEC-07-A escape hatch). It
# complements the version pin (which ebrecipe) and the primitive checksum suite
# (the numbers are unchanged); this layer asserts the fields are
# still there and correctly shaped.
#
# Field requirements are the gradepath-CONSUMED subset (a superset check would
# couple gradepath to ebrecipe internals it does not read). Verified live
# against ebrecipe 0.5.0:
#   eb_estimates : theta_hat, s, unit_id, ... (gradepath reads theta_hat, s)
#   eb_prior     : support, density, scale, ... (scale == "r")
#   posterior W  : a bare J x M row-stochastic matrix (rowSums == 1)
#
#' @keywords internal
#' @noRd
gp_assert_ebrecipe_boundary <- function(estimates,
                                        prior = NULL,
                                        weights = NULL,
                                        pairwise = NULL) {
  # (1) eb_estimates schema: required fields present + theta_hat/s length-matched.
  req_est <- c("theta_hat", "s", "unit_id")
  miss <- setdiff(req_est, names(estimates))
  if (length(miss) > 0L) {
    .gradepath_abort(
      "ebrecipe eb_estimates missing field(s): %s.",
      paste(miss, collapse = ", ")
    )
  }
  if (length(estimates$theta_hat) != length(estimates$s)) {
    .gradepath_abort(
      "ebrecipe eb_estimates: theta_hat (%d) and s (%d) must be length-matched.",
      length(estimates$theta_hat), length(estimates$s)
    )
  }

  # (2) eb_prior schema (when supplied): required fields + scale == 'r'
  #     (the standardized-residual scale gradepath's native core operates on).
  if (!is.null(prior)) {
    req_prior <- c("support", "density", "scale")
    miss <- setdiff(req_prior, names(prior))
    if (length(miss) > 0L) {
      .gradepath_abort(
        "ebrecipe eb_prior missing field(s): %s.",
        paste(miss, collapse = ", ")
      )
    }
    if (!identical(prior$scale, "r")) {
      .gradepath_abort(
        "ebrecipe eb_prior$scale must be 'r' (standardized residual scale); got %s.",
        deparse(prior$scale)
      )
    }
    if (length(prior$support) != length(prior$density)) {
      .gradepath_abort(
        "ebrecipe eb_prior: support (%d) and density (%d) must be length-matched.",
        length(prior$support), length(prior$density)
      )
    }
  }

  # (3) posterior weights W (when supplied): a J x M row-stochastic matrix
  #     (the .gp_eb_posterior_weights orientation; frozen invariant 5). A
  #     transposed / column-stochastic W from an upstream change fails here.
  if (!is.null(weights)) {
    if (!is.matrix(weights) || !is.numeric(weights)) {
      .gradepath_abort("ebrecipe posterior weights W must be a numeric matrix.")
    }
    rs <- rowSums(weights)
    if (any(abs(rs - 1) > 1e-6)) {
      .gradepath_abort(
        "ebrecipe posterior weights W must be row-stochastic (row sums == 1); max deviation %.2e.",
        max(abs(rs - 1))
      )
    }
  }

  # (4) pairwise Pi (when supplied): square J x J with diagonal 0.5
  #     (the gradepath diagonal-projection convention; frozen invariant 6).
  if (!is.null(pairwise)) {
    if (!is.matrix(pairwise) || nrow(pairwise) != ncol(pairwise)) {
      .gradepath_abort("pairwise Pi must be a square J x J matrix.")
    }
    if (any(abs(diag(pairwise) - 0.5) > 1e-8)) {
      .gradepath_abort("pairwise Pi diagonal must be 0.5.")
    }
  }

  invisible(TRUE)
}
