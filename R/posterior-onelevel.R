# =============================================================================
# posterior-onelevel.R  --  one-level posterior summaries (gp_posterior)
# -----------------------------------------------------------------------------
# Chapter 11 (one-level slice) + Chapter 9 typed gp_posterior. From the W-seam
# output (W [J x M] + the r-scale support, plus the per-unit reporting_support
# [M x J]) this builds the per-unit one-level posterior summaries:
#
#   r-scale posterior mean    pm_r[j]    = sum_m W[j, m] * support[m]
#   r-scale posterior sd      psd_r[j]   = sqrt(sum_m W[j, m] (support[m]-pm_r[j])^2)
#   central interval at level via the per-unit CDF cumsum(W[j, ]) +
#     findInterval(left.open = TRUE) on support (the SAME bit-exact convention
#     the pairwise stage uses).
#
# The VALIDATED gp_posterior carries the r-scale summaries (validate_gp_posterior
# requires scale == "r"). The per-unit REPORTING-scale (theta) summaries -- the
# same posterior mass ridden against each unit's back-transformed axis
# reporting_support[, j] -- ride in `metadata$reporting` (Chapter 8 W4: only the
# axis is back-transformed, never W).
#
# CROSS-CHECK (frozen invariant #8, one-level oracle): pm_r equals
# ebrecipe::eb_shrink(method = "nonparametric")$posterior$.posterior_mean to
# machine epsilon, because W is GP-W-EXACT and the nonparametric
# posterior mean is the SAME sum_m W[j,m] support[m] (ebrecipe's
# .eb_posterior_mean_np). This is a sanity gate, not the production estimand.
#
# GP-DEC-08-B: no ebrecipe ::: in production R/; the oracle is touched only in
# tests/testthat/test-posterior.R (via the .gp_eb_* seam wrapper / public
# ebrecipe::eb_shrink).
# =============================================================================


#' One-level per-unit posterior summaries from the W-seam output
#'
#' @description
#' Computes the one-level posterior mean, sd, and central interval for each unit
#' from the posterior weight matrix `W` and the r-scale support, and returns a
#' validated `gp_posterior`. The per-unit reporting-scale (theta) summaries are
#' attached in `metadata$reporting` when the reporting axis is available.
#'
#' @param weights Either the list returned by `gp_posterior_weights()` (carrying
#'   `$W`, `$support`, optionally `$reporting_support`) OR a bare `J x M`
#'   row-stochastic matrix (in which case `support` is taken from `prior`).
#' @param prior The prior with public field `support` (length M). Used only to
#'   obtain the r-scale support grid when `weights` is a bare matrix; if
#'   `weights` already carries `$support` that is used and `prior` may be `NULL`.
#' @param r_estimates The r-scale estimates (`theta_hat` = v_hat, `s` = s_v,
#'   optionally `original_s` and `id`/`label`); carried onto the object as the
#'   per-unit `estimate`/`se`.
#' @param precision_fit Optional KRW beta-GMM fit; only used to confirm a
#'   reporting axis can be built (the axis itself comes from `weights`).
#' @param control Optional `gp_control`; `control$interval_level` sets the
#'   central interval (default 0.90 when absent).
#' @param interval_level Optional explicit central-interval level (overrides
#'   `control$interval_level`).
#'
#' @return A validated `gp_posterior` (r-scale), with `metadata$reporting`
#'   carrying the per-unit theta-scale summaries when available.
#' @keywords internal
#' @noRd
gp_posterior_onelevel <- function(weights, prior = NULL, r_estimates,
                                  precision_fit = NULL, control = NULL,
                                  interval_level = NULL) {
  ## ---- unpack W + support + (optional) reporting_support --------------------
  reporting_support <- NULL
  if (is.list(weights) && !is.null(weights$W)) {
    W <- weights$W
    support <- if (!is.null(weights$support)) as.numeric(weights$support)
               else as.numeric(prior$support)
    reporting_support <- weights$reporting_support
  } else if (is.matrix(weights)) {
    W <- weights
    if (is.null(prior) || is.null(prior$support)) {
      .gradepath_abort(
        "`prior$support` is required when `weights` is a bare matrix.")
    }
    support <- as.numeric(prior$support)
  } else {
    .gradepath_abort(
      "`weights` must be a gp_posterior_weights() list or a J x M matrix.")
  }
  if (!is.matrix(W) || !is.numeric(W)) {
    .gradepath_abort("`weights$W` must be a numeric matrix.")
  }
  J <- nrow(W); M <- ncol(W)
  if (length(support) != M) {
    .gradepath_abort(
      "`support` (%d) must match ncol(W) (%d).", length(support), M)
  }

  ## ---- r-scale estimates carried onto the object ---------------------------
  if (!is.list(r_estimates)) {
    .gradepath_abort("`r_estimates` must be a list.")
  }
  theta_hat <- .gradepath_validate_numeric_vector(
    r_estimates$theta_hat, "r_estimates$theta_hat")
  s_v <- .gradepath_validate_numeric_vector(r_estimates$s, "r_estimates$s")
  if (length(theta_hat) != J || length(s_v) != J) {
    .gradepath_abort(
      "`r_estimates$theta_hat`/`s` (%d/%d) must have length J = %d.",
      length(theta_hat), length(s_v), J)
  }
  id <- if (!is.null(r_estimates$id)) as.character(r_estimates$id)
        else if (!is.null(r_estimates$unit_id)) as.character(r_estimates$unit_id)
        else as.character(seq_len(J))
  label <- if (!is.null(r_estimates$label)) as.character(r_estimates$label) else id

  ## ---- interval level ------------------------------------------------------
  level <- if (!is.null(interval_level)) interval_level
           else if (!is.null(control) && !is.null(control$interval_level)) control$interval_level
           else 0.90
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    .gradepath_abort("`interval_level` must be a scalar in (0, 1); got %s.",
                     deparse(level))
  }

  ## ---- r-scale per-unit summaries ------------------------------------------
  ## mean = W %*% support ; sd from the per-unit second moment.
  pm_r <- as.numeric(W %*% support)
  e2 <- as.numeric(W %*% (support^2))
  var_r <- pmax(e2 - pm_r^2, 0)
  psd_r <- sqrt(var_r)

  ## central interval via the per-unit CDF (bit-exact cumsum + findInterval,
  ## left.open = TRUE -- the same convention the pairwise stage uses).
  ci <- .gp_posterior_interval(W, support, level)

  ## ---- per-unit reporting-scale summaries (theta), W4-safe -----------------
  reporting <- NULL
  if (!is.null(reporting_support)) {
    rs <- reporting_support                       # M x J
    if (!is.matrix(rs) || nrow(rs) != M || ncol(rs) != J) {
      .gradepath_abort(
        "`reporting_support` must be M x J = %d x %d.", M, J)
    }
    pm_t <- numeric(J); psd_t <- numeric(J)
    lo_t <- numeric(J); up_t <- numeric(J)
    for (j in seq_len(J)) {
      axis <- rs[, j]
      wj <- W[j, ]
      mj <- sum(wj * axis)
      pm_t[j] <- mj
      psd_t[j] <- sqrt(max(sum(wj * (axis - mj)^2), 0))
      # BOUNDARY CONVENTION (CCR-20): .gp_cdf_index uses findInterval(left.open =
      # TRUE), i.e. the first m with cdf[m] >= p (weak inequality). KRW
      # get_posteriors.m:44-47 uses the strict cdf[m] > p. The two differ ONLY at
      # exact-equality boundaries cdf[m] == p (measure-zero on real data) and
      # affect ONLY this reporting-scale CI, never the posterior mean. No code
      # change for M1 (documented difference, not a parity break).
      cdf <- cumsum(wj)
      lo_t[j] <- axis[.gp_cdf_index(cdf, (1 - level) / 2)]
      up_t[j] <- axis[.gp_cdf_index(cdf, (1 + level) / 2)]
    }
    reporting <- list(
      posterior_mean = pm_t, posterior_sd = psd_t,
      lower = lo_t, upper = up_t, scale = "theta", level = level)
  }

  ## ---- assemble + validate (r-scale object) --------------------------------
  obj <- new_gp_posterior(
    estimate = theta_hat, se = s_v, id = id, label = label,
    posterior_mean = pm_r, posterior_sd = psd_r,
    lower = ci$lower, upper = ci$upper, scale = "r",
    metadata = list(level = level, interval_level = level,
                    reporting = reporting,
                    has_reporting = !is.null(reporting)))
  validate_gp_posterior(obj)
}


# ---------------------------------------------------------------------------
# .gp_cdf_index  --  bit-exact per-unit CDF quantile index
# ---------------------------------------------------------------------------
# Given a per-unit CDF (cumsum of row-stochastic weights) and a target
# probability p, return the support index of the p-quantile via
# findInterval(left.open = TRUE) + 1 (the smallest m with CDF[m] >= p), clamped
# to [1, M]. This is the SAME cumsum + findInterval(left.open = TRUE) convention
# the pairwise integral uses, so the posterior interval and the pairwise
# probabilities agree bit-for-bit on the grid.
#' @keywords internal
#' @noRd
.gp_cdf_index <- function(cdf, p) {
  M <- length(cdf)
  idx <- findInterval(p, cdf, left.open = TRUE) + 1L
  if (idx < 1L) idx <- 1L
  if (idx > M) idx <- M
  idx
}


# ---------------------------------------------------------------------------
# .gp_posterior_interval  --  r-scale central interval for every unit
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.gp_posterior_interval <- function(W, support, level) {
  J <- nrow(W)
  lo <- numeric(J); up <- numeric(J)
  p_lo <- (1 - level) / 2
  p_up <- (1 + level) / 2
  for (j in seq_len(J)) {
    cdf <- cumsum(W[j, ])
    lo[j] <- support[.gp_cdf_index(cdf, p_lo)]
    up[j] <- support[.gp_cdf_index(cdf, p_up)]
  }
  list(lower = lo, upper = up)
}
