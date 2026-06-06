# =============================================================================
# posterior-weights.R  --  The engine<->grading seam: posterior weight matrix W
#                          (Chapter 8)
# -----------------------------------------------------------------------------
# gp_posterior_weights(prior, r_estimates, precision_fit) recomputes the J x M
# row-stochastic posterior weight matrix W -- the single load-bearing object that
# crosses from gradepath's native deconvolution to the social-choice / grading
# layer (it is what gp_pairwise integrates into pi_ij, Chapter 4).
#
# DECISION GP-DEC-08-B (sec-ch08-acquire / sec-ch08-decisions): W is recomputed
# gradepath-NATIVE from the PUBLIC prior fields (prior$support, prior$density,
# prior$scale) plus the r-scale estimates (r_estimates$theta_hat == v_hat,
# r_estimates$s == s_v). gradepath NEVER calls ebrecipe:::.eb_posterior_weights
# (nor .eb_prior_grid_mass / .eb_safe_log / .eb_safe_normalize) in production R/.
# That internal is referenced ONLY in tests/testthat/test-w-seam-golden.R,
# through the seam wrapper .gp_eb_posterior_weights(), as a bit-exact CI tripwire
# (GP-W-EXACT: max|W_recompute - internal| == 0).
#
# RECIPE -- replicates ebrecipe 0.5.0 .eb_posterior_weights EXACTLY, against the
# public fields, so the golden master is 0 to the last bit (verified live on the
# preserved krw_firms race/multiplicative 97 x 1000 W). The deparsed oracle is:
#
#   mass      <- .eb_prior_grid_mass(support, density)        # scalar mean-spacing MASS
#   log_mass  <- .eb_safe_log(mass)                           # -Inf for <=0
#   log_lik   <- .eb_normal_mixture_matrix(theta_hat, s, support, log=TRUE)  # J x M
#   log_post  <- sweep(log_lik, 2L, log_mass, `+`)            # J x M
#   log_denom <- .eb_row_log_sum_exp(log_post)                # length J
#   weights   <- matrix(0, J, M)
#   valid     <- is.finite(log_denom)
#   weights[valid, ]  <- exp(log_post[valid, ] - log_denom[valid])   # @eq-w
#   weights[!valid, ] <- mass (recycled byrow)                       # W5 fallback
#   .eb_safe_normalize(weights, margin = 1L)                  # clamp rows to 1
#
# MASS, NOT DENSITY (sec-ch08-recipe @eq-mass, the uniform-spacing trap) -- the
# oracle forms the prior GRID MASS as the raw density scaled by the MEAN grid
# spacing, then safe-normalized to sum 1. Verified bit-for-bit against the
# deparsed ebrecipe 0.5.0 internal:
#     .eb_prior_grid_mass(support, density) =
#         .eb_safe_normalize(pmax(density, 0) * mean(diff(support)))
# i.e. a SINGLE scalar spacing = mean(diff(support)), NOT trapezoidal half-width
# endpoints. (An earlier trapezoidal guess was off by 4.03e-06 in the resulting
# W -- the half-width endpoints are wrong; the real engine uses the scalar mean
# spacing. On a uniform grid the constant spacing cancels under the normalize, so
# mass == density / sum(density), but we replicate the exact engine form so the
# golden master is 0 to the last bit.) gp_posterior_weights replicates this
# natively (.gp_prior_grid_mass), confirmed max|gp_mass - eb_mass| == 0.
#
# CONTRACT (sec-spec-w, invariants W1-W5, frozen invariant #5):
#   W1  dim(W) == c(J, M) ............ units are ROWS, grid points are COLUMNS
#   W2  all(abs(rowSums(W)-1) < 1e-12) each ROW is one unit's posterior over the grid
#   W3  W is on the r-scale (prior$scale == "r"); columns index prior$support
#   W4  W is NEVER back-transformed to theta; only the SUPPORT is, per unit
#   W5  all(W >= 0) && all(is.finite(W)); a degenerate row falls back to the mass
#
# ORIENTATION HAZARD (sec-ch08-orientation): v1 shipped W transposed (M x J,
# column-normalized). This file returns J x M, row-normalized, and the
# .gp_assert_w_row_stochastic(W, J, M, tol) guard converts any transpose into a
# hard, located error rather than a silently-wrong pi_ij.
# =============================================================================


#' Recompute the nonparametric posterior weight matrix W from a fitted prior
#'
#' Realizes the four-step posterior recipe (@eq-w / @eq-mass) against the PUBLIC
#' prior fields (`support`, `density`, `scale`) plus r-scale estimates. Does NOT
#' call `ebrecipe:::.eb_posterior_weights` or its mass/log/normalize helpers
#' (GP-DEC-08-B); it reuses only the seam primitives
#' `.gp_eb_normal_mixture_matrix()` and `.gp_eb_row_log_sum_exp()` and replicates
#' the scalar mean-spacing grid-mass / safe-log / safe-normalize arithmetic natively, so
#' the result is bit-identical to the ebrecipe internal on the one-level slice
#' (GP-W-EXACT, the test-only golden master).
#'
#' Returns W on the r-scale (contract sec-spec-w) plus the per-unit reporting
#' support for the pairwise stage (sec-ch08-reporting). W4: the reporting
#' back-transform is applied to the SUPPORT axis per unit, never to W.
#'
#' @param prior A validated prior with the public `eb_prior` surface:
#'   `prior$support` (length-M numeric, strictly increasing), `prior$density`
#'   (length-M numeric, >= 0), and `prior$scale` (the string `"r"`). Accepts
#'   ebrecipe's `eb_prior` AND gradepath's native `gp_prior` (both expose these
#'   three public fields); no KRW reporting coefficient is read from
#'   `prior$spline_info`.
#' @param r_estimates Standardized r-scale estimates supplying `theta_hat`
#'   (length-J, the r-scale point estimates `v_hat`), `s` (length-J, the r-scale
#'   noise scales `s_v = s^(1-beta)`), and -- when a reporting support is wanted --
#'   `original_s` (length-J, the pre-standardization reporting-scale standard
#'   errors). Assemble these from the FLAT `gp_estimation_core()` fit and the GMM
#'   input: `theta_hat = fit$v_hat`, `s = fit$s_v`, and
#'   `original_s = gp_krw_gmm_input()$s` (the pre-standardization SE). (The
#'   estimation-core fit is a flat list -- it carries no `$standardized`
#'   sub-list.)
#' @param precision_fit gradepath's KRW beta-GMM fit -- either the flat
#'   `gp_estimation_core()` list or the classed `gp_estimation_fit` returned by
#'   `gp_w_seam()`. The reporting support reads `beta` and (additive path only)
#'   `mu` from the fit's top level, wrapped parameters, or classed-fit provenance
#'   and derives the model form from `precision_fit$characteristic`: race =
#'   `"multiplicative"` (the preserved race path), gender = `"additive"` (adds the
#'   mu intercept). A legacy parameters-wrapped fit
#'   (`precision_fit$parameters$beta/mu/model_form`, e.g. a `gp_precision_fit`) is
#'   still honored when present, but `gp_estimation_core()` itself is flat.
#'   Optional: if NULL (or if `r_estimates$original_s` is absent) the reporting
#'   support is omitted and only the r-scale W + support are returned.
#'
#' @return A list with:
#'   * `W`                 : J x M numeric, row-stochastic, r-scale  (W1-W5)
#'   * `support`           : length-M r-scale grid (= `prior$support`)
#'   * `reporting_support` : M x J numeric, unit-specific theta-scale axis
#'                           (omitted -- not present -- when no reporting metadata)
#' @keywords internal
#' @noRd
gp_posterior_weights <- function(prior, r_estimates, precision_fit = NULL) {
  # ---- 0. validate + extract PUBLIC prior fields (no ::: ; no private state) --
  if (!is.list(prior)) {
    .gradepath_abort("`prior` must be a list with public fields support/density/scale.")
  }
  .gradepath_validate_required_keys(prior, c("support", "density", "scale"), "prior")
  if (!identical(prior$scale, "r")) {                        # W3: r-scale only
    .gradepath_abort(
      "`prior$scale` must be 'r' (standardized residual scale); got %s.",
      deparse(prior$scale)
    )
  }
  support <- as.numeric(prior$support)
  density <- as.numeric(prior$density)
  if (length(support) != length(density)) {
    .gradepath_abort(
      "`prior$support` (%d) and `prior$density` (%d) must be length-matched.",
      length(support), length(density)
    )
  }
  if (any(!is.finite(support)) || any(!is.finite(density))) {
    .gradepath_abort("`prior$support` and `prior$density` must be finite.")
  }

  # ---- r-scale estimates (the LIKELIHOOD inputs only -- never original_s here) -
  if (!is.list(r_estimates)) {
    .gradepath_abort("`r_estimates` must be a list with `theta_hat` and `s`.")
  }
  .gradepath_validate_required_keys(r_estimates, c("theta_hat", "s"), "r_estimates")
  theta_hat <- .gradepath_validate_numeric_vector(r_estimates$theta_hat, "r_estimates$theta_hat")
  s_v       <- .gradepath_validate_numeric_vector(r_estimates$s, "r_estimates$s")
  if (length(theta_hat) != length(s_v)) {
    .gradepath_abort(
      "`r_estimates$theta_hat` (%d) and `r_estimates$s` (%d) must be length-matched.",
      length(theta_hat), length(s_v)
    )
  }

  J <- length(theta_hat)
  M <- length(support)

  # ---- THE RECIPE -- mirrors ebrecipe 0.5.0 .eb_posterior_weights EXACTLY -----
  # (mass = density * mean(diff(support)), normalized; NOT raw density; see header.)

  # 1. prior GRID MASS = density * mean grid spacing, then renormalize to 1.
  mass     <- .gp_prior_grid_mass(support, density)        # @eq-mass, native
  log_mass <- .gp_safe_log(mass)                           # -Inf for mass <= 0

  # 2. log-likelihood J x M : row j uses unit j's own s_v[j] (seam primitive).
  log_lik  <- .gp_eb_normal_mixture_matrix(theta_hat, s_v, support, log = TRUE)

  # 3. log posterior (un-normalized): add column-indexed log prior mass to each row.
  log_post <- sweep(log_lik, 2L, log_mass, `+`)

  # 4. per-row log marginal likelihood (the denominator of @eq-w).
  log_denom <- .gp_eb_row_log_sum_exp(log_post)

  # 5. exponentiate the valid rows; degenerate rows fall back to the prior MASS
  #    (W5) -- exactly the oracle's valid_rows / !valid_rows branches.
  W <- matrix(0, nrow = J, ncol = M)
  valid_rows <- is.finite(log_denom)
  if (any(valid_rows)) {
    centered <- log_post[valid_rows, , drop = FALSE] - log_denom[valid_rows]
    W[valid_rows, ] <- exp(centered)
  }
  if (any(!valid_rows)) {                                   # W5 fallback = prior mass
    W[!valid_rows, ] <- matrix(rep(mass, sum(!valid_rows)),
                               nrow = sum(!valid_rows), byrow = TRUE)
  }

  # 6. final row-normalize (clamps W2 against drift; matches .eb_safe_normalize).
  W <- .gp_safe_normalize_rows(W)

  # ---- CONTRACT GUARD (W1 + W2 + W5) -- frozen invariant #5 ------------------
  # Turns the orientation transpose hazard and any drift into a hard, located
  # error. tol 1e-12 enforces the tighter W2 bound (the helper default is 1e-9).
  .gp_assert_w_row_stochastic(W, n_units = J, n_support = M, tol = 1e-12)

  # ---- per-unit reporting support (M x J): back-transform the AXIS, not W -----
  # W4: W is never mapped to theta; only the SUPPORT is, per unit, from the KRW
  # beta-GMM metadata (original_s, beta, mu, model_form) -- NOT prior$spline_info
  # (sec-ch08-psitrap, hazard 2). Computed only when the metadata is available.
  reporting_support <- .gp_reporting_support(support, r_estimates, precision_fit)

  out <- list(W = W, support = support)
  if (!is.null(reporting_support)) {
    out$reporting_support <- reporting_support
  }
  out
}


# ---------------------------------------------------------------------------
# .gp_prior_grid_mass  --  scalar mean-spacing grid mass (native, == .eb_prior_grid_mass)
# ---------------------------------------------------------------------------
#
# Native bit-for-bit replica of ebrecipe 0.5.0 .eb_prior_grid_mass (deparsed,
# base64-verified). The grid "mass" at each support point is the raw density
# scaled by a SINGLE SCALAR grid spacing = mean(diff(support)), then renormalized
# to sum 1:
#   spacing = mean(diff(support))            # one scalar, NOT per-point widths
#   mass    = pmax(density, 0) * spacing     # density * scalar spacing
#   mass    = mass / sum(mass)               # safe-normalize to sum 1
# This is NOT a trapezoidal half-width form (an earlier trapezoidal half-width
# guess was off by 4.03e-06 in W; the real engine uses the scalar mean spacing).
# On a uniform grid the constant spacing cancels under the normalize, so
# mass == density / sum(density). Recomputed natively (GP-DEC-08-B): NOT a call
# into ebrecipe. Degenerate total (<= 0 or non-finite) falls back to uniform
# 1/n, mirroring the engine.
#' @keywords internal
#' @noRd
.gp_prior_grid_mass <- function(support, density) {
  # NOTE: for a native gp_prior, `density` is ALREADY probability MASSES (sums to
  # 1; the Matlab g_xi/sum(g_xi) convention), NOT a per-point PDF like
  # ebrecipe::eb_prior$density -- see the gp_prior contract in
  # class-objects-output.R. On a uniform grid the scalar spacing below cancels
  # under the renormalize, so this returns the same masses back.
  # Bit-for-bit replica of ebrecipe 0.5.0 .eb_prior_grid_mass (deparsed, verified):
  # the grid mass is the raw density scaled by the MEAN grid spacing, then
  # safe-normalized to sum 1. (NOT the trapezoidal half-width endpoints -- the
  # real engine uses a single scalar spacing = mean(diff(support)); see NOTES /
  # dev/_massdiag.txt. On a uniform grid the constant spacing cancels under the
  # normalize, so mass == density / sum(density).)
  n <- length(support)
  if (n == 0L) {
    return(numeric(0))
  }
  if (n > 1L) {
    spacing <- mean(diff(support))
    if (!is.finite(spacing) || spacing <= 0) {
      .gradepath_abort("`prior$support` must be strictly increasing.")
    }
    return(.gp_safe_normalize_vec(pmax(density, 0) * spacing))
  }
  .gp_safe_normalize_vec(pmax(density, 0))
}

# Vector safe-normalize (native == .eb_safe_normalize on a vector): x / sum(x),
# or uniform 1/n when the total is non-finite or <= 0.
#' @keywords internal
#' @noRd
.gp_safe_normalize_vec <- function(x) {
  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    return(rep(1 / length(x), length(x)))
  }
  x / total
}


# ---------------------------------------------------------------------------
# .gp_safe_log  --  zero-safe log (native, == .eb_safe_log)
# ---------------------------------------------------------------------------
# log(x) for finite x > 0, -Inf elsewhere. Native replica of .eb_safe_log.
#' @keywords internal
#' @noRd
.gp_safe_log <- function(x, eps = .Machine$double.xmin) {
  # Bit-for-bit replica of ebrecipe 0.5.0 .eb_safe_log: log(pmax(x, eps)) with
  # eps = .Machine$double.xmin (NOT -Inf for x <= 0).
  log(pmax(x, eps))
}


# ---------------------------------------------------------------------------
# .gp_safe_normalize_rows  --  row-stochastic normalize (native, == .eb_safe_normalize margin=1)
# ---------------------------------------------------------------------------
# Row-normalize a numeric matrix; a row whose sum is non-finite or <= 0 is set to
# the uniform 1/ncol(x) so the result stays row-stochastic and finite (mirrors
# .eb_safe_normalize(x, margin = 1L) exactly).
#' @keywords internal
#' @noRd
.gp_safe_normalize_rows <- function(x) {
  row_sums <- rowSums(x)
  safe <- is.finite(row_sums) & row_sums > 0
  out <- x
  out[safe, ] <- x[safe, , drop = FALSE] / row_sums[safe]
  if (any(!safe)) {
    out[!safe, ] <- 1 / ncol(x)
  }
  out
}


# ---------------------------------------------------------------------------
# .gp_reporting_support  --  per-unit theta-scale axis (M x J), W4-safe
# ---------------------------------------------------------------------------
#
# For unit j with original reporting-scale SE original_s[j], KRW precision
# exponent beta, intercept mu, and r-scale grid point support[m]:
#   multiplicative (race) : theta_jm = (original_s[j])^beta * support[m]
#   additive (gender)     : theta_jm = mu + (original_s[j])^beta * support[m]
# Returns an M x J matrix (column j = unit j's grid on its own reporting scale)
# so it feeds the column-per-unit pairwise integral of Chapter 4 directly.
#
# Returns NULL when the reporting metadata is unavailable (no precision_fit, or
# no r_estimates$original_s) -- the r-scale W is still fully usable on its own.
#' @keywords internal
#' @noRd
.gp_reporting_support <- function(support, r_estimates, precision_fit) {
  if (is.null(precision_fit) || is.null(r_estimates$original_s)) {
    return(NULL)
  }
  original_s <- .gradepath_validate_numeric_vector(r_estimates$original_s, "r_estimates$original_s")
  if (length(original_s) != length(r_estimates$theta_hat)) {
    .gradepath_abort(
      "`r_estimates$original_s` (%d) must match `r_estimates$theta_hat` (%d).",
      length(original_s), length(r_estimates$theta_hat)
    )
  }
  # Read the KRW beta-GMM metadata from the fit. Accept BOTH the raw
  # gp_estimation_core() fit (flat $beta / $mu / $characteristic) AND a
  # parameters-wrapped fit ($parameters$beta / ...). model_form is derived from
  # the characteristic when not given explicitly: race = multiplicative,
  # gender = additive.
  pars <- if (!is.null(precision_fit$parameters)) precision_fit$parameters else precision_fit
  beta_val <- if (!is.null(pars$beta)) pars$beta else precision_fit$beta
  beta <- .gradepath_validate_scalar_numeric(beta_val, "precision_fit$beta")
  model_form <- pars$model_form
  if (is.null(model_form)) {
    ch <- if (!is.null(pars$characteristic)) pars$characteristic else precision_fit$characteristic
    model_form <- if (identical(ch, "gender")) "additive" else "multiplicative"
  }

  # M x J: outer(support, original_s^beta) gives reporting_support[m, j].
  reporting_support <- outer(support, original_s^beta, `*`)

  if (identical(model_form, "additive")) {
    mu_val <- if (!is.null(pars$mu)) {
      pars$mu
    } else if (!is.null(precision_fit$mu)) {
      precision_fit$mu
    } else {
      precision_fit$provenance$mu
    }
    mu <- .gradepath_validate_scalar_numeric(mu_val, "precision_fit$mu")
    reporting_support <- mu + reporting_support
  } else if (!identical(model_form, "multiplicative")) {
    .gradepath_abort(
      "Unsupported KRW reporting-support model form: %s.", deparse(model_form)
    )
  }
  reporting_support
}
