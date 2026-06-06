#' W-seam: native beta-GMM estimation core
#'
#' @description
#' The parity-fragile heart of `gradepath`: a native two-step efficient GMM
#' estimator for the heteroskedasticity exponent `beta`, with a sandwich
#' covariance, a single-N overidentification J-statistic, and the GMM-implied
#' coupling carriers (`m_hat`, `V_m`) consumed by the deconvolution (Chapter 11).
#'
#' This is a faithful port of the KRW (2024) Matlab reference
#' (`moment_conditions.m`, `gmm_obj.m`, `estimate_lsqnonlin.m`,
#' `get_moments_gmm.m`) for the **one-level** model (`group_fx == 0`). The
#' two-level group-effect branch (`group_fx == 1`, `sigma_psi`) is out of scope
#' here; the `group_fx` argument is retained only to keep the
#' signature stable.
#'
#' All functions are internal (not exported); `gp_estimation_core()` is called by
#' the orchestrator. See the design and the frozen invariants (#3 coupling,
#' #7 single-N).
#'
#' IMPORTANT (naming collision): the GMM weight matrix produced here,
#' `W_gmm = chol(Omega^-1)` (equivalently the quadratic-form weight `Omega^-1`),
#' is a DISTINCT object from the J x M *posterior* weight matrix `W`
#' (the deconvolution posterior weights). Same letter, different
#' object. Nothing in this file produces or touches the posterior W.
#'
#' @section PARITY (race beta = 0.510, gender = 1.255 -- REPRODUCED):
#' On KRW's *actual* GMM input -- the per-firm `(theta_hat, s)` estimates the
#' Matlab `estimate_lsqnonlin.m` consumed, loaded via `gp_krw_gmm_input()` from
#' `inst/extdata/krw-gmm-input/theta_estimates_matlab_<char>.csv` -- this faithful
#' two-step GMM reproduces KRW Table 3 to three decimals: **race beta = 0.5095**
#' (mu = 0.3073, sigma_xi = 0.2066, SE_beta = 0.190 -- vs published
#' 0.510/0.308/0.207/0.190) and **gender beta = 1.2554** (vs 1.255). No clamping,
#' no seed-pinning: the **two-step optimal weighting** (`W = chol(Omega1^-1)`) is
#' what makes 0.510 the robust optimum -- step 2 converges there from KRW's seed,
#' a naive start, and even a high-beta start.
#'
#' WHY THE INPUT MATTERS (a build-time discovery, flagged for the maintainer):
#' an earlier pass fed `ebrecipe::krw_firms` (the bundled example) and
#' got a spurious large-beta optimum (~2.1) -- because `ebrecipe::krw_firms` is a
#' *different* numeric series from KRW's real GMM input (sorted correlation ~0.94,
#' not a rescaling; theta sd ~6x and s sd ~12x smaller). On that wrong-scale input
#' the SEs are << 1, so `s^beta -> 0` and a degenerate large-beta basin opens up;
#' on KRW's real input that basin does not dominate and 0.510 is recovered
#' cleanly. The beta-GMM therefore uses KRW's real input as the parity ground
#' truth (Decision B: the Matlab is the authority). Whether `ebrecipe::krw_firms`
#' should itself be re-based to KRW's GMM input is an OPEN DATA-PROVENANCE QUESTION
#' for the maintainer (it does not block this core, which loads the real input
#' directly). The +-0.001 race gate is exercised in the parity tests.
#'
#' @keywords internal
#' @noRd
NULL

# ---------------------------------------------------------------------------
# eb_estimates -> moment-conditions data matrix
# ---------------------------------------------------------------------------

#' Build the Matlab `[d, theta_hat, s]` data matrix from an `eb_estimates`
#'
#' @description
#' The Matlab moment routines consume a numeric matrix `data = [d, theta_hat, s]`
#' (group, point estimate, standard error). gradepath's stage-1 container is an
#' `ebrecipe::eb_input()` object (an `eb_estimates` list with `theta_hat`, `s`,
#' `unit_id`, ... fields and **no** group column and **no** `characteristic`).
#' This helper assembles the matrix the moment routines expect. The group column
#' `d` is filled with `1`s: at `group_fx == 0` (the only model implemented here)
#' it is never read, and a constant placeholder keeps the column shape stable for
#' the eventual two-level extension that will replace it with the real
#' (integer-coded) industry grouping.
#'
#' @param estimates An `eb_estimates` (from `.gp_eb_input` / `ebrecipe::eb_input`),
#'   or any list carrying numeric `theta_hat` and `s` of equal length.
#' @return A numeric `N x 3` matrix with columns `c("d", "theta_hat", "s")`.
#' @keywords internal
#' @noRd
.gp_estimation_data_matrix <- function(estimates) {
  theta_hat <- as.numeric(estimates$theta_hat)
  s <- as.numeric(estimates$s)
  if (length(theta_hat) == 0L || length(s) == 0L) {
    .gradepath_abort(
      "`input` must carry non-empty `theta_hat` and `s` vectors.",
      class = "gradepath_validation_error"
    )
  }
  if (length(theta_hat) != length(s)) {
    .gradepath_abort(
      sprintf("`theta_hat` (%d) and `s` (%d) must be length-matched.",
              length(theta_hat), length(s)),
      class = "gradepath_validation_error"
    )
  }
  cbind(d = rep(1, length(theta_hat)), theta_hat = theta_hat, s = s)
}

# ---------------------------------------------------------------------------
# Jacobian backend flag (set by gp_jacobian on first use; observable for tests
# and for the build-report SUMMARY).
# ---------------------------------------------------------------------------

.gp_jacobian_backend <- local({
  used <- NA_character_
  list(
    get = function() used,
    set = function(x) used <<- x
  )
})

#' Numerical Jacobian of a vector-valued function
#'
#' @description
#' A self-contained base-R central-difference Jacobian (no external dependency;
#' the previous optional `numDeriv` fast-path was removed so the core has no
#' Suggested-package coupling). The chosen backend is recorded in
#' `.gp_jacobian_backend` (queried by tests / the build report).
#'
#' The Matlab reference obtains the moment Jacobian `G` and the moment-map
#' Jacobian `dm` via `lsqnonlin(..., 'MaxIter', 0)` (a zero-iteration solve that
#' returns the finite-difference Jacobian at the supplied point). The R
#' equivalent is a direct numerical Jacobian; this is a documented numerical
#' (not structural) difference and is parity-safe to within finite-difference
#' tolerance.
#'
#' @param fn A function mapping a length-`p` numeric vector to a length-`m`
#'   numeric vector.
#' @param delta The length-`p` point at which to differentiate.
#' @param eps Central-difference step (fallback path only).
#' @return An `m x p` numeric Jacobian matrix.
#' @keywords internal
#' @noRd
gp_jacobian <- function(fn, delta, eps = 1e-6) {
  .gp_jacobian_backend$set("central-diff")
  p <- length(delta)
  f0 <- fn(delta)
  m <- length(f0)
  J <- matrix(0, nrow = m, ncol = p)
  for (j in seq_len(p)) {
    h <- eps * max(1, abs(delta[j]))
    dp <- delta; dp[j] <- dp[j] + h
    dn <- delta; dn[j] <- dn[j] - h
    J[, j] <- (fn(dp) - fn(dn)) / (2 * h)
  }
  J
}

#' Moment conditions for the beta-GMM model (port of `moment_conditions.m`)
#'
#' @description
#' Forms the four standardized-residual moment conditions for the one-level
#' latent-signal model. RACE is multiplicative (location enters as
#' `mu = exp(mu_raw)`, residual `v_hat = theta_hat / s^beta`); GENDER is additive
#' (location subtracted before dividing, `v_hat = (theta_hat - mu) / s^beta`).
#' This structural asymmetry mirrors the Matlab exactly and must not be
#' "symmetrized".
#'
#' The parameter vector is `delta = c(mu_raw, log_sigma_xi, beta)` (3 params). The
#' 4th raw parameter (`log_sigma_psi`, the two-level group effect) is the
#' `group_fx == 1` branch and is **not** implemented here.
#'
#' @param delta Length-3 numeric `c(mu_raw, log_sigma_xi, beta)`.
#' @param data Numeric N x 3 matrix with columns `[group, theta_hat, se]`
#'   (Matlab `[d, theta_hat, s]`). `group` is unused when `group_fx == 0`.
#' @param characteristic `"race"` or `"gender"`.
#' @param group_fx Integer; must be `0` here (one-level). `1` is out of scope.
#' @param get_cov Logical; if `TRUE`, also return the 4 x 4 moment covariance
#'   `Omega = (1/N) sum_i g_i g_i'`.
#' @return A list with `g` (length-4), `g_i` (N x 4), `Omega` (4 x 4 or `NULL`),
#'   `v_hat`, `s_v` (length-N standardized-residual carriers), and `r` (length-N).
#' @keywords internal
#' @noRd
gp_moment_conditions <- function(delta, data, characteristic,
                                 group_fx = 0L, get_cov = FALSE) {
  if (!identical(as.integer(group_fx), 0L)) {
    .gradepath_abort(
      paste0(
        "`group_fx == 1` (the two-level group-effect / sigma_psi model) is ",
        "out of scope for the W-seam estimation core. Only the ",
        "one-level model (`group_fx == 0`) is implemented here."
      ),
      class = "gradepath_scope_error"
    )
  }

  ## Unpack data (Matlab column order [d, theta_hat, s]).
  theta_hat <- data[, 2L]
  s <- data[, 3L]
  N <- length(s)

  ## Unpack parameters: delta = [mu_raw, log_sigma_xi, beta].
  mu <- delta[1L]
  sigma_xi <- exp(delta[2L])
  beta <- delta[3L]
  sigma_psi <- 0            # group_fx == 0 (one-level)
  sigma <- sqrt(sigma_psi^2 + sigma_xi^2)

  if (identical(characteristic, "race")) {
    ## RACE: multiplicative location, reparametrized positive.
    mu <- exp(mu)
    ## Full Matlab sigma formula (the sigma_psi terms vanish at group_fx == 0,
    ## but we keep them verbatim for line-by-line parity with the reference).
    sigma <- sqrt((sigma_xi^2) * (sigma_psi^2) +
                    (mu^2) * (sigma_psi^2) +
                    (sigma_xi^2))
  }

  ## Form the standardized residual r.
  if (identical(characteristic, "race")) {
    v_hat <- theta_hat / (s^beta)
    s_v <- s^(1 - beta)
    r <- (v_hat - mu) / sqrt((sigma^2) + (s_v^2))
  } else if (identical(characteristic, "gender")) {
    v_hat <- (theta_hat - mu) / (s^beta)
    s_v <- s^(1 - beta)
    r <- v_hat / sqrt((sigma^2) + (s_v^2))
  } else {
    .gradepath_abort(
      sprintf("Unknown `characteristic`: '%s' (expected 'race' or 'gender').",
              characteristic),
      class = "gradepath_validation_error"
    )
  }

  ## Four moment conditions: g_i = [r, r*s, r^2 - 1, (r^2 - 1)*s].
  g_i <- cbind(r, r * s, r^2 - 1, (r^2 - 1) * s)
  colnames(g_i) <- c("r", "r_s", "r2m1", "r2m1_s")
  g <- colMeans(g_i)

  ## Moment covariance Omega = (1/N) sum_i g_i[i,] %o% g_i[i,] = crossprod(g_i)/N.
  ## crossprod(g_i) = t(g_i) %*% g_i = sum_i (g_i[i,] %o% g_i[i,]); identical to
  ## the Matlab loop, just vectorized.
  Omega <- NULL
  if (isTRUE(get_cov)) {
    Omega <- crossprod(g_i) / N
    dimnames(Omega) <- NULL
  }

  list(g = unname(g), g_i = unname(g_i), Omega = Omega,
       v_hat = v_hat, s_v = s_v, r = r)
}

#' Covariance of the moment conditions (port of the `get_cov` branch)
#'
#' @param g_i An N x 4 matrix of per-observation moment contributions.
#' @return The 4 x 4 covariance `(1/N) sum_i g_i g_i'`.
#' @keywords internal
#' @noRd
gp_moment_cov <- function(g_i) {
  N <- nrow(g_i)
  Omega <- crossprod(g_i) / N
  dimnames(Omega) <- NULL
  Omega
}

#' GMM start values (the Matlab reference starts)
#'
#' @description
#' RACE returns the exact Matlab one-step start `[-1.1789, -1.5761, 0.5099]`
#' (`estimate_lsqnonlin.m` L65). GENDER: the Matlab one-level model has **no**
#' published fixed start -- the reference draws 1000 `randn(3,1)` multi-starts
#' and keeps the best. We return a deterministic default `[0, -1.5, 0.5]`
#' (additive location near 0, moderate log-SD, beta ~ 0.5). (The 4-param start
#' `[-1.1598, -1.9752, 0.5193, -0.7595]` in the reference is the **two-level**
#' gender start, out of scope here.)
#'
#' @param characteristic `"race"` or `"gender"`.
#' @return A length-3 numeric start vector `c(mu_raw, log_sigma_xi, beta)`.
#' @keywords internal
#' @noRd
gp_gmm_start <- function(characteristic) {
  if (identical(characteristic, "race")) {
    return(c(-1.1789, -1.5761, 0.5099))
  }
  if (identical(characteristic, "gender")) {
    return(c(0, -1.5, 0.5))
  }
  .gradepath_abort(
    sprintf("Unknown `characteristic`: '%s'.", characteristic),
    class = "gradepath_validation_error"
  )
}

#' Minimize the GMM objective J = N * g' W g for one weighting (port of the
#' `lsqnonlin` minimization in `estimate_lsqnonlin.m` + `gmm_obj.m`)
#'
#' @description
#' The Matlab minimizes via `lsqnonlin` (a local trust-region Levenberg-Marquardt
#' method) on the weighted moment vector `chol(W) %*% g`, whose sum of squares is
#' `g' W g`, **seeded at KRW's own published estimate**. Here we minimize the
#' algebraically-equivalent scalar objective `J = N * g' W g` with `optim`.
#'
#' WEIGHT convention: `W` is the GMM weight in **quadratic-form** convention --
#' step 1 `W = diag(4)`; step 2 `W = solve(Omega1)`. (The Matlab passes
#' `chol(inv(Omega1))` to its residual form; `g' chol' chol g = g' solve(Omega1)
#' g`, the same quadratic form.)
#'
#' OPTIMIZER NOTE (see the file-header PARITY section): on KRW's real GMM input
#' the two-step optimal weighting makes `beta = 0.510` (race) the robust optimum,
#' reached by deterministic `"BFGS"` from the seed -- no clamping or special
#' damping needed. (On the wrong-scale `ebrecipe::krw_firms` the objective is
#' multi-modal with a spurious large-`beta` basin; that is an input-provenance
#' issue, not an optimizer one -- see the file header.)
#'
#' @param data Numeric N x 3 matrix `[group, theta_hat, se]`.
#' @param characteristic `"race"` or `"gender"`.
#' @param W The 4 x 4 GMM weight matrix (quadratic-form convention).
#' @param start Length-3 start `c(mu_raw, log_sigma_xi, beta)`.
#' @param optimizer `"BFGS"` (default), `"Nelder-Mead"`, `"CG"`, `"nlm"`, or
#'   `"L-BFGS-B"` (with `beta` boxed to `[0, beta_max]`).
#' @param reltol Relative convergence tolerance (`optim`) / `gradtol`+`steptol`
#'   for `nlm`.
#' @param maxit Max iterations.
#' @param beta_max Upper box bound on `beta` for the `"L-BFGS-B"` path only
#'   (default `Inf` == effectively unboxed).
#' @return A list with `delta`, `objective` (the minimized `J = N g' W g`),
#'   `convergence` (0 == success), `message`, and `optimizer`.
#' @keywords internal
#' @noRd
gp_gmm_min <- function(data, characteristic, W, start,
                       optimizer = "BFGS", reltol = 1e-10, maxit = 1000L,
                       beta_max = Inf) {
  N <- nrow(data)
  g_of <- function(delta) {
    gp_moment_conditions(delta, data, characteristic,
                         group_fx = 0L, get_cov = FALSE)$g
  }
  obj <- function(delta) as.numeric(N * crossprod(g_of(delta), W %*% g_of(delta)))

  if (identical(optimizer, "nlm")) {
    fit <- stats::nlm(f = obj, p = start, gradtol = reltol,
                      steptol = reltol, iterlim = maxit)
    conv <- if (fit$code %in% c(1L, 2L)) 0L else as.integer(fit$code)
    return(list(delta = fit$estimate, objective = fit$minimum,
                convergence = conv,
                message = sprintf("nlm code %d", fit$code),
                optimizer = "nlm"))
  }

  if (identical(optimizer, "L-BFGS-B")) {
    fit <- stats::optim(
      par = start, fn = obj, method = "L-BFGS-B",
      lower = c(-Inf, -Inf, 0), upper = c(Inf, Inf, beta_max),
      control = list(factr = 1e7, maxit = maxit)
    )
    return(list(delta = fit$par, objective = fit$value,
                convergence = as.integer(fit$convergence),
                message = if (is.null(fit$message)) NA_character_ else fit$message,
                optimizer = "L-BFGS-B"))
  }

  method <- if (optimizer %in% c("BFGS", "Nelder-Mead", "CG")) optimizer else "BFGS"
  fit <- stats::optim(par = start, fn = obj, method = method,
                      control = list(reltol = reltol, maxit = maxit))
  list(delta = fit$par, objective = fit$value,
       convergence = as.integer(fit$convergence),
       message = if (is.null(fit$message)) NA_character_ else fit$message,
       optimizer = method)
}

#' GMM-implied moments for the deconvolution coupling (port of
#' `get_moments_gmm.m`, `extra_moments = 0`, `group_fx = 0`)
#'
#' @description
#' Returns the natural-scale moments the deconvolution will match: RACE ->
#' `[exp(mu_raw), sigma_xi]` (length 2); GENDER -> `[sigma_xi]` (length 1). With
#' `extra_moments = 0` and `group_fx = 0` these depend only on `delta` (the
#' `data` argument is unused; it is retained for signature parity with the Matlab
#' and for the `extra_moments == 1` paper-reporting path). The full
#' `extra_moments == 1` path (`E_theta`, SD of theta, within-share) is included
#' for faithfulness but is NOT used by the estimation core (coupling uses
#' `extra = 0`, frozen invariant #3).
#'
#' @param delta Length-3 numeric `c(mu_raw, log_sigma_xi, beta)`.
#' @param characteristic `"race"` or `"gender"`.
#' @param data Optional N x 3 matrix `[group, theta_hat, se]`; only used when
#'   `extra_moments == 1`.
#' @param group_fx Integer; must be `0` here.
#' @param extra_moments Integer `0` (default; coupling) or `1` (paper reporting).
#' @return A numeric vector of GMM-implied moments.
#' @keywords internal
#' @noRd
gp_get_moments <- function(delta, characteristic, data = NULL,
                           group_fx = 0L, extra_moments = 0L) {
  if (!identical(as.integer(group_fx), 0L)) {
    .gradepath_abort(
      "`group_fx == 1` is out of scope.",
      class = "gradepath_scope_error"
    )
  }

  mu <- delta[1L]
  sigma_xi <- exp(delta[2L])
  beta <- delta[3L]
  sigma_psi <- 0
  sigma <- sqrt(sigma_psi^2 + sigma_xi^2)
  if (identical(characteristic, "race")) {
    mu <- exp(mu)
    sigma <- sqrt((sigma_xi^2) * (sigma_psi^2) +
                    (mu^2) * (sigma_psi^2) +
                    (sigma_xi^2))
  }

  ## Extra (paper-reporting) moments -- not used by the coupling.
  m_extra <- numeric(0)
  if (identical(as.integer(extra_moments), 1L)) {
    if (is.null(data)) {
      .gradepath_abort(
        "`data` is required when `extra_moments == 1`.",
        class = "gradepath_validation_error"
      )
    }
    s <- data[, 3L]
    if (identical(characteristic, "race")) {
      E_theta <- mu * mean(s^beta)
      E_theta_2 <- mean(s^(2 * beta)) * ((sigma^2) + (mu^2))
    } else {
      E_theta <- mu
      E_theta_2 <- (mu^2) + ((sigma^2) * mean(s^(2 * beta)))
    }
    m_extra <- c(E_theta, sqrt(E_theta_2 - (E_theta^2)))
  }

  ## Core moment vector (coupling carriers).
  m <- numeric(0)
  if (identical(characteristic, "race")) {
    m <- mu                       # exp(mu_raw)
  }
  m <- c(m, sigma_xi)
  m <- c(m, m_extra)
  unname(m)
}

#' Two-step efficient GMM driver (port of `estimate_lsqnonlin.m`, GMM section)
#'
#' @description
#' Runs the full two-step estimator on an `eb_estimates` input and returns the
#' fit object plus the deconvolution coupling carriers. Steps:
#' \enumerate{
#'   \item Step 1 (identity weighting): minimize `N g' I g` from the Matlab start.
#'   \item Form `Omega1` at the step-1 optimum; set the GMM weight
#'     `W_gmm = solve(Omega1)` (Matlab `chol(inv(Omega1))`). `chol(W_gmm)` is
#'     computed as a positive-definiteness guard (singular `Omega1` -> clean
#'     abort).
#'   \item Step 2 (efficient weighting): re-minimize `N g' W_gmm g` from delta1.
#'   \item Sandwich covariance `C = solve(t(G) %*% solve(Omega2) %*% G) / N`
#'     with **one** factor of N (frozen invariant #7); `SE = sqrt(diag(C))`.
#'   \item Single-N J-statistic `J_stat = N * g2' W_gmm g2` (frozen invariant #7,
#'     df = 1). The J-stat weight is `W_gmm = solve(Omega1)` -- the SAME weight
#'     used to *find* delta2 -- NOT `solve(Omega2)`; the sandwich uses
#'     `solve(Omega2)`. This Omega1-vs-Omega2 split is verbatim from the Matlab
#'     (`J_hat_2` from the step-2 `lsqnonlin` weighted by `chol(inv(Omega1))`;
#'     `C` uses `inv(Omega2)`).
#'   \item Coupling (frozen invariant #3, computed ONCE): `m_hat`, its Jacobian
#'     `dm`, and `V_m = dm %*% C %*% t(dm)`.
#' }
#' See the file-header PARITY section: on KRW's real GMM input
#' (`gp_krw_gmm_input()`) the returned `beta` reproduces KRW Table 3 (race
#' 0.5095 vs 0.510; gender 1.2554 vs 1.255). Feeding the wrong-scale
#' `ebrecipe::krw_firms` instead yields a spurious large-beta optimum -- an
#' input-provenance issue (open question for the maintainer), not an estimator
#' defect.
#'
#' @param input An `eb_estimates` (from `.gp_eb_input` / `ebrecipe::eb_input`)
#'   carrying the r-scale `theta_hat` and `s`. NOTE: `eb_estimates` carries **no**
#'   `characteristic` field, so `characteristic` is a required, explicit argument
#'   (not inferred from `input`).
#' @param control Optional `gp_control` (see `gp_control`). Only `control$seed`
#'   is read here (to seed the optional gender multi-start). The optimizer /
#'   tolerance / max-iter are NOT part of the `gp_control` surface; they are
#'   taken from the arguments below.
#' @param characteristic `"race"` or `"gender"` (required).
#' @param n_starts Integer; optional extra seeded random starts for GENDER (the
#'   Matlab one-level gender has no fixed start). Set `0L` to skip. Ignored for
#'   race. Default `8L`.
#' @param optimizer Passed to `gp_gmm_min`; default `"BFGS"`.
#' @param reltol Relative convergence tolerance for the inner minimizations.
#' @param maxit Max iterations for the inner minimizations.
#' @return A named list (the beta-GMM fit object) with: `beta`, `mu`,
#'   `sigma_xi`, `delta`, `SE`, `J_stat`, `df`, `m_hat`, `V_m`, `v_hat`, `s_v`,
#'   `characteristic`, plus diagnostics `Omega1`, `Omega2`, `W_gmm`, `W_chol`,
#'   `G`, `C`, `objective`, `convergence`, and a `provenance` stamp.
#' @keywords internal
#' @noRd
gp_estimation_core <- function(input, control = NULL, characteristic = NULL,
                               n_starts = 8L, optimizer = "BFGS",
                               reltol = 1e-10, maxit = 1000L) {
  if (is.null(characteristic)) {
    .gradepath_abort(
      paste0("`characteristic` is required ('race' or 'gender'): the ",
             "eb_estimates input carries no characteristic field."),
      class = "gradepath_validation_error"
    )
  }
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (is.null(control)) control <- gp_control()

  data <- .gp_estimation_data_matrix(input)
  N <- nrow(data)
  n_moments <- 4L
  n_params <- 3L
  df <- n_moments - n_params            # 4 - 3 = 1

  step_min <- function(W, st) {
    gp_gmm_min(data, characteristic, W, st,
               optimizer = optimizer, reltol = reltol, maxit = maxit)
  }

  ## ----- Step 1: identity weighting -----------------------------------------
  W1 <- diag(n_moments)
  start <- gp_gmm_start(characteristic)
  fit1 <- step_min(W1, start)

  ## GENDER has no canonical fixed Matlab start (the reference draws randn x1000).
  ## A small seeded multi-start hardens the gender step-1 fit. Race keeps the
  ## single published start. (Deterministic: seed from control$seed, else 1.)
  if (identical(characteristic, "gender") && n_starts > 0L) {
    if (!is.null(control$seed)) set.seed(control$seed) else set.seed(1L)
    for (k in seq_len(n_starts)) {
      st_k <- start + stats::rnorm(n_params, sd = c(0.5, 0.5, 0.2))
      fk <- tryCatch(step_min(W1, st_k), error = function(e) NULL)
      if (!is.null(fk) && is.finite(fk$objective) &&
          fk$objective < fit1$objective) {
        fit1 <- fk
      }
    }
  }
  delta1 <- fit1$delta

  ## ----- Optimal weight from Omega1 -----------------------------------------
  mc1 <- gp_moment_conditions(delta1, data, characteristic,
                              group_fx = 0L, get_cov = TRUE)
  Omega1 <- mc1$Omega

  ## PD guard: Matlab does chol(inv(Omega1)); a singular/indefinite Omega1 makes
  ## the GMM weight undefined. Surface a clear, data-pointing error.
  W_gmm <- tryCatch(
    solve(Omega1),
    error = function(e) {
      .gradepath_abort(
        paste0("Step-1 moment covariance Omega1 is singular; the optimal GMM ",
               "weight (Matlab `chol(inv(Omega1))`) is undefined. Check for too ",
               "few firms or collinear moments. (", conditionMessage(e), ")"),
        class = "gradepath_singular_error"
      )
    }
  )
  ## Cholesky of the weight = Matlab `chol(inv(Omega1))`; also confirms PD-ness.
  W_chol <- tryCatch(
    chol(W_gmm),
    error = function(e) {
      .gradepath_abort(
        paste0("inv(Omega1) is not positive definite; `chol(inv(Omega1))` ",
               "(the Matlab GMM weight factor) fails. Check the data/moments. (",
               conditionMessage(e), ")"),
        class = "gradepath_singular_error"
      )
    }
  )

  ## ----- Step 2: efficient weighting ----------------------------------------
  fit2 <- step_min(W_gmm, delta1)
  delta2 <- fit2$delta

  ## ----- Moments + covariance at the step-2 optimum -------------------------
  mc2 <- gp_moment_conditions(delta2, data, characteristic,
                              group_fx = 0L, get_cov = TRUE)
  g2 <- mc2$g
  Omega2 <- mc2$Omega

  ## ----- Sandwich covariance (single N; invariant #7) -----------------------
  ## G = d g / d delta at delta2 (4 x 3).
  g_fun <- function(d) {
    gp_moment_conditions(d, data, characteristic,
                         group_fx = 0L, get_cov = FALSE)$g
  }
  G <- gp_jacobian(g_fun, delta2)

  Omega2_inv <- tryCatch(
    solve(Omega2),
    error = function(e) {
      .gradepath_abort(
        paste0("Step-2 moment covariance Omega2 is singular; the sandwich ",
               "covariance is undefined. (", conditionMessage(e), ")"),
        class = "gradepath_singular_error"
      )
    }
  )
  bread <- t(G) %*% Omega2_inv %*% G          # G' Omega2^-1 G  (3 x 3)
  C <- solve(bread) / N                       # one N (invariant #7)
  C <- (C + t(C)) / 2                         # symmetrize tiny FP asymmetry
  SE <- sqrt(diag(C))

  ## ----- Single-N J-statistic (invariant #7) --------------------------------
  ## Uses W_gmm = solve(Omega1) -- the SAME weight that found delta2 (Matlab
  ## J_hat_2 from the step-2 lsqnonlin weighted by chol(inv(Omega1))), NOT
  ## solve(Omega2). df = n_moments - n_params = 1.
  J_stat <- as.numeric(N * crossprod(g2, W_gmm %*% g2))

  ## ----- Coupling carriers (invariant #3: computed ONCE) --------------------
  m_hat <- gp_get_moments(delta2, characteristic, data = data,
                          group_fx = 0L, extra_moments = 0L)
  m_fun <- function(d) {
    gp_get_moments(d, characteristic, data = data,
                   group_fx = 0L, extra_moments = 0L)
  }
  dm <- gp_jacobian(m_fun, delta2)            # length(m_hat) x 3
  if (is.null(dim(dm))) dm <- matrix(dm, nrow = length(m_hat))
  V_m <- dm %*% C %*% t(dm)
  V_m <- (V_m + t(V_m)) / 2                   # symmetrize

  ## ----- Unpack natural-scale parameters ------------------------------------
  mu_raw <- delta2[1L]
  sigma_xi <- exp(delta2[2L])
  beta <- delta2[3L]
  mu <- if (identical(characteristic, "race")) exp(mu_raw) else mu_raw

  ## ----- Provenance ----------------------------------------------------------
  ## Flat named list (every entry top-level) so consumers/tests read
  ## `provenance$step`, `provenance$n_moments`, etc. directly.
  provenance <- .gradepath_new_provenance(
    step = "w-seam:beta-GMM",
    optimizer = optimizer,
    reltol = reltol,
    maxit = maxit,
    n_moments = n_moments,
    n_params = n_params,
    df = df,
    characteristic = characteristic,
    jacobian_backend = .gp_jacobian_backend$get(),
    convergence_step1 = fit1$convergence,
    convergence_step2 = fit2$convergence
  )

  list(
    beta = beta,
    mu = mu,
    sigma_xi = sigma_xi,
    delta = unname(delta2),
    SE = unname(SE),
    J_stat = J_stat,
    df = df,
    m_hat = m_hat,
    V_m = V_m,
    v_hat = mc2$v_hat,
    s_v = mc2$s_v,
    characteristic = characteristic,
    ## diagnostics
    Omega1 = Omega1,
    Omega2 = Omega2,
    W_gmm = W_gmm,
    W_chol = W_chol,
    G = G,
    C = C,
    objective = fit2$objective,
    convergence = c(step1 = fit1$convergence, step2 = fit2$convergence),
    provenance = provenance
  )
}
