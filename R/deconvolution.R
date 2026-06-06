# =============================================================================
# deconvolution.R  --  native ONE-LEVEL deconvolution prior
# -----------------------------------------------------------------------------
# gradepath-NATIVE one-level deconvolution: a faithful base-R port of the
# group_fx == 0 branch of Kline-Rose-Walters' estimate_lsqnonlin.m (the spline /
# softmax / penalty machinery; the actual "decons" math lives in
# matlab_support/estimate_lsqnonlin.m + likelihood.m + get_model_moments.m, NOT
# in the top-level driver decons.m which is just the Table-3 export loop).
#
# This is NOT a call to ebrecipe::eb_deconvolve. It REUSES ebrecipe's validated
# numerical primitives (via the Step-2.1 .gp_eb_* seam wrappers) as internal
# helpers, and uses ebrecipe::eb_deconvolve only as a ONE-LEVEL cross-check
# ORACLE in the tests (frozen invariant #8: one-level cross-check only; the
# native TWO-level path -- group_fx == 1 -- is a LATER step and is deliberately
# NOT implemented here).
#
# It CONSUMES the Step-3.2 coupling carriers off the single `gp_estimation_fit`
# (frozen invariant #3 -- nothing is re-estimated): the standardized-residual
# carriers `v_hat` / `s_v`, the GMM coupling moments `m_hat`, the moment
# covariance `V_m`, the support caps `caps`, and `beta`. The penalty is selected
# by the coupling J = (m_hat - m_model)' V_m^{-1} (m_hat - m_model) (the
# GMM-moment fit of the prior-implied moments), exactly as in the Matlab
# `Js(comb) = (m_hat - m_hat_curr)' inv(V_m) (m_hat - m_hat_curr)` with
# selection `[J_curr, minidx] = min(Js)`.
#
# OUTPUT: a `gp_prior` (the Step-1.5 shell): support, density (probability
# masses summing to 1 on the support grid), mean, scale = "r", diagnostics,
# metadata.  See `new_gp_prior()` / `validate_gp_prior()` in
# R/class-objects-output.R.
#
# MATLAB GROUND-TRUTH LINE CITATIONS (estimate_lsqnonlin.m unless noted):
#   * supp_pts = 1000 for group_fx == 0 ......................... L14-L18
#   * race one-level penalty grid c_xi = 0.001:0.001:0.2 ........ L30-L33
#   * gender one-level penalty grid c_xi = 0.00025:0.00025:0.01 . L40-L43
#   * one-level support caps (race: lo=0; +/-5sigma capped 7sigma) L222-L266
#       (already ported in .gp_coupling_caps(); consumed here as fit$caps)
#   * support grid linspace(lo, hi, supp_pts) .................. L277
#   * one-level psi is degenerate: supp_psi=1 (race) / 0 (gender) L280-L287
#   * per-penalty likelihood max -> g_xi ...................... L366-L369
#   * model moments m_hat_curr = get_model_moments(supp,g_xi,0) . L370
#       (gender keeps only the 2nd moment: m_hat_curr(2)) ...... L371-L373
#   * J = (m_hat - m_hat_curr)' inv(V_m) (m_hat - m_hat_curr) ... L378
#   * penalty selection [J_curr, minidx] = min(Js) ............. L386
#   * final normalisation g_xi = g_xi / sum(g_xi) .............. L469
#   likelihood.m:
#   * softmax g_xi = exp(Q a - logsumexp(Q a)) ................. L43
#   * race one-level alpha = T_xi free coeffs (no mean constr.) . L25 (gender
#       subtracts one and fsolves mean_xi == 0 ................. L26-L29)
#   * penalised objective logL = -(loglik - c*sqrt(a'a)) ....... L95
#   get_model_moments.m: mean = sum(supp.*g)/sum(g);
#       sd = sqrt(sum((supp-mu)^2 g)/sum(g)) .................. L4-L9
#
# ebrecipe-primitive REUSE (via the .gp_eb_* seam; do NOT reimplement the math):
#   .gp_eb_spline_basis(support, n_knots)        -> Q (M x T natural-spline basis)
#   .gp_eb_normal_mixture_matrix(theta, s, supp) -> list(log_P)  (N x M)
#   .gp_eb_softmax_density(Q, alpha)             -> list(g, log_g)  (g sums to 1)
#   .gp_eb_penalized_loglik(alpha, Q, log_P, penalty) -> NEGATIVE penalised
#       log-likelihood (objective to MINIMISE); penalty term is
#       `penalty * ||alpha||_2` (L2 NORM, not squared) -- matches Matlab's
#       c_xi * sqrt(alpha' alpha) (empirically verified: the (obj - obj0) /
#       (c * ||alpha||) ratio is constant across c, while / (c * ||alpha||^2)
#       is not).
#   .gp_eb_full_alpha(alpha_free, Q, supp, target_mean) -> a FULL-length alpha
#       (length T) whose softmax density has mean EXACTLY target_mean (the
#       mean-constraint root-solve; used for GENDER one-level, target_mean = 0).
#
# Internal; gp_deconvolve() may be EXPORTED in a later step (it is the natural
# one-call native one-level deconvolution surface). For now @keywords internal.
# =============================================================================

#' Resolve the package abort helper whether sourced or namespaced
#'
#' Mirrors `.gp_coupling_abort` so this file works both under
#' `devtools::load_all()` (namespace) and when `source()`-d on top of it.
#'
#' @keywords internal
#' @noRd
.gp_decon_abort <- function(msg, class = "gradepath_error") {
  fn <- tryCatch(get(".gradepath_abort", envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (is.function(fn)) return(fn(msg, class = class))
  if (exists(".gradepath_abort", mode = "function")) {
    return(get(".gradepath_abort")(msg, class = class))
  }
  cnd <- structure(
    class = c(class, "gradepath_error", "error", "condition"),
    list(message = msg, call = NULL)
  )
  stop(cnd)
}

#' Resolve a package internal by name (namespace-first), erroring if absent
#' @keywords internal
#' @noRd
.gp_decon_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  .gp_decon_abort(
    sprintf("Required gradepath internal '%s' not found.", name),
    class = "gradepath_internal_error"
  )
}

# ---------------------------------------------------------------------------
# Support grid
# ---------------------------------------------------------------------------

#' One-level support grid from the coupling caps
#'
#' @description
#' Port of `supp_xi = (supp_xi_min : step : supp_xi_max)'` with
#' `supp_pts = 1000` (estimate_lsqnonlin.m L277, L17). The caps themselves are
#' the Step-3.2 `fit$caps` (already a faithful port of the L222-L266 support
#' block, via `.gp_coupling_caps`). We do NOT recompute the caps here
#' (invariant #3).
#'
#' @param caps Named numeric `c(lo = ., hi = .)` (the fit's `caps`).
#' @param supp_pts Integer number of grid points (1000 for one-level).
#' @return Numeric vector of length `supp_pts` from `lo` to `hi` inclusive.
#' @keywords internal
#' @noRd
.gp_decon_support <- function(caps, supp_pts = 1000L) {
  lo <- unname(caps[["lo"]])
  hi <- unname(caps[["hi"]])
  if (!is.finite(lo) || !is.finite(hi) || !(lo < hi)) {
    .gp_decon_abort(
      sprintf("Invalid support caps: lo = %g, hi = %g (need lo < hi).", lo, hi),
      class = "gradepath_validation_error"
    )
  }
  seq(lo, hi, length.out = as.integer(supp_pts))
}

# ---------------------------------------------------------------------------
# Prior-implied moments (get_model_moments.m, one-level)
# ---------------------------------------------------------------------------

#' Mean and standard deviation of a discrete prior on a grid
#'
#' @description
#' Port of `get_model_moments.m` with `higher_moments = 0`:
#' `mu = sum(supp .* g) / sum(g)`, `sigma = sqrt(sum((supp - mu)^2 g) / sum(g))`.
#' `g` is the probability-mass vector on `support` (sums to 1, but we divide by
#' `sum(g)` defensively exactly as the Matlab does).
#'
#' @param support Numeric grid.
#' @param g Numeric probability masses on `support`.
#' @return Numeric `c(mean = ., sd = .)`.
#' @keywords internal
#' @noRd
.gp_decon_model_moments <- function(support, g) {
  sg <- sum(g)
  mu <- sum(support * g) / sg
  sigma <- sqrt(sum(((support - mu)^2) * g) / sg)
  c(mean = mu, sd = sigma)
}

#' The coupling J for one prior fit
#'
#' @description
#' `J = (m_hat - m_model)' solve(V_m) (m_hat - m_model)` (estimate_lsqnonlin.m
#' L378). `m_model` is built from the prior-implied (mean, sd) to match the
#' SHAPE of the coupling `m_hat`:
#'   * race:   `m_hat = c(mu, sigma_xi)` -> `m_model = c(prior_mean, prior_sd)`.
#'   * gender: `m_hat = c(sigma_xi)`     -> `m_model = c(prior_sd)` (the Matlab
#'             keeps only `m_hat_curr(2)`, the std dev; L371-L373).
#'
#' @param m_hat Coupling moments (length 2 race / 1 gender).
#' @param V_m Moment covariance (same dimension).
#' @param mom Named `c(mean, sd)` from `.gp_decon_model_moments()`.
#' @param characteristic "race" or "gender".
#' @return Scalar J (>= 0 for PSD V_m).
#' @keywords internal
#' @noRd
.gp_decon_J <- function(m_hat, V_m, mom, characteristic) {
  m_model <- if (identical(characteristic, "race")) {
    c(unname(mom[["mean"]]), unname(mom[["sd"]]))
  } else {
    unname(mom[["sd"]])
  }
  d <- as.numeric(m_hat) - m_model
  as.numeric(crossprod(d, solve(V_m, d)))
}

# ---------------------------------------------------------------------------
# Single penalised log-spline fit
# ---------------------------------------------------------------------------

#' Fit the softmax-density prior for ONE penalty (one-level)
#'
#' @description
#' Maximises the penalised log-likelihood over the spline coefficients `alpha`
#' for a single penalty `c_xi`, returning the implied density `g` on the
#' support, the (negative) penalised objective value, and the fitted `alpha`.
#'
#' Faithful to likelihood.m (group_fx == 0):
#'   * RACE: all `T = ncol(Q)` coefficients are FREE (no mean constraint on xi
#'     in the one-level branch; the mean_target constraint applies only to psi
#'     in the two-level branch). We minimise
#'     `.gp_eb_penalized_loglik(alpha, Q, log_P, penalty = c_xi)` directly.
#'   * GENDER: the mean of xi is constrained to 0 (likelihood.m L26-L29 fsolves
#'     `mean_xi == 0`). We optimise over `T - 1` FREE coefficients and rebuild
#'     the full-length `alpha` with `.gp_eb_full_alpha(.., target_mean = 0)`,
#'     whose softmax density has mean exactly 0; the penalty is then
#'     `c_xi * ||full_alpha||` to match Matlab's `c_xi * sqrt(alpha' alpha)`
#'     (where `alpha` there is the FULL, constraint-completed vector).
#'
#' `.gp_eb_penalized_loglik` returns the NEGATIVE penalised log-likelihood
#' (objective to minimise), so we minimise it with `stats::optim` (BFGS), the
#' base-R analogue of Matlab's `fminunc`.
#'
#' @param c_xi Scalar penalty.
#' @param Q Spline basis (M x T).
#' @param log_P Log normal-mixture matrix (N x M) from
#'   `.gp_eb_normal_mixture_matrix`.
#' @param support The support grid (length M).
#' @param characteristic "race" or "gender".
#' @param start Numeric starting value for the FREE coefficients.
#' @param max_iter `stats::optim` iteration cap.
#' @param tol `stats::optim` `reltol`.
#' @return list(g, alpha, objective, converged, convergence_code,
#'   convergence_message).
#' @keywords internal
#' @noRd
.gp_decon_fit_one <- function(c_xi, Q, log_P, support, characteristic,
                              start, max_iter = 5000L, tol = 1e-8) {
  eb_softmax <- .gp_decon_get(".gp_eb_softmax_density")
  eb_pll     <- .gp_decon_get(".gp_eb_penalized_loglik")
  is_gender  <- identical(characteristic, "gender")

  if (is_gender) {
    eb_full <- .gp_decon_get(".gp_eb_full_alpha")
    ## objective over the T-1 FREE coefficients; rebuild full alpha (mean == 0).
    obj <- function(a) {
      full <- eb_full(a, Q, support, target_mean = 0)
      eb_pll(full, Q, log_P, penalty = c_xi)
    }
  } else {
    ## race one-level: alpha is fully free (T coefficients).
    obj <- function(a) eb_pll(a, Q, log_P, penalty = c_xi)
  }

  opt <- stats::optim(
    par = start, fn = obj, method = "BFGS",
    control = list(maxit = as.integer(max_iter), reltol = tol)
  )

  alpha_full <- if (is_gender) {
    .gp_decon_get(".gp_eb_full_alpha")(opt$par, Q, support, target_mean = 0)
  } else {
    opt$par
  }
  sd <- eb_softmax(Q, alpha_full)
  g  <- if (is.list(sd)) sd[["g"]] else sd

  list(g = g, alpha = alpha_full, objective = opt$value,
       converged = isTRUE(opt$convergence == 0L),
       convergence_code = as.integer(opt$convergence),
       convergence_message = if (is.null(opt$message)) NA_character_ else as.character(opt$message))
}

# ---------------------------------------------------------------------------
# Penalty grid selection (the coupling-J rule)
# ---------------------------------------------------------------------------

#' Default one-level penalty grid for a characteristic
#'
#' @description
#' Matches the Matlab one-level grids (estimate_lsqnonlin.m L30-L33 / L40-L43):
#'   race   `seq(0.001, 0.2, by = 0.001)`   (200 nodes)
#'   gender `seq(0.00025, 0.01, by = 0.00025)` (40 nodes)
#' Overridable through `control$deconv_penalty_grid_race/_gender`.
#'
#' @keywords internal
#' @noRd
.gp_decon_penalty_grid <- function(characteristic, control = NULL) {
  if (!is.null(control)) {
    g <- if (identical(characteristic, "race")) {
      control$deconv_penalty_grid_race
    } else {
      control$deconv_penalty_grid_gender
    }
    if (is.numeric(g) && length(g) >= 1L) return(g)
  }
  if (identical(characteristic, "race")) {
    seq(0.001, 0.2, by = 0.001)
  } else {
    seq(0.00025, 0.01, by = 0.00025)
  }
}

# ---------------------------------------------------------------------------
# The native one-level deconvolution
# ---------------------------------------------------------------------------

#' Native one-level deconvolution prior from a coupling fit
#'
#' @description
#' Takes the `gp_estimation_fit` (the single
#' coupling object) and produces the one-level mixing-density prior `gp_prior`
#' on the r/xi-scale, by penalised log-spline maximum likelihood with the
#' penalty selected to minimise the coupling J (the GMM-moment fit of the
#' prior-implied moments). Faithful base-R port of the `group_fx == 0` branch of
#' `estimate_lsqnonlin.m` (see file header for line citations).
#'
#' Invariants: #8 (one-level only -- the native two-level path is a later step;
#' this function refuses anything but a one-level coupling fit). #3 (consume the
#' carriers off the ONE fit -- `v_hat`, `s_v`, `m_hat`, `V_m`, `caps`, `beta`;
#' nothing is re-estimated).
#'
#' @param fit A validated `gp_estimation_fit` (from `gp_w_seam()` /
#'   `gp_estimation_coupling()`).
#' @param control Optional `gp_control`; reads `deconv_supp_pts`,
#'   `deconv_n_knots`, `deconv_penalty_grid_race/_gender`, `deconv_max_iter`,
#'   `tol`, and (for the opt-in multistart) `seed`. Falls back to the Matlab
#'   defaults when absent.
#' @param n_starts Positive integer; number of starts per penalty node.
#'   `n_starts = 1L` (the DEFAULT) uses the single deterministic zeros-start and
#'   reproduces the historical result bit-for-bit. `n_starts > 1L` is an OPT-IN
#'   seeded multistart that, per penalty node, ALSO fits from `n_starts - 1`
#'   additional `rnorm` starts and keeps the lowest-objective fit for that node
#'   (mirroring KRW's per-penalty `randn` restart, `estimate_lsqnonlin.m:356`).
#'   Fully deterministic given the seed (`control$seed`, else a fixed default).
#' @return A validated object of class `"gp_prior"` with slots
#'   `support`, `density` (probability masses summing to 1), `mean`,
#'   `scale = "r"`, `diagnostics` (selected penalty, J, log-likelihood,
#'   convergence, the full penalty/J grid, prior-implied moments vs `m_hat`,
#'   the support caps), `metadata` (characteristic, spline info, beta,
#'   coupling carriers length).
#' @keywords internal
#' @noRd
gp_deconvolve <- function(fit, control = NULL, n_starts = 1L) {
  ## --- validate n_starts (positive integer) -------------------------------
  if (length(n_starts) != 1L || !is.numeric(n_starts) ||
      !is.finite(n_starts) || n_starts < 1 ||
      n_starts != as.integer(n_starts)) {
    .gp_decon_abort(
      sprintf("`n_starts` must be a single positive integer (>= 1); got %s.",
              paste(format(n_starts), collapse = ", ")),
      class = "gradepath_validation_error"
    )
  }
  n_starts <- as.integer(n_starts)

  ## --- validate the coupling fit (invariant #3 contract) ------------------
  validate_fit <- tryCatch(.gp_decon_get("validate_gp_estimation_fit"),
                           error = function(e) NULL)
  if (!is.null(validate_fit)) validate_fit(fit)
  if (!inherits(fit, "gp_estimation_fit")) {
    .gp_decon_abort(
      "`fit` must be a `gp_estimation_fit` (the Step-3.2 coupling object).",
      class = "gradepath_validation_error"
    )
  }

  characteristic <- fit$characteristic
  if (!characteristic %in% c("race", "gender")) {
    .gp_decon_abort(
      sprintf("Unknown `characteristic`: '%s'.", characteristic),
      class = "gradepath_validation_error"
    )
  }

  ## --- one-level only (invariant #8) --------------------------------------
  ## The coupling fit is one-level iff its m_hat has the one-level shape
  ## (race length 2 = c(mu, sigma_xi); gender length 1 = c(sigma_xi)). A
  ## two-level fit would carry an extra sigma_psi moment; refuse it here.
  expected_m <- if (identical(characteristic, "race")) 2L else 1L
  if (length(fit$m_hat) != expected_m) {
    .gp_decon_abort(
      sprintf(paste0("gp_deconvolve() is one-level only (invariant #8): ",
                     "expected m_hat length %d for %s, got %d. The native ",
                     "two-level deconvolution is a later step."),
              expected_m, characteristic, length(fit$m_hat)),
      class = "gradepath_validation_error"
    )
  }

  ## --- control knobs (Matlab defaults when absent) ------------------------
  supp_pts <- if (!is.null(control$deconv_supp_pts)) {
    as.integer(control$deconv_supp_pts)
  } else 1000L
  n_knots <- if (!is.null(control$deconv_n_knots)) {
    as.integer(control$deconv_n_knots)
  } else 5L
  max_iter <- if (!is.null(control$deconv_max_iter)) {
    as.integer(control$deconv_max_iter)
  } else 5000L
  opt_tol <- if (!is.null(control$tol)) control$tol else 1e-8

  ## --- carriers off the ONE fit (invariant #3) ----------------------------
  v_hat <- fit$v_hat
  s_v   <- fit$s_v
  m_hat <- fit$m_hat
  V_m   <- fit$V_m
  caps  <- fit$caps

  ## --- support grid + spline basis ----------------------------------------
  support <- .gp_decon_support(caps, supp_pts = supp_pts)
  eb_basis  <- .gp_decon_get(".gp_eb_spline_basis")
  eb_mixmat <- .gp_decon_get(".gp_eb_normal_mixture_matrix")
  Q <- eb_basis(support, n_knots = n_knots)
  Tdim <- ncol(Q)

  ## normal-mixture matrix on the carriers (the standardized residuals on the
  ## r-scale and their scaled standard errors). likelihood.m one-level reduces
  ## to P_tilde_i = sum_m g_m * normpdf((v_i - psi*supp_m)/s_v_i)/s_v_i with
  ## psi degenerate (=1 race, =0 gender) -- i.e. the standard normal-mixture
  ## convolution of the prior with the per-firm Gaussian noise. That is exactly
  ## the `.gp_eb_normal_mixture_matrix(theta_hat, s, support)` log-likelihood
  ## kernel.
  mm <- eb_mixmat(v_hat, s_v, support, log = TRUE)
  log_P <- if (is.list(mm)) mm[["log_P"]] else mm

  ## --- penalty grid -------------------------------------------------------
  pen_grid <- .gp_decon_penalty_grid(characteristic, control)
  is_gender <- identical(characteristic, "gender")
  n_free <- if (is_gender) Tdim - 1L else Tdim

  ## Starts for the per-penalty optimisation.
  ##   * DEFAULT (`n_starts == 1L`): a SINGLE deterministic zeros-start. The
  ##     penalised log-spline objective is non-convex; KRW's reference draws a
  ##     FRESH `randn` start per penalty node (estimate_lsqnonlin.m:356, inside
  ##     `parfor` + `fminunc`). We substitute the zeros-start here purely for
  ##     replication determinism. This is the historical behaviour and is what
  ##     `gp_deconvolve(fit)` (no `n_starts`) reproduces bit-for-bit.
  ##   * OPT-IN (`n_starts > 1L`): start #1 is STILL the zeros vector, plus
  ##     `n_starts - 1` additional SEEDED `rnorm` starts per penalty node,
  ##     mirroring KRW's per-penalty `randn` restart. For each node we keep the
  ##     lowest-objective fit, then run the existing coupling-J penalty
  ##     selection unchanged. Because start #1 is always zeros, the multistart
  ##     can only find an equal-or-lower objective at each node -- it never
  ##     degrades the default. The RNG is seeded ONCE here (from `control$seed`
  ##     when present, else a fixed default of 1L) so the whole multistart is
  ##     fully reproducible. The draws are precomputed up front so the seeding
  ##     is independent of how many penalty nodes happen to be skipped.
  start <- rep(0, n_free)
  extra_starts <- if (n_starts > 1L) {
    seed <- if (!is.null(control$seed)) as.integer(control$seed) else 1L
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    ## A list (length nP) of matrices, each (n_starts - 1) x n_free, of rnorm
    ## restart offsets -- one independent block of draws per penalty node, so a
    ## skipped/failed node does not shift the draws used by later nodes.
    lapply(seq_len(length(pen_grid)), function(.) {
      matrix(stats::rnorm((n_starts - 1L) * n_free), nrow = n_starts - 1L,
             ncol = n_free)
    })
  } else NULL

  ## --- loop over penalties, score by coupling J ---------------------------
  nP <- length(pen_grid)
  J_vec    <- rep(NA_real_, nP)
  obj_vec  <- rep(NA_real_, nP)
  conv_vec <- rep(FALSE, nP)
  conv_code_vec <- rep(NA_integer_, nP)
  conv_msg_vec  <- rep(NA_character_, nP)
  best <- list(J = Inf, idx = NA_integer_, g = NULL, alpha = NULL,
               objective = NA_real_, mom = NULL, converged = FALSE,
               convergence_code = NA_integer_,
               convergence_message = NA_character_)

  fit_node <- function(c_xi, k) {
    ## start #1: the deterministic zeros vector (always present).
    f <- tryCatch(
      .gp_decon_fit_one(c_xi, Q, log_P, support, characteristic,
                        start = start, max_iter = max_iter, tol = opt_tol),
      error = function(e) NULL
    )
    if (!is.null(f) && any(!is.finite(f$g))) f <- NULL
    ## opt-in: additional seeded rnorm starts; keep the lowest objective.
    if (!is.null(extra_starts)) {
      offs <- extra_starts[[k]]
      for (j in seq_len(nrow(offs))) {
        st_j <- start + offs[j, ]
        fj <- tryCatch(
          .gp_decon_fit_one(c_xi, Q, log_P, support, characteristic,
                            start = st_j, max_iter = max_iter, tol = opt_tol),
          error = function(e) NULL
        )
        if (is.null(fj) || any(!is.finite(fj$g))) next
        if (is.null(f) || (is.finite(fj$objective) &&
                           fj$objective < f$objective)) {
          f <- fj
        }
      }
    }
    f
  }

  for (k in seq_len(nP)) {
    c_xi <- pen_grid[k]
    f <- fit_node(c_xi, k)
    if (is.null(f) || any(!is.finite(f$g))) next
    mom <- .gp_decon_model_moments(support, f$g)
    Jk  <- .gp_decon_J(m_hat, V_m, mom, characteristic)
    J_vec[k]    <- Jk
    obj_vec[k]  <- f$objective
    conv_vec[k] <- f$converged
    conv_code_vec[k] <- f$convergence_code
    conv_msg_vec[k]  <- f$convergence_message
    if (is.finite(Jk) && Jk < best$J) {
      best <- list(J = Jk, idx = k, g = f$g, alpha = f$alpha,
                   objective = f$objective, mom = mom,
                   converged = f$converged,
                   convergence_code = f$convergence_code,
                   convergence_message = f$convergence_message)
    }
  }

  if (is.na(best$idx)) {
    .gp_decon_abort(
      "Deconvolution failed: no penalty on the grid produced a finite fit.",
      class = "gradepath_estimation_error"
    )
  }
  if (!isTRUE(best$converged)) {
    warning(sprintf(
      "Selected deconvolution fit did not converge (optim convergence code %s): %s",
      best$convergence_code,
      if (is.na(best$convergence_message)) "no optimizer message" else best$convergence_message
    ), call. = FALSE)
  }

  ## --- final density: normalise to probability masses (g_xi / sum) --------
  ## estimate_lsqnonlin.m L469: g_xi = g_xi / sum(g_xi). The softmax already
  ## sums to 1, but we renormalise defensively to satisfy the gp_prior contract
  ## exactly (validator requires |sum - 1| <= 1e-6).
  density <- best$g / sum(best$g)
  prior_mean <- sum(support * density)

  ## --- assemble the gp_prior ----------------------------------------------
  new_prior      <- .gp_decon_get("new_gp_prior")
  validate_prior <- .gp_decon_get("validate_gp_prior")

  diagnostics <- list(
    method            = "native-one-level-spline-softmax",
    penalty           = pen_grid[best$idx],
    penalty_index     = best$idx,
    J                 = best$J,
    log_likelihood    = -best$objective,   # objective is NEGATIVE penalised LL
    objective         = best$objective,
    converged         = best$converged,
    convergence_code  = best$convergence_code,
    convergence_message = best$convergence_message,
    penalty_grid      = pen_grid,
    J_grid            = J_vec,
    objective_grid    = obj_vec,
    converged_grid    = conv_vec,
    convergence_code_grid = conv_code_vec,
    convergence_message_grid = conv_msg_vec,
    model_moments     = best$mom,          # c(mean, sd) of the prior
    coupling_m_hat    = m_hat,
    caps              = caps,
    n_knots           = n_knots,
    spline_dim        = Tdim
  )
  metadata <- list(
    characteristic = characteristic,
    beta           = fit$beta,
    n_carriers     = length(v_hat),
    supp_pts       = supp_pts,
    alpha          = best$alpha
  )

  prior <- new_prior(
    support     = support,
    density     = density,
    mean        = prior_mean,
    scale       = "r",
    diagnostics = diagnostics,
    metadata    = metadata
  )
  validate_prior(prior)
}
