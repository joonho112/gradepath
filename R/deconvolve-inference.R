# =============================================================================
# deconvolve-inference.R -- deconvolution inference primitives
# -----------------------------------------------------------------------------
# Native KRW inference layer around the deconvolution spline coefficients:
# finite-difference Hessians, sandwich covariance, delta-method SEs, and a seeded
# parametric-bootstrap harness. The binding Matlab lines are estimate_lsqnonlin.m
# 425-549 and delta_method.m.
# =============================================================================

#' Resolve package abort helper whether sourced or namespaced
#' @keywords internal
#' @noRd
.gp_inf_abort <- function(msg, class = "gradepath_error") {
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

#' Resolve a package internal by name
#' @keywords internal
#' @noRd
.gp_inf_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  .gp_inf_abort(
    sprintf("Required gradepath internal '%s' not found.", name),
    class = "gradepath_internal_error"
  )
}

#' @keywords internal
#' @noRd
.gp_inf_or <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
#' @noRd
.gp_inf_par <- function(par, name = "par") {
  par <- as.numeric(par)
  if (length(par) < 1L || any(!is.finite(par))) {
    .gp_inf_abort(sprintf("`%s` must be a non-empty finite numeric vector.", name),
                  class = "gradepath_validation_error")
  }
  par
}

#' @keywords internal
#' @noRd
.gp_inf_step <- function(par, step = NULL, name = "step") {
  if (is.null(step)) {
    return((.Machine$double.eps)^(1 / 4) * pmax(abs(par), 1))
  }
  step <- as.numeric(step)
  if (length(step) == 1L) step <- rep(step, length(par))
  if (length(step) != length(par) || any(!is.finite(step)) || any(step <= 0)) {
    .gp_inf_abort(
      sprintf("`%s` must be a positive finite scalar or length(par) vector.", name),
      class = "gradepath_validation_error"
    )
  }
  step
}

#' @keywords internal
#' @noRd
.gp_inf_eval_scalar <- function(fn, par, ...) {
  val <- fn(par, ...)
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
    .gp_inf_abort("Inference objective functions must return one finite numeric value.",
                  class = "gradepath_validation_error")
  }
  as.numeric(val)
}

#' @keywords internal
#' @noRd
.gp_inf_eval_vector <- function(fn, par, ...) {
  val <- fn(par, ...)
  if (!is.numeric(val) || length(val) < 1L || any(!is.finite(val))) {
    .gp_inf_abort("Delta-method target functions must return a finite numeric vector.",
                  class = "gradepath_validation_error")
  }
  out <- as.numeric(val)
  names(out) <- names(val)
  out
}

#' Finite-difference Hessian for a scalar objective
#'
#' Central finite differences at `par`, used to mirror Matlab's `fminunc(...,
#' MaxIter=0)` Hessian extraction for penalized and unpenalized objectives.
#' @keywords internal
#' @noRd
gp_hessian <- function(object, par = NULL, input = NULL, ..., step = NULL,
                       fit = NULL, penalized = TRUE) {
  if (!is.function(object)) {
    fit_arg <- if (is.null(fit)) par else fit
    return(gp_deconvolution_hessian(object, fit_arg, input,
                                    penalized = penalized, step = step))
  }
  fn <- object
  if (!is.function(fn)) {
    .gp_inf_abort("`fn` must be a function.", class = "gradepath_validation_error")
  }
  par <- .gp_inf_par(par)
  h <- .gp_inf_step(par, step)
  p <- length(par)
  H <- matrix(0, nrow = p, ncol = p)
  f0 <- .gp_inf_eval_scalar(fn, par, ...)

  for (i in seq_len(p)) {
    ei <- numeric(p)
    ei[i] <- h[i]

    f_plus <- .gp_inf_eval_scalar(fn, par + ei, ...)
    f_minus <- .gp_inf_eval_scalar(fn, par - ei, ...)
    H[i, i] <- (f_plus - 2 * f0 + f_minus) / (h[i]^2)

    if (i < p) {
      for (j in (i + 1L):p) {
        ej <- numeric(p)
        ej[j] <- h[j]
        f_pp <- .gp_inf_eval_scalar(fn, par + ei + ej, ...)
        f_pm <- .gp_inf_eval_scalar(fn, par + ei - ej, ...)
        f_mp <- .gp_inf_eval_scalar(fn, par - ei + ej, ...)
        f_mm <- .gp_inf_eval_scalar(fn, par - ei - ej, ...)
        H[i, j] <- (f_pp - f_pm - f_mp + f_mm) / (4 * h[i] * h[j])
        H[j, i] <- H[i, j]
      }
    }
  }

  dimnames(H) <- list(names(par), names(par))
  H
}

#' @keywords internal
#' @noRd
.gp_jacobian <- function(fn, par, step = NULL, ...) {
  par <- .gp_inf_par(par)
  h <- .gp_inf_step(par, step)
  f0 <- .gp_inf_eval_vector(fn, par, ...)
  p <- length(par)
  J <- matrix(NA_real_, nrow = length(f0), ncol = p)

  for (j in seq_len(p)) {
    ej <- numeric(p)
    ej[j] <- h[j]
    fp <- .gp_inf_eval_vector(fn, par + ej, ...)
    fm <- .gp_inf_eval_vector(fn, par - ej, ...)
    if (length(fp) != length(f0) || length(fm) != length(f0)) {
      .gp_inf_abort("Delta-method target length changed under finite differences.",
                    class = "gradepath_validation_error")
    }
    J[, j] <- (fp - fm) / (2 * h[j])
  }

  rownames(J) <- names(f0)
  colnames(J) <- names(par)
  J
}

#' Penalized/unpenalized sandwich covariance for deconvolution coefficients
#'
#' Computes `H1 = Hessian(penalized_fn)`, `H2 = Hessian(unpenalized_fn)`, and
#' `V = solve(H1) %*% H2 %*% solve(H1)`, matching estimate_lsqnonlin.m:425-434.
#' @keywords internal
#' @noRd
gp_sandwich_vcov <- function(par, penalized_fn, unpenalized_fn, step = NULL) {
  if (!is.function(penalized_fn) || !is.function(unpenalized_fn)) {
    .gp_inf_abort("`penalized_fn` and `unpenalized_fn` must both be functions.",
                  class = "gradepath_validation_error")
  }
  par <- .gp_inf_par(par)
  H1 <- gp_hessian(penalized_fn, par, step = step)
  H2 <- gp_hessian(unpenalized_fn, par, step = step)
  rc <- tryCatch(rcond(H1), error = function(e) 0)
  if (!is.finite(rc) || rc < 1e-12) {
    .gp_inf_abort(sprintf("Penalized Hessian is numerically singular (rcond = %.3g).", rc),
                  class = "gradepath_singular_error")
  }
  H1_inv <- solve(H1)
  V <- H1_inv %*% H2 %*% H1_inv
  V <- (V + t(V)) / 2
  dimnames(V) <- dimnames(H1)
  list(H1 = H1, H2 = H2, vcov = V, rcond_H1 = rc)
}

#' Delta-method standard errors for functions of deconvolution coefficients
#'
#' Computes `f_hat`, the finite-difference Jacobian `df`, and
#' `sqrt(diag(df %*% vcov %*% t(df)))`, matching delta_method.m plus
#' estimate_lsqnonlin.m:444-455.
#' @keywords internal
#' @noRd
gp_delta_method_se <- function(object, fn = NULL, vcov = NULL, input = NULL,
                               ..., step = NULL, f_hat = NULL, fit = NULL,
                               supp_pts_theta = 250L) {
  if (!is.numeric(object)) {
    fit_arg <- if (is.null(fit)) vcov else fit
    return(gp_deconvolution_delta_method_se(
      prior = object,
      vcov = fn,
      fit = fit_arg,
      input = input,
      step = step,
      supp_pts_theta = supp_pts_theta
    ))
  }
  par <- object
  if (!is.function(fn)) {
    .gp_inf_abort("`fn` must be a function.", class = "gradepath_validation_error")
  }
  par <- .gp_inf_par(par)
  if (!is.matrix(vcov) || nrow(vcov) != ncol(vcov) || nrow(vcov) != length(par) ||
      any(!is.finite(vcov))) {
    .gp_inf_abort("`vcov` must be a finite square matrix with dimension length(par).",
                  class = "gradepath_validation_error")
  }
  f0 <- if (is.null(f_hat)) .gp_inf_eval_vector(fn, par, ...) else as.numeric(f_hat)
  if (length(f0) < 1L || any(!is.finite(f0))) {
    .gp_inf_abort("`f_hat` must be a finite numeric vector when supplied.",
                  class = "gradepath_validation_error")
  }
  J <- .gp_jacobian(fn, par, step = step, ...)
  if (nrow(J) != length(f0)) {
    .gp_inf_abort("Jacobian row count does not match `f_hat` length.",
                  class = "gradepath_validation_error")
  }
  cov_f <- J %*% vcov %*% t(J)
  cov_f <- (cov_f + t(cov_f)) / 2
  se <- sqrt(pmax(0, diag(cov_f)))
  names(se) <- names(f0)
  list(value = f0, se = se, covariance = cov_f, jacobian = J)
}

#' @keywords internal
#' @noRd
.gp_inf_distribution_moments <- function(support, g, prefix) {
  support <- as.numeric(support)
  g <- as.numeric(g)
  if (length(support) != length(g) || length(g) < 2L || any(!is.finite(support)) ||
      any(!is.finite(g)) || any(g < 0) || sum(g) <= 0) {
    .gp_inf_abort(sprintf("Invalid support/mass pair for `%s` moments.", prefix),
                  class = "gradepath_validation_error")
  }
  g <- g / sum(g)
  mu <- sum(support * g)
  sd <- sqrt(sum(((support - mu)^2) * g))
  if (!is.finite(sd) || sd <= 0) {
    .gp_inf_abort(sprintf("`%s` distribution has zero variance.", prefix),
                  class = "gradepath_validation_error")
  }
  z <- (support - mu) / sd
  out <- c(
    mean = mu,
    sd = sd,
    skew = sum((z^3) * g),
    excess_kurt = sum((z^4) * g) - 3
  )
  names(out) <- paste(prefix, names(out), sep = "_")
  out
}

#' @keywords internal
#' @noRd
.gp_inf_split_twolevel_free <- function(par, Q_xi, Q_eta, characteristic) {
  par <- .gp_inf_par(par)
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (!is.matrix(Q_xi) || !is.matrix(Q_eta)) {
    .gp_inf_abort("`Q_xi` and `Q_eta` must be matrices.",
                  class = "gradepath_validation_error")
  }
  T_xi <- ncol(Q_xi)
  T_eta <- ncol(Q_eta)
  n_xi_free <- if (identical(characteristic, "race")) T_xi else T_xi - 1L
  n_eta_free <- T_eta - 1L
  n_free <- n_xi_free + n_eta_free
  if (length(par) != n_free) {
    .gp_inf_abort(sprintf("Free parameter length %d != expected two-level length %d.",
                          length(par), n_free),
                  class = "gradepath_validation_error")
  }
  list(
    alpha_xi_free = par[seq_len(n_xi_free)],
    alpha_eta_free = par[n_xi_free + seq_len(n_eta_free)],
    n_xi_free = n_xi_free,
    n_eta_free = n_eta_free
  )
}

#' KRW delta-method targets from a packed two-level free-parameter vector
#'
#' Returns the `delta_method.m` functions of interest for `group_fx == 1`: mean,
#' standard deviation, skewness, and excess kurtosis for xi, eta/psi, and theta.
#' @keywords internal
#' @noRd
gp_delta_targets <- function(par, Q_xi, supp_xi, Q_eta, supp_eta,
                             s, mu = 0, beta, characteristic,
                             supp_pts_theta = 250L) {
  characteristic <- match.arg(characteristic, c("race", "gender"))
  split <- .gp_inf_split_twolevel_free(par, Q_xi, Q_eta, characteristic)
  if (nrow(Q_xi) != length(supp_xi) || nrow(Q_eta) != length(supp_eta)) {
    .gp_inf_abort("Support lengths must match the row counts of `Q_xi` and `Q_eta`.",
                  class = "gradepath_validation_error")
  }
  s <- as.numeric(s)
  if (length(s) < 1L || any(!is.finite(s)) || any(s <= 0)) {
    .gp_inf_abort("`s` must be a non-empty positive finite numeric vector.",
                  class = "gradepath_validation_error")
  }
  if (length(beta) != 1L || !is.finite(beta) || length(mu) != 1L || !is.finite(mu)) {
    .gp_inf_abort("`mu` and `beta` must be finite scalars.",
                  class = "gradepath_validation_error")
  }

  eb_softmax <- .gp_inf_get(".gp_eb_softmax_density")
  eb_full <- .gp_inf_get(".gp_eb_full_alpha")

  alpha_xi <- if (identical(characteristic, "race")) {
    split$alpha_xi_free
  } else {
    eb_full(split$alpha_xi_free, Q_xi, supp_xi, target_mean = 0)
  }
  eta_target <- if (identical(characteristic, "race")) 1 else 0
  alpha_eta <- eb_full(split$alpha_eta_free, Q_eta, supp_eta,
                       target_mean = eta_target)

  sx <- eb_softmax(Q_xi, alpha_xi)
  se <- eb_softmax(Q_eta, alpha_eta)
  g_xi <- if (is.list(sx)) sx$g else sx
  g_eta <- if (is.list(se)) se$g else se
  theta <- .gp_inf_get("gp_pushforward_theta")(
    supp_xi = supp_xi,
    g_xi = g_xi,
    supp_eta = supp_eta,
    g_eta = g_eta,
    s = s,
    mu = mu,
    beta = beta,
    characteristic = characteristic,
    supp_pts_theta = supp_pts_theta
  )

  c(
    .gp_inf_distribution_moments(supp_xi, g_xi, "xi"),
    .gp_inf_distribution_moments(supp_eta, g_eta, "eta"),
    .gp_inf_distribution_moments(theta$support, theta$g, "theta")
  )
}

#' @keywords internal
#' @noRd
.gp_inf_decon_context <- function(prior, fit, input, supp_pts_theta = 250L) {
  if (is.null(input) || is.null(input$s)) {
    .gp_inf_abort("`input` with original standard errors `s` is required.",
                  class = "gradepath_validation_error")
  }
  characteristic <- .gp_inf_or(fit$characteristic, prior$metadata$characteristic)
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (is.null(prior$diagnostics$support_eta) || is.null(prior$diagnostics$g_eta) ||
      is.null(prior$diagnostics$c_xi) || is.null(prior$diagnostics$c_eta) ||
      is.null(prior$metadata$alpha_xi_free) || is.null(prior$metadata$alpha_eta_free)) {
    .gp_inf_abort("`prior` must be a Step-7.2 two-level deconvolution prior.",
                  class = "gradepath_validation_error")
  }
  if (is.null(fit$v_hat) || is.null(fit$s_v) || is.null(fit$beta)) {
    .gp_inf_abort("`fit` must carry `v_hat`, `s_v`, and `beta`.",
                  class = "gradepath_validation_error")
  }
  raw_industry <- .gp_inf_or(input$industry, fit$industry)
  if (is.null(raw_industry)) {
    .gp_inf_abort("`input$industry` or `fit$industry` is required.",
                  class = "gradepath_validation_error")
  }
  industry_idx <- as.integer(match(raw_industry, sort(unique(raw_industry), method = "radix")))
  if (length(industry_idx) != length(fit$v_hat) || length(input$s) != length(fit$v_hat)) {
    .gp_inf_abort("`input$s`, `industry`, and `fit$v_hat` must be length-matched.",
                  class = "gradepath_validation_error")
  }

  n_knots <- .gp_inf_or(prior$diagnostics$n_knots, 5L)
  eb_basis <- .gp_inf_get(".gp_eb_spline_basis")
  Q_xi <- eb_basis(prior$support, n_knots = n_knots)
  Q_eta <- eb_basis(prior$diagnostics$support_eta, n_knots = n_knots)
  par <- c(prior$metadata$alpha_xi_free, prior$metadata$alpha_eta_free)
  split <- .gp_inf_split_twolevel_free(par, Q_xi, Q_eta, characteristic)
  mu <- .gp_inf_or(fit$mu, if (identical(characteristic, "race") &&
                              !is.null(fit$m_hat)) fit$m_hat[1L] else 0)

  make_obj <- function(c_xi, c_eta) {
    force(c_xi)
    force(c_eta)
    function(p) {
      sp <- .gp_inf_split_twolevel_free(p, Q_xi, Q_eta, characteristic)
      gp_two_level_likelihood(
        sp$alpha_xi_free, sp$alpha_eta_free,
        r = as.numeric(fit$v_hat),
        s_v = as.numeric(fit$s_v),
        industry_idx = industry_idx,
        Q_xi = Q_xi,
        supp_xi = prior$support,
        Q_eta = Q_eta,
        supp_eta = prior$diagnostics$support_eta,
        c_xi = c_xi,
        c_eta = c_eta,
        characteristic = characteristic,
        return_densities = FALSE
      )
    }
  }
  target_fn <- function(p) {
    gp_delta_targets(
      p, Q_xi = Q_xi, supp_xi = prior$support,
      Q_eta = Q_eta, supp_eta = prior$diagnostics$support_eta,
      s = input$s,
      mu = mu,
      beta = fit$beta,
      characteristic = characteristic,
      supp_pts_theta = supp_pts_theta
    )
  }

  list(
    par = par,
    split = split,
    characteristic = characteristic,
    c_xi = prior$diagnostics$c_xi,
    c_eta = prior$diagnostics$c_eta,
    objective_penalized = make_obj(prior$diagnostics$c_xi, prior$diagnostics$c_eta),
    objective_unpenalized = make_obj(0, 0),
    target_fn = target_fn
  )
}

#' Contract-shaped Hessian for a two-level deconvolution prior
#' @keywords internal
#' @noRd
gp_deconvolution_hessian <- function(prior, fit, input, penalized = TRUE,
                                     step = NULL) {
  ctx <- .gp_inf_decon_context(prior, fit, input)
  obj <- if (isTRUE(penalized)) ctx$objective_penalized else ctx$objective_unpenalized
  gp_hessian(obj, ctx$par, step = step)
}

#' Contract-shaped sandwich covariance for a two-level deconvolution prior
#' @keywords internal
#' @noRd
gp_deconvolution_sandwich <- function(prior, fit, input, step = NULL) {
  ctx <- .gp_inf_decon_context(prior, fit, input)
  gp_sandwich_vcov(
    par = ctx$par,
    penalized_fn = ctx$objective_penalized,
    unpenalized_fn = ctx$objective_unpenalized,
    step = step
  )
}

#' Contract-shaped delta-method SEs for a two-level deconvolution prior
#' @keywords internal
#' @noRd
gp_deconvolution_delta_method_se <- function(prior, vcov, fit, input,
                                             step = NULL,
                                             supp_pts_theta = 250L) {
  ctx <- .gp_inf_decon_context(prior, fit, input,
                               supp_pts_theta = supp_pts_theta)
  gp_delta_method_se(ctx$par, ctx$target_fn, vcov, step = step)
}

#' Deconvolution inference bundle for a two-level prior and GMM fit
#'
#' Reconstructs the selected deconvolution objective at `alpha_hat`, computes the
#' penalized/unpenalized sandwich covariance, then applies the KRW delta targets.
#' @keywords internal
#' @noRd
gp_deconvolution_inference <- function(prior, fit, input,
                                       hessian_step = NULL,
                                       jacobian_step = NULL,
                                       supp_pts_theta = 250L) {
  if (is.null(input) || is.null(input$s)) {
    .gp_inf_abort("`input` with original standard errors `s` is required.",
                  class = "gradepath_validation_error")
  }
  characteristic <- .gp_inf_or(fit$characteristic, prior$metadata$characteristic)
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (is.null(prior$diagnostics$support_eta) || is.null(prior$diagnostics$g_eta) ||
      is.null(prior$diagnostics$c_xi) || is.null(prior$diagnostics$c_eta) ||
      is.null(prior$metadata$alpha_xi_free) || is.null(prior$metadata$alpha_eta_free)) {
    .gp_inf_abort("`prior` must be a Step-7.2 two-level deconvolution prior.",
                  class = "gradepath_validation_error")
  }
  if (is.null(fit$v_hat) || is.null(fit$s_v) || is.null(fit$beta)) {
    .gp_inf_abort("`fit` must carry `v_hat`, `s_v`, and `beta`.",
                  class = "gradepath_validation_error")
  }
  raw_industry <- .gp_inf_or(input$industry, fit$industry)
  if (is.null(raw_industry)) {
    .gp_inf_abort("`input$industry` or `fit$industry` is required.",
                  class = "gradepath_validation_error")
  }
  industry_idx <- as.integer(match(raw_industry, sort(unique(raw_industry), method = "radix")))
  if (length(industry_idx) != length(fit$v_hat) || length(input$s) != length(fit$v_hat)) {
    .gp_inf_abort("`input$s`, `industry`, and `fit$v_hat` must be length-matched.",
                  class = "gradepath_validation_error")
  }

  n_knots <- .gp_inf_or(prior$diagnostics$n_knots, 5L)
  eb_basis <- .gp_inf_get(".gp_eb_spline_basis")
  Q_xi <- eb_basis(prior$support, n_knots = n_knots)
  Q_eta <- eb_basis(prior$diagnostics$support_eta, n_knots = n_knots)
  par <- c(prior$metadata$alpha_xi_free, prior$metadata$alpha_eta_free)
  split <- .gp_inf_split_twolevel_free(par, Q_xi, Q_eta, characteristic)

  make_obj <- function(c_xi, c_eta) {
    force(c_xi)
    force(c_eta)
    function(p) {
      sp <- .gp_inf_split_twolevel_free(p, Q_xi, Q_eta, characteristic)
      gp_two_level_likelihood(
        sp$alpha_xi_free, sp$alpha_eta_free,
        r = as.numeric(fit$v_hat),
        s_v = as.numeric(fit$s_v),
        industry_idx = industry_idx,
        Q_xi = Q_xi,
        supp_xi = prior$support,
        Q_eta = Q_eta,
        supp_eta = prior$diagnostics$support_eta,
        c_xi = c_xi,
        c_eta = c_eta,
        characteristic = characteristic,
        return_densities = FALSE
      )
    }
  }

  sand <- gp_sandwich_vcov(
    par = par,
    penalized_fn = make_obj(prior$diagnostics$c_xi, prior$diagnostics$c_eta),
    unpenalized_fn = make_obj(0, 0),
    step = hessian_step
  )
  target_fn <- function(p) {
    gp_delta_targets(
      p, Q_xi = Q_xi, supp_xi = prior$support,
      Q_eta = Q_eta, supp_eta = prior$diagnostics$support_eta,
      s = input$s,
      mu = .gp_inf_or(fit$mu, if (identical(characteristic, "race")) fit$m_hat[1L] else 0),
      beta = fit$beta,
      characteristic = characteristic,
      supp_pts_theta = supp_pts_theta
    )
  }
  delta <- gp_delta_method_se(par, target_fn, sand$vcov, step = jacobian_step)
  structure(
    list(
      parameters = par,
      H1 = sand$H1,
      H2 = sand$H2,
      vcov = sand$vcov,
      value = delta$value,
      se = delta$se,
      covariance = delta$covariance,
      jacobian = delta$jacobian,
      diagnostics = list(
        method = "native-deconvolution-inference",
        characteristic = characteristic,
        c_xi = prior$diagnostics$c_xi,
        c_eta = prior$diagnostics$c_eta,
        rcond_H1 = sand$rcond_H1,
        n_xi_free = split$n_xi_free,
        n_eta_free = split$n_eta_free
      )
    ),
    class = c("gp_deconvolution_inference", "list")
  )
}

#' @keywords internal
#' @noRd
.gp_param_bootstrap_loop <- function(simulate, estimate, statistic,
                                     n_rep = 500L, seed = 15238L,
                                     keep_draws = TRUE) {
  if (!is.function(simulate) || !is.function(estimate) || !is.function(statistic)) {
    .gp_inf_abort("`simulate`, `estimate`, and `statistic` must be functions.",
                  class = "gradepath_validation_error")
  }
  if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1 ||
      n_rep != as.integer(n_rep)) {
    .gp_inf_abort("`n_rep` must be a positive integer.",
                  class = "gradepath_validation_error")
  }
  n_rep <- as.integer(n_rep)
  if (length(seed) != 1L || !is.finite(seed)) {
    .gp_inf_abort("`seed` must be a finite scalar.",
                  class = "gradepath_validation_error")
  }
  seed <- as.integer(seed)

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

  draws <- NULL
  stat_names <- NULL
  for (b in seq_len(n_rep)) {
    sim <- simulate(b)
    est <- estimate(sim, b)
    stat_raw <- statistic(est, b, sim)
    stat_names_b <- names(stat_raw)
    stat <- as.numeric(stat_raw)
    if (length(stat) < 1L || any(!is.finite(stat))) {
      .gp_inf_abort("`statistic` must return a finite numeric vector.",
                    class = "gradepath_validation_error")
    }
    if (is.null(draws)) {
      draws <- matrix(NA_real_, nrow = n_rep, ncol = length(stat))
      stat_names <- stat_names_b
      if (!is.null(stat_names)) colnames(draws) <- stat_names
    } else if (ncol(draws) != length(stat)) {
      .gp_inf_abort("Bootstrap statistic length changed across replications.",
                    class = "gradepath_validation_error")
    }
    draws[b, ] <- stat
  }

  se <- if (n_rep == 1L) rep(NA_real_, ncol(draws)) else apply(draws, 2L, stats::sd)
  mean <- colMeans(draws)
  if (!is.null(stat_names)) {
    names(mean) <- stat_names
    names(se) <- stat_names
  }
  structure(
    list(
      mean = mean,
      se = se,
      draws = if (isTRUE(keep_draws)) draws else NULL,
      n_rep = n_rep,
      seed = seed,
      diagnostics = list(method = "seeded-parametric-bootstrap")
    ),
    class = c("gp_param_bootstrap", "list")
  )
}

#' Seeded two-level parametric bootstrap
#'
#' Reuses the selected support, basis, and penalties from the two-level prior,
#' simulates eta by industry and xi by firm, refits only the selected penalty
#' node, and returns bootstrap standard errors for `gp_delta_targets()`. A generic
#' closure-driven mode is retained for small tests by supplying `simulate`,
#' `estimate`, and `statistic`.
#' @keywords internal
#' @noRd
gp_param_bootstrap <- function(prior = NULL, fit = NULL, input = NULL,
                               n_rep = 500L, B = NULL, seed = 15238L,
                               keep_draws = TRUE,
                               max_iter = 500L, tol = 1e-6,
                               supp_pts_theta = 250L,
                               simulate = NULL, estimate = NULL,
                               statistic = NULL) {
  if (!is.null(B)) {
    if (!missing(n_rep) && !identical(as.numeric(n_rep), as.numeric(B))) {
      .gp_inf_abort("Specify only one of `n_rep` or `B`, or give them the same value.",
                    class = "gradepath_validation_error")
    }
    n_rep <- B
  }
  if (is.function(simulate) || is.function(estimate) || is.function(statistic)) {
    return(.gp_param_bootstrap_loop(simulate, estimate, statistic,
                                    n_rep = n_rep, seed = seed,
                                    keep_draws = keep_draws))
  }
  if (is.null(prior) || is.null(fit) || is.null(input)) {
    .gp_inf_abort("`prior`, `fit`, and `input` are required for two-level bootstrap mode.",
                  class = "gradepath_validation_error")
  }
  characteristic <- match.arg(.gp_inf_or(fit$characteristic, prior$metadata$characteristic),
                              c("race", "gender"))
  if (is.null(input$s) || is.null(fit$s_v) || is.null(fit$v_hat) ||
      is.null(prior$diagnostics$support_eta) || is.null(prior$diagnostics$g_eta) ||
      is.null(prior$diagnostics$c_xi) || is.null(prior$diagnostics$c_eta)) {
    .gp_inf_abort("Two-level bootstrap needs input$s, fit$v_hat/s_v, and Step-7.2 prior diagnostics.",
                  class = "gradepath_validation_error")
  }
  raw_industry <- .gp_inf_or(input$industry, fit$industry)
  if (is.null(raw_industry)) {
    .gp_inf_abort("Two-level bootstrap needs industry membership.",
                  class = "gradepath_validation_error")
  }
  industry_idx <- as.integer(match(raw_industry, sort(unique(raw_industry), method = "radix")))
  N <- length(fit$v_hat)
  if (length(industry_idx) != N || length(fit$s_v) != N || length(input$s) != N) {
    .gp_inf_abort("`input$s`, `industry`, `fit$v_hat`, and `fit$s_v` must be length-matched.",
                  class = "gradepath_validation_error")
  }

  n_knots <- .gp_inf_or(prior$diagnostics$n_knots, 5L)
  eb_basis <- .gp_inf_get(".gp_eb_spline_basis")
  Q_xi <- eb_basis(prior$support, n_knots = n_knots)
  Q_eta <- eb_basis(prior$diagnostics$support_eta, n_knots = n_knots)
  par_hat <- c(prior$metadata$alpha_xi_free, prior$metadata$alpha_eta_free)
  split_hat <- .gp_inf_split_twolevel_free(par_hat, Q_xi, Q_eta, characteristic)
  n_free <- length(par_hat)

  g_xi <- as.numeric(prior$density)
  g_xi <- g_xi / sum(g_xi)
  g_eta <- as.numeric(prior$diagnostics$g_eta)
  g_eta <- g_eta / sum(g_eta)
  supp_xi <- as.numeric(prior$support)
  supp_eta <- as.numeric(prior$diagnostics$support_eta)
  J <- max(industry_idx)

  simulate_one <- function(b) {
    eta_idx <- sample.int(length(supp_eta), J, replace = TRUE, prob = g_eta)
    eta_by_firm <- supp_eta[eta_idx][industry_idx]
    xi <- supp_xi[sample.int(length(supp_xi), N, replace = TRUE, prob = g_xi)]
    v_latent <- if (identical(characteristic, "race")) eta_by_firm * xi else eta_by_firm + xi
    v_hat_bs <- v_latent + as.numeric(fit$s_v) * stats::rnorm(N)
    list(v_hat = v_hat_bs, start = par_hat + stats::rnorm(n_free))
  }
  estimate_one <- function(sim, b) {
    sp <- .gp_inf_split_twolevel_free(sim$start, Q_xi, Q_eta, characteristic)
    .gp_2l_fit_node(
      c_xi = prior$diagnostics$c_xi,
      c_eta = prior$diagnostics$c_eta,
      r = sim$v_hat,
      s_v = as.numeric(fit$s_v),
      industry_idx = industry_idx,
      Q_xi = Q_xi,
      supp_xi = supp_xi,
      Q_eta = Q_eta,
      supp_eta = supp_eta,
      characteristic = characteristic,
      start = c(sp$alpha_xi_free, sp$alpha_eta_free),
      max_iter = max_iter,
      tol = tol
    )
  }
  statistic_one <- function(est, b, sim) {
    gp_delta_targets(
      c(est$alpha_xi_free, est$alpha_eta_free),
      Q_xi = Q_xi,
      supp_xi = supp_xi,
      Q_eta = Q_eta,
      supp_eta = supp_eta,
      s = input$s,
      mu = .gp_inf_or(fit$mu, if (identical(characteristic, "race")) fit$m_hat[1L] else 0),
      beta = fit$beta,
      characteristic = characteristic,
      supp_pts_theta = supp_pts_theta
    )
  }

  out <- .gp_param_bootstrap_loop(
    simulate = simulate_one,
    estimate = estimate_one,
    statistic = statistic_one,
    n_rep = n_rep,
    seed = seed,
    keep_draws = keep_draws
  )
  out$diagnostics$method <- "two-level-parametric-bootstrap"
  out$diagnostics$characteristic <- characteristic
  out$diagnostics$c_xi <- prior$diagnostics$c_xi
  out$diagnostics$c_eta <- prior$diagnostics$c_eta
  out$diagnostics$n_xi_free <- split_hat$n_xi_free
  out$diagnostics$n_eta_free <- split_hat$n_eta_free
  out
}
