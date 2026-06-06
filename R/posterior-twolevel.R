# =============================================================================
# posterior-twolevel.R  --  native two-level posterior summaries
# -----------------------------------------------------------------------------
# Ports the `group_fx == 1` posterior-summary contract in KRW's
# matlab_support/get_posteriors.m: industry likelihoods are shared within an
# industry, then posterior means and weighted-percentile intervals are computed
# over the latent xi x eta support.  This implementation uses the deterministic
# factorization
#
#   m_i(a)   = sum_b g_xi(b) ell_i(a,b)
#   Z_k(a)   = g_eta(a) prod_{i in k} m_i(a)
#   q_i(b|a) = g_xi(b) ell_i(a,b) / m_i(a)
#
# rather than materializing all xi combinations in an industry.  The validated
# gp_posterior shell stays on the r/v scale because validate_gp_posterior()
# currently accepts only scale == "r"; theta-scale KRW posterior summaries ride
# in metadata$reporting, following the one-level posterior pattern.
# =============================================================================

#' Resolve package abort helper
#' @keywords internal
#' @noRd
.gp_post2_abort <- function(msg, ..., class = "gradepath_error") {
  args <- list(...)
  if (length(args) > 0L) {
    msg <- do.call(sprintf, c(list(msg), args))
  }
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

#' Validate and normalize a discrete support/mass pair
#' @keywords internal
#' @noRd
.gp_post2_mass <- function(support, mass, name) {
  support <- as.numeric(support)
  mass <- as.numeric(mass)
  if (length(support) == 0L || length(support) != length(mass)) {
    .gp_post2_abort(
      "`%s` support and mass must be non-empty and length-matched.",
      name,
      class = "gradepath_validation_error"
    )
  }
  if (any(!is.finite(support)) || any(!is.finite(mass)) ||
      any(mass < 0) || !is.finite(sum(mass)) || sum(mass) <= 0) {
    .gp_post2_abort(
      "`%s` support/mass must be finite with non-negative positive total mass.",
      name,
      class = "gradepath_validation_error"
    )
  }
  if (any(diff(support) <= 0)) {
    .gp_post2_abort(
      "`%s` support must be strictly increasing.",
      name,
      class = "gradepath_validation_error"
    )
  }
  list(support = support, mass = mass / sum(mass))
}

#' Weighted percentile with KRW wprctile.m default type 5 convention
#'
#' Matlab wprctile.m rescales weights to sum n, then uses
#' pk = (cumsum(sortedW) - sortedW / 2) / n and linear interpolation over
#' c(0, pk, 1).  The rescaling cancels algebraically, so this implementation
#' uses normalized cumulative weights directly.
#' @keywords internal
#' @noRd
.gp_weighted_percentile <- function(x, p, w = NULL) {
  x <- as.numeric(x)
  p <- as.numeric(p)
  if (is.null(w)) w <- rep(1, length(x))
  w <- as.numeric(w)
  if (length(x) == 0L || length(w) != length(x)) {
    .gp_post2_abort("`x` and `w` must be non-empty and length-matched.",
                    class = "gradepath_validation_error")
  }
  if (length(p) == 0L || any(!is.finite(p)) || any(p < 0 | p > 100)) {
    .gp_post2_abort("`p` must contain values in [0, 100].",
                    class = "gradepath_validation_error")
  }
  if (any(!is.finite(x)) || any(!is.finite(w)) || any(w < 0) ||
      sum(w) <= 0 || !is.finite(sum(w))) {
    .gp_post2_abort("`x`/`w` must be finite and `w` must have positive non-negative mass.",
                    class = "gradepath_validation_error")
  }

  ord <- order(x)
  xs <- x[ord]
  ws <- w[ord]
  pk <- (cumsum(ws) - 0.5 * ws) / sum(ws)
  stats::approx(
    x = c(0, pk, 1),
    y = c(xs[1L], xs, xs[length(xs)]),
    xout = p / 100,
    ties = "ordered",
    rule = 2
  )$y
}

#' Pull the xi/eta support masses from the two-level prior object
#' @keywords internal
#' @noRd
.gp_post2_prior_parts <- function(prior) {
  if (!is.list(prior)) {
    .gp_post2_abort("`prior` must be a list with two-level support fields.",
                    class = "gradepath_validation_error")
  }
  if (is.null(prior$support) || is.null(prior$density)) {
    .gp_post2_abort("`prior$support` and `prior$density` are required.",
                    class = "gradepath_validation_error")
  }
  if (is.null(prior$diagnostics$support_eta) ||
      is.null(prior$diagnostics$g_eta)) {
    .gp_post2_abort(
      "`prior$diagnostics$support_eta` and `prior$diagnostics$g_eta` are required.",
      class = "gradepath_validation_error"
    )
  }
  xi <- .gp_post2_mass(prior$support, prior$density, "prior$xi")
  eta <- .gp_post2_mass(prior$diagnostics$support_eta,
                        prior$diagnostics$g_eta,
                        "prior$eta")
  list(supp_xi = xi$support, g_xi = xi$mass,
       supp_eta = eta$support, g_eta = eta$mass)
}

#' Resolve characteristic/beta/mu and carrier vectors
#' @keywords internal
#' @noRd
.gp_post2_carriers <- function(input, prior, fit) {
  if (!is.list(input)) {
    .gp_post2_abort("`input` must be a list carrying theta_hat, s, and industry.",
                    class = "gradepath_validation_error")
  }
  theta_hat <- as.numeric(input$theta_hat)
  s <- as.numeric(input$s)
  if (length(theta_hat) == 0L || length(theta_hat) != length(s) ||
      any(!is.finite(theta_hat)) || any(!is.finite(s)) || any(s <= 0)) {
    .gp_post2_abort(
      "`input$theta_hat` and `input$s` must be finite, positive-SE, length-matched vectors.",
      class = "gradepath_validation_error"
    )
  }
  N <- length(theta_hat)

  characteristic <- fit$characteristic %gp_or%
    prior$metadata$characteristic %gp_or% input$characteristic
  if (is.null(characteristic)) {
    .gp_post2_abort("`fit$characteristic` (race/gender) is required.",
                    class = "gradepath_validation_error")
  }
  characteristic <- match.arg(as.character(characteristic), c("race", "gender"))

  beta <- fit$beta %gp_or% prior$metadata$beta
  beta <- as.numeric(beta)
  if (length(beta) != 1L || !is.finite(beta)) {
    .gp_post2_abort("`fit$beta` must be a finite scalar.",
                    class = "gradepath_validation_error")
  }
  mu <- fit$mu %gp_or% prior$metadata$mu %gp_or% 0
  mu <- as.numeric(mu)
  if (length(mu) != 1L || !is.finite(mu)) {
    .gp_post2_abort("`fit$mu`/`prior$metadata$mu` must be a finite scalar.",
                    class = "gradepath_validation_error")
  }

  v_hat <- fit$v_hat
  if (is.null(v_hat)) {
    v_hat <- if (identical(characteristic, "race")) {
      theta_hat / (s^beta)
    } else {
      (theta_hat - mu) / (s^beta)
    }
  }
  v_hat <- as.numeric(v_hat)
  s_v <- fit$s_v
  if (is.null(s_v)) s_v <- s^(1 - beta)
  s_v <- as.numeric(s_v)
  if (length(v_hat) != N || length(s_v) != N ||
      any(!is.finite(v_hat)) || any(!is.finite(s_v)) || any(s_v <= 0)) {
    .gp_post2_abort("`fit$v_hat`/`fit$s_v` must be finite positive-SE length-N vectors.",
                    class = "gradepath_validation_error")
  }

  industry <- fit$industry %gp_or% input$industry
  if (is.null(industry) || length(industry) != N || any(is.na(industry))) {
    .gp_post2_abort("`industry` must be a non-missing length-N vector.",
                    class = "gradepath_validation_error")
  }
  if (is.numeric(industry) || is.integer(industry)) {
    raw_levels <- sort(unique(as.numeric(industry)), method = "radix")
    industry_idx <- match(as.numeric(industry), raw_levels)
    levels <- as.character(raw_levels)
  } else {
    levels <- sort(unique(as.character(industry)), method = "radix")
    industry_idx <- match(as.character(industry), levels)
  }
  if (length(levels) < 1L) {
    .gp_post2_abort("`industry` must contain at least one group.",
                    class = "gradepath_validation_error")
  }

  id <- input$unit_id %gp_or% input$id %gp_or% input$ids %gp_or% seq_len(N)
  id <- as.character(id)
  if (length(id) != N || any(is.na(id)) || any(duplicated(id))) {
    id <- as.character(seq_len(N))
  }
  label <- input$label %gp_or% input$labels %gp_or% id
  label <- as.character(label)
  if (length(label) != N || any(is.na(label))) label <- id

  list(
    theta_hat = theta_hat,
    s = s,
    v_hat = v_hat,
    s_v = s_v,
    industry = industry_idx,
    industry_levels = levels,
    characteristic = characteristic,
    beta = beta,
    mu = mu,
    id = id,
    label = label
  )
}

#' Build log kernels log(g_xi * likelihood) for every firm and eta node
#' @keywords internal
#' @noRd
.gp_post2_log_kernels <- function(parts, carriers) {
  supp_xi <- parts$supp_xi
  supp_eta <- parts$supp_eta
  log_g_xi <- log(parts$g_xi)
  N <- length(carriers$v_hat)
  E <- length(supp_eta)
  M <- length(supp_xi)

  log_kernel <- vector("list", N)
  log_m <- matrix(NA_real_, nrow = N, ncol = E)
  for (i in seq_len(N)) {
    means <- if (identical(carriers$characteristic, "race")) {
      outer(supp_eta, supp_xi, `*`)
    } else {
      outer(supp_eta, supp_xi, `+`)
    }
    ll <- matrix(stats::dnorm(
      carriers$v_hat[i],
      mean = means,
      sd = carriers$s_v[i],
      log = TRUE
    ), nrow = E, ncol = M)
    K <- sweep(ll, 2L, log_g_xi, `+`)
    log_kernel[[i]] <- K
    log_m[i, ] <- .gp_row_logsumexp(K)
  }
  list(log_kernel = log_kernel, log_m = log_m, E = E, M = M)
}

#' Native two-level per-industry posterior summaries
#'
#' @description
#' Computes two-level posterior summaries for the `group_fx == 1` KRW model using
#' deterministic xi/eta quadrature.  The validated `gp_posterior` object carries
#' r/v-scale summaries; `metadata$reporting` carries the theta-scale KRW
#' `posteriors_groupfx1`-style table (`posterior_mean`,
#' `posterior_second_moment`, `lower`, `upper`).
#'
#' @param input A two-level input list with `theta_hat`, original `s`, and
#'   `industry` membership (e.g. `gp_krw_gmm_input()`).
#' @param prior A two-level `gp_prior` from `gp_deconvolve_groups()`, with
#'   `support`/`density` for xi and `diagnostics$support_eta`/`g_eta` for eta.
#' @param fit A two-level fit from `gp_two_level_gmm()`.
#' @param control Optional control list; `control$interval_level` is used when
#'   `interval_level` is absent.
#' @param interval_level Central interval level.  Defaults to 0.90; use 0.95 for
#'   archive 2.5/97.5 percentile parity reads.
#'
#' @return A validated `gp_posterior`.
#' @keywords internal
#' @noRd
gp_posterior_twolevel <- function(input, prior, fit, control = NULL,
                                  interval_level = NULL) {
  parts <- .gp_post2_prior_parts(prior)
  carriers <- .gp_post2_carriers(input, prior, fit)
  N <- length(carriers$v_hat)

  level <- interval_level %gp_or%
    (if (!is.null(control) && !is.null(control$interval_level)) {
      control$interval_level
    } else {
      0.90
    })
  level <- as.numeric(level)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    .gp_post2_abort("`interval_level` must be a scalar in (0, 1).",
                    class = "gradepath_validation_error")
  }
  pct <- c((1 - level) / 2, (1 + level) / 2) * 100

  kernels <- .gp_post2_log_kernels(parts, carriers)
  log_g_eta <- log(parts$g_eta)
  K <- max(carriers$industry)
  E <- kernels$E

  log_z <- matrix(NA_real_, nrow = K, ncol = E)
  log_den <- numeric(K)
  eta_posterior <- matrix(NA_real_, nrow = K, ncol = E)
  for (k in seq_len(K)) {
    idx <- which(carriers$industry == k)
    log_z[k, ] <- log_g_eta + colSums(kernels$log_m[idx, , drop = FALSE])
    log_den[k] <- .gp_logsumexp(log_z[k, ])
    if (!is.finite(log_den[k])) {
      .gp_post2_abort("Degenerate industry posterior mass for industry %d.", k,
                      class = "gradepath_validation_error")
    }
    eta_posterior[k, ] <- exp(log_z[k, ] - log_den[k])
  }

  pm_r <- psd_r <- e2_r <- lo_r <- up_r <- numeric(N)
  pm_t <- psd_t <- e2_t <- lo_t <- up_t <- numeric(N)
  e_xi <- e_eta <- e_sbar_eta <- lo_sbar_eta <- up_sbar_eta <- numeric(N)
  s_beta <- carriers$s^carriers$beta
  sbar_by_industry <- numeric(K)
  for (k in seq_len(K)) {
    idx <- which(carriers$industry == k)
    sbar_by_industry[k] <- mean(s_beta[idx])
  }
  sbar <- sbar_by_industry[carriers$industry]

  for (i in seq_len(N)) {
    k <- carriers$industry[i]
    log_w <- sweep(
      kernels$log_kernel[[i]],
      1L,
      log_z[k, ] - log_den[k] - kernels$log_m[i, ],
      `+`
    )
    w <- exp(as.vector(log_w))
    total <- sum(w)
    if (!is.finite(total) || total <= 0) {
      .gp_post2_abort("Degenerate posterior mass for unit %d.", i,
                      class = "gradepath_validation_error")
    }
    w <- w / total

    eta_grid <- rep(parts$supp_eta, times = length(parts$supp_xi))
    xi_grid <- rep(parts$supp_xi, each = length(parts$supp_eta))
    # as.vector(matrix) is column-major: eta varies fastest within each xi column.
    if (identical(carriers$characteristic, "race")) {
      r_grid <- eta_grid * xi_grid
      theta_grid <- s_beta[i] * r_grid
    } else {
      r_grid <- eta_grid + xi_grid
      theta_grid <- carriers$mu + s_beta[i] * r_grid
    }
    sbar_eta_grid <- sbar[i] * eta_grid

    pm_r[i] <- sum(w * r_grid)
    e2_r[i] <- sum(w * r_grid^2)
    psd_r[i] <- sqrt(max(e2_r[i] - pm_r[i]^2, 0))
    q_r <- .gp_weighted_percentile(r_grid, pct, w)
    lo_r[i] <- min(q_r[1L], pm_r[i])
    up_r[i] <- max(q_r[2L], pm_r[i])

    pm_t[i] <- sum(w * theta_grid)
    e2_t[i] <- sum(w * theta_grid^2)
    psd_t[i] <- sqrt(max(e2_t[i] - pm_t[i]^2, 0))
    q_t <- .gp_weighted_percentile(theta_grid, pct, w)
    lo_t[i] <- q_t[1L]
    up_t[i] <- q_t[2L]

    e_xi[i] <- sum(w * xi_grid)
    e_eta[i] <- sum(w * eta_grid)
    e_sbar_eta[i] <- sum(w * sbar_eta_grid)
    q_sbar <- .gp_weighted_percentile(sbar_eta_grid, pct, w)
    lo_sbar_eta[i] <- q_sbar[1L]
    up_sbar_eta[i] <- q_sbar[2L]
  }

  posteriors <- data.frame(
    posterior_mean = pm_t,
    posterior_second_moment = e2_t,
    lower = lo_t,
    upper = up_t
  )
  reporting <- list(
    posterior_mean = pm_t,
    posterior_sd = psd_t,
    posterior_second_moment = e2_t,
    lower = lo_t,
    upper = up_t,
    posteriors = posteriors,
    scale = "theta",
    level = level,
    percentile_convention = "collapsed_marginal_wprctile_type5"
  )
  two_level <- list(
    characteristic = carriers$characteristic,
    beta = carriers$beta,
    mu = carriers$mu,
    industry = carriers$industry,
    industry_levels = carriers$industry_levels,
    support_xi = parts$supp_xi,
    g_xi = parts$g_xi,
    support_eta = parts$supp_eta,
    g_eta = parts$g_eta,
    log_m = kernels$log_m,
    log_z = log_z,
    log_denominator = log_den,
    eta_posterior = eta_posterior,
    log_kernel = kernels$log_kernel,
    sbar = sbar,
    s_beta = s_beta,
    original_s = carriers$s,
    percentile_convention = "collapsed_marginal_wprctile_type5",
    posterior_components = data.frame(
      E_r = pm_r,
      E_r2 = e2_r,
      E_theta = pm_t,
      E_theta2 = e2_t,
      E_xi = e_xi,
      E_eta = e_eta,
      E_sbar_eta = e_sbar_eta,
      lower_sbar_eta = lo_sbar_eta,
      upper_sbar_eta = up_sbar_eta
    )
  )

  obj <- new_gp_posterior(
    estimate = carriers$v_hat,
    se = carriers$s_v,
    id = carriers$id,
    label = carriers$label,
    posterior_mean = pm_r,
    posterior_sd = psd_r,
    lower = lo_r,
    upper = up_r,
    scale = "r",
    metadata = list(
      level = level,
      interval_level = level,
      reporting = reporting,
      has_reporting = TRUE,
      r_second_moment = e2_r,
      two_level = two_level
    )
  )
  validate_gp_posterior(obj)
}
