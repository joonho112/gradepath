# =============================================================================
# two-level-pushforward.R  --  native TWO-level theta pushforward
# -----------------------------------------------------------------------------
# Ports KRW's matlab_support/get_g_theta.m for the two-support group_fx == 1 path.
# The one-level `.gp_eb_pushforward_theta` remains a cross-check only; it cannot
# integrate over the second eta/psi support.
# =============================================================================

#' Resolve package abort helper
#' @keywords internal
#' @noRd
.gp_pushfwd_abort <- function(msg, class = "gradepath_error") {
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

#' Validate and normalize a discrete support/density pair
#' @keywords internal
#' @noRd
.gp_pushfwd_density <- function(support, density, name) {
  support <- as.numeric(support)
  density <- as.numeric(density)
  if (length(support) == 0L || length(support) != length(density)) {
    .gp_pushfwd_abort(
      sprintf("`%s` support and density must be non-empty and length-matched.", name),
      class = "gradepath_validation_error"
    )
  }
  if (any(!is.finite(support)) || any(!is.finite(density)) ||
      any(density < 0) || !is.finite(sum(density)) || sum(density) <= 0) {
    .gp_pushfwd_abort(
      sprintf("`%s` support/density must be finite with non-negative positive-mass density.", name),
      class = "gradepath_validation_error"
    )
  }
  if (any(diff(support) <= 0)) {
    .gp_pushfwd_abort(
      sprintf("`%s` support must be strictly increasing.", name),
      class = "gradepath_validation_error"
    )
  }
  list(support = support, density = density / sum(density))
}

#' Theta values on the full N x xi x eta support
#' @keywords internal
#' @noRd
.gp_theta_values <- function(supp_xi, supp_eta, s, mu, beta, characteristic) {
  vals <- vector("list", length(s))
  xi_eta <- if (identical(characteristic, "race")) {
    as.vector(outer(supp_xi, supp_eta, `*`))
  } else {
    as.vector(outer(supp_xi, supp_eta, `+`))
  }
  for (i in seq_along(s)) {
    vals[[i]] <- if (identical(characteristic, "race")) {
      (s[i]^beta) * xi_eta
    } else {
      mu + (s[i]^beta) * xi_eta
    }
  }
  vals
}

#' Deposit one support value's mass on the nearest theta-grid node(s)
#'
#' Mirrors get_g_theta.m:39-41: `diffs=abs(val-supp_theta); mindiff=min(diffs);
#' G(:,t)=G(:,t)+mass*(diffs==mindiff)`. On an exact tie, Matlab adds the full
#' mass to every tied node; the final `g_theta/sum(g_theta)` normalization then
#' makes this behave as the reference implementation.
#'
#' @keywords internal
#' @noRd
.gp_pushfwd_deposit <- function(G_col, val, mass, supp_theta, step) {
  M <- length(supp_theta)
  if (M == 1L || step == 0) {
    G_col[1L] <- G_col[1L] + mass
    return(G_col)
  }
  pos <- (val - supp_theta[1L]) / step + 1
  lo <- max(1L, min(M, floor(pos)))
  hi <- max(1L, min(M, ceiling(pos)))
  if (lo == hi) {
    G_col[lo] <- G_col[lo] + mass
    return(G_col)
  }
  d_lo <- abs(val - supp_theta[lo])
  d_hi <- abs(val - supp_theta[hi])
  mindiff <- min(d_lo, d_hi)
  if (d_lo == mindiff) G_col[lo] <- G_col[lo] + mass
  if (d_hi == mindiff) G_col[hi] <- G_col[hi] + mass
  G_col
}

#' Native two-level theta pushforward
#'
#' @description
#' Faithful port of `matlab_support/get_g_theta.m`. Given `G_xi`, `G_eta`, standard
#' errors `s`, location `mu`, precision exponent `beta`, and `characteristic`,
#' builds the full transformed theta support:
#'
#' - race: `theta = s^beta * xi * eta`
#' - gender: `theta = mu + s^beta * (xi + eta)`
#'
#' It then creates the fixed 250-point theta grid, deposits each `xi x eta` mass at
#' every nearest grid node (`diffs == min(diffs)`, not `which.min`), averages the
#' conditional distributions over firms, and mass-normalizes.
#'
#' @param supp_xi,g_xi Within-unit support and probability masses.
#' @param supp_eta,g_eta Between-industry support and probability masses.
#' @param s Length-N original standard errors.
#' @param mu Additive location for gender; ignored for race except for signature
#'   symmetry with Matlab.
#' @param beta Precision exponent.
#' @param characteristic `"race"` or `"gender"`.
#' @param supp_pts_theta Fixed theta-grid size; defaults to 250 to match Matlab.
#' @return A list with `support`, `g` (normalized theta probability mass), and
#'   `density` (`g` divided by theta-grid width, matching archive CSV scale),
#'   plus lightweight diagnostics.
#' @keywords internal
#' @noRd
gp_pushforward_theta <- function(supp_xi, g_xi, supp_eta, g_eta,
                                 s, mu = 0, beta, characteristic,
                                 supp_pts_theta = 250L) {
  characteristic <- match.arg(characteristic, c("race", "gender"))
  xi <- .gp_pushfwd_density(supp_xi, g_xi, "xi")
  eta <- .gp_pushfwd_density(supp_eta, g_eta, "eta")
  s <- as.numeric(s)
  if (length(s) == 0L || any(!is.finite(s)) || any(s <= 0)) {
    .gp_pushfwd_abort("`s` must be a non-empty positive finite numeric vector.",
                      class = "gradepath_validation_error")
  }
  if (length(beta) != 1L || !is.finite(beta)) {
    .gp_pushfwd_abort("`beta` must be a finite scalar.",
                      class = "gradepath_validation_error")
  }
  if (length(mu) != 1L || !is.finite(mu)) {
    .gp_pushfwd_abort("`mu` must be a finite scalar.",
                      class = "gradepath_validation_error")
  }
  supp_pts_theta <- as.integer(supp_pts_theta)
  if (supp_pts_theta < 2L) {
    .gp_pushfwd_abort("`supp_pts_theta` must be at least 2.",
                      class = "gradepath_validation_error")
  }

  values_by_s <- .gp_theta_values(
    xi$support, eta$support, s, mu = mu, beta = beta,
    characteristic = characteristic
  )
  all_values <- unlist(values_by_s, use.names = FALSE)
  theta_min <- min(all_values)
  theta_max <- max(all_values)
  if (!is.finite(theta_min) || !is.finite(theta_max) || theta_max <= theta_min) {
    .gp_pushfwd_abort(
      "`theta` transformed support must have positive width.",
      class = "gradepath_validation_error"
    )
  }
  supp_theta <- seq(theta_min, theta_max, length.out = supp_pts_theta)
  step <- if (length(supp_theta) > 1L) supp_theta[2L] - supp_theta[1L] else 0

  weights <- as.vector(outer(xi$density, eta$density, `*`))
  G <- matrix(0, nrow = supp_pts_theta, ncol = length(s))
  for (t in seq_along(s)) {
    vals <- values_by_s[[t]]
    col <- numeric(supp_pts_theta)
    for (k in seq_along(vals)) {
      col <- .gp_pushfwd_deposit(col, vals[k], weights[k], supp_theta, step)
    }
    G[, t] <- col
  }

  g_theta <- rowMeans(G)
  g_theta <- g_theta / sum(g_theta)
  density_theta <- g_theta / step
  list(
    support = supp_theta,
    g = g_theta,
    density = density_theta,
    diagnostics = list(
      method = "native-two-level-get_g_theta",
      characteristic = characteristic,
      supp_pts_theta = supp_pts_theta,
      grid_width = step,
      n_carriers = length(s),
      n_xi = length(xi$support),
      n_eta = length(eta$support)
    )
  )
}
