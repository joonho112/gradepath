# =============================================================================
# two-level.R  --  same-industry override + five Pi matrices
# -----------------------------------------------------------------------------
# Consumes the deterministic posterior blocks produced by
# gp_posterior_twolevel() and assembles the KRW group_fx == 1 pairwise outputs:
# Pi_theta, Pi_sq_theta, Pi_xi, Pi_psi, and Pi_sbar_psi (archive Pi_bar).
#
# The load-bearing contract is get_posteriors.m:188-199:
#
#   L_ij = L_i * L_j              for different industries
#   L_ij = L_i                    for same industries, counted once
#
# In deterministic quadrature form, same-industry pairs integrate over one
# shared eta posterior; cross-industry pairs integrate over the product of two
# independent eta posteriors.  The raw Matlab-style matrices are preserved, and
# the theta probability matrix is also wrapped as a cleaned gp_pairwise object.
# =============================================================================

#' Resolve package abort helper
#' @keywords internal
#' @noRd
.gp_tl_abort <- function(msg, ..., class = "gradepath_error") {
  args <- list(...)
  if (length(args) > 0L) msg <- do.call(sprintf, c(list(msg), args))
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

#' Normalize a finite non-negative mass vector
#' @keywords internal
#' @noRd
.gp_tl_normalize <- function(w, name = "w") {
  w <- as.numeric(w)
  if (length(w) == 0L || any(!is.finite(w)) || any(w < 0) ||
      !is.finite(sum(w)) || sum(w) <= 0) {
    .gp_tl_abort("`%s` must be finite non-negative weights with positive mass.",
                 name, class = "gradepath_validation_error")
  }
  w / sum(w)
}

#' Weighted strict ordering probability and positive squared gap
#' @keywords internal
#' @noRd
.gp_tl_prepare_weighted <- function(y, wy) {
  y <- as.numeric(y)
  wy <- .gp_tl_normalize(wy, "wy")
  if (length(y) != length(wy) || any(!is.finite(y))) {
    .gp_tl_abort("Prepared values and weights must be finite and length-matched.",
                 class = "gradepath_validation_error")
  }
  ord <- order(y)
  ys <- y[ord]
  wys <- wy[ord]
  list(
    y = ys,
    cw0 = c(0, cumsum(wys)),
    cw1 = c(0, cumsum(wys * ys)),
    cw2 = c(0, cumsum(wys * ys^2))
  )
}

#' Weighted strict ordering probability using a prepared rhs distribution
#' @keywords internal
#' @noRd
.gp_tl_weighted_pair_prepared <- function(x, wx, prepared_y, squared = FALSE) {
  x <- as.numeric(x)
  wx <- .gp_tl_normalize(wx, "wx")
  if (length(x) != length(wx) || any(!is.finite(x)) ||
      !is.list(prepared_y) || is.null(prepared_y$y) ||
      is.null(prepared_y$cw0)) {
    .gp_tl_abort("Values and prepared weights must be finite and length-matched.",
                 class = "gradepath_validation_error")
  }
  idx <- findInterval(x, prepared_y$y, left.open = TRUE)
  if (!squared) {
    return(sum(wx * prepared_y$cw0[idx + 1L]))
  }

  inner <- x^2 * prepared_y$cw0[idx + 1L] -
    2 * x * prepared_y$cw1[idx + 1L] +
    prepared_y$cw2[idx + 1L]
  sum(wx * inner)
}

#' Weighted strict ordering probability and positive squared gap
#' @keywords internal
#' @noRd
.gp_tl_weighted_pair <- function(x, wx, y, wy, squared = FALSE) {
  .gp_tl_weighted_pair_prepared(
    x,
    wx,
    .gp_tl_prepare_weighted(y, wy),
    squared = squared
  )
}

#' Extract and validate two-level posterior metadata
#' @keywords internal
#' @noRd
.gp_tl_extract <- function(posterior) {
  if (!inherits(posterior, "gp_posterior")) {
    .gp_tl_abort("`posterior` must be a gp_posterior from gp_posterior_twolevel().",
                 class = "gradepath_validation_error")
  }
  posterior <- validate_gp_posterior(posterior)
  tl <- posterior$metadata$two_level
  req <- c("characteristic", "beta", "mu", "industry", "support_xi", "g_xi",
           "support_eta", "g_eta", "log_m", "log_z", "log_denominator",
           "eta_posterior", "log_kernel", "sbar")
  if (!is.list(tl) || any(!req %in% names(tl))) {
    .gp_tl_abort(
      "`posterior$metadata$two_level` is missing two-level likelihood blocks.",
      class = "gradepath_validation_error"
    )
  }

  N <- length(posterior$id)
  industry <- as.integer(tl$industry)
  if (length(industry) != N || anyNA(industry) || any(industry < 1L)) {
    .gp_tl_abort("`two_level$industry` must be a positive integer length-N vector.",
                 class = "gradepath_validation_error")
  }
  K <- max(industry)
  if (!identical(sort(unique(industry)), seq_len(K))) {
    .gp_tl_abort("`two_level$industry` must be dense and cover 1:K.",
                 class = "gradepath_validation_error")
  }
  supp_xi <- as.numeric(tl$support_xi)
  supp_eta <- as.numeric(tl$support_eta)
  M <- length(supp_xi)
  E <- length(supp_eta)
  if (!is.matrix(tl$log_m) || !identical(dim(tl$log_m), c(N, E)) ||
      !is.matrix(tl$log_z) || !identical(dim(tl$log_z), c(K, E)) ||
      !is.matrix(tl$eta_posterior) || !identical(dim(tl$eta_posterior), c(K, E)) ||
      length(tl$log_kernel) != N) {
    .gp_tl_abort("Two-level likelihood block dimensions are inconsistent.",
                 class = "gradepath_validation_error")
  }
  if (any(!is.finite(tl$eta_posterior)) || any(tl$eta_posterior < 0) ||
      any(abs(rowSums(tl$eta_posterior) - 1) > 1e-8)) {
    .gp_tl_abort("`two_level$eta_posterior` must be row-stochastic.",
                 class = "gradepath_validation_error")
  }
  for (i in seq_len(N)) {
    if (!is.matrix(tl$log_kernel[[i]]) ||
        !identical(dim(tl$log_kernel[[i]]), c(E, M))) {
      .gp_tl_abort("`two_level$log_kernel[[%d]]` must be E x M.", i,
                   class = "gradepath_validation_error")
    }
  }
  if (length(tl$log_denominator) != K || length(tl$sbar) != N ||
      any(!is.finite(tl$sbar))) {
    .gp_tl_abort("Two-level denominator/sbar lengths are inconsistent.",
                 class = "gradepath_validation_error")
  }
  for (k in seq_len(K)) {
    idx <- industry == k
    if (any(abs(tl$sbar[idx] - tl$sbar[idx][1L]) > 1e-10)) {
      .gp_tl_abort("`two_level$sbar` must be constant within each industry.",
                   class = "gradepath_validation_error")
    }
  }
  tl$ids <- posterior$id
  tl$N <- N; tl$K <- K; tl$M <- M; tl$E <- E
  tl$beta <- as.numeric(tl$beta)
  tl$mu <- as.numeric(tl$mu)
  if (length(tl$beta) != 1L || !is.finite(tl$beta) ||
      length(tl$mu) != 1L || !is.finite(tl$mu)) {
    .gp_tl_abort("`two_level$beta` and `two_level$mu` must be finite scalars.",
                 class = "gradepath_validation_error")
  }
  if (!identical(tl$characteristic, "race") &&
      !identical(tl$characteristic, "gender")) {
    .gp_tl_abort("`two_level$characteristic` must be 'race' or 'gender'.",
                 class = "gradepath_validation_error")
  }
  if (is.null(tl$s_beta)) {
    if (abs(1 - tl$beta) < 1e-8) {
      .gp_tl_abort(
        "`two_level$s_beta` is required when beta is too close to 1 to reconstruct from s_v.",
        class = "gradepath_validation_error"
      )
    }
    tl$s_beta <- as.numeric(posterior$se)^(tl$beta / (1 - tl$beta))
  } else {
    tl$s_beta <- as.numeric(tl$s_beta)
  }
  if (length(tl$s_beta) != N || any(!is.finite(tl$s_beta)) ||
      any(tl$s_beta <= 0)) {
    .gp_tl_abort("`two_level$s_beta` must be a finite positive length-N vector.",
                 class = "gradepath_validation_error")
  }
  tl
}

#' Conditional xi weights q_i(xi | eta, data)
#' @keywords internal
#' @noRd
.gp_tl_qxi <- function(tl, i, a) {
  .gp_tl_normalize(exp(tl$log_kernel[[i]][a, ] - tl$log_m[i, a]), "q_xi")
}

#' Theta support for one unit at one eta node
#' @keywords internal
#' @noRd
.gp_tl_theta_values <- function(tl, i, eta, xi) {
  s_beta <- tl$s_beta
  if (is.null(s_beta) || length(s_beta) < i || !is.finite(s_beta[i])) {
    .gp_tl_abort(
      "`posterior$metadata$two_level$s_beta` is required for theta pairwise matrices.",
      class = "gradepath_validation_error"
    )
  }
  if (identical(tl$characteristic, "race")) {
    s_beta[i] * eta * xi
  } else {
    tl$mu + s_beta[i] * (eta + xi)
  }
}

#' Build marginal distributions for a unit over all eta x xi nodes
#' @keywords internal
#' @noRd
.gp_tl_unit_marginals <- function(tl, i) {
  k <- tl$industry[i]
  vals_theta <- vals_xi <- weights <- numeric(0)
  for (a in seq_len(tl$E)) {
    q <- .gp_tl_qxi(tl, i, a)
    wa <- tl$eta_posterior[k, a] * q
    vals_theta <- c(vals_theta, .gp_tl_theta_values(tl, i, tl$support_eta[a], tl$support_xi))
    vals_xi <- c(vals_xi, tl$support_xi)
    weights <- c(weights, wa)
  }
  list(
    theta = list(x = vals_theta, w = .gp_tl_normalize(weights, "theta weights")),
    xi = list(x = vals_xi, w = .gp_tl_normalize(weights, "xi weights")),
    psi = list(x = tl$support_eta,
               w = .gp_tl_normalize(tl$eta_posterior[k, ], "psi weights")),
    sbar_psi = list(x = tl$sbar[i] * tl$support_eta,
                    w = .gp_tl_normalize(tl$eta_posterior[k, ], "sbar weights"))
  )
}

#' Same-industry pair integral: one shared eta posterior
#' @keywords internal
#' @noRd
.gp_tl_pair_same <- function(tl, i, j) {
  k <- tl$industry[i]
  out <- c(Pi_theta = 0, Pi_sq_theta = 0, Pi_xi = 0, Pi_psi = 0,
           Pi_sbar_psi = 0)
  for (a in seq_len(tl$E)) {
    wi <- .gp_tl_qxi(tl, i, a)
    wj <- .gp_tl_qxi(tl, j, a)
    eta <- tl$support_eta[a]
    theta_i <- .gp_tl_theta_values(tl, i, eta, tl$support_xi)
    theta_j <- .gp_tl_theta_values(tl, j, eta, tl$support_xi)
    ea <- tl$eta_posterior[k, a]
    out["Pi_theta"] <- out["Pi_theta"] +
      ea * .gp_tl_weighted_pair(theta_i, wi, theta_j, wj)
    out["Pi_sq_theta"] <- out["Pi_sq_theta"] +
      ea * .gp_tl_weighted_pair(theta_i, wi, theta_j, wj, squared = TRUE)
    out["Pi_xi"] <- out["Pi_xi"] +
      ea * .gp_tl_weighted_pair(tl$support_xi, wi, tl$support_xi, wj)
    # psi and sbar_psi share the same eta inside an industry, so strict
    # same-industry orderings are false by construction.
  }
  out
}

#' Cross-industry pair integral: product of independent eta posteriors
#' @keywords internal
#' @noRd
.gp_tl_prepare_marginal <- function(d) {
  list(
    theta = .gp_tl_prepare_weighted(d$theta$x, d$theta$w),
    xi = .gp_tl_prepare_weighted(d$xi$x, d$xi$w),
    psi = .gp_tl_prepare_weighted(d$psi$x, d$psi$w),
    sbar_psi = .gp_tl_prepare_weighted(d$sbar_psi$x, d$sbar_psi$w)
  )
}

#' Cross-industry pair integral: product of independent eta posteriors
#' @keywords internal
#' @noRd
.gp_tl_pair_cross <- function(di, dj, prepared_j = NULL) {
  if (is.null(prepared_j)) prepared_j <- .gp_tl_prepare_marginal(dj)
  c(
    Pi_theta = .gp_tl_weighted_pair_prepared(di$theta$x, di$theta$w,
                                             prepared_j$theta),
    Pi_sq_theta = .gp_tl_weighted_pair_prepared(di$theta$x, di$theta$w,
                                                prepared_j$theta,
                                                squared = TRUE),
    Pi_xi = .gp_tl_weighted_pair_prepared(di$xi$x, di$xi$w, prepared_j$xi),
    Pi_psi = .gp_tl_weighted_pair_prepared(di$psi$x, di$psi$w,
                                           prepared_j$psi),
    Pi_sbar_psi = .gp_tl_weighted_pair_prepared(di$sbar_psi$x,
                                                di$sbar_psi$w,
                                                prepared_j$sbar_psi)
  )
}

#' Clean a raw strict-upper probability matrix into a gp_pairwise matrix
#' @keywords internal
#' @noRd
.gp_tl_clean_probability <- function(raw) {
  N <- nrow(raw)
  upper <- matrix(0, N, N)
  if (N >= 2L) {
    for (i in seq_len(N - 1L)) {
      for (j in seq.int(i + 1L, N)) {
        upper[i, j] <- raw[i, j]
      }
    }
  }
  .gp_pairwise_cleanup_matrix(upper)
}

#' Wrap a cleaned probability matrix in a gp_pairwise object
#' @keywords internal
#' @noRd
.gp_tl_new_pairwise <- function(ids, matrix, control, n_industries,
                                provenance_step) {
  dimnames(matrix) <- list(ids, ids)
  new_gp_pairwise(
    ids = ids,
    matrix = matrix,
    power = 0L,
    cleanup = list(antisymmetry = TRUE,
                   diagonal = .gp_pairwise_diagonal,
                   zero_floor = .gp_pairwise_zero_floor),
    source = list(stage = "posterior",
                  rule = "groupfx1_archive_matrix",
                  assumption = "grouped_industry_dependence"),
    control = control,
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      step = provenance_step,
      n_units = length(ids),
      n_industries = n_industries,
      rule = "same_industry_likelihood_once",
      cleanup_order = "antisymmetry_diagonal_zero_floor"
    ),
    warnings = character(0)
  )
}

#' Native two-level pairwise Pi matrices
#'
#' @description
#' Builds the five `group_fx == 1` Pi matrices from a two-level posterior.  Raw
#' matrices follow the Matlab strict-inequality convention (diagonal zero).  The
#' returned `$pairwise` slot is the cleaned `gp_pairwise` object for `Pi_theta`
#' (diagonal 0.5, antisymmetry, zero-floor) used by the grading layer.
#'
#' @param posterior A `gp_posterior` returned by `gp_posterior_twolevel()`.
#' @param ids Optional ids. Defaults to `posterior$id`.
#' @param control Optional `gp_control`.
#'
#' @return A list of class `gp_two_level_pi`.
#' @keywords internal
#' @noRd
gp_twolevel_pi <- function(posterior, ids = NULL, control = NULL) {
  tl <- .gp_tl_extract(posterior)
  N <- tl$N
  ids <- as.character(ids %gp_or% tl$ids)
  if (length(ids) != N || anyNA(ids) || any(duplicated(ids))) {
    .gp_tl_abort("`ids` must be a unique character vector of length N.",
                 class = "gradepath_validation_error")
  }
  ctrl <- if (is.null(control)) gp_control() else validate_gp_control(control)

  raw <- list(
    Pi_theta = matrix(0, N, N),
    Pi_sq_theta = matrix(0, N, N),
    Pi_xi = matrix(0, N, N),
    Pi_psi = matrix(0, N, N),
    Pi_sbar_psi = matrix(0, N, N)
  )
  marginals <- lapply(seq_len(N), function(i) .gp_tl_unit_marginals(tl, i))
  prepared_marginals <- lapply(marginals, .gp_tl_prepare_marginal)

  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      if (i == j) next
      vals <- if (tl$industry[i] == tl$industry[j]) {
        .gp_tl_pair_same(tl, i, j)
      } else {
        .gp_tl_pair_cross(marginals[[i]], marginals[[j]],
                          prepared_marginals[[j]])
      }
      for (nm in names(vals)) raw[[nm]][i, j] <- vals[[nm]]
    }
  }
  raw$Pi_bar <- raw$Pi_sbar_psi
  for (nm in names(raw)) dimnames(raw[[nm]]) <- list(ids, ids)

  reps <- vapply(seq_len(tl$K), function(k) which(tl$industry == k)[1L], integer(1))
  bar_ids <- tl$industry_levels %gp_or% as.character(seq_len(tl$K))
  bar_ids <- as.character(bar_ids)
  if (length(bar_ids) != tl$K || anyNA(bar_ids) || any(duplicated(bar_ids))) {
    bar_ids <- paste0("industry_", seq_len(tl$K))
  }
  Pi_sbar_psi_raw <- raw$Pi_sbar_psi[reps, reps, drop = FALSE]
  dimnames(Pi_sbar_psi_raw) <- list(bar_ids, bar_ids)

  pi_theta_clean <- .gp_tl_clean_probability(raw$Pi_theta)
  dimnames(pi_theta_clean) <- list(ids, ids)
  pi_bar_clean <- .gp_tl_clean_probability(Pi_sbar_psi_raw)
  dimnames(pi_bar_clean) <- list(bar_ids, bar_ids)

  pairwise_theta <- .gp_tl_new_pairwise(
    ids = ids,
    matrix = pi_theta_clean,
    control = ctrl,
    n_industries = tl$K,
    provenance_step = "two-level-pairwise-theta"
  )
  pairwise_bar <- .gp_tl_new_pairwise(
    ids = bar_ids,
    matrix = pi_bar_clean,
    control = ctrl,
    n_industries = tl$K,
    provenance_step = "two-level-pairwise-sbar-psi"
  )

  out <- list(
    ids = ids,
    industry = tl$industry,
    industry_levels = bar_ids,
    industry_representatives = reps,
    raw = list(
      Pi_theta = raw$Pi_theta,
      Pi_sbar_psi = raw$Pi_sbar_psi,
      Pi_bar = raw$Pi_bar,
      Pi_sq_theta = raw$Pi_sq_theta,
      Pi_xi = raw$Pi_xi,
      Pi_psi = raw$Pi_psi
    ),
    Pi_theta = pi_theta_clean,
    Pi_sbar_psi = raw$Pi_sbar_psi,
    Pi_bar = raw$Pi_bar,
    Pi_sbar_psi_industry = pi_bar_clean,
    Pi_bar_industry = pi_bar_clean,
    Pi_sq_theta = raw$Pi_sq_theta,
    Pi_xi = raw$Pi_xi,
    Pi_psi = raw$Pi_psi,
    pairwise = pairwise_theta,
    pairwise_theta = pairwise_theta,
    pairwise_bar = pairwise_bar,
    source = list(
      group_fx = 1L,
      rule = "same_industry_likelihood_once",
      same_industry = "shared_eta_integral",
      cross_industry = "independent_eta_product"
    ),
    metadata = list(
      industry = tl$industry,
      industry_levels = bar_ids,
      raw_Pi_sbar_psi_industry = Pi_sbar_psi_raw,
      characteristic = tl$characteristic,
      percentile_convention = tl$percentile_convention %gp_or% NA_character_
    ),
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      step = "two-level-pi",
      n_units = N,
      n_industries = tl$K,
      rule = "same_industry_likelihood_once"
    ),
    warnings = character(0)
  )
  validate_gp_twolevel_pi(structure(out, class = c("gp_twolevel_pi", "list")))
}

#' Backward-compatible alias for the two-level Pi builder
#' @keywords internal
#' @noRd
gp_pairwise_twolevel <- function(posterior, ids = NULL, control = NULL) {
  gp_twolevel_pi(posterior, ids = ids, control = control)
}

#' Validate a two-level Pi object
#' @keywords internal
#' @noRd
validate_gp_twolevel_pi <- function(x) {
  if (!inherits(x, "gp_twolevel_pi")) {
    .gp_tl_abort("Expected a gp_twolevel_pi object.",
                 class = "gradepath_validation_error")
  }
  req <- c("ids", "industry", "industry_levels", "industry_representatives",
           "raw", "Pi_theta", "Pi_sbar_psi", "Pi_bar", "Pi_sq_theta",
           "Pi_xi", "Pi_psi", "Pi_sbar_psi_industry", "Pi_bar_industry",
           "pairwise", "pairwise_theta", "pairwise_bar", "source",
           "metadata", "schema_version", "provenance", "warnings")
  if (any(!req %in% names(x))) {
    .gp_tl_abort("`gp_twolevel_pi` is missing required fields.",
                 class = "gradepath_validation_error")
  }
  ids <- as.character(x$ids)
  N <- length(ids)
  K <- length(x$industry_levels)
  raw_req <- c("Pi_theta", "Pi_sbar_psi", "Pi_bar", "Pi_sq_theta", "Pi_xi", "Pi_psi")
  if (any(!raw_req %in% names(x$raw))) {
    .gp_tl_abort("`gp_twolevel_pi$raw` is missing one or more Pi matrices.",
                 class = "gradepath_validation_error")
  }
  for (nm in raw_req) {
    mat <- x$raw[[nm]]
    if (!is.matrix(mat) || !is.numeric(mat) || !identical(dim(mat), c(N, N)) ||
        any(!is.finite(mat))) {
      .gp_tl_abort("`gp_twolevel_pi$raw$%s` must be a finite N x N matrix.", nm,
                   class = "gradepath_validation_error")
    }
  }
  for (nm in c("Pi_theta", "Pi_sbar_psi", "Pi_bar", "Pi_sq_theta", "Pi_xi", "Pi_psi")) {
    mat <- x[[nm]]
    if (!is.matrix(mat) || !is.numeric(mat) || !identical(dim(mat), c(N, N)) ||
        any(!is.finite(mat))) {
      .gp_tl_abort("`gp_twolevel_pi$%s` must be a finite N x N matrix.", nm,
                   class = "gradepath_validation_error")
    }
  }
  for (nm in c("Pi_sbar_psi_industry", "Pi_bar_industry")) {
    mat <- x[[nm]]
    if (!is.matrix(mat) || !is.numeric(mat) || !identical(dim(mat), c(K, K)) ||
        any(!is.finite(mat))) {
      .gp_tl_abort("`gp_twolevel_pi$%s` must be a finite K x K matrix.", nm,
                   class = "gradepath_validation_error")
    }
  }
  if (!isTRUE(all.equal(x$raw$Pi_bar, x$raw$Pi_sbar_psi, tolerance = 0)) ||
      !isTRUE(all.equal(x$Pi_bar, x$Pi_sbar_psi, tolerance = 0)) ||
      !isTRUE(all.equal(x$Pi_bar_industry, x$Pi_sbar_psi_industry, tolerance = 0))) {
    .gp_tl_abort("`Pi_bar` fields must alias their matching `Pi_sbar_psi` fields.",
                 class = "gradepath_validation_error")
  }
  prob_raw <- c("Pi_theta", "Pi_sbar_psi", "Pi_bar", "Pi_xi", "Pi_psi")
  for (nm in prob_raw) {
    mat <- x$raw[[nm]]
    if (any(mat < 0 | mat > 1) || any(abs(diag(mat)) > 1e-12)) {
      .gp_tl_abort("`gp_twolevel_pi$raw$%s` must be a raw probability matrix with zero diagonal.",
                   nm, class = "gradepath_validation_error")
    }
  }
  if (any(x$raw$Pi_sq_theta < 0) || any(abs(diag(x$raw$Pi_sq_theta)) > 1e-12)) {
    .gp_tl_abort("`gp_twolevel_pi$raw$Pi_sq_theta` must be non-negative with zero diagonal.",
                 class = "gradepath_validation_error")
  }
  if (!isTRUE(all.equal(x$Pi_sbar_psi, x$raw$Pi_sbar_psi, tolerance = 0)) ||
      !isTRUE(all.equal(x$Pi_bar, x$raw$Pi_bar, tolerance = 0)) ||
      !isTRUE(all.equal(x$Pi_sq_theta, x$raw$Pi_sq_theta, tolerance = 0)) ||
      !isTRUE(all.equal(x$Pi_xi, x$raw$Pi_xi, tolerance = 0)) ||
      !isTRUE(all.equal(x$Pi_psi, x$raw$Pi_psi, tolerance = 0))) {
    .gp_tl_abort("Top-level raw diagnostic Pi matrices must agree with `raw`.",
                 class = "gradepath_validation_error")
  }
  validate_gp_pairwise(x$pairwise_theta)
  validate_gp_pairwise(x$pairwise)
  validate_gp_pairwise(x$pairwise_bar)
  if (!identical(x$pairwise, x$pairwise_theta)) {
    .gp_tl_abort("`gp_twolevel_pi$pairwise` must alias `pairwise_theta`.",
                 class = "gradepath_validation_error")
  }
  if (!isTRUE(all.equal(x$Pi_theta, x$pairwise_theta$matrix, tolerance = 0)) ||
      !isTRUE(all.equal(x$Pi_bar_industry, x$pairwise_bar$matrix, tolerance = 0))) {
    .gp_tl_abort("Cleaned Pi matrices must agree with their gp_pairwise wrappers.",
                 class = "gradepath_validation_error")
  }
  x
}
