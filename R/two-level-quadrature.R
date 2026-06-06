# =============================================================================
# two-level-quadrature.R  --  deterministic nested quadrature wrapper
# -----------------------------------------------------------------------------
# Formalizes the preferred M2 posterior path: outer eta quadrature, inner
# per-firm xi quadrature, and pairwise domains that encode the same-industry
# override structurally.  The low-level integrations live in the underlying routines
# (`gp_posterior_twolevel()` and `gp_twolevel_pi()`); this file provides the
# explicit quadrature orchestration and diagnostics needed by the
# fixture gate.
# =============================================================================

#' Resolve package abort helper
#' @keywords internal
#' @noRd
.gp_tlq_abort <- function(msg, ..., class = "gradepath_error") {
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

#' Numerically safe log10(sum(exp(log_terms)))
#' @keywords internal
#' @noRd
.gp_tlq_log10_sum_exp <- function(log_terms) {
  if (length(log_terms) == 0L) return(-Inf)
  m <- max(log_terms)
  if (!is.finite(m)) return(m / log(10))
  (m + log(sum(exp(log_terms - m)))) / log(10)
}

#' Complexity/state-count diagnostics for the nested quadrature path
#' @keywords internal
#' @noRd
.gp_tlq_complexity <- function(posterior) {
  tl <- .gp_tl_extract(posterior)
  N <- tl$N
  K <- tl$K
  E <- tl$E
  M <- tl$M
  sizes <- tabulate(tl$industry, nbins = K)
  same_ordered <- sum(sizes * pmax(sizes - 1L, 0L))
  total_ordered <- N * pmax(N - 1L, 0L)
  cross_ordered <- total_ordered - same_ordered

  naive_logs <- log(E) + sizes * log(M)
  naive_log10 <- .gp_tlq_log10_sum_exp(naive_logs)
  naive_states <- if (naive_log10 < 15) round(10^naive_log10) else Inf

  list(
    method = "deterministic_nested_quadrature",
    materializes_full_grid = FALSE,
    full_grid_pattern_avoided = "N_x_N_x_eta_x_eta_x_xi_x_xi_and_industry_xi_power",
    eta_nodes = E,
    xi_nodes = M,
    n_units = N,
    n_industries = K,
    industry_sizes = sizes,
    max_industry_size = max(sizes),
    posterior_unit_cells = N * E * M,
    materialized_kernel_cells = N * E * M,
    same_industry_ordered_pairs = same_ordered,
    cross_industry_ordered_pairs = cross_ordered,
    same_industry_pair_cells = same_ordered * E * M^2,
    cross_industry_pair_cells = cross_ordered * E^2 * M^2,
    naive_industry_full_grid_states = naive_states,
    naive_industry_full_grid_log10_states = naive_log10,
    forbidden_full_xi_grid_cells = naive_states,
    forbidden_full_xi_grid_log10_cells = naive_log10,
    same_industry_eta_domain = "single_shared_eta",
    cross_industry_eta_domain = "independent_eta_product",
    inner_integral = "per_firm_xi"
  )
}

#' Deterministic two-level nested quadrature fit
#'
#' @description
#' Runs or reuses the posterior, then builds the five Pi
#' matrices through deterministic nested quadrature.  This is the quadrature
#' first path described in Chapter 11: no simulation draws and no full joint grid
#' over all firms' xi states are materialized.
#'
#' @param input,prior,fit Inputs for `gp_posterior_twolevel()`. Required unless
#'   `posterior` is supplied.
#' @param posterior Optional precomputed `gp_posterior_twolevel()` result.
#' @param control Optional control list.
#' @param interval_level Optional posterior interval level passed to the posterior.
#' @param ids Optional ids passed to `gp_twolevel_pi()`.
#' @param include_g_theta Logical; when `TRUE`, also attach the
#'   two-support theta pushforward artifact.
#' @param supp_pts_theta Theta-grid size for `gp_pushforward_theta()`.
#'
#' @return A validated `gp_twolevel_quadrature` object.
#' @keywords internal
#' @noRd
gp_twolevel_quadrature <- function(input = NULL, prior = NULL, fit = NULL,
                                   posterior = NULL, control = NULL,
                                   interval_level = NULL, ids = NULL,
                                   include_g_theta = TRUE,
                                   supp_pts_theta = 250L) {
  if (is.null(posterior)) {
    if (is.null(input) || is.null(prior) || is.null(fit)) {
      .gp_tlq_abort(
        "`input`, `prior`, and `fit` are required when `posterior` is not supplied.",
        class = "gradepath_validation_error"
      )
    }
    posterior <- gp_posterior_twolevel(
      input = input,
      prior = prior,
      fit = fit,
      control = control,
      interval_level = interval_level
    )
  } else {
    posterior <- validate_gp_posterior(posterior)
  }

  pi <- gp_twolevel_pi(posterior, ids = ids, control = control)
  diagnostics <- .gp_tlq_complexity(posterior)
  diagnostics$posterior_percentile_convention <-
    posterior$metadata$two_level$percentile_convention %gp_or% NA_character_
  diagnostics$include_g_theta <- isTRUE(include_g_theta)
  diagnostics$supp_pts_theta <- as.integer(supp_pts_theta)

  g_theta <- NULL
  if (isTRUE(include_g_theta)) {
    tl <- .gp_tl_extract(posterior)
    original_s <- posterior$metadata$two_level$original_s %gp_or%
      (if (!is.null(input)) input$s else NULL)
    if (is.null(original_s)) {
      .gp_tlq_abort(
        "`posterior$metadata$two_level$original_s` or `input$s` is required to build `g_theta`.",
        class = "gradepath_validation_error"
      )
    }
    original_s <- as.numeric(original_s)
    if (length(original_s) != tl$N || any(!is.finite(original_s)) ||
        any(original_s <= 0)) {
      .gp_tlq_abort(
        "`original_s` must be a finite positive vector with one entry per unit.",
        class = "gradepath_validation_error"
      )
    }
    g_theta <- gp_pushforward_theta(
      supp_xi = tl$support_xi,
      g_xi = tl$g_xi,
      supp_eta = tl$support_eta,
      g_eta = tl$g_eta,
      s = original_s,
      mu = tl$mu,
      beta = tl$beta,
      characteristic = tl$characteristic,
      supp_pts_theta = supp_pts_theta
    )
  }

  artifacts <- list(
    posteriors = posterior$metadata$reporting$posteriors,
    Pi_theta = pi$raw$Pi_theta,
    Pi_sq_theta = pi$raw$Pi_sq_theta,
    Pi_xi = pi$raw$Pi_xi,
    Pi_psi = pi$raw$Pi_psi,
    Pi_sbar_psi = pi$raw$Pi_sbar_psi,
    Pi_bar = pi$raw$Pi_bar,
    Pi_bar_industry = pi$Pi_bar_industry,
    g_theta = g_theta
  )

  out <- list(
    posterior = posterior,
    pi = pi,
    g_theta = g_theta,
    artifacts = artifacts,
    pairwise_theta = pi$pairwise_theta,
    pairwise_bar = pi$pairwise_bar,
    method = "quadrature",
    diagnostics = diagnostics,
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      step = "two-level-quadrature",
      method = "deterministic_nested_quadrature",
      n_units = diagnostics$n_units,
      n_industries = diagnostics$n_industries,
      materializes_full_grid = FALSE
    ),
    warnings = character(0)
  )
  validate_gp_twolevel_quadrature(
    structure(out, class = c("gp_twolevel_quadrature", "list"))
  )
}

#' Validate a two-level quadrature object
#' @keywords internal
#' @noRd
validate_gp_twolevel_quadrature <- function(x) {
  if (!inherits(x, "gp_twolevel_quadrature")) {
    .gp_tlq_abort("Expected a gp_twolevel_quadrature object.",
                  class = "gradepath_validation_error")
  }
  req <- c("posterior", "pi", "g_theta", "artifacts", "pairwise_theta",
           "pairwise_bar", "method", "diagnostics", "schema_version",
           "provenance", "warnings")
  if (any(!req %in% names(x))) {
    .gp_tlq_abort("`gp_twolevel_quadrature` is missing required fields.",
                  class = "gradepath_validation_error")
  }
  validate_gp_posterior(x$posterior)
  validate_gp_twolevel_pi(x$pi)
  validate_gp_pairwise(x$pairwise_theta)
  validate_gp_pairwise(x$pairwise_bar)
  if (!is.null(x$g_theta)) {
    gt <- x$g_theta
    if (!is.list(gt) || any(!c("support", "g", "density", "diagnostics") %in% names(gt)) ||
        length(gt$support) != length(gt$g) ||
        length(gt$support) != length(gt$density) ||
        any(!is.finite(gt$support)) || any(!is.finite(gt$g)) ||
        any(!is.finite(gt$density)) || any(gt$g < 0) ||
        abs(sum(gt$g) - 1) > 1e-8 ||
        !identical(gt$diagnostics$n_carriers, x$diagnostics$n_units)) {
      .gp_tlq_abort("`gp_twolevel_quadrature$g_theta` is malformed.",
                    class = "gradepath_validation_error")
    }
  }
  art_req <- c("posteriors", "Pi_theta", "Pi_sq_theta", "Pi_xi", "Pi_psi",
               "Pi_sbar_psi", "Pi_bar", "Pi_bar_industry", "g_theta")
  if (!is.list(x$artifacts) || any(!art_req %in% names(x$artifacts))) {
    .gp_tlq_abort("`gp_twolevel_quadrature$artifacts` is incomplete.",
                  class = "gradepath_validation_error")
  }
  if (!isTRUE(all.equal(x$artifacts$posteriors,
                        x$posterior$metadata$reporting$posteriors,
                        tolerance = 0)) ||
      !isTRUE(all.equal(x$artifacts$Pi_theta, x$pi$raw$Pi_theta, tolerance = 0)) ||
      !isTRUE(all.equal(x$artifacts$Pi_sq_theta, x$pi$raw$Pi_sq_theta,
                        tolerance = 0)) ||
      !isTRUE(all.equal(x$artifacts$Pi_xi, x$pi$raw$Pi_xi, tolerance = 0)) ||
      !isTRUE(all.equal(x$artifacts$Pi_psi, x$pi$raw$Pi_psi, tolerance = 0)) ||
      !isTRUE(all.equal(x$artifacts$Pi_sbar_psi, x$pi$raw$Pi_sbar_psi,
                        tolerance = 0)) ||
      !isTRUE(all.equal(x$artifacts$Pi_bar, x$pi$raw$Pi_bar, tolerance = 0)) ||
      !isTRUE(all.equal(x$artifacts$Pi_bar_industry, x$pi$Pi_bar_industry,
                        tolerance = 0))) {
    .gp_tlq_abort("Quadrature artifacts must alias posterior/pi outputs.",
                  class = "gradepath_validation_error")
  }
  if (!identical(x$artifacts$g_theta, x$g_theta)) {
    .gp_tlq_abort("`artifacts$g_theta` must alias `g_theta`.",
                  class = "gradepath_validation_error")
  }
  if (!identical(x$pairwise_theta, x$pi$pairwise_theta) ||
      !identical(x$pairwise_bar, x$pi$pairwise_bar)) {
    .gp_tlq_abort("Quadrature pairwise aliases must match `pi` pairwise objects.",
                  class = "gradepath_validation_error")
  }
  if (!identical(x$method, "quadrature")) {
    .gp_tlq_abort("`gp_twolevel_quadrature$method` must be 'quadrature'.",
                  class = "gradepath_validation_error")
  }
  d <- x$diagnostics
  needed <- c("method", "materializes_full_grid", "eta_nodes", "xi_nodes",
              "n_units", "n_industries", "industry_sizes",
              "same_industry_eta_domain", "cross_industry_eta_domain",
              "inner_integral")
  if (!is.list(d) || any(!needed %in% names(d))) {
    .gp_tlq_abort("`gp_twolevel_quadrature$diagnostics` is incomplete.",
                  class = "gradepath_validation_error")
  }
  if (!identical(d$method, "deterministic_nested_quadrature") ||
      !identical(d$materializes_full_grid, FALSE) ||
      !identical(d$same_industry_eta_domain, "single_shared_eta") ||
      !identical(d$cross_industry_eta_domain, "independent_eta_product") ||
      !identical(d$inner_integral, "per_firm_xi")) {
    .gp_tlq_abort("Quadrature diagnostics do not satisfy the nested-domain contract.",
                  class = "gradepath_validation_error")
  }
  if (length(d$industry_sizes) != d$n_industries ||
      sum(d$industry_sizes) != d$n_units) {
    .gp_tlq_abort("Quadrature industry-size diagnostics are inconsistent.",
                  class = "gradepath_validation_error")
  }
  x
}
