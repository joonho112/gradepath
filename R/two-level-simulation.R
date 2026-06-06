# =============================================================================
# two-level-simulation.R  --  seeded simulation fallback
# -----------------------------------------------------------------------------
# Retains the KRW get_posteriors.m:93-199 simulation path as an explicit
# approximate fallback for the two-level posterior.  This path intentionally does
# not claim Matlab RNG parity: it uses R's seeded categorical sampler in the same
# mnrnd-style dimensions, records the seed in build metadata, and marks the
# producer status APPROXIMATE_OK.  The default claim-bearing path remains the
# deterministic quadrature facade.
# =============================================================================

#' Resolve package abort helper
#' @keywords internal
#' @noRd
.gp_tls_abort <- function(msg, ..., class = "gradepath_error") {
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

#' Normalize and validate simulation controls
#' @keywords internal
#' @noRd
.gp_tls_controls <- function(control = NULL, n_draws = NULL, seed = NULL,
                             interval_level = NULL) {
  ctrl <- if (is.null(control)) gp_control() else validate_gp_control(control)
  if (is.null(n_draws)) {
    n_draws <- ctrl$sims %gp_or% ctrl$twolevel_sims %gp_or% 10000L
  }
  n_draws <- as.numeric(n_draws)
  if (length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws < 1 || abs(n_draws - round(n_draws)) > 1e-8) {
    .gp_tls_abort("`n_draws` must be a positive integer.",
                  class = "gradepath_validation_error")
  }
  n_draws <- as.integer(round(n_draws))

  if (is.null(seed)) seed <- ctrl$seed %gp_or% 1234L
  seed <- as.numeric(seed)
  if (length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || abs(seed - round(seed)) > 1e-8) {
    .gp_tls_abort("`seed` must be a non-negative integer.",
                  class = "gradepath_validation_error")
  }
  seed <- as.integer(round(seed))

  level <- interval_level %gp_or% ctrl$interval_level %gp_or% 0.90
  level <- as.numeric(level)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    .gp_tls_abort("`interval_level` must be a scalar in (0, 1).",
                  class = "gradepath_validation_error")
  }

  list(control = ctrl, n_draws = n_draws, seed = seed,
       interval_level = level)
}

#' Preserve caller RNG state while running a seeded simulation
#' @keywords internal
#' @noRd
.gp_tls_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

#' Softmax-normalized likelihood weights from log likelihoods
#' @keywords internal
#' @noRd
.gp_tls_log_weights <- function(log_w, name = "log_w") {
  log_w <- as.numeric(log_w)
  if (length(log_w) == 0L || any(!is.finite(log_w))) {
    .gp_tls_abort("`%s` must contain finite log weights.", name,
                  class = "gradepath_validation_error")
  }
  m <- max(log_w)
  w <- exp(log_w - m)
  total <- sum(w)
  if (!is.finite(total) || total <= 0) {
    .gp_tls_abort("`%s` has degenerate total mass.", name,
                  class = "gradepath_validation_error")
  }
  w / total
}

#' Seeded mnrnd-style draws and broadcast industry likelihoods
#' @keywords internal
#' @noRd
.gp_tls_draws <- function(posterior, n_draws, seed) {
  tl <- .gp_tl_extract(posterior)
  v_hat <- as.numeric(posterior$estimate)
  s_v <- as.numeric(posterior$se)
  if (length(v_hat) != tl$N || length(s_v) != tl$N ||
      any(!is.finite(v_hat)) || any(!is.finite(s_v)) || any(s_v <= 0)) {
    .gp_tls_abort("`posterior` must carry finite v_hat/s_v vectors.",
                  class = "gradepath_validation_error")
  }

  .gp_tls_with_seed(seed, {
    draw_seeds <- sample.int(n_draws, n_draws, replace = TRUE)
    xi <- psi <- theta <- r_effect <- log_L <- matrix(NA_real_, tl$N, n_draws)
    s_bar <- matrix(rep(tl$sbar, n_draws), nrow = tl$N, ncol = n_draws)

    for (r in seq_len(n_draws)) {
      set.seed(draw_seeds[r])
      xi_idx <- sample.int(tl$M, tl$N, replace = TRUE, prob = tl$g_xi)
      eta_idx <- sample.int(tl$E, tl$K, replace = TRUE, prob = tl$g_eta)
      xi_r <- tl$support_xi[xi_idx]
      eta_by_industry <- tl$support_eta[eta_idx]
      psi_r <- eta_by_industry[tl$industry]

      xi[, r] <- xi_r
      psi[, r] <- psi_r
      if (identical(tl$characteristic, "race")) {
        r_effect[, r] <- psi_r * xi_r
        theta[, r] <- tl$s_beta * r_effect[, r]
      } else {
        r_effect[, r] <- psi_r + xi_r
        theta[, r] <- tl$mu + tl$s_beta * r_effect[, r]
      }

      for (k in seq_len(tl$K)) {
        idx <- which(tl$industry == k)
        mean_k <- if (identical(tl$characteristic, "race")) {
          psi_r[idx] * xi_r[idx]
        } else {
          psi_r[idx] + xi_r[idx]
        }
        lk <- sum(stats::dnorm(v_hat[idx], mean = mean_k, sd = s_v[idx],
                               log = TRUE))
        log_L[idx, r] <- lk
      }
    }

    list(
      xi = xi,
      psi = psi,
      s_bar = s_bar,
      r_effect = r_effect,
      theta = theta,
      log_L = log_L,
      draw_seeds = draw_seeds
    )
  })
}

#' Summarize simulation posterior draws unit by unit
#' @keywords internal
#' @noRd
.gp_tls_posteriors <- function(posterior, draws, interval_level) {
  N <- nrow(draws$theta)
  pct <- c((1 - interval_level) / 2, (1 + interval_level) / 2) * 100
  pm_t <- e2_t <- psd_t <- lo_t <- up_t <- numeric(N)
  pm_r <- e2_r <- psd_r <- lo_r <- up_r <- numeric(N)
  e_xi <- e_eta <- e_sbar_eta <- lo_sbar_eta <- up_sbar_eta <- numeric(N)
  for (i in seq_len(N)) {
    w <- .gp_tls_log_weights(draws$log_L[i, ], "unit log likelihood")
    theta_i <- draws$theta[i, ]
    r_i <- draws$r_effect[i, ]
    sbar_eta_i <- draws$s_bar[i, ] * draws$psi[i, ]

    pm_t[i] <- sum(w * theta_i)
    e2_t[i] <- sum(w * theta_i^2)
    psd_t[i] <- sqrt(max(e2_t[i] - pm_t[i]^2, 0))
    q_t <- .gp_weighted_percentile(theta_i, pct, w)
    lo_t[i] <- q_t[1L]
    up_t[i] <- q_t[2L]

    pm_r[i] <- sum(w * r_i)
    e2_r[i] <- sum(w * r_i^2)
    psd_r[i] <- sqrt(max(e2_r[i] - pm_r[i]^2, 0))
    q_r <- .gp_weighted_percentile(r_i, pct, w)
    lo_r[i] <- min(q_r[1L], pm_r[i])
    up_r[i] <- max(q_r[2L], pm_r[i])

    e_xi[i] <- sum(w * draws$xi[i, ])
    e_eta[i] <- sum(w * draws$psi[i, ])
    e_sbar_eta[i] <- sum(w * sbar_eta_i)
    q_sbar <- .gp_weighted_percentile(sbar_eta_i, pct, w)
    lo_sbar_eta[i] <- q_sbar[1L]
    up_sbar_eta[i] <- q_sbar[2L]
  }

  list(
    reporting = list(
      posterior_mean = pm_t,
      posterior_sd = psd_t,
      posterior_second_moment = e2_t,
      lower = lo_t,
      upper = up_t,
      posteriors = data.frame(
        posterior_mean = pm_t,
        posterior_second_moment = e2_t,
        lower = lo_t,
        upper = up_t
      ),
      scale = "theta",
      level = interval_level,
      percentile_convention = "simulation_wprctile_type5"
    ),
    components = data.frame(
      E_r = pm_r,
      E_r2 = e2_r,
      E_theta = pm_t,
      E_theta2 = e2_t,
      E_xi = e_xi,
      E_eta = e_eta,
      E_sbar_eta = e_sbar_eta,
      lower_sbar_eta = lo_sbar_eta,
      upper_sbar_eta = up_sbar_eta
    ),
    shell = list(
      pm_r = pm_r,
      psd_r = psd_r,
      lo_r = lo_r,
      up_r = up_r,
      e2_r = e2_r
    )
  )
}

#' Build simulation Pi matrices from stored draws
#' @keywords internal
#' @noRd
.gp_tls_pi <- function(base_posterior, draws, ids, control, seed, n_draws) {
  tl <- .gp_tl_extract(base_posterior)
  raw <- list(
    Pi_theta = matrix(0, tl$N, tl$N),
    Pi_sbar_psi = matrix(0, tl$N, tl$N),
    Pi_bar = matrix(0, tl$N, tl$N),
    Pi_sq_theta = matrix(0, tl$N, tl$N),
    Pi_xi = matrix(0, tl$N, tl$N),
    Pi_psi = matrix(0, tl$N, tl$N)
  )

  for (i in seq_len(tl$N)) {
    for (j in seq_len(tl$N)) {
      if (i == j) next
      log_w <- if (tl$industry[i] == tl$industry[j]) {
        draws$log_L[i, ]
      } else {
        draws$log_L[i, ] + draws$log_L[j, ]
      }
      w <- .gp_tls_log_weights(log_w, "pair log likelihood")
      theta_gap <- draws$theta[i, ] - draws$theta[j, ]
      raw$Pi_theta[i, j] <- sum(w * (theta_gap > 0))
      raw$Pi_sq_theta[i, j] <- sum(w * pmax(theta_gap, 0)^2)
      raw$Pi_xi[i, j] <- sum(w * (draws$xi[i, ] > draws$xi[j, ]))
      raw$Pi_psi[i, j] <- sum(w * (draws$psi[i, ] > draws$psi[j, ]))
      raw$Pi_sbar_psi[i, j] <- sum(
        w * ((draws$s_bar[i, ] * draws$psi[i, ]) >
               (draws$s_bar[j, ] * draws$psi[j, ]))
      )
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
    control = control,
    n_industries = tl$K,
    provenance_step = "two-level-simulation-pairwise-theta"
  )
  pairwise_bar <- .gp_tl_new_pairwise(
    ids = bar_ids,
    matrix = pi_bar_clean,
    control = control,
    n_industries = tl$K,
    provenance_step = "two-level-simulation-pairwise-sbar-psi"
  )

  obj <- list(
    ids = ids,
    industry = tl$industry,
    industry_levels = bar_ids,
    industry_representatives = reps,
    raw = raw,
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
      same_industry = "seeded_simulation_L_i",
      cross_industry = "seeded_simulation_L_i_times_L_j",
      producer_status = "APPROXIMATE_OK"
    ),
    metadata = list(
      industry = tl$industry,
      industry_levels = bar_ids,
      raw_Pi_sbar_psi_industry = Pi_sbar_psi_raw,
      characteristic = tl$characteristic,
      percentile_convention = "simulation_wprctile_type5",
      simulation = list(seed = seed, n_draws = n_draws)
    ),
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      step = "two-level-simulation-pi",
      n_units = tl$N,
      n_industries = tl$K,
      seed = seed,
      n_draws = n_draws,
      rule = "same_industry_likelihood_once"
    ),
    warnings = "APPROXIMATE_OK: seeded R simulation fallback, not Matlab RNG parity"
  )
  validate_gp_twolevel_pi(structure(obj, class = c("gp_twolevel_pi", "list")))
}

#' Build metadata for the simulation fallback
#' @keywords internal
#' @noRd
.gp_tls_build_metadata <- function(seed, n_draws, extra = list()) {
  list(
    key = "twolevel-simulation-fallback",
    seed = seed,
    built_at = Sys.time(),
    r_version = paste0(R.version$major, ".", R.version$minor),
    gradepath_version = .gradepath_package_version(),
    backend = NA_character_,
    solver_metadata = list(
      backend_env = NA_character_,
      gurobi_cl_available = NA,
      gurobi_version = NA_character_,
      gurobi_r_package_available = NA,
      highs_package_available = NA,
      roi_package_available = NA
    ),
    source_hash = NA_character_,
    extra = c(list(
      method = "simulate",
      producer_status = "APPROXIMATE_OK",
      tolerance_class = "approximate",
      n_draws = n_draws,
      rng_kind = RNGkind(),
      seed_stream = "R_sample_int_draw_seeds_then_categorical_samples",
      matlab_rng_parity = FALSE
    ), extra)
  )
}

#' Seeded two-level simulation fallback
#'
#' @description
#' Runs the KRW `get_posteriors.m:93-199` simulation layout in R: draw xi for
#' every unit and eta for every industry, broadcast each industry's likelihood to
#' its member rows, compute weighted posterior summaries, and assemble the same
#' five Pi matrices with the same-industry `L_ij = L_i` override.  The result is
#' marked `APPROXIMATE_OK` because R cannot bit-match Matlab's `mnrnd`/twister
#' stream.
#'
#' @return A validated `gp_twolevel_simulation` object.
#' @keywords internal
#' @noRd
gp_twolevel_simulation <- function(input = NULL, prior = NULL, fit = NULL,
                                   posterior = NULL, control = NULL,
                                   interval_level = NULL, ids = NULL,
                                   n_draws = NULL, seed = NULL,
                                   include_g_theta = TRUE,
                                   supp_pts_theta = 250L,
                                   keep_draws = FALSE) {
  ctl <- .gp_tls_controls(
    control = control,
    n_draws = n_draws,
    seed = seed,
    interval_level = interval_level
  )
  if (is.null(posterior)) {
    if (is.null(input) || is.null(prior) || is.null(fit)) {
      .gp_tls_abort(
        "`input`, `prior`, and `fit` are required when `posterior` is not supplied.",
        class = "gradepath_validation_error"
      )
    }
    posterior <- gp_posterior_twolevel(
      input = input,
      prior = prior,
      fit = fit,
      control = ctl$control,
      interval_level = ctl$interval_level
    )
  } else {
    posterior <- validate_gp_posterior(posterior)
  }
  tl <- .gp_tl_extract(posterior)
  ids <- as.character(ids %gp_or% posterior$id)
  if (length(ids) != tl$N || anyNA(ids) || any(duplicated(ids))) {
    .gp_tls_abort("`ids` must be a unique character vector of length N.",
                  class = "gradepath_validation_error")
  }

  draws <- .gp_tls_draws(posterior, ctl$n_draws, ctl$seed)
  summaries <- .gp_tls_posteriors(posterior, draws, ctl$interval_level)
  sim_tl <- posterior$metadata$two_level
  sim_tl$posterior_components <- summaries$components
  sim_tl$percentile_convention <- "simulation_wprctile_type5"
  sim_tl$simulation <- list(
    seed = ctl$seed,
    n_draws = ctl$n_draws,
    producer_status = "APPROXIMATE_OK",
    tolerance_class = "approximate",
    matlab_rng_parity = FALSE,
    likelihood = "log_likelihood_softmax"
  )

  sim_posterior <- new_gp_posterior(
    estimate = posterior$estimate,
    se = posterior$se,
    id = ids,
    label = posterior$label,
    posterior_mean = summaries$shell$pm_r,
    posterior_sd = summaries$shell$psd_r,
    lower = summaries$shell$lo_r,
    upper = summaries$shell$up_r,
    scale = "r",
    metadata = list(
      level = ctl$interval_level,
      interval_level = ctl$interval_level,
      reporting = summaries$reporting,
      has_reporting = TRUE,
      r_second_moment = summaries$shell$e2_r,
      two_level = sim_tl
    )
  )
  sim_posterior <- validate_gp_posterior(sim_posterior)

  pi <- .gp_tls_pi(posterior, draws, ids = ids, control = ctl$control,
                   seed = ctl$seed, n_draws = ctl$n_draws)

  g_theta <- NULL
  if (isTRUE(include_g_theta)) {
    original_s <- posterior$metadata$two_level$original_s %gp_or%
      (if (!is.null(input)) input$s else NULL)
    if (is.null(original_s)) {
      .gp_tls_abort(
        "`posterior$metadata$two_level$original_s` or `input$s` is required to build `g_theta`.",
        class = "gradepath_validation_error"
      )
    }
    original_s <- as.numeric(original_s)
    if (length(original_s) != tl$N || any(!is.finite(original_s)) ||
        any(original_s <= 0)) {
      .gp_tls_abort(
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

  build_metadata <- .gp_tls_build_metadata(
    seed = ctl$seed,
    n_draws = ctl$n_draws,
    extra = list(
      interval_level = ctl$interval_level,
      n_units = tl$N,
      n_industries = tl$K,
      characteristic = tl$characteristic
    )
  )
  artifacts <- list(
    posteriors = sim_posterior$metadata$reporting$posteriors,
    Pi_theta = pi$raw$Pi_theta,
    Pi_sq_theta = pi$raw$Pi_sq_theta,
    Pi_xi = pi$raw$Pi_xi,
    Pi_psi = pi$raw$Pi_psi,
    Pi_sbar_psi = pi$raw$Pi_sbar_psi,
    Pi_bar = pi$raw$Pi_bar,
    Pi_bar_industry = pi$Pi_bar_industry,
    g_theta = g_theta
  )
  diagnostics <- list(
    method = "seeded_simulation_fallback",
    simulation_engine = "R_multinomial_importance",
    producer_status = "APPROXIMATE_OK",
    tolerance_class = "approximate",
    n_draws = ctl$n_draws,
    seed = ctl$seed,
    rng_kind = RNGkind(),
    n_units = tl$N,
    n_industries = tl$K,
    industry_sizes = tabulate(tl$industry, nbins = tl$K),
    eta_nodes = tl$E,
    xi_nodes = tl$M,
    materializes_full_grid = FALSE,
    materializes_draw_matrices = TRUE,
    likelihood_scale = "log_softmax",
    same_industry_likelihood = "L_i",
    cross_industry_likelihood = "L_i_times_L_j",
    percentile_convention = "simulation_wprctile_type5",
    include_g_theta = isTRUE(include_g_theta),
    supp_pts_theta = as.integer(supp_pts_theta),
    matlab_rng_parity = FALSE
  )

  out <- list(
    posterior = sim_posterior,
    pi = pi,
    g_theta = g_theta,
    artifacts = artifacts,
    pairwise_theta = pi$pairwise_theta,
    pairwise_bar = pi$pairwise_bar,
    method = "simulate",
    diagnostics = diagnostics,
    metadata = list(
      producer_status = "APPROXIMATE_OK",
      tolerance_class = "approximate",
      build_metadata = build_metadata,
      simulation = sim_tl$simulation
    ),
    draws = if (isTRUE(keep_draws)) draws else NULL,
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      step = "two-level-simulation",
      method = "seeded_simulation_fallback",
      seed = ctl$seed,
      n_draws = ctl$n_draws,
      n_units = tl$N,
      n_industries = tl$K
    ),
    warnings = "APPROXIMATE_OK: seeded R simulation fallback, not Matlab RNG parity"
  )
  out <- validate_gp_twolevel_simulation(
    structure(out, class = c("gp_twolevel_simulation", "list"))
  )
  attr(out, "build_metadata") <- build_metadata
  out
}

#' Backward-compatible verb alias for the simulation fallback
#' @keywords internal
#' @noRd
gp_twolevel_simulate <- function(...) {
  gp_twolevel_simulation(...)
}

#' Validate a two-level simulation object
#' @keywords internal
#' @noRd
validate_gp_twolevel_simulation <- function(x) {
  if (!inherits(x, "gp_twolevel_simulation")) {
    .gp_tls_abort("Expected a gp_twolevel_simulation object.",
                  class = "gradepath_validation_error")
  }
  req <- c("posterior", "pi", "g_theta", "artifacts", "pairwise_theta",
           "pairwise_bar", "method", "diagnostics", "metadata", "draws",
           "schema_version", "provenance", "warnings")
  if (any(!req %in% names(x))) {
    .gp_tls_abort("`gp_twolevel_simulation` is missing required fields.",
                  class = "gradepath_validation_error")
  }
  validate_gp_posterior(x$posterior)
  validate_gp_twolevel_pi(x$pi)
  validate_gp_pairwise(x$pairwise_theta)
  validate_gp_pairwise(x$pairwise_bar)
  if (!identical(x$pairwise_theta, x$pi$pairwise_theta) ||
      !identical(x$pairwise_bar, x$pi$pairwise_bar)) {
    .gp_tls_abort("Simulation pairwise aliases must match `pi` pairwise objects.",
                  class = "gradepath_validation_error")
  }
  if (!identical(x$method, "simulate")) {
    .gp_tls_abort("`gp_twolevel_simulation$method` must be 'simulate'.",
                  class = "gradepath_validation_error")
  }
  d <- x$diagnostics
  needed <- c("method", "producer_status", "tolerance_class", "n_draws",
              "seed", "n_units", "n_industries", "industry_sizes",
              "materializes_full_grid",
              "same_industry_likelihood", "cross_industry_likelihood",
              "matlab_rng_parity")
  if (!is.list(d) || any(!needed %in% names(d))) {
    .gp_tls_abort("`gp_twolevel_simulation$diagnostics` is incomplete.",
                  class = "gradepath_validation_error")
  }
  if (!identical(d$method, "seeded_simulation_fallback") ||
      !identical(d$producer_status, "APPROXIMATE_OK") ||
      !identical(d$tolerance_class, "approximate") ||
      !identical(d$materializes_full_grid, FALSE) ||
      !identical(d$same_industry_likelihood, "L_i") ||
      !identical(d$cross_industry_likelihood, "L_i_times_L_j") ||
      !identical(d$matlab_rng_parity, FALSE)) {
    .gp_tls_abort("Simulation diagnostics do not satisfy the fallback contract.",
                  class = "gradepath_validation_error")
  }
  if (length(d$industry_sizes) != d$n_industries ||
      sum(d$industry_sizes) != d$n_units ||
      d$n_draws < 1L || d$seed < 0L) {
    .gp_tls_abort("Simulation diagnostics have inconsistent sizes or controls.",
                  class = "gradepath_validation_error")
  }
  art_req <- c("posteriors", "Pi_theta", "Pi_sq_theta", "Pi_xi", "Pi_psi",
               "Pi_sbar_psi", "Pi_bar", "Pi_bar_industry", "g_theta")
  if (!is.list(x$artifacts) || any(!art_req %in% names(x$artifacts))) {
    .gp_tls_abort("`gp_twolevel_simulation$artifacts` is incomplete.",
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
    .gp_tls_abort("Simulation artifacts must alias posterior/pi outputs.",
                  class = "gradepath_validation_error")
  }
  if (!identical(x$artifacts$g_theta, x$g_theta)) {
    .gp_tls_abort("`artifacts$g_theta` must alias `g_theta`.",
                  class = "gradepath_validation_error")
  }
  if (!is.null(x$g_theta)) {
    gt <- x$g_theta
    if (!is.list(gt) || any(!c("support", "g", "density", "diagnostics") %in% names(gt)) ||
        length(gt$support) != length(gt$g) ||
        length(gt$support) != length(gt$density) ||
        any(!is.finite(gt$support)) || any(!is.finite(gt$g)) ||
        any(!is.finite(gt$density)) || any(gt$g < 0) ||
        abs(sum(gt$g) - 1) > 1e-8 ||
        !identical(gt$diagnostics$n_carriers, d$n_units)) {
      .gp_tls_abort("`gp_twolevel_simulation$g_theta` is malformed.",
                    class = "gradepath_validation_error")
    }
  }
  bm <- x$metadata$build_metadata %gp_or% attr(x, "build_metadata", exact = TRUE)
  if (!is.list(bm) || !identical(bm$seed, d$seed) ||
      !identical(bm$extra$producer_status, "APPROXIMATE_OK") ||
      !identical(bm$extra$tolerance_class, "approximate") ||
      !identical(bm$extra$matlab_rng_parity, FALSE)) {
    .gp_tls_abort("Simulation build metadata must record seed and APPROXIMATE_OK status.",
                  class = "gradepath_validation_error")
  }
  if (!is.null(x$draws)) {
    draw_req <- c("xi", "psi", "s_bar", "r_effect", "theta", "log_L", "draw_seeds")
    if (!is.list(x$draws) || any(!draw_req %in% names(x$draws))) {
      .gp_tls_abort("`gp_twolevel_simulation$draws` is malformed.",
                    class = "gradepath_validation_error")
    }
    for (nm in setdiff(draw_req, "draw_seeds")) {
      mat <- x$draws[[nm]]
      if (!is.matrix(mat) || !identical(dim(mat), c(d$n_units, d$n_draws)) ||
          any(!is.finite(mat))) {
        .gp_tls_abort("`gp_twolevel_simulation$draws$%s` has invalid dimensions.",
                      nm, class = "gradepath_validation_error")
      }
    }
    if (length(x$draws$draw_seeds) != d$n_draws ||
        any(!is.finite(x$draws$draw_seeds))) {
      .gp_tls_abort("`gp_twolevel_simulation$draws$draw_seeds` is invalid.",
                    class = "gradepath_validation_error")
    }
  }
  x
}
