# =============================================================================
# test-sim-fallback.R  --  seeded simulation fallback
# -----------------------------------------------------------------------------
# Pins the approximate fallback path for the two-level posterior.  The preferred
# path is deterministic quadrature; this file verifies the retained mnrnd-style
# seeded simulation path, including RNG hygiene, build metadata, raw draw
# percentiles, and the same-industry likelihood override.
# =============================================================================

.tls_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}

.sim2 <- function(...) .tls_get("gp_twolevel_simulation")(...)
.valid_sim2 <- function(...) .tls_get("validate_gp_twolevel_simulation")(...)
.quad2 <- function(...) .tls_get("gp_twolevel_quadrature")(...)
.push2 <- function(...) .tls_get("gp_pushforward_theta")(...)

.tiny_prior_sim <- function(characteristic = "race") {
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (identical(characteristic, "race")) {
    supp_xi <- c(0.45, 1.25)
    g_xi <- c(0.45, 0.55)
    supp_eta <- c(0.70, 1.45)
    g_eta <- c(0.40, 0.60)
  } else {
    supp_xi <- c(-0.50, 0.40)
    g_xi <- c(0.35, 0.65)
    supp_eta <- c(-0.30, 0.50)
    g_eta <- c(0.55, 0.45)
  }
  structure(
    list(
      support = supp_xi,
      density = g_xi,
      mean = sum(supp_xi * g_xi),
      scale = "r",
      diagnostics = list(group_fx = 1L, support_eta = supp_eta, g_eta = g_eta),
      metadata = list(characteristic = characteristic, beta = 0.55, mu = 0)
    ),
    class = c("gp_prior", "list")
  )
}

.tiny_case_sim <- function(characteristic = "race") {
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (identical(characteristic, "race")) {
    s <- c(0.80, 1.35, 1.05)
    beta <- 0.55
    mu <- 0
    v_hat <- c(0.72, 1.34, 1.08)
    theta_hat <- v_hat * s^beta
  } else {
    s <- c(0.90, 1.20, 1.10)
    beta <- 0.65
    mu <- 0.20
    v_hat <- c(-0.10, 0.25, 0.05)
    theta_hat <- mu + v_hat * s^beta
  }
  s_v <- c(0.30, 0.34, 0.28)
  prior <- .tiny_prior_sim(characteristic)
  prior$metadata$beta <- beta
  prior$metadata$mu <- mu
  list(
    input = list(
      theta_hat = theta_hat,
      s = s,
      industry = c(1, 1, 2),
      unit_id = paste0(characteristic, "_", 1:3),
      label = paste(characteristic, 1:3)
    ),
    prior = prior,
    fit = list(
      characteristic = characteristic,
      beta = beta,
      mu = mu,
      v_hat = v_hat,
      s_v = s_v,
      industry = c(1, 1, 2)
    )
  )
}

.wlog <- function(log_w) {
  w <- exp(log_w - max(log_w))
  w / sum(w)
}

.manual_sim_pair <- function(sim, i, j) {
  draws <- sim$draws
  log_w <- if (sim$pi$industry[i] == sim$pi$industry[j]) {
    draws$log_L[i, ]
  } else {
    draws$log_L[i, ] + draws$log_L[j, ]
  }
  w <- .wlog(log_w)
  theta_gap <- draws$theta[i, ] - draws$theta[j, ]
  c(
    Pi_theta = sum(w * (theta_gap > 0)),
    Pi_sq_theta = sum(w * pmax(theta_gap, 0)^2),
    Pi_xi = sum(w * (draws$xi[i, ] > draws$xi[j, ])),
    Pi_psi = sum(w * (draws$psi[i, ] > draws$psi[j, ])),
    Pi_sbar_psi = sum(
      w * ((draws$s_bar[i, ] * draws$psi[i, ]) >
             (draws$s_bar[j, ] * draws$psi[j, ]))
    )
  )
}

test_that("seeded simulation fallback is reproducible and preserves caller RNG", {
  b <- .tiny_case_sim("race")
  set.seed(8675309)
  before <- .Random.seed
  s1 <- .sim2(b$input, b$prior, b$fit, n_draws = 80L, seed = 1234L,
              keep_draws = TRUE)
  after <- .Random.seed
  expect_identical(after, before)

  s2 <- .sim2(b$input, b$prior, b$fit, n_draws = 80L, seed = 1234L,
              keep_draws = TRUE)
  s3 <- .sim2(b$input, b$prior, b$fit, n_draws = 80L, seed = 4321L,
              keep_draws = TRUE)

  expect_s3_class(s1, "gp_twolevel_simulation")
  expect_silent(.valid_sim2(s1))
  expect_equal(s1$artifacts$posteriors, s2$artifacts$posteriors, tolerance = 0)
  expect_equal(s1$artifacts$Pi_theta, s2$artifacts$Pi_theta, tolerance = 0)
  expect_false(isTRUE(all.equal(s1$draws$draw_seeds, s3$draws$draw_seeds,
                                tolerance = 0)))
})

test_that("seeded simulation fallback preserves an absent RNG state", {
  b <- .tiny_case_sim("race")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  sim <- .sim2(b$input, b$prior, b$fit, n_draws = 20L, seed = 1234L,
               include_g_theta = FALSE)
  expect_silent(.valid_sim2(sim))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("simulation metadata records approximate status and build seed", {
  b <- .tiny_case_sim("gender")
  sim <- .sim2(b$input, b$prior, b$fit, n_draws = 70L, seed = 99L,
               include_g_theta = FALSE)
  bm <- attr(sim, "build_metadata")

  expect_identical(sim$method, "simulate")
  expect_identical(sim$diagnostics$method, "seeded_simulation_fallback")
  expect_identical(sim$diagnostics$producer_status, "APPROXIMATE_OK")
  expect_identical(sim$diagnostics$tolerance_class, "approximate")
  expect_identical(sim$diagnostics$simulation_engine, "R_multinomial_importance")
  expect_identical(sim$diagnostics$seed, 99L)
  expect_identical(sim$diagnostics$n_draws, 70L)
  expect_identical(sim$diagnostics$matlab_rng_parity, FALSE)
  expect_identical(sim$metadata$producer_status, "APPROXIMATE_OK")
  expect_identical(sim$metadata$tolerance_class, "approximate")
  expect_identical(bm$seed, 99L)
  expect_identical(bm$extra$producer_status, "APPROXIMATE_OK")
  expect_identical(bm$extra$tolerance_class, "approximate")
  expect_identical(bm$extra$matlab_rng_parity, FALSE)
  expect_null(sim$g_theta)
  expect_null(sim$artifacts$g_theta)
})

test_that("simulation uses raw length-R draw percentiles and posterior weights", {
  b <- .tiny_case_sim("race")
  sim <- .sim2(b$input, b$prior, b$fit, n_draws = 90L, seed = 111L,
               interval_level = 0.95, keep_draws = TRUE)
  wpr <- .tls_get(".gp_weighted_percentile")
  theta <- sim$draws$theta[1, ]
  w <- .wlog(sim$draws$log_L[1, ])
  q <- wpr(theta, c(2.5, 97.5), w)

  expect_identical(sim$posterior$metadata$reporting$percentile_convention,
                   "simulation_wprctile_type5")
  expect_equal(sim$posterior$metadata$reporting$posterior_mean[1],
               sum(theta * w), tolerance = 1e-12)
  expect_equal(sim$posterior$metadata$reporting$posterior_second_moment[1],
               sum(theta^2 * w), tolerance = 1e-12)
  expect_equal(sim$artifacts$posteriors$lower[1], q[1], tolerance = 1e-12)
  expect_equal(sim$artifacts$posteriors$upper[1], q[2], tolerance = 1e-12)
})

test_that("simulation five Pi matrices use the same-industry L_i override", {
  b <- .tiny_case_sim("race")
  sim <- .sim2(b$input, b$prior, b$fit, n_draws = 120L, seed = 222L,
               keep_draws = TRUE)
  same <- .manual_sim_pair(sim, 1, 2)
  cross <- .manual_sim_pair(sim, 1, 3)

  for (nm in names(same)) {
    expect_equal(sim$pi$raw[[nm]][1, 2], unname(same[nm]), tolerance = 1e-12)
    expect_equal(sim$pi$raw[[nm]][1, 3], unname(cross[nm]), tolerance = 1e-12)
  }
  expect_equal(sim$pi$raw$Pi_psi[1, 2], 0, tolerance = 0)
  expect_equal(sim$pi$raw$Pi_sbar_psi[1, 2], 0, tolerance = 0)

  wrong_w <- .wlog(sim$draws$log_L[1, ] + sim$draws$log_L[2, ])
  wrong <- sum(wrong_w * (sim$draws$theta[1, ] > sim$draws$theta[2, ]))
  expect_gt(abs(sim$pi$raw$Pi_theta[1, 2] - wrong), 1e-3)
})

test_that("simulation artifacts are shape-compatible with quadrature facade", {
  b <- .tiny_case_sim("gender")
  sim <- .sim2(b$input, b$prior, b$fit, n_draws = 75L, seed = 333L,
               supp_pts_theta = 31L)
  quad <- .quad2(b$input, b$prior, b$fit, supp_pts_theta = 31L)

  expect_named(sim$artifacts, names(quad$artifacts))
  expect_equal(dim(sim$artifacts$posteriors), dim(quad$artifacts$posteriors))
  expect_equal(dim(sim$artifacts$Pi_theta), dim(quad$artifacts$Pi_theta))
  expect_equal(dim(sim$artifacts$Pi_bar_industry),
               dim(quad$artifacts$Pi_bar_industry))
  expect_equal(length(sim$artifacts$g_theta$support),
               length(quad$artifacts$g_theta$support))
  expect_equal(sim$artifacts$g_theta$g, quad$artifacts$g_theta$g, tolerance = 0)
  expect_silent(.tls_get("validate_gp_pairwise")(sim$pairwise_theta))
  expect_silent(.tls_get("validate_gp_pairwise")(sim$pairwise_bar))
})

test_that("simulation validator catches malformed objects and bad controls", {
  b <- .tiny_case_sim("race")
  expect_error(.sim2(b$input, b$prior, b$fit, n_draws = 0),
               class = "gradepath_validation_error")
  expect_error(.sim2(b$input, b$prior, b$fit, n_draws = 10, seed = -1),
               class = "gradepath_validation_error")

  sim <- .sim2(b$input, b$prior, b$fit, n_draws = 50L, seed = 444L)
  bad <- sim
  bad$diagnostics$producer_status <- "OK"
  expect_error(.valid_sim2(bad), class = "gradepath_validation_error")

  bad2 <- sim
  bad2$artifacts$Pi_xi <- bad2$artifacts$Pi_xi + 0.01
  expect_error(.valid_sim2(bad2), class = "gradepath_validation_error")

  bad3 <- sim
  bad3$metadata$build_metadata$seed <- 1L
  expect_error(.valid_sim2(bad3), class = "gradepath_validation_error")
})
