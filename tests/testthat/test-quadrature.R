# =============================================================================
# test-quadrature.R  --  deterministic nested quadrature wrapper
# -----------------------------------------------------------------------------
# The math is exercised elsewhere.  This file pins the formal quadrature
# path: it orchestrates posterior + Pi construction, records nested-domain
# diagnostics, and refuses to masquerade as a full-grid or simulation path.
# =============================================================================

.tlq_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}

.post2q <- function(...) .tlq_get("gp_posterior_twolevel")(...)
.pi2q <- function(...) .tlq_get("gp_twolevel_pi")(...)
.quad <- function(...) .tlq_get("gp_twolevel_quadrature")(...)
.valid_quad <- function(...) .tlq_get("validate_gp_twolevel_quadrature")(...)
.valid_pi <- function(...) .tlq_get("validate_gp_twolevel_pi")(...)
.push <- function(...) .tlq_get("gp_pushforward_theta")(...)

.tiny_prior_q <- function(characteristic = "race") {
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

.tiny_case_q <- function(characteristic = "race") {
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
  prior <- .tiny_prior_q(characteristic)
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

test_that("quadrature wrapper reproduces direct outputs", {
  b <- .tiny_case_q("race")
  post <- .post2q(b$input, b$prior, b$fit, interval_level = 0.95)
  pi <- .pi2q(post)
  q <- .quad(input = b$input, prior = b$prior, fit = b$fit,
             interval_level = 0.95)

  expect_s3_class(q, "gp_twolevel_quadrature")
  expect_silent(.valid_quad(q))
  expect_equal(q$posterior$metadata$reporting$posteriors,
               post$metadata$reporting$posteriors,
               tolerance = 0)
  expect_equal(q$pi$raw$Pi_theta, pi$raw$Pi_theta, tolerance = 0)
  expect_equal(q$pi$raw$Pi_sbar_psi, pi$raw$Pi_sbar_psi, tolerance = 0)
  expect_equal(q$pi$Pi_theta, pi$Pi_theta, tolerance = 0)
  expect_identical(q$pairwise_theta, q$pi$pairwise_theta)
  expect_identical(q$pairwise_bar, q$pi$pairwise_bar)
  expect_equal(q$artifacts$posteriors,
               post$metadata$reporting$posteriors,
               tolerance = 0)
  expect_equal(q$artifacts$Pi_theta, pi$raw$Pi_theta, tolerance = 0)
  expect_equal(q$artifacts$Pi_bar, pi$raw$Pi_bar, tolerance = 0)
})

test_that("quadrature can consume a precomputed posterior", {
  b <- .tiny_case_q("gender")
  post <- .post2q(b$input, b$prior, b$fit)
  q <- .quad(posterior = post)

  expect_silent(.valid_quad(q))
  expect_silent(.valid_pi(q$pi))
  expect_equal(q$posterior$metadata$reporting$posterior_mean,
               post$metadata$reporting$posterior_mean,
               tolerance = 0)
  expect_identical(q$method, "quadrature")
  expect_false(is.null(q$g_theta))
})

test_that("quadrature diagnostics pin nested domains and avoid full grids", {
  b <- .tiny_case_q("race")
  q <- .quad(input = b$input, prior = b$prior, fit = b$fit)
  d <- q$diagnostics

  expect_identical(d$method, "deterministic_nested_quadrature")
  expect_false(d$materializes_full_grid)
  expect_identical(d$same_industry_eta_domain, "single_shared_eta")
  expect_identical(d$cross_industry_eta_domain, "independent_eta_product")
  expect_identical(d$inner_integral, "per_firm_xi")
  expect_equal(d$n_units, 3L)
  expect_equal(d$n_industries, 2L)
  expect_equal(d$eta_nodes, 2L)
  expect_equal(d$xi_nodes, 2L)
  expect_equal(d$industry_sizes, c(2L, 1L))
  expect_equal(d$posterior_unit_cells, 12)
  expect_equal(d$materialized_kernel_cells, 12)
  expect_equal(d$same_industry_ordered_pairs, 2)
  expect_equal(d$cross_industry_ordered_pairs, 4)
  expect_equal(d$same_industry_pair_cells, 16)
  expect_equal(d$cross_industry_pair_cells, 64)
  expect_equal(d$naive_industry_full_grid_states, 12)
  expect_equal(d$forbidden_full_xi_grid_cells, 12)
})

test_that("quadrature g_theta artifact equals direct two-support pushforward", {
  b <- .tiny_case_q("race")
  q <- .quad(input = b$input, prior = b$prior, fit = b$fit,
             supp_pts_theta = 31L)
  tl <- q$posterior$metadata$two_level
  direct <- .push(
    supp_xi = tl$support_xi,
    g_xi = tl$g_xi,
    supp_eta = tl$support_eta,
    g_eta = tl$g_eta,
    s = b$input$s,
    mu = tl$mu,
    beta = tl$beta,
    characteristic = tl$characteristic,
    supp_pts_theta = 31L
  )

  expect_equal(q$g_theta$support, direct$support, tolerance = 0)
  expect_equal(q$g_theta$g, direct$g, tolerance = 0)
  expect_equal(q$g_theta$density, direct$density, tolerance = 0)
  expect_identical(q$artifacts$g_theta, q$g_theta)
})

test_that("quadrature can skip g_theta artifact explicitly", {
  b <- .tiny_case_q("race")
  q <- .quad(input = b$input, prior = b$prior, fit = b$fit,
             include_g_theta = FALSE)
  expect_null(q$g_theta)
  expect_null(q$artifacts$g_theta)
  expect_silent(.valid_quad(q))
})

test_that("quadrature refuses incomplete inputs and malformed objects", {
  b <- .tiny_case_q("race")
  expect_error(.quad(input = b$input, prior = b$prior),
               class = "gradepath_validation_error")

  q <- .quad(input = b$input, prior = b$prior, fit = b$fit)
  bad <- q
  bad$diagnostics$materializes_full_grid <- TRUE
  expect_error(.valid_quad(bad), class = "gradepath_validation_error")

  bad2 <- q
  bad2$pairwise_theta <- q$pairwise_bar
  expect_error(.valid_quad(bad2), class = "gradepath_validation_error")

  bad3 <- q
  bad3$artifacts$Pi_xi <- bad3$artifacts$Pi_xi + 0.01
  expect_error(.valid_quad(bad3), class = "gradepath_validation_error")

  bad4 <- q
  bad4$posterior$metadata$two_level$original_s <- 1
  expect_error(.quad(posterior = bad4$posterior),
               class = "gradepath_validation_error")
})

test_that("quadrature scales on a large same-industry smoke without xi^n enumeration", {
  N <- 14L
  b <- .tiny_case_q("race")
  input <- list(
    theta_hat = rep(b$input$theta_hat[1], N) + seq_len(N) * 0.001,
    s = rep(b$input$s[1], N) + seq_len(N) * 0.002,
    industry = rep(1L, N),
    unit_id = paste0("big_", seq_len(N)),
    label = paste0("big ", seq_len(N))
  )
  fit <- list(
    characteristic = "race",
    beta = b$fit$beta,
    mu = 0,
    v_hat = rep(b$fit$v_hat[1], N) + seq_len(N) * 0.001,
    s_v = rep(b$fit$s_v[1], N),
    industry = input$industry
  )
  prior <- b$prior
  q <- .quad(input = input, prior = prior, fit = fit,
             include_g_theta = FALSE)

  expect_silent(.valid_quad(q))
  expect_false(q$diagnostics$materializes_full_grid)
  expect_equal(q$diagnostics$n_units, N)
  expect_equal(q$diagnostics$n_industries, 1L)
  expect_equal(q$diagnostics$materialized_kernel_cells, N * 2L * 2L)
  expect_equal(q$diagnostics$forbidden_full_xi_grid_cells, 2 * (2^N))
  expect_gt(q$diagnostics$forbidden_full_xi_grid_cells,
            q$diagnostics$materialized_kernel_cells)
})
