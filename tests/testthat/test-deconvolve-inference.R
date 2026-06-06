# =============================================================================
# test-deconvolve-inference.R -- deconvolution inference layer
# =============================================================================

.gp_inf_quad <- function(A) {
  function(x) 0.5 * as.numeric(crossprod(x, A %*% x))
}

.gp_inf_small_twolevel <- function() {
  supp_xi <- seq(0.2, 2.0, length.out = 9)
  supp_eta <- seq(0.5, 1.5, length.out = 9)
  Q_xi <- .gp_eb_spline_basis(supp_xi, n_knots = 3)
  Q_eta <- .gp_eb_spline_basis(supp_eta, n_knots = 3)

  alpha_xi <- c(0.10, -0.05, 0.02)
  alpha_eta_free <- c(0.03, -0.02)
  g_xi <- .gp_eb_softmax_density(Q_xi, alpha_xi)$g
  alpha_eta <- .gp_eb_full_alpha(alpha_eta_free, Q_eta, supp_eta,
                                 target_mean = 1)
  g_eta <- .gp_eb_softmax_density(Q_eta, alpha_eta)$g

  prior <- new_gp_prior(
    support = supp_xi,
    density = g_xi,
    mean = sum(supp_xi * g_xi),
    scale = "r",
    diagnostics = list(
      support_eta = supp_eta,
      g_eta = g_eta,
      c_xi = 0.001,
      c_eta = 0.001,
      n_knots = 3
    ),
    metadata = list(
      characteristic = "race",
      alpha_xi_free = alpha_xi,
      alpha_eta_free = alpha_eta_free,
      beta = 0.5
    )
  )
  input <- list(s = c(0.4, 0.6, 0.5, 0.7), industry = c(1, 1, 2, 2))
  fit <- list(
    characteristic = "race",
    v_hat = c(0.9, 1.1, 1.0, 1.2),
    s_v = input$s^(1 - 0.5),
    beta = 0.5,
    mu = 1,
    m_hat = c(1, 0.5, 0.2),
    industry = input$industry
  )
  list(prior = prior, fit = fit, input = input, Q_xi = Q_xi, Q_eta = Q_eta)
}

test_that("gp_hessian matches a quadratic objective and is symmetric", {
  A <- matrix(c(4, 1, 1, 3), 2, 2)
  H <- gp_hessian(.gp_inf_quad(A), c(a = 0.7, b = -0.4), step = 1e-4)
  expect_equal(unname(H), A, tolerance = 1e-6)
  expect_equal(H, t(H), tolerance = 1e-12)
})

test_that("sandwich covariance uses penalized and unpenalized Hessians", {
  H1 <- matrix(c(4, 1, 1, 3), 2, 2)
  H2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  out <- gp_sandwich_vcov(
    par = c(0.2, -0.3),
    penalized_fn = .gp_inf_quad(H1),
    unpenalized_fn = .gp_inf_quad(H2),
    step = 1e-4
  )
  ref <- solve(H1) %*% H2 %*% solve(H1)
  expect_equal(unname(out$H1), H1, tolerance = 1e-6)
  expect_equal(unname(out$H2), H2, tolerance = 1e-6)
  expect_equal(unname(out$vcov), ref, tolerance = 1e-6)
  expect_equal(out$vcov, t(out$vcov), tolerance = 1e-12)

  expect_error(
    gp_sandwich_vcov(
      par = c(1, 1),
      penalized_fn = function(x) 1,
      unpenalized_fn = .gp_inf_quad(diag(2)),
      step = 1e-4
    ),
    class = "gradepath_singular_error"
  )
})

test_that("delta-method SEs use the finite-difference Jacobian", {
  fn <- function(x) c(sum = x[1] + 2 * x[2], diff = x[1] - x[2])
  V <- matrix(c(0.25, 0.05, 0.05, 0.16), 2, 2)
  out <- gp_delta_method_se(c(1, 2), fn, V, step = 1e-5)
  J <- matrix(c(1, 2, 1, -1), nrow = 2, byrow = TRUE)
  expect_equal(out$value, fn(c(1, 2)), tolerance = 1e-12)
  expect_equal(unname(out$jacobian), J, tolerance = 1e-8)
  expect_equal(unname(out$se), sqrt(diag(J %*% V %*% t(J))), tolerance = 1e-8)

  out4 <- gp_delta_method_se(c(1, 2), fn, 4 * V, step = 1e-5)
  expect_equal(out4$se, 2 * out$se, tolerance = 1e-8)
  out0 <- gp_delta_method_se(c(1, 2), fn, 0 * V, step = 1e-5)
  expect_equal(out0$se, c(sum = 0, diff = 0), tolerance = 1e-12)
})

test_that("gp_delta_targets returns KRW xi/eta/theta moment functionals", {
  obj <- .gp_inf_small_twolevel()
  par <- c(obj$prior$metadata$alpha_xi_free,
           obj$prior$metadata$alpha_eta_free)

  f1 <- gp_delta_targets(
    par, obj$Q_xi, obj$prior$support,
    obj$Q_eta, obj$prior$diagnostics$support_eta,
    s = c(0.5, 1.0), mu = 1, beta = 0.5, characteristic = "race",
    supp_pts_theta = 31L
  )
  f2 <- gp_delta_targets(
    par, obj$Q_xi, obj$prior$support,
    obj$Q_eta, obj$prior$diagnostics$support_eta,
    s = c(0.8, 1.6), mu = 1, beta = 0.5, characteristic = "race",
    supp_pts_theta = 31L
  )

  expect_length(f1, 12L)
  expect_named(
    f1,
    c("xi_mean", "xi_sd", "xi_skew", "xi_excess_kurt",
      "eta_mean", "eta_sd", "eta_skew", "eta_excess_kurt",
      "theta_mean", "theta_sd", "theta_skew", "theta_excess_kurt")
  )
  expect_equal(f1[1:8], f2[1:8], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(f1[9:12], f2[9:12])))
})

test_that("generic parametric bootstrap is seeded and preserves caller RNG", {
  set.seed(987)
  before <- .Random.seed
  out1 <- gp_param_bootstrap(
    simulate = function(b) stats::rnorm(2),
    estimate = function(sim, b) sim + b / 10,
    statistic = function(est, b, sim) c(mean = mean(est), first = est[1]),
    B = 6,
    seed = 123
  )
  after <- .Random.seed

  out2 <- gp_param_bootstrap(
    simulate = function(b) stats::rnorm(2),
    estimate = function(sim, b) sim + b / 10,
    statistic = function(est, b, sim) c(mean = mean(est), first = est[1]),
    n_rep = 6,
    seed = 123
  )
  out3 <- gp_param_bootstrap(
    simulate = function(b) stats::rnorm(2),
    estimate = function(sim, b) sim + b / 10,
    statistic = function(est, b, sim) c(mean = mean(est), first = est[1]),
    n_rep = 6,
    seed = 456
  )

  expect_identical(after, before)
  expect_equal(out1$draws, out2$draws, tolerance = 0)
  expect_false(isTRUE(all.equal(out1$draws, out3$draws)))
  expect_equal(out1$se, apply(out1$draws, 2, stats::sd), tolerance = 1e-12)
  expect_named(out1$se, c("mean", "first"))
  expect_error(
    gp_param_bootstrap(
      simulate = function(b) b,
      estimate = function(sim, b) sim,
      statistic = function(est, b, sim) est,
      n_rep = 2,
      B = 3
    ),
    class = "gradepath_validation_error"
  )
})

test_that("two-level parametric bootstrap reuses selected supports and penalties", {
  obj <- .gp_inf_small_twolevel()
  set.seed(654)
  before <- .Random.seed
  out <- gp_param_bootstrap(
    obj$prior,
    obj$fit,
    obj$input,
    B = 2,
    seed = 44,
    max_iter = 1,
    tol = 1e-3,
    supp_pts_theta = 31L
  )
  after <- .Random.seed

  expect_identical(after, before)
  expect_equal(dim(out$draws), c(2L, 12L))
  expect_named(out$se, colnames(out$draws))
  expect_equal(out$se, apply(out$draws, 2, stats::sd), tolerance = 1e-12)
  expect_identical(out$diagnostics$method, "two-level-parametric-bootstrap")
  expect_equal(out$diagnostics$c_xi, obj$prior$diagnostics$c_xi)
  expect_equal(out$diagnostics$c_eta, obj$prior$diagnostics$c_eta)
})

test_that("deconvolution inference wrapper returns a 12-target SE bundle", {
  obj <- .gp_inf_small_twolevel()
  out <- gp_deconvolution_inference(
    obj$prior,
    obj$fit,
    obj$input,
    hessian_step = 1e-3,
    jacobian_step = 1e-4,
    supp_pts_theta = 31L
  )

  expect_s3_class(out, "gp_deconvolution_inference")
  expect_length(out$value, 12L)
  expect_length(out$se, 12L)
  expect_true(all(is.finite(out$se)))
  expect_equal(dim(out$vcov), c(5L, 5L))
  expect_equal(out$vcov, t(out$vcov), tolerance = 1e-10)
  expect_identical(out$diagnostics$method, "native-deconvolution-inference")

  H1 <- gp_hessian(obj$prior, obj$fit, obj$input, penalized = TRUE,
                   step = 1e-3)
  H2 <- gp_hessian(obj$prior, obj$fit, obj$input, penalized = FALSE,
                   step = 1e-3)
  expect_equal(H1, out$H1, tolerance = 1e-10)
  expect_equal(H2, out$H2, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(H1, H2)))

  sw <- gp_deconvolution_sandwich(obj$prior, obj$fit, obj$input, step = 1e-3)
  dm <- gp_delta_method_se(obj$prior, sw$vcov, obj$fit, obj$input,
                           step = 1e-4, supp_pts_theta = 31L)
  expect_equal(sw$vcov, out$vcov, tolerance = 1e-10)
  expect_equal(dm$se, out$se, tolerance = 1e-10)
})

test_that("deconvolution inference wrapper rejects malformed inputs loudly", {
  obj <- .gp_inf_small_twolevel()
  bad_input <- obj$input
  bad_input$s <- bad_input$s[-1]
  expect_error(
    gp_deconvolution_inference(obj$prior, obj$fit, bad_input),
    class = "gradepath_validation_error"
  )
})
