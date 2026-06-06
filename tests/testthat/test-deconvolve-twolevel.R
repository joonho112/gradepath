# =============================================================================
# test-deconvolve-twolevel.R  --  native two-level deconvolution tests
# -----------------------------------------------------------------------------
# Structural and algebraic tests for the native group_fx == 1 deconvolution. The
# expensive real-data 72/108-node penalty grid is slow-gated; default tests use
# tiny synthetic carriers and fixed penalties.
# =============================================================================

skip_if_no_2l_decon <- function() {
  skip_if_not_installed("ebrecipe")
  skip_if_not(exists("gp_deconvolve_groups", mode = "function"),
              "two-level deconvolution source (gp_deconvolve_groups) not loaded")
}

.gp_2l_test_fit <- function(characteristic = "race") {
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (identical(characteristic, "race")) {
    v_hat <- c(0.8, 1.2, 0.6, 1.4)
    s_v <- c(0.18, 0.20, 0.16, 0.22)
    industry <- c(1L, 1L, 2L, 2L)
    list(
      characteristic = "race",
      m_hat = c(1.0, 0.30, 0.20),
      V_m = diag(c(0.20, 0.10, 0.08)),
      v_hat = v_hat,
      s_v = s_v,
      industry = industry,
      beta = 0.5
    )
  } else {
    v_hat <- c(-0.20, 0.10, 0.40, -0.30)
    s_v <- c(0.18, 0.20, 0.16, 0.22)
    industry <- c(1L, 1L, 2L, 2L)
    list(
      characteristic = "gender",
      m_hat = c(0.25, 0.15),
      V_m = diag(c(0.10, 0.08)),
      v_hat = v_hat,
      s_v = s_v,
      industry = industry,
      beta = 1.1
    )
  }
}

.gp_2l_test_control <- function() {
  ctl <- gp_control()
  ctl$deconv2l_supp_pts <- 11L
  ctl$deconv2l_n_knots <- 3L
  ctl$deconv2l_max_iter <- 60L
  ctl$tol <- 1e-6
  ctl
}

.gp_2l_literal_loglik <- function(r, s_v, industry, supp_xi, g_xi,
                                  supp_eta, g_eta, characteristic) {
  out <- 0
  for (k in sort(unique(industry))) {
    idx <- industry == k
    y <- r[idx]
    s <- s_v[idx]
    p_hat <- numeric(length(supp_eta))
    for (l in seq_along(supp_eta)) {
      psi <- supp_eta[l]
      P <- matrix(0, nrow = length(y), ncol = length(supp_xi))
      for (i in seq_along(y)) {
        mean_i <- if (identical(characteristic, "race")) {
          psi * supp_xi
        } else {
          psi + supp_xi
        }
        P[i, ] <- stats::dnorm(y[i], mean = mean_i, sd = s[i])
      }
      p_hat[l] <- prod(P %*% g_xi)
    }
    out <- out + log(sum(p_hat * g_eta))
  }
  out
}

test_that("two-level empirical components and support caps match Matlab formulas", {
  skip_if_no_2l_decon()

  fit_r <- .gp_2l_test_fit("race")
  comp_r <- .gp_2l_components(fit_r$v_hat, fit_r$industry, "race", mu = fit_r$m_hat[1])
  expect_equal(comp_r$psi_hat, rep(c(1, 1), each = 2), tolerance = 1e-12)
  expect_equal(comp_r$xi_hat, fit_r$v_hat, tolerance = 1e-12)
  caps_r <- .gp_2l_support_caps(fit_r$m_hat, comp_r$psi_hat, comp_r$xi_hat, "race")
  expect_equal(caps_r$xi, c(lo = 0, hi = 2.5), tolerance = 1e-12)
  expect_equal(caps_r$eta, c(lo = 0, hi = 2.0), tolerance = 1e-12)

  fit_g <- .gp_2l_test_fit("gender")
  comp_g <- .gp_2l_components(fit_g$v_hat, fit_g$industry, "gender")
  expect_equal(comp_g$psi_hat, rep(c(-0.05, 0.05), each = 2), tolerance = 1e-12)
  expect_equal(comp_g$xi_hat, fit_g$v_hat - comp_g$psi_hat, tolerance = 1e-12)
  caps_g <- .gp_2l_support_caps(fit_g$m_hat, comp_g$psi_hat, comp_g$xi_hat, "gender")
  expect_equal(caps_g$xi, c(lo = -1.25, hi = 1.25), tolerance = 1e-12)
  expect_equal(caps_g$eta, c(lo = -0.75, hi = 0.75), tolerance = 1e-12)
})

test_that("default two-level penalty grids match the Matlab combvec shape", {
  skip_if_no_2l_decon()
  gr <- .gp_2l_penalty_grid("race")
  gg <- .gp_2l_penalty_grid("gender")
  expect_equal(nrow(gr), 72L)
  expect_equal(nrow(gg), 108L)
  expect_equal(unname(unlist(gr[1, ])), c(0.08, 0.0025), tolerance = 1e-12)
  expect_equal(gr[2, "c_xi"], 0.085, tolerance = 1e-12)
  expect_equal(unname(unlist(gr[10, ])), c(0.08, 0.005), tolerance = 1e-12)
  expect_equal(unname(unlist(gg[1, ])), c(0.01, 0.0025), tolerance = 1e-12)
})

test_that("log-space likelihood matches the literal Matlab loop on tiny race/gender problems", {
  skip_if_no_2l_decon()
  for (ch in c("race", "gender")) {
    fit <- .gp_2l_test_fit(ch)
    comp <- .gp_2l_components(fit$v_hat, fit$industry, ch,
                              mu = if (ch == "race") fit$m_hat[1] else NULL)
    caps <- .gp_2l_support_caps(fit$m_hat, comp$psi_hat, comp$xi_hat, ch)
    supp_xi <- .gp_2l_grid(caps$xi, supp_pts = 9L)
    supp_eta <- .gp_2l_grid(caps$eta, supp_pts = 9L)
    Q_xi <- .gp_eb_spline_basis(supp_xi, n_knots = 3L)
    Q_eta <- .gp_eb_spline_basis(supp_eta, n_knots = 3L)
    ax <- rep(0, ncol(Q_xi) - if (ch == "gender") 1L else 0L)
    ae <- rep(0, ncol(Q_eta) - 1L)

    out <- gp_two_level_likelihood(
      ax, ae, fit$v_hat, fit$s_v, fit$industry,
      Q_xi, supp_xi, Q_eta, supp_eta,
      c_xi = 0.01, c_eta = 0.01, characteristic = ch,
      return_densities = TRUE
    )
    literal <- .gp_2l_literal_loglik(
      fit$v_hat, fit$s_v, fit$industry, supp_xi, out$g_xi,
      supp_eta, out$g_eta, ch
    )
    expect_equal(out$loglik, literal, tolerance = 1e-10)
  }
})

test_that("mean pins land on eta for race and xi/eta for gender", {
  skip_if_no_2l_decon()
  for (ch in c("race", "gender")) {
    fit <- .gp_2l_test_fit(ch)
    comp <- .gp_2l_components(fit$v_hat, fit$industry, ch,
                              mu = if (ch == "race") fit$m_hat[1] else NULL)
    caps <- .gp_2l_support_caps(fit$m_hat, comp$psi_hat, comp$xi_hat, ch)
    supp_xi <- .gp_2l_grid(caps$xi, supp_pts = 11L)
    supp_eta <- .gp_2l_grid(caps$eta, supp_pts = 11L)
    Q_xi <- .gp_eb_spline_basis(supp_xi, n_knots = 3L)
    Q_eta <- .gp_eb_spline_basis(supp_eta, n_knots = 3L)
    ax <- rep(0, ncol(Q_xi) - if (ch == "gender") 1L else 0L)
    ae <- rep(0, ncol(Q_eta) - 1L)
    out <- gp_two_level_likelihood(
      ax, ae, fit$v_hat, fit$s_v, fit$industry,
      Q_xi, supp_xi, Q_eta, supp_eta,
      c_xi = 0.01, c_eta = 0.01, characteristic = ch,
      return_densities = TRUE
    )
    expect_equal(sum(supp_eta * out$g_eta), if (ch == "race") 1 else 0,
                 tolerance = 1e-8)
    if (identical(ch, "gender")) {
      expect_equal(sum(supp_xi * out$g_xi), 0, tolerance = 1e-8)
    }
  }
})

test_that("gp_deconvolve_groups returns a valid gp_prior-shaped two-level object", {
  skip_if_no_2l_decon()
  fit <- .gp_2l_test_fit("race")
  prior <- suppressWarnings(gp_deconvolve_groups(
    fit,
    control = .gp_2l_test_control(),
    supp_pts = 11L,
    penalties = list(c_xi = 0.01, c_eta = 0.01)
  ))
  expect_s3_class(prior, "gp_prior")
  expect_silent(validate_gp_prior(prior))
  expect_identical(prior$scale, "r")
  expect_equal(sum(prior$density), 1, tolerance = 1e-8)
  expect_equal(sum(prior$diagnostics$g_eta), 1, tolerance = 1e-8)
  expect_identical(prior$diagnostics$group_fx, 1L)
  expect_identical(prior$diagnostics$grid_is_fixed, TRUE)
  expect_equal(prior$diagnostics$c_xi, 0.01)
  expect_equal(prior$diagnostics$c_eta, 0.01)
  expect_length(prior$diagnostics$model_moments, 3L)
  expect_true(is.finite(prior$diagnostics$J))
  expect_null(prior$diagnostics$g_theta)
})

test_that("two-level multistart preserves the caller RNG state", {
  skip_if_no_2l_decon()
  fit <- .gp_2l_test_fit("race")
  set.seed(20260702)
  seed_before <- .Random.seed
  prior <- suppressWarnings(gp_deconvolve_groups(
    fit,
    control = .gp_2l_test_control(),
    n_starts = 2L,
    supp_pts = 11L,
    penalties = list(c_xi = 0.01, c_eta = 0.01)
  ))
  expect_s3_class(prior, "gp_prior")
  expect_identical(.Random.seed, seed_before)
})

test_that("V_m inverse and industry-order guards fail loudly", {
  skip_if_no_2l_decon()
  bad <- .gp_2l_test_fit("race")
  bad$V_m <- matrix(0, 3, 3)
  expect_error(
    gp_deconvolve_groups(bad, control = .gp_2l_test_control(),
                         penalties = list(c_xi = 0.01, c_eta = 0.01)),
    class = "gradepath_singular_error"
  )

  fit <- .gp_2l_test_fit("race")
  expect_error(
    gp_deconvolve_groups(
      fit,
      input = list(industry = c(1L, 2L, 1L, 2L)),
      control = .gp_2l_test_control(),
      penalties = list(c_xi = 0.01, c_eta = 0.01)
    ),
    class = "gradepath_validation_error"
  )
})

test_that("real-data fixed-penalty two-level solve is available behind the slow gate", {
  skip_if_no_2l_decon()
  skip_if_not(identical(Sys.getenv("GRADEPATH_RUN_SLOW_TESTS"), "1"),
              "set GRADEPATH_RUN_SLOW_TESTS=1 for the real-data fixed-penalty two-level solve")
  fit <- gp_two_level_gmm(gp_krw_gmm_input("race"), characteristic = "race",
                          n_starts = 5L)
  ctl <- gp_control()
  ctl$deconv2l_supp_pts <- 25L
  ctl$deconv2l_n_knots <- 5L
  ctl$deconv2l_max_iter <- 80L
  prior <- suppressWarnings(gp_deconvolve_groups(
    fit,
    control = ctl,
    penalties = list(c_xi = 0.105, c_eta = 0.0025)
  ))
  expect_s3_class(prior, "gp_prior")
  expect_silent(validate_gp_prior(prior))
  expect_equal(prior$diagnostics$c_xi, 0.105)
  expect_equal(prior$diagnostics$c_eta, 0.0025)
})
