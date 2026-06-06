# =============================================================================
# test-same-industry-override.R  --  same-industry override + 5 Pi
# -----------------------------------------------------------------------------
# Pins the group_fx == 1 pairwise contract from get_posteriors.m:188-199:
# same-industry pairs count the shared industry likelihood once; cross-industry
# pairs factorize across independent eta effects.  Tiny supports keep the oracle
# fully explicit.
# =============================================================================

.tl_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}

.post2 <- function(...) .tl_get("gp_posterior_twolevel")(...)
.tlpi <- function(...) .tl_get("gp_twolevel_pi")(...)
.valid_tlpi <- function(...) .tl_get("validate_gp_twolevel_pi")(...)
.valid_pair <- function(...) .tl_get("validate_gp_pairwise")(...)
.rowlse2 <- function(...) .tl_get(".gp_row_logsumexp")(...)

.tiny_prior_82 <- function() {
  structure(
    list(
      support = c(0.45, 1.25),
      density = c(0.45, 0.55),
      mean = 0.45 * 0.45 + 1.25 * 0.55,
      scale = "r",
      diagnostics = list(
        group_fx = 1L,
        support_eta = c(0.70, 1.45),
        g_eta = c(0.40, 0.60)
      ),
      metadata = list(characteristic = "race", beta = 0.55, mu = 0)
    ),
    class = c("gp_prior", "list")
  )
}

.tiny_case_82 <- function() {
  s <- c(0.80, 1.35, 1.05)
  beta <- 0.55
  v_hat <- c(0.72, 1.34, 1.08)
  s_v <- c(0.30, 0.34, 0.28)
  list(
    input = list(
      theta_hat = v_hat * s^beta,
      s = s,
      industry = c(1, 1, 2),
      unit_id = paste0("u", 1:3),
      label = paste0("unit ", 1:3)
    ),
    prior = .tiny_prior_82(),
    fit = list(
      characteristic = "race",
      beta = beta,
      mu = 0,
      v_hat = v_hat,
      s_v = s_v,
      industry = c(1, 1, 2)
    )
  )
}

.qxi <- function(tl, i, a) {
  w <- exp(tl$log_kernel[[i]][a, ] - tl$log_m[i, a])
  w / sum(w)
}

.theta_vals <- function(tl, i, a) {
  tl$s_beta[i] * tl$support_eta[a] * tl$support_xi
}

.strict_prob <- function(x, wx, y, wy) {
  sum(outer(wx, wy, `*`) * outer(x, y, `>`))
}

.sq_gap <- function(x, wx, y, wy) {
  gap <- outer(x, y, `-`)
  sum(outer(wx, wy, `*`) * (pmax(gap, 0)^2))
}

.same_manual <- function(post, i, j) {
  tl <- post$metadata$two_level
  k <- tl$industry[i]
  out <- c(Pi_theta = 0, Pi_sq_theta = 0, Pi_xi = 0,
           Pi_psi = 0, Pi_sbar_psi = 0)
  for (a in seq_along(tl$support_eta)) {
    qi <- .qxi(tl, i, a)
    qj <- .qxi(tl, j, a)
    ti <- .theta_vals(tl, i, a)
    tj <- .theta_vals(tl, j, a)
    ea <- tl$eta_posterior[k, a]
    out["Pi_theta"] <- out["Pi_theta"] + ea * .strict_prob(ti, qi, tj, qj)
    out["Pi_sq_theta"] <- out["Pi_sq_theta"] + ea * .sq_gap(ti, qi, tj, qj)
    out["Pi_xi"] <- out["Pi_xi"] + ea * .strict_prob(tl$support_xi, qi,
                                                     tl$support_xi, qj)
  }
  out
}

.cross_manual <- function(post, i, j) {
  tl <- post$metadata$two_level
  dist <- function(unit) {
    k <- tl$industry[unit]
    w <- x_theta <- x_xi <- numeric(0)
    for (a in seq_along(tl$support_eta)) {
      q <- .qxi(tl, unit, a)
      wa <- tl$eta_posterior[k, a] * q
      w <- c(w, wa)
      x_theta <- c(x_theta, .theta_vals(tl, unit, a))
      x_xi <- c(x_xi, tl$support_xi)
    }
    list(theta = x_theta, xi = x_xi, w = w / sum(w),
         eta = tl$support_eta,
         w_eta = tl$eta_posterior[k, ] / sum(tl$eta_posterior[k, ]),
         sbar_eta = tl$sbar[unit] * tl$support_eta)
  }
  di <- dist(i); dj <- dist(j)
  c(
    Pi_theta = .strict_prob(di$theta, di$w, dj$theta, dj$w),
    Pi_sq_theta = .sq_gap(di$theta, di$w, dj$theta, dj$w),
    Pi_xi = .strict_prob(di$xi, di$w, dj$xi, dj$w),
    Pi_psi = .strict_prob(di$eta, di$w_eta, dj$eta, dj$w_eta),
    Pi_sbar_psi = .strict_prob(di$sbar_eta, di$w_eta,
                               dj$sbar_eta, dj$w_eta)
  )
}

.same_squared_wrong <- function(b, i = 1L, j = 2L) {
  p <- b$prior
  f <- b$fit
  supp_xi <- p$support
  g_xi <- p$density / sum(p$density)
  supp_eta <- p$diagnostics$support_eta
  g_eta <- p$diagnostics$g_eta / sum(p$diagnostics$g_eta)
  num <- den <- 0
  for (a in seq_along(supp_eta)) {
    eta <- supp_eta[a]
    for (bi in seq_along(supp_xi)) {
      for (bj in seq_along(supp_xi)) {
        li <- stats::dnorm(f$v_hat[i], mean = eta * supp_xi[bi], sd = f$s_v[i])
        lj <- stats::dnorm(f$v_hat[j], mean = eta * supp_xi[bj], sd = f$s_v[j])
        w_wrong <- g_eta[a] * g_xi[bi] * g_xi[bj] * (li * lj)^2
        ti <- b$input$s[i]^f$beta * eta * supp_xi[bi]
        tj <- b$input$s[j]^f$beta * eta * supp_xi[bj]
        num <- num + (ti > tj) * w_wrong
        den <- den + w_wrong
      }
    }
  }
  num / den
}

test_that("same-industry override uses one shared eta likelihood, not squared evidence", {
  b <- .tiny_case_82()
  post <- .post2(b$input, b$prior, b$fit)
  pi <- .tlpi(post)
  man <- .same_manual(post, 1, 2)

  expect_equal(pi$raw$Pi_theta[1, 2], unname(man["Pi_theta"]), tolerance = 1e-12)
  expect_equal(pi$raw$Pi_sq_theta[1, 2], unname(man["Pi_sq_theta"]), tolerance = 1e-12)
  expect_equal(pi$raw$Pi_xi[1, 2], unname(man["Pi_xi"]), tolerance = 1e-12)
  expect_equal(pi$raw$Pi_psi[1, 2], 0, tolerance = 0)
  expect_equal(pi$raw$Pi_sbar_psi[1, 2], 0, tolerance = 0)

  wrong <- .same_squared_wrong(b, 1, 2)
  expect_gt(abs(pi$raw$Pi_theta[1, 2] - wrong), 1e-3)
})

test_that("cross-industry pairs factorize across independent eta effects", {
  b <- .tiny_case_82()
  post <- .post2(b$input, b$prior, b$fit)
  pi <- .tlpi(post)
  man <- .cross_manual(post, 1, 3)

  for (nm in names(man)) {
    expect_equal(pi$raw[[nm]][1, 3], unname(man[nm]), tolerance = 1e-12)
  }
})

test_that("five raw matrices and cleaned pairwise wrappers have expected shape", {
  b <- .tiny_case_82()
  pi <- .tlpi(.post2(b$input, b$prior, b$fit))

  expect_s3_class(pi, "gp_twolevel_pi")
  expect_silent(.valid_tlpi(pi))
  for (nm in c("Pi_theta", "Pi_sbar_psi", "Pi_sq_theta", "Pi_xi", "Pi_psi")) {
    expect_equal(dim(pi$raw[[nm]]), c(3L, 3L))
    expect_true(all(diag(pi$raw[[nm]]) == 0))
    expect_true(all(is.finite(pi$raw[[nm]])))
  }
  expect_true(all(pi$raw$Pi_sq_theta >= 0))

  expect_equal(dim(pi$Pi_theta), c(3L, 3L))
  expect_equal(dim(pi$Pi_sbar_psi), c(3L, 3L))
  expect_equal(dim(pi$Pi_sbar_psi_industry), c(2L, 2L))
  expect_equal(pi$Pi_bar, pi$Pi_sbar_psi, tolerance = 0)
  expect_equal(pi$Pi_bar_industry, pi$Pi_sbar_psi_industry, tolerance = 0)
  expect_silent(.valid_pair(pi$pairwise_theta))
  expect_silent(.valid_pair(pi$pairwise_bar))
  expect_true(all(diag(pi$pairwise_theta$matrix) == 0.5))
  expect_true(all(diag(pi$pairwise_bar$matrix) == 0.5))
  expect_gte(min(pi$pairwise_theta$matrix[row(pi$pairwise_theta$matrix) !=
                                           col(pi$pairwise_theta$matrix)]), 1e-7)
})

test_that("Pi_bar uses stable numeric industry order and representatives", {
  b <- .tiny_case_82()
  b$input$industry <- c(10, 10, 2)
  b$fit$industry <- b$input$industry
  post <- .post2(b$input, b$prior, b$fit)
  pi <- .tlpi(post)

  expect_identical(pi$industry_levels, c("2", "10"))
  expect_equal(pi$industry_representatives, c(3L, 1L))
  expected <- pi$raw$Pi_sbar_psi[pi$industry_representatives,
                                 pi$industry_representatives]
  dimnames(expected) <- dimnames(pi$metadata$raw_Pi_sbar_psi_industry)
  expect_equal(pi$metadata$raw_Pi_sbar_psi_industry, expected, tolerance = 0)
})

test_that("strict > ties contribute zero", {
  expect_equal(.tl_get(".gp_tl_weighted_pair")(c(1, 1), c(0.4, 0.6),
                                               c(1, 1), c(0.2, 0.8)),
               0, tolerance = 0)
})

test_that("malformed metadata and beta-near-one reconstruction fail clearly", {
  b <- .tiny_case_82()
  post <- .post2(b$input, b$prior, b$fit)
  bad <- post
  bad$metadata$two_level$log_m <- bad$metadata$two_level$log_m[, 1, drop = FALSE]
  expect_error(.tlpi(bad), class = "gradepath_validation_error")

  bad_eta <- post
  bad_eta$metadata$two_level$eta_posterior[1, ] <- c(0.2, 0.2)
  expect_error(.tlpi(bad_eta), class = "gradepath_validation_error")

  bad_sbar <- post
  bad_sbar$metadata$two_level$sbar[1] <- bad_sbar$metadata$two_level$sbar[1] + 0.1
  expect_error(.tlpi(bad_sbar), class = "gradepath_validation_error")

  old <- post
  old$metadata$two_level$s_beta <- NULL
  old$metadata$two_level$beta <- 1
  expect_error(.tlpi(old), class = "gradepath_validation_error")
})
