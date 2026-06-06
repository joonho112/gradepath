# =============================================================================
# test-posterior-twolevel.R  --  two-level posterior summaries
# -----------------------------------------------------------------------------
# Tiny synthetic parity checks for the deterministic factorization in
# gp_posterior_twolevel().  The oracle below literally enumerates eta x xi^n
# states inside each industry to check the group likelihood, then collapses each
# unit's duplicated marginal atoms to the eta x xi_i support consumed by the
# deterministic quadrature.
# =============================================================================

.gp2post_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}

.post2 <- function(...) .gp2post_get("gp_posterior_twolevel")(...)
.valid_post <- function(...) .gp2post_get("validate_gp_posterior")(...)
.wpr <- function(...) .gp2post_get(".gp_weighted_percentile")(...)
.rowlse <- function(...) .gp2post_get(".gp_row_logsumexp")(...)

.wpr_literal <- function(x, p, w = NULL) {
  x <- as.numeric(x)
  if (is.null(w)) w <- rep(1, length(x))
  w <- as.numeric(w)
  ord <- order(x)
  xs <- x[ord]
  ws <- w[ord]
  pk <- (cumsum(ws) - 0.5 * ws) / sum(ws)
  stats::approx(c(0, pk, 1), c(xs[1], xs, xs[length(xs)]),
                xout = as.numeric(p) / 100, ties = "ordered", rule = 2)$y
}

.collapse_atoms <- function(x, w) {
  key <- sprintf("%.17g", x)
  keep <- !duplicated(key)
  list(x = x[keep], w = as.numeric(rowsum(w, group = key, reorder = FALSE)[, 1L]))
}

.tiny_prior <- function(characteristic = "race") {
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (identical(characteristic, "race")) {
    supp_xi <- c(0.65, 1.15)
    g_xi <- c(0.35, 0.65)
    supp_eta <- c(0.75, 1.30)
    g_eta <- c(0.55, 0.45)
  } else {
    supp_xi <- c(-0.55, 0.35)
    g_xi <- c(0.40, 0.60)
    supp_eta <- c(-0.25, 0.45)
    g_eta <- c(0.30, 0.70)
  }
  structure(
    list(
      support = supp_xi,
      density = g_xi,
      mean = sum(supp_xi * g_xi),
      scale = "r",
      diagnostics = list(group_fx = 1L, support_eta = supp_eta, g_eta = g_eta),
      metadata = list(characteristic = characteristic, beta = 0.5)
    ),
    class = c("gp_prior", "list")
  )
}

.tiny_case <- function(characteristic = "race") {
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (identical(characteristic, "race")) {
    s <- c(0.80, 1.40, 1.05)
    beta <- 0.55
    mu <- 0
    v_hat <- c(0.92, 1.18, 0.82)
    s_v <- c(0.31, 0.36, 0.29)
    theta_hat <- v_hat * s^beta
  } else {
    s <- c(0.90, 1.30, 1.10)
    beta <- 0.65
    mu <- 0.18
    v_hat <- c(-0.15, 0.32, 0.08)
    s_v <- c(0.27, 0.33, 0.30)
    theta_hat <- mu + v_hat * s^beta
  }
  input <- list(
    theta_hat = theta_hat,
    s = s,
    industry = c("A", "A", "B"),
    unit_id = paste0(characteristic, "_", seq_along(s)),
    label = paste(characteristic, seq_along(s))
  )
  fit <- list(
    characteristic = characteristic,
    beta = beta,
    mu = mu,
    v_hat = v_hat,
    s_v = s_v,
    industry = input$industry
  )
  prior <- .tiny_prior(characteristic)
  prior$metadata$beta <- beta
  prior$metadata$mu <- mu
  list(input = input, prior = prior, fit = fit)
}

.oracle_twolevel <- function(input, prior, fit, interval_level = 0.90) {
  supp_xi <- prior$support
  g_xi <- prior$density / sum(prior$density)
  supp_eta <- prior$diagnostics$support_eta
  g_eta <- prior$diagnostics$g_eta / sum(prior$diagnostics$g_eta)
  industry <- match(as.character(fit$industry),
                    sort(unique(as.character(fit$industry)), method = "radix"))
  N <- length(fit$v_hat)
  pct <- c((1 - interval_level) / 2, (1 + interval_level) / 2) * 100
  vals_r <- vals_t <- weights <- vector("list", N)
  for (i in seq_len(N)) {
    vals_r[[i]] <- vals_t[[i]] <- weights[[i]] <- numeric(0)
  }

  s_beta <- input$s^fit$beta
  for (k in sort(unique(industry))) {
    idx <- which(industry == k)
    grids <- expand.grid(rep(list(seq_along(supp_xi)), length(idx)),
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    for (a in seq_along(supp_eta)) {
      eta <- supp_eta[a]
      for (row in seq_len(nrow(grids))) {
        xi_idx <- as.integer(grids[row, ])
        xi <- supp_xi[xi_idx]
        w <- g_eta[a] * prod(g_xi[xi_idx])
        for (h in seq_along(idx)) {
          i <- idx[h]
          mean_i <- if (identical(fit$characteristic, "race")) eta * xi[h]
                    else eta + xi[h]
          w <- w * stats::dnorm(fit$v_hat[i], mean = mean_i, sd = fit$s_v[i])
        }
        for (h in seq_along(idx)) {
          i <- idx[h]
          r <- if (identical(fit$characteristic, "race")) eta * xi[h]
               else eta + xi[h]
          theta <- if (identical(fit$characteristic, "race")) s_beta[i] * r
                   else fit$mu + s_beta[i] * r
          vals_r[[i]] <- c(vals_r[[i]], r)
          vals_t[[i]] <- c(vals_t[[i]], theta)
          weights[[i]] <- c(weights[[i]], w)
        }
      }
    }
  }

  out <- data.frame(
    E_r = numeric(N), E_r2 = numeric(N), lo_r = numeric(N), up_r = numeric(N),
    E_theta = numeric(N), E_theta2 = numeric(N),
    lo_theta = numeric(N), up_theta = numeric(N)
  )
  for (i in seq_len(N)) {
    w <- weights[[i]] / sum(weights[[i]])
    out$E_r[i] <- sum(w * vals_r[[i]])
    out$E_r2[i] <- sum(w * vals_r[[i]]^2)
    atoms_r <- .collapse_atoms(vals_r[[i]], w)
    qr <- .wpr_literal(atoms_r$x, pct, atoms_r$w)
    out$lo_r[i] <- min(qr[1], out$E_r[i])
    out$up_r[i] <- max(qr[2], out$E_r[i])
    out$E_theta[i] <- sum(w * vals_t[[i]])
    out$E_theta2[i] <- sum(w * vals_t[[i]]^2)
    atoms_t <- .collapse_atoms(vals_t[[i]], w)
    qt <- .wpr_literal(atoms_t$x, pct, atoms_t$w)
    out$lo_theta[i] <- qt[1]
    out$up_theta[i] <- qt[2]
  }
  out
}

test_that("weighted percentile matches Matlab wprctile type-5 convention", {
  x <- c(4, 1, 8, 2, 6)
  p <- c(0, 25, 50, 75, 100)
  expect_equal(.wpr(x, p), as.numeric(stats::quantile(x, p / 100, type = 5)),
               tolerance = 1e-12)

  w <- c(2, 0, 1, 3, 4)
  expect_equal(.wpr(x, c(5, 50, 95), w),
               .wpr_literal(x, c(5, 50, 95), w),
               tolerance = 1e-12)
  expect_error(.wpr(x, 50, c(0, 0, 0, 0, 0)),
               class = "gradepath_validation_error")
  expect_error(.wpr(x, 50, c(1, -1, 1, 1, 1)),
               class = "gradepath_validation_error")
})

test_that("race posterior factorization matches literal eta x xi^n enumeration", {
  b <- .tiny_case("race")
  post <- .post2(b$input, b$prior, b$fit, interval_level = 0.90)
  oracle <- .oracle_twolevel(b$input, b$prior, b$fit, interval_level = 0.90)

  expect_s3_class(post, "gp_posterior")
  expect_silent(.valid_post(post))
  expect_identical(post$scale, "r")
  expect_equal(post$posterior_mean, oracle$E_r, tolerance = 1e-12)
  expect_equal(post$posterior_sd, sqrt(pmax(oracle$E_r2 - oracle$E_r^2, 0)),
               tolerance = 1e-12)
  expect_equal(post$lower, oracle$lo_r, tolerance = 1e-12)
  expect_equal(post$upper, oracle$up_r, tolerance = 1e-12)

  rep <- post$metadata$reporting
  expect_identical(rep$scale, "theta")
  expect_equal(rep$posterior_mean, oracle$E_theta, tolerance = 1e-12)
  expect_equal(rep$posterior_second_moment, oracle$E_theta2, tolerance = 1e-12)
  expect_equal(rep$lower, oracle$lo_theta, tolerance = 1e-12)
  expect_equal(rep$upper, oracle$up_theta, tolerance = 1e-12)
  expect_equal(rep$posteriors$posterior_mean, rep$posterior_mean, tolerance = 0)
  expect_equal(rep$posteriors$posterior_second_moment,
               rep$posterior_second_moment, tolerance = 0)
  expect_identical(rep$percentile_convention,
                   "collapsed_marginal_wprctile_type5")
})

test_that("gender posterior uses the additive mu + s^beta * (eta + xi) transform", {
  b <- .tiny_case("gender")
  post <- .post2(b$input, b$prior, b$fit, interval_level = 0.95)
  oracle <- .oracle_twolevel(b$input, b$prior, b$fit, interval_level = 0.95)

  expect_silent(.valid_post(post))
  expect_equal(post$posterior_mean, oracle$E_r, tolerance = 1e-12)
  expect_equal(post$metadata$reporting$posterior_mean,
               oracle$E_theta, tolerance = 1e-12)
  expect_equal(post$metadata$reporting$posterior_second_moment,
               oracle$E_theta2, tolerance = 1e-12)
  expect_equal(post$metadata$reporting$lower,
               oracle$lo_theta, tolerance = 1e-12)
  expect_equal(post$metadata$reporting$upper,
               oracle$up_theta, tolerance = 1e-12)
})

test_that("interval level controls weighted theta percentiles", {
  b <- .tiny_case("race")
  post90 <- .post2(b$input, b$prior, b$fit)
  post95 <- .post2(b$input, b$prior, b$fit, interval_level = 0.95)
  oracle95 <- .oracle_twolevel(b$input, b$prior, b$fit, interval_level = 0.95)

  expect_equal(post90$metadata$interval_level, 0.90, tolerance = 0)
  expect_equal(post95$metadata$interval_level, 0.95, tolerance = 0)
  expect_equal(post95$metadata$reporting$lower, oracle95$lo_theta,
               tolerance = 1e-12)
  expect_equal(post95$metadata$reporting$upper, oracle95$up_theta,
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(post90$metadata$reporting$lower,
                                post95$metadata$reporting$lower)))
})

test_that("theta summaries use original input$s while r-scale likelihood carriers stay fixed", {
  b <- .tiny_case("race")
  post_a <- .post2(b$input, b$prior, b$fit)
  b$input$s <- b$input$s * c(1.4, 0.8, 1.7)
  post_b <- .post2(b$input, b$prior, b$fit)

  expect_equal(post_a$posterior_mean, post_b$posterior_mean, tolerance = 0)
  expect_equal(post_a$posterior_sd, post_b$posterior_sd, tolerance = 0)
  expect_false(isTRUE(all.equal(post_a$metadata$reporting$posterior_mean,
                                post_b$metadata$reporting$posterior_mean)))
})

test_that("two-level metadata keeps likelihood blocks", {
  b <- .tiny_case("race")
  post <- .post2(b$input, b$prior, b$fit)
  tl <- post$metadata$two_level
  expect_equal(dim(tl$log_m), c(3L, 2L))
  expect_equal(dim(tl$log_z), c(2L, 2L))
  expect_equal(dim(tl$eta_posterior), c(2L, 2L))
  expect_length(tl$log_kernel, 3L)
  expect_true(all(abs(rowSums(tl$eta_posterior) - 1) < 1e-12))
  expect_equal(tl$posterior_components$E_theta,
               post$metadata$reporting$posterior_mean,
               tolerance = 0)
  expect_identical(tl$percentile_convention,
                   "collapsed_marginal_wprctile_type5")
})

test_that("metadata reconstructs same-industry and cross-industry eta weighting", {
  b <- .tiny_case("race")
  post <- .post2(b$input, b$prior, b$fit)
  tl <- post$metadata$two_level

  # Same industry: units 1 and 2 share the single industry likelihood, so the
  # pairwise eta marginal is the one industry posterior eta distribution.
  qi_mass <- exp(.rowlse(tl$log_kernel[[1]]) - tl$log_m[1, ])
  qj_mass <- exp(.rowlse(tl$log_kernel[[2]]) - tl$log_m[2, ])
  same_eta_marginal <- exp(tl$log_z[1, ] - tl$log_denominator[1]) *
    qi_mass * qj_mass
  expect_equal(qi_mass, rep(1, length(qi_mass)), tolerance = 1e-12)
  expect_equal(qj_mass, rep(1, length(qj_mass)), tolerance = 1e-12)
  expect_equal(same_eta_marginal, tl$eta_posterior[1, ],
               tolerance = 1e-12)

  # Cross industry: unit 1 and unit 3 have independent industry effects, so the
  # eta joint marginal is the outer product of the two industry eta posteriors.
  qh_mass <- exp(.rowlse(tl$log_kernel[[3]]) - tl$log_m[3, ])
  cross_eta_joint <- outer(
    exp(tl$log_z[1, ] - tl$log_denominator[1]) * qi_mass,
    exp(tl$log_z[2, ] - tl$log_denominator[2]) * qh_mass,
    `*`
  )
  expect_equal(cross_eta_joint,
               outer(tl$eta_posterior[1, ], tl$eta_posterior[2, ], `*`),
               tolerance = 1e-12)
})

test_that("numeric industry labels preserve numeric ordering in metadata", {
  b <- .tiny_case("race")
  b$input$industry <- c(1, 10, 2)
  b$fit$industry <- b$input$industry
  post <- .post2(b$input, b$prior, b$fit)

  expect_identical(post$metadata$two_level$industry_levels, c("1", "2", "10"))
  expect_equal(post$metadata$two_level$industry, c(1L, 3L, 2L))
})
