# =============================================================================
# test-pushforward-theta.R  --  native two-level theta pushforward
# =============================================================================

skip_if_no_pushforward <- function() {
  skip_if_not(exists("gp_pushforward_theta", mode = "function"),
              "two-level pushforward source (gp_pushforward_theta) not loaded")
}

.gp_pushfwd_literal <- function(supp_xi, g_xi, supp_eta, g_eta,
                                s, mu, beta, characteristic,
                                supp_pts_theta = 250L) {
  vals_all <- numeric()
  vals_by_t <- vector("list", length(s))
  for (t in seq_along(s)) {
    vals <- numeric()
    for (m in seq_along(supp_xi)) {
      for (l in seq_along(supp_eta)) {
        vals <- c(vals, if (identical(characteristic, "race")) {
          (s[t]^beta) * supp_xi[m] * supp_eta[l]
        } else {
          mu + (s[t]^beta) * (supp_xi[m] + supp_eta[l])
        })
      }
    }
    vals_by_t[[t]] <- vals
    vals_all <- c(vals_all, vals)
  }
  supp_theta <- seq(min(vals_all), max(vals_all), length.out = supp_pts_theta)
  G <- matrix(0, nrow = supp_pts_theta, ncol = length(s))
  for (t in seq_along(s)) {
    k <- 1L
    for (m in seq_along(supp_xi)) {
      for (l in seq_along(supp_eta)) {
        val <- vals_by_t[[t]][k]
        diffs <- abs(val - supp_theta)
        mindiff <- min(diffs)
        G[, t] <- G[, t] + g_xi[m] * g_eta[l] * (diffs == mindiff)
        k <- k + 1L
      }
    }
  }
  g <- rowMeans(G)
  g <- g / sum(g)
  list(support = supp_theta, g = g, density = g / mean(diff(supp_theta)))
}

.gp_pushfwd_fixture_path <- function(file) {
  p <- system.file("extdata", "fixtures", file, package = "gradepath")
  if (nzchar(p) && file.exists(p)) {
    return(p)
  }
  file.path(.gp_pushfwd_repo_root(), "inst", "extdata", "fixtures", file)
}

.gp_pushfwd_repo_root <- function() {
  candidates <- c(".", "..", "../..", "../../..")
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) &&
        dir.exists(file.path(cand, "inst", "extdata"))) {
      return(normalizePath(cand))
    }
  }
  normalizePath(".")
}

.gp_pushfwd_slow <- function() {
  tolower(Sys.getenv("GRADEPATH_RUN_SLOW_TESTS")) %in% c("1", "true", "yes")
}

test_that("gp_pushforward_theta returns a 250-point normalized theta distribution", {
  skip_if_no_pushforward()
  out <- gp_pushforward_theta(
    supp_xi = c(0, 1, 2),
    g_xi = c(0.2, 0.5, 0.3),
    supp_eta = c(0.5, 1.5),
    g_eta = c(0.4, 0.6),
    s = c(0.5, 1.0),
    mu = 0,
    beta = 0.5,
    characteristic = "race"
  )
  expect_length(out$support, 250L)
  expect_length(out$g, 250L)
  expect_length(out$density, 250L)
  expect_equal(sum(out$g), 1, tolerance = 1e-12)
  expect_true(all(out$density >= 0))
  expect_equal(sum(out$density) * out$diagnostics$grid_width, 1,
               tolerance = 1e-12)
  expect_equal(min(out$support), 0, tolerance = 1e-12)
})

test_that("native pushforward matches the literal get_g_theta loop for race and gender", {
  skip_if_no_pushforward()
  supp_xi <- c(-0.5, 0.25, 1.0)
  g_xi <- c(0.2, 0.3, 0.5)
  supp_eta <- c(-0.25, 0.75)
  g_eta <- c(0.4, 0.6)
  s <- c(0.5, 1.25)
  for (ch in c("race", "gender")) {
    sx <- if (identical(ch, "race")) c(0.25, 0.5, 1.0) else supp_xi
    se <- if (identical(ch, "race")) c(0.5, 1.5) else supp_eta
    out <- gp_pushforward_theta(
      sx, g_xi, se, g_eta, s = s, mu = 0.1, beta = 0.7,
      characteristic = ch, supp_pts_theta = 31L
    )
    ref <- .gp_pushfwd_literal(
      sx, g_xi, se, g_eta, s = s, mu = 0.1, beta = 0.7,
      characteristic = ch, supp_pts_theta = 31L
    )
    expect_equal(out$support, ref$support, tolerance = 1e-12)
    expect_equal(out$g, ref$g, tolerance = 1e-12)
    expect_equal(out$density, ref$density, tolerance = 1e-12)
  }
})

test_that("nearest-grid tie behavior mirrors diffs == min(diffs)", {
  skip_if_no_pushforward()
  col <- .gp_pushfwd_deposit(
    G_col = c(0, 0),
    val = 0.5,
    mass = 1,
    supp_theta = c(0, 1),
    step = 1
  )
  expect_equal(col, c(1, 1))
})

test_that("gender additive transform includes mu and eta + xi", {
  skip_if_no_pushforward()
  out <- gp_pushforward_theta(
    supp_xi = c(-1, 1),
    g_xi = c(0.5, 0.5),
    supp_eta = c(-0.5, 0.5),
    g_eta = c(0.5, 0.5),
    s = 1,
    mu = 2,
    beta = 1,
    characteristic = "gender",
    supp_pts_theta = 5L
  )
  expect_equal(out$support, seq(0.5, 3.5, length.out = 5), tolerance = 1e-12)
})

test_that("degenerate eta slice matches the ebrecipe one-level pushforward", {
  skip_if_no_pushforward()
  skip_if_not_installed("ebrecipe")
  support <- c(0.2, 0.8, 1.5, 2.1)
  g <- c(0.1, 0.2, 0.3, 0.4)
  s <- c(0.45, 0.9)
  psi_1 <- 0.1
  psi_2 <- 0.5

  eb_race <- .gp_eb_pushforward_theta(
    support, g, s = s, psi_1 = psi_1, psi_2 = psi_2,
    characteristic = "white"
  )
  gp_race <- gp_pushforward_theta(
    support, g, supp_eta = exp(psi_1), g_eta = 1,
    s = s, mu = 0, beta = psi_2, characteristic = "race",
    supp_pts_theta = length(support)
  )
  expect_equal(gp_race$support, eb_race$support, tolerance = 1e-12)
  expect_equal(gp_race$g, eb_race$g, tolerance = 1e-12)
  expect_equal(gp_race$density, eb_race$density, tolerance = 1e-12)

  eb_gender <- .gp_eb_pushforward_theta(
    support, g, s = s, psi_1 = psi_1, psi_2 = psi_2,
    characteristic = "male"
  )
  gp_gender <- gp_pushforward_theta(
    support, g, supp_eta = 0, g_eta = 1,
    s = s, mu = psi_1, beta = psi_2, characteristic = "gender",
    supp_pts_theta = length(support)
  )
  expect_equal(gp_gender$support, eb_gender$support, tolerance = 1e-12)
  expect_equal(gp_gender$g, eb_gender$g, tolerance = 1e-12)
  expect_equal(gp_gender$density, eb_gender$density, tolerance = 1e-12)
})

test_that("groupfx1 theta fixtures are support-density tables", {
  for (ch in c("race", "gender")) {
    path <- .gp_pushfwd_fixture_path(sprintf("g_theta_groupfx1_%s.csv", ch))
    skip_if_not(file.exists(path), sprintf("missing %s g_theta fixture", ch))
    x <- as.matrix(utils::read.csv(path, header = FALSE))
    expect_identical(dim(x), c(250L, 2L))
    expect_true(all(is.finite(x)))
    expect_true(all(diff(x[, 1]) > 0))
    expect_true(all(x[, 2] >= 0))
    expect_equal(sum(x[, 2]) * mean(diff(x[, 1])), 1, tolerance = 2e-5)
  }
})

test_that("slow archive supports reproduce groupfx1 theta fixtures approximately", {
  skip_if_not(.gp_pushfwd_slow(),
              "Set GRADEPATH_RUN_SLOW_TESTS=1 for archive g_theta fixture parity.")
  root <- .gp_pushfwd_repo_root()
  dump_root <- file.path(
    root, "log", "029_external-review-prep", "materials", "06_reference-code",
    "krw-companion-public", "dump"
  )

  for (ch in c("race", "gender")) {
    xi_path <- file.path(dump_root, sprintf("g_xi_groupfx1_%s.csv", ch))
    eta_path <- file.path(dump_root, sprintf("g_psi_%s.csv", ch))
    theta_path <- .gp_pushfwd_fixture_path(sprintf("g_theta_groupfx1_%s.csv", ch))
    input_path <- file.path(
      root, "inst", "extdata", "krw-gmm-input",
      sprintf("theta_estimates_matlab_%s.csv", ch)
    )
    skip_if_not(all(file.exists(c(xi_path, eta_path, theta_path, input_path))),
                sprintf("missing archive fixtures for %s", ch))

    xi <- utils::read.csv(xi_path, header = FALSE)
    eta <- utils::read.csv(eta_path, header = FALSE)
    theta <- utils::read.csv(theta_path, header = FALSE)
    input <- utils::read.csv(input_path, header = FALSE)
    mu <- as.numeric(xi[[3]][1])
    beta <- as.numeric(xi[[4]][1])

    out <- gp_pushforward_theta(
      supp_xi = as.numeric(xi[[1]]),
      g_xi = as.numeric(xi[[2]]),
      supp_eta = as.numeric(eta[[1]]),
      g_eta = as.numeric(eta[[2]]),
      s = as.numeric(input[[3]]),
      mu = mu,
      beta = beta,
      characteristic = ch
    )

    expect_equal(out$support, as.numeric(theta[[1]]), tolerance = 5e-4)
    expect_equal(out$density, as.numeric(theta[[2]]), tolerance = 5e-3)
    expect_equal(sum(out$g), 1, tolerance = 1e-12)
  }
})

test_that("invalid pushforward inputs fail with gradepath validation errors", {
  skip_if_no_pushforward()
  expect_error(
    gp_pushforward_theta(c(0, 1), c(1), c(1), c(1), s = 1, beta = 1,
                         characteristic = "race"),
    class = "gradepath_validation_error"
  )
  expect_error(
    gp_pushforward_theta(c(0, 1), c(0.5, 0.5), c(1), c(1), s = -1, beta = 1,
                         characteristic = "race"),
    class = "gradepath_validation_error"
  )
  expect_error(
    gp_pushforward_theta(c(1), c(1), c(1), c(1), s = 1, beta = 1,
                         characteristic = "race"),
    class = "gradepath_validation_error"
  )
})
