# =============================================================================
# test-twolevel-deconv-parity.R -- Matlab parity gate
# =============================================================================

.gp_75_slow <- function() {
  tolower(Sys.getenv("GRADEPATH_RUN_SLOW_TESTS")) %in% c("1", "true", "yes")
}

.gp_75_heavy <- function() {
  .gp_75_slow() &&
    tolower(Sys.getenv("GRADEPATH_RUN_HEAVY_TESTS")) %in% c("1", "true", "yes")
}

.gp_75_root <- function() {
  candidates <- c(".", "..", "../..", "../../..")
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) &&
        dir.exists(file.path(cand, "inst", "extdata"))) {
      return(normalizePath(cand))
    }
  }
  normalizePath(".")
}

.gp_75_dump_root <- function() {
  env <- Sys.getenv("GRADEPATH_KRW_ARCHIVE", unset = "")
  candidates <- c(
    if (nzchar(env)) {
      c(env, file.path(env, "dump"), file.path(env, "replication", "dump"))
    } else character(),
    file.path(.gp_75_root(), "..", "dev", "codes",
              "KRW-2024-replication-archive", "replication", "dump"),
    file.path(.gp_75_root(), "..", "dev", "codes",
              "KRW-2024-replication-archive", "replication", "code", "dump"),
    file.path(.gp_75_root(), "log", "029_external-review-prep",
              "materials", "06_reference-code", "krw-companion-public", "dump")
  )
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  hit <- candidates[dir.exists(candidates)][1L]
  if (length(hit) == 0L || is.na(hit)) {
    testthat::skip("KRW archive dump not available.")
  }
  hit
}

.gp_75_fixture <- function(file) {
  p <- system.file("extdata", "fixtures", file, package = "gradepath")
  if (nzchar(p) && file.exists(p)) return(p)
  file.path(.gp_75_root(), "inst", "extdata", "fixtures", file)
}

.gp_75_read <- function(path, ncol, label) {
  testthat::skip_if_not(file.exists(path), sprintf("missing %s", label))
  x <- utils::read.csv(path, header = FALSE)
  testthat::expect_identical(ncol(x), as.integer(ncol))
  testthat::expect_true(all(is.finite(as.matrix(x))))
  x
}

.gp_75_density_mass <- function(support, density) {
  mass <- as.numeric(density) * mean(diff(as.numeric(support)))
  mass / sum(mass)
}

.gp_75_expect_regular_grid <- function(support, tolerance = 2e-4) {
  dx <- diff(as.numeric(support))
  testthat::expect_lte(max(dx) - min(dx), tolerance)
}

.gp_75_moments <- function(support, mass) {
  support <- as.numeric(support)
  mass <- as.numeric(mass) / sum(mass)
  mu <- sum(support * mass)
  sd <- sqrt(sum(((support - mu)^2) * mass))
  c(mean = mu, sd = sd)
}

.gp_75_expect_moments_close <- function(actual, expected, tolerance) {
  testthat::expect_lte(max(abs(as.numeric(actual) - as.numeric(expected))),
                       tolerance)
}

.gp_75_cdf_gap <- function(support_a, mass_a, support_b, mass_b) {
  support_a <- as.numeric(support_a)
  support_b <- as.numeric(support_b)
  mass_a <- as.numeric(mass_a) / sum(mass_a)
  mass_b <- as.numeric(mass_b) / sum(mass_b)
  oa <- order(support_a)
  ob <- order(support_b)
  support_a <- support_a[oa]
  support_b <- support_b[ob]
  mass_a <- mass_a[oa]
  mass_b <- mass_b[ob]
  grid <- sort(unique(c(support_a, support_b)), method = "radix")
  cdf_a <- stats::approxfun(support_a, cumsum(mass_a), method = "constant",
                            f = 0, yleft = 0, yright = 1, ties = "ordered")
  cdf_b <- stats::approxfun(support_b, cumsum(mass_b), method = "constant",
                            f = 0, yleft = 0, yright = 1, ties = "ordered")
  max(abs(cdf_a(grid) - cdf_b(grid)))
}

.gp_75_selected_penalty <- function(ch) {
  grid <- .gp_75_read(
    file.path(.gp_75_dump_root(), sprintf("grid_results_groupfx1_%s.csv", ch)),
    4L,
    sprintf("%s grid_results", ch)
  )
  names(grid) <- c("c_xi", "c_eta", "J", "logL")
  grid[which.min(grid$J), c("c_xi", "c_eta")]
}

.gp_75_archive <- function(ch) {
  dump <- .gp_75_dump_root()
  list(
    xi = .gp_75_read(file.path(dump, sprintf("g_xi_groupfx1_%s.csv", ch)),
                     4L, sprintf("%s g_xi", ch)),
    eta = .gp_75_read(file.path(dump, sprintf("g_psi_%s.csv", ch)),
                      3L, sprintf("%s g_psi", ch)),
    theta = .gp_75_read(file.path(dump, sprintf("g_theta_groupfx1_%s.csv", ch)),
                        2L, sprintf("%s dump g_theta", ch)),
    theta_fixture = .gp_75_read(.gp_75_fixture(sprintf("g_theta_groupfx1_%s.csv", ch)),
                                2L, sprintf("%s package g_theta", ch)),
    penalty = .gp_75_selected_penalty(ch)
  )
}

test_that("archive groupfx1 deconvolution fixtures have the expected structure", {
  skip_if_not(.gp_75_slow(),
              "Set GRADEPATH_RUN_SLOW_TESTS=1 for archive parity.")
  for (ch in c("race", "gender")) {
    a <- .gp_75_archive(ch)
    expect_identical(dim(a$xi), c(200L, 4L))
    expect_identical(dim(a$eta), c(200L, 3L))
    expect_identical(dim(a$theta), c(250L, 2L))
    expect_identical(dim(a$theta_fixture), c(250L, 2L))
    expect_equal(a$theta_fixture, a$theta, tolerance = 0)
    expect_true(all(diff(a$xi[[1]]) > 0))
    expect_true(all(diff(a$eta[[1]]) > 0))
    expect_true(all(diff(a$theta[[1]]) > 0))
    .gp_75_expect_regular_grid(a$xi[[1]])
    .gp_75_expect_regular_grid(a$eta[[1]])
    .gp_75_expect_regular_grid(a$theta[[1]])
    expect_equal(sum(a$xi[[2]]) * mean(diff(a$xi[[1]])), 1, tolerance = 5e-4)
    expect_equal(sum(a$eta[[2]]) * mean(diff(a$eta[[1]])), 1, tolerance = 5e-4)
    expect_equal(sum(a$theta[[2]]) * mean(diff(a$theta[[1]])), 1,
                 tolerance = 5e-4)
  }
})

test_that("archive grid selects the documented two-level penalties", {
  skip_if_not(.gp_75_slow(),
              "Set GRADEPATH_RUN_SLOW_TESTS=1 for archive parity.")
  expect_equal(.gp_75_selected_penalty("race")$c_xi, 0.105, tolerance = 1e-12)
  expect_equal(.gp_75_selected_penalty("race")$c_eta, 0.0025, tolerance = 1e-12)
  expect_equal(.gp_75_selected_penalty("gender")$c_xi, 0.025, tolerance = 1e-12)
  expect_equal(.gp_75_selected_penalty("gender")$c_eta, 0.0025, tolerance = 1e-12)
})

test_that("archive supports reproduce bundled groupfx1 theta fixtures", {
  skip_if_not(.gp_75_slow(),
              "Set GRADEPATH_RUN_SLOW_TESTS=1 for archive parity.")
  for (ch in c("race", "gender")) {
    a <- .gp_75_archive(ch)
    input <- gp_krw_gmm_input(ch)
    xi_mass <- .gp_75_density_mass(a$xi[[1]], a$xi[[2]])
    eta_mass <- .gp_75_density_mass(a$eta[[1]], a$eta[[2]])
    out <- gp_pushforward_theta(
      supp_xi = a$xi[[1]],
      g_xi = xi_mass,
      supp_eta = a$eta[[1]],
      g_eta = eta_mass,
      s = input$s,
      mu = a$xi[[3]][1],
      beta = a$xi[[4]][1],
      characteristic = ch
    )
    expect_equal(out$support, a$theta_fixture[[1]], tolerance = 5e-4)
    expect_equal(out$density, a$theta_fixture[[2]], tolerance = 5e-3)
    expect_equal(sum(out$g), 1, tolerance = 1e-12)
  }
})

test_that("native fixed-penalty two-level deconvolution matches archive moments", {
  skip_if_not(.gp_75_slow(),
              "Set GRADEPATH_RUN_SLOW_TESTS=1 for native parity.")
  for (ch in c("race", "gender")) {
    a <- .gp_75_archive(ch)
    input <- gp_krw_gmm_input(ch)
    fit <- gp_two_level_gmm(input, characteristic = ch)
    native_supp_pts <- if (identical(ch, "race")) 100L else 50L
    prior <- gp_deconvolve_groups(
      fit,
      input,
      control = list(deconv2l_n_knots = 5L, deconv2l_max_iter = 500L,
                     tol = 1e-7),
      n_starts = 1L,
      supp_pts = native_supp_pts,
      penalties = list(c_xi = a$penalty$c_xi, c_eta = a$penalty$c_eta)
    )

    xi_mass <- .gp_75_density_mass(a$xi[[1]], a$xi[[2]])
    eta_mass <- .gp_75_density_mass(a$eta[[1]], a$eta[[2]])
    native_theta <- gp_pushforward_theta(
      prior$support,
      prior$density,
      prior$diagnostics$support_eta,
      prior$diagnostics$g_eta,
      s = input$s,
      mu = a$xi[[3]][1],
      beta = fit$beta,
      characteristic = ch
    )
    theta_mass <- .gp_75_density_mass(a$theta_fixture[[1]],
                                      a$theta_fixture[[2]])

    expect_equal(range(prior$support), range(a$xi[[1]]), tolerance = 5e-3)
    expect_equal(range(prior$diagnostics$support_eta), range(a$eta[[1]]),
                 tolerance = 5e-3)
    expect_equal(range(native_theta$support), range(a$theta_fixture[[1]]),
                 tolerance = 5e-3)
    expect_equal(sum(prior$density), 1, tolerance = 1e-10)
    expect_equal(sum(prior$diagnostics$g_eta), 1, tolerance = 1e-10)
    expect_equal(sum(native_theta$g), 1, tolerance = 1e-10)

    .gp_75_expect_moments_close(
      .gp_75_moments(prior$support, prior$density),
      .gp_75_moments(a$xi[[1]], xi_mass),
      tolerance = 0.025
    )
    .gp_75_expect_moments_close(
      .gp_75_moments(prior$diagnostics$support_eta, prior$diagnostics$g_eta),
      .gp_75_moments(a$eta[[1]], eta_mass),
      tolerance = 0.025
    )
    .gp_75_expect_moments_close(
      .gp_75_moments(native_theta$support, native_theta$g),
      .gp_75_moments(a$theta_fixture[[1]], theta_mass),
      tolerance = 0.025
    )
  }
})
