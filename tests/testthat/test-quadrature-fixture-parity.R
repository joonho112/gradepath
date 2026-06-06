# =============================================================================
# test-quadrature-fixture-parity.R  --  fixture parity gate
# -----------------------------------------------------------------------------
# Pins the M2 decision input: deterministic quadrature can promote L02 continuous
# targets only if produced Pi_theta, posteriors, and g_theta match the archive
# group_fx == 1 fixtures inside the L02 band.
# =============================================================================

.tlfg_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}

.quad85 <- function(...) .tlfg_get("gp_twolevel_quadrature")(...)
.gate85 <- function(...) .tlfg_get("gp_twolevel_fixture_gate")(...)
.valid_gate85 <- function(...) .tlfg_get("validate_gp_twolevel_fixture_gate")(...)

.gp_85_slow <- function() {
  tolower(Sys.getenv("GRADEPATH_RUN_SLOW_TESTS")) %in% c("1", "true", "yes")
}

.gp_85_heavy <- function() {
  .gp_85_slow() &&
    tolower(Sys.getenv("GRADEPATH_RUN_HEAVY_TESTS")) %in% c("1", "true", "yes")
}

.gp_85_root <- function() {
  candidates <- c(".", "..", "../..", "../../..")
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) &&
        dir.exists(file.path(cand, "inst", "extdata"))) {
      return(normalizePath(cand))
    }
  }
  normalizePath(".")
}

.gp_85_dump_root <- function() {
  root <- .gp_85_root()
  env <- Sys.getenv("GRADEPATH_KRW_ARCHIVE", unset = "")
  candidates <- c(
    if (nzchar(env)) {
      c(env, file.path(env, "dump"), file.path(env, "replication", "dump"))
    } else character(),
    file.path(root, "log", "029_external-review-prep",
              "materials", "06_reference-code", "krw-companion-public",
              "dump"),
    file.path(root, "log", "029_external-review-prep",
              "materials", "06_reference-code", "krw-replication-archive",
              "dump")
  )
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  hit <- candidates[dir.exists(candidates)][1L]
  if (length(hit) == 0L || is.na(hit)) {
    testthat::skip("KRW archive dump not available.")
  }
  hit
}

.gp_85_fixture_dir <- function() {
  p <- system.file("extdata", "fixtures", package = "gradepath")
  if (nzchar(p) && dir.exists(p)) return(p)
  file.path(.gp_85_root(), "inst", "extdata", "fixtures")
}

.gp_85_density_mass <- function(support, density) {
  mass <- as.numeric(density) * mean(diff(as.numeric(support)))
  mass / sum(mass)
}

.gp_85_check_gap <- function(gate, artifact) {
  vals <- gate$checks$max_abs[gate$checks$artifact == artifact]
  if (length(vals) != 1L) {
    stop(sprintf("Expected one `%s` gap.", artifact), call. = FALSE)
  }
  vals[[1L]]
}

.tiny_prior_85 <- function() {
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

.tiny_case_85 <- function() {
  s <- c(0.80, 1.35, 1.05)
  beta <- 0.55
  v_hat <- c(0.72, 1.34, 1.08)
  list(
    input = list(
      theta_hat = v_hat * s^beta,
      s = s,
      industry = c(1, 1, 2),
      unit_id = paste0("u", 1:3),
      label = paste0("unit ", 1:3)
    ),
    prior = .tiny_prior_85(),
    fit = list(
      characteristic = "race",
      beta = beta,
      mu = 0,
      v_hat = v_hat,
      s_v = c(0.30, 0.34, 0.28),
      industry = c(1, 1, 2)
    )
  )
}

.write_gate_fixtures <- function(dir, characteristic, artifacts) {
  utils::write.table(artifacts$Pi_theta,
                     file.path(dir, sprintf("Pi_groupfx1_%s.csv", characteristic)),
                     sep = ",", row.names = FALSE, col.names = FALSE)
  utils::write.table(artifacts$posteriors,
                     file.path(dir, sprintf("posteriors_groupfx1_%s.csv", characteristic)),
                     sep = ",", row.names = FALSE, col.names = FALSE)
  g <- cbind(artifacts$g_theta$support, artifacts$g_theta$density)
  utils::write.table(g,
                     file.path(dir, sprintf("g_theta_groupfx1_%s.csv", characteristic)),
                     sep = ",", row.names = FALSE, col.names = FALSE)
}

.archive_artifacts_85 <- function(ch) {
  dump <- .gp_85_dump_root()
  xi <- utils::read.csv(file.path(dump, sprintf("g_xi_groupfx1_%s.csv", ch)),
                        header = FALSE)
  eta <- utils::read.csv(file.path(dump, sprintf("g_psi_%s.csv", ch)),
                         header = FALSE)
  input <- gp_krw_gmm_input(ch)
  beta <- as.numeric(xi[[4]][1])
  mu <- as.numeric(xi[[3]][1])
  prior <- structure(
    list(
      support = as.numeric(xi[[1]]),
      density = .gp_85_density_mass(xi[[1]], xi[[2]]),
      mean = NA_real_,
      scale = "r",
      diagnostics = list(
        group_fx = 1L,
        support_eta = as.numeric(eta[[1]]),
        g_eta = .gp_85_density_mass(eta[[1]], eta[[2]])
      ),
      metadata = list(characteristic = ch, beta = beta, mu = mu)
    ),
    class = c("gp_prior", "list")
  )
  fit <- list(characteristic = ch, beta = beta, mu = mu,
              industry = input$industry)
  post <- gp_posterior_twolevel(input, prior, fit, interval_level = 0.95)
  g_theta <- gp_pushforward_theta(
    supp_xi = prior$support,
    g_xi = prior$density,
    supp_eta = prior$diagnostics$support_eta,
    g_eta = prior$diagnostics$g_eta,
    s = input$s,
    mu = mu,
    beta = beta,
    characteristic = ch
  )
  list(
    input = input,
    prior = prior,
    fit = fit,
    artifacts = list(
      posteriors = post$metadata$reporting$posteriors,
      Pi_theta = NULL,
      g_theta = g_theta
    )
  )
}

test_that("fixture gate passes a quadrature object against matching fixtures", {
  b <- .tiny_case_85()
  q <- .quad85(b$input, b$prior, b$fit, supp_pts_theta = 31L)
  dir <- tempfile("gp-tlfg-")
  dir.create(dir)
  .write_gate_fixtures(dir, "race", q$artifacts)

  gate <- .gate85(q, characteristic = "race", fixture_dir = dir)
  expect_s3_class(gate, "gp_twolevel_fixture_gate")
  expect_silent(.valid_gate85(gate))
  expect_true(gate$pass)
  expect_identical(gate$class_decision, "banded_candidate")
  expect_identical(gate$producer_status, "OK")
  expect_true(all(gate$checks$status == "PASS"))
  expect_lte(max(gate$checks$max_abs), 1e-12)
})

test_that("fixture gate keeps L02 approximate on failed or skipped checks", {
  b <- .tiny_case_85()
  q <- .quad85(b$input, b$prior, b$fit, supp_pts_theta = 31L)
  dir <- tempfile("gp-tlfg-")
  dir.create(dir)
  bad_artifacts <- q$artifacts
  bad_artifacts$posteriors$posterior_mean[1] <-
    bad_artifacts$posteriors$posterior_mean[1] + 0.25
  .write_gate_fixtures(dir, "race", bad_artifacts)

  gate_fail <- .gate85(q, characteristic = "race", fixture_dir = dir)
  expect_false(gate_fail$pass)
  expect_identical(gate_fail$class_decision, "approximate")
  expect_identical(gate_fail$producer_status, "APPROXIMATE_OK")
  expect_true(any(gate_fail$checks$status == "FAIL"))

  gate_skip <- .gate85(q, characteristic = "race", fixture_dir = dir,
                       require_pi = FALSE)
  expect_false(gate_skip$pass)
  expect_identical(gate_skip$class_decision, "approximate")
  expect_identical(gate_skip$checks$status[gate_skip$checks$artifact == "Pi_theta"],
                   "SKIP")
})

test_that("bundled groupfx1 fixtures have gate shapes", {
  dir <- .gp_85_fixture_dir()
  for (ch in c("race", "gender")) {
    pi <- as.matrix(utils::read.csv(file.path(dir, sprintf("Pi_groupfx1_%s.csv", ch)),
                                    header = FALSE))
    post <- as.matrix(utils::read.csv(file.path(dir, sprintf("posteriors_groupfx1_%s.csv", ch)),
                                      header = FALSE))
    gt <- as.matrix(utils::read.csv(file.path(dir, sprintf("g_theta_groupfx1_%s.csv", ch)),
                                    header = FALSE))
    expect_identical(dim(pi), c(97L, 97L))
    expect_identical(dim(post), c(97L, 4L))
    expect_identical(dim(gt), c(250L, 2L))
    expect_true(all(is.finite(pi)))
    expect_true(all(is.finite(post)))
    expect_true(all(is.finite(gt)))
    expect_equal(diag(pi), rep(0, 97L), tolerance = 0)
    off <- row(pi) != col(pi)
    expect_lte(max(abs((pi + t(pi))[off] - 1)), 2e-4)
    expect_true(all(diff(gt[, 1]) > 0))
    expect_true(all(gt[, 2] >= 0))
    expect_equal(sum(gt[, 2]) * mean(diff(gt[, 1])), 1, tolerance = 5e-4)
  }
})

test_that("archive posterior and g_theta gaps match recorded M2 fixture snapshot", {
  skip_if_not(.gp_85_slow(),
              "Set GRADEPATH_RUN_SLOW_TESTS=1 for archive artifact parity.")
  recorded <- gp_m2_acceptance()$fixture_gates
  for (ch in c("race", "gender")) {
    a <- .archive_artifacts_85(ch)
    gate <- .gate85(
      artifacts = a$artifacts,
      characteristic = ch,
      fixture_dir = .gp_85_fixture_dir(),
      require_pi = FALSE
    )
    expect_silent(.valid_gate85(gate))
    checked <- gate$checks[gate$checks$status != "SKIP", ]
    expect_true(all(checked$status == "PASS"))
    expect_lte(checked$max_abs[checked$artifact == "posteriors"],
               gate$tolerances$posteriors)
    expect_lte(checked$max_abs[checked$artifact == "g_theta_support"],
               gate$tolerances$g_theta_support)
    expect_lte(checked$max_abs[checked$artifact == "g_theta_density"],
               gate$tolerances$g_theta_density)
    expect_false(gate$pass)
    expect_identical(gate$class_decision, "approximate")

    snap <- recorded[recorded$characteristic == ch, , drop = FALSE]
    expect_equal(nrow(snap), 1L)
    expect_identical(
      gate$checks$status[gate$checks$artifact == "Pi_theta"],
      "SKIP"
    )
    expect_lt(abs(.gp_85_check_gap(gate, "posteriors") -
                    snap$posteriors_max_abs), 5e-10)
    expect_lt(abs(.gp_85_check_gap(gate, "g_theta_support") -
                    snap$g_theta_support_max_abs), 5e-10)
    expect_lt(abs(.gp_85_check_gap(gate, "g_theta_density") -
                    snap$g_theta_density_max_abs), 5e-10)
  }
})

test_that("heavy full quadrature gate can compare Pi_theta against fixtures", {
  skip_if_not(.gp_85_heavy(),
              "Set GRADEPATH_RUN_SLOW_TESTS=1 and GRADEPATH_RUN_HEAVY_TESTS=1 for full Pi gate.")
  recorded <- gp_m2_acceptance()$fixture_gates
  for (ch in c("race", "gender")) {
    a <- .archive_artifacts_85(ch)
    q <- .quad85(
      input = a$input,
      prior = a$prior,
      fit = a$fit,
      interval_level = 0.95
    )
    gate <- .gate85(q, characteristic = ch, fixture_dir = .gp_85_fixture_dir())
    expect_silent(.valid_gate85(gate))
    expect_false(any(gate$checks$status == "SKIP"))
    expect_true(all(is.finite(gate$checks$max_abs)))

    pi_gap <- .gp_85_check_gap(gate, "Pi_theta")
    pi_status <- gate$checks$status[gate$checks$artifact == "Pi_theta"]
    snap <- recorded[recorded$characteristic == ch, , drop = FALSE]
    expect_equal(nrow(snap), 1L)
    expect_lt(abs(pi_gap - snap$pi_theta_max_abs), 5e-10)
    expect_equal(gate$tolerances$Pi_theta, 0.01, tolerance = 0)

    if (identical(ch, "race")) {
      expect_gt(pi_gap, 0.01)
      expect_identical(pi_status, "FAIL")
      expect_false(gate$pass)
      expect_identical(gate$class_decision, "approximate")
      expect_identical(gate$producer_status, "APPROXIMATE_OK")
    } else {
      expect_identical(ch, "gender")
      expect_lte(pi_gap, 0.01 + 1e-12)
      expect_identical(pi_status, "PASS")
      expect_true(gate$pass)
      expect_identical(gate$class_decision, "banded_candidate")
      expect_identical(gate$producer_status, "OK")
      expect_true(all(gate$checks$status == "PASS"))
    }
  }
})
