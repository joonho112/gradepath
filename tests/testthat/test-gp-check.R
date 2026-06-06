test_that("gp_check() applies the registry tolerance boundary", {
  reg <- data.frame(
    id = "x",
    paper_value = "10",
    unit = "count",
    tolerance = "0",
    class = "exact",
    milestone = "M1",
    quantity = "toy",
    stringsAsFactors = FALSE
  )

  expect_identical(gp_check("x", 10, registry = reg)$status, "PASS")
  expect_identical(gp_check("x", 10 + 0.5e-9, registry = reg)$status, "PASS")
  expect_identical(gp_check("x", 10 + 2e-9, registry = reg)$status, "FAIL")
  expect_equal(gp_check("x", 9, registry = reg)$delta, -1)
})

test_that("gp_check() reports NA, no-tolerance, unknown-id, and producer status", {
  reg <- data.frame(
    id = c("missing_paper", "missing_tol", "ok"),
    paper_value = c(NA_character_, "1", "1"),
    unit = c("count", "other", "other"),
    tolerance = c("0", NA_character_, "0.1"),
    class = c("exact", "banded", "banded"),
    milestone = c("M1", "M1", "M1"),
    quantity = c("missing paper", "missing tol", "ok"),
    stringsAsFactors = FALSE
  )

  expect_identical(gp_check("missing_paper", 1, registry = reg)$status, "n/a")
  expect_identical(gp_check("ok", NA_real_, registry = reg)$status, "n/a")
  expect_identical(gp_check("missing_tol", 1, registry = reg)$status, "no-tol")

  unavailable <- gp_check(
    "ok",
    1,
    producer_status = "SOLVER_BACKEND_UNAVAILABLE",
    registry = reg
  )
  expect_identical(unavailable$status, "UNVERIFIED")
  expect_identical(unavailable$reason, "SOLVER_BACKEND_UNAVAILABLE")
  expect_identical(unavailable$producer_status, "SOLVER_BACKEND_UNAVAILABLE")

  expect_error(gp_check("not-present", 1, registry = reg), regexp = "not found")
})

test_that("gp_check() carries unit, milestone, and registry class separately from status", {
  row <- gp_check("race_baseline_dr", 3.9)
  expect_identical(row$unit, "percent")
  expect_identical(row$class, "banded")
  expect_identical(row$milestone, "M1")
  expect_identical(row$status, "PASS")

  # Unit conversion is caller-owned: a proportion is not silently multiplied.
  expect_identical(gp_check("race_baseline_dr", 0.039)$status, "FAIL")
})

test_that("harness adapters read headerless CSV and map arcsine to probability", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(c("1,2,3", "4,5,6"), tmp)
  mat <- .gp_read_headerless_csv(tmp)
  expect_identical(dim(mat), c(2L, 3L))
  expect_identical(names(mat), c("V1", "V2", "V3"))
  expect_equal(as.numeric(mat[2, 3]), 6)

  x <- c(0, pi / 6, pi / 2)
  expect_equal(.gp_asin2p(x), sin(x)^2)
})

test_that("gp_r2() matches a tiny hand-computable fixture", {
  dat <- data.frame(
    estimate = c(-1, -0.5, 0.5, 1),
    se = rep(0.1, 4),
    grade = c(1, 1, 2, 2)
  )
  total_var <- stats::var(dat$estimate) - mean(dat$se^2)
  group_mean <- tapply(dat$estimate, dat$grade, mean)
  group_n <- as.numeric(tabulate(dat$grade))
  overall <- stats::weighted.mean(group_mean, group_n)
  between <- sum(group_n * (group_mean - overall)^2) / (sum(group_n) - 1)
  expect_equal(gp_r2(dat, scale = "proportion"), between / total_var)
  expect_equal(gp_r2(dat), 100 * between / total_var)

  expect_true(is.na(gp_r2(data.frame(
    estimate = c(1, 1, 1),
    se = c(0.1, 0.1, 0.1),
    grade = c(1, 1, 2)
  ))))
})

test_that("gp_krw_r2() ports the KRW second-moment R2 recipe", {
  theta <- data.frame(
    log_dif = c(-1, -0.5, 0.5, 1),
    log_dif_se = rep(0.1, 4),
    stringsAsFactors = FALSE
  )
  post_mean_beta <- c(-1, -0.5, 0.5, 1)
  ranking <- data.frame(
    obs_idx = 1:4,
    post_mean_beta = post_mean_beta,
    pmean_sq = post_mean_beta^2 + 0.04,
    check.names = FALSE
  )
  ranking[["grades_lamb0.25"]] <- c(1, 1, 2, 2)

  out <- gp_krw_r2(theta, ranking)
  denominator <- stats::var(theta$log_dif) - mean(theta$log_dif_se^2)
  between <- 0.5725
  expect_equal(out$n, 4)
  expect_equal(out$grade_count, 2)
  expect_identical(out$grade_col, "grades_lamb0.25")
  expect_equal(out$overall_variance, denominator, tolerance = 1e-12)
  expect_equal(out$between_variance, between, tolerance = 1e-12)
  expect_equal(out$between_sd, sqrt(between), tolerance = 1e-12)
  expect_equal(out$r2_proportion, between / denominator, tolerance = 1e-12)
  expect_equal(out$r2, 100 * between / denominator, tolerance = 1e-12)
  expect_identical(out$scale, "percent")

  ranking_stata <- ranking
  ranking_stata[["grades_lamb025"]] <- ranking_stata[["grades_lamb0.25"]]
  ranking_stata[["grades_lamb0.25"]] <- NULL
  prop <- gp_krw_r2(theta, ranking_stata, scale = "proportion")
  expect_identical(prop$grade_col, "grades_lamb025")
  expect_equal(prop$r2, out$r2_proportion, tolerance = 1e-12)

  expect_error(
    gp_krw_r2(theta, ranking[, setdiff(names(ranking), "pmean_sq")]),
    regexp = "pmean_sq"
  )
  bad_theta <- theta
  bad_theta$log_dif <- rep(1, 4)
  expect_error(gp_krw_r2(bad_theta, ranking), regexp = "denominator")
})

test_that("gp_krw_r2() matches KRW companion rsquared.do golden outputs", {
  root_candidates <- c(
    file.path("..", "KRW-2024-companion-public"),
    file.path("..", "..", "KRW-2024-companion-public"),
    file.path("..", "..", "..", "KRW-2024-companion-public")
  )
  root_ok <- vapply(
    root_candidates,
    function(path) {
      file.exists(file.path(path, "code", "rsquared.do")) &&
        file.exists(file.path(path, "data", "theta_estimates_race.csv")) &&
        file.exists(file.path(path, "dump", "ranking_results_log_dif_binary_race.csv"))
    },
    logical(1)
  )
  skip_if_not(
    any(root_ok),
    "KRW companion archive is not available"
  )
  root <- root_candidates[which(root_ok)[[1L]]]

  read_ref <- function(path) {
    utils::read.csv(file.path(root, path), stringsAsFactors = FALSE, check.names = FALSE)
  }
  calc <- function(demographic, model) {
    gp_krw_r2(
      read_ref(file.path("data", sprintf("theta_estimates_%s.csv", demographic))),
      read_ref(file.path("dump", sprintf("ranking_results_%s_%s.csv", model, demographic)))
    )
  }

  race_base <- calc("race", "log_dif_binary")
  gender_base <- calc("gender", "log_dif_binary")
  race_industry <- calc("race", "industry_rfe_binary")
  gender_industry <- calc("gender", "industry_rfe_binary")

  expect_equal(race_base$r2, 24.5972214331222, tolerance = 1e-10)
  expect_equal(race_base$between_sd, 0.0339817975178244, tolerance = 1e-12)
  expect_equal(gender_base$r2, 43.7177861484058, tolerance = 1e-10)
  expect_equal(race_industry$r2, 70.1403256675552, tolerance = 1e-10)
  expect_equal(gender_industry$r2, 37.7055522895268, tolerance = 1e-10)

  expect_identical(gp_check("race_baseline_r2", race_base$r2)$status, "PASS")
  expect_identical(
    gp_check("race_baseline_betweengrade_sd", race_base$between_sd)$status,
    "PASS"
  )
  expect_identical(gp_check("gender_baseline_r2", gender_base$r2)$status, "PASS")
  expect_identical(gp_check("race_industry_r2", race_industry$r2)$status, "PASS")
  expect_identical(gp_check("gender_industry_r2", gender_industry$r2)$status, "PASS")

  race_prop <- gp_krw_r2(
    read_ref(file.path("data", "theta_estimates_race.csv")),
    read_ref(file.path("dump", "ranking_results_log_dif_binary_race.csv")),
    scale = "proportion"
  )
  expect_identical(gp_check("race_baseline_r2", race_prop$r2)$status, "FAIL")
})
