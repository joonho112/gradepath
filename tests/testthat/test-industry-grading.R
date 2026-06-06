# =============================================================================
# test-industry-grading.R -- two-level grading/report-card route
# =============================================================================

.tlg_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}

.tlg_pairwise <- function(...) .tlg_get("gp_twolevel_pairwise")(...)
.tlg_grade <- function(...) .tlg_get("gp_twolevel_grade")(...)
.tlg_card <- function(...) .tlg_get("gp_twolevel_report_card")(...)
.tlg_valid_bundle <- function(...) .tlg_get("validate_gp_twolevel_pairwise_bundle")(...)
.tlg_valid_grade <- function(...) .tlg_get("validate_gp_twolevel_grade")(...)
.tlg_quad <- function(...) .tlg_get("gp_twolevel_quadrature")(...)

.tlg_control <- function() {
  gp_control(backend = "highs")
}

.tlg_strict_matrix <- function(n, prefix = "u") {
  ids <- paste0(prefix, seq_len(n))
  m <- matrix(0, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      m[i, j] <- if (i < j) 0.95 else 0.05
    }
  }
  m
}

.tlg_tiny_prior <- function(characteristic = "race") {
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (identical(characteristic, "race")) {
    supp_xi <- c(0.45, 1.25)
    g_xi <- c(0.45, 0.55)
    supp_eta <- c(0.70, 1.45)
    g_eta <- c(0.40, 0.60)
    beta <- 0.55
    mu <- 0
  } else {
    supp_xi <- c(-0.50, 0.40)
    g_xi <- c(0.35, 0.65)
    supp_eta <- c(-0.30, 0.50)
    g_eta <- c(0.55, 0.45)
    beta <- 0.65
    mu <- 0.20
  }
  structure(
    list(
      support = supp_xi,
      density = g_xi,
      mean = sum(supp_xi * g_xi),
      scale = "r",
      diagnostics = list(group_fx = 1L, support_eta = supp_eta, g_eta = g_eta),
      metadata = list(characteristic = characteristic, beta = beta, mu = mu)
    ),
    class = c("gp_prior", "list")
  )
}

.tlg_tiny_case <- function(characteristic = "race") {
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
  prior <- .tlg_tiny_prior(characteristic)
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
      s_v = c(0.30, 0.34, 0.28),
      industry = c(1, 1, 2)
    )
  )
}

test_that("two-level pairwise dispatch applies the grading cleanup contract", {
  bundle <- .tlg_pairwise(
    Pi_theta = .tlg_strict_matrix(4L, "firm_"),
    Pi_bar = .tlg_strict_matrix(3L, "industry_"),
    control = .tlg_control()
  )

  expect_s3_class(bundle, "gp_twolevel_pairwise_bundle")
  expect_silent(.tlg_valid_bundle(bundle))
  expect_s3_class(bundle$pairwise_theta, "gp_pairwise")
  expect_s3_class(bundle$pairwise_bar, "gp_pairwise")
  expect_identical(bundle$pairwise_theta$ids, paste0("firm_", 1:4))
  expect_identical(bundle$pairwise_bar$ids, paste0("industry_", 1:3))
  expect_equal(unname(diag(bundle$Pi_theta)), rep(0.5, 4L), tolerance = 0)
  expect_equal(unname(diag(bundle$Pi_bar)), rep(0.5, 3L), tolerance = 0)
  expect_equal(unname(bundle$Pi_theta + t(bundle$Pi_theta)),
               matrix(1, 4L, 4L), tolerance = 1e-12)
  expect_gte(min(bundle$Pi_theta[row(bundle$Pi_theta) !=
                                  col(bundle$Pi_theta)]),
             .tlg_get(".gp_pairwise_zero_floor"))
  expect_identical(bundle$cleanup$antisymmetry, TRUE)
  expect_identical(bundle$cleanup$diagonal, 0.5)
})

test_that("two-level grade routes Pi_theta to industry_rfe and Pi_bar to btwn", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  bundle <- .tlg_pairwise(
    Pi_theta = .tlg_strict_matrix(4L, "firm_"),
    Pi_bar = .tlg_strict_matrix(3L, "industry_"),
    control = .tlg_control()
  )

  graded <- .tlg_grade(
    bundle,
    control = .tlg_control(),
    lambda_grid = c(0.25, 1),
    build_report_cards = FALSE
  )

  expect_s3_class(graded, "gp_twolevel_grade")
  expect_silent(.tlg_valid_grade(graded))
  expect_identical(graded$industry_rfe$model, "industry_rfe")
  expect_identical(graded$btwn$model, "btwn")
  expect_identical(graded$industry_rfe$pairwise, bundle$pairwise_theta)
  expect_identical(graded$btwn$pairwise, bundle$pairwise_bar)
  expect_identical(graded$industry_rfe$grade_count, 4L)
  expect_identical(graded$btwn$grade_count, 3L)
  expect_identical(graded$industry_rfe$producer_status, "OK")
  expect_identical(graded$btwn$producer_status, "OK")
  expect_identical(graded$producer_status, "OK")
  expect_equal(graded$lambda_grid, c(0.25, 1), tolerance = 0)
  expect_identical(graded$provenance$route,
                   "industry_rfe_Pi_theta__btwn_Pi_bar")
  expect_true(isTRUE(graded$provenance$m1_safe))
})

test_that("two-level grade assembles firm and between-industry report cards when posterior-backed", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  b <- .tlg_tiny_case("race")
  q <- .tlg_quad(
    input = b$input,
    prior = b$prior,
    fit = b$fit,
    include_g_theta = FALSE,
    control = .tlg_control()
  )

  graded <- .tlg_grade(
    q,
    control = .tlg_control(),
    lambda_grid = c(0.25, 1)
  )
  firm_card <- .tlg_card(graded, model = "industry_rfe")
  industry_card <- .tlg_card(graded, model = "btwn")

  expect_silent(.tlg_valid_grade(graded))
  expect_s3_class(firm_card, "gp_report_card")
  expect_s3_class(industry_card, "gp_report_card")
  expect_identical(sort(firm_card$ids), sort(q$pairwise_theta$ids))
  expect_identical(sort(industry_card$ids), sort(q$pairwise_bar$ids))
  expect_identical(nrow(firm_card$table), length(q$pairwise_theta$ids))
  expect_identical(nrow(industry_card$table), length(q$pairwise_bar$ids))
  expect_equal(firm_card$table$selected_lambda,
               rep(0.25, nrow(firm_card$table)), tolerance = 0)
  expect_equal(industry_card$table$selected_lambda,
               rep(0.25, nrow(industry_card$table)), tolerance = 0)
})

test_that("two-level report-card selector refuses grades without posterior cards", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  bundle <- .tlg_pairwise(
    Pi_theta = .tlg_strict_matrix(4L, "firm_"),
    Pi_bar = .tlg_strict_matrix(3L, "industry_"),
    control = .tlg_control()
  )
  graded <- .tlg_grade(
    bundle,
    control = .tlg_control(),
    lambda_grid = c(0.25, 1),
    build_report_cards = FALSE
  )

  expect_error(
    .tlg_card(graded, model = "btwn"),
    class = "gradepath_validation_error"
  )
})
