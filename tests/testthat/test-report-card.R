report_ids <- c("u1", "u2", "u3", "u4")

report_control <- function() {
  gp_control(backend = "highs")
}

report_grade_fit <- function(lambda, grades, ids = report_ids, control = report_control()) {
  new_gp_grade_fit(
    ids = ids,
    lambda = lambda,
    assignment = data.frame(id = ids, grade = as.integer(grades), stringsAsFactors = FALSE),
    summary = list(grade_count = as.integer(length(unique(grades)))),
    objective = list(value = 0, raw = 0, canonical = 0),
    backend = list(name = control$backend, path = "fixture"),
    control = control
  )
}

report_grade_path <- function(ids = report_ids, control = report_control()) {
  selected <- report_grade_fit(0.25, c(2L, 1L, 2L, 3L), ids, control)
  endpoint <- report_grade_fit(1, c(2L, 1L, 1L, 2L), ids, control)
  fits <- list(selected, endpoint)
  lambda_grid <- c(0.25, 1)
  new_gp_grade_path(
    ids = ids,
    lambda_grid = lambda_grid,
    fits = fits,
    summary = data.frame(
      lambda = lambda_grid,
      grade_count = vapply(fits, function(fit) fit$summary$grade_count, integer(1))
    ),
    backend = list(name = control$backend),
    selection = list(
      selected_lambda = 0.25,
      selection_rule = "baseline_lambda_0.25",
      endpoint_lambda = 1
    ),
    control = control
  )
}

report_grade_path_nonbaseline <- function(ids = report_ids, control = report_control()) {
  selected <- report_grade_fit(0.5, c(2L, 1L, 2L, 3L), ids, control)
  endpoint <- report_grade_fit(1, c(2L, 1L, 1L, 2L), ids, control)
  fits <- list(selected, endpoint)
  lambda_grid <- c(0.5, 1)
  new_gp_grade_path(
    ids = ids,
    lambda_grid = lambda_grid,
    fits = fits,
    summary = data.frame(
      lambda = lambda_grid,
      grade_count = vapply(fits, function(fit) fit$summary$grade_count, integer(1))
    ),
    backend = list(name = control$backend),
    selection = list(
      selected_lambda = 0.5,
      selection_rule = "nonbaseline_fixture",
      endpoint_lambda = 1
    ),
    control = control
  )
}

report_posterior <- function(ids = report_ids, labels = paste("Posterior", ids)) {
  pm <- c(0.4, 0.1, 0.2, 0.3)
  new_gp_posterior(
    estimate = c(4, 1, 2, 3),
    se = rep(0.2, length(ids)),
    id = ids,
    label = labels,
    posterior_mean = pm,
    posterior_sd = rep(0.05, length(ids)),
    lower = pm - 0.01,
    upper = pm + 0.01,
    scale = "r"
  )
}

report_estimates <- function(ids = report_ids, labels = paste("Label", ids)) {
  list(
    unit_id = ids,
    label = labels,
    theta_hat = c(10, 20, 30, 40),
    s = rep(0.4, length(ids))
  )
}

report_pairwise <- function(ids = report_ids, control = report_control()) {
  m <- matrix(0.5, length(ids), length(ids), dimnames = list(ids, ids))
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      if (i != j) {
        m[i, j] <- if (i < j) 0.6 else 0.4
      }
    }
  }
  new_gp_pairwise(
    ids = ids,
    matrix = m,
    power = 0L,
    cleanup = list(
      antisymmetry = TRUE,
      diagonal = .gp_pairwise_diagonal,
      zero_floor = .gp_pairwise_zero_floor
    ),
    source = list(
      stage = "posterior",
      rule = "outer_product",
      assumption = "one_level_independence"
    ),
    control = control
  )
}

report_prior <- function() {
  support <- seq(-1, 1, length.out = 11)
  density <- rep(1 / length(support), length(support))
  new_gp_prior(
    support = support,
    density = density,
    mean = sum(support * density),
    scale = "r"
  )
}

report_fit <- function() {
  path <- report_grade_path()
  posterior <- report_posterior()
  card <- gp_report_card(
    report_estimates(),
    posterior = posterior,
    selected_grade = path$fits[[1L]],
    grade_path = path
  )
  new_gp_fit(
    ids = report_ids,
    estimates = ebrecipe::eb_input(
      theta_hat = c(10, 20, 30, 40),
      s = rep(0.4, length(report_ids)),
      unit_id = report_ids
    ),
    prior = report_prior(),
    posterior = posterior,
    precision_fit = NULL,
    pairwise = report_pairwise(),
    grade_path = path,
    selected_grade = path$fits[[1L]],
    report_card = card,
    control = report_control()
  )
}

test_that("gp_report_card() assembles endpoint-sorted rows with selected grades", {
  path <- report_grade_path()
  card <- gp_report_card(
    report_estimates(),
    posterior = report_posterior(),
    selected_grade = path$fits[[1L]],
    grade_path = path
  )

  expect_s3_class(card, "gp_report_card")
  expect_identical(names(card$table), .gp_report_card_table_columns)
  expect_identical(card$table$id, c("u2", "u3", "u1", "u4"))
  expect_identical(card$table$sort_rank, seq_len(4L))
  expect_identical(card$table$grade, c(1L, 2L, 2L, 3L))
  expect_equal(card$table$selected_lambda, rep(0.25, 4L))
  expect_identical(card$table$label, paste("Label", card$table$id))
  expect_equal(card$table$posterior_mean, c(0.1, 0.2, 0.4, 0.3))
  expect_equal(card$table$estimate, c(20, 30, 10, 40))
})

test_that("endpoint ties use id as the secondary sort key", {
  path <- report_grade_path()
  card <- gp_report_card(
    report_estimates(),
    posterior = report_posterior(),
    selected_grade = path$fits[[1L]],
    grade_path = path
  )

  expect_identical(card$table$id[1:2], c("u2", "u3"))
})

test_that("gp_report_card() requires the selected grade to be the stored path member", {
  path <- report_grade_path()
  selected <- path$fits[[1L]]
  selected$backend$path <- "same_grades_but_not_the_stored_fit"

  expect_error(
    gp_report_card(
      report_estimates(),
      posterior = report_posterior(),
      selected_grade = selected,
      grade_path = path
    ),
    regexp = "stored path member",
    class = "gradepath_error"
  )
})

test_that("gp_report_card() rejects id drift across stage-wise inputs", {
  path <- report_grade_path()
  posterior <- report_posterior(ids = rev(report_ids))

  expect_error(
    gp_report_card(
      report_estimates(),
      posterior = posterior,
      selected_grade = path$fits[[1L]],
      grade_path = path
    ),
    regexp = "align exactly",
    class = "gradepath_error"
  )

  expect_error(
    gp_report_card(
      report_estimates(ids = rev(report_ids)),
      posterior = report_posterior(),
      selected_grade = path$fits[[1L]],
      grade_path = path
    ),
    regexp = "estimates.*align",
    class = "gradepath_error"
  )
})

test_that("gp_report_card() uses reporting-scale posterior summaries when available", {
  path <- report_grade_path()
  posterior <- report_posterior()
  posterior$metadata$reporting <- list(
    posterior_mean = c(40, 10, 20, 30),
    posterior_sd = rep(1, 4),
    lower = c(39, 9, 19, 29),
    upper = c(41, 11, 21, 31),
    scale = "theta",
    level = 0.9
  )

  card <- gp_report_card(
    report_estimates(),
    posterior = posterior,
    selected_grade = path$fits[[1L]],
    grade_path = path
  )

  expect_equal(card$table$posterior_mean, c(10, 20, 40, 30))
  expect_equal(card$table$lower, c(9, 19, 39, 29))
  expect_equal(card$table$upper, c(11, 21, 41, 31))
})

test_that("reporting-scale cards validate inside gp_fit", {
  fit <- report_fit()
  fit$posterior$metadata$reporting <- list(
    posterior_mean = c(40, 10, 20, 30),
    posterior_sd = rep(1, 4),
    lower = c(39, 9, 19, 29),
    upper = c(41, 11, 21, 31),
    scale = "theta",
    level = 0.9
  )
  fit$report_card <- gp_report_card(
    report_estimates(),
    posterior = fit$posterior,
    selected_grade = fit$selected_grade,
    grade_path = fit$grade_path
  )

  expect_silent(validate_gp_fit(fit))
})

test_that("gp_report_card(path, posterior=) uses the path selected lambda", {
  path <- report_grade_path_nonbaseline()
  card <- gp_report_card(path, posterior = report_posterior())

  expect_equal(card$selected_lambda, 0.5)
  expect_equal(card$table$selected_lambda, rep(0.5, 4L))
})

test_that("gp_report_card() can omit estimates and falls back to posterior columns", {
  path <- report_grade_path()
  card <- gp_report_card(
    posterior = report_posterior(labels = paste("Fallback", report_ids)),
    selected_grade = path$fits[[1L]],
    grade_path = path
  )

  expect_identical(card$table$label, paste("Fallback", card$table$id))
  expect_equal(card$table$estimate, c(1, 2, 4, 3))
  expect_equal(card$table$se, rep(0.2, 4L))
})

test_that("gp_report_card() reads a stored report card from gp_fit", {
  fit <- report_fit()
  expect_silent(validate_gp_fit(fit))

  card <- gp_report_card(fit)

  expect_identical(card, fit$report_card)
  expect_identical(as.data.frame(fit), fit$report_card$table)
})

test_that("validate_gp_fit() enforces selected fit identity and estimate/se agreement", {
  fit <- report_fit()
  cloned <- fit$selected_grade
  cloned$backend$path <- "same_grades_different_fit"
  fit$selected_grade <- cloned
  expect_error(validate_gp_fit(fit), "stored path member", class = "gradepath_error")

  fit2 <- report_fit()
  fit2$report_card$table$estimate[1] <- fit2$report_card$table$estimate[1] + 1
  expect_error(validate_gp_fit(fit2), "estimate must agree", class = "gradepath_error")

  fit3 <- report_fit()
  fit3$report_card$table$se[1] <- fit3$report_card$table$se[1] + 1
  expect_error(validate_gp_fit(fit3), "se must agree", class = "gradepath_error")
})

test_that("KRW report-card source-table adapter checks registry units and producer status", {
  f5 <- data.frame(
    firm_name = c("Genuine Parts (Napa Auto)", "Charter / Spectrum"),
    log_dif = c(0.3303573, -0.0322733),
    post_mean_beta = c(0.25011, -0.043934),
    ind_post_mean_beta = c(0.23469, -0.027109),
    grade1 = c(1, 97),
    grade1_ind = c(4, 73),
    stringsAsFactors = FALSE
  )
  f6 <- data.frame(
    firm_name = c("Builders FirstSource", "Ascena (Ann Taylor / Loft)"),
    log_dif = c(1.568249, -0.6610414),
    post_mean_beta = c(0.89603, -0.23168),
    ind_post_mean_beta = c(0.67345, -0.18021),
    grade1 = c(1, 96),
    grade1_ind = c(1, 90),
    stringsAsFactors = FALSE
  )
  reg <- data.frame(
    id = c("f5_genuineparts_theta", "f5_charterspectrum_condrank_baseline"),
    paper_value = c("0.33", "97"),
    unit = c("log_diff", "count"),
    tolerance = c("0.01", "0"),
    class = c("banded", "exact"),
    milestone = c("M1", "M1"),
    quantity = c("theta", "rank"),
    stringsAsFactors = FALSE
  )

  out <- .gp_krw_report_card_targets(
    f5,
    f6,
    targets = c("f5_genuineparts_theta", "f5_charterspectrum_condrank_baseline"),
    registry = reg
  )
  expect_equal(out$replicated, c(0.3303573, 97))
  expect_identical(out$status, c("PASS", "PASS"))
  expect_identical(out$unit, c("log_diff", "count"))

  unavailable <- .gp_krw_report_card_targets(
    f5,
    f6,
    targets = "f5_genuineparts_theta",
    producer_status = "SOLVER_GAP",
    registry = reg
  )
  expect_identical(unavailable$status, "UNVERIFIED")
  expect_identical(unavailable$reason, "SOLVER_GAP")

  expect_error(
    .gp_krw_report_card_targets(f5[-1, ], f6, targets = "f5_genuineparts_theta", registry = reg),
    regexp = "matches"
  )
  expect_error(
    .gp_krw_report_card_targets(
      rbind(f5[1, ], f5[1, ]),
      f6,
      targets = "f5_genuineparts_theta",
      registry = reg
    ),
    regexp = "Expected exactly one firm"
  )
  expect_error(
    .gp_krw_report_card_targets(f5, f6, targets = "not_a_registered_f5_id", registry = reg),
    regexp = "Unknown"
  )
})
