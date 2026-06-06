make_api_monolith_fixture <- function() {
  ids <- paste0("u", seq_len(4))
  data <- data.frame(
    theta_hat = c(-0.45, -0.10, 0.20, 0.55),
    s = c(0.20, 0.22, 0.24, 0.26),
    unit_id = ids,
    label = paste("Unit", seq_len(4)),
    stringsAsFactors = FALSE
  )

  coupling <- new_gp_estimation_fit(
    beta = 0.4,
    m_hat = c(0.0, 1.0),
    V_m = diag(c(0.10, 0.20)),
    v_hat = c(-1.1, -0.25, 0.35, 1.0),
    s_v = c(0.55, 0.50, 0.48, 0.52),
    J = 0,
    df = 1L,
    p_value = 1,
    report = list(E_theta = 0, SD_theta = 1),
    caps = c(lo = -2, hi = 2),
    characteristic = "race",
    provenance = list(mu = 0, characteristic = "race")
  )

  support <- c(-1, 0, 1)
  prior <- new_gp_prior(
    support = support,
    density = c(0.25, 0.50, 0.25),
    mean = 0,
    scale = "r",
    diagnostics = list(method = "test-fixture"),
    metadata = list(characteristic = "race")
  )

  W <- matrix(
    c(
      0.80, 0.18, 0.02,
      0.25, 0.55, 0.20,
      0.10, 0.40, 0.50,
      0.02, 0.18, 0.80
    ),
    nrow = 4,
    byrow = TRUE
  )
  weights <- list(
    W = W,
    support = support,
    reporting_support = matrix(rep(support, times = 4), nrow = 3),
    ids = ids
  )

  list(data = data, coupling = coupling, prior = prior, weights = weights)
}

with_mocked_monolith_front_half <- function(parts) {
  testthat::local_mocked_bindings(
    gp_w_seam = function(input, characteristic, control = NULL) {
      validate_gp_estimation_fit(parts$coupling)
      parts$coupling
    },
    gp_deconvolve = function(fit, control = NULL) {
      validate_gp_prior(parts$prior)
      parts$prior
    },
    gp_posterior_weights = function(prior, r_estimates, precision_fit = NULL) {
      parts$weights
    },
    .package = "gradepath",
    .env = parent.frame()
  )
}

make_api_grade_fit <- function(ids,
                               lambda,
                               grades,
                               control,
                               status = "optimal",
                               warnings = character()) {
  new_gp_grade_fit(
    ids = ids,
    lambda = lambda,
    assignment = data.frame(
      id = ids,
      grade = as.integer(grades),
      stringsAsFactors = FALSE
    ),
    summary = list(
      grade_count = as.integer(length(unique(grades))),
      status = status,
      n_units = length(ids)
    ),
    objective = list(value = 0, raw = 0, canonical = 0),
    backend = list(
      name = control$backend,
      path = "fixture",
      status = status
    ),
    control = control,
    warnings = warnings
  )
}

make_api_status_grade_path <- function(ids,
                                       control,
                                       selected_lambda = 0.25,
                                       selection_rule = "baseline_lambda_0.25") {
  warning <- "Solver returned status `gap_reached` at selected lambda."
  selected <- make_api_grade_fit(
    ids = ids,
    lambda = selected_lambda,
    grades = c(2L, 1L, 2L, 3L),
    control = control,
    status = "gap_reached",
    warnings = warning
  )
  endpoint <- make_api_grade_fit(
    ids = ids,
    lambda = 1,
    grades = c(2L, 1L, 1L, 2L),
    control = control,
    status = "optimal"
  )
  fits <- list(selected, endpoint)
  lambda_grid <- c(selected_lambda, 1)

  new_gp_grade_path(
    ids = ids,
    lambda_grid = lambda_grid,
    fits = fits,
    summary = data.frame(
      lambda = lambda_grid,
      grade_count = vapply(fits, function(fit) fit$summary$grade_count, integer(1)),
      status = vapply(fits, function(fit) fit$summary$status, character(1)),
      stringsAsFactors = FALSE
    ),
    backend = list(name = control$backend, path = "fixture"),
    selection = list(
      selected_lambda = selected_lambda,
      selection_rule = selection_rule,
      endpoint_lambda = 1
    ),
    control = control,
    warnings = warning
  )
}

test_that("public exports are present", {
  exports <- getNamespaceExports("gradepath")
  expect_true(all(c(
    "gp_pairwise",
    "gp_preview",
    "krw_report_card",
    "gradepath"
  ) %in% exports))
  expect_identical(gradepath, krw_report_card)
})

test_that("gp_preview validates input and reports the no-solve plan", {
  parts <- make_api_monolith_fixture()
  testthat::local_mocked_bindings(
    gp_w_seam = function(...) stop("GMM should not run in gp_preview()"),
    gp_deconvolve = function(...) stop("deconvolution should not run in gp_preview()"),
    gp_posterior_weights = function(...) stop("posterior weights should not run in gp_preview()"),
    gp_posterior_onelevel = function(...) stop("posterior should not run in gp_preview()"),
    gp_pairwise = function(...) stop("pairwise should not run in gp_preview()"),
    gp_grade_path = function(...) stop("solver should not run in gp_preview()"),
    .package = "gradepath"
  )

  preview <- gp_preview(parts$data, "race", control = gp_control())
  expect_s3_class(preview, "gp_preview")
  expect_silent(validate_gp_preview(preview))
  expect_identical(preview$status$status, "OK")
  expect_identical(preview$workflow, "one_level_independence")
  expect_identical(preview$estimated_solves, length(.gp_grade_operational_default_grid()))
  expect_identical(preview$ids, parts$data$unit_id)

  grouped <- gp_preview(parts$data, "race", groups = c("a", "a", "b", "b"))
  expect_identical(grouped$status$status, "GROUPS_ERROR")
  expect_identical(grouped$workflow, "grouped_pending")

  bad <- parts$data
  bad$s[2] <- NA_real_
  expect_error(
    gp_preview(bad, "race"),
    class = "gp_status_input_error"
  )
  expect_error(
    gp_preview(parts$data, "race", control = list()),
    class = "gp_status_input_error"
  )

  bad_industry <- as.list(parts$data)
  bad_industry$industry <- c("x", "y")
  expect_error(
    gp_preview(bad_industry, "race"),
    class = "gp_status_input_error"
  )
})

test_that("API input rejects ambiguous generic and demographic-specific columns", {
  ids <- paste0("u", seq_len(4))
  demographic_specific <- data.frame(
    theta_hat_race = c(-0.4, -0.1, 0.2, 0.5),
    se_race = c(0.20, 0.21, 0.22, 0.23),
    theta_hat_gender = c(1.1, 1.2, 1.3, 1.4),
    se_gender = c(0.30, 0.31, 0.32, 0.33),
    unit_id = ids,
    label = paste("Unit", seq_len(4)),
    stringsAsFactors = FALSE
  )

  race <- .gp_api_estimates(demographic_specific, "race")
  expect_equal(race$theta_hat, demographic_specific$theta_hat_race)
  expect_equal(race$s, demographic_specific$se_race)
  gender <- .gp_api_estimates(demographic_specific, "gender")
  expect_equal(gender$theta_hat, demographic_specific$theta_hat_gender)
  expect_equal(gender$s, demographic_specific$se_gender)

  generic <- data.frame(
    theta_hat = c(-0.4, -0.1, 0.2, 0.5),
    s = c(0.20, 0.21, 0.22, 0.23),
    unit_id = ids,
    stringsAsFactors = FALSE
  )
  generic_input <- .gp_api_estimates(generic, "race")
  expect_equal(generic_input$theta_hat, generic$theta_hat)
  expect_equal(generic_input$s, generic$s)

  ambiguous_theta <- demographic_specific
  ambiguous_theta$theta_hat <- c(9, 9, 9, 9)
  expect_error(
    .gp_api_estimates(ambiguous_theta, "race"),
    regexp = "both demographic-specific",
    class = "gp_status_input_error"
  )
  expect_error(
    gp_preview(ambiguous_theta, "race"),
    regexp = "both demographic-specific",
    class = "gp_status_input_error"
  )

  ambiguous_estimate <- demographic_specific
  ambiguous_estimate$estimate <- c(9, 9, 9, 9)
  expect_error(
    .gp_api_estimates(ambiguous_estimate, "race"),
    regexp = "both demographic-specific",
    class = "gp_status_input_error"
  )

  ambiguous_se <- demographic_specific
  ambiguous_se$s <- c(9, 9, 9, 9)
  expect_error(
    .gp_api_estimates(ambiguous_se, "race"),
    regexp = "both demographic-specific",
    class = "gp_status_input_error"
  )

  ambiguous_se_alias <- demographic_specific
  ambiguous_se_alias$se <- c(9, 9, 9, 9)
  expect_error(
    .gp_api_estimates(ambiguous_se_alias, "race"),
    regexp = "both demographic-specific",
    class = "gp_status_input_error"
  )
})

test_that("wrong demographic-specific column with no generic fallback errors", {
  # Only theta_hat_gender / se_gender are present (no generic theta_hat/s, no
  # race-specific columns). Asking for `demographic = "race"` resolves NO race
  # estimate: conflict_specific matches a gender column but generic does not, so
  # the conflict guard is correctly silent and the missing required estimate is
  # what must error. This pins the path against a future narrowing of the
  # conflict_specific logic. The code is already correct; this is a guard.
  ids <- paste0("u", seq_len(4))
  gender_only <- data.frame(
    theta_hat_gender = c(1.1, 1.2, 1.3, 1.4),
    se_gender = c(0.30, 0.31, 0.32, 0.33),
    unit_id = ids,
    label = paste("Unit", seq_len(4)),
    stringsAsFactors = FALSE
  )

  # The classed input error is gp_status_input_error (status INPUT_ERROR), with
  # the message naming the missing `theta_hat` estimate.
  err <- expect_error(
    .gp_api_estimates(gender_only, "race"),
    regexp = "theta_hat",
    class = "gp_status_input_error"
  )
  expect_identical(err$status, "INPUT_ERROR")

  expect_error(
    krw_report_card(gender_only, "race"),
    regexp = "theta_hat",
    class = "gp_status_input_error"
  )
  expect_error(
    gp_preview(gender_only, "race"),
    regexp = "theta_hat",
    class = "gp_status_input_error"
  )
})

test_that("krw_report_card() builds a valid gp_fit and equals the piped surface", {
  skip_if_not_installed("highs")
  parts <- make_api_monolith_fixture()
  with_mocked_monolith_front_half(parts)
  control <- gp_control(
    lambda_grid = c(0.25, 1),
    backend = "highs",
    precision_rule = "krw_gmm",
    time_limit = 5,
    mip_gap = 0.05
  )

  fit <- krw_report_card(parts$data, "race", control = control, lambda = 0.25)
  expect_s3_class(fit, "gp_fit")
  expect_silent(validate_gp_fit(fit))
  expect_identical(fit$provenance$producer_status, "OK")
  expect_equal(fit$estimates$original_s, parts$data$s)
  expect_s3_class(fit$precision_fit, "gp_precision_fit")
  expect_identical(fit$precision_fit$parameters$model_form, "multiplicative")

  expect_identical(gp_pairwise(fit), fit$pairwise)
  expect_s3_class(
    gp_pairwise(weights = parts$weights$W,
                reporting_support = parts$weights$reporting_support,
                ids = parts$weights$ids,
                control = control),
    "gp_pairwise"
  )
  expect_error(
    gp_pairwise(get_posterior(fit)),
    regexp = "cannot rebuild",
    class = "gradepath_error"
  )

  path <- gp_grade_path(gp_pairwise(fit), control = get_control(fit),
                        selected_lambda = 0.25)
  selected <- gp_select_grade(path, lambda = 0.25)
  card <- gp_report_card(
    fit$estimates,
    posterior = get_posterior(fit),
    selected_grade = selected,
    grade_path = path
  )
  expect_equal(selected$assignment, fit$selected_grade$assignment)
  expect_equal(card$table, fit$report_card$table)

  bad <- fit
  bad$estimates$original_s <- NULL
  expect_error(
    validate_gp_fit(bad),
    regexp = "original_s"
  )
})

test_that("krw_report_card() propagates selected non-optimal solver status", {
  parts <- make_api_monolith_fixture()
  with_mocked_monolith_front_half(parts)
  control <- gp_control(
    lambda_grid = c(0.25, 1),
    backend = "highs",
    precision_rule = "krw_gmm",
    time_limit = 5,
    mip_gap = 0.05
  )

  testthat::local_mocked_bindings(
    gp_grade_path = function(pairwise,
                             lambda_grid = NULL,
                             control = NULL,
                             selected_lambda = 0.25,
                             selection_rule = "baseline_lambda_0.25",
                             ...) {
      validate_gp_pairwise(pairwise)
      make_api_status_grade_path(
        ids = pairwise$ids,
        control = validate_gp_control(control),
        selected_lambda = selected_lambda,
        selection_rule = selection_rule
      )
    },
    .package = "gradepath"
  )

  fit <- krw_report_card(parts$data, "race", control = control, lambda = 0.25)
  expect_silent(validate_gp_fit(fit))
  expect_identical(fit$selected_grade$backend$status, "gap_reached")
  expect_identical(fit$selected_grade$summary$status, "gap_reached")
  expect_identical(fit$provenance$producer_status, "SOLVER_GAP")
  expect_true(any(grepl("gap_reached", fit$warnings, fixed = TRUE)))
})

test_that("krw_report_card() preserves backend correctness status taxonomy", {
  parts <- make_api_monolith_fixture()
  with_mocked_monolith_front_half(parts)
  control <- gp_control(
    lambda_grid = c(0.25, 1),
    backend = "highs",
    precision_rule = "krw_gmm",
    time_limit = 5,
    mip_gap = 0.05
  )

  testthat::local_mocked_bindings(
    gp_grade_path = function(...) {
      .gradepath_abort_solver_objective_mismatch("mock objective mismatch")
    },
    .package = "gradepath"
  )
  objective_error <- tryCatch(
    krw_report_card(parts$data, "race", control = control, lambda = 0.25),
    error = function(e) e
  )
  expect_s3_class(objective_error, "gp_status_solver_objective_mismatch")
  expect_identical(objective_error$status, "SOLVER_OBJECTIVE_MISMATCH")
  expect_false(identical(objective_error$status, "SOLVER_BACKEND_UNAVAILABLE"))

  testthat::local_mocked_bindings(
    gp_grade_path = function(...) {
      empty_detect <- list(
        gurobi = FALSE,
        callr = FALSE,
        Matrix = FALSE,
        jsonlite = FALSE,
        gurobi_cl = FALSE,
        gurobi_cl_path = "",
        gurobi_cl_smoke = FALSE,
        binding_smoke = FALSE
      )
      .gp_gurobi_abort_unavailable(empty_detect)
    },
    .package = "gradepath"
  )
  unavailable <- tryCatch(
    krw_report_card(parts$data, "race", control = control, lambda = 0.25),
    error = function(e) e
  )
  expect_s3_class(unavailable, "gp_status_solver_backend_unavailable")
  expect_identical(unavailable$status, "SOLVER_BACKEND_UNAVAILABLE")
})

test_that("krw_report_card() rejects unsupported M1 groups and bad control", {
  parts <- make_api_monolith_fixture()
  expect_error(
    krw_report_card(parts$data, "race", groups = c("a", "a", "b", "b")),
    class = "gp_status_groups_error"
  )
  expect_error(
    krw_report_card(parts$data, "race", control = gp_control(backend = "highs")),
    regexp = "precision_rule",
    class = "gp_status_input_error"
  )
  expect_error(
    krw_report_card(
      parts$data,
      "race",
      control = gp_control(
        lambda_grid = c(0, 0.25, 1),
        backend = "highs",
        precision_rule = "krw_gmm"
      ),
      lambda = 0.5
    ),
    regexp = "lambda",
    class = "gp_status_input_error"
  )
})

test_that("krw_report_card() threads acceptance_mode into gp_grade_path()", {
  parts <- make_api_monolith_fixture()
  with_mocked_monolith_front_half(parts)
  control <- gp_control(
    lambda_grid = c(0.25, 1),
    backend = "highs",
    precision_rule = "krw_gmm",
    time_limit = 5,
    mip_gap = 0.05
  )
  seen <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    gp_grade_path = function(pairwise,
                             lambda_grid = NULL,
                             control = NULL,
                             selected_lambda = 0.25,
                             selection_rule = "baseline_lambda_0.25",
                             acceptance_mode = FALSE,
                             ...) {
      validate_gp_pairwise(pairwise)
      seen$acceptance_mode <- acceptance_mode
      make_api_status_grade_path(
        ids = pairwise$ids,
        control = validate_gp_control(control),
        selected_lambda = selected_lambda,
        selection_rule = selection_rule
      )
    },
    .package = "gradepath"
  )

  krw_report_card(parts$data, "race", control = control, lambda = 0.25)
  expect_false(seen$acceptance_mode)

  krw_report_card(
    parts$data, "race",
    control = control, lambda = 0.25, acceptance_mode = TRUE
  )
  expect_true(seen$acceptance_mode)
})

test_that("krw_report_card() validates acceptance_mode as a scalar logical", {
  parts <- make_api_monolith_fixture()
  expect_error(
    krw_report_card(parts$data, "race", acceptance_mode = "yes"),
    regexp = "acceptance_mode",
    class = "gradepath_error"
  )
})

test_that("status helpers normalize solver and API statuses", {
  expect_identical(.gp_status_normalize("optimal"), "OK")
  expect_identical(.gp_status_normalize("time_limit"), "SOLVER_TIME_LIMIT")
  expect_identical(.gp_status_normalize("gap_reached"), "SOLVER_GAP")
  expect_identical(.gp_status_normalize("not-a-status"), "UNVERIFIED")
  err <- tryCatch(
    .gp_status_abort("INPUT_ERROR", "bad input"),
    error = identity
  )
  expect_identical(.gp_status_from_condition(err), "INPUT_ERROR")
})
