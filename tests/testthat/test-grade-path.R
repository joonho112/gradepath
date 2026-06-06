test_that("gp_grade() solves one lambda using the pairwise control backend", {
  skip_if_not_installed("digest")
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  fixture <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")

  fit <- gp_grade(fixture$pairwise, lambda = 1)

  expect_s3_class(fit, "gp_grade_fit")
  expect_identical(fit$lambda, 1)
  expect_identical(fit$backend$name, "highs")
  expect_identical(fit$assignment$grade, c(1L, 2L, 3L))
  expect_identical(fit$summary$grade_count, 3L)
  expect_equal(fit$objective$raw, -2.68, tolerance = 1e-8)
  expect_equal(fit$objective$value, -5.36, tolerance = 1e-8)
  expect_equal(fit$objective$value, 2 * fit$objective$raw, tolerance = 1e-8)
})

test_that("gp_grade() requires the requested lambda to be on the active grid", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")$pairwise

  expect_error(
    gp_grade(pairwise, lambda = 0.333),
    regexp = "Available values",
    class = "gradepath_grade_error"
  )
})

test_that("gp_grade_path() uses the operational default grid instead of 101 solves", {
  skip_if_not_installed("digest")
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")$pairwise

  path <- gp_grade_path(pairwise)

  expect_s3_class(path, "gp_grade_path")
  expect_identical(path$backend$name, "highs")
  expect_identical(path$provenance$warm_start_strategy, "none")
  expect_identical(path$lambda_grid, .gp_grade_operational_default_grid())
  expect_lt(length(path$lambda_grid), length(.gp_control_default_lambda_grid))
  expect_true(0.25 %in% path$lambda_grid)
  expect_true(1 %in% path$lambda_grid)
  expect_length(path$fits, length(path$lambda_grid))
  expect_identical(path$summary$lambda, path$lambda_grid)
  expect_identical(
    path$summary$grade_count,
    vapply(path$fits, function(fit) fit$summary$grade_count, integer(1))
  )
  expect_identical(path$selection$selected_lambda, 0.25)
  expect_identical(path$selection$endpoint_lambda, 1)
  expect_identical(
    as.numeric(path$provenance$solve_order),
    c(1, rev(setdiff(path$lambda_grid, 1)))
  )
  expect_true(all(vapply(
    path$fits,
    function(fit) identical(fit$backend$warmstart, "none") &&
      is.na(fit$backend$warm_start_from_lambda) &&
      identical(fit$backend$warm_start_used, FALSE),
    logical(1)
  )))
})

test_that("gp_grade_path() stores fits in ascending grid order after endpoint-first solving", {
  skip_if_not_installed("digest")
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")$pairwise
  grid <- c(0, 0.25, 0.5, 1)

  path <- gp_grade_path(pairwise, lambda_grid = grid)

  expect_identical(path$lambda_grid, grid)
  expect_identical(
    vapply(path$fits, function(fit) fit$lambda, numeric(1)),
    grid
  )
  expect_identical(path$summary$lambda, grid)
  expect_identical(as.numeric(path$provenance$solve_order), c(1, 0.5, 0.25, 0))
})

test_that("gp_grade_path() validates selected lambda and endpoint membership", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")$pairwise

  expect_error(
    gp_grade_path(pairwise, lambda_grid = c(0, 0.5, 1)),
    regexp = "parity anchor 0.25",
    class = "gp_control_error"
  )
  expect_error(
    gp_grade_path(pairwise, lambda_grid = c(0, 0.25, 0.5), selected_lambda = 0.25),
    regexp = "parity anchor 1",
    class = "gp_control_error"
  )
  expect_error(
    gp_grade_path(pairwise, lambda_grid = c(0, 0.25, 0.5, 1), selected_lambda = 0.75),
    regexp = "selected_lambda",
    class = "gradepath_grade_error"
  )
})

test_that("backend adapters reject mismatched control backends", {
  fixture <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")

  expect_error(
    gp_open_solve(
      fixture$problem,
      backend = "highs",
      control = gp_control(backend = "glpk")
    ),
    regexp = "must match",
    class = "gradepath_backend_error"
  )
  expect_error(
    gp_gurobi_robust(
      fixture$problem,
      control = gp_control(backend = "highs"),
      force_path = "gurobi_cl_bigM"
    ),
    regexp = "must be `gurobi`",
    class = "gradepath_backend_error"
  )
})

test_that("gp_select_grade() returns a stored fit and never re-solves", {
  skip_if_not_installed("digest")
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")$pairwise
  path <- gp_grade_path(
    pairwise,
    lambda_grid = c(0, 0.25, 0.5, 1),
    selected_lambda = 0.5,
    selection_rule = "manual"
  )

  expect_identical(gp_select_grade(path), path$fits[[2L]])
  expect_identical(gp_select_grade(path, lambda = 0.5), path$fits[[3L]])
  expect_error(
    gp_select_grade(path, lambda = 0.75),
    regexp = "Available values",
    class = "gradepath_grade_error"
  )
})

test_that("gp_grade_path validator rejects ambiguous near-duplicate selections", {
  skip_if_not_installed("digest")
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")$pairwise
  path <- gp_grade_path(
    pairwise,
    lambda_grid = c(0, 0.25, 0.5, 1)
  )
  near_grid <- c(0, 0.25, 0.250000005, 1)
  fit_template <- path$fits[[2L]]
  near_fit <- new_gp_grade_fit(
    ids = fit_template$ids,
    lambda = near_grid[[3L]],
    assignment = fit_template$assignment,
    summary = fit_template$summary,
    objective = fit_template$objective,
    backend = fit_template$backend,
    control = fit_template$control,
    provenance = fit_template$provenance,
    warnings = fit_template$warnings
  )

  path$lambda_grid <- near_grid
  path$fits[[3L]] <- near_fit
  path$summary$lambda <- near_grid

  expect_error(
    validate_gp_grade_path(path),
    regexp = "exactly one member",
    class = "gradepath_error"
  )
})

test_that("Gurobi grade paths pass warm-start metadata when available", {
  skip_if_not_installed("digest")
  skip_if(!.gp_gurobi_cl_smoke_ok(), "gurobi_cl cannot solve a licensed smoke model")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "gurobi")$pairwise

  path <- gp_grade_path(pairwise, lambda_grid = c(0.25, 1))
  selected <- gp_select_grade(path)
  endpoint <- gp_select_grade(path, lambda = 1)

  expect_identical(endpoint$backend$warmstart, "none")
  expect_identical(endpoint$backend$warm_start_from_status, NA_character_)
  expect_equal(selected$backend$warm_start_from_lambda, 1)
  expect_identical(selected$backend$warm_start_from_status, endpoint$backend$status)
  expect_true(isTRUE(selected$backend$warm_start_used))
})

test_that("Gurobi warm-start metadata records non-optimal source status", {
  skip_if_not_installed("digest")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "gurobi")$pairwise
  ids <- pairwise$ids
  grid <- c(0.25, 0.5, 1)
  solve_calls <- list()

  mock_result <- function(lambda, status, warm_start) {
    assignment <- data.frame(
      id = ids,
      grade = c(1L, 2L, 3L),
      stringsAsFactors = FALSE
    )
    list(
      primal = rep(lambda, 12),
      objval = -lambda,
      status = status,
      solver_meta = list(
        backend = "gurobi",
        path = "mock",
        encoding = "indicator",
        objbound = -lambda,
        mipgap = if (identical(status, "gap_reached")) 0.01 else 0,
        runtime = 0,
        warmstart = if (is.null(warm_start)) "none" else "vector",
        problem_hash = paste0("problem-hash-", lambda),
        problem_signature = "problem-signature"
      ),
      validation = list(
        assignment = assignment,
        canon = -2 * lambda
      )
    )
  }

  testthat::local_mocked_bindings(
    .gp_grade_solve_backend = function(pairwise,
                                       lambda,
                                       control,
                                       matrix_opt,
                                       warm_start = NULL,
                                       acceptance_mode = FALSE) {
      solve_calls[[length(solve_calls) + 1L]] <<- list(
        lambda = lambda,
        warm_start = warm_start
      )
      status <- if (identical(lambda, 1)) "gap_reached" else "optimal"
      mock_result(lambda, status, warm_start)
    },
    .package = "gradepath"
  )

  path <- gp_grade_path(pairwise, lambda_grid = grid)
  endpoint <- gp_select_grade(path, lambda = 1)
  middle <- gp_select_grade(path, lambda = 0.5)
  selected <- gp_select_grade(path, lambda = 0.25)

  expect_identical(
    vapply(solve_calls, `[[`, numeric(1), "lambda"),
    c(1, 0.5, 0.25)
  )
  expect_null(solve_calls[[1L]]$warm_start)
  expect_identical(solve_calls[[2L]]$warm_start, rep(1, 12))
  expect_identical(solve_calls[[3L]]$warm_start, rep(0.5, 12))

  expect_false(endpoint$backend$warm_start_used)
  expect_true(is.na(endpoint$backend$warm_start_from_lambda))
  expect_identical(endpoint$backend$warm_start_from_status, NA_character_)

  expect_true(middle$backend$warm_start_used)
  expect_equal(middle$backend$warm_start_from_lambda, 1)
  expect_identical(middle$backend$warm_start_from_status, "gap_reached")
  expect_false(middle$backend$warm_start_from_acceptance_ready)
  expect_identical(
    middle$provenance$warm_start_from_status,
    "gap_reached"
  )
  expect_false(middle$provenance$warm_start_from_acceptance_ready)

  expect_true(selected$backend$warm_start_used)
  expect_equal(selected$backend$warm_start_from_lambda, 0.5)
  expect_identical(selected$backend$warm_start_from_status, "optimal")
  expect_true(selected$backend$warm_start_from_acceptance_ready)
  expect_identical(
    selected$provenance$warm_start_from_status,
    "optimal"
  )
  expect_true(selected$provenance$warm_start_from_acceptance_ready)
  expect_identical(path$provenance$warm_start_nonoptimal_policy, "allowed_labeled")

  expect_identical(
    path$provenance$warm_start_sources$warm_start_from_status,
    c("optimal", "gap_reached", NA_character_)
  )
  expect_identical(
    path$provenance$warm_start_sources$warm_start_from_acceptance_ready,
    c(TRUE, FALSE, NA)
  )
  expect_identical(
    path$provenance$warm_start_sources$warm_start_used,
    c(TRUE, TRUE, FALSE)
  )
})

# --- CCR-01: acceptance_mode arg-threading ----------------------------------
# These tests confirm the opt-in solver acceptance mode is reachable end-to-end
# through the public grade verbs WITHOUT changing default behaviour. They mock
# `gp_gurobi_robust()` (no licensed solver required) and observe whether
# `gp_grade_path()` / `gp_grade()` thread `acceptance_mode` down to it and
# record it back on the fit's `backend$acceptance_mode`. The mock mirrors the
# real honesty policy: with `acceptance_mode = FALSE` the first (gap_reached)
# path is kept; with `TRUE` the cascade proves `optimal`. The actual
# path1=gap_reached -> path2=optimal cascade inside `gp_gurobi_robust()` is
# covered directly in test-gurobi-robust.R.

ccr01_mock_robust_result <- function(ids, acceptance_mode) {
  status <- if (isTRUE(acceptance_mode)) "optimal" else "gap_reached"
  path <- if (isTRUE(acceptance_mode)) "subprocess_R_binding" else "in_process"
  list(
    primal = c(0, 1, 1, 0, 0, 1, 0, 0, 0, 3, 2, 1),
    objval = -1,
    status = status,
    solver_meta = list(
      backend = "gurobi",
      path = path,
      encoding = "indicator",
      objbound = -1,
      mipgap = if (isTRUE(acceptance_mode)) 0 else 0.01,
      runtime = 0,
      acceptance_mode = isTRUE(acceptance_mode),
      warmstart = "none",
      problem_hash = "problem-hash",
      problem_signature = "problem-signature"
    ),
    validation = list(
      assignment = data.frame(
        id = ids,
        grade = c(1L, 2L, 3L),
        stringsAsFactors = FALSE
      ),
      canon = -2
    )
  )
}

test_that("gp_grade_path() threads acceptance_mode to gp_gurobi_robust and records it", {
  pairwise <- step45_pairwise(step45_strict3_matrix(), backend = "gurobi")
  ctl <- gp_control(backend = "gurobi", time_limit = 10, mip_gap = 0.05)
  seen <- new.env(parent = emptyenv())
  seen$acceptance_mode <- logical()

  testthat::local_mocked_bindings(
    gp_gurobi_robust = function(problem,
                                control = NULL,
                                warm_start = NULL,
                                matrix_opt = NULL,
                                force_path = NULL,
                                acceptance_mode = FALSE) {
      seen$acceptance_mode <- c(seen$acceptance_mode, isTRUE(acceptance_mode))
      ccr01_mock_robust_result(pairwise$ids, acceptance_mode)
    },
    .package = "gradepath"
  )

  # Default: invocation-fallback only -> the arg arrives FALSE and the gap
  # solve is reported honestly with `gap_reached`.
  default_path <- gp_grade_path(
    pairwise,
    lambda_grid = c(0.25, 1),
    control = ctl
  )
  default_fit <- gp_select_grade(default_path, lambda = 1)
  expect_true(all(seen$acceptance_mode == FALSE))
  expect_false(default_fit$backend$acceptance_mode)
  expect_identical(default_fit$backend$status, "gap_reached")

  # Opt-in: the arg arrives TRUE; the solver proves optimal and the flag is
  # recorded back on every fit.
  seen$acceptance_mode <- logical()
  accepted_path <- gp_grade_path(
    pairwise,
    lambda_grid = c(0.25, 1),
    control = ctl,
    acceptance_mode = TRUE
  )
  accepted_fit <- gp_select_grade(accepted_path, lambda = 1)
  expect_true(all(seen$acceptance_mode == TRUE))
  expect_true(accepted_fit$backend$acceptance_mode)
  expect_identical(accepted_fit$backend$status, "optimal")
  expect_true(all(vapply(
    accepted_path$fits,
    function(fit) isTRUE(fit$backend$acceptance_mode),
    logical(1)
  )))
})

test_that("gp_grade_path() validates acceptance_mode as a scalar logical", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")$pairwise
  expect_error(
    gp_grade_path(pairwise, acceptance_mode = "yes"),
    regexp = "acceptance_mode",
    class = "gradepath_error"
  )
})

test_that("gp_grade() threads acceptance_mode to the backend solver (default FALSE)", {
  pairwise <- step45_pairwise(step45_strict3_matrix(), backend = "gurobi")
  ctl <- gp_control(backend = "gurobi", time_limit = 10)
  seen <- new.env(parent = emptyenv())

  mock_result <- function(acceptance_mode) {
    list(
      primal = c(0, 1, 1, 0, 0, 1, 0, 0, 0, 3, 2, 1),
      objval = -1,
      status = "optimal",
      solver_meta = list(
        backend = "gurobi",
        path = "mock",
        encoding = "indicator",
        objbound = -1,
        mipgap = 0,
        runtime = 0,
        acceptance_mode = isTRUE(acceptance_mode),
        warmstart = "none",
        problem_hash = "problem-hash",
        problem_signature = "problem-signature"
      ),
      validation = list(
        assignment = data.frame(
          id = pairwise$ids,
          grade = c(1L, 2L, 3L),
          stringsAsFactors = FALSE
        ),
        canon = -2
      )
    )
  }

  testthat::local_mocked_bindings(
    .gp_grade_solve_backend = function(pairwise,
                                       lambda,
                                       control,
                                       matrix_opt,
                                       warm_start = NULL,
                                       acceptance_mode = FALSE) {
      seen$acceptance_mode <- acceptance_mode
      mock_result(acceptance_mode)
    },
    .package = "gradepath"
  )

  default_fit <- gp_grade(pairwise, lambda = 1, control = ctl)
  expect_false(seen$acceptance_mode)
  expect_false(default_fit$backend$acceptance_mode)

  accepted_fit <- gp_grade(pairwise, lambda = 1, control = ctl, acceptance_mode = TRUE)
  expect_true(seen$acceptance_mode)
  expect_true(accepted_fit$backend$acceptance_mode)
})
