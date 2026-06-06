test_that("producer-status contract maps solver statuses conservatively", {
  contract <- .gp_status_contract()

  row <- function(source_status) {
    contract[identical_na(contract$source_status, source_status), , drop = FALSE]
  }
  identical_na <- function(x, y) {
    if (is.na(y)) {
      is.na(x)
    } else {
      identical(x, y) | (!is.na(x) & x == y)
    }
  }

  expect_identical(row("optimal")$producer_status, "OK")
  expect_identical(row("OK")$producer_status, "OK")
  expect_identical(row("gap_reached")$producer_status, "SOLVER_GAP")
  expect_identical(row("suboptimal")$producer_status, "SOLVER_GAP")
  expect_identical(row("time_limit")$producer_status, "SOLVER_TIME_LIMIT")
  expect_identical(row("infeasible")$producer_status, "SOLVER_INFEASIBLE")
  expect_identical(row("output_invalid")$producer_status, "SOLVER_OUTPUT_INVALID")
  expect_identical(row("objective_mismatch")$producer_status, "SOLVER_OBJECTIVE_MISMATCH")
  expect_identical(row("canonical_mismatch")$producer_status, "SOLVER_CANONICAL_MISMATCH")
  expect_identical(row("unavailable")$producer_status, "SOLVER_BACKEND_UNAVAILABLE")
  expect_identical(row(NA_character_)$producer_status, "UNVERIFIED")

  expect_true(.gp_status_acceptance_ready("optimal"))
  expect_true(.gp_status_acceptance_ready("OK"))
  expect_false(.gp_status_acceptance_ready("APPROXIMATE_OK"))
  expect_true(.gp_status_requires_unverified("gap_reached"))
  expect_true(.gp_status_requires_unverified("time_limit"))

  non_ok <- setdiff(.gp_status_values, "OK")
  expect_false(any(vapply(non_ok, .gp_status_acceptance_ready, logical(1))))
})

test_that("selected grade status uses the selected solver status", {
  ids <- c("a", "b", "c")
  ctl <- gp_control(lambda_grid = c(0.25, 1), backend = "highs")
  assignment <- data.frame(
    id = ids,
    grade = c(1L, 2L, 3L),
    score = c(2L, 1L, 0L),
    stringsAsFactors = FALSE
  )
  fit <- new_gp_grade_fit(
    ids = ids,
    lambda = 0.25,
    assignment = assignment,
    summary = list(grade_count = 3L, status = "gap_reached", n_units = 3L),
    objective = list(value = -2, raw = -1, canonical = -2),
    backend = list(
      name = "gurobi",
      path = "gurobi_cl_indicator",
      encoding = "indicator",
      status = "gap_reached",
      objbound = -2.1,
      mipgap = 0.01,
      runtime = 1,
      warmstart = "none",
      warm_start_from_lambda = NA_real_,
      warm_start_used = FALSE,
      problem_hash = "toy"
    ),
    control = ctl,
    provenance = list(),
    warnings = "Solver returned status `gap_reached`."
  )

  expect_identical(.gp_selected_grade_solver_status(fit), "SOLVER_GAP")
  expect_identical(.gp_producer_status_from_selected_grade(fit), "SOLVER_GAP")
})

test_that("backend condition taxonomy preserves solver correctness failures", {
  condition_from <- function(expr) {
    tryCatch(
      expr,
      error = function(e) e
    )
  }

  unavailable <- condition_from(.gradepath_abort_backend_unavailable("missing backend"))
  expect_identical(.gp_status_from_condition(unavailable), "SOLVER_BACKEND_UNAVAILABLE")

  infeasible <- condition_from(.gradepath_abort_solver_infeasible("infeasible"))
  expect_identical(.gp_status_from_condition(infeasible), "SOLVER_INFEASIBLE")

  output_invalid <- condition_from(.gradepath_abort_solver_output_invalid("bad output"))
  expect_identical(.gp_status_from_condition(output_invalid), "SOLVER_OUTPUT_INVALID")

  objective <- condition_from(.gradepath_abort_solver_objective_mismatch("bad objective"))
  expect_identical(.gp_status_from_condition(objective), "SOLVER_OBJECTIVE_MISMATCH")
  expect_false(identical(
    .gp_status_from_condition(objective),
    "SOLVER_BACKEND_UNAVAILABLE"
  ))

  canonical <- condition_from(.gradepath_abort_solver_canonical_mismatch("bad canonical"))
  expect_identical(.gp_status_from_condition(canonical), "SOLVER_CANONICAL_MISMATCH")

  generic_backend <- condition_from(.gradepath_abort_backend("generic backend bug"))
  expect_identical(.gp_status_from_condition(generic_backend), "INTERNAL_ERROR")

  timeout <- condition_from(
    .gp_status_abort(
      "SOLVER_TIME_LIMIT",
      "timeout",
      class = c("gradepath_gurobi_cli_timeout", "gradepath_backend_error")
    )
  )
  expect_identical(.gp_status_from_condition(timeout), "SOLVER_TIME_LIMIT")
})

test_that("gp_check consumer boundary routes non-OK producers to UNVERIFIED", {
  reg <- data.frame(
    id = "toy",
    paper_value = "1",
    unit = "count",
    tolerance = "0",
    class = "exact",
    milestone = "M1",
    quantity = "toy target",
    stringsAsFactors = FALSE
  )
  statuses <- c(
    gap_reached = "SOLVER_GAP",
    time_limit = "SOLVER_TIME_LIMIT",
    infeasible = "SOLVER_INFEASIBLE",
    unavailable = "SOLVER_BACKEND_UNAVAILABLE",
    output_invalid = "SOLVER_OUTPUT_INVALID",
    objective_mismatch = "SOLVER_OBJECTIVE_MISMATCH",
    canonical_mismatch = "SOLVER_CANONICAL_MISMATCH",
    APPROXIMATE_OK = "APPROXIMATE_OK",
    CACHE_STALE = "CACHE_STALE",
    unknown_status = "UNVERIFIED"
  )
  for (status in names(statuses)) {
    check <- gp_check(
      "toy",
      replicated = 1,
      producer_status = status,
      registry = reg
    )
    expect_identical(check$status, "UNVERIFIED", info = status)
    expect_identical(check$reason, unname(statuses[[status]]), info = status)
    expect_identical(check$producer_status, unname(statuses[[status]]), info = status)
    expect_true(is.na(check$replicated), info = status)
  }

  expect_identical(
    gp_check("toy", replicated = 1, producer_status = "optimal", registry = reg)$status,
    "PASS"
  )
})
