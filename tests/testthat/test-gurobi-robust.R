make_step43_pairwise <- function(P, ids = rownames(P), control = gp_control(time_limit = 10)) {
  new_gp_pairwise(
    ids = ids,
    matrix = P,
    power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = 0.5, zero_floor = 1e-7),
    source = list(
      stage = "posterior",
      rule = "outer_product",
      assumption = "one_level_independence"
    ),
    control = control,
    provenance = .gradepath_new_provenance(producer = "test")
  )
}

make_step43_problem <- function(lambda = 1) {
  ids <- c("a", "b", "c")
  P <- matrix(
    c(
      0.50, 0.95, 0.97,
      0.05, 0.50, 0.92,
      0.03, 0.08, 0.50
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  list(
    P = P,
    matrix_opt = {
      P0 <- P
      diag(P0) <- 0
      P0
    },
    problem = gp_grade_problem(make_step43_pairwise(P), lambda = lambda, encoding = "indicator")
  )
}

make_step43_gurobi_raw <- function(problem,
                                   status = "OPTIMAL",
                                   mipgap = 0,
                                   runtime = 0.1) {
  D <- matrix(
    c(
      0, 1, 1,
      0, 0, 1,
      0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  primal <- c(as.vector(t(D)), c(3, 2, 1))
  objval <- sum(problem$objective * primal)
  list(
    x = primal,
    objval = objval,
    status = status,
    objbound = objval,
    mipgap = mipgap,
    runtime = runtime
  )
}

make_fake_gurobi_cl_dir <- function() {
  dir <- tempfile("fake-gurobi-cl-")
  dir.create(dir)
  cli <- file.path(dir, "gurobi_cl")
  writeLines(c("#!/bin/sh", "exit 0"), cli, useBytes = TRUE)
  Sys.chmod(cli, mode = "0755")
  dir
}

test_that("Gurobi parameter mapping pins KRW tolerances and rejects wrong aliases", {
  ctl <- gp_control(
    time_limit = 12,
    mip_gap = 0.01,
    solver_options = list(Threads = 1)
  )
  params <- .gp_gurobi_params(control = ctl)

  expect_identical(params$OutputFlag, 0)
  expect_identical(params$FeasibilityTol, 1e-9)
  expect_identical(params$IntFeasTol, 1e-9)
  expect_identical(params$OptimalityTol, 1e-9)
  expect_equal(params$TimeLimit, 12)
  expect_equal(params$MIPGap, 0.01)
  expect_equal(params$Threads, 1)

  expect_error(
    .gp_gurobi_params(control = gp_control(solver_options = list(time_limit = 1, TimeLimit = 2))),
    regexp = "time_limit.*TimeLimit",
    class = "gp_control_error"
  )
  expect_error(
    .gp_gurobi_params(control = gp_control(solver_options = list(mip_gap = 0.01, MIPGap = 0.02))),
    regexp = "mip_gap.*MIPGap",
    class = "gp_control_error"
  )
  expect_error(
    .gp_gurobi_params(control = gp_control(solver_options = list(max_time = 2))),
    regexp = "Gurobi",
    class = "gp_control_error"
  )
  expect_error(
    .gp_gurobi_params(control = gp_control(solver_options = list(mip_rel_gap = 0.02))),
    regexp = "Gurobi",
    class = "gp_control_error"
  )
})

test_that("indicator model uses native genconind constraints, not big-M rows", {
  skip_if_not_installed("Matrix")
  fixture <- make_step43_problem(lambda = 1)
  problem <- fixture$problem

  model <- .gp_gurobi_model_indicator(problem)

  expect_s4_class(model$A, "dgCMatrix")
  expect_identical(dim(model$A), c(0L, problem$constraint_ncol))
  expect_length(model$genconind, 2L * problem$n_units * (problem$n_units - 1L))
  expect_null(model$constrnames)
  expect_identical(model$sense, character())
  expect_equal(model$ub[problem$diagonal_positions], c(0, 0, 0))
  expect_identical(model$varnames, .gp_gurobi_varnames(problem))

  first_one <- model$genconind[[1L]]
  first_zero <- model$genconind[[problem$n_units * (problem$n_units - 1L) + 1L]]
  expect_identical(first_one$binvar, problem$d_indices[1, 2])
  expect_identical(first_one$binval, 1L)
  expect_identical(first_one$sense, ">")
  expect_equal(first_one$rhs, 1)
  expect_identical(first_zero$binvar, problem$d_indices[1, 2])
  expect_identical(first_zero$binval, 0L)
  expect_identical(first_zero$sense, ">")
  expect_equal(first_zero$rhs, 0)
})

test_that("big-M Gurobi model translates the canonical rows and keeps the same hash", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("digest")
  fixture <- make_step43_problem(lambda = 0.25)
  problem_indicator <- fixture$problem
  problem_bigM <- problem_indicator
  problem_bigM$encoding <- "bigM"
  problem_bigM <- validate_gp_grade_problem(problem_bigM)

  model <- .gp_gurobi_model_bigM(problem_bigM)

  expect_s4_class(model$A, "dgCMatrix")
  expect_identical(dim(model$A), c(problem_bigM$constraint_nrow, problem_bigM$constraint_ncol))
  expect_identical(model$sense, rep("<", problem_bigM$constraint_nrow))
  expect_identical(model$rhs, problem_bigM$rhs)
  expect_null(model$genconind)
  expect_identical(gp_problem_hash(problem_indicator), gp_problem_hash(problem_bigM))
})

test_that("LP writers preserve indicator and big-M encodings explicitly", {
  fixture <- make_step43_problem(lambda = 1)
  problem <- fixture$problem
  ind_file <- tempfile(fileext = ".lp")
  bigm_file <- tempfile(fileext = ".lp")

  .gp_gurobi_write_lp_indicator(problem, ind_file)
  .gp_gurobi_write_lp_bigM(problem, bigm_file)
  ind_lines <- readLines(ind_file, warn = FALSE)
  bigm_lines <- readLines(bigm_file, warn = FALSE)

  expect_true(any(grepl("->", ind_lines, fixed = TRUE)))
  expect_true(any(grepl("ind_D_1_2_one", ind_lines, fixed = TRUE)))
  expect_true(any(grepl("0 <= D_1_1 <= 0", ind_lines, fixed = TRUE)))
  expect_false(any(grepl("->", bigm_lines, fixed = TRUE)))
  expect_true(any(grepl("bigM_1:", bigm_lines, fixed = TRUE)))
  expect_true(any(grepl("Binaries", ind_lines, fixed = TRUE)))
  expect_true(any(grepl("Generals", bigm_lines, fixed = TRUE)))
})

test_that("solution parsers rebuild primals by variable name", {
  skip_if_not_installed("jsonlite")
  fixture <- make_step43_problem(lambda = 1)
  varnames <- .gp_gurobi_varnames(fixture$problem)

  json_file <- tempfile(fileext = ".json")
  jsonlite::write_json(
    list(
      SolutionInfo = list(Status = 2, ObjVal = -1.25, ObjBound = -1.25, MIPGap = 0, Runtime = 0.1),
      Vars = list(
        list(VarName = varnames[[3]], X = 1),
        list(VarName = varnames[[1]], X = 0),
        list(VarName = varnames[[10]], X = 3)
      )
    ),
    json_file,
    auto_unbox = TRUE
  )
  parsed_json <- .gp_gurobi_parse_json(json_file, varnames)
  expect_equal(parsed_json$x[[3]], 1)
  expect_equal(parsed_json$x[[2]], 0)
  expect_equal(parsed_json$x[[10]], 3)
  expect_equal(parsed_json$objval, -1.25)
  expect_identical(parsed_json$status, "OPTIMAL")

  sol_file <- tempfile(fileext = ".sol")
  writeLines(
    c(
      "# Solution for model obj",
      "# Objective value = -2.5",
      paste(varnames[[11]], 2),
      paste(varnames[[2]], 1)
    ),
    sol_file
  )
  parsed_sol <- .gp_gurobi_parse_sol(sol_file, varnames)
  expect_equal(parsed_sol$x[[11]], 2)
  expect_equal(parsed_sol$x[[2]], 1)
  expect_equal(parsed_sol$x[[3]], 0)
  expect_equal(parsed_sol$objval, -2.5)
})

test_that("base JSON parser preserves structured CLI metadata without jsonlite", {
  fixture <- make_step43_problem(lambda = 1)
  varnames <- .gp_gurobi_varnames(fixture$problem)
  json_file <- tempfile(fileext = ".json")
  writeLines(
    paste0(
      '{ "SolutionInfo": { "Status": 2, "Runtime": 0.25, ',
      '"ObjVal": -2.68, "ObjBound": -2.68, "MIPGap": 0 }, ',
      '"Vars": [ { "VarName": "', varnames[[2]], '", "X": 1 }, ',
      '{ "VarName": "', varnames[[10]], '", "X": 3 } ] }'
    ),
    json_file
  )

  sol <- .gp_gurobi_parse_json_base(json_file)

  expect_equal(sol$SolutionInfo$Status, 2)
  expect_equal(sol$SolutionInfo$Runtime, 0.25)
  expect_equal(sol$SolutionInfo$ObjVal, -2.68)
  expect_equal(sol$SolutionInfo$ObjBound, -2.68)
  expect_equal(sol$SolutionInfo$MIPGap, 0)
  expect_equal(sol$Vars[[1]]$VarName, varnames[[2]])
  expect_equal(sol$Vars[[1]]$X, 1)
  expect_equal(sol$Vars[[2]]$VarName, varnames[[10]])
  expect_equal(sol$Vars[[2]]$X, 3)
})

test_that("normalized Gurobi results validate bounds, grades, hash, and 2x objective", {
  skip_if_not_installed("digest")
  fixture <- make_step43_problem(lambda = 1)
  problem <- fixture$problem
  D <- matrix(
    c(
      0, 1, 1,
      0, 0, 1,
      0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  primal <- c(as.vector(t(D)), c(3, 2, 1))
  objval <- sum(problem$objective * primal)

  out <- .gp_gurobi_normalize_result(
    result = list(
      x = primal,
      objval = objval,
      status = "OPTIMAL",
      objbound = objval,
      mipgap = 0,
      runtime = 0
    ),
    problem = problem,
    params = .gp_gurobi_params(control = gp_control(time_limit = 10)),
    path = "gurobi_cl_bigM",
    encoding = "bigM",
    matrix_opt = fixture$matrix_opt
  )

  expect_s3_class(out, "gp_gurobi_result")
  expect_identical(out$status, "optimal")
  expect_identical(out$validation$grades, c(1L, 2L, 3L))
  expect_identical(out$validation$problem_hash, gp_problem_hash(problem))
  expect_identical(
    out$validation$problem_signature,
    gp_problem_feasible_signature(problem)
  )
  expect_identical(
    out$solver_meta$problem_signature,
    gp_problem_feasible_signature(problem)
  )
  expect_equal(out$validation$canon, 2 * out$objval)
})

test_that("Gurobi normalization classifies solver correctness failures", {
  fixture <- make_step43_problem(lambda = 1)
  problem <- fixture$problem
  raw <- make_step43_gurobi_raw(problem, status = "OPTIMAL")
  params <- .gp_gurobi_params(control = gp_control(time_limit = 10))

  objective_error <- tryCatch(
    .gp_gurobi_normalize_result(
      result = modifyList(raw, list(objval = raw$objval + 1)),
      problem = problem,
      params = params,
      path = "gurobi_cl_bigM",
      encoding = "bigM",
      matrix_opt = fixture$matrix_opt
    ),
    error = function(e) e
  )
  expect_s3_class(objective_error, "gradepath_solver_objective_mismatch")
  expect_identical(
    .gp_status_from_condition(objective_error),
    "SOLVER_OBJECTIVE_MISMATCH"
  )
  expect_false(identical(
    .gp_status_from_condition(objective_error),
    "SOLVER_BACKEND_UNAVAILABLE"
  ))

  invalid_output <- tryCatch(
    .gp_gurobi_normalize_result(
      result = modifyList(raw, list(x = raw$x[-1])),
      problem = problem,
      params = params,
      path = "gurobi_cl_bigM",
      encoding = "bigM",
      matrix_opt = fixture$matrix_opt
    ),
    error = function(e) e
  )
  expect_s3_class(invalid_output, "gradepath_solver_output_invalid")
  expect_identical(
    .gp_status_from_condition(invalid_output),
    "SOLVER_OUTPUT_INVALID"
  )

  infeasible <- tryCatch(
    .gp_gurobi_normalize_result(
      result = modifyList(raw, list(status = "INFEASIBLE")),
      problem = problem,
      params = params,
      path = "gurobi_cl_bigM",
      encoding = "bigM",
      matrix_opt = fixture$matrix_opt
    ),
    error = function(e) e
  )
  expect_s3_class(infeasible, "gradepath_solver_infeasible")
  expect_identical(.gp_status_from_condition(infeasible), "SOLVER_INFEASIBLE")
})

test_that(".gp_gurobi_status only demotes OPTIMAL to gap_reached inside the commanded band", {
  # Rationale: Gurobi reports OPTIMAL only after a solution is within its own
  # working tolerances. `.gp_gurobi_status()` adds a *conservative* demotion:
  # a strictly positive MIPGap that still sits inside a user-commanded `MIPGap`
  # limit is surfaced as `gap_reached` so acceptance-mode can prefer a strict
  # zero-gap path. The demotion deliberately fires ONLY inside the acceptable
  # band -- two corners must therefore stay `optimal`.

  # (a) existing acceptable-gap case: positive gap within the commanded limit
  # is demoted to `gap_reached`.
  expect_identical(
    .gp_gurobi_status("OPTIMAL", mipgap = 0.01, params = list(MIPGap = 0.05)),
    "gap_reached"
  )

  # (b) existing zero-gap case: a (near-)zero gap stays `optimal` even with a
  # commanded limit, because mipgap <= 1e-9 falls below the demotion threshold.
  expect_identical(
    .gp_gurobi_status("OPTIMAL", mipgap = 0, params = list(MIPGap = 0.05)),
    "optimal"
  )

  # (i) UNTESTED CORNER: OPTIMAL with mipgap > gap_limit -> stays `optimal`.
  # Gurobi only reports OPTIMAL within tolerance, so the conservative
  # gap-demotion only fires INSIDE the acceptable band (mipgap <= limit). A gap
  # above the commanded limit is therefore left as `optimal`, never demoted.
  expect_identical(
    .gp_gurobi_status("OPTIMAL", mipgap = 0.10, params = list(MIPGap = 0.05)),
    "optimal"
  )

  # (ii) UNTESTED CORNER: OPTIMAL with a positive mipgap but NO commanded MIPGap
  # (gap_limit NULL) -> stays `optimal`. With no user-supplied limit the gap is
  # within Gurobi's own default tolerance, so there is nothing to demote.
  expect_identical(
    .gp_gurobi_status("OPTIMAL", mipgap = 0.10, params = list()),
    "optimal"
  )
})

test_that("acceptance mode continues past non-optimal Gurobi paths", {
  fixture <- make_step43_problem(lambda = 1)
  problem <- fixture$problem
  detect <- list(
    gurobi = TRUE,
    callr = TRUE,
    Matrix = TRUE,
    jsonlite = TRUE,
    gurobi_cl = FALSE,
    gurobi_cl_path = "",
    gurobi_cl_smoke = FALSE,
    binding_smoke = TRUE
  )
  gap_result <- make_step43_gurobi_raw(problem, status = "OPTIMAL", mipgap = 0.01)
  optimal_result <- make_step43_gurobi_raw(problem, status = "OPTIMAL", mipgap = 0)

  testthat::local_mocked_bindings(
    .gp_gurobi_detect = function(smoke = TRUE) detect,
    .gp_gurobi_path_available = function(path, detect) {
      path %in% c("in_process", "subprocess_R_binding")
    },
    .gp_gurobi_solve_inprocess = function(...) gap_result,
    .gp_gurobi_solve_subprocess = function(...) optimal_result,
    .package = "gradepath"
  )

  ctl <- gp_control(time_limit = 10, mip_gap = 0.05)
  ordinary <- gp_gurobi_robust(
    problem,
    control = ctl,
    matrix_opt = fixture$matrix_opt
  )
  expect_identical(ordinary$status, "gap_reached")
  expect_identical(ordinary$solver_meta$path, "in_process")
  expect_false(ordinary$solver_meta$acceptance_mode)
  expect_identical(ordinary$solver_meta$attempts$path, "in_process")
  expect_identical(ordinary$solver_meta$attempts$status, "gap_reached")
  expect_identical(
    ordinary$solver_meta$attempts$reason,
    "non_acceptance_ready_default_policy"
  )

  accepted <- gp_gurobi_robust(
    problem,
    control = ctl,
    matrix_opt = fixture$matrix_opt,
    acceptance_mode = TRUE
  )
  expect_identical(accepted$status, "optimal")
  expect_identical(accepted$solver_meta$path, "subprocess_R_binding")
  expect_true(accepted$solver_meta$acceptance_mode)
  expect_identical(
    accepted$solver_meta$attempts$path,
    c("in_process", "subprocess_R_binding")
  )
  expect_identical(
    accepted$solver_meta$attempts$status,
    c("gap_reached", "optimal")
  )
  expect_identical(
    accepted$solver_meta$attempts$reason[[1]],
    "non_acceptance_ready_in_acceptance_mode"
  )
})

test_that("vector force_path drives a genuine two-path acceptance cascade", {
  # CCR-19: force_path now accepts a character VECTOR, so a real two-leg
  # fall-through can be exercised: leg 1 returns gap_reached (not
  # acceptance-ready), the cascade continues, leg 2 returns optimal and is
  # returned. A scalar force_path stays a length-1 vector and is unaffected.
  fixture <- make_step43_problem(lambda = 1)
  problem <- fixture$problem
  detect <- list(
    gurobi = TRUE,
    callr = TRUE,
    Matrix = TRUE,
    jsonlite = TRUE,
    gurobi_cl = TRUE,
    gurobi_cl_path = "/fake/gurobi_cl",
    gurobi_cl_smoke = TRUE,
    binding_smoke = TRUE
  )
  gap_result <- make_step43_gurobi_raw(problem, status = "OPTIMAL", mipgap = 0.01)
  optimal_result <- make_step43_gurobi_raw(problem, status = "OPTIMAL", mipgap = 0)
  cli_calls <- new.env(parent = emptyenv())
  cli_calls$encoding <- character()

  testthat::local_mocked_bindings(
    .gp_gurobi_detect = function(smoke = TRUE) detect,
    .gp_gurobi_path_available = function(path, detect) TRUE,
    .gp_gurobi_cli_solve = function(problem, params, warm_start = NULL, encoding) {
      cli_calls$encoding <- c(cli_calls$encoding, encoding)
      if (identical(encoding, "indicator")) gap_result else optimal_result
    },
    .package = "gradepath"
  )

  out <- gp_gurobi_robust(
    problem,
    control = gp_control(time_limit = 10, mip_gap = 0.05),
    matrix_opt = fixture$matrix_opt,
    force_path = c("gurobi_cl_indicator", "gurobi_cl_bigM"),
    acceptance_mode = TRUE
  )

  # Path 1 (indicator) returned gap_reached; the cascade fell through to path 2
  # (bigM) which returned optimal and is the accepted result.
  expect_identical(out$status, "optimal")
  expect_identical(out$solver_meta$path, "gurobi_cl_bigM")
  expect_true(out$solver_meta$acceptance_mode)
  expect_identical(cli_calls$encoding, c("indicator", "bigM"))
  expect_identical(
    out$solver_meta$attempts$path,
    c("gurobi_cl_indicator", "gurobi_cl_bigM")
  )
  expect_identical(
    out$solver_meta$attempts$status,
    c("gap_reached", "optimal")
  )
  expect_identical(
    out$solver_meta$attempts$reason[[1]],
    "non_acceptance_ready_in_acceptance_mode"
  )
})

test_that("gurobi_cl availability is smoke-gated", {
  detect <- list(
    gurobi = FALSE,
    callr = FALSE,
    Matrix = FALSE,
    jsonlite = TRUE,
    gurobi_cl = TRUE,
    gurobi_cl_path = "/fake/gurobi_cl",
    gurobi_cl_smoke = FALSE,
    binding_smoke = FALSE
  )

  expect_false(.gp_gurobi_path_available("gurobi_cl_indicator", detect))
  expect_false(.gp_gurobi_path_available("gurobi_cl_bigM", detect))

  detect$gurobi_cl_smoke <- TRUE
  expect_true(.gp_gurobi_path_available("gurobi_cl_indicator", detect))
  expect_true(.gp_gurobi_path_available("gurobi_cl_bigM", detect))

  detect$gurobi_cl_smoke <- FALSE
  testthat::local_mocked_bindings(
    .gp_gurobi_detect = function(smoke = TRUE) detect,
    .package = "gradepath"
  )
  paths <- .gp_available_gurobi_paths(smoke = TRUE)
  expect_false(any(grepl("^gurobi_cl", paths)))
})

test_that("gurobi_cl detection records failed smoke tests as unavailable", {
  fake_dir <- make_fake_gurobi_cl_dir()
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(fake_dir, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  testthat::local_mocked_bindings(
    .gp_gurobi_cl_smoke_ok = function(timeout = 10) FALSE,
    .gp_gurobi_binding_smoke_ok = function(timeout = 10) FALSE,
    .package = "gradepath"
  )

  detect <- .gp_gurobi_detect(smoke = TRUE)
  expect_true(detect$gurobi_cl)
  expect_false(detect$gurobi_cl_smoke)
  expect_false(.gp_gurobi_path_available("gurobi_cl_bigM", detect))
})

test_that("gurobi_cl solve passes process timeout and reports timeout explicitly", {
  fake_dir <- make_fake_gurobi_cl_dir()
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(fake_dir, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  fixture <- make_step43_problem(lambda = 1)
  params <- .gp_gurobi_params(control = gp_control(time_limit = 2))
  seen <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    .gp_system2 = function(command, args, stdout, stderr, timeout) {
      seen$timeout <- timeout
      structure("mock process timed out", status = 124L, timeout = TRUE)
    },
    .package = "gradepath"
  )

  expect_error(
    .gp_gurobi_cli_solve(
      fixture$problem,
      params = params,
      encoding = "indicator"
    ),
    class = "gradepath_gurobi_cli_timeout"
  )
  expect_equal(seen$timeout, .gp_gurobi_cli_timeout(params))
})

test_that("gurobi_cl solve reports nonzero exits and missing JSON explicitly", {
  fake_dir <- make_fake_gurobi_cl_dir()
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(fake_dir, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  fixture <- make_step43_problem(lambda = 1)
  params <- .gp_gurobi_params(control = gp_control(time_limit = 2))

  testthat::local_mocked_bindings(
    .gp_system2 = function(command, args, stdout, stderr, timeout) {
      structure("mock license failure", status = 1L)
    },
    .package = "gradepath"
  )
  expect_error(
    .gp_gurobi_cli_solve(
      fixture$problem,
      params = params,
      encoding = "indicator"
    ),
    class = "gradepath_gurobi_cli_error"
  )

  testthat::local_mocked_bindings(
    .gp_system2 = function(command, args, stdout, stderr, timeout) {
      structure(character(), status = 0L)
    },
    .package = "gradepath"
  )
  expect_error(
    .gp_gurobi_cli_solve(
      fixture$problem,
      params = params,
      encoding = "bigM"
    ),
    class = "gradepath_gurobi_cli_missing_json"
  )
})

test_that("robust cascade records CLI timeout attempts explicitly", {
  fixture <- make_step43_problem(lambda = 1)
  problem <- fixture$problem
  detect <- list(
    gurobi = FALSE,
    callr = FALSE,
    Matrix = FALSE,
    jsonlite = TRUE,
    gurobi_cl = TRUE,
    gurobi_cl_path = "/fake/gurobi_cl",
    gurobi_cl_smoke = TRUE,
    binding_smoke = FALSE
  )
  calls <- new.env(parent = emptyenv())
  calls$encoding <- character()

  testthat::local_mocked_bindings(
    .gp_gurobi_detect = function(smoke = TRUE) detect,
    .gp_gurobi_cli_solve = function(problem, params, warm_start = NULL, encoding) {
      calls$encoding <- c(calls$encoding, encoding)
      if (identical(encoding, "indicator")) {
        .gp_status_abort(
          "SOLVER_TIME_LIMIT",
          "mock gurobi_cl timeout",
          class = c("gradepath_gurobi_cli_timeout", "gradepath_backend_error")
        )
      }
      make_step43_gurobi_raw(problem, status = "OPTIMAL")
    },
    .package = "gradepath"
  )

  out <- gp_gurobi_robust(
    problem,
    control = gp_control(time_limit = 2),
    matrix_opt = fixture$matrix_opt
  )
  expect_identical(out$status, "optimal")
  expect_identical(out$solver_meta$path, "gurobi_cl_bigM")
  expect_identical(calls$encoding, c("indicator", "bigM"))

  attempts <- out$solver_meta$attempts
  timeout_row <- attempts[attempts$path == "gurobi_cl_indicator", , drop = FALSE]
  expect_identical(timeout_row$status, "time_limit")
  expect_match(timeout_row$reason, "mock gurobi_cl timeout")
})

test_that("unavailable error points to Gurobi installation and open fallback", {
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
  expect_error(
    .gp_gurobi_abort_unavailable(empty_detect),
    regexp = "gp_control\\(backend = \"highs\"\\)",
    class = "gradepath_gurobi_unavailable"
  )
})

test_that("gurobi_cl indicator and big-M paths solve the strict 3-unit fixture equivalently", {
  skip_if_not_installed("digest")
  skip_if(Sys.which("gurobi_cl") == "", "gurobi_cl not on PATH")
  skip_if(!.gp_gurobi_cl_smoke_ok(), "gurobi_cl cannot solve a licensed smoke model")
  fixture <- make_step43_problem(lambda = 1)
  problem <- fixture$problem
  ctl <- gp_control(time_limit = 10)

  indicator <- gp_gurobi_robust(
    problem,
    control = ctl,
    matrix_opt = fixture$matrix_opt,
    force_path = "gurobi_cl_indicator"
  )
  bigm <- gp_gurobi_robust(
    problem,
    control = ctl,
    matrix_opt = fixture$matrix_opt,
    force_path = "gurobi_cl_bigM"
  )

  expect_identical(indicator$validation$problem_hash, bigm$validation$problem_hash)
  expect_identical(
    indicator$validation$problem_signature,
    bigm$validation$problem_signature
  )
  expect_identical(indicator$validation$grades, c(1L, 2L, 3L))
  expect_identical(bigm$validation$grades, c(1L, 2L, 3L))
  expect_equal(indicator$objval, bigm$objval, tolerance = 1e-9)
  expect_equal(indicator$validation$canon, 2 * indicator$objval)
  expect_equal(bigm$validation$canon, 2 * bigm$objval)
  expect_identical(indicator$solver_meta$encoding, "indicator")
  expect_identical(bigm$solver_meta$encoding, "bigM")
})
