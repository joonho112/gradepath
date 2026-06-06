make_step44_pairwise <- function(P, ids = rownames(P), backend = "highs") {
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
    control = gp_control(backend = backend),
    provenance = .gradepath_new_provenance(producer = "test")
  )
}

make_step44_strict_problem <- function(n = 3, lambda = 1) {
  ids <- paste0("u", seq_len(n))
  P <- matrix(0.5, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i < j) P[i, j] <- 0.9 + (0.01 * (j - i))
      if (i > j) P[i, j] <- 1 - P[j, i]
    }
  }
  diag(P) <- 0.5
  matrix_opt <- P
  diag(matrix_opt) <- 0
  list(
    P = P,
    matrix_opt = matrix_opt,
    problem = gp_grade_problem(make_step44_pairwise(P), lambda = lambda, encoding = "bigM")
  )
}

make_step44_seeded_problem <- function(n = 8, lambda = 0.25, seed = 1) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  score <- sort(stats::runif(n), decreasing = TRUE)
  P <- matrix(0.5, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) {
        P[i, j] <- stats::plogis(5 * (score[i] - score[j]))
      }
    }
  }
  diag(P) <- 0.5
  matrix_opt <- P
  diag(matrix_opt) <- 0
  list(
    P = P,
    matrix_opt = matrix_opt,
    problem = gp_grade_problem(make_step44_pairwise(P), lambda = lambda, encoding = "bigM")
  )
}

test_that("open backend availability reports only installed engines", {
  available <- .gp_available_open_backends()

  expect_true(all(available %in% .gp_open_backends))
  expect_setequal(.gp_open_backend_packages("highs"), c("highs", "Matrix"))
  for (backend in available) {
    expect_true(.gp_open_backend_available(backend))
  }
})

test_that("open backend control mapping is strict and backend-specific", {
  highs_controls <- .gp_open_solver_controls(
    solver_options = list(time_limit = 3, mip_gap = 0.05, threads = 1),
    backend = "highs"
  )
  expect_equal(highs_controls$time_limit, 3)
  expect_equal(highs_controls$mip_rel_gap, 0.05)
  expect_equal(highs_controls$threads, 1)
  expect_false("max_time" %in% names(highs_controls))

  highs_native <- .gp_open_solver_controls(
    solver_options = list(max_time = 4, mip_rel_gap = 0.02),
    backend = "highs"
  )
  expect_equal(highs_native$time_limit, 4)
  expect_equal(highs_native$mip_rel_gap, 0.02)

  if (requireNamespace("highs", quietly = TRUE)) {
    highs_runtime <- do.call(highs::highs_control, highs_controls)
    expect_equal(highs_runtime$time_limit, 3)
    expect_false("max_time" %in% names(highs_runtime))
  }

  symphony_controls <- .gp_open_solver_controls(
    solver_options = list(time_limit = 5, mip_gap = 0.01),
    backend = "symphony"
  )
  expect_equal(symphony_controls$time_limit, 5)
  expect_equal(symphony_controls$gap_limit, 0.01)

  glpk_controls <- .gp_open_solver_controls(
    solver_options = list(time_limit = 6),
    backend = "glpk"
  )
  expect_equal(glpk_controls$tm_limit, 6000)

  expect_error(
    .gp_open_solver_controls(solver_options = list(mip_gap = 0.1), backend = "glpk"),
    regexp = "does not expose",
    class = "gp_control_error"
  )
  expect_error(
    .gp_open_solver_controls(solver_options = list(time_limit = 1, max_time = 2), backend = "highs"),
    regexp = "multiple time",
    class = "gp_control_error"
  )
  expect_error(
    .gp_open_solver_controls(solver_options = list(mip_gap = 0.1, mip_rel_gap = 0.2), backend = "highs"),
    regexp = "multiple MIP-gap",
    class = "gp_control_error"
  )
})

test_that("ROI OP construction uses v2 bounds and canonical big-M rows", {
  skip_if_not_installed("ROI")
  skip_if_not_installed("slam")
  fixture <- make_step44_strict_problem(n = 3, lambda = 1)
  problem <- fixture$problem

  op <- .gp_roi_op(problem)
  bounds <- ROI::bounds(op)
  lower_full <- rep(0, problem$constraint_ncol)
  lower_full[bounds$lower$ind] <- bounds$lower$val
  upper_full <- rep(Inf, problem$constraint_ncol)
  upper_full[bounds$upper$ind] <- bounds$upper$val

  expect_s3_class(op, "OP")
  expect_equal(lower_full, problem$lower)
  expect_equal(upper_full, problem$upper)
  expect_equal(upper_full[problem$diagonal_positions], c(0, 0, 0))
  expect_equal(lower_full[problem$grade_indices], c(1, 1, 1))
  expect_equal(upper_full[problem$grade_indices], c(3, 3, 3))
})

test_that("open result normalization validates objective, grades, and hash", {
  skip_if_not_installed("digest")
  fixture <- make_step44_strict_problem(n = 3, lambda = 1)
  problem <- fixture$problem
  D <- matrix(c(0, 1, 1, 0, 0, 1, 0, 0, 0), 3, byrow = TRUE)
  primal <- c(as.vector(t(D)), c(3, 2, 1))
  objval <- sum(problem$objective * primal)

  out <- .gp_open_normalize_result(
    result = list(
      x = primal,
      objval = objval,
      status = "OPTIMAL",
      message = "Optimal",
      objbound = objval,
      mipgap = 0,
      runtime = 0
    ),
    problem = problem,
    params = list(),
    backend = "highs",
    path = "highs",
    matrix_opt = fixture$matrix_opt
  )

  expect_s3_class(out, "gp_open_backend_result")
  expect_identical(out$status, "optimal")
  expect_identical(out$validation$grades, c(1L, 2L, 3L))
  expect_identical(out$validation$problem_hash, gp_problem_hash(problem))
  expect_equal(out$validation$canon, 2 * out$objval)

  bad <- primal
  bad[1] <- 1
  expect_error(
    .gp_open_normalize_result(
      result = list(x = bad, objval = sum(problem$objective * bad), status = "OPTIMAL"),
      problem = problem,
      params = list(),
      backend = "highs",
      path = "highs"
    ),
    regexp = "bounds",
    class = "gradepath_error"
  )
  expect_error(
    .gp_open_normalize_result(
      result = list(x = primal[-1], objval = objval, status = "OPTIMAL"),
      problem = problem,
      params = list(),
      backend = "highs",
      path = "highs"
    ),
    regexp = "wrong length",
    class = "gradepath_error"
  )
  expect_error(
    .gp_open_normalize_result(
      result = list(x = primal, objval = objval + 1, status = "OPTIMAL"),
      problem = problem,
      params = list(),
      backend = "highs",
      path = "highs"
    ),
    regexp = "c'x",
    class = "gradepath_error"
  )
})

test_that("open MILP backends solve the strict chain identically", {
  skip_if_not_installed("digest")
  fixture <- make_step44_strict_problem(n = 3, lambda = 1)
  available <- .gp_available_open_backends()
  skip_if(length(available) == 0L, "no open MILP backend installed")

  results <- lapply(available, function(backend) {
    gp_open_solve(
      fixture$problem,
      backend,
      control = gp_control(backend = backend, time_limit = 10),
      matrix_opt = fixture$matrix_opt
    )
  })
  names(results) <- available

  for (backend in names(results)) {
    expect_identical(results[[backend]]$validation$grades, c(1L, 2L, 3L))
    expect_equal(results[[backend]]$validation$canon, 2 * results[[backend]]$objval)
    expect_identical(results[[backend]]$validation$problem_hash, gp_problem_hash(fixture$problem))
  }

  ref <- results[[1L]]
  for (backend in names(results)[-1L]) {
    expect_identical(results[[backend]]$validation$grades, ref$validation$grades)
    expect_equal(results[[backend]]$objval, ref$objval, tolerance = 1e-6)
  }
})

test_that("open MILP backends agree on a seeded antisymmetric battery", {
  skip_if_not_installed("digest")
  available <- .gp_available_open_backends()
  skip_if(length(available) < 2L, "need at least two open MILP backends")

  fixtures <- list(
    make_step44_seeded_problem(n = 8, lambda = 0, seed = 11),
    make_step44_seeded_problem(n = 9, lambda = 0.25, seed = 12),
    make_step44_seeded_problem(n = 10, lambda = 0.5, seed = 13),
    make_step44_seeded_problem(n = 8, lambda = 1, seed = 14)
  )

  for (fixture in fixtures) {
    results <- lapply(available, function(backend) {
      gp_open_solve(
        fixture$problem,
        backend,
        control = gp_control(backend = backend, time_limit = 15),
        matrix_opt = fixture$matrix_opt
      )
    })
    names(results) <- available
    ref <- results[[1L]]
    for (backend in names(results)[-1L]) {
      expect_identical(results[[backend]]$validation$problem_hash, ref$validation$problem_hash)
      expect_identical(results[[backend]]$validation$grades, ref$validation$grades)
      expect_equal(results[[backend]]$objval, ref$objval, tolerance = 1e-6)
    }
  }
})

test_that("Alabama DP is a lambda-zero pooled oracle only", {
  skip_if_not_installed("digest")
  fixture <- make_step44_seeded_problem(n = 8, lambda = 0, seed = 22)
  available <- .gp_available_open_backends()
  skip_if(length(available) == 0L, "no open MILP backend installed")

  dp <- gp_alabama_dp(fixture$problem, matrix_opt = fixture$matrix_opt)
  milp <- gp_open_solve(
    fixture$problem,
    available[[1L]],
    control = gp_control(backend = available[[1L]], time_limit = 10),
    matrix_opt = fixture$matrix_opt
  )

  expect_identical(dp$validation$grades, rep(1L, fixture$problem$n_units))
  expect_identical(milp$validation$grades, dp$validation$grades)
  expect_equal(dp$objval, 0)
  expect_equal(milp$objval, 0)
  expect_equal(dp$validation$canon, 0)

  lambda_positive <- make_step44_seeded_problem(n = 8, lambda = 0.25, seed = 22)
  expect_error(
    gp_alabama_dp(lambda_positive$problem, matrix_opt = lambda_positive$matrix_opt),
    regexp = "lambda = 0",
    class = "gradepath_error"
  )
})
