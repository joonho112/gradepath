test_that("public pairwise objects reject contaminated diagonals", {
  skip_if_not_installed("digest")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 0.25)$pairwise
  bad <- pairwise
  diag(bad$matrix) <- c(0.123, 0.456, 0.789)

  expect_error(
    validate_gp_pairwise(bad),
    regexp = "diagonal",
    class = "gradepath_error"
  )
})

test_that("grade problem projects any incoming pairwise diagonal out of the IP", {
  skip_if_not_installed("digest")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 0.25)$pairwise
  perturbed <- pairwise
  diag(perturbed$matrix) <- c(0.123, 0.456, 0.789)

  reference <- gp_grade_problem(pairwise, lambda = 0.25, encoding = "bigM")
  projected <- gp_grade_problem(perturbed, lambda = 0.25, encoding = "bigM")

  expect_identical(gp_problem_hash(projected), gp_problem_hash(reference))
  expect_equal(projected$lower[projected$diagonal_positions], rep(0, 3))
  expect_equal(projected$upper[projected$diagonal_positions], rep(0, 3))
  expect_equal(projected$objective[projected$diagonal_positions], rep(0, 3))
  expect_identical(
    projected$metadata$diagonal_projection,
    "upper_bound_zero_and_objective_zero"
  )
})

test_that("grade problem diagonal projection does not mask off-diagonal corruption", {
  skip_if_not_installed("digest")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 0.25)$pairwise
  bad <- pairwise
  diag(bad$matrix) <- c(0.123, 0.456, 0.789)
  bad$matrix[1, 2] <- 0.40

  expect_error(
    gp_grade_problem(bad, lambda = 0.25, encoding = "bigM"),
    regexp = "antisymmetry",
    class = "gradepath_error"
  )
})

test_that("solver results preserve projected zero diagonal decisions", {
  skip_if_not_installed("digest")
  fixture <- step45_problem(step45_strict3_matrix(), lambda = 1)
  results <- step45_solver_results(fixture$problem, fixture$matrix_opt)

  step45_skip_without_solver(results)

  for (name in names(results)) {
    out <- results[[name]]
    n_units <- nrow(out$validation$decision_matrix)
    expect_true(out$validation$diagonal_zero)
    expect_lte(max(abs(out$primal[gp_diag_positions(n_units)])), 1e-6)
    expect_equal(out$validation$canon, 2 * out$objval, tolerance = 1e-6)
    expect_identical(out$validation$grades, c(1L, 2L, 3L))
  }
})

test_that("solving is invariant to pairwise diagonal perturbations", {
  skip_if_not_installed("digest")
  clean <- step45_problem(step45_strict3_matrix(), lambda = 1)
  perturbed_pairwise <- clean$pairwise
  diag(perturbed_pairwise$matrix) <- c(0.123, 0.456, 0.789)
  perturbed_problem <- gp_grade_problem(
    perturbed_pairwise,
    lambda = 1,
    encoding = "bigM"
  )

  clean_results <- step45_solver_results(clean$problem, clean$matrix_opt)
  perturbed_results <- step45_solver_results(perturbed_problem, clean$matrix_opt)
  step45_skip_without_solver(clean_results)
  step45_skip_without_solver(perturbed_results)

  expect_setequal(names(perturbed_results), names(clean_results))
  for (name in names(clean_results)) {
    ref <- clean_results[[name]]
    other <- perturbed_results[[name]]
    expect_identical(other$validation$problem_hash, ref$validation$problem_hash)
    expect_identical(other$validation$grades, ref$validation$grades)
    expect_equal(other$validation$canon, ref$validation$canon, tolerance = 1e-8)
    expect_equal(other$objval, ref$objval, tolerance = 1e-8)
    expect_lte(max(abs(other$primal[gp_diag_positions(3L)])), 1e-6)
  }
})

test_that("diagonal contamination of problem bounds/objective is rejected before solve", {
  skip_if_not_installed("digest")
  fixture <- step45_problem(step45_strict3_matrix(), lambda = 1)

  bad_objective <- fixture$problem
  bad_objective$objective[bad_objective$diagonal_positions[1L]] <- 99
  expect_error(
    gp_open_solve(bad_objective, backend = "highs"),
    regexp = "diagonal",
    class = "gradepath_error"
  )
  expect_error(
    gp_gurobi_robust(bad_objective, force_path = "gurobi_cl_bigM"),
    regexp = "diagonal",
    class = "gradepath_error"
  )

  bad_upper <- fixture$problem
  bad_upper$upper[bad_upper$diagonal_positions[1L]] <- 1
  expect_error(
    gp_open_solve(bad_upper, backend = "highs"),
    regexp = "diagonal",
    class = "gradepath_error"
  )
  expect_error(
    gp_gurobi_robust(bad_upper, force_path = "gurobi_cl_bigM"),
    regexp = "diagonal",
    class = "gradepath_error"
  )

  bad_matrix_opt <- fixture$matrix_opt
  diag(bad_matrix_opt) <- c(1, 2, 3)
  if ("highs" %in% .gp_available_open_backends()) {
    expect_error(
      gp_open_solve(fixture$problem, backend = "highs", matrix_opt = bad_matrix_opt),
      regexp = "projected zero",
      class = "gradepath_backend_error"
    )
  }
  expect_error(
    gp_gurobi_robust(
      fixture$problem,
      force_path = "gurobi_cl_bigM",
      matrix_opt = bad_matrix_opt
    ),
    regexp = "projected zero",
    class = "gradepath_backend_error"
  )
})
