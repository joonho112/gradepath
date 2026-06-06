test_that("installed solver paths solve the same strict n=3 canonical problem", {
  skip_if_not_installed("digest")
  fixture <- step45_problem(step45_strict3_matrix(), lambda = 1)
  problem <- fixture$problem
  results <- step45_solver_results(problem, fixture$matrix_opt)

  step45_skip_without_solver(results, min_results = 2L)

  for (name in names(results)) {
    res <- results[[name]]
    expect_length(res$primal, 3L^2 + 3L)
    expect_identical(res$validation$grades, c(1L, 2L, 3L))
    expect_equal(res$objval, -2.68, tolerance = 1e-8)
    expect_equal(res$validation$canon, -5.36, tolerance = 1e-8)
    expect_equal(res$validation$canon, 2 * res$objval, tolerance = 1e-8)
    expect_true(res$validation$diagonal_zero)
    expect_lte(max(abs(res$primal[gp_diag_positions(3L)])), 1e-6)
    expect_identical(res$validation$problem_hash, gp_problem_hash(problem))
    expect_identical(res$solver_meta$problem_hash, gp_problem_hash(problem))
    expect_identical(
      res$validation$problem_signature,
      gp_problem_feasible_signature(problem)
    )
    expect_identical(
      res$solver_meta$problem_signature,
      gp_problem_feasible_signature(problem)
    )
  }

  ref <- results[[1L]]
  for (name in names(results)[-1L]) {
    other <- results[[name]]
    expect_silent(gp_assert_same_problem(ref, other, n = 3L))
    expect_silent(gp_assert_same_answer(ref, other))
    expect_equal(other$validation$decision_matrix, ref$validation$decision_matrix,
                 tolerance = 1e-8)
    expect_equal(other$objval, ref$objval, tolerance = 1e-6)
  }
})

test_that("same-problem and same-answer assertions reject drift", {
  skip_if_not_installed("digest")
  fixture <- step45_problem(step45_strict3_matrix(), lambda = 1)
  D <- matrix(
    c(0, 1, 1,
      0, 0, 1,
      0, 0, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(fixture$problem$ids, fixture$problem$ids)
  )
  primal <- numeric(fixture$problem$constraint_ncol)
  primal[seq_len(9)] <- as.vector(t(D))
  primal[fixture$problem$grade_indices] <- c(3, 2, 1)
  raw <- sum(fixture$problem$objective * primal)
  ref <- .gp_open_normalize_result(
    result = list(
      x = primal,
      objval = raw,
      status = "OPTIMAL",
      message = "fixture",
      objbound = raw,
      mipgap = 0,
      runtime = 0
    ),
    problem = fixture$problem,
    params = list(),
    backend = "highs",
    path = "fixture",
    encoding = "bigM",
    matrix_opt = fixture$matrix_opt
  )
  same <- ref

  expect_silent(gp_assert_same_problem(ref, same, n = 3L))
  expect_silent(gp_assert_same_answer(ref, same))

  bad_problem <- same
  bad_problem$validation$problem_hash <- paste0(bad_problem$validation$problem_hash, "x")
  expect_error(
    gp_assert_same_problem(ref, bad_problem, n = 3L),
    regexp = "same canonical",
    class = "gradepath_backend_error"
  )

  bad_signature <- same
  bad_signature$validation$problem_signature <- paste0(
    bad_signature$validation$problem_signature,
    "x"
  )
  expect_error(
    gp_assert_same_problem(ref, bad_signature, n = 3L),
    regexp = "same feasible",
    class = "gradepath_backend_error"
  )

  mutated_problem <- fixture$problem
  mutated_problem$constraint_x[[1L]] <- mutated_problem$constraint_x[[1L]] + 0.25
  mutated_problem <- validate_gp_grade_problem(mutated_problem)
  mutated <- .gp_open_normalize_result(
    result = list(
      x = primal,
      objval = raw,
      status = "OPTIMAL",
      message = "fixture",
      objbound = raw,
      mipgap = 0,
      runtime = 0
    ),
    problem = mutated_problem,
    params = list(),
    backend = "highs",
    path = "fixture",
    encoding = "bigM",
    matrix_opt = fixture$matrix_opt
  )
  expect_false(identical(
    mutated$validation$problem_signature,
    ref$validation$problem_signature
  ))
  expect_error(
    gp_assert_same_problem(ref, mutated, n = 3L),
    regexp = "same canonical",
    class = "gradepath_backend_error"
  )

  bad_answer <- same
  bad_answer$validation$grades <- rev(bad_answer$validation$grades)
  expect_error(
    gp_assert_same_answer(ref, bad_answer),
    regexp = "different grade",
    class = "gradepath_backend_error"
  )
})
