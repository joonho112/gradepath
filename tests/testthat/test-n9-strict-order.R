test_that("N9 strict-order fixture has grade-label range [1, n]", {
  skip_if_not_installed("digest")
  fixture <- step45_problem(step45_n9_matrix(), lambda = 1)
  problem <- fixture$problem

  expect_equal(problem$lower[problem$grade_indices], rep(1, 4))
  expect_equal(problem$upper[problem$grade_indices], rep(4, 4))
  expect_equal(problem$m_big, 4)
  expect_identical(problem$metadata$label_bounds, "[1,n]")

  bad_upper <- problem
  bad_upper$upper[bad_upper$grade_indices] <- 3
  expect_error(
    validate_gp_grade_problem(bad_upper),
    regexp = "\\[1, n\\]",
    class = "gradepath_grade_error"
  )

  bad_m <- problem
  bad_m$m_big <- 3
  expect_error(
    validate_gp_grade_problem(bad_m),
    regexp = "must equal n",
    class = "gradepath_grade_error"
  )
})

test_that("N9 strict-order fixture solves to four distinct grades", {
  skip_if_not_installed("digest")
  fixture <- step45_problem(step45_n9_matrix(), lambda = 1)
  problem <- fixture$problem
  results <- step45_solver_results(problem, fixture$matrix_opt)

  step45_skip_without_solver(results)

  for (name in names(results)) {
    res <- results[[name]]
    expect_identical(res$status, "optimal")
    expect_identical(res$validation$grades, 1:4)
    expect_identical(length(unique(res$validation$grades)), 4L)
    expect_equal(unname(rowSums(res$validation$decision_matrix)), c(3, 2, 1, 0),
                 tolerance = 1e-8)
    expect_equal(res$objval, -5.56, tolerance = 1e-8)
    expect_equal(res$validation$canon, -11.12, tolerance = 1e-8)
    expect_equal(res$validation$canon, 2 * res$objval, tolerance = 1e-8)
    expect_identical(res$validation$problem_hash, gp_problem_hash(problem))
    expect_equal(
      sort(as.integer(round(res$primal[problem$grade_indices]))),
      1:4,
      tolerance = 1e-8
    )
  }
})

test_that("N9 registry target is present as an exact M1 count gate", {
  skip_if_not_installed("digest")
  row <- gp_registry[gp_registry$id == "n9_strict4_ngrades", , drop = FALSE]

  expect_equal(nrow(row), 1L)
  expect_equal(as.numeric(row$paper_value), 4)
  expect_identical(row$unit, "count")
  expect_equal(as.numeric(row$tolerance), 0)
  expect_identical(row$class, "exact")
  expect_identical(row$milestone, "M1")
})
