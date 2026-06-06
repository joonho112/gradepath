make_step42_pairwise <- function(P, ids = rownames(P)) {
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
    control = gp_control(),
    provenance = .gradepath_new_provenance(producer = "test")
  )
}

test_that("decision matrix extraction uses row-major D and ignores grade auxiliaries", {
  ids <- c("a", "b", "c")
  P <- matrix(
    c(
      0.5, 0.9, 0.8,
      0.1, 0.5, 0.7,
      0.2, 0.3, 0.5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  problem <- gp_grade_problem(make_step42_pairwise(P), lambda = 1)
  D <- matrix(
    as.integer(c(
      0, 1, 1,
      0, 0, 1,
      0, 0, 0
    )),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  primal <- c(as.vector(t(D)), c(99, 88, 77))

  out <- .gradepath_grade_decision_matrix(problem, primal)

  expect_identical(out, D)
})

test_that("decision matrix extraction enforces length, binarity, and diagonal zero", {
  ids <- c("a", "b")
  P <- matrix(c(0.5, 0.9, 0.1, 0.5), 2, byrow = TRUE, dimnames = list(ids, ids))
  problem <- gp_grade_problem(make_step42_pairwise(P), lambda = 1)
  D <- matrix(c(0, 1, 0, 0), 2, byrow = TRUE)
  primal <- c(as.vector(t(D)), c(2, 1))

  expect_no_error(.gradepath_grade_decision_matrix(problem, primal + c(rep(5e-7, 4), 0, 0)))
  expect_error(.gradepath_grade_decision_matrix(problem, primal[-1]),
               regexp = "length", class = "gradepath_error")
  expect_error(.gradepath_grade_decision_matrix(problem, c(primal, 1)),
               regexp = "length", class = "gradepath_error")

  expect_no_error(.gradepath_grade_decision_matrix(
    problem,
    primal + c(0, 5e-7, -5e-7, 0, 0, 0)
  ))

  bad_nonbinary <- primal
  bad_nonbinary[2] <- 0.2
  expect_error(.gradepath_grade_decision_matrix(problem, bad_nonbinary),
               regexp = "non-binary", class = "gradepath_error")

  bad_above <- primal
  bad_above[2] <- 1 + 2e-6
  expect_error(.gradepath_grade_decision_matrix(problem, bad_above),
               regexp = "non-binary", class = "gradepath_error")

  bad_below <- primal
  bad_below[3] <- -2e-6
  expect_error(.gradepath_grade_decision_matrix(problem, bad_below),
               regexp = "non-binary", class = "gradepath_error")

  bad_diag <- primal
  bad_diag[1] <- 1
  expect_error(.gradepath_grade_decision_matrix(problem, bad_diag),
               regexp = "diagonal", class = "gradepath_error")
})

test_that("grade assignment is dense rank of row sums descending", {
  ids <- c("a", "b", "c", "d")
  D <- matrix(
    c(
      0, 1, 1, 1,
      0, 0, 1, 1,
      0, 0, 0, 0,
      0, 0, 1, 0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )

  assignment <- .gradepath_grade_assignment(D, ids)

  expect_identical(assignment$id, ids)
  expect_identical(assignment$grade, c(1L, 2L, 4L, 3L))
  expect_false(identical(
    assignment$grade,
    as.integer(match(colSums(D), sort(unique(colSums(D)), decreasing = TRUE)))
  ))

  diag_bad <- D
  diag_bad[1, 1] <- 1
  expect_error(.gradepath_grade_assignment(diag_bad, ids),
               regexp = "zero diagonal", class = "gradepath_error")
})

test_that("grade assignment rounds tolerated near-binary entries before row sums", {
  ids <- c("a", "b", "c")
  D <- matrix(
    c(
      0, 0.999999995, 0,
      1, 0, 0,
      0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )

  assignment <- .gradepath_grade_assignment(D, ids)

  expect_identical(assignment$grade, c(1L, 1L, 2L))
})

test_that("strict n=4 chain exports grades 1, 2, 3, 4", {
  ids <- as.character(1:4)
  D <- outer(seq_len(4), seq_len(4), FUN = "<") * 1
  dimnames(D) <- list(ids, ids)

  assignment <- .gradepath_grade_assignment(D, ids)

  expect_identical(rowSums(D), c("1" = 3, "2" = 2, "3" = 1, "4" = 0))
  expect_identical(assignment$grade, 1:4)
})

test_that("canonical objective equals twice the raw MILP objective across lambdas", {
  ids <- c("a", "b", "c")
  P <- matrix(
    c(
      0.5, 0.9, 0.8,
      0.1, 0.5, 0.6,
      0.2, 0.4, 0.5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  pair <- make_step42_pairwise(P)
  D <- matrix(
    c(
      0, 1, 1,
      0, 0, 1,
      0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  primal <- c(as.vector(t(D)), c(3, 2, 1))
  matrix_opt <- .gp_grade_matrix_for_optimization(pair)

  for (lambda in c(0, 0.25, 0.5, 1)) {
    problem <- gp_grade_problem(pair, lambda = lambda)
    canon <- .gradepath_grade_objective_value(matrix_opt, D, lambda = lambda)
    raw <- sum(problem$objective * primal)

    expect_equal(canon, 2 * raw, tolerance = 1e-12)
  }
})

test_that("two-unit threshold sign matches the grade-IP objective", {
  ids <- c("hi", "lo")
  P <- matrix(c(0.5, 0.81, 0.19, 0.5), 2, byrow = TRUE, dimnames = list(ids, ids))
  pair <- make_step42_pairwise(P)
  problem <- gp_grade_problem(pair, lambda = 0.25)
  separate <- c(0, 1, 0, 0, 2, 1)
  tie <- c(0, 0, 0, 0, 1, 1)
  reversed <- c(0, 0, 1, 0, 1, 2)
  D_sep <- .gradepath_grade_decision_matrix(problem, separate)
  assignment <- .gradepath_grade_assignment(D_sep, ids)
  canon_sep <- .gradepath_grade_objective_value(
    .gp_grade_matrix_for_optimization(pair),
    D_sep,
    lambda = 0.25
  )
  raw_sep <- sum(problem$objective * separate)

  expect_identical(assignment$grade, c(1L, 2L))
  expect_equal(raw_sep, -0.0125, tolerance = 1e-12)
  expect_equal(canon_sep, -0.025, tolerance = 1e-12)
  expect_equal(canon_sep, 2 * raw_sep, tolerance = 1e-12)
  expect_equal(sum(problem$objective * tie), 0, tolerance = 1e-12)
  expect_lt(sum(problem$objective * separate), sum(problem$objective * tie))
  expect_gt(sum(problem$objective * reversed), sum(problem$objective * tie))

  P2 <- matrix(c(0.5, 0.79, 0.21, 0.5), 2, byrow = TRUE, dimnames = list(ids, ids))
  pair2 <- make_step42_pairwise(P2)
  problem2 <- gp_grade_problem(pair2, lambda = 0.25)

  expect_gt(sum(problem2$objective * separate), sum(problem2$objective * tie))
})

test_that("two-unit at-threshold separation ties the no-separation objective", {
  ids <- c("hi", "lo")
  P <- matrix(c(0.5, 0.80, 0.20, 0.5), 2, byrow = TRUE, dimnames = list(ids, ids))
  pair <- make_step42_pairwise(P)
  problem <- gp_grade_problem(pair, lambda = 0.25)
  separate <- c(0, 1, 0, 0, 2, 1)
  tie <- c(0, 0, 0, 0, 1, 1)

  expect_equal(sum(problem$objective * separate), 0, tolerance = 1e-12)
  expect_equal(sum(problem$objective * tie), 0, tolerance = 1e-12)
})

test_that("objective helper validates dimensions and lambda", {
  D <- diag(2)
  P <- matrix(c(0, 0.8, 0.2, 0), 2, byrow = TRUE)

  expect_error(.gradepath_grade_objective_value(P, D[1, , drop = FALSE], 0.25),
               regexp = "same dimensions", class = "gradepath_error")
  expect_error(.gradepath_grade_objective_value(P, D, 1.5),
               regexp = "lambda", class = "gradepath_error")
  expect_error(.gradepath_grade_objective_value(P + diag(0.5, 2), matrix(0, 2, 2), 0.25),
               regexp = "zero diagonals", class = "gradepath_error")
})
