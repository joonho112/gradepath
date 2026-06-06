test_that("gp_grade_problem builds the row-major grade-IP core", {
  ids <- c("a", "b", "c")
  P <- matrix(
    c(
      0.5, 0.8, 0.6,
      0.2, 0.5, 0.7,
      0.4, 0.3, 0.5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  pair <- new_gp_pairwise(
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

  problem <- gp_grade_problem(pair, lambda = 0.25)

  expect_s3_class(problem, "gp_grade_problem")
  expect_identical(names(problem), .gp_grade_problem_fields)
  expect_identical(problem$n_units, 3L)
  expect_identical(problem$d_indices, matrix(seq_len(9), 3, 3, byrow = TRUE))
  expect_identical(problem$grade_indices, 10:12)
  expect_identical(problem$constraint_ncol, 12L)
  expect_identical(problem$constraint_nrow, 12L)
  expect_identical(problem$m_big, 3)

  expect_identical(problem$diagonal_positions, gp_diag_positions(3))
  expect_equal(problem$lower[problem$diagonal_positions], c(0, 0, 0))
  expect_equal(problem$upper[problem$diagonal_positions], c(0, 0, 0))
  expect_equal(problem$objective[problem$diagonal_positions], c(0, 0, 0))

  grade_slots <- problem$grade_indices
  expect_equal(problem$lower[grade_slots], c(1, 1, 1))
  expect_equal(problem$upper[grade_slots], c(3, 3, 3))
  expect_identical(problem$types[1:9], rep("B", 9))
  expect_identical(problem$types[10:12], rep("I", 3))
})

test_that("gp_grade_problem objective follows row-major P_ji - lambda P_ij", {
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
  pair <- new_gp_pairwise(
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

  problem <- gp_grade_problem(pair, lambda = 0.25)
  P0 <- P
  diag(P0) <- 0
  expected <- c(as.vector(t(t(P0) - 0.25 * P0)), rep(0, 3))
  expected[gp_diag_positions(3)] <- 0

  expect_equal(problem$objective, expected)
  expect_equal(problem$objective[2], P0[2, 1] - 0.25 * P0[1, 2])
  expect_equal(problem$objective[4], P0[1, 2] - 0.25 * P0[2, 1])
})

test_that("gp_grade_problem encodes the big-M off-diagonal constraints", {
  ids <- c("a", "b")
  P <- matrix(
    c(0.5, 0.9, 0.1, 0.5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  pair <- new_gp_pairwise(
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

  problem <- gp_grade_problem(pair, lambda = 1)
  expect_identical(problem$constraint_nrow, 4L)
  expect_identical(problem$constraint_i, c(1L, 2L, 1L, 2L, 1L, 2L,
                                           3L, 4L, 3L, 4L, 3L, 4L))
  expect_identical(problem$constraint_j, c(2L, 3L, 5L, 6L, 6L, 5L,
                                           2L, 3L, 5L, 6L, 6L, 5L))
  expect_equal(problem$constraint_x, c(-2, -2, 1, 1, -1, -1, 2, 2, -1, -1, 1, 1))
  expect_equal(problem$rhs, c(0, 0, 1, 1))
})

test_that("gp_problem_hash is stable and encoding-invariant", {
  skip_if_not_installed("digest")

  ids <- c("a", "b", "c")
  P <- matrix(
    c(
      0.5, 0.8, 0.6,
      0.2, 0.5, 0.7,
      0.4, 0.3, 0.5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  pair <- new_gp_pairwise(
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

  big_m <- gp_grade_problem(pair, lambda = 0.25, encoding = "bigM")
  indicator <- gp_grade_problem(pair, lambda = 0.25, encoding = "indicator")
  other_lambda <- gp_grade_problem(pair, lambda = 1, encoding = "bigM")

  expect_identical(gp_problem_hash(big_m), gp_problem_hash(big_m))
  expect_identical(gp_problem_hash(big_m), gp_problem_hash(indicator))
  expect_identical(
    gp_problem_feasible_signature(big_m),
    gp_problem_feasible_signature(indicator)
  )
  expect_identical(
    gp_problem_feasible_signature(big_m),
    gp_problem_feasible_signature(other_lambda)
  )
  expect_false(identical(gp_problem_hash(big_m), gp_problem_hash(other_lambda)))
})

test_that("gp_problem_feasible_signature is constraint-complete and triplet-order-invariant", {
  skip_if_not_installed("digest")

  ids <- c("a", "b", "c")
  P <- matrix(
    c(
      0.5, 0.8, 0.6,
      0.2, 0.5, 0.7,
      0.4, 0.3, 0.5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  pair <- new_gp_pairwise(
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
  problem <- gp_grade_problem(pair, lambda = 0.25, encoding = "bigM")
  signature <- gp_problem_feasible_signature(problem)

  shuffled <- problem
  ord <- rev(seq_along(shuffled$constraint_i))
  shuffled$constraint_i <- shuffled$constraint_i[ord]
  shuffled$constraint_j <- shuffled$constraint_j[ord]
  shuffled$constraint_x <- shuffled$constraint_x[ord]
  shuffled <- validate_gp_grade_problem(shuffled)
  expect_identical(gp_problem_feasible_signature(shuffled), signature)

  changed_coef <- problem
  changed_coef$constraint_x[[1L]] <- changed_coef$constraint_x[[1L]] + 0.25
  changed_coef <- validate_gp_grade_problem(changed_coef)
  expect_false(identical(
    gp_problem_feasible_signature(changed_coef),
    signature
  ))
  expect_false(identical(gp_problem_hash(changed_coef), gp_problem_hash(problem)))

  changed_rhs <- problem
  changed_rhs$rhs[[1L]] <- changed_rhs$rhs[[1L]] + 1
  changed_rhs <- validate_gp_grade_problem(changed_rhs)
  expect_false(identical(
    gp_problem_feasible_signature(changed_rhs),
    signature
  ))
  expect_false(identical(gp_problem_hash(changed_rhs), gp_problem_hash(problem)))

  core <- .gp_problem_feasible_core(problem)
  changed_bound <- problem
  changed_bound$upper[changed_bound$grade_indices[[1L]]] <- changed_bound$n_units - 1
  expect_false(identical(.gp_problem_feasible_core(changed_bound)$upper, core$upper))
  expect_error(
    gp_problem_feasible_signature(changed_bound),
    regexp = "Grade-label",
    class = "gradepath_error"
  )

  changed_type <- problem
  changed_type$types[changed_type$grade_indices[[1L]]] <- "C"
  expect_false(identical(.gp_problem_feasible_core(changed_type)$types, core$types))
  expect_error(
    gp_problem_feasible_signature(changed_type),
    regexp = "binary D then integer",
    class = "gradepath_error"
  )

  # CCR-14: the feasible signature digests constraint DIRECTION (sense), but no
  # prior test proves it because validate_gp_grade_problem() rejects any
  # non-"<=" direction before the signature is ever computed. Mutate one row's
  # direction on a COPY and call the digest core DIRECTLY (bypassing the
  # validator): both the per-row constraint_rows digest AND the full
  # .gp_problem_feasible_core() result must change, confirming sense is bound.
  changed_direction <- problem
  changed_direction$direction[[1L]] <- ">="
  expect_false(identical(
    .gp_problem_feasible_core(changed_direction)$constraint_rows,
    core$constraint_rows
  ))
  expect_false(identical(
    .gp_problem_feasible_core(changed_direction),
    core
  ))
  # And the validator still guards the contract that the signature relies on.
  expect_error(
    validate_gp_grade_problem(changed_direction),
    regexp = "<=",
    class = "gradepath_error"
  )
})

test_that("validate_gp_grade_problem rejects diagonal and grade-bound drift", {
  ids <- c("a", "b")
  P <- matrix(
    c(0.5, 0.9, 0.1, 0.5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  pair <- new_gp_pairwise(
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
  problem <- gp_grade_problem(pair, lambda = 0.25)

  bad_diag <- problem
  bad_diag$upper[gp_diag_positions(2)[1]] <- 1
  expect_error(validate_gp_grade_problem(bad_diag), regexp = "diagonal",
               class = "gradepath_error")

  bad_bound <- problem
  bad_bound$upper[bad_bound$grade_indices[1]] <- 1
  expect_error(validate_gp_grade_problem(bad_bound), regexp = "Grade-label",
               class = "gradepath_error")

  bad_d <- problem
  bad_d$upper[setdiff(seq_len(4), gp_diag_positions(2))[1]] <- 0
  expect_error(validate_gp_grade_problem(bad_d), regexp = "Off-diagonal",
               class = "gradepath_error")
})

test_that("gp_grade_problem rejects singleton pairwise inputs", {
  pair <- new_gp_pairwise(
    ids = "a",
    matrix = matrix(0.5, nrow = 1, dimnames = list("a", "a")),
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

  expect_error(gp_grade_problem(pair, lambda = 0.25), regexp = "at least two",
               class = "gradepath_error")
})

test_that("the strict n=4 chain is feasible under M=n and grade bounds [1,n]", {
  ids <- as.character(1:4)
  P <- matrix(
    c(
      0.5, 0.95, 0.97, 0.99,
      0.05, 0.5, 0.95, 0.97,
      0.03, 0.05, 0.5, 0.95,
      0.01, 0.03, 0.05, 0.5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  pair <- new_gp_pairwise(
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
  problem <- gp_grade_problem(pair, lambda = 1)

  D <- outer(seq_len(4), seq_len(4), FUN = "<") * 1
  primal <- c(as.vector(t(D)), c(4, 3, 2, 1))
  lhs <- numeric(problem$constraint_nrow)
  for (k in seq_along(problem$constraint_x)) {
    lhs[problem$constraint_i[k]] <- lhs[problem$constraint_i[k]] +
      problem$constraint_x[k] * primal[problem$constraint_j[k]]
  }

  expect_true(all(primal[problem$diagonal_positions] == 0))
  expect_true(all(primal >= problem$lower))
  expect_true(all(primal <= problem$upper))
  expect_true(all(lhs <= problem$rhs + 1e-12))
})
