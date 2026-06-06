frontier_pairwise <- function(pairwise_matrix = NULL,
                              ids = NULL,
                              control = gp_control(backend = "highs")) {
  if (is.null(pairwise_matrix)) {
    pairwise_matrix <- matrix(c(0.5, 0.9, 0.1, 0.5), 2, byrow = TRUE)
  }
  if (is.null(ids)) {
    ids <- paste0("unit_", seq_len(nrow(pairwise_matrix)))
  }
  pairwise_matrix <- as.matrix(pairwise_matrix)
  storage.mode(pairwise_matrix) <- "double"
  dimnames(pairwise_matrix) <- list(ids, ids)

  new_gp_pairwise(
    ids = ids,
    matrix = pairwise_matrix,
    power = 0L,
    cleanup = list(
      antisymmetry = TRUE,
      diagonal = .gp_pairwise_diagonal,
      zero_floor = .gp_pairwise_zero_floor
    ),
    source = list(
      stage = "posterior",
      rule = "outer_product",
      assumption = "one_level_independence"
    ),
    control = control
  )
}

frontier_fit <- function(ids, lambda, grades, control) {
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
      status = "optimal",
      n_units = length(ids)
    ),
    objective = list(
      value = 0,
      raw = 0,
      canonical = 0
    ),
    backend = list(
      name = control$backend,
      path = "fixture",
      encoding = "fixture",
      status = "optimal",
      objbound = NA_real_,
      mipgap = NA_real_,
      runtime = NA_real_,
      warmstart = "none",
      warm_start_from_lambda = NA_real_,
      warm_start_used = FALSE,
      problem_hash = "fixture"
    ),
    control = control
  )
}

frontier_path <- function(pairwise,
                          fit_specs = list(
                            list(lambda = 0.25, grades = c(1L, 1L)),
                            list(lambda = 1, grades = c(1L, 2L))
                          ),
                          enriched = TRUE,
                          carry_pairwise = FALSE) {
  control <- pairwise$control
  fits <- lapply(
    fit_specs,
    function(spec) frontier_fit(
      ids = pairwise$ids,
      lambda = spec$lambda,
      grades = spec$grades,
      control = control
    )
  )
  lambda_grid <- vapply(fit_specs, function(spec) spec$lambda, numeric(1))
  summary <- if (isTRUE(enriched)) {
    frontier_table(pairwise, fits)
  } else {
    data.frame(
      lambda = lambda_grid,
      grade_count = vapply(fits, function(fit) fit$summary$grade_count, integer(1))
    )
  }

  new_gp_grade_path(
    ids = pairwise$ids,
    lambda_grid = lambda_grid,
    fits = fits,
    summary = summary,
    backend = list(name = control$backend),
    selection = list(
      selected_lambda = 0.25,
      selection_rule = "baseline_lambda_0.25",
      endpoint_lambda = 1
    ),
    control = control,
    provenance = if (isTRUE(carry_pairwise)) list(pairwise = pairwise) else list()
  )
}

frontier_strict_matrix <- function(n = 4, high = 0.9) {
  out <- matrix(0.5, nrow = n, ncol = n)
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      out[i, j] <- high
      out[j, i] <- 1 - high
    }
  }
  out
}

test_that("frontier metrics port the v1 blockwise two-unit example", {
  pairwise <- frontier_pairwise()
  path <- frontier_path(pairwise)

  frontier <- gp_frontier(path, pairwise = pairwise)

  expect_s3_class(frontier, "gp_frontier")
  expect_equal(
    frontier$table$discordance_rate,
    c(0, 0.1),
    tolerance = 1e-8
  )
  expect_equal(frontier$table$reliability, c(1, 0.9), tolerance = 1e-8)
  expect_equal(frontier$table$tau_bar, c(0, 0.8), tolerance = 1e-8)
  expect_equal(
    frontier$table$reliability,
    1 - frontier$table$discordance_rate,
    tolerance = 1e-8
  )

  selected_matrix <- frontier$dr_matrix
  expect_identical(dim(selected_matrix), c(1L, 1L))
  expect_true(is.na(selected_matrix[1, 1]))

  endpoint <- gp_frontier(path, pairwise = pairwise, selected_lambda = 1)
  expect_identical(dim(endpoint$dr_matrix), c(2L, 2L))
  expect_equal(endpoint$dr_matrix[2, 1], 0.1, tolerance = 1e-8)
  expect_true(is.na(endpoint$dr_matrix[1, 2]))
  expect_true(is.na(endpoint$dr_matrix[1, 1]))
})

test_that("gp_frontier() can recompute from an explicit pairwise profile", {
  pairwise <- frontier_pairwise()
  path <- frontier_path(pairwise, enriched = FALSE, carry_pairwise = FALSE)

  expect_error(
    gp_frontier(path),
    regexp = "must include columns",
    class = "gradepath_error"
  )

  frontier <- gp_frontier(path, pairwise = pairwise)

  expect_equal(frontier$table$discordance_rate, c(0, 0.1), tolerance = 1e-8)
  expect_equal(frontier$table$tau_bar, c(0, 0.8), tolerance = 1e-8)
})

test_that("gp_frontier() scores the five naive benchmark rules when supplied", {
  pairwise <- frontier_pairwise(
    pairwise_matrix = frontier_strict_matrix(4),
    ids = paste0("unit_", 1:4)
  )
  path <- frontier_path(
    pairwise,
    fit_specs = list(
      list(lambda = 0.25, grades = c(1L, 2L, 3L, 4L)),
      list(lambda = 1, grades = c(1L, 2L, 3L, 4L))
    )
  )

  frontier <- gp_frontier(
    path,
    pairwise = pairwise,
    benchmark_scores = list(
      raw_estimate = stats::setNames(c(1, 2, 3, 4), rev(pairwise$ids)),
      posterior_mean = c(4, 3, 2, 1),
      linear_shrinkage = c(4, 3, 1, 2)
    )
  )

  expect_identical(
    frontier$benchmarks$benchmark,
    c(
      "raw_estimate_rank",
      "posterior_mean_rank",
      "linear_shrinkage_rank",
      "posterior_mean_decile",
      "posterior_mean_quartile"
    )
  )
  expect_equal(
    frontier$benchmarks$reliability,
    1 - frontier$benchmarks$discordance_rate,
    tolerance = 1e-8
  )
  expect_identical(frontier$benchmarks$grade_count, rep(4L, 5L))
  expect_equal(
    frontier$benchmarks$tau_bar[[1L]],
    frontier$table$tau_bar[[1L]],
    tolerance = 1e-8
  )
  expect_error(
    gp_frontier(
      path,
      pairwise = pairwise,
      benchmark_scores = list(
        raw_estimate = stats::setNames(c(4, 3, 2, 1), paste0("bad_", 1:4)),
        posterior_mean = c(4, 3, 2, 1),
        linear_shrinkage = c(4, 3, 1, 2)
      )
    ),
    regexp = "names must match",
    class = "gradepath_error"
  )
  expect_error(
    gp_frontier(
      path,
      benchmark_scores = list(
        raw_estimate = c(4, 3, 2, 1),
        posterior_mean = c(4, 3, 2, 1),
        linear_shrinkage = c(4, 3, 1, 2)
      )
    ),
    regexp = "require an explicit `pairwise`",
    class = "gradepath_error"
  )
})

test_that("gp_frontier() validates pairwise alignment", {
  pairwise <- frontier_pairwise()
  path <- frontier_path(pairwise)
  mismatched <- frontier_pairwise(ids = c("other_1", "other_2"))

  expect_error(
    gp_frontier(path, pairwise = mismatched),
    regexp = "must match",
    class = "gradepath_error"
  )
})

test_that("frontier_table() validates every fit id vector", {
  pairwise <- frontier_pairwise()
  control <- pairwise$control
  good <- frontier_fit(
    ids = pairwise$ids,
    lambda = 0.25,
    grades = c(1L, 1L),
    control = control
  )
  bad_second <- frontier_fit(
    ids = c(pairwise$ids[[1L]], "other_2"),
    lambda = 1,
    grades = c(1L, 2L),
    control = control
  )

  expect_error(
    frontier_table(pairwise, list(good, bad_second)),
    regexp = "mismatch at fit\\(s\\): 2",
    class = "gradepath_error"
  )
})

test_that("validate_gp_frontier() enforces reliability complement", {
  pairwise <- frontier_pairwise()
  frontier <- gp_frontier(frontier_path(pairwise), pairwise = pairwise)
  frontier$table$reliability[[2L]] <- 0.5

  expect_error(
    validate_gp_frontier(frontier),
    regexp = "1 - discordance_rate",
    class = "gradepath_error"
  )

  frontier <- gp_frontier(frontier_path(pairwise), pairwise = pairwise, selected_lambda = 1)
  frontier$dr_matrix[1, 2] <- 0.2
  expect_error(
    validate_gp_frontier(frontier),
    regexp = "lower triangle",
    class = "gradepath_error"
  )
})

test_that("validate_gp_frontier() enforces lower-triangular DR matrix storage", {
  pairwise <- frontier_pairwise()
  frontier <- gp_frontier(
    frontier_path(pairwise),
    pairwise = pairwise,
    selected_lambda = 1
  )

  bad_upper <- frontier
  bad_upper$dr_matrix[1, 2] <- 0.1
  expect_error(
    validate_gp_frontier(bad_upper),
    regexp = "lower triangle",
    class = "gradepath_error"
  )

  bad_diag <- frontier
  bad_diag$dr_matrix[1, 1] <- 0
  expect_error(
    validate_gp_frontier(bad_diag),
    regexp = "lower triangle",
    class = "gradepath_error"
  )
})

test_that("gp_grade_path() stores an enriched frontier summary", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  pairwise <- step45_problem(step45_strict3_matrix(), lambda = 1, backend = "highs")$pairwise

  path <- gp_grade_path(pairwise, lambda_grid = c(0, 0.25, 1))
  expected <- frontier_table(pairwise, path$fits)

  expect_true(all(.gp_frontier_table_columns %in% names(path$summary)))
  expect_equal(
    path$summary[, names(expected)],
    expected,
    tolerance = 1e-8
  )
  expect_null(path$provenance$pairwise)
})
