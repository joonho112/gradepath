# Tests for gp_grade_fit and gp_grade_path.
#
# Covers: valid gp_grade_fit / gp_grade_path validate; the cross-checks fire --
#   misaligned fit lambda vs lambda_grid; selected_lambda not in lambda_grid;
#   fit ids != path ids; summary grade_count vs stored fit mismatch; plus the
#   per-object guards (grade-count vs assignment; non-contiguous / non-integer
#   grades; ids length mismatch; lambda out of [0,1]; wrong slot order; non-
#   gp_control control; and confirmation that the v1 replication_mode lock is
#   GONE -- a selected_lambda != 0.25 path with a krw-ish setup still validates).

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

gpB_ids <- function() c("a", "b", "c", "d")

# A valid gp_grade_fit at a given lambda with a chosen grade vector. `grade`
# defaults to a 2-grade split of 4 ids.
gpB_valid_grade_fit <- function(lambda = 0.25,
                                ids = gpB_ids(),
                                grade = c(1L, 1L, 2L, 2L),
                                control = gp_control()) {
  k <- length(unique(grade))
  new_gp_grade_fit(
    ids = ids,
    lambda = lambda,
    assignment = data.frame(id = ids, grade = grade, stringsAsFactors = FALSE),
    summary = list(grade_count = k),
    objective = list(value = -0.123),
    backend = list(name = "highs"),
    control = control,
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(producer = "test"),
    warnings = character()
  )
}

# A valid gp_grade_path over a small grid that contains the parity anchors
# {0.25, 1.00} (required by the shared lambda-grid validator). Each fit's lambda
# aligns to the grid; grade_count varies along the grid.
gpB_grid <- function() c(0, 0.25, 0.5, 1)

gpB_valid_grade_path <- function(ids = gpB_ids(),
                                 grid = gpB_grid(),
                                 selected_lambda = 0.25,
                                 control = gp_control()) {
  # Distinct grade vectors so grade_count differs across the grid (coarser as
  # lambda grows is the usual shape; here just make them distinct + valid).
  grade_by_lambda <- list(
    c(1L, 2L, 3L, 4L),  # lambda 0.00 -> 4 grades
    c(1L, 1L, 2L, 2L),  # lambda 0.25 -> 2 grades
    c(1L, 1L, 2L, 2L),  # lambda 0.50 -> 2 grades
    c(1L, 1L, 1L, 1L)   # lambda 1.00 -> 1 grade
  )
  fits <- lapply(seq_along(grid), function(i) {
    gpB_valid_grade_fit(
      lambda = grid[[i]],
      ids = ids,
      grade = grade_by_lambda[[i]],
      control = control
    )
  })
  grade_counts <- vapply(fits, function(f) f$summary$grade_count, integer(1))
  summary <- data.frame(lambda = grid, grade_count = grade_counts)

  new_gp_grade_path(
    ids = ids,
    lambda_grid = grid,
    fits = fits,
    summary = summary,
    backend = list(name = "highs"),
    selection = list(
      selected_lambda = selected_lambda,
      selection_rule = "first_at_target_count",
      endpoint_lambda = 1
    ),
    control = control,
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(producer = "test"),
    warnings = character()
  )
}

# ===========================================================================
# gp_grade_fit
# ===========================================================================

test_that("a valid gp_grade_fit constructs and validates", {
  fit <- gpB_valid_grade_fit()
  expect_s3_class(fit, "gp_grade_fit")
  expect_identical(names(fit), .gp_grade_fit_fields)
  expect_identical(fit$summary$grade_count, 2L)
  expect_identical(fit$assignment$id, gpB_ids())
  expect_type(fit$assignment$grade, "integer")
  expect_no_error(validate_gp_grade_fit(fit))
})

test_that("gp_grade_fit validates optional warm-start metadata when present", {
  fit <- gpB_valid_grade_fit()
  expect_no_error(validate_gp_grade_fit(fit))

  fit$backend$warm_start_from_lambda <- 1
  fit$backend$warm_start_from_status <- "gap_reached"
  fit$backend$warm_start_from_acceptance_ready <- FALSE
  fit$backend$warm_start_used <- TRUE
  expect_no_error(validate_gp_grade_fit(fit))

  bad_lambda <- fit
  bad_lambda$backend$warm_start_from_lambda <- "1"
  expect_error(
    validate_gp_grade_fit(bad_lambda),
    regexp = "warm_start_from_lambda",
    class = "gradepath_error"
  )

  bad_status <- fit
  bad_status$backend$warm_start_from_status <- ""
  expect_error(
    validate_gp_grade_fit(bad_status),
    regexp = "warm_start_from_status",
    class = "gradepath_error"
  )

  bad_ready <- fit
  bad_ready$backend$warm_start_from_acceptance_ready <- "FALSE"
  expect_error(
    validate_gp_grade_fit(bad_ready),
    regexp = "warm_start_from_acceptance_ready",
    class = "gradepath_error"
  )
})

test_that("gp_grade_fit normalizes / accepts contiguous integer grades", {
  # Already-contiguous {1,2,3} across 4 ids with a repeat.
  fit <- gpB_valid_grade_fit(grade = c(1L, 2L, 2L, 3L))
  expect_identical(fit$summary$grade_count, 3L)
})

test_that("gp_grade_fit rejects non-contiguous grade labels", {
  # {1, 3} is not contiguous from 1.
  fit_raw <- gpB_ids()
  expect_error(
    new_gp_grade_fit(
      ids = fit_raw,
      lambda = 0.25,
      assignment = data.frame(id = fit_raw, grade = c(1L, 1L, 3L, 3L)),
      summary = list(grade_count = 2L),
      objective = list(value = 0),
      backend = list(name = "highs"),
      control = gp_control()
    ),
    regexp = "contiguous",
    class = "gradepath_error"
  )
})

test_that("gp_grade_fit rejects non-integer grades", {
  ids <- gpB_ids()
  expect_error(
    new_gp_grade_fit(
      ids = ids,
      lambda = 0.25,
      assignment = data.frame(id = ids, grade = c(1, 1.5, 2, 2)),
      summary = list(grade_count = 3L),
      objective = list(value = 0),
      backend = list(name = "highs"),
      control = gp_control()
    ),
    regexp = "integer",
    class = "gradepath_error"
  )
})

test_that("gp_grade_fit rejects summary$grade_count != realized count", {
  ids <- gpB_ids()
  expect_error(
    new_gp_grade_fit(
      ids = ids,
      lambda = 0.25,
      assignment = data.frame(id = ids, grade = c(1L, 1L, 2L, 2L)),  # 2 grades
      summary = list(grade_count = 3L),                              # claims 3
      objective = list(value = 0),
      backend = list(name = "highs"),
      control = gp_control()
    ),
    regexp = "grade_count",
    class = "gradepath_error"
  )
})

test_that("gp_grade_fit rejects assignment$id != ids", {
  ids <- gpB_ids()
  expect_error(
    new_gp_grade_fit(
      ids = ids,
      lambda = 0.25,
      assignment = data.frame(id = c("a", "b", "c", "Z"), grade = c(1L, 1L, 2L, 2L)),
      summary = list(grade_count = 2L),
      objective = list(value = 0),
      backend = list(name = "highs"),
      control = gp_control()
    ),
    regexp = "must exactly match",
    class = "gradepath_error"
  )
})

test_that("gp_grade_fit rejects assignment row count != length(ids)", {
  ids <- gpB_ids()
  expect_error(
    new_gp_grade_fit(
      ids = ids,
      lambda = 0.25,
      assignment = data.frame(id = ids[1:3], grade = c(1L, 1L, 2L)),
      summary = list(grade_count = 2L),
      objective = list(value = 0),
      backend = list(name = "highs"),
      control = gp_control()
    ),
    regexp = "one row per id",
    class = "gradepath_error"
  )
})

test_that("gp_grade_fit rejects lambda outside [0, 1]", {
  expect_error(gpB_valid_grade_fit(lambda = 1.5), regexp = "lambda",
               class = "gradepath_error")
  expect_error(gpB_valid_grade_fit(lambda = -0.1), regexp = "lambda",
               class = "gradepath_error")
})

test_that("gp_grade_fit rejects wrong slot order", {
  fit <- gpB_valid_grade_fit()
  reordered <- fit[c("lambda", "ids", "assignment", "summary", "objective",
                     "backend", "control", "schema_version", "provenance",
                     "warnings")]
  class(reordered) <- c("gp_grade_fit", "list")
  expect_error(validate_gp_grade_fit(reordered), regexp = "canonical order",
               class = "gradepath_error")
})

test_that("gp_grade_fit rejects a non-gp_control control slot", {
  fit <- gpB_valid_grade_fit()
  fit$control <- list(backend = "highs")
  expect_error(validate_gp_grade_fit(fit), class = "gradepath_error")
})

# ===========================================================================
# gp_grade_path
# ===========================================================================

test_that("a valid gp_grade_path constructs and validates", {
  path <- gpB_valid_grade_path()
  expect_s3_class(path, "gp_grade_path")
  expect_identical(names(path), .gp_grade_path_fields)
  expect_length(path$fits, length(gpB_grid()))
  expect_true(all(vapply(path$fits, inherits, logical(1), "gp_grade_fit")))
  expect_identical(path$summary$grade_count, c(4L, 2L, 2L, 1L))
  expect_no_error(validate_gp_grade_path(path))
})

# --- cross-check: misaligned fit lambda vs lambda_grid ---------------------

test_that("gp_grade_path rejects a fit whose lambda != lambda_grid[i]", {
  path <- gpB_valid_grade_path()
  # Break the 2nd fit's lambda (0.25 -> 0.30) without touching the grid.
  path$fits[[2]]$lambda <- 0.30
  expect_error(
    validate_gp_grade_path(path),
    regexp = "align to `gp_grade_path\\$lambda_grid`|align",
    class = "gradepath_error"
  )
})

# --- cross-check: selected_lambda not in lambda_grid -----------------------

test_that("gp_grade_path rejects selected_lambda not in lambda_grid", {
  path <- gpB_valid_grade_path()
  path$selection$selected_lambda <- 0.37   # not on {0, 0.25, 0.5, 1}
  expect_error(
    validate_gp_grade_path(path),
    regexp = "exactly one member",
    class = "gradepath_error"
  )
})

test_that("gp_grade_path rejects endpoint_lambda not in lambda_grid", {
  path <- gpB_valid_grade_path()
  path$selection$endpoint_lambda <- 0.9    # not on the grid
  expect_error(
    validate_gp_grade_path(path),
    regexp = "exactly one member",
    class = "gradepath_error"
  )
})

# --- cross-check: fit ids != path ids --------------------------------------

test_that("gp_grade_path rejects a fit whose ids != path ids", {
  path <- gpB_valid_grade_path()
  # Rebuild fit 3 on a DIFFERENT id set (same length so the grade vector fits).
  bad_ids <- c("a", "b", "c", "ZZ")
  path$fits[[3]] <- gpB_valid_grade_fit(
    lambda = gpB_grid()[[3]],
    ids = bad_ids,
    grade = c(1L, 1L, 2L, 2L)
  )
  expect_error(
    validate_gp_grade_path(path),
    regexp = "same ids",
    class = "gradepath_error"
  )
})

# --- cross-check: summary grade_count vs stored fit mismatch ---------------

test_that("gp_grade_path rejects summary grade_count != stored fits' counts", {
  path <- gpB_valid_grade_path()
  # Corrupt one summary count so it disagrees with the fit it summarizes.
  path$summary$grade_count[2] <- 3L   # fit[[2]] actually has 2 grades
  expect_error(
    validate_gp_grade_path(path),
    regexp = "align to the stored fits",
    class = "gradepath_error"
  )
})

test_that("gp_grade_path rejects summary$lambda misaligned to the grid", {
  path <- gpB_valid_grade_path()
  path$summary$lambda[3] <- 0.55   # grid has 0.5 there
  expect_error(
    validate_gp_grade_path(path),
    regexp = "align exactly",
    class = "gradepath_error"
  )
})

# --- structural guards -----------------------------------------------------

test_that("gp_grade_path rejects fits length != grid length", {
  path <- gpB_valid_grade_path()
  path$fits <- path$fits[1:3]   # one short
  expect_error(
    validate_gp_grade_path(path),
    regexp = "one fit per lambda",
    class = "gradepath_error"
  )
})

test_that("gp_grade_path rejects a non-strictly-increasing lambda_grid", {
  # Build a path then corrupt the grid ordering (and matching summary/fits) so
  # the lambda-grid validator is what trips.
  path <- gpB_valid_grade_path()
  path$lambda_grid <- c(0, 0.5, 0.25, 1)   # not increasing
  expect_error(
    validate_gp_grade_path(path),
    regexp = "increasing",
    class = "gradepath_error"
  )
})

test_that("gp_grade_path rejects wrong slot order", {
  path <- gpB_valid_grade_path()
  reordered <- path[c("lambda_grid", "ids", "fits", "summary", "backend",
                      "selection", "control", "schema_version", "provenance",
                      "warnings")]
  class(reordered) <- c("gp_grade_path", "list")
  expect_error(validate_gp_grade_path(reordered), regexp = "canonical order",
               class = "gradepath_error")
})

# --- v2: the replication_mode lock is GONE ---------------------------------

test_that("gp_grade_path accepts a non-0.25 selected_lambda (no replication lock)", {
  # v1 forced selected_lambda = 0.25 under replication_mode. v2 has no such
  # field/branch: selecting lambda = 0.5 (a grid member) must validate cleanly.
  path <- gpB_valid_grade_path(selected_lambda = 0.5)
  expect_s3_class(path, "gp_grade_path")
  expect_identical(path$selection$selected_lambda, 0.5)
})

test_that("gp_grade_path accepts an arbitrary first-class selection_rule string", {
  path <- gpB_valid_grade_path()
  path$selection$selection_rule <- "elbow_of_DR_curve"
  expect_no_error(validate_gp_grade_path(path))
})
