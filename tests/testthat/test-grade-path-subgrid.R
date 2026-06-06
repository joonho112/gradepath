# Regression — a gp_grade_path solved on a NON-anchor sub-grid is valid.
#
# Per the design, a gp_grade_path's `lambda_grid` slot must be strictly
# increasing / unique / in [0,1], but NEED NOT contain the parity anchors
# {0.25, 1.00} -- a user may solve a sub-grid such as c(0.1, 0.5, 0.9). The
# {0.25, 1.00} requirement lives only on gp_control's reference grid. This pins
# that validate_gp_grade_path uses the anchor-FREE core check, not the
# anchor-requiring control validator (QA-1.4 should-fix).
#
# Constructors are exercised with the canonical field shapes (summary /
# objective / backend are lists; the three audit slots are explicit), matching
# the working builders in test-class-grade.R.

test_that("a gp_grade_path on a non-anchor sub-grid validates", {
  ids <- c("a", "b", "c")
  sub_grid <- c(0.1, 0.5, 0.9)              # no 0.25, no 1.00

  grade_by_lambda <- list(
    c(1L, 2L, 3L),   # 0.1 -> 3 grades
    c(1L, 1L, 2L),   # 0.5 -> 2 grades
    c(1L, 1L, 1L)    # 0.9 -> 1 grade
  )

  fits <- lapply(seq_along(sub_grid), function(i) {
    grade <- grade_by_lambda[[i]]
    new_gp_grade_fit(
      ids        = ids,
      lambda     = sub_grid[[i]],
      assignment = data.frame(id = ids, grade = grade, stringsAsFactors = FALSE),
      summary    = list(grade_count = length(unique(grade))),
      objective  = list(value = -0.1 * i),
      backend    = list(name = "gurobi"),
      control    = gp_control(),
      schema_version = .gradepath_schema_version,
      provenance = .gradepath_new_provenance(producer = "test"),
      warnings   = character()
    )
  })

  grade_counts <- vapply(fits, function(f) f$summary$grade_count, integer(1))

  path <- new_gp_grade_path(
    ids         = ids,
    lambda_grid = sub_grid,
    fits        = fits,
    summary     = data.frame(lambda = sub_grid, grade_count = grade_counts),
    backend     = list(name = "gurobi"),
    selection   = list(selected_lambda = 0.5,            # in the sub-grid
                       selection_rule  = "manual",
                       endpoint_lambda = 0.9),           # in the sub-grid
    control     = gp_control(),
    schema_version = .gradepath_schema_version,
    provenance  = .gradepath_new_provenance(producer = "test"),
    warnings    = character()
  )

  expect_s3_class(path, "gp_grade_path")
  expect_no_error(validate_gp_grade_path(path))
})

test_that("gp_control still requires the parity anchors on its reference grid", {
  # The anchor requirement must remain on gp_control (control error), so the
  # relaxation above did not weaken the reference-grid contract.
  expect_error(
    gp_control(lambda_grid = c(0.1, 0.5, 0.9)),
    class = "gradepath_error"
  )
})
