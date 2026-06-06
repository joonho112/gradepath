# Tests for the gp_report_card constructor and validator guards.
# Run with testthat::test_file() after sourcing the package R/ files.

# ---- helpers ---------------------------------------------------------------

rc_table <- function(J = 4L, selected_lambda = 0.25) {
  data.frame(
    id             = paste0("u", seq_len(J)),
    label          = paste0("Unit ", seq_len(J)),
    grade          = c(1L, 1L, 2L, 3L)[seq_len(J)],
    sort_rank      = seq_len(J),
    selected_lambda = rep(selected_lambda, J),
    posterior_mean = seq_len(J) + 0.0,
    lower          = seq_len(J) - 0.5,
    upper          = seq_len(J) + 0.5,
    estimate       = seq_len(J) + 0.1,
    se             = rep(0.3, J),
    stringsAsFactors = FALSE
  )
}

make_rc <- function(J = 4L, table = rc_table(J, selected_lambda), grades = c(1L, 1L, 2L, 3L)[seq_len(J)],
                    selected_lambda = 0.25) {
  new_gp_report_card(
    ids             = paste0("u", seq_len(J)),
    table           = table,
    selected_lambda = selected_lambda,
    grades          = grades,
    control         = gp_control()
  )
}

# ---- valid object passes ---------------------------------------------------

test_that("a well-formed gp_report_card validates", {
  rc <- make_rc()
  expect_silent(validate_gp_report_card(rc))
  expect_identical(validate_gp_report_card(rc), rc)
  expect_s3_class(rc, "gp_report_card")
  expect_identical(names(rc), .gp_report_card_fields)
  expect_identical(names(rc$table), .gp_report_card_table_columns)
})

test_that("constructor does not validate (cheap), per conventions §3", {
  # bad grades length, but new_ should still build the structure
  bad <- new_gp_report_card(
    ids = c("a", "b"), table = rc_table(2L),
    selected_lambda = 0.25, grades = c(1L, 2L, 3L),  # wrong length
    control = gp_control()
  )
  expect_s3_class(bad, "gp_report_card")
  expect_error(validate_gp_report_card(bad))          # validator catches it
})

# ---- field / order guards --------------------------------------------------

test_that("wrong field order is rejected", {
  rc <- make_rc()
  # rebuild as a plain list in the wrong order (avoid S3-subset coercion warning)
  reordered <- structure(
    unclass(rc)[c("table", "ids", "selected_lambda", "grades", "control",
                  "schema_version", "provenance", "warnings")],
    class = c("gp_report_card", "list")
  )
  expect_error(validate_gp_report_card(reordered), "canonical order")
})

# ---- ids guards ------------------------------------------------------------

test_that("non-unique ids are rejected", {
  rc <- make_rc()
  rc$ids <- c("u1", "u1", "u3", "u4")
  expect_error(validate_gp_report_card(rc), "unique")
})

test_that("table id column inconsistent with ids is rejected", {
  rc <- make_rc()
  rc$table$id <- c("x1", "x2", "x3", "x4")     # disagree with rc$ids
  expect_error(validate_gp_report_card(rc), "must equal")
})

# ---- table row-count guard -------------------------------------------------

test_that("table with wrong row count (!= J) is rejected", {
  rc <- make_rc()
  rc$table <- rc$table[1:3, ]                   # 3 rows, J = 4
  expect_error(validate_gp_report_card(rc), "J-row")
})

# ---- grades guards ---------------------------------------------------------

test_that("grades not integer length J is rejected", {
  rc <- make_rc()
  rc$grades <- c(1.5, 2.0, 3.0, 4.0)            # non-integer
  expect_error(validate_gp_report_card(rc), "integer")

  rc2 <- make_rc()
  rc2$grades <- c(1L, 2L, 3L)                   # wrong length
  expect_error(validate_gp_report_card(rc2), "length")
})

test_that("table grade column inconsistent with grades slot is rejected", {
  rc <- make_rc()
  rc$table$grade <- c(3L, 2L, 1L, 1L)           # disagree with rc$grades
  expect_error(validate_gp_report_card(rc), "must equal")
})

# ---- payload guards ---------------------------------------------------------

test_that("lower <= posterior_mean <= upper guard fires (posterior_mean col)", {
  rc <- make_rc()
  rc$table$lower[2] <- rc$table$posterior_mean[2] + 1   # lower > mean
  expect_error(validate_gp_report_card(rc), "lower")

  rc2 <- make_rc()
  rc2$table$upper[3] <- rc2$table$posterior_mean[3] - 1 # upper < mean
  expect_error(validate_gp_report_card(rc2), "upper|lower")
})

test_that("payload columns, sort_rank, selected_lambda, and se are required", {
  rc <- make_rc()
  rc$table$posterior_mean <- NULL
  expect_error(validate_gp_report_card(rc), "must include columns")

  rc2 <- make_rc()
  rc2$table$sort_rank <- c(1L, 2L, 4L, 3L)
  expect_error(validate_gp_report_card(rc2), "sort_rank")

  rc3 <- make_rc()
  rc3$table$selected_lambda <- rep(0.5, 4L)
  expect_error(validate_gp_report_card(rc3), "selected_lambda")

  rc4 <- make_rc()
  rc4$table$se[1] <- 0
  expect_error(validate_gp_report_card(rc4), "positive")
})

test_that("reduced tables are rejected by the payload schema", {
  tab <- rc_table()[, c("id", "label", "grade", "estimate", "se")]
  rc <- make_rc(table = tab)
  expect_error(validate_gp_report_card(rc), "must include columns")
})

# ---- selected_lambda guard -------------------------------------------------

test_that("selected_lambda out of [0,1] is rejected", {
  rc <- make_rc(selected_lambda = 1.5)
  expect_error(validate_gp_report_card(rc), "selected_lambda")
})
