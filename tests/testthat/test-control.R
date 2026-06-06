# Tests for the pruned gp_control() surface + the validation/math utils.
# testthat 3e. These pin: the gp_control class + its six fields + defaults; the
# validation rules (lambda / backend / interval / seed); the ergonomic
# time_limit / mip_gap normalization and double-specify rejection; and the three
# structural/numeric helpers monotone_relabel / gp_diag_positions /
# .gp_row_logsumexp.

# ---------------------------------------------------------------------------
# gp_control(): class + the six fields + defaults
# ---------------------------------------------------------------------------

test_that("gp_control() returns a gp_control with exactly six user fields + audit slots", {
  ctl <- gp_control()

  expect_s3_class(ctl, "gp_control")

  # Exact canonical field order: six user fields then the two audit slots.
  expect_identical(
    names(ctl),
    c("lambda_grid", "backend", "precision_rule", "interval_level",
      "solver_options", "seed", "schema_version", "provenance")
  )

  # The six user-facing fields.
  expect_true(all(
    c("lambda_grid", "backend", "precision_rule", "interval_level",
      "solver_options", "seed") %in% names(ctl)
  ))

  # Dropped v1 fields must NOT be present.
  expect_false(any(
    c("workflow", "replication_mode", "target_bundle", "selection_rule",
      "backend_fallback", "grade_backend", "claim_bearing", "claim_scope",
      "backend_claim_tier") %in% names(ctl)
  ))
})

test_that("gp_control() defaults match the spec", {
  ctl <- gp_control()

  # interval_level default 0.90.
  expect_equal(ctl$interval_level, 0.90)

  # backend default "gurobi" (via getOption("gradepath.backend", "gurobi")).
  expect_identical(ctl$backend, "gurobi")

  # precision_rule default "none".
  expect_identical(ctl$precision_rule, "none")

  # lambda_grid default = seq(0,1,0.01) and contains the parity anchors.
  expect_equal(ctl$lambda_grid, seq(0, 1, by = 0.01))
  expect_true(any(abs(ctl$lambda_grid - 0.25) < 1e-8))
  expect_true(any(abs(ctl$lambda_grid - 1.00) < 1e-8))

  # solver_options default empty named list; seed default NULL.
  expect_identical(ctl$solver_options, list())
  expect_null(ctl$seed)

  # schema tag is v2.
  expect_identical(ctl$schema_version, "v2")

  # provenance is a non-empty named list (audit stamp).
  expect_type(ctl$provenance, "list")
  expect_true(length(ctl$provenance) > 0L)
  expect_true(!is.null(names(ctl$provenance)) && all(nzchar(names(ctl$provenance))))
})

test_that("gp_control() honors the gradepath.backend option for the default", {
  old <- options(gradepath.backend = "highs")
  on.exit(options(old), add = TRUE)
  expect_identical(gp_control()$backend, "highs")
})

# ---------------------------------------------------------------------------
# gp_control(): lambda_grid validation
# ---------------------------------------------------------------------------

test_that("lambda_grid validation rejects bad grids", {
  # Out of [0, 1].
  expect_error(gp_control(lambda_grid = c(-0.1, 0.25, 1)), class = "gp_control_error")
  expect_error(gp_control(lambda_grid = c(0.25, 1, 1.2)), class = "gp_control_error")
  # Not strictly increasing.
  expect_error(gp_control(lambda_grid = c(1, 0.25, 0)), class = "gp_control_error")
  # Duplicated.
  expect_error(gp_control(lambda_grid = c(0.25, 0.25, 1)), class = "gp_control_error")
  # Missing a parity anchor (no 0.25).
  expect_error(gp_control(lambda_grid = c(0, 0.5, 1)), class = "gp_control_error")
  # Missing the other anchor (no 1.00).
  expect_error(gp_control(lambda_grid = c(0, 0.25, 0.5)), class = "gp_control_error")
  # Non-finite.
  expect_error(gp_control(lambda_grid = c(0.25, NA, 1)), class = "gradepath_error")
})

test_that("lambda_grid accepts a valid narrow grid containing the anchors", {
  ctl <- gp_control(lambda_grid = c(0, 0.25, 0.5, 1))
  expect_equal(ctl$lambda_grid, c(0, 0.25, 0.5, 1))
  # NULL falls back to the reference grid.
  expect_equal(gp_control(lambda_grid = NULL)$lambda_grid, seq(0, 1, by = 0.01))
})

# ---------------------------------------------------------------------------
# gp_control(): backend validation + alias normalization
# ---------------------------------------------------------------------------

test_that("backend validation accepts the four canonical backends", {
  expect_identical(gp_control(backend = "gurobi")$backend, "gurobi")
  expect_identical(gp_control(backend = "highs")$backend, "highs")
  expect_identical(gp_control(backend = "glpk")$backend, "glpk")
  expect_identical(gp_control(backend = "symphony")$backend, "symphony")
})

test_that("backend roi_* aliases normalize to canonical short names", {
  expect_identical(gp_control(backend = "roi_highs")$backend, "highs")
  expect_identical(gp_control(backend = "roi_glpk")$backend, "glpk")
  expect_identical(gp_control(backend = "roi_symphony")$backend, "symphony")
})

test_that("backend validation rejects unknown backends and 'auto'", {
  expect_error(gp_control(backend = "auto"), class = "gp_control_error")
  expect_error(gp_control(backend = "cbc"), class = "gp_control_error")
  expect_error(gp_control(backend = 1L), class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# gp_control(): interval_level + precision_rule + seed validation
# ---------------------------------------------------------------------------

test_that("interval_level must be in the open (0, 1)", {
  expect_error(gp_control(interval_level = 0), class = "gradepath_error")
  expect_error(gp_control(interval_level = 1), class = "gradepath_error")
  expect_error(gp_control(interval_level = 1.5), class = "gradepath_error")
  expect_equal(gp_control(interval_level = 0.95)$interval_level, 0.95)
})

test_that("precision_rule is matched against {none, krw_gmm}", {
  expect_identical(gp_control(precision_rule = "none")$precision_rule, "none")
  expect_identical(gp_control(precision_rule = "krw_gmm")$precision_rule, "krw_gmm")

  # The no-arg call resolves the formal default c("none", "krw_gmm") to "none".
  expect_identical(gp_control()$precision_rule, "none")

  # A bad precision_rule is rejected with the package's classed control error
  # (a subclass of gradepath_error), NOT match.arg()'s base simpleError.
  err <- expect_error(
    gp_control(precision_rule = "log_linear"),
    class = "gp_control_error"
  )
  expect_s3_class(err, "gradepath_error")
  expect_match(conditionMessage(err), "none")
  expect_match(conditionMessage(err), "krw_gmm")

  err2 <- expect_error(
    gp_control(precision_rule = "bad"),
    class = "gp_control_error"
  )
  expect_s3_class(err2, "gradepath_error")

  # tryCatch on the v2-named class catches it (the discriminating contract).
  expect_identical(
    tryCatch(gp_control(precision_rule = "bad"),
             gp_control_error = function(e) "caught"),
    "caught"
  )
})

test_that("seed must be a non-negative integerish value or NULL", {
  expect_identical(gp_control(seed = 7)$seed, 7L)
  expect_identical(gp_control(seed = 0)$seed, 0L)
  expect_null(gp_control(seed = NULL)$seed)
  expect_error(gp_control(seed = -1), class = "gradepath_error")
  expect_error(gp_control(seed = 1.5), class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# gp_control(): solver-option aliases -> solver_options + double-specify
# ---------------------------------------------------------------------------

test_that("time_limit / mip_gap ergonomic aliases land in solver_options", {
  ctl <- gp_control(time_limit = 30, mip_gap = 0.01)
  expect_equal(ctl$solver_options$time_limit, 30)
  expect_equal(ctl$solver_options$mip_gap, 0.01)
})

test_that("solver_options passed directly are validated and preserved", {
  ctl <- gp_control(solver_options = list(time_limit = 10, mip_gap = 0.05, threads = 4))
  expect_equal(ctl$solver_options$time_limit, 10)
  expect_equal(ctl$solver_options$mip_gap, 0.05)
  expect_equal(ctl$solver_options$threads, 4)   # unknown key passes through
})

test_that("HiGHS-native option spellings are accepted and range-checked", {
  ctl <- gp_control(solver_options = list(max_time = 12, mip_rel_gap = 0.02))
  expect_equal(ctl$solver_options$max_time, 12)
  expect_equal(ctl$solver_options$mip_rel_gap, 0.02)
  # Bad HiGHS-native values are rejected.
  expect_error(
    gp_control(solver_options = list(mip_rel_gap = 2)),
    class = "gp_control_error"
  )
})

test_that("double-specifying a knob errors (canonical and HiGHS alias)", {
  # ergonomic alias + same key inside solver_options.
  expect_error(
    gp_control(time_limit = 30, solver_options = list(time_limit = 10)),
    class = "gp_control_error"
  )
  expect_error(
    gp_control(mip_gap = 0.01, solver_options = list(mip_gap = 0.02)),
    class = "gp_control_error"
  )
  # ergonomic alias + the backend-native spelling of the same knob.
  expect_error(
    gp_control(time_limit = 30, solver_options = list(max_time = 10)),
    class = "gp_control_error"
  )
  expect_error(
    gp_control(mip_gap = 0.01, solver_options = list(mip_rel_gap = 0.02)),
    class = "gp_control_error"
  )
})

test_that("bad solver_options shapes are rejected", {
  expect_error(gp_control(solver_options = list(1, 2)), class = "gp_control_error")  # unnamed
  expect_error(gp_control(solver_options = c(1, 2)), class = "gp_control_error")     # not a list
  expect_error(gp_control(time_limit = -5), class = "gradepath_error")              # nonpositive
})

# ---------------------------------------------------------------------------
# validate_gp_control(): round-trips a well-formed object, rejects tampering
# ---------------------------------------------------------------------------

test_that("validate_gp_control() is idempotent and order-sensitive", {
  ctl <- gp_control()
  expect_identical(validate_gp_control(ctl), invisible(ctl))

  # Reordering fields breaks the canonical-order assertion.
  scrambled <- ctl[c("backend", "lambda_grid", "precision_rule", "interval_level",
                     "solver_options", "seed", "schema_version", "provenance")]
  class(scrambled) <- c("gp_control", "list")
  expect_error(validate_gp_control(scrambled), class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# monotone_relabel()
# ---------------------------------------------------------------------------

test_that("monotone_relabel() maps to contiguous 1..k by sorted-unique value", {
  expect_identical(monotone_relabel(c(5, 5, 2, 9)), c(2L, 2L, 1L, 3L))
  expect_identical(monotone_relabel(c(3, 1, 2)), c(3L, 1L, 2L))
  expect_identical(monotone_relabel(c(7, 7, 7)), c(1L, 1L, 1L))
  expect_identical(monotone_relabel(c(2, 1)), c(2L, 1L))
  # Already-canonical input is a fixed point.
  expect_identical(monotone_relabel(c(1L, 2L, 3L)), c(1L, 2L, 3L))
  # Returns integer.
  expect_type(monotone_relabel(c(5, 5, 2, 9)), "integer")
})

test_that("monotone_relabel() rejects NA / empty", {
  expect_error(monotone_relabel(c(1, NA, 2)), class = "gradepath_error")
  expect_error(monotone_relabel(numeric(0)), class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# gp_diag_positions()
# ---------------------------------------------------------------------------

test_that("gp_diag_positions() returns the column-major diagonal indices", {
  expect_identical(gp_diag_positions(3), c(1L, 5L, 9L))
  expect_identical(gp_diag_positions(1), 1L)
  expect_identical(gp_diag_positions(4), c(1L, 6L, 11L, 16L))
  # Cross-check against base which() on a real diagonal mask.
  n <- 5L
  expect_equal(gp_diag_positions(n), which(diag(n) == 1))
  # n == 0 -> empty.
  expect_identical(gp_diag_positions(0), integer(0))
})

test_that("gp_diag_positions() rejects bad n", {
  expect_error(gp_diag_positions(-1), class = "gradepath_error")
  expect_error(gp_diag_positions(2.5), class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# .gp_row_logsumexp()  (and the scalar / normalize siblings)
# ---------------------------------------------------------------------------

test_that(".gp_row_logsumexp() matches log(rowSums(exp(M)))", {
  set.seed(1)
  M <- matrix(rnorm(5 * 4), nrow = 5)
  expect_equal(.gp_row_logsumexp(M), log(rowSums(exp(M))))

  # A single row works (vector promoted to 1-row matrix).
  v <- c(-1, 0, 2, 3)
  expect_equal(.gp_row_logsumexp(matrix(v, nrow = 1)), log(sum(exp(v))))
})

test_that(".gp_row_logsumexp() is overflow-stable and handles -Inf rows", {
  # Large magnitudes: naive exp() overflows to Inf; stable version stays finite.
  big <- matrix(c(1000, 1001, 1002, 1, 2, 3), nrow = 2, byrow = TRUE)
  out <- .gp_row_logsumexp(big)
  expect_true(all(is.finite(out)))
  expect_equal(out[1], 1002 + log(exp(-2) + exp(-1) + exp(0)))

  # An all -Inf row is a clean "zero mass" -> -Inf, not NaN.
  z <- matrix(c(-Inf, -Inf, -Inf, 0, 0, 0), nrow = 2, byrow = TRUE)
  out2 <- .gp_row_logsumexp(z)
  expect_identical(out2[1], -Inf)
  expect_equal(out2[2], log(3))

  # NA is rejected.
  expect_error(.gp_row_logsumexp(matrix(c(1, NA, 2, 3), nrow = 2)), class = "gradepath_error")
})

test_that(".gp_logsumexp() / .gp_softmax() / .gp_safe_normalize() behave", {
  x <- c(-1, 0, 2)
  expect_equal(.gp_logsumexp(x), log(sum(exp(x))))
  expect_equal(sum(.gp_softmax(x)), 1)
  expect_equal(.gp_softmax(x), exp(x) / sum(exp(x)))

  # safe-normalize divides by total; degenerate (all-zero) -> uniform.
  expect_equal(.gp_safe_normalize(c(1, 1, 2)), c(0.25, 0.25, 0.5))
  expect_equal(.gp_safe_normalize(c(0, 0, 0)), rep(1 / 3, 3))
  # negatives rejected.
  expect_error(.gp_safe_normalize(c(1, -1)), class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# provenance helper
# ---------------------------------------------------------------------------

test_that(".gradepath_new_provenance() requires names and drops NULLs", {
  p <- .gradepath_new_provenance(a = 1, b = NULL, c = "x")
  expect_identical(names(p), c("a", "c"))
  expect_error(.gradepath_new_provenance(1), "named")
  expect_identical(.gradepath_new_provenance(), list())
})
