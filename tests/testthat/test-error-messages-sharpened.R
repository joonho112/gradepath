# Sharpened runtime error messages (issues 05 + 04, runtime halves).

test_that("gp_report_card() names the missing stage-wise argument(s) + a recipe", {
  # All three absent -> names all three + the fit-slots recipe.
  e <- tryCatch(gp_report_card(), error = function(err) conditionMessage(err))
  expect_match(e, "missing required argument")
  expect_match(e, "posterior")
  expect_match(e, "selected_grade")
  expect_match(e, "grade_path")
  expect_match(e, "get_posterior(fit)", fixed = TRUE)
  expect_false(grepl("\\btier\\b", e))
  expect_false(grepl("claim", e, ignore.case = TRUE))
})

test_that("gp_twolevel_report_card() explains the posterior-backed / M2 scope", {
  skip_if(!requireNamespace("highs", quietly = TRUE), "highs backend unavailable")
  ctl <- gp_control(backend = "highs")
  strict <- function(n, p) {
    ids <- paste0(p, seq_len(n)); m <- matrix(0, n, n, dimnames = list(ids, ids))
    for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) m[i, j] <- if (i < j) 0.95 else 0.05
    m
  }
  tlb <- gp_twolevel_pairwise(strict(4, "f"), strict(3, "i"), control = ctl)
  tlg <- gp_twolevel_grade(tlb, control = ctl, lambda_grid = c(0.25, 1))
  e <- tryCatch(gp_twolevel_report_card(tlg, model = "industry_rfe"),
                error = function(err) conditionMessage(err))
  expect_match(e, "posterior-backed")
  expect_match(e, "PARTIAL_ACCEPTED")
  expect_match(e, "gp_m2_status")
  expect_false(grepl("quadrature", e))   # no internal function name leaked
  expect_false(grepl("\\btier\\b", e))
})
