# gp_calibrate() budget heads-up. gp_calibrate re-solves a grade
# integer program once per draw, so a large n_sim and/or fit is slow; the guard
# warns before a likely-slow run and stays quiet for small/quick checks. It is a
# pure informational message (suppressible) and never changes the result.

test_that("the budget message fires for a likely-slow run and is quiet for small ones", {
  bm <- gradepath:::.gp_cal_budget_message

  # Quiet for quick checks (small n_sim and/or small fit), and for bad input.
  expect_null(bm(2L, 24L))
  expect_null(bm(10L, 24L))
  expect_null(bm(20L, 24L))            # 480 <= 600
  expect_null(bm(5L, 97L))             # 485 <= 600
  expect_null(bm(200L, NA_real_))      # non-finite n_units -> no message

  # Fires for a genuinely large workload.
  expect_type(bm(200L, 97L), "character")   # the vignette's old default: hours
  expect_type(bm(8L, 97L), "character")      # a large fit even at modest n_sim

  msg <- bm(200L, 97L)
  expect_match(msg, "integer program")
  expect_match(msg, "200 refit\\+grade solves on 97 units")
  expect_match(msg, "gp_preview")
  # GP-DEC-14-A house style.
  expect_false(grepl("\\btier\\b", msg))
  expect_false(grepl("claim", msg, ignore.case = TRUE))
})
