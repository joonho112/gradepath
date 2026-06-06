# The runtime solver messages must steer to Gurobi (the strongly
# recommended default backend) and name HiGHS as the explicit license-free
# last-resort fallback. Wording-only change; no behavior change.

test_that("the gurobi-unavailable error recommends Gurobi and names the HiGHS last resort", {
  ns <- asNamespace("gradepath")
  detect <- ns$.gp_gurobi_detect()
  for (k in c("gurobi", "callr", "Matrix", "jsonlite", "gurobi_cl",
              "gurobi_cl_smoke", "binding_smoke")) {
    if (k %in% names(detect)) detect[[k]] <- FALSE
  }
  if ("gurobi_cl_path" %in% names(detect)) detect$gurobi_cl_path <- ""

  msg <- tryCatch(ns$.gp_gurobi_abort_unavailable(detect),
                  error = function(e) conditionMessage(e))

  expect_match(msg, "strongly recommended", fixed = TRUE)
  expect_match(msg, "Gurobi", fixed = TRUE)
  expect_match(msg, "highs", fixed = TRUE)
  expect_match(msg, "gurobi.com/academia")
  expect_match(msg, "last resort", fixed = TRUE)
  # GP-DEC-14-A house style: no claim/tier vocabulary in user-facing output.
  expect_false(grepl("\\btier\\b", msg))
  expect_false(grepl("claim", msg, ignore.case = TRUE))
})

test_that("the .onAttach Gurobi notice uses the recommended-default / last-resort framing", {
  ns <- asNamespace("gradepath")
  src <- paste(deparse(body(ns$.onAttach)), collapse = " ")
  expect_match(src, "strongly recommended", fixed = TRUE)
  expect_match(src, "last resort", fixed = TRUE)
  expect_match(src, "highs", fixed = TRUE)
  expect_match(src, "academia")
})
