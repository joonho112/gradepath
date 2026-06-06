# =============================================================================
# test-m2-status.R -- M2-only public status surface
# =============================================================================

test_that("gp_m2_status summarizes the current partial M2 status", {
  status <- gp_m2_status()

  expect_s3_class(status, "gp_m2_status")
  expect_s3_class(status, "data.frame")
  expect_identical(
    names(status),
    c("component", "status", "detail", "evidence", "producer_status", "source")
  )
  expect_identical(
    status$component,
    c("formal", "L01_industry_grade_counts", "L02_race_continuous",
      "L02_gender_continuous", "N10_support_guards", "M1_status_boundary")
  )
  expect_identical(attr(status, "m2_formal_status"), "PARTIAL_ACCEPTED")
  expect_false(attr(status, "m1_status_recalculated"))
  expect_s3_class(attr(status, "m2_acceptance"), "gp_m2_acceptance")
  expect_true("gp_m2_status" %in% getNamespaceExports("gradepath"))

  expect_identical(status$status[status$component == "formal"],
                   "PARTIAL_ACCEPTED")
  expect_identical(status$status[status$component == "L01_industry_grade_counts"],
                   "PASS")
  expect_identical(status$status[status$component == "L02_race_continuous"],
                   "APPROXIMATE_OK")
  expect_match(
    status$detail[status$component == "L02_race_continuous"],
    "0\\.0121756787 > 0\\.01"
  )
  expect_identical(status$status[status$component == "L02_gender_continuous"],
                   "PROMOTED")
  expect_match(
    status$detail[status$component == "L02_gender_continuous"],
    "0\\.0066335820 <= 0\\.01"
  )
  expect_match(
    status$evidence[status$component == "L02_gender_continuous"],
    "not direct paper-value reproduction"
  )
  expect_identical(status$status[status$component == "N10_support_guards"],
                   "EVIDENCE_OK")
  expect_match(
    status$evidence[status$component == "N10_support_guards"],
    "synthetic guard evidence"
  )
  expect_identical(status$status[status$component == "M1_status_boundary"],
                   "NOT_RECALCULATED")
  expect_match(
    status$evidence[status$component == "M1_status_boundary"],
    "does not edit or recompute R/status\\.R"
  )
  expect_identical(status$source[status$component == "M1_status_boundary"],
                   "gp_m2_status")
})

test_that("gp_m2_status accepts an explicit acceptance object", {
  gates <- data.frame(
    characteristic = c("race", "gender"),
    pass = c(TRUE, TRUE),
    class_decision = c("banded_candidate", "banded_candidate"),
    producer_status = c("OK", "OK"),
    reason = c("toy_pass", "toy_pass"),
    pi_theta_max_abs = c(0.009, 0.006),
    posteriors_max_abs = c(0.001, 0.001),
    g_theta_support_max_abs = c(1e-5, 1e-5),
    g_theta_density_max_abs = c(1e-4, 1e-4),
    stringsAsFactors = FALSE
  )
  accepted <- gp_m2_acceptance(fixture_gates = gates)
  status <- gp_m2_status(accepted)

  expect_identical(attr(status, "m2_formal_status"), "ACCEPTED")
  expect_false(attr(status, "m1_status_recalculated"))
  expect_identical(status$status[status$component == "formal"], "ACCEPTED")
  expect_identical(status$status[status$component == "L02_race_continuous"],
                   "PROMOTED")
  expect_match(
    status$detail[status$component == "L02_race_continuous"],
    "0\\.0090000000 <= 0\\.01"
  )
  expect_identical(status$status[status$component == "L02_gender_continuous"],
                   "PROMOTED")
})
