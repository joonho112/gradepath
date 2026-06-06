# =============================================================================
# test-m2-acceptance.R -- M2 acceptance and promotion gate
# =============================================================================

.m2_acceptance_file <- function(file) {
  path <- system.file("extdata", "acceptance", file, package = "gradepath")
  if (!nzchar(path) || !file.exists(path)) {
    path <- file.path("inst", "extdata", "acceptance", file)
  }
  path
}

test_that("M2 effective registry promotes L02 only when the panel fixture gate passes", {
  original <- gp_registry[gp_registry$id %in% c(
    "race_industry_dr",
    "race_industry_tau",
    "race_industry_r2",
    "gender_industry_dr",
    "gender_industry_tau",
    "gender_industry_r2"
  ), c("id", "class")]
  expect_true(all(original$class == "approximate"))

  effective <- gp_m2_promoted_registry()
  eff <- effective[match(original$id, effective$id), c("id", "class")]

  expect_identical(eff$class[grepl("^race_", eff$id)],
                   rep("approximate", 3L))
  expect_identical(eff$class[grepl("^gender_", eff$id)],
                   rep("banded", 3L))

  gates <- attr(effective, "m2_fixture_gates")
  expect_equal(nrow(gates), 2L)
  expect_identical(gates$characteristic, c("race", "gender"))
  expect_identical(gates$producer_status, c("APPROXIMATE_OK", "OK"))
})

test_that("M2 recorded fixture snapshot pins measurements", {
  gates <- gp_m2_acceptance()$fixture_gates
  gates <- gates[match(c("race", "gender"), gates$characteristic), ]

  expect_identical(gates$characteristic, c("race", "gender"))
  expect_identical(gates$pass, c(FALSE, TRUE))
  expect_identical(gates$class_decision,
                   c("approximate", "banded_candidate"))
  expect_identical(gates$producer_status, c("APPROXIMATE_OK", "OK"))
  expect_equal(gates$pi_theta_max_abs,
               c(0.0121756787, 0.0066335820), tolerance = 0)
  expect_equal(gates$posteriors_max_abs,
               c(0.0021779164, 0.0020897030), tolerance = 0)
  expect_equal(gates$g_theta_support_max_abs,
               c(0.0001578607, 0.0000888980), tolerance = 0)
  expect_equal(gates$g_theta_density_max_abs,
               c(0.0015353847, 0.0007618801), tolerance = 0)
  expect_gt(gates$pi_theta_max_abs[gates$characteristic == "race"], 0.01)
  expect_lte(gates$pi_theta_max_abs[gates$characteristic == "gender"],
             0.01 + 1e-12)
  expect_lte(max(gates$posteriors_max_abs), 0.01 + 1e-12)
  expect_lte(max(gates$g_theta_support_max_abs), 5e-4 + 1e-12)
  expect_lte(max(gates$g_theta_density_max_abs), 5e-3 + 1e-12)
})

test_that("M2 acceptance scorecard records L01, N10, override, and panel-specific L02 decisions", {
  sc <- gp_m2_acceptance()
  tab <- sc$table

  expect_s3_class(sc, "gp_m2_acceptance")
  expect_identical(
    tab$status[tab$id == "m2_acceptance"],
    "PARTIAL_ACCEPTED"
  )
  expect_identical(
    tab$status[tab$id %in% c("race_industry_ngrades",
                             "gender_industry_ngrades")],
    c("PASS", "PASS")
  )
  expect_identical(
    tab$status[tab$id == "same_industry_override"],
    "EVIDENCE_OK"
  )
  expect_true(all(tab$status[tab$gate == "N10"] == "EVIDENCE_OK"))
  expect_identical(
    tab$status[tab$id %in% c("race_industry_dr",
                             "race_industry_tau",
                             "race_industry_r2")],
    rep("APPROXIMATE_OK", 3L)
  )
  expect_identical(
    tab$effective_class[tab$id %in% c("race_industry_dr",
                                      "race_industry_tau",
                                      "race_industry_r2")],
    rep("approximate", 3L)
  )
  expect_identical(
    tab$status[tab$id %in% c("gender_industry_dr",
                             "gender_industry_tau",
                             "gender_industry_r2")],
    rep("PROMOTED", 3L)
  )
  expect_identical(
    tab$effective_class[tab$id %in% c("gender_industry_dr",
                                      "gender_industry_tau",
                                      "gender_industry_r2")],
    rep("banded", 3L)
  )
  expect_false(any(tab$status == "PASS" & tab$producer_status != "OK"))
  expect_equal(sc$summary$fail, 0)
})

test_that("M2 promotion is data-driven and can promote all panels when both gates pass", {
  gates <- data.frame(
    characteristic = c("race", "gender"),
    pass = c(TRUE, TRUE),
    class_decision = c("banded_candidate", "banded_candidate"),
    producer_status = c("OK", "OK"),
    reason = c("toy_pass", "toy_pass"),
    stringsAsFactors = FALSE
  )
  sc <- gp_m2_acceptance(fixture_gates = gates)
  tab <- sc$table

  expect_identical(tab$status[tab$id == "m2_acceptance"], "ACCEPTED")
  expect_true(all(tab$status[tab$gate == "L02"] == "PROMOTED"))
  expect_true(all(tab$effective_class[tab$gate == "L02"] == "banded"))
  eff <- sc$effective_registry
  expect_true(all(eff$class[eff$id %in% c("race_industry_dr",
                                          "race_industry_tau",
                                          "race_industry_r2",
                                          "gender_industry_dr",
                                          "gender_industry_tau",
                                          "gender_industry_r2")] == "banded"))
})

test_that("M2 acceptance scorecard artifact is explicit and not globally promoted", {
  path <- .m2_acceptance_file("m2-scorecard.csv")
  expect_true(file.exists(path), info = "missing M2 acceptance artifact")
  artifact <- utils::read.csv(path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  required <- c("gate", "layer", "target", "id", "group", "status", "reason",
                "producer_status", "registry_class", "effective_class",
                "paper", "replicated", "delta", "tol", "unit",
                "fixture_decision", "source", "notes")
  expect_true(all(required %in% names(artifact)))

  expect_identical(
    artifact$status[artifact$id == "m2_acceptance"],
    "PARTIAL_ACCEPTED"
  )
  expect_true(all(artifact$status[artifact$gate == "N10"] == "EVIDENCE_OK"))
  expect_true(all(artifact$status[artifact$id %in% c("race_industry_dr",
                                                     "race_industry_tau",
                                                     "race_industry_r2")] ==
                    "APPROXIMATE_OK"))
  expect_true(all(artifact$status[artifact$id %in% c("gender_industry_dr",
                                                     "gender_industry_tau",
                                                     "gender_industry_r2")] ==
                    "PROMOTED"))
  expect_false(all(artifact$effective_class[artifact$gate == "L02"] == "banded"))
  expect_false(any(artifact$status == "PASS" &
                     artifact$producer_status != "OK"))
})
