.m1_acceptance_file <- function(file) {
  path <- system.file("extdata", "acceptance", file, package = "gradepath")
  if (!nzchar(path) || !file.exists(path)) {
    path <- file.path("inst", "extdata", "acceptance", file)
  }
  path
}

.m1_cached_file <- function(file) {
  path <- system.file("extdata", "cached", file, package = "gradepath")
  if (!nzchar(path) || !file.exists(path)) {
    path <- file.path("inst", "extdata", "cached", file)
  }
  path
}

.m1_read_csv <- function(file) {
  path <- .m1_acceptance_file(file)
  expect_true(file.exists(path), info = paste("missing M1 artifact:", file))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

test_that("M1 acceptance scorecard artifact is explicit and not silently green", {
  scorecard <- .m1_read_csv("m1-scorecard.csv")
  required <- c("gate", "layer", "target", "id", "group", "status", "reason", "notes")
  expect_true(all(required %in% names(scorecard)))
  expect_true(all(c(
    "M1_ACCEPTANCE",
    "K01",
    "K02",
    "K03-K04",
    "K05-K07",
    "K08",
    "K09",
    "K12-K13",
    "K14-K15",
    "NAMES",
    "N9",
    "SOLVER"
  ) %in% scorecard$gate))

  overall <- scorecard[scorecard$gate == "M1_ACCEPTANCE", , drop = FALSE]
  expect_equal(nrow(overall), 1L)
  expect_identical(overall$layer, "formal_gate")
  # the one-level race+gender M1 core is ACCEPTED, backed by
  # proven-optimal (mipgap=0) selected solves at the published distributions. The
  # only non-PASS registered items are DEFERRED two-level/names rows (within-share,
  # names M1_NAMES, industry M2); no UNVERIFIED rows remain.
  expect_identical(overall$status, "ACCEPTED")
  expect_true(any(scorecard$status == "PASS"))
  expect_true(any(scorecard$status == "EVIDENCE_OK"))
  expect_true(any(scorecard$status == "DEFERRED"))
  expect_false(any(scorecard$status == "UNVERIFIED"))
  expect_false(all(scorecard$status == "PASS"))

  pass_rows <- scorecard[scorecard$status == "PASS", , drop = FALSE]
  expect_true(all(pass_rows$layer %in% c("registered_target", "formal_gate", "target_pending_registry")))
  evidence <- scorecard[scorecard$layer %in% c("component_evidence", "diagnostic_evidence"), , drop = FALSE]
  expect_true(nrow(evidence) > 0L)
  expect_true(all(evidence$status == "EVIDENCE_OK"))
  expect_false(any(evidence$status == "PASS"))
  expect_setequal(
    evidence$id,
    c(
      "GP-W-EXACT",
      "Pi_groupfx0",
      "race_n97_controlled_selected",
      "gender_deconv_source_truth",
      "step45_strict3",
      "step45_diagonal"
    )
  )

  pending_registry <- scorecard[scorecard$layer == "target_pending_registry", , drop = FALSE]
  expect_setequal(pending_registry$id, "race_baseline_distribution")
  # proven-optimal incumbent is 2/81/14, matching the published
  # distribution (same proven solve as race_baseline_ngrades).
  expect_true(all(pending_registry$status == "PASS"))
})

test_that("M1 registry scorecard covers every M1 id without false acceptance", {
  skip_if_not_installed("ebrecipe")

  m1_ids <- gp_registry$id[gp_registry$milestone == "M1"]
  race_beta <- unname(gp_w_seam(gp_krw_gmm_input("race"), "race")$beta)
  krw_r2 <- .m1_read_csv("m1-krw-r2-diagnostic.csv")
  krw_r2_m1 <- krw_r2[krw_r2$milestone == "M1" & krw_r2$status == "PASS", , drop = FALSE]
  report_targets <- .m1_read_csv("m1-report-card-targets.csv")
  report_m1 <- report_targets[
    report_targets$milestone == "M1" & report_targets$status == "PASS",
    ,
    drop = FALSE
  ]
  replicated_values <- c(
    nrow(krw_firms),
    length(unique(krw_firms$industry)),
    race_beta,
    4,
    as.numeric(krw_r2_m1$replicated),
    as.numeric(report_m1$replicated)
  )
  replicated <- data.frame(
    id = c(
      "scale_n_firms_graded",
      "race_n_industries",
      "t3_race_ni_beta",
      "n9_strict4_ngrades",
      krw_r2_m1$target_id,
      report_m1$id
    ),
    replicated = replicated_values,
    producer_status = if (is.na(race_beta)) {
      c("OK", "OK", "GMM_NONCONVERGED", "OK",
        rep("OK", nrow(krw_r2_m1) + nrow(report_m1)))
    } else {
      rep("OK", 4 + nrow(krw_r2_m1) + nrow(report_m1))
    },
    group = c(
      "Design",
      "Design",
      "Precision",
      "Solver fixtures",
      rep("KRW R2", nrow(krw_r2_m1)),
      rep("Report card", nrow(report_m1))
    ),
    stringsAsFactors = FALSE
  )
  # the race/gender selected-grade (lambda=0.25) and endpoint
  # (lambda=1) targets are now PROVEN OPTIMAL (mipgap=0) at the published
  # distributions under the predeclared acceptance policy (acceptance_mode=TRUE,
  # time_limit=120, mip_gap=0); see m1-solver-evidence.csv. Their proven verdicts
  # (DR rows in the registry's own units: percent for the selected DR, proportion
  # for dr_l1) are supplied with producer_status="OK" so the registry scorecard
  # scores them. A live-solve re-derivation is the slow-gated test below.
  proven_grades <- data.frame(
    id = c(
      "race_baseline_ngrades", "race_baseline_dr", "race_baseline_tau",
      "race_baseline_worst_n", "race_baseline_best_n",
      "race_baseline_dr_l1", "race_baseline_tau_l1",
      "gender_baseline_ngrades", "gender_baseline_dr", "gender_baseline_tau"
    ),
    replicated = c(
      3, 3.8692210, 0.2069798, 2, 14,
      0.26815392, 0.4636922,
      4, 1.8066945, 0.1208678
    ),
    producer_status = "OK",
    group = "Grades",
    stringsAsFactors = FALSE
  )
  replicated <- rbind(replicated, proven_grades)
  sc <- gp_validate_targets(replicated, targets = m1_ids)

  expect_s3_class(sc, "gp_scorecard")
  expect_setequal(sc$checks$id, m1_ids)
  expect_false(any(sc$checks$status == "FAIL"))
  expect_setequal(
    sc$checks$id[sc$checks$status == "PASS"],
    replicated$id
  )
  expect_setequal(
    sc$checks$id[sc$checks$status == "UNVERIFIED"],
    setdiff(m1_ids, replicated$id)
  )
  expect_true(is.na(sc$pass_rate) || sc$pass_rate == 1)
})

test_that("runtime M1 scorecard upholds the PASS=>OK producer-status invariant", {
  # CCR-13 regression: no PASS row may ever carry a non-OK producer status,
  # asserted on a RUNTIME gp_scorecard object (not just the bundled CSV).
  # Mix acceptance-ready (OK) and not-acceptance-ready (SOLVER_GAP) producers so
  # the comparator must actively suppress PASS for the non-OK rows.
  replicated <- data.frame(
    id = c(
      "scale_n_firms_graded",  # OK + on-target => PASS
      "race_n_industries",     # OK + on-target => PASS
      "race_baseline_ngrades", # non-OK => must NOT be PASS
      "gender_baseline_ngrades" # non-OK => must NOT be PASS
    ),
    replicated = c(97, length(unique(krw_firms$industry)), 3, 4),
    producer_status = c("OK", "OK", "SOLVER_GAP", "SOLVER_TIME_LIMIT"),
    group = "M1",
    stringsAsFactors = FALSE
  )
  sc <- gp_validate_targets(replicated, targets = replicated$id)

  expect_s3_class(sc, "gp_scorecard")
  checks <- sc$checks
  # Core cross-column invariant.
  expect_identical(
    sum(checks$status == "PASS" & checks$producer_status != "OK"),
    0L
  )
  # The OK rows do pass; the non-OK rows are routed to UNVERIFIED.
  expect_true(any(checks$status == "PASS"))
  non_ok <- checks[checks$producer_status != "OK", , drop = FALSE]
  expect_true(nrow(non_ok) > 0L)
  expect_true(all(non_ok$status == "UNVERIFIED"))
})

test_that("KRW R2 diagnostic artifact ports rsquared.do and keeps units explicit", {
  diagnostic <- .m1_read_csv("m1-krw-r2-diagnostic.csv")
  required <- c(
    "diagnostic_id",
    "demographic",
    "model",
    "target_id",
    "quantity",
    "value",
    "producer",
    "source_script",
    "source_script_sha256",
    "source_theta",
    "source_ranking",
    "grade_col",
    "n",
    "grade_count",
    "overall_sd",
    "between_sd",
    "r2_proportion",
    "r2_percent",
    "replicated",
    "paper",
    "delta",
    "tolerance",
    "unit",
    "class",
    "milestone",
    "status",
    "reason",
    "notes"
  )
  expect_true(all(required %in% names(diagnostic)))
  expect_setequal(
    diagnostic$target_id,
    c(
      "race_baseline_r2",
      "race_baseline_betweengrade_sd",
      "gender_baseline_r2",
      "race_industry_r2",
      "gender_industry_r2"
    )
  )
  expect_true(all(diagnostic$producer == "gp_krw_r2"))
  expect_true(all(diagnostic$status == "PASS"))
  expect_true(all(nzchar(diagnostic$source_script_sha256)))

  row <- function(id) diagnostic[diagnostic$target_id == id, , drop = FALSE]
  expect_equal(as.numeric(row("race_baseline_r2")$r2_percent), 24.5972214331222, tolerance = 1e-10)
  expect_equal(as.numeric(row("race_baseline_betweengrade_sd")$between_sd), 0.0339817975178244, tolerance = 1e-12)
  expect_equal(as.numeric(row("gender_baseline_r2")$r2_percent), 43.7177861484058, tolerance = 1e-10)
  expect_equal(as.numeric(row("race_industry_r2")$r2_percent), 70.1403256675552, tolerance = 1e-10)
  expect_equal(as.numeric(row("gender_industry_r2")$r2_percent), 37.7055522895268, tolerance = 1e-10)
  expect_identical(row("race_baseline_r2")$unit, "percent")
  expect_identical(row("race_baseline_betweengrade_sd")$unit, "sd")
  expect_setequal(
    diagnostic$target_id[diagnostic$milestone == "M1"],
    c("race_baseline_r2", "race_baseline_betweengrade_sd", "gender_baseline_r2")
  )
})

test_that("Table F5/F6 report-card target artifact has registered rows and M2 scope", {
  diagnostic <- .m1_read_csv("m1-report-card-targets.csv")
  required <- c(
    "id",
    "source_table",
    "source_table_path",
    "source_table_sha256",
    "firm_pattern",
    "column",
    "producer",
    "replicated",
    "producer_status",
    "group",
    "quantity",
    "paper",
    "delta",
    "tol",
    "unit",
    "class",
    "milestone",
    "status",
    "reason",
    "notes"
  )
  expect_true(all(required %in% names(diagnostic)))
  expect_equal(nrow(diagnostic), 16L)
  expect_true(all(diagnostic$producer == "krw_report_card_table_cells"))
  expect_true(all(diagnostic$producer_status == "OK"))
  expect_true(all(diagnostic$status == "PASS"))
  expect_true(all(nzchar(diagnostic$source_table_sha256)))

  expect_equal(sum(diagnostic$milestone == "M1"), 13L)
  expect_setequal(
    diagnostic$id[diagnostic$milestone == "M2"],
    c(
      "f5_genuineparts_postmean_industry",
      "f5_autonation_condrank_industry",
      "f6_buildersfirstsource_postmean_industry"
    )
  )

  row <- function(id) diagnostic[diagnostic$id == id, , drop = FALSE]
  expect_equal(as.numeric(row("f5_genuineparts_theta")$replicated), 0.3303573, tolerance = 1e-7)
  expect_equal(as.numeric(row("f5_genuineparts_postmean_baseline")$replicated), 0.25011, tolerance = 1e-8)
  expect_equal(as.numeric(row("f5_charterspectrum_condrank_baseline")$replicated), 97)
  expect_equal(as.numeric(row("f6_buildersfirstsource_theta")$replicated), 1.568249, tolerance = 1e-7)
  expect_equal(as.numeric(row("f6_ascena_condrank_baseline")$replicated), 96)
  expect_identical(row("f5_charterspectrum_condrank_baseline")$unit, "count")
  expect_identical(row("f6_buildersfirstsource_theta")$unit, "log_diff")
})

test_that("cached M1 scorecard artifact has build metadata", {
  path <- .m1_cached_file("m1-acceptance-scorecard.rds")
  expect_true(file.exists(path), info = "missing cached M1 scorecard artifact")
  artifact <- readRDS(path)
  meta <- attr(artifact, "build_metadata", exact = TRUE)

  expect_type(artifact, "list")
  expect_s3_class(artifact$scorecard, "gp_scorecard")
  expect_true(all(gp_registry$id[gp_registry$milestone == "M1"] %in%
                    artifact$scorecard$checks$id))
  expect_true(all(c(
    "scorecard",
    "human_scorecard",
    "hard_pass_ids",
    "deferred_ids",
    "component_gates",
    "scale_evidence",
    "race_selected_diagnostic",
    "gender_deconvolution_diagnostic",
    "names_scope_decision",
    "krw_r2_diagnostic",
    "report_card_targets",
    "provenance"
  ) %in% names(artifact)))
  expect_true(is.list(meta))
  expect_true(all(c(
    "key",
    "seed",
    "built_at",
    "r_version",
    "gradepath_version",
    "backend",
    "solver_metadata",
    "source_hash",
    "extra"
  ) %in% names(meta)))
  expect_identical(meta$key, "m1-acceptance-scorecard")
  expect_match(meta$r_version, "^\\d+\\.\\d+\\.\\d+")
  expect_named(meta$solver_metadata, c(
    "backend_env",
    "gurobi_cl_available",
    "gurobi_version",
    "gurobi_r_package_available",
    "highs_package_available",
    "roi_package_available"
  ))
  expect_identical(meta$source_hash, meta$extra$source_hash)
  expect_true("dependency_manifest" %in% names(meta$extra))
  expect_equal(meta$extra$dependency_count, nrow(meta$extra$dependency_manifest))
  expect_true(all(c("path", "category", "reason") %in% names(meta$extra$dependency_manifest)))
  expect_true("layer" %in% names(artifact$human_scorecard))
  expect_true(all(artifact$component_gates$layer == "component_evidence"))
  expect_true(all(artifact$component_gates$status == "EVIDENCE_OK"))
  expect_true(any(artifact$scorecard$checks$status == "UNVERIFIED"))

  # the cached (machine-generated) scorecard must AGREE with
  # the ACCEPTED gate -- the race/gender selected-grade core is live-solved by the
  # builder under the acceptance policy and scores PASS here, not UNVERIFIED.
  core <- c("race_baseline_ngrades", "gender_baseline_ngrades",
            "race_baseline_dr", "race_baseline_tau",
            "gender_baseline_dr", "gender_baseline_tau",
            "race_baseline_worst_n", "race_baseline_best_n",
            "race_baseline_dr_l1", "race_baseline_tau_l1")
  core_rows <- artifact$scorecard$checks[artifact$scorecard$checks$id %in% core, , drop = FALSE]
  expect_setequal(core_rows$id, core)
  expect_true(all(core_rows$status == "PASS"))
  expect_true(all(core_rows$producer_status == "OK"))
  # The ONLY UNVERIFIED M1 id in the cached checks is the two-level within-share.
  uv <- artifact$scorecard$checks$id[artifact$scorecard$checks$status == "UNVERIFIED"]
  expect_setequal(uv, "t3b_race_withinshare")
  # No false PASS in the cached scorecard.
  expect_identical(
    sum(artifact$scorecard$checks$status == "PASS" &
          artifact$scorecard$checks$producer_status != "OK"),
    0L
  )
  # Provenance reflects acceptance, not the pre-repair "not accepted" note.
  expect_match(artifact$provenance$note, "ACCEPTED")

  # freshness guard -- the
  # cached artifact must match the CURRENT source hash (catches a stale build).
  # Resolve the build script via system.file() rather than a wd-relative path:
  # under devtools::test() the working directory is tests/testthat/, so the old
  # "inst/scripts/..." relative path silently no-op'd even in the dev tree. The
  # check now runs wherever the manifest sources are present, and cleanly skips
  # ONLY in a stripped tarball (where data-raw/ is .Rbuildignore'd, so the cache
  # is frozen and its hash cannot be -- and need not be -- recomputed).
  build_script <- system.file("scripts", "build-cached-assets.R",
                              package = "gradepath")
  # A binary install may ship extdata/cached/ but strip inst/scripts/; if the
  # build script is absent there is nothing to recompute against, so skip cleanly
  # (fail-safe) rather than letting sys.source("") error.
  skip_if_not(nzchar(build_script),
              "build script not available (binary install without inst/scripts)")
  bs_env <- new.env(parent = baseenv())
  sys.source(build_script, envir = bs_env)
  # Resolve the source root that holds the most manifest files (mirrors the
  # CCR-02 root-resolution in test-run-all.R), so the hash is computed against
  # the live source regardless of the test working directory.
  manifest_paths <- bs_env$cache_m1_dependency_paths()
  candidate_roots <- c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", ".."),
    file.path(system.file(package = "gradepath"), ".."),
    dirname(dirname(build_script))
  )
  candidate_roots <- unique(candidate_roots[file.exists(candidate_roots)])
  candidate_roots <- normalizePath(candidate_roots, winslash = "/", mustWork = TRUE)
  hit_counts <- vapply(
    candidate_roots,
    function(root) sum(file.exists(file.path(root, manifest_paths))),
    integer(1)
  )
  package_root <- candidate_roots[[which.max(hit_counts)]]
  # Need every manifest file present to recompute the full hash; otherwise this
  # is a stripped tarball with a frozen cache -- skip cleanly.
  all_present <- all(file.exists(file.path(package_root, manifest_paths)))
  skip_if_not(all_present,
              "manifest sources not fully available (frozen tarball cache)")
  expect_identical(meta$source_hash, bs_env$cache_m1_source_hash(root = package_root))
})

test_that("M1 quick gates are wired to the registry comparator", {
  skip_if_not_installed("ebrecipe")

  expect_identical(gp_check("scale_n_firms_graded", 97)$status, "PASS")
  expect_identical(
    gp_check("race_n_industries", length(unique(krw_firms$industry)))$status,
    "PASS"
  )
  expect_identical(gp_check("n9_strict4_ngrades", 4)$status, "PASS")

  race_coupling <- gp_w_seam(gp_krw_gmm_input("race"), "race")
  race_beta <- unname(race_coupling$beta)
  race_beta_check <- gp_check("t3_race_ni_beta", race_beta)
  expect_identical(race_beta_check$status, "PASS")
  expect_lte(abs(race_beta_check$delta), race_beta_check$tol + 1e-9)
})

test_that("n=97 M1 solver evidence records the proven-optimal acceptance solve", {
  evidence <- .m1_read_csv("m1-solver-evidence.csv")
  required <- c(
    "demographic",
    "workflow",
    "acceptance_status",
    "reason",
    "elapsed_sec",
    "backend",
    "path",
    "time_limit_sec",
    "mip_gap_target",
    "selected_lambda",
    "solver_status",
    "incumbent_grade_count",
    "expected_grade_count",
    "incumbent_distribution",
    "expected_distribution",
    "incumbent_objective",
    "best_bound",
    "mipgap",
    "solver_runtime_sec",
    "problem_hash"
  )
  expect_true(all(required %in% names(evidence)))

  # the recorded evidence is now the PROVEN-OPTIMAL solve
  # under the predeclared acceptance policy (acceptance_mode, time_limit=120,
  # mip_gap=0), replacing the earlier gap_reached/time_limit evidence.
  race <- evidence[evidence$demographic == "race", , drop = FALSE]
  expect_equal(nrow(race), 1L)
  expect_identical(race$acceptance_status, "ACCEPTED")
  expect_identical(race$reason, "PROVEN_OPTIMAL")
  expect_identical(race$backend, "gurobi")
  expect_identical(race$path, "gurobi_cl_bigM")
  expect_equal(as.numeric(race$time_limit_sec), 120)
  expect_equal(as.numeric(race$mip_gap_target), 0)
  expect_equal(as.numeric(race$selected_lambda), 0.25)
  expect_identical(race$solver_status, "optimal")
  expect_equal(as.numeric(race$incumbent_grade_count), 3)
  expect_equal(as.numeric(race$expected_grade_count), 3)
  expect_identical(race$incumbent_distribution, "1:2;2:81;3:14")
  expect_identical(race$expected_distribution, "1:2;2:81;3:14")
  expect_equal(as.numeric(race$mipgap), 0)
  expect_true(is.finite(as.numeric(race$incumbent_objective)))
  expect_true(is.finite(as.numeric(race$best_bound)))
  expect_true(is.finite(as.numeric(race$solver_runtime_sec)))
  expect_true(nzchar(race$problem_hash))

  gender <- evidence[evidence$demographic == "gender", , drop = FALSE]
  expect_equal(nrow(gender), 1L)
  expect_identical(gender$acceptance_status, "ACCEPTED")
  expect_identical(gender$reason, "PROVEN_OPTIMAL")
  expect_identical(gender$backend, "gurobi")
  expect_identical(gender$path, "gurobi_cl_bigM")
  expect_equal(as.numeric(gender$time_limit_sec), 120)
  expect_equal(as.numeric(gender$mip_gap_target), 0)
  expect_equal(as.numeric(gender$selected_lambda), 0.25)
  expect_identical(gender$solver_status, "optimal")
  expect_equal(as.numeric(gender$incumbent_grade_count), 4)
  expect_equal(as.numeric(gender$expected_grade_count), 4)
  expect_identical(gender$incumbent_distribution, "1:1;2:3;3:89;4:4")
  expect_identical(gender$expected_distribution, "1:1;2:3;3:89;4:4")
  expect_equal(as.numeric(gender$mipgap), 0)
  expect_true(is.finite(as.numeric(gender$incumbent_objective)))
  expect_true(is.finite(as.numeric(gender$best_bound)))
  expect_true(is.finite(as.numeric(gender$solver_runtime_sec)))
  expect_true(nzchar(gender$problem_hash))
})

test_that("race n=97 selected-solve diagnostic separates controls from upstream drift", {
  diagnostic <- .m1_read_csv("m1-race-selected-diagnostic.csv")
  required <- c(
    "diagnostic_id",
    "input_source",
    "workflow",
    "control_time_limit_sec",
    "control_mip_gap_target",
    "selected_lambda",
    "producer_status",
    "acceptance_ready",
    "selected_solver_status",
    "selected_grade_count",
    "selected_distribution",
    "expected_grade_count",
    "expected_distribution",
    "published_grade_count_match",
    "published_distribution_match",
    "published_distribution_accepted",
    "mipgap",
    "solver_runtime_sec",
    "pairwise_reference_status",
    "pairwise_reference_cor_offdiag",
    "pairwise_reference_max_abs_diff",
    "pairwise_threshold_disagreement_0_5",
    "diagnosis"
  )
  expect_true(all(required %in% names(diagnostic)))
  expect_setequal(
    diagnostic$diagnostic_id,
    c(
      "initial_probe",
      "archive_pi_controlled_probe",
      "native_pi_controlled_probe"
    )
  )

  initial <- diagnostic[diagnostic$diagnostic_id == "initial_probe", , drop = FALSE]
  expect_identical(initial$producer_status, "SOLVER_GAP")
  expect_identical(initial$selected_solver_status, "gap_reached")
  expect_equal(as.numeric(initial$selected_grade_count), 6)
  expect_identical(initial$selected_distribution, "1:1;2:13;3:40;4:34;5:5;6:4")
  expect_false(as.logical(initial$published_distribution_accepted))

  archive <- diagnostic[diagnostic$diagnostic_id == "archive_pi_controlled_probe", , drop = FALSE]
  expect_identical(archive$selected_solver_status, "time_limit")
  expect_equal(as.numeric(archive$selected_grade_count), 3)
  expect_identical(archive$selected_distribution, "1:2;2:81;3:14")
  expect_true(as.logical(archive$published_distribution_match))
  expect_false(as.logical(archive$published_distribution_accepted))

  native <- diagnostic[diagnostic$diagnostic_id == "native_pi_controlled_probe", , drop = FALSE]
  expect_identical(native$producer_status, "OK")
  expect_identical(native$selected_solver_status, "optimal")
  expect_identical(native$selected_distribution, "1:2;2:81;3:14")
  expect_true(as.logical(native$published_distribution_accepted))
  expect_identical(native$diagnosis, "ACCEPTED")
  expect_gt(as.numeric(native$pairwise_reference_cor_offdiag), 0.999)
  expect_lt(as.numeric(native$pairwise_reference_max_abs_diff), 0.001)
  expect_equal(as.numeric(native$pairwise_threshold_disagreement_0_5), 0)
})

test_that("gender deconvolution diagnostic separates source-truth input from bundled example", {
  diagnostic <- .m1_read_csv("m1-gender-deconvolution-diagnostic.csv")
  required <- c(
    "diagnostic_id",
    "input_source",
    "source_role",
    "producer_status",
    "deconv_status",
    "beta",
    "m_hat_sigma_xi",
    "caps_lo",
    "caps_hi",
    "deconv_penalty",
    "solver_status",
    "selected_grade_count",
    "selected_distribution",
    "expected_grade_count",
    "expected_distribution",
    "grade_count_match",
    "distribution_match",
    "accepted",
    "diagnosis"
  )
  expect_true(all(required %in% names(diagnostic)))
  expect_setequal(
    diagnostic$diagnostic_id,
    c("krw_gmm_input_controlled", "krw_firms_public_example")
  )

  source_truth <- diagnostic[
    diagnostic$diagnostic_id == "krw_gmm_input_controlled",
    ,
    drop = FALSE
  ]
  expect_identical(source_truth$source_role, "source_truth")
  expect_identical(source_truth$deconv_status, "OK")
  expect_identical(source_truth$producer_status, "SOLVER_TIME_LIMIT")
  expect_equal(as.numeric(source_truth$beta), 1.255381, tolerance = 1e-6)
  expect_equal(as.numeric(source_truth$m_hat_sigma_xi), 1.233586, tolerance = 1e-6)
  expect_equal(as.numeric(source_truth$caps_lo), -6.16793, tolerance = 1e-5)
  expect_equal(as.numeric(source_truth$caps_hi), 6.16793, tolerance = 1e-5)
  expect_true(is.finite(as.numeric(source_truth$deconv_penalty)))
  expect_identical(source_truth$solver_status, "time_limit")
  expect_equal(as.numeric(source_truth$selected_grade_count), 4)
  expect_identical(source_truth$selected_distribution, "1:1;2:3;3:89;4:4")
  expect_equal(as.numeric(source_truth$expected_grade_count), 4)
  expect_identical(source_truth$expected_distribution, "1:1;2:3;3:89;4:4")
  expect_true(as.logical(source_truth$grade_count_match))
  expect_true(as.logical(source_truth$distribution_match))
  expect_false(as.logical(source_truth$accepted))
  expect_identical(
    source_truth$diagnosis,
    "DECONV_FIXED_SOLVER_CERTIFICATE_PENDING"
  )

  public_example <- diagnostic[
    diagnostic$diagnostic_id == "krw_firms_public_example",
    ,
    drop = FALSE
  ]
  expect_identical(public_example$source_role, "bundled_example")
  expect_identical(public_example$deconv_status, "FAILED")
  expect_identical(public_example$producer_status, "DECONV_BOUNDARY_ERROR")
  expect_gt(as.numeric(public_example$beta), 2.9)
  expect_gt(as.numeric(public_example$m_hat_sigma_xi), 800)
  expect_lt(as.numeric(public_example$caps_lo), -4000)
  expect_gt(as.numeric(public_example$caps_hi), 6000)
  expect_identical(
    public_example$diagnosis,
    "INPUT_SOURCE_MISMATCH_BOUNDARY_REPRODUCER"
  )
})

test_that("names pipeline is explicitly deferred out of the current M1 gate", {
  scope <- .m1_read_csv("m1-names-scope.csv")
  required <- c("record_type", "name", "status", "path_or_ids", "role", "decision", "notes")
  expect_true(all(required %in% names(scope)))

  decision <- scope[scope$record_type == "decision", , drop = FALSE]
  expect_equal(nrow(decision), 1L)
  expect_identical(decision$name, "names_pipeline")
  expect_identical(decision$status, "DEFERRED_TO_M1_NAMES")
  expect_identical(decision$decision, "DEFER")

  names_ids <- gp_registry$id[
    grepl("^names_", gp_registry$id) | gp_registry$id == "scale_n_names"
  ]
  expect_length(names_ids, 17L)
  expect_true(all(gp_registry$milestone[match(names_ids, gp_registry$id)] == "M1_NAMES"))
  expect_false(any(names_ids %in% gp_registry$id[gp_registry$milestone == "M1"]))

  scorecard <- .m1_read_csv("m1-scorecard.csv")
  names_row <- scorecard[scorecard$id == "names_ngrades_l025", , drop = FALSE]
  expect_equal(nrow(names_row), 1L)
  expect_identical(names_row$gate, "NAMES")
  expect_identical(names_row$layer, "deferred_milestone")
  expect_identical(names_row$status, "DEFERRED")
  expect_identical(names_row$reason, "M1_NAMES_SCOPE")
})

test_that("M1 race+gender targets are accepted; two-level/names items stay deferred", {
  scorecard <- .m1_read_csv("m1-scorecard.csv")

  row <- function(id) scorecard[scorecard$id == id, , drop = FALSE]

  # the race/gender selected-grade targets are now PROVEN
  # OPTIMAL (mipgap=0) at the published distributions, so they are accepted PASSes.
  expect_identical(row("race_baseline_ngrades")$status, "PASS")
  expect_identical(row("race_baseline_ngrades")$layer, "registered_target")
  expect_identical(row("gender_baseline_ngrades")$status, "PASS")
  expect_identical(row("race_baseline_worst_n")$status, "PASS")
  expect_identical(row("race_baseline_best_n")$status, "PASS")
  expect_identical(row("race_baseline_dr")$status, "PASS")
  # Within-industry share is a TWO-LEVEL variance decomposition; deferred (not a
  # one-level blocker), like the industry rows.
  expect_identical(row("t3b_race_withinshare")$status, "DEFERRED")
  expect_identical(row("t3b_race_withinshare")$reason, "TWO_LEVEL_WITHIN_SHARE")
  # Names + industry/two-level report cards remain separate deferred milestones.
  expect_identical(row("names_ngrades_l025")$reason, "M1_NAMES_SCOPE")
  expect_identical(row("names_ngrades_l025")$layer, "deferred_milestone")
  expect_identical(row("race_baseline_r2")$status, "PASS")
  expect_identical(row("race_baseline_r2")$reason, "")
  expect_identical(row("race_baseline_betweengrade_sd")$status, "PASS")
  expect_identical(row("gender_baseline_r2")$status, "PASS")
  expect_identical(row("f5_genuineparts_theta")$status, "PASS")
  expect_identical(row("f5_charterspectrum_condrank_baseline")$status, "PASS")
  expect_identical(row("f6_buildersfirstsource_theta")$status, "PASS")
  expect_identical(row("f6_ascena_condrank_baseline")$status, "PASS")
  expect_identical(
    row("f5_f6_industry_report_card_cells")$layer,
    "deferred_milestone"
  )
  expect_identical(
    row("f5_f6_industry_report_card_cells")$reason,
    "M2_INDUSTRY_SCOPE"
  )
})

test_that("M1 race+gender selected solves prove optimal at the published distributions (slow, acceptance-mode)", {
  skip_if_not_installed("ebrecipe")
  skip_if_not_installed("gurobi")
  if (!nzchar(Sys.getenv("GRADEPATH_RUN_SLOW_TESTS"))) {
    skip("slow n=97 acceptance solve (set GRADEPATH_RUN_SLOW_TESTS=1 to run)")
  }
  # the live proof behind the ACCEPTED M1 one-level gate.
  # Under the predeclared acceptance policy the selected solves must prove OPTIMAL
  # (mipgap 0) at the published distributions, with producer_status OK and every
  # auto-extracted registry target passing.
  ctrl <- gp_control(lambda_grid = c(0.25, 1), backend = "gurobi",
                     precision_rule = "krw_gmm", time_limit = 120, mip_gap = 0)
  expected <- list(race = c(2L, 81L, 14L), gender = c(1L, 3L, 89L, 4L))
  for (dem in c("race", "gender")) {
    inp <- gp_krw_gmm_input(dem)
    dat <- data.frame(theta_hat = inp$theta_hat, s = inp$s, unit_id = inp$unit_id,
                      stringsAsFactors = FALSE)
    fit <- krw_report_card(dat, demographic = dem, control = ctrl, lambda = 0.25,
                           acceptance_mode = TRUE)
    expect_identical(fit$provenance$producer_status, "OK")
    expect_identical(fit$selected_grade$backend$status, "optimal")
    expect_lt(fit$selected_grade$backend$mipgap, 1e-6)
    dist <- as.integer(table(fit$selected_grade$assignment$grade))
    expect_identical(sort(dist, decreasing = TRUE),
                     sort(expected[[dem]], decreasing = TRUE))
    sc <- gp_validate_targets(fit)
    expect_false(any(sc$checks$status == "FAIL"))
    expect_true(all(sc$checks$status[sc$checks$status != "n/a"] == "PASS"))
    expect_true(all(sc$checks$producer_status == "OK"))
  }
})
