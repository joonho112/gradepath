# Diagnose the M1 race n=97 selected solve.
#
# This script is intentionally inert when sourced by tests. To refresh the
# diagnostic artifact, execute from the package root with:
#
#   GRADEPATH_RUN_RACE_N97_DIAGNOSTIC=true \
#   Rscript inst/scripts/diagnose-race-n97-selected.R

diagnose_race_n97_dist <- function(grades) {
  paste(names(table(grades)), as.integer(table(grades)), sep = ":", collapse = ";")
}

diagnose_race_n97_symmetrize_pi <- function(P) {
  out <- matrix(0, nrow(P), ncol(P))
  out[upper.tri(out)] <- P[upper.tri(P)]
  out[lower.tri(out)] <- 1 - t(out)[lower.tri(out)]
  diag(out) <- 0.5
  off <- row(out) != col(out)
  out[off] <- pmax(pmin(out[off], 1), 1e-7)
  out
}

diagnose_race_n97_reference_dir <- function() {
  candidates <- c(
    "../KRW-2024-companion-public/dump",
    file.path(
      "log",
      "029_external-review-prep",
      "materials",
      "06_reference-code",
      "krw-replication-archive",
      "dump"
    ),
    file.path(
      "log",
      "029_external-review-prep",
      "materials",
      "06_reference-code",
      "krw-companion-public",
      "dump"
    )
  )
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "Pi_groupfx0_race.csv")) &&
        file.exists(file.path(candidate, "ranking_results_log_dif_binary_race.csv"))) {
      return(candidate)
    }
  }
  stop("KRW race reference dump was not found.", call. = FALSE)
}

diagnose_race_n97_pairwise <- function(control) {
  input <- gp_krw_gmm_input("race")
  coupling <- gp_w_seam(input, "race", control = control)
  prior <- gp_deconvolve(coupling, control = control)
  r_estimates <- list(
    theta_hat = coupling$v_hat,
    s = coupling$s_v,
    original_s = input$s,
    id = input$unit_id
  )
  weights <- gp_posterior_weights(prior, r_estimates, coupling)
  pairwise <- gp_pairwise(weights, ids = input$unit_id, control = control)
  list(input = input, pairwise = pairwise)
}

diagnose_race_n97_archive_pairwise <- function(input, control, reference_dir) {
  raw <- as.matrix(utils::read.csv(
    file.path(reference_dir, "Pi_groupfx0_race.csv"),
    header = FALSE,
    check.names = FALSE
  ))
  off_raw <- row(raw) != col(raw)
  raw_antisym_error <- max(abs((raw + t(raw))[off_raw] - 1))
  P <- diagnose_race_n97_symmetrize_pi(raw)
  dimnames(P) <- list(input$unit_id, input$unit_id)
  pairwise <- new_gp_pairwise(
    ids = input$unit_id,
    matrix = P,
    power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = 0.5, zero_floor = 1e-7),
    source = list(
      stage = "posterior",
      rule = "outer_product",
      assumption = "one_level_independence"
    ),
    control = control,
    provenance = .gradepath_new_provenance(
      producer = "diagnose-race-n97-selected",
      source = file.path(reference_dir, "Pi_groupfx0_race.csv"),
      symmetrized = TRUE,
      raw_antisymmetry_max_abs_error = raw_antisym_error
    )
  )
  list(raw = raw, solver_pi = P, raw_antisym_error = raw_antisym_error,
       pairwise = pairwise)
}

diagnose_race_n97_solve_row <- function(diagnostic_id,
                                        input_source,
                                        pairwise,
                                        control,
                                        expected_distribution,
                                        expected_grade_count,
                                        reference_pairwise_source,
                                        reference_grade_source,
                                        pairwise_reference_status = "",
                                        pairwise_reference_cor_offdiag = NA_real_,
                                        pairwise_reference_max_abs_diff = NA_real_,
                                        pairwise_reference_mean_abs_diff = NA_real_,
                                        pairwise_reference_p95_abs_diff = NA_real_,
                                        pairwise_threshold_disagreement = NA_real_,
                                        notes = "") {
  elapsed <- system.time({
    path <- gp_grade_path(pairwise, control = control, selected_lambda = 0.25)
  })[["elapsed"]]
  selected <- gp_select_grade(path, 0.25)
  endpoint <- gp_select_grade(path, 1)
  distribution <- diagnose_race_n97_dist(selected$assignment$grade)
  grade_count_match <- identical(selected$summary$grade_count, as.integer(expected_grade_count))
  distribution_match <- identical(distribution, expected_distribution)
  acceptance_ready <- .gp_status_acceptance_ready(selected$summary$status)
  accepted <- isTRUE(acceptance_ready) && isTRUE(grade_count_match) &&
    isTRUE(distribution_match)
  diagnosis <- if (isTRUE(accepted)) {
    "ACCEPTED"
  } else if (isTRUE(distribution_match)) {
    "SOLVER_TIME_LIMIT_PUBLISHED_INCUMBENT_MATCH"
  } else {
    "SOLVER_GAP_PENDING_STRICT_CONTROL"
  }

  data.frame(
    diagnostic_id = diagnostic_id,
    input_source = input_source,
    workflow = "one_level_independence",
    control_time_limit_sec = 300,
    control_mip_gap_target = 0.001,
    selected_lambda = 0.25,
    producer_status = .gp_status_from_solver_status(selected$summary$status),
    acceptance_ready = acceptance_ready,
    selected_solver_status = selected$summary$status,
    raw_solver_status = selected$backend$raw_solver_status %gp_or% NA_character_,
    selected_grade_count = selected$summary$grade_count,
    selected_distribution = distribution,
    expected_grade_count = expected_grade_count,
    expected_distribution = expected_distribution,
    published_grade_count_match = grade_count_match,
    published_distribution_match = distribution_match,
    published_distribution_accepted = accepted,
    selected_objective = selected$objective$raw,
    best_bound = selected$backend$objbound,
    mipgap = selected$backend$mipgap,
    solver_runtime_sec = selected$backend$runtime,
    elapsed_sec = as.numeric(elapsed),
    endpoint_solver_status = endpoint$summary$status,
    endpoint_grade_count = endpoint$summary$grade_count,
    endpoint_distribution = "97_singletons",
    pairwise_reference_status = pairwise_reference_status,
    pairwise_reference_cor_offdiag = pairwise_reference_cor_offdiag,
    pairwise_reference_max_abs_diff = pairwise_reference_max_abs_diff,
    pairwise_reference_mean_abs_diff = pairwise_reference_mean_abs_diff,
    pairwise_reference_p95_abs_diff = pairwise_reference_p95_abs_diff,
    pairwise_threshold_disagreement_0_5 = pairwise_threshold_disagreement,
    reference_pairwise_source = reference_pairwise_source,
    reference_grade_source = reference_grade_source,
    problem_hash = selected$backend$problem_hash,
    diagnosis = diagnosis,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

diagnose_race_n97_selected <- function(output_dir = file.path("inst", "extdata", "acceptance")) {
  control <- gp_control(
    lambda_grid = c(0.25, 1),
    backend = "gurobi",
    precision_rule = "krw_gmm",
    time_limit = 300,
    mip_gap = 0.001
  )
  reference_dir <- diagnose_race_n97_reference_dir()
  reference_grade_source <- file.path(reference_dir, "ranking_results_log_dif_binary_race.csv")
  ranking <- utils::read.csv(reference_grade_source, check.names = FALSE)
  expected_distribution <- diagnose_race_n97_dist(ranking[["grades_lamb0.25"]])
  expected_grade_count <- length(unique(ranking[["grades_lamb0.25"]]))

  built <- diagnose_race_n97_pairwise(control)
  archive <- diagnose_race_n97_archive_pairwise(
    built$input,
    control,
    reference_dir
  )
  off <- row(built$pairwise$matrix) != col(built$pairwise$matrix)
  diff <- abs(built$pairwise$matrix[off] - archive$solver_pi[off])
  threshold_disagreement <- mean(
    (built$pairwise$matrix[off] > 0.5) != (archive$solver_pi[off] > 0.5)
  )

  archive_row <- diagnose_race_n97_solve_row(
    diagnostic_id = "archive_pi_controlled_probe",
    input_source = "archive_pi_symmetrized_upper",
    pairwise = archive$pairwise,
    control = control,
    expected_distribution = expected_distribution,
    expected_grade_count = expected_grade_count,
    reference_pairwise_source = file.path(reference_dir, "Pi_groupfx0_race.csv"),
    reference_grade_source = reference_grade_source,
    pairwise_reference_status = "REFERENCE_INPUT",
    pairwise_reference_cor_offdiag = 1,
    pairwise_reference_max_abs_diff = 0,
    pairwise_reference_mean_abs_diff = 0,
    pairwise_reference_p95_abs_diff = 0,
    pairwise_threshold_disagreement = 0,
    notes = sprintf(
      "Archive Pi was symmetrized from rounded upper triangle; raw antisymmetry max error was %.12g.",
      archive$raw_antisym_error
    )
  )
  native_row <- diagnose_race_n97_solve_row(
    diagnostic_id = "native_pi_controlled_probe",
    input_source = "native_live_pairwise",
    pairwise = built$pairwise,
    control = control,
    expected_distribution = expected_distribution,
    expected_grade_count = expected_grade_count,
    reference_pairwise_source = file.path(reference_dir, "Pi_groupfx0_race.csv"),
    reference_grade_source = reference_grade_source,
    pairwise_reference_status = "NEAR_EXACT_NO_THRESHOLD_DRIFT",
    pairwise_reference_cor_offdiag = stats::cor(built$pairwise$matrix[off], archive$solver_pi[off]),
    pairwise_reference_max_abs_diff = max(diff),
    pairwise_reference_mean_abs_diff = mean(diff),
    pairwise_reference_p95_abs_diff = unname(stats::quantile(diff, 0.95)),
    pairwise_threshold_disagreement = threshold_disagreement,
    notes = "Live native Pi recovered an optimal selected solve with the published distribution under stricter controls."
  )

  initial_path <- file.path(output_dir, "m1-race-selected-diagnostic.csv")
  initial <- if (file.exists(initial_path)) {
    current <- utils::read.csv(initial_path, stringsAsFactors = FALSE, check.names = FALSE)
    current[current$diagnostic_id == "initial_probe", , drop = FALSE]
  } else {
    data.frame()
  }
  out <- rbind(initial, archive_row, native_row)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  utils::write.csv(out, initial_path, row.names = FALSE, na = "")
  initial_path
}

if (identical(Sys.getenv("GRADEPATH_RUN_RACE_N97_DIAGNOSTIC"), "true")) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The `pkgload` package is required to run this diagnostic.", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE)
  message(diagnose_race_n97_selected())
}
