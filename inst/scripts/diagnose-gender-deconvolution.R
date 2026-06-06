# Diagnose the M1 gender deconvolution boundary failure.
#
# This script is intentionally inert when sourced by tests. To refresh the
# diagnostic artifact, execute from the package root with:
#
#   GRADEPATH_RUN_GENDER_DECONV_DIAGNOSTIC=true \
#   Rscript inst/scripts/diagnose-gender-deconvolution.R

diagnose_gender_dist <- function(grades) {
  paste(names(table(grades)), as.integer(table(grades)), sep = ":", collapse = ";")
}

diagnose_gender_reference_dir <- function() {
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
    if (file.exists(file.path(candidate, "ranking_results_log_dif_binary_gender.csv"))) {
      return(candidate)
    }
  }
  stop("KRW gender reference dump was not found.", call. = FALSE)
}

diagnose_gender_expected <- function() {
  reference_dir <- diagnose_gender_reference_dir()
  ranking <- utils::read.csv(
    file.path(reference_dir, "ranking_results_log_dif_binary_gender.csv"),
    check.names = FALSE
  )
  grades <- ranking[["grades_lamb0.25"]]
  list(
    source = file.path(reference_dir, "ranking_results_log_dif_binary_gender.csv"),
    grade_count = length(unique(grades)),
    distribution = diagnose_gender_dist(grades)
  )
}

diagnose_gender_fit_summary <- function(data, source_role, input_source, control) {
  api_input <- .gp_api_estimates(data, "gender")
  fit <- gp_w_seam(api_input$estimates, "gender", control = control)
  list(
    api_input = api_input,
    fit = fit,
    source_role = source_role,
    input_source = input_source
  )
}

diagnose_gender_deconvolve <- function(fit, control) {
  elapsed <- system.time({
    prior <- tryCatch(gp_deconvolve(fit, control = control), error = function(e) e)
  })[["elapsed"]]
  list(prior = prior, elapsed = as.numeric(elapsed))
}

diagnose_gender_controlled_row <- function(expected, control) {
  data <- gp_krw_gmm_input("gender")
  prepared <- diagnose_gender_fit_summary(
    data,
    source_role = "source_truth",
    input_source = "inst/extdata/krw-gmm-input/theta_estimates_matlab_gender.csv",
    control = control
  )
  deconv <- diagnose_gender_deconvolve(prepared$fit, control)
  if (inherits(deconv$prior, "error")) {
    stop("Source-truth gender deconvolution unexpectedly failed.", call. = FALSE)
  }

  elapsed <- system.time({
    full_fit <- krw_report_card(data, "gender", control = control, lambda = 0.25)
  })[["elapsed"]]
  selected <- full_fit$selected_grade
  distribution <- diagnose_gender_dist(selected$assignment$grade)
  grade_count_match <- identical(selected$summary$grade_count, as.integer(expected$grade_count))
  distribution_match <- identical(distribution, expected$distribution)
  accepted <- isTRUE(.gp_status_acceptance_ready(selected$summary$status)) &&
    isTRUE(grade_count_match) &&
    isTRUE(distribution_match)
  f <- prepared$fit
  prior <- deconv$prior

  data.frame(
    diagnostic_id = "krw_gmm_input_controlled",
    input_source = prepared$input_source,
    workflow = "one_level_independence",
    source_role = prepared$source_role,
    producer_status = .gp_producer_status_from_selected_grade(selected),
    deconv_status = "OK",
    deconv_message = "",
    beta = f$beta,
    m_hat_sigma_xi = f$m_hat[1L],
    V_m = f$V_m[1L, 1L],
    caps_lo = f$caps[["lo"]],
    caps_hi = f$caps[["hi"]],
    v_min = min(f$v_hat),
    v_max = max(f$v_hat),
    s_min = min(f$s_v),
    s_max = max(f$s_v),
    theta_min = min(prepared$api_input$theta_hat),
    theta_max = max(prepared$api_input$theta_hat),
    se_min = min(prepared$api_input$s),
    se_max = max(prepared$api_input$s),
    deconv_penalty = prior$diagnostics$penalty,
    deconv_J = prior$diagnostics$J,
    deconv_converged = prior$diagnostics$converged,
    deconv_elapsed_sec = deconv$elapsed,
    selected_lambda = 0.25,
    solver_status = selected$summary$status,
    raw_solver_status = selected$backend$raw_solver_status %gp_or% NA_character_,
    selected_grade_count = selected$summary$grade_count,
    selected_distribution = distribution,
    expected_grade_count = expected$grade_count,
    expected_distribution = expected$distribution,
    grade_count_match = grade_count_match,
    distribution_match = distribution_match,
    accepted = accepted,
    selected_objective = selected$objective$raw,
    best_bound = selected$backend$objbound,
    mipgap = selected$backend$mipgap,
    solver_runtime_sec = selected$backend$runtime,
    elapsed_sec = as.numeric(elapsed),
    problem_hash = selected$backend$problem_hash,
    diagnosis = if (isTRUE(accepted)) {
      "ACCEPTED"
    } else {
      "DECONV_FIXED_SOLVER_CERTIFICATE_PENDING"
    },
    notes = paste(
      "KRW Matlab GMM input deconvolves successfully and matches the",
      "published gender lambda 0.25 distribution but the selected solve",
      "stopped before an acceptance-ready solver certificate."
    ),
    stringsAsFactors = FALSE
  )
}

diagnose_gender_bad_input_row <- function(expected, control) {
  prepared <- diagnose_gender_fit_summary(
    krw_firms,
    source_role = "bundled_example",
    input_source = "data/krw_firms.rda",
    control = control
  )
  deconv <- diagnose_gender_deconvolve(prepared$fit, control)
  failed <- inherits(deconv$prior, "error")
  f <- prepared$fit

  data.frame(
    diagnostic_id = "krw_firms_public_example",
    input_source = prepared$input_source,
    workflow = "one_level_independence",
    source_role = prepared$source_role,
    producer_status = if (failed) "DECONV_BOUNDARY_ERROR" else "OK",
    deconv_status = if (failed) "FAILED" else "OK",
    deconv_message = if (failed) conditionMessage(deconv$prior) else "",
    beta = f$beta,
    m_hat_sigma_xi = f$m_hat[1L],
    V_m = f$V_m[1L, 1L],
    caps_lo = f$caps[["lo"]],
    caps_hi = f$caps[["hi"]],
    v_min = min(f$v_hat),
    v_max = max(f$v_hat),
    s_min = min(f$s_v),
    s_max = max(f$s_v),
    theta_min = min(prepared$api_input$theta_hat),
    theta_max = max(prepared$api_input$theta_hat),
    se_min = min(prepared$api_input$s),
    se_max = max(prepared$api_input$s),
    deconv_penalty = NA_real_,
    deconv_J = NA_real_,
    deconv_converged = NA,
    deconv_elapsed_sec = deconv$elapsed,
    selected_lambda = 0.25,
    solver_status = "",
    raw_solver_status = "",
    selected_grade_count = NA_integer_,
    selected_distribution = "",
    expected_grade_count = expected$grade_count,
    expected_distribution = expected$distribution,
    grade_count_match = FALSE,
    distribution_match = FALSE,
    accepted = FALSE,
    selected_objective = NA_real_,
    best_bound = NA_real_,
    mipgap = NA_real_,
    solver_runtime_sec = NA_real_,
    elapsed_sec = deconv$elapsed,
    problem_hash = "",
    diagnosis = "INPUT_SOURCE_MISMATCH_BOUNDARY_REPRODUCER",
    notes = paste(
      "Bundled krw_firms is a public example input and not the KRW Matlab GMM",
      "source-truth series used for M1 precision parity."
    ),
    stringsAsFactors = FALSE
  )
}

diagnose_gender_deconvolution <- function(output_dir = file.path("inst", "extdata", "acceptance")) {
  control <- gp_control(
    lambda_grid = c(0.25, 1),
    backend = "gurobi",
    precision_rule = "krw_gmm",
    time_limit = 60,
    mip_gap = 0.01
  )
  expected <- diagnose_gender_expected()
  out <- rbind(
    diagnose_gender_controlled_row(expected, control),
    diagnose_gender_bad_input_row(expected, control)
  )
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  path <- file.path(output_dir, "m1-gender-deconvolution-diagnostic.csv")
  utils::write.csv(out, path, row.names = FALSE, na = "")
  path
}

if (identical(Sys.getenv("GRADEPATH_RUN_GENDER_DECONV_DIAGNOSTIC"), "true")) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The `pkgload` package is required to run this diagnostic.", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE)
  message(diagnose_gender_deconvolution())
}
