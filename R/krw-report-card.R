# =============================================================================
# krw-report-card.R -- one-call monolith
# -----------------------------------------------------------------------------
# Runs the M1 one-level KRW chain:
#   eb_input -> beta-GMM W seam -> native deconvolution -> W -> posterior ->
#   gp_pairwise -> gp_grade_path -> gp_select_grade -> gp_report_card -> gp_fit.
# =============================================================================

.gp_api_with_status <- function(status, expr) {
  tryCatch(
    expr,
    error = function(e) {
      mapped <- .gp_status_from_condition(e)
      if (mapped %in% c(
        "INPUT_ERROR",
        "GROUPS_ERROR",
        "SOLVER_BACKEND_UNAVAILABLE",
        "SOLVER_INFEASIBLE",
        "SOLVER_TIME_LIMIT",
        "SOLVER_GAP",
        "SOLVER_OUTPUT_INVALID",
        "SOLVER_OBJECTIVE_MISMATCH",
        "SOLVER_CANONICAL_MISMATCH"
      )) {
        status <- mapped
      }
      .gp_status_abort(
        status,
        "%s",
        conditionMessage(e),
        class = setdiff(class(e), c("error", "condition"))
      )
    }
  )
}

.gp_precision_fit_from_coupling <- function(coupling) {
  validate_gp_estimation_fit(coupling)
  model_form <- if (identical(coupling$characteristic, "gender")) {
    "additive"
  } else {
    "multiplicative"
  }
  mu <- coupling$provenance$mu %gp_or% coupling$report$E_theta
  if (is.null(mu) && identical(model_form, "multiplicative")) {
    mu <- 0
  }

  validate_gp_precision_fit(new_gp_precision_fit(
    parameters = list(
      model_form = model_form,
      beta = as.numeric(coupling$beta),
      mu = as.numeric(mu),
      characteristic = coupling$characteristic
    ),
    moments = list(
      m_hat = coupling$m_hat,
      V_m = coupling$V_m,
      J = coupling$J,
      df = coupling$df,
      p_value = coupling$p_value
    ),
    diagnostics = list(
      report = coupling$report,
      caps = coupling$caps
    ),
    scale = "r",
    provenance = .gradepath_new_provenance(
      producer = "krw_report_card",
      step = "precision-fit-shell",
      source = "gp_estimation_fit"
    )
  ))
}

.gp_r_estimates_from_coupling <- function(coupling, input) {
  list(
    theta_hat = coupling$v_hat,
    s = coupling$s_v,
    original_s = input$original_s,
    id = input$ids,
    label = input$label
  )
}

#' Run the one-level KRW report-card pipeline in a single call
#'
#' `krw_report_card()` is the one-call entry point for the M1 one-level
#' Kline-Rose-Walters grading workflow. Give it firm-level estimates and a
#' demographic; it builds the stage-1 input through the ebrecipe seam, runs
#' gradepath's native one-level KRW core (beta-GMM precision seam -> deconvolution
#' -> posterior weights), then chains the public social-choice verbs
#' ([gp_pairwise()] -> [gp_grade_path()] -> [gp_select_grade()] ->
#' [gp_report_card()]) into one validated `gp_fit` you can print, plot, and pass to
#' the accessors. `gradepath()` is an alias.
#'
#' @param data Data frame, list, or `ebrecipe::eb_estimates` carrying `theta_hat`/`s`
#'   (or the demographic-specific `theta_hat_race`/`se_race`,
#'   `theta_hat_gender`/`se_gender`) plus an optional `unit_id`/`firm_id`/`label`. Do
#'   not supply both a generic and a demographic-specific column for one quantity.
#' @param demographic Character scalar, `"race"` or `"gender"`; selects the estimate
#'   columns and the precision-dependence form (multiplicative for race, additive for
#'   gender). Defaults to `"race"`.
#' @param groups Optional grouping vector. Grouped/two-level execution is not
#'   implemented in this monolith; pass `NULL` (the default). A non-`NULL` value
#'   raises a `GROUPS_ERROR`. See [gp_twolevel_report_card()] for the two-level path.
#' @param control A [gp_control] object. The monolith requires
#'   `precision_rule = "krw_gmm"` and defaults to that rule.
#' @param lambda Numeric scalar in `[0, 1]`; the frontier penalty at which the final
#'   grade is selected from the path. Defaults to `0.25` (KRW's published selection).
#' @param acceptance_mode Logical; solver solution-quality policy threaded to the
#'   internal [gp_grade_path()] (gurobi backend only; ignored for open backends).
#'   `FALSE` (default) reports a `gap_reached`/`time_limit` solve honestly with that
#'   status (it lands `UNVERIFIED` downstream); `TRUE` additionally attempts later
#'   Gurobi paths to prove the optimum and never relabels a non-optimal solve as
#'   `optimal`.
#' @param ... Unused; an error is raised if any argument is passed.
#'
#' @return A validated `gp_fit` object (a list of class `c("gp_fit", "list")`) with
#'   the public slots: \describe{
#'   \item{`ids`}{Character vector; the one canonical unit-id order shared by every
#'     downstream slot (an enforced invariant).}
#'   \item{`estimates`}{The stage-1 `ebrecipe::eb_estimates` input.}
#'   \item{`prior`}{A `gp_prior`; the native deconvolved prior on the r-scale.}
#'   \item{`posterior`}{A `gp_posterior`; per-unit posterior summaries
#'     (`posterior_mean`, `lower`, `upper`, `estimate`, `se`).}
#'   \item{`precision_fit`}{A `gp_precision_fit`; the beta-GMM precision fit and its
#'     J-statistic moments.}
#'   \item{`pairwise`}{A `gp_pairwise`; the J x J posterior outranking structure.}
#'   \item{`grade_path`}{A `gp_grade_path`; the solved frontier of grade assignments
#'     across the penalty grid, with per-solve solver status.}
#'   \item{`selected_grade`}{A `gp_grade_fit`; the assignment at the selected
#'     `lambda` (exactly the `grade_path` member at that penalty).}
#'   \item{`report_card`}{A `gp_report_card`; the per-unit grade-label table with
#'     posterior summaries.}
#'   \item{`control`}{The validated [gp_control] used for the run.}
#'   \item{`provenance`}{Named list: producer, demographic, workflow, and the honest
#'     `producer_status` (e.g. `ACCEPTED` / `UNVERIFIED`).}
#'   \item{`warnings`, `schema_version`}{Internal audit slots.}
#' }
#'
#' @details
#' The chain runs in order: stage-1 input (ebrecipe seam) -> beta-GMM precision seam
#' -> native deconvolution -> posterior weights -> one-level posterior ->
#' [gp_pairwise()] -> [gp_grade_path()] -> [gp_select_grade()] -> [gp_report_card()].
#' Grade labels are integers in `{1, ..., n}` and carry no ranking-superiority
#' statement of any kind. Solver honesty: a gurobi solve that stops at a gap or time
#' limit keeps that status (surfaced as an `UNVERIFIED` producer status);
#' `acceptance_mode = TRUE` only adds optimization attempts, never a relabel.
#'
#' @note Replicating KRW: do NOT pass the bundled [krw_firms] dataset here. It is a
#'   public example on a different numeric scale and will NOT reproduce the published
#'   KRW (2024) results -- the beta-GMM lands on a spurious large-beta optimum (race
#'   beta ~ 2.1 / gender beta ~ 3.0 versus the published 0.51 / 1.26) and the gender
#'   path errors in deconvolution with `DECONV_BOUNDARY_ERROR`. For replication read
#'   KRW's real Matlab GMM series (shipped under `inst/extdata/krw-gmm-input/`) directly,
#'   e.g. `read.csv(system.file("extdata/krw-gmm-input/theta_estimates_matlab_race.csv",
#'   package = "gradepath"), header = FALSE)` (column 2 = `theta_hat`, column 3 = `s`);
#'   see [krw_firms] and the applied vignettes.
#'
#' @examples
#' # A small example input bundled with the package (24 firms; an example subset,
#' # not the parity input). The runnable solve uses the open HiGHS backend.
#' inp <- readRDS(system.file("extdata/examples/tiny_input.rds", package = "gradepath"))
#'
#' \donttest{
#' fit <- krw_report_card(inp, demographic = "race",
#'                        control = gp_control(backend = "highs",
#'                                             precision_rule = "krw_gmm"))
#' fit                       # result-first console print
#' get_report_card(fit)      # the per-unit grade-label report card
#' }
#'
#' @seealso [gp_control()], [gp_report_card()], [gp_grade_path()],
#'   [gp_select_grade()], [get_report_card()], [gp_plot_report_card()],
#'   [gp_twolevel_report_card()], [krw_firms]
#' @family gradepath-pipeline
#' @export
krw_report_card <- function(data,
                            demographic = c("race", "gender"),
                            groups = NULL,
                            control = gp_control(precision_rule = "krw_gmm"),
                            lambda = 0.25,
                            acceptance_mode = FALSE,
                            ...) {
  .gp_api_check_empty_dots("krw_report_card", ...)
  demographic <- .gp_api_match_demographic(demographic)
  acceptance_mode <- .gradepath_validate_scalar_logical(
    acceptance_mode,
    "acceptance_mode"
  )
  if (!is.null(groups)) {
    .gp_status_abort(
      "GROUPS_ERROR",
      "`groups` is not implemented in the M1 monolith; use `groups = NULL` for the one-level path."
    )
  }
  control <- .gp_api_monolith_control(control)
  lambda <- .gp_api_report_lambda(lambda, control)
  input <- .gp_api_estimates(data, demographic)

  coupling <- .gp_api_with_status(
    "GMM_NONCONVERGED",
    gp_w_seam(input$estimates, demographic, control = control)
  )
  precision_fit <- .gp_precision_fit_from_coupling(coupling)
  r_estimates <- .gp_r_estimates_from_coupling(coupling, input)

  prior <- .gp_api_with_status(
    "DECONV_BOUNDARY_ERROR",
    gp_deconvolve(coupling, control = control)
  )
  weights <- .gp_api_with_status(
    "WEIGHT_DEGENERATE",
    gp_posterior_weights(prior, r_estimates, precision_fit)
  )
  posterior <- .gp_api_with_status(
    "WEIGHT_DEGENERATE",
    gp_posterior_onelevel(
      weights,
      prior = prior,
      r_estimates = r_estimates,
      precision_fit = precision_fit,
      control = control
    )
  )
  pairwise <- .gp_api_with_status(
    "PAIRWISE_INVARIANT_ERROR",
    gp_pairwise(weights, ids = input$ids, control = control)
  )
  grade_path <- .gp_api_with_status(
    "SOLVER_BACKEND_UNAVAILABLE",
    gp_grade_path(
      pairwise,
      control = control,
      selected_lambda = lambda,
      acceptance_mode = acceptance_mode
    )
  )
  selected_grade <- gp_select_grade(grade_path, lambda = lambda)
  producer_status <- .gp_producer_status_from_selected_grade(selected_grade)
  report_card <- gp_report_card(
    input$estimates,
    posterior = posterior,
    selected_grade = selected_grade,
    grade_path = grade_path
  )

  validate_gp_fit(new_gp_fit(
    ids = input$ids,
    estimates = input$estimates,
    prior = prior,
    posterior = posterior,
    precision_fit = precision_fit,
    pairwise = pairwise,
    grade_path = grade_path,
    selected_grade = selected_grade,
    report_card = report_card,
    control = control,
    provenance = .gradepath_new_provenance(
      producer = "krw_report_card",
      producer_status = producer_status,
      demographic = demographic,
      workflow = "one_level_independence",
      selected_lambda = lambda
    ),
    warnings = unique(c(
      grade_path$warnings,
      selected_grade$warnings,
      report_card$warnings
    ))
  ))
}

#' @rdname krw_report_card
#' @export
gradepath <- krw_report_card
