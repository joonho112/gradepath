# =============================================================================
# preview.R -- dry-run preview
# -----------------------------------------------------------------------------
# Validate user input and report the planned one-level workflow and grade-IP
# solve count without running GMM, deconvolution, posterior weights, pairwise,
# or any solver backend.
# =============================================================================

.gp_preview_fields <- c(
  "ids",
  "n_units",
  "demographic",
  "groups",
  "workflow",
  "lambda_grid",
  "estimated_solves",
  "backend",
  "control",
  "status",
  "schema_version",
  "provenance",
  "warnings"
)

new_gp_preview <- function(ids, n_units, demographic, groups, workflow,
                           lambda_grid, estimated_solves, backend, control,
                           status = .gp_status_ok(), warnings = character()) {
  structure(
    list(
      ids = ids,
      n_units = n_units,
      demographic = demographic,
      groups = groups,
      workflow = workflow,
      lambda_grid = lambda_grid,
      estimated_solves = estimated_solves,
      backend = backend,
      control = control,
      status = status,
      schema_version = .gradepath_schema_version,
      provenance = .gradepath_new_provenance(
        producer = "gp_preview",
        no_solve = TRUE
      ),
      warnings = warnings
    ),
    class = c("gp_preview", "list")
  )
}

validate_gp_preview <- function(x) {
  .gradepath_validate_list_class(x, "gp_preview")
  .gradepath_validate_named_fields(x, .gp_preview_fields, "gp_preview")

  ids <- .gradepath_validate_character_vector(x$ids, "gp_preview$ids", unique = TRUE)
  n_units <- .gradepath_validate_integerish(x$n_units, "gp_preview$n_units", min = 1L)
  if (!identical(n_units, length(ids))) {
    .gradepath_abort("`gp_preview$n_units` must equal length(ids).")
  }
  .gradepath_validate_scalar_character(
    x$demographic,
    "gp_preview$demographic",
    allowed = c("race", "gender")
  )
  if (!is.null(x$groups)) {
    .gradepath_validate_character_vector(x$groups, "gp_preview$groups")
    if (length(x$groups) != n_units) {
      .gradepath_abort("`gp_preview$groups` must have length n_units.")
    }
  }
  .gradepath_validate_scalar_character(
    x$workflow,
    "gp_preview$workflow",
    allowed = c("one_level_independence", "grouped_pending")
  )
  .gp_validate_lambda_grid_core(x$lambda_grid, "gp_preview$lambda_grid")
  .gradepath_validate_integerish(
    x$estimated_solves,
    "gp_preview$estimated_solves",
    min = 0L
  )
  .gradepath_validate_scalar_character(
    x$backend,
    "gp_preview$backend",
    allowed = .gp_control_backends_canonical
  )
  validate_gp_control(x$control)
  validate_gp_status(x$status)
  if (!is.character(x$warnings) || anyNA(x$warnings)) {
    .gradepath_abort("`gp_preview$warnings` must be a character vector.")
  }
  .gradepath_validate_scalar_character(x$schema_version, "gp_preview$schema_version")
  .gradepath_validate_named_list(x$provenance, "gp_preview$provenance")

  x
}

#' Preview a gradepath run without solving
#'
#' `gp_preview()` is the solve-free dry run for the one-level KRW workflow. Give
#' it the same firm-level input and demographic you would hand to
#' [krw_report_card()]; it validates the input grammar and reports the inferred
#' plan -- the unit ids, the inferred workflow, the chosen backend, the
#' operational frontier penalty grid, and the number of grade-IP solves the run
#' would perform -- then returns a `gp_preview` object WITHOUT running the
#' beta-GMM, deconvolution, posterior, pairwise, or any solver stage. Reach for
#' it to size and sanity-check a run before committing to a full solve.
#'
#' @param data Data frame, list, or `ebrecipe::eb_estimates` object carrying
#'   `theta_hat`/`s` (or the demographic-specific `theta_hat_race`/`se_race`,
#'   `theta_hat_gender`/`se_gender`) plus an optional `unit_id`/`firm_id`/`label`.
#'   Do not supply both a generic and a demographic-specific column for the same
#'   quantity in one input.
#' @param demographic Character scalar, `"race"` or `"gender"`; selects the
#'   estimate columns the input grammar reads. Defaults to `"race"`.
#' @param groups Optional grouping vector, or `NULL`. Grouped/two-level execution
#'   is not implemented in the M1 monolith; a non-`NULL` value is previewed as a
#'   pending grouped workflow (with an explanatory warning and a `GROUPS_ERROR`
#'   status), not run. Defaults to `NULL` (the one-level workflow).
#' @param control A [gp_control] object, or `NULL` to use the defaults. Supplies
#'   the backend and the lambda grid the preview reports. Defaults to
#'   [gp_control()].
#' @param ... Unused; an error is raised if any argument is passed.
#'
#' @return A validated S3 object of class `c("gp_preview", "list")`: the planned
#'   run with no results. \describe{
#'   \item{`ids`}{Character vector; the canonical unit-id order parsed from the
#'     input.}
#'   \item{`n_units`}{Integer scalar; the number of units (equals
#'     `length(ids)`).}
#'   \item{`demographic`}{Character scalar; `"race"` or `"gender"`.}
#'   \item{`groups`}{Character vector of length `n_units`, or `NULL`; the parsed
#'     grouping when one was supplied.}
#'   \item{`workflow`}{Character scalar; the inferred workflow, either
#'     `"one_level_independence"` or `"grouped_pending"`.}
#'   \item{`lambda_grid`}{Numeric vector; the operational frontier penalty grid
#'     the run would sweep.}
#'   \item{`estimated_solves`}{Integer scalar; the number of grade-IP solves the
#'     run would perform (the length of `lambda_grid`).}
#'   \item{`backend`}{Character scalar; the grade-IP backend the run would use.}
#'   \item{`control`}{The validated [gp_control] the preview was built from.}
#'   \item{`status`}{A `gp_status`; `OK` for a one-level plan, or `GROUPS_ERROR`
#'     when a grouping was supplied.}
#'   \item{`schema_version`, `provenance`, `warnings`}{Internal audit slots; the
#'     `provenance` stamp records that no solve was performed.}
#' }
#'
#' @examples
#' # Solve-free: load the bundled tiny example input and preview the planned run
#' # (no Gurobi, no solve, instant).
#' inp <- readRDS(system.file("extdata/examples/tiny_input.rds", package = "gradepath"))
#' pv <- gp_preview(inp, demographic = "race")
#' pv                  # result-first console print of the planned run
#' pv$n_units          # number of units the run would grade
#' pv$estimated_solves # number of grade-IP solves the run would perform
#'
#' @seealso [krw_report_card()], [gp_grade_path()], [gp_control()]
#' @family gradepath-grade
#' @export
gp_preview <- function(data,
                       demographic = c("race", "gender"),
                       groups = NULL,
                       control = gp_control(),
                       ...) {
  .gp_api_check_empty_dots("gp_preview", ...)
  demographic <- .gp_api_match_demographic(demographic)
  control <- if (is.null(control)) gp_control() else .gp_api_validate_control(control)
  input <- .gp_api_estimates(data, demographic)

  warnings <- character()
  workflow <- "one_level_independence"
  status <- .gp_status_ok("Preview OK")
  group_values <- NULL
  if (!is.null(groups)) {
    group_values <- .gp_api_input_vector(
      groups,
      "groups",
      n = length(input$ids)
    )
    workflow <- "grouped_pending"
    status <- new_gp_status(
      "GROUPS_ERROR",
      message = "Grouped/two-level execution is not implemented in the M1 preview."
    )
    warnings <- c(
      warnings,
      "Grouped/two-level execution is planned for a later phase; the M1 monolith runs one-level only."
    )
  }

  grid <- .gp_api_operational_lambda_grid(control)
  validate_gp_preview(new_gp_preview(
    ids = input$ids,
    n_units = length(input$ids),
    demographic = demographic,
    groups = group_values,
    workflow = workflow,
    lambda_grid = grid,
    estimated_solves = length(grid),
    backend = control$backend,
    control = control,
    status = status,
    warnings = warnings
  ))
}
