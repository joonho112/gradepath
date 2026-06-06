# =============================================================================
# status.R -- public-runtime status contract
# -----------------------------------------------------------------------------
# Status values are runtime producer states, not tolerance classes. They are kept
# internal in M1 so the Phase-5 harness can consume one shared vocabulary without
# turning every stage into a new exported API.
# =============================================================================

.gp_status_values <- c(
  "OK",
  "INPUT_ERROR",
  "GROUPS_ERROR",
  "GMM_NONCONVERGED",
  "DECONV_BOUNDARY_ERROR",
  "WEIGHT_DEGENERATE",
  "PAIRWISE_INVARIANT_ERROR",
  "SOLVER_BACKEND_UNAVAILABLE",
  "SOLVER_INFEASIBLE",
  "SOLVER_TIME_LIMIT",
  "SOLVER_GAP",
  "SOLVER_OUTPUT_INVALID",
  "SOLVER_OBJECTIVE_MISMATCH",
  "SOLVER_CANONICAL_MISMATCH",
  "CACHE_STALE",
  "APPROXIMATE_OK",
  "UNVERIFIED",
  "INTERNAL_ERROR"
)

.gp_status_match <- function(status, name = "status") {
  status <- .gradepath_validate_scalar_character(status, name)
  status <- toupper(status)
  if (!status %in% .gp_status_values) {
    .gradepath_abort(
      "`%s` must be one of: %s.",
      name,
      paste(.gp_status_values, collapse = ", ")
    )
  }
  status
}

new_gp_status <- function(status = "OK", message = NULL, detail = list()) {
  status <- .gp_status_match(status)
  if (!is.null(message)) {
    message <- .gradepath_validate_scalar_character(message, "message")
  }
  .gradepath_validate_named_list(detail, "detail")
  structure(
    list(status = status, message = message, detail = detail),
    class = c("gp_status", "list")
  )
}

validate_gp_status <- function(x) {
  .gradepath_validate_list_class(x, "gp_status")
  .gradepath_validate_named_fields(x, c("status", "message", "detail"), "gp_status")
  .gp_status_match(x$status)
  if (!is.null(x$message)) {
    .gradepath_validate_scalar_character(x$message, "gp_status$message")
  }
  .gradepath_validate_named_list(x$detail, "gp_status$detail")
  x
}

.gp_status_ok <- function(message = "OK", detail = list()) {
  new_gp_status("OK", message = message, detail = detail)
}

.gp_status_abort <- function(status, message, ..., class = NULL, call = FALSE) {
  status <- .gp_status_match(status)
  if (!is.character(message) || length(message) != 1L ||
      is.na(message) || !nzchar(message)) {
    stop("`message` must be a length-1 non-empty character string.", call. = FALSE)
  }
  if (!is.null(class) &&
      (!is.character(class) || anyNA(class) || any(!nzchar(class)))) {
    stop("`class` must be NULL or a non-empty character vector.", call. = FALSE)
  }

  condition <- structure(
    list(
      message = sprintf(message, ...),
      call = if (isTRUE(call)) sys.call(-1L) else NULL,
      status = status
    ),
    class = unique(c(
      paste0("gp_status_", tolower(status)),
      class,
      "gradepath_error",
      "error",
      "condition"
    ))
  )
  stop(condition)
}

.gp_status_normalize <- function(status) {
  if (is.null(status) || length(status) != 1L || is.na(status)) {
    return("UNVERIFIED")
  }
  key <- tolower(as.character(status))
  mapped <- switch(
    key,
    ok = "OK",
    optimal = "OK",
    input_error = "INPUT_ERROR",
    groups_error = "GROUPS_ERROR",
    infeasible = "SOLVER_INFEASIBLE",
    solver_infeasible = "SOLVER_INFEASIBLE",
    time_limit = "SOLVER_TIME_LIMIT",
    solver_time_limit = "SOLVER_TIME_LIMIT",
    gap = "SOLVER_GAP",
    gap_reached = "SOLVER_GAP",
    suboptimal = "SOLVER_GAP",
    solver_gap = "SOLVER_GAP",
    output_invalid = "SOLVER_OUTPUT_INVALID",
    solver_output_invalid = "SOLVER_OUTPUT_INVALID",
    objective_mismatch = "SOLVER_OBJECTIVE_MISMATCH",
    solver_objective_mismatch = "SOLVER_OBJECTIVE_MISMATCH",
    canonical_mismatch = "SOLVER_CANONICAL_MISMATCH",
    solver_canonical_mismatch = "SOLVER_CANONICAL_MISMATCH",
    unavailable = "SOLVER_BACKEND_UNAVAILABLE",
    solver_backend_unavailable = "SOLVER_BACKEND_UNAVAILABLE",
    unverified = "UNVERIFIED",
    internal_error = "INTERNAL_ERROR",
    toupper(as.character(status))
  )
  if (mapped %in% .gp_status_values) mapped else "UNVERIFIED"
}

# Producer-status contract for the M1 acceptance harness.
#
# Runtime solver statuses are lower-level facts. Producer statuses are the
# package-level states consumed by gp_check()/gp_validate_targets(). Only `OK`
# is acceptance-ready: every other producer status must route registered target
# rows to UNVERIFIED, even if the numeric value happens to match a paper value.
# This is the M1 acceptance posture: SOLVER_GAP, SOLVER_TIME_LIMIT, backend
# unavailability, and internal errors are evidence states, not claim-ready
# passes.
#
# For a report-card fit, the selected lambda controls the fit-level producer
# status. Non-selected lambda failures remain path-level facts until a future
# path-wide target explicitly depends on them.
.gp_status_from_solver_status <- function(status) {
  .gp_status_normalize(status)
}

.gp_status_acceptance_ready <- function(status) {
  identical(.gp_status_normalize(status), "OK")
}

.gp_status_requires_unverified <- function(status) {
  !.gp_status_acceptance_ready(status)
}

.gp_selected_grade_solver_status <- function(selected_grade) {
  selected_grade <- validate_gp_grade_fit(selected_grade)
  status <- NULL
  if (!is.null(selected_grade$backend$status)) {
    status <- selected_grade$backend$status
  } else if (!is.null(selected_grade$summary$status)) {
    status <- selected_grade$summary$status
  }
  .gp_status_from_solver_status(status)
}

.gp_producer_status_from_selected_grade <- function(selected_grade) {
  .gp_selected_grade_solver_status(selected_grade)
}

.gp_status_contract <- function() {
  source_status <- c(
    "optimal",
    "OK",
    "APPROXIMATE_OK",
    "gap_reached",
    "suboptimal",
    "time_limit",
    "infeasible",
    "output_invalid",
    "objective_mismatch",
    "canonical_mismatch",
    "unavailable",
    "INTERNAL_ERROR",
    NA_character_
  )
  producer_status <- unname(vapply(
    source_status,
    .gp_status_from_solver_status,
    character(1)
  ))
  acceptance_ready <- unname(vapply(
    producer_status,
    .gp_status_acceptance_ready,
    logical(1)
  ))
  data.frame(
    source_status = source_status,
    producer_status = producer_status,
    acceptance_ready = acceptance_ready,
    target_routing = ifelse(
      acceptance_ready,
      "eligible_for_PASS_FAIL",
      "UNVERIFIED"
    ),
    stringsAsFactors = FALSE
  )
}

.gp_status_from_condition <- function(cnd) {
  if (!is.null(cnd$status)) {
    return(.gp_status_match(cnd$status))
  }
  cls <- class(cnd)
  if ("gradepath_gurobi_unavailable" %in% cls ||
      "gradepath_backend_unavailable" %in% cls) {
    return("SOLVER_BACKEND_UNAVAILABLE")
  }
  if ("gradepath_solver_infeasible" %in% cls) {
    return("SOLVER_INFEASIBLE")
  }
  if ("gradepath_solver_output_invalid" %in% cls) {
    return("SOLVER_OUTPUT_INVALID")
  }
  if ("gradepath_solver_objective_mismatch" %in% cls) {
    return("SOLVER_OBJECTIVE_MISMATCH")
  }
  if ("gradepath_solver_canonical_mismatch" %in% cls) {
    return("SOLVER_CANONICAL_MISMATCH")
  }
  if ("gp_status_input_error" %in% cls) {
    return("INPUT_ERROR")
  }
  if ("gp_status_groups_error" %in% cls) {
    return("GROUPS_ERROR")
  }
  if ("gradepath_validation_error" %in% cls) {
    return("INPUT_ERROR")
  }
  if ("gradepath_backend_error" %in% cls) {
    return("INTERNAL_ERROR")
  }
  if ("gradepath_error" %in% cls) {
    return("INTERNAL_ERROR")
  }
  "UNVERIFIED"
}
