# =============================================================================
# api-input.R -- shared public API input adapters
# -----------------------------------------------------------------------------
# The monolith and preview both enter through user data. This file keeps the
# column grammar and ebrecipe::eb_input boundary in one place.
# =============================================================================

.gp_api_match_demographic <- function(demographic) {
  match.arg(demographic, c("race", "gender"))
}

.gp_api_check_empty_dots <- function(fn, ...) {
  dots <- list(...)
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }
  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("<unnamed>", length(dots))
  } else {
    dot_names[is.na(dot_names) | !nzchar(dot_names)] <- "<unnamed>"
  }
  .gradepath_abort(
    "`%s()` does not accept unused argument(s): %s.",
    fn,
    paste(dot_names, collapse = ", ")
  )
}

.gp_api_first_field <- function(data, names) {
  for (name in names) {
    if (!is.null(data[[name]])) {
      return(data[[name]])
    }
  }
  NULL
}

.gp_api_available_fields <- function(data, names) {
  names[vapply(names, function(name) !is.null(data[[name]]), logical(1))]
}

.gp_api_demographic_field <- function(data,
                                      demographic,
                                      specific,
                                      conflict_specific = specific,
                                      generic,
                                      label,
                                      data_name) {
  specific_hit <- .gp_api_available_fields(data, specific)
  conflict_specific_hit <- .gp_api_available_fields(data, conflict_specific)
  generic_hit <- .gp_api_available_fields(data, generic)
  if (length(conflict_specific_hit) > 0L && length(generic_hit) > 0L) {
    .gp_status_abort(
      "INPUT_ERROR",
      "`%s` contains both demographic-specific column(s) %s and generic column(s) %s for `%s` with demographic `%s`; provide only one form.",
      data_name,
      paste(sprintf("`%s`", conflict_specific_hit), collapse = ", "),
      paste(sprintf("`%s`", generic_hit), collapse = ", "),
      label,
      demographic
    )
  }
  if (length(specific_hit) > 0L) {
    return(data[[specific_hit[[1L]]]])
  }
  .gp_api_first_field(data, generic)
}

.gp_api_input_vector <- function(x, name, n = NULL, numeric = FALSE,
                                 positive = FALSE) {
  if (is.null(x)) {
    .gp_status_abort("INPUT_ERROR", "`%s` is required.", name)
  }
  if (isTRUE(numeric)) {
    out <- suppressWarnings(as.numeric(x))
    if (length(out) != length(x) || anyNA(out) || any(!is.finite(out))) {
      .gp_status_abort("INPUT_ERROR", "`%s` must be finite numeric.", name)
    }
    if (isTRUE(positive) && any(out <= 0)) {
      .gp_status_abort("INPUT_ERROR", "`%s` must be strictly positive.", name)
    }
  } else {
    out <- as.character(x)
    if (length(out) != length(x) || anyNA(out) || any(!nzchar(out))) {
      .gp_status_abort("INPUT_ERROR", "`%s` must be non-empty character values.", name)
    }
  }
  if (!is.null(n) && length(out) != n) {
    .gp_status_abort("INPUT_ERROR", "`%s` must have length %d.", name, as.integer(n))
  }
  out
}

.gp_api_validate_control <- function(control, name = "control") {
  tryCatch(
    validate_gp_control(control),
    error = function(e) {
      .gp_status_abort(
        "INPUT_ERROR",
        "`%s` must be a valid gp_control object: %s",
        name,
        conditionMessage(e),
        class = setdiff(class(e), c("error", "condition"))
      )
    }
  )
}

.gp_api_labels_from_estimates <- function(estimates, ids) {
  labels <- estimates$label %gp_or% estimates$labels
  if (is.null(labels) && is.data.frame(estimates$covariates) &&
      "label" %in% names(estimates$covariates)) {
    labels <- estimates$covariates$label
  }
  if (is.null(labels)) {
    return(ids)
  }
  labels <- as.character(labels)
  if (!is.null(names(labels))) {
    labels <- labels[ids]
  }
  .gp_api_input_vector(labels, "label", n = length(ids))
}

.gp_api_input_columns <- function(data, demographic, data_name = "data") {
  if (!is.list(data) && !is.data.frame(data)) {
    .gp_status_abort(
      "INPUT_ERROR",
      "`%s` must be a data frame, list, or ebrecipe::eb_estimates object.",
      data_name
    )
  }
  theta <- .gp_api_demographic_field(
    data,
    demographic = demographic,
    specific = paste0("theta_hat_", demographic),
    conflict_specific = c("theta_hat_race", "theta_hat_gender"),
    generic = c("theta_hat", "estimate"),
    label = "theta_hat",
    data_name = data_name
  )
  se <- .gp_api_demographic_field(
    data,
    demographic = demographic,
    specific = paste0("se_", demographic),
    conflict_specific = c("se_race", "se_gender"),
    generic = c("s", "se"),
    label = "s",
    data_name = data_name
  )
  ids <- .gp_api_first_field(
    data,
    c("unit_id", "firm_id", "id", "ids")
  )

  theta <- .gp_api_input_vector(theta, "theta_hat", numeric = TRUE)
  n <- length(theta)
  se <- .gp_api_input_vector(se, "s", n = n, numeric = TRUE, positive = TRUE)
  if (is.null(ids)) {
    ids <- as.character(seq_len(n))
  } else {
    ids <- .gp_api_input_vector(ids, "unit_id", n = n)
  }
  if (anyDuplicated(ids)) {
    .gp_status_abort("INPUT_ERROR", "`unit_id` values must be unique.")
  }

  label <- .gp_api_first_field(data, c("label", "labels"))
  if (is.null(label)) {
    label <- ids
  } else {
    label <- .gp_api_input_vector(label, "label", n = n)
  }

  covariates <- data.frame(label = label, stringsAsFactors = FALSE)
  industry <- .gp_api_first_field(data, c("industry", "group", "group_id"))
  if (!is.null(industry)) {
    covariates$industry <- .gp_api_input_vector(industry, "industry", n = n)
  }

  list(
    theta_hat = theta,
    s = se,
    ids = ids,
    label = label,
    covariates = covariates
  )
}

.gp_api_estimates <- function(data, demographic, data_name = "data") {
  demographic <- .gp_api_match_demographic(demographic)

  if (inherits(data, "eb_estimates")) {
    tryCatch(
      .gp_eb_validate_estimates(data),
      error = function(e) {
        .gp_status_abort(
          "INPUT_ERROR",
          "`%s` is not a valid ebrecipe::eb_estimates object: %s",
          data_name,
          conditionMessage(e)
        )
      }
    )
    theta <- .gp_api_input_vector(.gp_estimates_theta(data), "theta_hat", numeric = TRUE)
    n <- length(theta)
    se <- .gp_api_input_vector(.gp_estimates_se(data), "s", n = n,
                               numeric = TRUE, positive = TRUE)
    ids <- .gp_estimates_id(data)
    if (is.null(ids)) {
      ids <- as.character(seq_len(n))
    } else {
      ids <- .gp_api_input_vector(ids, "unit_id", n = n)
    }
    if (anyDuplicated(ids)) {
      .gp_status_abort("INPUT_ERROR", "`unit_id` values must be unique.")
    }
    original_s <- .gp_estimates_original_s(data)
    if (is.null(original_s)) {
      data$original_s <- se
      original_s <- se
    } else {
      original_s <- .gp_api_input_vector(original_s, "original_s", n = n,
                                         numeric = TRUE, positive = TRUE)
      data$original_s <- original_s
    }
    label <- .gp_api_labels_from_estimates(data, ids)
    return(list(
      estimates = data,
      ids = ids,
      label = label,
      theta_hat = theta,
      s = se,
      original_s = original_s,
      demographic = demographic
    ))
  }

  cols <- .gp_api_input_columns(data, demographic, data_name = data_name)
  estimates <- tryCatch(
    .gp_eb_input(
      theta_hat = cols$theta_hat,
      s = cols$s,
      unit_id = cols$ids,
      covariates = cols$covariates
    ),
    error = function(e) {
      .gp_status_abort(
        "INPUT_ERROR",
        "Unable to create ebrecipe input from `%s`: %s",
        data_name,
        conditionMessage(e)
      )
    }
  )
  estimates$original_s <- cols$s

  list(
    estimates = estimates,
    ids = cols$ids,
    label = cols$label,
    theta_hat = cols$theta_hat,
    s = cols$s,
    original_s = cols$s,
    demographic = demographic
  )
}

.gp_api_monolith_control <- function(control) {
  if (is.null(control)) {
    control <- gp_control(precision_rule = "krw_gmm")
  }
  control <- .gp_api_validate_control(control)
  if (!identical(control$precision_rule, "krw_gmm")) {
    .gp_status_abort(
      "INPUT_ERROR",
      "`control$precision_rule` must be \"krw_gmm\" for `krw_report_card()`."
    )
  }
  control
}

.gp_api_operational_lambda_grid <- function(control) {
  control <- .gp_api_validate_control(control)
  if (identical(as.numeric(control$lambda_grid),
                as.numeric(.gp_control_default_lambda_grid))) {
    .gp_grade_operational_default_grid()
  } else {
    control$lambda_grid
  }
}

.gp_api_report_lambda <- function(lambda, control) {
  lambda <- tryCatch(
    .gp_grade_validate_lambda(lambda, "lambda"),
    error = function(e) {
      .gp_status_abort(
        "INPUT_ERROR",
        "%s",
        conditionMessage(e),
        class = setdiff(class(e), c("error", "condition"))
      )
    }
  )

  grid <- .gp_api_operational_lambda_grid(control)
  if (sum(abs(grid - lambda) < 1e-8) != 1L) {
    .gp_status_abort(
      "INPUT_ERROR",
      "`lambda` must match exactly one solved lambda in `control$lambda_grid`. Available values: %s.",
      paste(format(grid), collapse = ", ")
    )
  }
  lambda
}
