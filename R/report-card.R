# =============================================================================
# report-card.R -- report-card assembly
# -----------------------------------------------------------------------------
# The report card is a unit-level payload. Rows are ordered by the endpoint
# (Condorcet, lambda = 1) fit; the displayed grade comes from the selected fit.
# Frontier quality metrics stay in gp_frontier / gp_grade_path.
# =============================================================================

.gp_report_card_check_empty_dots <- function(...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    if (is.null(dot_names)) {
      dot_names <- rep("<unnamed>", length(dots))
    } else {
      dot_names[is.na(dot_names) | !nzchar(dot_names)] <- "<unnamed>"
    }
    .gradepath_abort(
      "`gp_report_card()` does not accept unused arguments: %s.",
      paste(dot_names, collapse = ", ")
    )
  }
  invisible(NULL)
}

.gp_report_card_resolve_inputs <- function(estimates,
                                           posterior,
                                           selected_grade,
                                           grade_path) {
  if (inherits(estimates, "gp_fit")) {
    if (!is.null(posterior) || !is.null(selected_grade) || !is.null(grade_path)) {
      .gradepath_abort(
        "`gp_report_card()` with a `gp_fit` input does not accept stage-wise arguments."
      )
    }
    fit <- validate_gp_fit(estimates)
    return(list(fit = fit))
  }

  if (inherits(estimates, "gp_grade_path")) {
    if (!is.null(grade_path)) {
      .gradepath_abort("`grade_path` was supplied twice.")
    }
    grade_path <- estimates
    estimates <- NULL
    if (is.null(selected_grade)) {
      selected_grade <- gp_select_grade(
        grade_path,
        lambda = grade_path$selection$selected_lambda
      )
    }
  } else if (inherits(estimates, "gp_grade_fit")) {
    if (!is.null(selected_grade)) {
      .gradepath_abort("`selected_grade` was supplied twice.")
    }
    selected_grade <- estimates
    estimates <- NULL
  }

  .missing_args <- c(
    if (is.null(posterior)) "posterior",
    if (is.null(selected_grade)) "selected_grade",
    if (is.null(grade_path)) "grade_path"
  )
  if (length(.missing_args) > 0L) {
    .gradepath_abort(paste0(
      "`gp_report_card()` stage-wise assembly is missing required argument(s): ",
      paste(.missing_args, collapse = ", "),
      ". Supply `posterior`, `selected_grade`, and `grade_path` together -- e.g. ",
      "from a fit's slots: gp_report_card(estimates = fit$estimates, ",
      "posterior = get_posterior(fit), selected_grade = fit$selected_grade, ",
      "grade_path = fit$grade_path)."
    ))
  }

  list(
    estimates = estimates,
    posterior = posterior,
    selected_grade = selected_grade,
    grade_path = grade_path
  )
}

.gp_report_card_vector <- function(x, name, n, numeric = FALSE) {
  if (is.null(x)) {
    .gradepath_abort("`%s` must not be NULL.", name)
  }
  if (isTRUE(numeric)) {
    out <- .gradepath_validate_numeric_vector(x, name)
  } else {
    out <- as.character(x)
    if (length(out) != length(x) || anyNA(out)) {
      .gradepath_abort("`%s` must be a character vector.", name)
    }
  }
  if (length(out) != n) {
    .gradepath_abort("`%s` must have length %d.", name, as.integer(n))
  }
  out
}

.gp_report_card_estimates_id <- function(estimates, n) {
  value <- estimates$unit_id %gp_or% estimates$id %gp_or% estimates$ids
  if (is.null(value)) {
    return(NULL)
  }
  .gp_report_card_vector(value, "estimates id", n)
}

.gp_report_card_estimates_label <- function(estimates, ids, fallback) {
  n <- length(ids)
  label <- estimates$label %gp_or% estimates$labels
  if (is.null(label) && !is.null(estimates$covariates) &&
      is.data.frame(estimates$covariates) &&
      "label" %in% names(estimates$covariates)) {
    label <- estimates$covariates$label
  }
  if (is.null(label)) {
    return(fallback)
  }
  if (!is.null(names(label))) {
    label <- label[ids]
  }
  .gp_report_card_vector(label, "estimates label", n)
}

.gp_report_card_estimates_frame <- function(estimates, posterior) {
  n <- length(posterior$id)
  ids <- as.character(posterior$id)

  if (is.null(estimates)) {
    return(data.frame(
      id = ids,
      label = as.character(posterior$label),
      estimate = as.numeric(posterior$estimate),
      se = as.numeric(posterior$se),
      stringsAsFactors = FALSE
    ))
  }

  if (inherits(estimates, "eb_estimates") && requireNamespace("ebrecipe", quietly = TRUE)) {
    .gp_eb_validate_estimates(estimates)
  }

  estimate_id <- .gp_report_card_estimates_id(estimates, n)
  if (!is.null(estimate_id) && !identical(estimate_id, ids)) {
    .gradepath_abort("`estimates` ids must align exactly with `posterior$id`.")
  }

  estimate <- estimates$theta_hat %gp_or% estimates$estimate
  se <- estimates$s %gp_or% estimates$se
  label <- .gp_report_card_estimates_label(
    estimates,
    ids = ids,
    fallback = as.character(posterior$label)
  )

  data.frame(
    id = ids,
    label = label,
    estimate = .gp_report_card_vector(estimate, "estimates estimate", n, numeric = TRUE),
    se = .gp_report_card_vector(se, "estimates se", n, numeric = TRUE),
    stringsAsFactors = FALSE
  )
}

.gp_report_card_posterior_frame <- function(posterior) {
  posterior <- validate_gp_posterior(posterior)
  ids <- .gradepath_validate_character_vector(
    as.character(posterior$id),
    "gp_posterior$id",
    unique = TRUE
  )
  n <- length(ids)

  values <- .gp_posterior_report_card_values(posterior)
  if (!identical(values$posterior_mean, as.numeric(posterior$posterior_mean)) ||
      !identical(values$lower, as.numeric(posterior$lower)) ||
      !identical(values$upper, as.numeric(posterior$upper))) {
    mean <- .gp_report_card_vector(
      values$posterior_mean,
      "gp_posterior$metadata$reporting$posterior_mean",
      n,
      numeric = TRUE
    )
    lower <- .gp_report_card_vector(
      values$lower,
      "gp_posterior$metadata$reporting$lower",
      n,
      numeric = TRUE
    )
    upper <- .gp_report_card_vector(
      values$upper,
      "gp_posterior$metadata$reporting$upper",
      n,
      numeric = TRUE
    )
  } else {
    mean <- posterior$posterior_mean
    lower <- posterior$lower
    upper <- posterior$upper
  }

  data.frame(
    id = ids,
    label = as.character(posterior$label),
    posterior_mean = as.numeric(mean),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    stringsAsFactors = FALSE
  )
}

.gp_report_card_fit_index <- function(grade_path, lambda, label) {
  match <- .gp_grade_exact_lambda_match(
    grade_path$lambda_grid,
    lambda,
    label
  )
  match$index
}

.gp_report_card_endpoint_fit <- function(grade_path) {
  grade_path$fits[[.gp_report_card_fit_index(
    grade_path,
    lambda = grade_path$selection$endpoint_lambda,
    label = "endpoint_lambda"
  )]]
}

.gp_report_card_selected_fit <- function(grade_path) {
  grade_path$fits[[.gp_report_card_fit_index(
    grade_path,
    lambda = grade_path$selection$selected_lambda,
    label = "selected_lambda"
  )]]
}

.gp_report_card_sort_index <- function(endpoint_fit) {
  endpoint_fit <- validate_gp_grade_fit(endpoint_fit)
  order(endpoint_fit$assignment$grade, endpoint_fit$assignment$id, method = "radix")
}

.gp_report_card_validate_stagewise <- function(estimates,
                                               posterior,
                                               selected_grade,
                                               grade_path) {
  posterior <- validate_gp_posterior(posterior)
  selected_grade <- validate_gp_grade_fit(selected_grade)
  grade_path <- validate_gp_grade_path(grade_path)

  ids <- as.character(posterior$id)
  if (!identical(ids, grade_path$ids) || !identical(ids, selected_grade$ids)) {
    .gradepath_abort(
      "`posterior$id`, `selected_grade$ids`, and `grade_path$ids` must align exactly."
    )
  }
  if (!identical(selected_grade$control, grade_path$control)) {
    .gradepath_abort(
      "`selected_grade` and `grade_path` must share identical controls."
    )
  }

  selected_fit <- .gp_report_card_selected_fit(grade_path)
  if (!identical(selected_grade, selected_fit)) {
    .gradepath_abort(
      "`selected_grade` must be the stored path member at `grade_path$selection$selected_lambda`."
    )
  }

  list(
    estimates = .gp_report_card_estimates_frame(estimates, posterior),
    posterior = .gp_report_card_posterior_frame(posterior),
    selected_grade = selected_grade,
    grade_path = grade_path,
    endpoint_fit = .gp_report_card_endpoint_fit(grade_path)
  )
}

.gp_report_card_data <- function(inputs) {
  sort_index <- .gp_report_card_sort_index(inputs$endpoint_fit)
  ordered_ids <- inputs$endpoint_fit$assignment$id[sort_index]
  selected_index <- match(ordered_ids, inputs$selected_grade$assignment$id)
  posterior_index <- match(ordered_ids, inputs$posterior$id)
  estimate_index <- match(ordered_ids, inputs$estimates$id)

  if (anyNA(selected_index) || anyNA(posterior_index) || anyNA(estimate_index)) {
    .gradepath_abort_internal(
      "Validated report-card inputs could not be aligned by id."
    )
  }

  data.frame(
    id = ordered_ids,
    label = inputs$estimates$label[estimate_index],
    grade = as.integer(inputs$selected_grade$assignment$grade[selected_index]),
    sort_rank = seq_along(ordered_ids),
    selected_lambda = rep(inputs$selected_grade$lambda, length(ordered_ids)),
    posterior_mean = inputs$posterior$posterior_mean[posterior_index],
    lower = inputs$posterior$lower[posterior_index],
    upper = inputs$posterior$upper[posterior_index],
    estimate = inputs$estimates$estimate[estimate_index],
    se = inputs$estimates$se[estimate_index],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Assemble a unit-level report card
#'
#' `gp_report_card()` assembles the per-unit report-card table: each unit with
#' its grade label, posterior summary, and original estimate, sorted from
#' most-extreme to least-extreme. It does not solve a grade problem -- the row
#' order comes from the endpoint fit cached in `grade_path`, while the displayed
#' `grade` column comes from `selected_grade`. Pass a finished `gp_fit` to read
#' back its stored card, or the stage objects to assemble one yourself; the
#' result prints, plots with [gp_plot_report_card()], and is what [get_report_card()]
#' returns.
#'
#' @param estimates One of: an `ebrecipe::eb_estimates`-like input container, a
#'   `gp_fit` (its stored report card is returned as-is), a `gp_grade_path`, or a
#'   `gp_grade_fit`. For stage-wise assembly it supplies the report-card `label`,
#'   `estimate`, and `se`; if omitted (`NULL`, the default) those columns are read
#'   from `posterior`.
#' @param posterior A `gp_posterior` (from [get_posterior()]) giving each unit's
#'   `posterior_mean`, `lower`, and `upper`. Required for stage-wise assembly;
#'   ignored when `estimates` is a `gp_fit`.
#' @param selected_grade A `gp_grade_fit`, usually from [gp_select_grade()];
#'   supplies the displayed `grade` column and the selected penalty. Required for
#'   stage-wise assembly.
#' @param grade_path A solved `gp_grade_path`; its endpoint fit fixes the row
#'   sort order (most-extreme grade first, ties broken by id). Required for
#'   stage-wise assembly.
#' @param ... Unused; reserved for future arguments. Passing any value raises an
#'   error.
#'
#' @return A validated `gp_report_card` object (a list of class
#'   `c("gp_report_card", "list")`) with the public slots: \describe{
#'   \item{`ids`}{Character vector; the unit-id order of the table rows
#'     (most-extreme grade first).}
#'   \item{`table`}{Data frame, one row per unit, with `id`, `label`, the integer
#'     `grade`, `sort_rank`, `selected_lambda`, the posterior summary
#'     `posterior_mean`/`lower`/`upper`, and the original `estimate`/`se`.}
#'   \item{`selected_lambda`}{Numeric; the penalty the displayed grades were taken
#'     from.}
#'   \item{`grades`}{Integer vector of the per-unit grade labels, aligned to `ids`.}
#'   \item{`control`}{The [gp_control] carried from `grade_path`.}
#'   \item{`provenance`, `warnings`, `schema_version`}{Producer metadata and
#'     internal audit slots.}
#' }
#'
#' @details Grade labels are integers in `{1, ..., n}` and carry no
#'   ranking-superiority statement of any kind: grade 1 is simply the most-extreme
#'   theta block.
#'
#' @note Welfare orientation (frozen INVARIANT 2). Grade 1 is the most-extreme
#'   theta, and its welfare reading FLIPS with the application: for firms the
#'   most-extreme block is the most discriminatory, while for names it is the
#'   best-treated. The label is therefore not hard-coded to one application; the
#'   console print and [format_gp_report_card_cli()] surface this same note. See
#'   Appendix A.
#'
#' @examples
#' # Read the stored card from the bundled fit (no Gurobi, no solve).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' gp_report_card(fit)              # returns fit$report_card
#'
#' # Re-assemble the same card from the fit's stage slots (solve-free).
#' rc <- gp_report_card(
#'   estimates      = fit$estimates,
#'   posterior      = fit$posterior,
#'   selected_grade = fit$selected_grade,
#'   grade_path     = fit$grade_path
#' )
#' head(rc$table)                   # per-unit grade-label table
#'
#' @seealso [gp_select_grade()], [gp_grade_path()], [get_report_card()],
#'   [gp_plot_report_card()], [format_gp_report_card_cli()]
#' @family gradepath-report-card
#' @export
gp_report_card <- function(estimates = NULL,
                           posterior = NULL,
                           selected_grade = NULL,
                           grade_path = NULL,
                           ...) {
  .gp_report_card_check_empty_dots(...)
  inputs <- .gp_report_card_resolve_inputs(
    estimates = estimates,
    posterior = posterior,
    selected_grade = selected_grade,
    grade_path = grade_path
  )
  if (!is.null(inputs$fit)) {
    return(validate_gp_report_card(inputs$fit$report_card))
  }

  inputs <- .gp_report_card_validate_stagewise(
    estimates = inputs$estimates,
    posterior = inputs$posterior,
    selected_grade = inputs$selected_grade,
    grade_path = inputs$grade_path
  )
  table <- .gp_report_card_data(inputs)

  validate_gp_report_card(new_gp_report_card(
    ids = table$id,
    table = table,
    selected_lambda = inputs$selected_grade$lambda,
    grades = table$grade,
    control = inputs$grade_path$control,
    provenance = .gradepath_new_provenance(
      producer = "gp_report_card",
      selected_lambda = inputs$selected_grade$lambda,
      endpoint_lambda = inputs$grade_path$selection$endpoint_lambda,
      sort_key = "endpoint_grade_then_id"
    ),
    warnings = unique(c(
      inputs$selected_grade$warnings,
      inputs$grade_path$warnings
    ))
  ))
}
