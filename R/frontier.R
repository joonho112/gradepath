# =============================================================================
# frontier.R -- frontier metrics for solved grade paths
# -----------------------------------------------------------------------------
# Ports the v1 blockwise frontier metric recipe verbatim in spirit:
# DR / reliability / tau_bar are computed over between-grade blocks and
# normalized by 2 * choose(n, 2).
# =============================================================================

.gp_frontier_fields <- c(
  "ids",
  "table",
  "benchmarks",
  "dr_matrix",
  "selection",
  "control",
  "schema_version",
  "provenance",
  "warnings"
)

.gp_frontier_table_columns <- c(
  "lambda",
  "grade_count",
  "discordance_rate",
  "reliability",
  "tau_bar"
)

.gp_frontier_benchmark_columns <- c(
  "benchmark",
  "grade_count",
  "discordance_rate",
  "reliability",
  "tau_bar"
)

.gradepath_frontier_metrics <- function(grades, pairwise_matrix) {
  grades <- .gradepath_validate_numeric_vector(grades, "grades")
  pairwise_matrix <- .gradepath_validate_numeric_matrix(
    pairwise_matrix,
    "pairwise_matrix"
  )

  if (length(grades) != nrow(pairwise_matrix) ||
      length(grades) != ncol(pairwise_matrix)) {
    .gradepath_abort(
      "Frontier metrics require one grade label per row and column of `pairwise_matrix`."
    )
  }
  if (!identical(nrow(pairwise_matrix), ncol(pairwise_matrix))) {
    .gradepath_abort("`pairwise_matrix` must be square.")
  }
  if (any(!is.finite(pairwise_matrix)) ||
      any(pairwise_matrix < 0) ||
      any(pairwise_matrix > 1)) {
    .gradepath_abort("`pairwise_matrix` entries must be finite values in [0, 1].")
  }
  if (any(abs(grades - round(grades)) > 1e-8) || any(grades < 1)) {
    .gradepath_abort("`grades` must be positive integer labels.")
  }

  grades <- as.integer(round(grades))
  if (length(grades) < 2L || length(unique(grades)) < 2L) {
    return(list(
      discordance_rate = 0,
      reliability = 1,
      tau_bar = 0
    ))
  }

  ordered_grades <- sort(unique(grades))
  normfactor <- choose(length(grades), 2L) * 2L
  discordance_rate <- 0
  tau_bar <- 0
  pairwise_matrix_t <- t(pairwise_matrix)

  for (grade_i in ordered_grades) {
    index_i <- grades == grade_i

    for (grade_j in ordered_grades) {
      index_j <- grades == grade_j

      if (grade_i > grade_j) {
        discordance_rate <- discordance_rate +
          sum(pairwise_matrix[index_i, index_j, drop = FALSE])
        tau_bar <- tau_bar +
          sum(pairwise_matrix_t[index_i, index_j, drop = FALSE]) -
          sum(pairwise_matrix[index_i, index_j, drop = FALSE])
      }

      if (grade_i < grade_j) {
        discordance_rate <- discordance_rate +
          sum(pairwise_matrix_t[index_i, index_j, drop = FALSE])
        tau_bar <- tau_bar +
          sum(pairwise_matrix[index_i, index_j, drop = FALSE]) -
          sum(pairwise_matrix_t[index_i, index_j, drop = FALSE])
      }
    }
  }

  discordance_rate <- discordance_rate / normfactor
  tau_bar <- tau_bar / normfactor

  list(
    discordance_rate = discordance_rate,
    reliability = 1 - discordance_rate,
    tau_bar = tau_bar
  )
}

.gp_frontier_metric_scalar <- function(fit, name) {
  value <- fit$summary[[name]]
  if (is.null(value)) {
    return(NA_real_)
  }
  as.numeric(value)
}

.gp_frontier_fit_metric_list <- function(fit, pairwise_matrix) {
  fit <- validate_gp_grade_fit(fit)
  .gradepath_frontier_metrics(
    grades = fit$assignment$grade,
    pairwise_matrix = pairwise_matrix
  )
}

frontier_table <- function(pairwise, fits) {
  pairwise <- validate_gp_pairwise(pairwise)
  if (!is.list(fits) || length(fits) < 1L) {
    .gradepath_abort("`fits` must be a non-empty list of `gp_grade_fit` objects.")
  }
  fits <- lapply(fits, validate_gp_grade_fit)
  mismatched <- which(vapply(
    fits,
    function(fit) !identical(pairwise$ids, fit$ids),
    logical(1)
  ))
  if (length(mismatched) > 0L) {
    .gradepath_abort(
      "`pairwise$ids` must match the stored fit ids for every fit; mismatch at fit(s): %s.",
      paste(mismatched, collapse = ", ")
    )
  }

  metrics <- lapply(
    fits,
    .gp_frontier_fit_metric_list,
    pairwise_matrix = pairwise$matrix
  )

  data.frame(
    lambda = vapply(fits, function(fit) fit$lambda, numeric(1)),
    grade_count = vapply(fits, function(fit) fit$summary$grade_count, integer(1)),
    discordance_rate = vapply(metrics, function(x) x$discordance_rate, numeric(1)),
    reliability = vapply(metrics, function(x) x$reliability, numeric(1)),
    tau_bar = vapply(metrics, function(x) x$tau_bar, numeric(1)),
    objective = vapply(fits, function(fit) fit$objective$value, numeric(1)),
    raw_objective = vapply(fits, function(fit) fit$objective$raw, numeric(1)),
    status = vapply(fits, function(fit) fit$summary$status, character(1)),
    row.names = NULL,
    check.names = FALSE
  )
}

.gp_frontier_empty_benchmarks <- function() {
  data.frame(
    benchmark = character(),
    grade_count = integer(),
    discordance_rate = numeric(),
    reliability = numeric(),
    tau_bar = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.gp_frontier_score_grades <- function(scores, n_bins = NULL) {
  scores <- .gradepath_validate_numeric_vector(scores, "scores")
  if (any(!is.finite(scores))) {
    .gradepath_abort("Benchmark scores must be finite.")
  }

  order_scores <- order(scores, decreasing = TRUE, na.last = NA)
  if (length(order_scores) != length(scores)) {
    .gradepath_abort("Benchmark scores must not contain missing values.")
  }

  if (is.null(n_bins)) {
    levels <- sort(unique(scores), decreasing = TRUE)
    return(as.integer(match(scores, levels)))
  }

  n_bins <- .gradepath_validate_integerish(n_bins, "n_bins", min = 1L)
  n <- length(scores)
  grade_sorted <- pmin(n_bins, ceiling(seq_len(n) * n_bins / n))
  grades <- integer(n)
  grades[order_scores] <- grade_sorted
  as.integer(grades)
}

.gp_frontier_align_benchmark_vector <- function(value, name, ids) {
  value_names <- names(value)
  value <- .gradepath_validate_numeric_vector(
    value,
    paste0("benchmark_scores$", name)
  )
  if (length(value) != length(ids)) {
    .gradepath_abort(
      "`benchmark_scores$%s` must have one value per unit.",
      name
    )
  }
  if (!is.null(value_names)) {
    if (anyNA(value_names) || any(!nzchar(value_names)) ||
        !setequal(value_names, ids)) {
      .gradepath_abort(
        "`benchmark_scores$%s` names must match `pairwise$ids`.",
        name
      )
    }
    names(value) <- value_names
    value <- value[ids]
  }
  as.numeric(value)
}

.gp_frontier_validate_benchmark_scores <- function(benchmark_scores, ids) {
  if (is.null(benchmark_scores)) {
    return(NULL)
  }
  if (!is.list(benchmark_scores)) {
    .gradepath_abort("`benchmark_scores` must be a named list or data frame.")
  }
  benchmark_scores <- as.list(benchmark_scores)
  required <- c("raw_estimate", "posterior_mean", "linear_shrinkage")
  missing <- setdiff(required, names(benchmark_scores))
  if (length(missing) > 0L) {
    .gradepath_abort(
      "`benchmark_scores` must include: %s.",
      paste(required, collapse = ", ")
    )
  }
  lapply(required, function(name) {
    .gp_frontier_align_benchmark_vector(
      benchmark_scores[[name]],
      name = name,
      ids = ids
    )
  }) |>
    stats::setNames(required)
}

.gp_frontier_benchmarks <- function(pairwise, benchmark_scores = NULL) {
  pairwise <- validate_gp_pairwise(pairwise)
  benchmark_scores <- .gp_frontier_validate_benchmark_scores(
    benchmark_scores,
    pairwise$ids
  )
  if (is.null(benchmark_scores)) {
    return(.gp_frontier_empty_benchmarks())
  }

  specs <- list(
    raw_estimate_rank = .gp_frontier_score_grades(benchmark_scores$raw_estimate),
    posterior_mean_rank = .gp_frontier_score_grades(benchmark_scores$posterior_mean),
    linear_shrinkage_rank = .gp_frontier_score_grades(benchmark_scores$linear_shrinkage),
    posterior_mean_decile = .gp_frontier_score_grades(
      benchmark_scores$posterior_mean,
      n_bins = min(10L, length(pairwise$ids))
    ),
    posterior_mean_quartile = .gp_frontier_score_grades(
      benchmark_scores$posterior_mean,
      n_bins = min(4L, length(pairwise$ids))
    )
  )

  rows <- lapply(names(specs), function(name) {
    grades <- specs[[name]]
    metrics <- .gradepath_frontier_metrics(grades, pairwise$matrix)
    data.frame(
      benchmark = name,
      grade_count = as.integer(length(unique(grades))),
      discordance_rate = metrics$discordance_rate,
      reliability = metrics$reliability,
      tau_bar = metrics$tau_bar,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  do.call(rbind, rows)
}

.gp_frontier_dr_matrix <- function(grades, pairwise_matrix) {
  grades <- .gradepath_validate_numeric_vector(grades, "grades")
  pairwise_matrix <- .gradepath_validate_numeric_matrix(
    pairwise_matrix,
    "pairwise_matrix"
  )
  if (length(grades) != nrow(pairwise_matrix) ||
      length(grades) != ncol(pairwise_matrix)) {
    .gradepath_abort(
      "Conditional DR matrix requires one grade label per row and column of `pairwise_matrix`."
    )
  }
  if (any(abs(grades - round(grades)) > 1e-8) || any(grades < 1)) {
    .gradepath_abort("`grades` must be positive integer labels.")
  }

  grades <- as.integer(round(grades))
  ordered_grades <- sort(unique(grades))
  labels <- paste0("grade_", ordered_grades)
  out <- matrix(
    NA_real_,
    nrow = length(ordered_grades),
    ncol = length(ordered_grades),
    dimnames = list(labels, labels)
  )
  if (length(ordered_grades) < 2L) {
    return(out)
  }

  for (a in seq_along(ordered_grades)) {
    for (b in seq_along(ordered_grades)) {
      if (a >= b) {
        next
      }
      more_extreme <- min(ordered_grades[[a]], ordered_grades[[b]])
      less_extreme <- max(ordered_grades[[a]], ordered_grades[[b]])
      idx_more <- grades == more_extreme
      idx_less <- grades == less_extreme
      concordant <- sum(pairwise_matrix[idx_more, idx_less, drop = FALSE])
      discordant <- sum(pairwise_matrix[idx_less, idx_more, drop = FALSE])
      denom <- concordant + discordant
      row_index <- which(ordered_grades == less_extreme)
      col_index <- which(ordered_grades == more_extreme)
      out[row_index, col_index] <- if (denom > 0) discordant / denom else NA_real_
    }
  }

  out
}

.gp_frontier_required_table <- function(table, object = "frontier table") {
  table <- .gradepath_validate_data_frame(table, object)
  missing <- setdiff(.gp_frontier_table_columns, names(table))
  if (length(missing) > 0L) {
    .gradepath_abort(
      "`%s` must include columns: %s.",
      object,
      paste(.gp_frontier_table_columns, collapse = ", ")
    )
  }
  table
}

new_gp_frontier <- function(ids,
                            table,
                            benchmarks,
                            dr_matrix,
                            selection,
                            control,
                            schema_version = .gradepath_schema_version,
                            provenance = list(),
                            warnings = character()) {
  x <- structure(
    list(
      ids = ids,
      table = table,
      benchmarks = benchmarks,
      dr_matrix = dr_matrix,
      selection = selection,
      control = control,
      schema_version = schema_version,
      provenance = provenance,
      warnings = warnings
    ),
    class = c("gp_frontier", "list")
  )

  validate_gp_frontier(x)
}

validate_gp_frontier <- function(x) {
  .gradepath_validate_list_class(x, "gp_frontier")
  .gradepath_validate_named_fields(x, .gp_frontier_fields, "gp_frontier")

  ids <- .gradepath_validate_character_vector(
    x$ids,
    "gp_frontier$ids",
    unique = TRUE
  )
  table <- .gp_frontier_required_table(x$table, "gp_frontier$table")
  if (nrow(table) < 1L) {
    .gradepath_abort("`gp_frontier$table` must contain at least one row.")
  }
  for (name in c("lambda", "discordance_rate", "reliability", "tau_bar")) {
    table[[name]] <- .gradepath_validate_numeric_vector(
      table[[name]],
      paste0("gp_frontier$table$", name)
    )
  }
  table$grade_count <- .gradepath_validate_numeric_vector(
    table$grade_count,
    "gp_frontier$table$grade_count"
  )
  if (any(table$lambda < 0 | table$lambda > 1) ||
      any(table$discordance_rate < 0 | table$discordance_rate > 1) ||
      any(table$reliability < 0 | table$reliability > 1) ||
      any(table$tau_bar < -1 | table$tau_bar > 1) ||
      any(abs(table$grade_count - round(table$grade_count)) > 1e-8) ||
      any(table$grade_count < 1)) {
    .gradepath_abort("`gp_frontier$table` contains invalid frontier values.")
  }
  if (!isTRUE(all.equal(
    table$reliability,
    1 - table$discordance_rate,
    tolerance = 1e-8
  ))) {
    .gradepath_abort("`gp_frontier$table$reliability` must equal `1 - discordance_rate`.")
  }
  table$grade_count <- as.integer(round(table$grade_count))

  benchmarks <- .gradepath_validate_data_frame(
    x$benchmarks,
    "gp_frontier$benchmarks"
  )
  missing_benchmark <- setdiff(.gp_frontier_benchmark_columns, names(benchmarks))
  if (length(missing_benchmark) > 0L) {
    .gradepath_abort(
      "`gp_frontier$benchmarks` must include columns: %s.",
      paste(.gp_frontier_benchmark_columns, collapse = ", ")
    )
  }
  if (nrow(benchmarks) > 0L) {
    benchmarks$benchmark <- .gradepath_validate_character_vector(
      benchmarks$benchmark,
      "gp_frontier$benchmarks$benchmark",
      unique = TRUE
    )
    for (name in c("discordance_rate", "reliability", "tau_bar")) {
      benchmarks[[name]] <- .gradepath_validate_numeric_vector(
        benchmarks[[name]],
        paste0("gp_frontier$benchmarks$", name)
      )
    }
    benchmarks$grade_count <- .gradepath_validate_numeric_vector(
      benchmarks$grade_count,
      "gp_frontier$benchmarks$grade_count"
    )
    if (any(benchmarks$discordance_rate < 0 | benchmarks$discordance_rate > 1) ||
        any(benchmarks$reliability < 0 | benchmarks$reliability > 1) ||
        any(benchmarks$tau_bar < -1 | benchmarks$tau_bar > 1) ||
        any(abs(benchmarks$grade_count - round(benchmarks$grade_count)) > 1e-8) ||
        any(benchmarks$grade_count < 1)) {
      .gradepath_abort("`gp_frontier$benchmarks` contains invalid frontier values.")
    }
    if (!isTRUE(all.equal(
      benchmarks$reliability,
      1 - benchmarks$discordance_rate,
      tolerance = 1e-8
    ))) {
      .gradepath_abort("`gp_frontier$benchmarks$reliability` must equal `1 - discordance_rate`.")
    }
    benchmarks$grade_count <- as.integer(round(benchmarks$grade_count))
  }

  dr_matrix <- x$dr_matrix
  if (!is.matrix(dr_matrix) || !is.numeric(dr_matrix)) {
    .gradepath_abort("`gp_frontier$dr_matrix` must be a numeric matrix.")
  }
  if (!identical(nrow(dr_matrix), ncol(dr_matrix))) {
    .gradepath_abort("`gp_frontier$dr_matrix` must be square.")
  }
  finite_dr <- dr_matrix[is.finite(dr_matrix)]
  if (length(finite_dr) > 0L && any(finite_dr < 0 | finite_dr > 1)) {
    .gradepath_abort("`gp_frontier$dr_matrix` finite entries must lie in [0, 1].")
  }
  if (length(dr_matrix) > 0L) {
    if (any(!is.na(diag(dr_matrix))) ||
        any(!is.na(dr_matrix[upper.tri(dr_matrix)]))) {
      .gradepath_abort("`gp_frontier$dr_matrix` must store conditional DR values only in the lower triangle.")
    }
  }
  if (nrow(dr_matrix) > 0L) {
    if (any(!is.na(diag(dr_matrix)))) {
      .gradepath_abort("`gp_frontier$dr_matrix` diagonal entries must be `NA`.")
    }
    if (any(!is.na(dr_matrix[upper.tri(dr_matrix)]))) {
      .gradepath_abort("`gp_frontier$dr_matrix` must store conditional DR values in the lower triangle only.")
    }
  }

  selection <- .gradepath_validate_named_list(
    x$selection,
    "gp_frontier$selection"
  )
  .gradepath_validate_required_keys(
    selection,
    "selected_lambda",
    "gp_frontier$selection"
  )
  selection$selected_lambda <- .gp_grade_validate_lambda(
    selection$selected_lambda,
    "gp_frontier$selection$selected_lambda"
  )

  control <- validate_gp_control(x$control)
  schema_version <- .gradepath_validate_scalar_character(
    x$schema_version,
    "gp_frontier$schema_version",
    allowed = .gradepath_schema_version
  )
  provenance <- .gradepath_validate_named_list(
    x$provenance,
    "gp_frontier$provenance"
  )
  warnings <- .gradepath_validate_warning_vector(
    x$warnings,
    "gp_frontier$warnings"
  )

  structure(
    list(
      ids = ids,
      table = table,
      benchmarks = benchmarks,
      dr_matrix = dr_matrix,
      selection = selection,
      control = control,
      schema_version = schema_version,
      provenance = provenance,
      warnings = warnings
    ),
    class = c("gp_frontier", "list")
  )
}

#' Summarize the frontier of a solved grade path
#'
#' `gp_frontier()` collects the frontier metrics for a solved grade path: the
#' per-penalty table of grade count, discordance rate, reliability, and average
#' rank distance, and -- when you also supply the pairwise outranking profile --
#' the five naive benchmark scorings and the selected-penalty conditional
#' discordance matrix. Reach for it after [gp_grade_path()] when you want the
#' frontier as a tidy object to print, plot with [gp_plot_frontier()], or compare
#' grades against simple baselines.
#'
#' @param grade_path A solved `gp_grade_path` (from [gp_grade_path()] or the
#'   `grade_path` slot of a `gp_fit`). Supplies the penalty grid and the per-penalty
#'   grade assignments the frontier is built over.
#' @param pairwise Optional `gp_pairwise` (from [gp_pairwise()] or [get_pairwise()]).
#'   When `NULL` (default) `gp_frontier()` returns only the frontier table cached by
#'   [gp_grade_path()]; the benchmarks and conditional discordance matrix are left
#'   empty. Supplying it must match `grade_path$ids` and unlocks both. Its `ids` must
#'   equal `grade_path$ids`.
#' @param selected_lambda Optional numeric penalty at which the conditional
#'   discordance matrix is evaluated. `NULL` (default) uses
#'   `grade_path$selection$selected_lambda`; any value supplied must match a penalty
#'   on the grid exactly.
#' @param benchmark_scores Optional named list (or one-row-per-unit data frame) of
#'   the per-unit scores the naive baselines are built from. Requires three numeric
#'   vectors named `raw_estimate`, `posterior_mean`, and `linear_shrinkage`, each one
#'   value per unit and named by (or column-aligned to) the unit ids in
#'   `pairwise$ids`; they are reordered to that id order. From these, five baseline
#'   grade scorings are formed -- a rank cut on each of the three vectors plus decile
#'   and quartile bins of `posterior_mean` -- and scored against the pairwise profile.
#'   `NULL` (default) leaves the benchmark table empty. Requires `pairwise`; passing
#'   scores without it raises an error.
#'
#' @return A validated `gp_frontier` object (a list of class
#'   `c("gp_frontier", "list")`) with the public slots: \describe{
#'   \item{`ids`}{Character vector; the canonical unit-id order, equal to
#'     `grade_path$ids`.}
#'   \item{`table`}{Data frame, one row per penalty on the grid, with `lambda`,
#'     `grade_count`, `discordance_rate`, `reliability`, and `tau_bar` (the average
#'     between-grade rank distance), plus the solve `objective`/`status` when a
#'     `pairwise` profile is supplied.}
#'   \item{`benchmarks`}{Data frame of the five naive baseline scorings
#'     (`raw_estimate_rank`, `posterior_mean_rank`, `linear_shrinkage_rank`,
#'     `posterior_mean_decile`, `posterior_mean_quartile`), each with its
#'     `grade_count`, `discordance_rate`, `reliability`, and `tau_bar`; empty when
#'     `benchmark_scores` (and `pairwise`) were not supplied.}
#'   \item{`dr_matrix`}{Lower-triangular matrix of conditional discordance rates
#'     between grade pairs at the selected penalty (row = less-extreme grade,
#'     column = more-extreme grade; diagonal and upper triangle are `NA`); a
#'     0 x 0 matrix when no `pairwise` profile was supplied.}
#'   \item{`selection`}{Named list recording `selected_lambda`, the penalty the
#'     conditional discordance matrix was evaluated at.}
#'   \item{`control`}{The [gp_control] carried from `grade_path`.}
#'   \item{`provenance`, `warnings`, `schema_version`}{Producer metadata and
#'     internal audit slots; `warnings` notes when benchmarks and the conditional
#'     discordance matrix were skipped for lack of a `pairwise` profile.}
#' }
#'
#' @details The frontier metrics are computed over between-grade blocks of the
#'   pairwise outranking profile and normalized by `2 * choose(n, 2)`. The
#'   `dr_matrix` reports, for each pair of grades, the share of cross-grade
#'   comparisons that run against the grade ordering (discordant), conditional on
#'   that pair. Grade labels are integers and carry no ranking-superiority statement;
#'   the benchmarks are descriptive baselines, not a contest.
#'
#' @examples
#' # Inspect the frontier of the bundled pre-solved fit (no Gurobi, no solve).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' post <- get_posterior(fit)
#' ids <- get_pairwise(fit)$ids
#'
#' fr <- gp_frontier(
#'   fit$grade_path,
#'   pairwise = get_pairwise(fit),
#'   benchmark_scores = list(
#'     raw_estimate     = setNames(post$estimate, ids),
#'     posterior_mean   = setNames(post$posterior_mean, ids),
#'     linear_shrinkage = setNames(0.5 * post$estimate, ids)
#'   )
#' )
#' fr$table        # per-penalty frontier metrics
#' fr$benchmarks   # the five naive baseline scorings
#'
#' @seealso [gp_grade_path()], [gp_pairwise()], [gp_plot_frontier()],
#'   [gp_plot_discordance()], [get_pairwise()]
#' @family gradepath-frontier
#' @export
gp_frontier <- function(grade_path,
                        pairwise = NULL,
                        selected_lambda = NULL,
                        benchmark_scores = NULL) {
  grade_path <- validate_gp_grade_path(grade_path)
  if (is.null(selected_lambda)) {
    selected_lambda <- grade_path$selection$selected_lambda
  }
  selected <- .gp_grade_exact_lambda_match(
    grade_path$lambda_grid,
    selected_lambda,
    "selected_lambda"
  )
  pairwise <- if (is.null(pairwise)) NULL else validate_gp_pairwise(pairwise)

  table <- if (is.null(pairwise)) {
    .gp_frontier_required_table(grade_path$summary, "grade_path$summary")
  } else {
    frontier_table(pairwise, grade_path$fits)
  }
  table <- table[, intersect(names(table), c(
    .gp_frontier_table_columns,
    "objective",
    "raw_objective",
    "status"
  )), drop = FALSE]

  if (is.null(pairwise)) {
    if (!is.null(benchmark_scores)) {
      .gradepath_abort("`benchmark_scores` require an explicit `pairwise` profile.")
    }
    benchmarks <- .gp_frontier_empty_benchmarks()
    dr_matrix <- matrix(NA_real_, nrow = 0L, ncol = 0L)
    warnings <- "Pairwise profile unavailable; benchmarks and conditional DR matrix were not computed."
  } else {
    if (!identical(pairwise$ids, grade_path$ids)) {
      .gradepath_abort("`pairwise$ids` must match `grade_path$ids`.")
    }
    benchmarks <- .gp_frontier_benchmarks(pairwise, benchmark_scores)
    dr_matrix <- .gp_frontier_dr_matrix(
      grades = grade_path$fits[[selected$index]]$assignment$grade,
      pairwise_matrix = pairwise$matrix
    )
    warnings <- character()
  }

  new_gp_frontier(
    ids = grade_path$ids,
    table = table,
    benchmarks = benchmarks,
    dr_matrix = dr_matrix,
    selection = list(selected_lambda = selected$value),
    control = grade_path$control,
    provenance = .gradepath_new_provenance(
      producer = "gp_frontier",
      pairwise_available = !is.null(pairwise)
    ),
    warnings = warnings
  )
}
