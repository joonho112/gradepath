# =============================================================================
# grade.R -- grade IP problem object and grade extraction helpers
# -----------------------------------------------------------------------------
# The grade IP uses a row-major decision-vector layout:
#   d_11, d_12, ..., d_1n, d_21, ..., d_nn, g_1, ..., g_n.
# This mirrors the KRW/v1 solver layout and keeps the MILP raw objective equal
# to one half of the canonical reported risk.
# =============================================================================

.gp_grade_problem_fields <- c(
  "ids",
  "n_units",
  "lambda",
  "d_indices",
  "grade_indices",
  "objective",
  "constraint_i",
  "constraint_j",
  "constraint_x",
  "constraint_nrow",
  "constraint_ncol",
  "direction",
  "rhs",
  "types",
  "lower",
  "upper",
  "diagonal_positions",
  "m_big",
  "encoding",
  "metadata"
)

.gp_grade_problem_encodings <- c("bigM", "indicator")

.gradepath_abort_grade <- function(message, ..., call = FALSE) {
  .gradepath_abort(
    message,
    ...,
    class = c("gradepath_grade_error", "gradepath_error"),
    call = call
  )
}

.gp_grade_validate_pairwise <- function(pairwise) {
  # The public gp_pairwise validator pins diag(pi) = 0.5. The grade IP ignores
  # the diagonal, so normalize only this boundary before validation to make the
  # solver fixture explicitly diagonal-invariant.
  if (is.list(pairwise) &&
      !is.null(pairwise$matrix) &&
      is.matrix(pairwise$matrix)) {
    pairwise$matrix <- as.matrix(pairwise$matrix)
    diag(pairwise$matrix) <- .gp_pairwise_diagonal
  }
  validate_gp_pairwise(pairwise)
}

.gp_grade_validate_lambda <- function(lambda, name = "lambda") {
  .gradepath_validate_scalar_numeric(
    lambda,
    name,
    lower = 0,
    upper = 1,
    include_lower = TRUE,
    include_upper = TRUE
  )
}

.gp_grade_matrix_for_optimization <- function(pairwise) {
  pairwise <- .gp_grade_validate_pairwise(pairwise)
  matrix_opt <- pairwise$matrix
  diag(matrix_opt) <- 0
  storage.mode(matrix_opt) <- "double"
  matrix_opt
}

.gp_grade_problem_indices <- function(n_units) {
  n_units <- .gradepath_validate_integerish(n_units, "n_units", min = 1L)
  d_indices <- matrix(
    seq_len(n_units * n_units),
    nrow = n_units,
    ncol = n_units,
    byrow = TRUE
  )
  grade_indices <- (n_units * n_units) + seq_len(n_units)

  list(d_indices = d_indices, grade_indices = grade_indices)
}

#' Build the grade IP problem object
#'
#' Internal constructor for the solver-facing grade problem. It projects
#' the pairwise diagonal out by setting both the diagonal objective coefficients
#' and diagonal upper bounds to zero, and uses integer grade labels in `[1, n]`
#' with big-M `M = n`.
#'
#' @keywords internal
#' @noRd
gp_grade_problem <- function(pairwise, lambda, encoding = "bigM") {
  pairwise <- .gp_grade_validate_pairwise(pairwise)
  lambda <- .gp_grade_validate_lambda(lambda)
  encoding <- .gradepath_validate_scalar_character(
    encoding,
    "encoding",
    allowed = .gp_grade_problem_encodings
  )

  matrix_opt <- .gp_grade_matrix_for_optimization(pairwise)
  ids <- pairwise$ids
  n_units <- length(ids)
  if (n_units < 2L) {
    .gradepath_abort_grade("`gp_grade_problem()` requires at least two ids.")
  }
  index <- .gp_grade_problem_indices(n_units)
  d_indices <- index$d_indices
  grade_indices <- index$grade_indices
  n_decision <- n_units * n_units
  n_total <- n_decision + n_units
  diag_pos <- gp_diag_positions(n_units)
  m_big <- as.numeric(n_units)

  pair_i <- rep(seq_len(n_units), each = n_units)
  pair_j <- rep(seq_len(n_units), times = n_units)
  off_diagonal <- pair_i != pair_j
  pair_i <- pair_i[off_diagonal]
  pair_j <- pair_j[off_diagonal]
  pair_d_index <- d_indices[cbind(pair_i, pair_j)]
  pair_count <- length(pair_d_index)

  row_upper <- seq_len(pair_count)
  row_lower <- pair_count + seq_len(pair_count)

  constraint_i <- c(
    row_upper,
    row_upper,
    row_upper,
    row_lower,
    row_lower,
    row_lower
  )
  constraint_j <- c(
    pair_d_index,
    grade_indices[pair_i],
    grade_indices[pair_j],
    pair_d_index,
    grade_indices[pair_i],
    grade_indices[pair_j]
  )
  constraint_x <- c(
    rep(-m_big, pair_count),
    rep(1, pair_count),
    rep(-1, pair_count),
    rep(m_big, pair_count),
    rep(-1, pair_count),
    rep(1, pair_count)
  )

  objective <- c(
    as.vector(t(t(matrix_opt) - (lambda * matrix_opt))),
    rep(0, n_units)
  )
  objective[diag_pos] <- 0

  lower <- c(rep(0, n_decision), rep(1, n_units))
  upper <- c(rep(1, n_decision), rep(n_units, n_units))
  upper[diag_pos] <- 0

  problem <- list(
    ids = ids,
    n_units = n_units,
    lambda = lambda,
    d_indices = d_indices,
    grade_indices = grade_indices,
    objective = objective,
    constraint_i = constraint_i,
    constraint_j = constraint_j,
    constraint_x = constraint_x,
    constraint_nrow = 2L * pair_count,
    constraint_ncol = n_total,
    direction = rep("<=", 2L * pair_count),
    rhs = c(rep(0, pair_count), rep(m_big - 1, pair_count)),
    types = c(rep("B", n_decision), rep("I", n_units)),
    lower = lower,
    upper = upper,
    diagonal_positions = diag_pos,
    m_big = m_big,
    encoding = encoding,
    metadata = .gradepath_new_provenance(
      producer = "gp_grade_problem",
      objective_layout = "row_major_D_then_grade",
      constraint_encoding = "bigM_canonical_rows",
      requested_solver_encoding = encoding,
      diagonal_projection = "upper_bound_zero_and_objective_zero",
      label_bounds = "[1,n]"
    )
  )

  validate_gp_grade_problem(problem)
}

validate_gp_grade_problem <- function(problem) {
  if (!is.list(problem)) {
    .gradepath_abort_grade("`problem` must be a list.")
  }
  .gradepath_validate_named_fields(
    problem,
    .gp_grade_problem_fields,
    "gp_grade_problem"
  )

  ids <- .gradepath_validate_character_vector(
    problem$ids,
    "gp_grade_problem$ids",
    unique = TRUE
  )
  n_units <- .gradepath_validate_integerish(
    problem$n_units,
    "gp_grade_problem$n_units",
    min = 1L
  )
  if (!identical(n_units, length(ids))) {
    .gradepath_abort_grade("`gp_grade_problem$n_units` must equal `length(ids)`.")
  }
  if (n_units < 2L) {
    .gradepath_abort_grade("`gp_grade_problem$n_units` must be at least 2.")
  }
  lambda <- .gp_grade_validate_lambda(problem$lambda, "gp_grade_problem$lambda")
  d_indices <- .gradepath_validate_numeric_matrix(
    problem$d_indices,
    "gp_grade_problem$d_indices"
  )
  if (!identical(dim(d_indices), c(n_units, n_units))) {
    .gradepath_abort_grade("`gp_grade_problem$d_indices` must be n x n.")
  }
  expected_d <- matrix(seq_len(n_units * n_units), n_units, n_units, byrow = TRUE)
  if (!identical(d_indices, expected_d)) {
    .gradepath_abort_grade("`gp_grade_problem$d_indices` must be row-major.")
  }

  grade_indices <- .gradepath_validate_numeric_vector(
    problem$grade_indices,
    "gp_grade_problem$grade_indices"
  )
  expected_grade <- (n_units * n_units) + seq_len(n_units)
  if (!identical(as.integer(grade_indices), expected_grade)) {
    .gradepath_abort_grade("`gp_grade_problem$grade_indices` are not contiguous after D.")
  }

  n_total <- n_units * n_units + n_units
  objective <- .gradepath_validate_numeric_vector(
    problem$objective,
    "gp_grade_problem$objective"
  )
  if (length(objective) != n_total) {
    .gradepath_abort_grade("`gp_grade_problem$objective` has the wrong length.")
  }

  constraint_i <- .gradepath_validate_numeric_vector(
    problem$constraint_i,
    "gp_grade_problem$constraint_i"
  )
  constraint_j <- .gradepath_validate_numeric_vector(
    problem$constraint_j,
    "gp_grade_problem$constraint_j"
  )
  constraint_x <- .gradepath_validate_numeric_vector(
    problem$constraint_x,
    "gp_grade_problem$constraint_x"
  )
  if (!identical(length(constraint_i), length(constraint_j)) ||
      !identical(length(constraint_i), length(constraint_x))) {
    .gradepath_abort_grade("Constraint triplets must have equal length.")
  }
  constraint_nrow <- .gradepath_validate_integerish(
    problem$constraint_nrow,
    "gp_grade_problem$constraint_nrow",
    min = 0L
  )
  constraint_ncol <- .gradepath_validate_integerish(
    problem$constraint_ncol,
    "gp_grade_problem$constraint_ncol",
    min = 1L
  )
  if (!identical(constraint_ncol, n_total)) {
    .gradepath_abort_grade("`gp_grade_problem$constraint_ncol` must equal n^2 + n.")
  }
  if (any(constraint_i < 1 | constraint_i > constraint_nrow) ||
      any(constraint_j < 1 | constraint_j > constraint_ncol)) {
    .gradepath_abort_grade("Constraint triplet indices are out of bounds.")
  }

  direction <- .gradepath_validate_character_vector(
    problem$direction,
    "gp_grade_problem$direction"
  )
  if (length(direction) != constraint_nrow || any(direction != "<=")) {
    .gradepath_abort_grade("`gp_grade_problem$direction` must be '<=' for every row.")
  }
  rhs <- .gradepath_validate_numeric_vector(problem$rhs, "gp_grade_problem$rhs")
  if (length(rhs) != constraint_nrow) {
    .gradepath_abort_grade("`gp_grade_problem$rhs` has the wrong length.")
  }

  types <- .gradepath_validate_character_vector(
    problem$types,
    "gp_grade_problem$types"
  )
  if (length(types) != n_total ||
      !identical(types[seq_len(n_units * n_units)], rep("B", n_units * n_units)) ||
      !identical(types[(n_units * n_units + 1L):n_total], rep("I", n_units))) {
    .gradepath_abort_grade("`gp_grade_problem$types` must be binary D then integer grades.")
  }

  lower <- .gradepath_validate_numeric_vector(
    problem$lower,
    "gp_grade_problem$lower"
  )
  upper <- .gradepath_validate_numeric_vector(
    problem$upper,
    "gp_grade_problem$upper"
  )
  if (length(lower) != n_total || length(upper) != n_total) {
    .gradepath_abort_grade("Problem bounds have the wrong length.")
  }

  diagonal_positions <- .gradepath_validate_numeric_vector(
    problem$diagonal_positions,
    "gp_grade_problem$diagonal_positions"
  )
  if (!identical(as.integer(diagonal_positions), gp_diag_positions(n_units))) {
    .gradepath_abort_grade("`gp_grade_problem$diagonal_positions` are invalid.")
  }
  diagonal_positions <- as.integer(diagonal_positions)
  if (any(lower[diagonal_positions] != 0) ||
      any(upper[diagonal_positions] != 0) ||
      any(objective[diagonal_positions] != 0)) {
    .gradepath_abort_grade("The diagonal decision variables must be projected to zero.")
  }

  d_slots <- seq_len(n_units * n_units)
  off_diagonal_slots <- setdiff(d_slots, diagonal_positions)
  if (any(lower[d_slots] != 0) || any(upper[off_diagonal_slots] != 1)) {
    .gradepath_abort_grade("Off-diagonal D bounds must be [0, 1] and diagonal D bounds must be zero.")
  }

  grade_slots <- (n_units * n_units + 1L):n_total
  if (any(lower[grade_slots] != 1) ||
      any(upper[grade_slots] != n_units)) {
    .gradepath_abort_grade("Grade-label bounds must be [1, n].")
  }

  m_big <- .gradepath_validate_scalar_numeric(
    problem$m_big,
    "gp_grade_problem$m_big",
    finite = TRUE,
    lower = 1,
    include_lower = TRUE
  )
  if (abs(m_big - n_units) > 1e-8) {
    .gradepath_abort_grade("`gp_grade_problem$m_big` must equal n.")
  }
  encoding <- .gradepath_validate_scalar_character(
    problem$encoding,
    "gp_grade_problem$encoding",
    allowed = .gp_grade_problem_encodings
  )
  metadata <- .gradepath_validate_named_list(
    problem$metadata,
    "gp_grade_problem$metadata"
  )

  structure(
    list(
      ids = ids,
      n_units = n_units,
      lambda = lambda,
      d_indices = d_indices,
      grade_indices = as.integer(grade_indices),
      objective = objective,
      constraint_i = as.integer(constraint_i),
      constraint_j = as.integer(constraint_j),
      constraint_x = constraint_x,
      constraint_nrow = constraint_nrow,
      constraint_ncol = constraint_ncol,
      direction = direction,
      rhs = rhs,
      types = types,
      lower = lower,
      upper = upper,
      diagonal_positions = diagonal_positions,
      m_big = m_big,
      encoding = encoding,
      metadata = metadata
    ),
    class = c("gp_grade_problem", "list")
  )
}

.gp_problem_require_digest <- function(caller) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    .gradepath_abort_grade("`%s()` requires the package `digest`.", caller)
  }
}

.gp_problem_num_key <- function(x) {
  sprintf("%.12g", round(as.numeric(x), 12L))
}

.gp_problem_variable_index <- function(problem) {
  data.frame(
    kind = c(rep("D", problem$n_units * problem$n_units), rep("g", problem$n_units)),
    i = c(rep(seq_len(problem$n_units), each = problem$n_units), seq_len(problem$n_units)),
    j = c(rep(seq_len(problem$n_units), times = problem$n_units), rep(NA_integer_, problem$n_units)),
    index = seq_len(problem$constraint_ncol)
  )
}

.gp_problem_constraint_row_key <- function(problem, row) {
  idx <- which(problem$constraint_i == row)
  if (length(idx) > 0L) {
    ord <- order(problem$constraint_j[idx], problem$constraint_x[idx])
    idx <- idx[ord]
  }
  terms <- if (length(idx) == 0L) {
    ""
  } else {
    paste(
      paste(problem$constraint_j[idx], .gp_problem_num_key(problem$constraint_x[idx]), sep = ":"),
      collapse = ","
    )
  }
  paste(
    problem$direction[[row]],
    .gp_problem_num_key(problem$rhs[[row]]),
    terms,
    sep = "|"
  )
}

.gp_problem_feasible_core <- function(problem) {
  list(
    n_units = problem$n_units,
    n_constraints = problem$constraint_nrow,
    n_variables = problem$constraint_ncol,
    variable_index = .gp_problem_variable_index(problem),
    d_indices = as.integer(problem$d_indices),
    grade_indices = as.integer(problem$grade_indices),
    constraint_rows = sort(vapply(
      seq_len(problem$constraint_nrow),
      function(row) .gp_problem_constraint_row_key(problem, row),
      character(1)
    )),
    lower = .gp_problem_num_key(problem$lower),
    upper = .gp_problem_num_key(problem$upper),
    types = problem$types,
    diagonal_positions = as.integer(problem$diagonal_positions)
  )
}

gp_problem_feasible_signature <- function(problem) {
  problem <- validate_gp_grade_problem(problem)
  .gp_problem_require_digest("gp_problem_feasible_signature")
  digest::digest(.gp_problem_feasible_core(problem), algo = "sha256")
}

gp_problem_hash <- function(problem) {
  problem <- validate_gp_grade_problem(problem)
  .gp_problem_require_digest("gp_problem_hash")

  core <- list(
    feasible_signature = gp_problem_feasible_signature(problem),
    lambda = .gp_problem_num_key(problem$lambda),
    objective = .gp_problem_num_key(problem$objective)
  )

  digest::digest(core, algo = "sha256")
}

.gradepath_grade_decision_matrix <- function(problem, primal, tol = 1e-6) {
  problem <- validate_gp_grade_problem(problem)
  primal <- .gradepath_validate_numeric_vector(primal, "primal")
  tol <- .gradepath_validate_scalar_numeric(
    tol,
    "tol",
    finite = TRUE,
    lower = 0,
    include_lower = FALSE
  )

  if (length(primal) != problem$constraint_ncol) {
    .gradepath_abort_grade(
      "Solver output must have length n^2 + n for the grade problem."
    )
  }

  d_length <- problem$n_units * problem$n_units
  d_values <- primal[seq_len(d_length)]
  if (any(abs(d_values - round(d_values)) > tol) ||
      any(round(d_values) < 0 | round(d_values) > 1)) {
    .gradepath_abort_grade(
      "Solver returned non-binary decision values for the grading relation."
    )
  }

  decision_matrix <- matrix(
    as.integer(round(d_values)),
    nrow = problem$n_units,
    ncol = problem$n_units,
    byrow = TRUE,
    dimnames = list(problem$ids, problem$ids)
  )

  if (any(diag(decision_matrix) != 0L)) {
    .gradepath_abort_grade("Solver returned nonzero diagonal decisions.")
  }

  decision_matrix
}

.gradepath_grade_assignment <- function(decision_matrix, ids) {
  decision_matrix <- .gradepath_validate_numeric_matrix(
    decision_matrix,
    "decision_matrix"
  )
  ids <- .gradepath_validate_character_vector(ids, "ids", unique = TRUE)

  if (!identical(dim(decision_matrix), c(length(ids), length(ids)))) {
    .gradepath_abort_grade("`decision_matrix` must be square with one row and column per id.")
  }

  if (any(abs(decision_matrix - round(decision_matrix)) > 1e-8) ||
      any(decision_matrix < 0 | decision_matrix > 1)) {
    .gradepath_abort_grade("`decision_matrix` must contain binary decisions.")
  }
  decision_matrix <- matrix(
    as.integer(round(decision_matrix)),
    nrow = nrow(decision_matrix),
    ncol = ncol(decision_matrix),
    dimnames = dimnames(decision_matrix)
  )
  if (any(diag(decision_matrix) != 0)) {
    .gradepath_abort_grade("`decision_matrix` must have zero diagonal decisions.")
  }

  row_score <- rowSums(decision_matrix)
  score_levels <- sort(unique(row_score), decreasing = TRUE)
  grade <- as.integer(match(row_score, score_levels))

  data.frame(
    id = ids,
    grade = grade,
    stringsAsFactors = FALSE
  )
}

.gradepath_grade_objective_value <- function(matrix_opt, decision_matrix, lambda) {
  matrix_opt <- .gradepath_validate_numeric_matrix(matrix_opt, "matrix_opt")
  decision_matrix <- .gradepath_validate_numeric_matrix(
    decision_matrix,
    "decision_matrix"
  )
  lambda <- .gp_grade_validate_lambda(lambda)

  if (!identical(dim(matrix_opt), dim(decision_matrix))) {
    .gradepath_abort_grade("`matrix_opt` and `decision_matrix` must have the same dimensions.")
  }
  if (!identical(nrow(matrix_opt), ncol(matrix_opt))) {
    .gradepath_abort_grade("`matrix_opt` and `decision_matrix` must be square.")
  }
  if (any(abs(decision_matrix - round(decision_matrix)) > 1e-8) ||
      any(decision_matrix < 0 | decision_matrix > 1)) {
    .gradepath_abort_grade("`decision_matrix` must contain binary decisions.")
  }
  decision_matrix <- matrix(
    as.integer(round(decision_matrix)),
    nrow = nrow(decision_matrix),
    ncol = ncol(decision_matrix),
    dimnames = dimnames(decision_matrix)
  )
  if (any(diag(decision_matrix) != 0) || any(diag(matrix_opt) != 0)) {
    .gradepath_abort_grade("The canonical objective requires projected zero diagonals.")
  }

  dp_value <- sum((matrix_opt * t(decision_matrix)) + (t(matrix_opt) * decision_matrix))
  tau_value <- sum(
    (matrix_opt * decision_matrix) +
      (t(matrix_opt) * t(decision_matrix)) -
      (matrix_opt * t(decision_matrix)) -
      (t(matrix_opt) * decision_matrix)
  )

  ((1 - lambda) * dp_value) - (lambda * tau_value)
}

.gp_grade_operational_default_grid <- function() {
  sort(unique(c(seq(0, 1, by = 0.1), .gp_control_required_lambda)))
}

.gp_grade_exact_lambda_match <- function(lambda_grid,
                                         lambda,
                                         lambda_name = "lambda") {
  lambda <- .gp_grade_validate_lambda(lambda, lambda_name)
  distance <- abs(lambda_grid - lambda)
  hits <- which(distance < 1e-8)
  if (length(hits) != 1L) {
    .gradepath_abort_grade(
      "`%s` must match exactly one solved lambda. Available values: %s.",
      lambda_name,
      paste(format(lambda_grid), collapse = ", ")
    )
  }
  list(index = hits[[1L]], value = as.numeric(lambda_grid[[hits[[1L]]]]))
}

.gp_grade_resolve_control <- function(pairwise, control = NULL, lambda_grid = NULL) {
  pairwise <- validate_gp_pairwise(pairwise)
  control <- if (is.null(control)) pairwise$control else validate_gp_control(control)

  grid <- if (is.null(lambda_grid)) {
    if (identical(as.numeric(control$lambda_grid), as.numeric(.gp_control_default_lambda_grid))) {
      .gp_grade_operational_default_grid()
    } else {
      control$lambda_grid
    }
  } else {
    .gp_control_validate_lambda_grid(lambda_grid, "lambda_grid")
  }

  control$lambda_grid <- grid
  validate_gp_control(control)
}

.gp_grade_solve_order <- function(lambda_grid) {
  endpoint_index <- .gp_grade_exact_lambda_match(
    lambda_grid,
    1,
    "endpoint_lambda"
  )$index
  c(endpoint_index, rev(setdiff(seq_along(lambda_grid), endpoint_index)))
}

.gp_grade_solve_backend <- function(pairwise,
                                    lambda,
                                    control,
                                    matrix_opt,
                                    warm_start = NULL,
                                    acceptance_mode = FALSE) {
  problem <- gp_grade_problem(
    pairwise,
    lambda = lambda,
    encoding = if (identical(control$backend, "gurobi")) "indicator" else "bigM"
  )

  if (identical(control$backend, "gurobi")) {
    gp_gurobi_robust(
      problem,
      control = control,
      warm_start = warm_start,
      matrix_opt = matrix_opt,
      acceptance_mode = acceptance_mode
    )
  } else {
    # Open backends solve a single invocation directly; there is no path
    # cascade to continue, so `acceptance_mode` does not apply to them.
    gp_open_solve(
      problem,
      backend = control$backend,
      control = control,
      matrix_opt = matrix_opt
    )
  }
}

.gp_grade_fit_from_result <- function(result,
                                      ids,
                                      lambda,
                                      control,
                                      warm_start_from = NA_real_,
                                      warm_start_from_status = NA_character_) {
  assignment <- result$validation$assignment
  grade_count <- as.integer(length(unique(assignment$grade)))
  solver_meta <- result$solver_meta
  if (!is.numeric(warm_start_from) || length(warm_start_from) != 1L) {
    .gradepath_abort_grade("`warm_start_from` must be a length-1 numeric value or NA.")
  }
  warm_start_from <- if (is.na(warm_start_from)) {
    NA_real_
  } else {
    .gradepath_validate_scalar_numeric(
      warm_start_from,
      "warm_start_from",
      lower = 0,
      upper = 1,
      include_lower = TRUE,
      include_upper = TRUE
    )
  }
  warm_start_used <- !is.na(warm_start_from) &&
    !identical(solver_meta$warmstart, "none")
  warm_start_from_status <- if (isTRUE(warm_start_used)) {
    .gradepath_validate_scalar_character(
      warm_start_from_status,
      "warm_start_from_status"
    )
  } else {
    NA_character_
  }
  warm_start_from_acceptance_ready <- if (isTRUE(warm_start_used)) {
    .gp_status_acceptance_ready(warm_start_from_status)
  } else {
    NA
  }
  warnings <- character()
  if (!identical(result$status, "optimal")) {
    warnings <- sprintf(
      "Solver returned status `%s` at lambda %s.",
      result$status,
      format(lambda)
    )
  }

  new_gp_grade_fit(
    ids = ids,
    lambda = lambda,
    assignment = assignment,
    summary = list(
      grade_count = grade_count,
      status = result$status,
      n_units = length(ids)
    ),
    objective = list(
      value = result$validation$canon,
      raw = result$objval,
      canonical = result$validation$canon
    ),
    backend = list(
      name = solver_meta$backend,
      path = solver_meta$path,
      encoding = solver_meta$encoding,
      status = result$status,
      raw_solver_status = solver_meta$solver_status %gp_or% NA_character_,
      objbound = solver_meta$objbound,
      mipgap = solver_meta$mipgap,
      runtime = solver_meta$runtime,
      acceptance_mode = solver_meta$acceptance_mode %gp_or% FALSE,
      solver_attempts = solver_meta$attempts %gp_or% NULL,
      warmstart = solver_meta$warmstart,
      warm_start_from_lambda = warm_start_from,
      warm_start_from_status = warm_start_from_status,
      warm_start_from_acceptance_ready = warm_start_from_acceptance_ready,
      warm_start_used = warm_start_used,
      problem_hash = solver_meta$problem_hash,
      problem_signature = solver_meta$problem_signature
    ),
    control = control,
    provenance = .gradepath_new_provenance(
      producer = "gp_grade",
      lambda = lambda,
      problem_hash = solver_meta$problem_hash,
      problem_signature = solver_meta$problem_signature,
      warm_start_from_lambda = warm_start_from,
      warm_start_from_status = warm_start_from_status,
      warm_start_from_acceptance_ready = warm_start_from_acceptance_ready,
      warm_start_used = warm_start_used
    ),
    warnings = warnings
  )
}

.gp_grade_path_warm_start_sources <- function(lambda_grid, fits) {
  data.frame(
    lambda = lambda_grid,
    warm_start_from_lambda = vapply(
      fits,
      function(fit) fit$backend$warm_start_from_lambda,
      numeric(1)
    ),
    warm_start_from_status = vapply(
      fits,
      function(fit) {
        status <- fit$backend$warm_start_from_status
        if (is.null(status) || length(status) == 0L || is.na(status[[1L]])) {
          NA_character_
        } else {
          as.character(status[[1L]])
        }
      },
      character(1)
    ),
    warm_start_from_acceptance_ready = vapply(
      fits,
      function(fit) {
        ready <- fit$backend$warm_start_from_acceptance_ready
        if (is.null(ready) || length(ready) == 0L || is.na(ready[[1L]])) {
          NA
        } else {
          isTRUE(ready[[1L]])
        }
      },
      logical(1)
    ),
    warm_start_used = vapply(
      fits,
      function(fit) isTRUE(fit$backend$warm_start_used),
      logical(1)
    ),
    stringsAsFactors = FALSE
  )
}

.gp_grade_path_summary <- function(lambda_grid, fits) {
  data.frame(
    lambda = lambda_grid,
    grade_count = vapply(fits, function(fit) fit$summary$grade_count, integer(1)),
    objective = vapply(fits, function(fit) fit$objective$value, numeric(1)),
    raw_objective = vapply(fits, function(fit) fit$objective$raw, numeric(1)),
    status = vapply(fits, function(fit) fit$summary$status, character(1)),
    stringsAsFactors = FALSE
  )
}

.gp_grade_path_backend <- function(control, fits) {
  list(
    name = control$backend,
    paths = unique(vapply(fits, function(fit) fit$backend$path, character(1))),
    encodings = unique(vapply(fits, function(fit) fit$backend$encoding, character(1)))
  )
}

.gp_grade_path_warnings <- function(fits) {
  warnings <- unique(unlist(lapply(fits, `[[`, "warnings"), use.names = FALSE))
  if (is.null(warnings)) character() else warnings
}

#' Solve the grade path over a lambda grid
#'
#' `gp_grade_path()` solves one grade integer program per penalty value on a
#' lambda grid, tracing the frontier between fewer grades and lower discordance.
#' It is the engine the one-level pipeline calls; reach for it directly to obtain
#' the whole solved path, then pull a single assignment with [gp_select_grade()]
#' (solve-free) or solve one penalty in isolation with [gp_grade()]. The
#' Condorcet endpoint (`lambda = 1`) is solved first as a warm-start anchor, the
#' remaining penalties are solved in descending order, and fits are stored back
#' in ascending grid order.
#'
#' @param pairwise A [gp_pairwise] object: the `J x J` posterior outranking
#'   structure the integer program scores.
#' @param lambda_grid Numeric vector of penalties in `[0, 1]`, or `NULL`. `NULL`
#'   (default) uses the operational small grid unless `pairwise` or `control`
#'   already carries a custom grid. An explicit grid must include the parity
#'   anchors `0.25` and `1`.
#' @param control Optional [gp_control]; the backend, solver tolerances, and
#'   active lambda grid. Defaults to `pairwise$control`.
#' @param selected_lambda Numeric scalar in the grid; the penalty recorded as the
#'   default reportable selection on the returned path. Defaults to `0.25` (KRW's
#'   published selection). Recording it never triggers a solve.
#' @param selection_rule Character scalar; a label stored alongside
#'   `selected_lambda` in the path's selection metadata. Defaults to
#'   `"baseline_lambda_0.25"`.
#' @param acceptance_mode Logical; the solver solution-quality policy, honored by
#'   the gurobi backend only (open backends solve a single invocation, so it is
#'   ignored there). `FALSE` (default) keeps the first Gurobi path that returns a
#'   normalized result, so a selected solve that stops at `gap_reached` /
#'   `time_limit` is reported honestly with that status (it lands `UNVERIFIED`
#'   downstream). `TRUE` additionally attempts later Gurobi paths to prove the
#'   optimum when an earlier path is not acceptance-ready; it stays honest and
#'   never relabels a non-optimal solve as `optimal`. The chosen value is recorded
#'   on each fit's `backend$acceptance_mode`.
#'
#' @return A validated [gp_grade_path] object (a list of class
#'   `c("gp_grade_path", "list")`) with the public slots: \describe{
#'   \item{`ids`}{Character vector; the canonical unit-id order shared by every
#'     fit on the path.}
#'   \item{`lambda_grid`}{Numeric vector; the strictly increasing penalty grid
#'     that was solved, including the `0.25` and `1` anchors.}
#'   \item{`fits`}{List of `gp_grade_fit`, one per `lambda_grid` value in
#'     ascending order; the per-penalty grade assignments and solver status.}
#'   \item{`summary`}{Data frame, one row per penalty, with `lambda`,
#'     `grade_count`, `discordance_rate`, `reliability`, `tau_bar`, `objective`,
#'     `raw_objective`, and `status` -- the frontier table.}
#'   \item{`backend`}{Named list; the resolved solver name, paths, and encodings
#'     used across the solves.}
#'   \item{`selection`}{Named list (`selected_lambda`, `selection_rule`,
#'     `endpoint_lambda`); the recorded default selection, not a re-solve.}
#'   \item{`control`}{The validated [gp_control] used for the run.}
#'   \item{`warnings`, `schema_version`, `provenance`}{Internal audit slots
#'     (solve order, warm-start strategy, accumulated warnings).}
#' }
#'
#' @details
#' Grade labels are contiguous integers in `{1, ..., k}` and carry no
#' ranking-superiority statement of any kind. Solving descends from the
#' `lambda = 1` endpoint so each gurobi solve can warm-start from the adjacent
#' penalty; open backends solve each penalty independently. Selecting a fit later
#' (with [gp_select_grade()]) only reads the stored list and never re-solves. See
#' KRW (2024, Section 4) and the grade-engine method vignette for the integer
#' program and the frontier penalty.
#'
#' @examples
#' # Instant: read the solved path from the bundled tiny fit (24 firms).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' fit$grade_path                 # the pre-solved gp_grade_path
#' fit$grade_path$summary         # the frontier table over the penalty grid
#'
#' # Live solve from the stored pairwise structure; open HiGHS backend, small grid.
#' \donttest{
#' path <- gp_grade_path(
#'   fit$pairwise,
#'   control = gp_control(backend = "highs", lambda_grid = c(0.25, 1))
#' )
#' gp_select_grade(path, lambda = 0.25)
#' }
#'
#' @seealso [gp_grade()], [gp_select_grade()], [krw_report_card()], [gp_pairwise()]
#' @family gradepath-grade
#' @export
gp_grade_path <- function(pairwise,
                          lambda_grid = NULL,
                          control = NULL,
                          selected_lambda = 0.25,
                          selection_rule = "baseline_lambda_0.25",
                          acceptance_mode = FALSE) {
  pairwise <- validate_gp_pairwise(pairwise)
  acceptance_mode <- .gradepath_validate_scalar_logical(
    acceptance_mode,
    "acceptance_mode"
  )
  control <- .gp_grade_resolve_control(
    pairwise,
    control = control,
    lambda_grid = lambda_grid
  )
  lambda_grid <- control$lambda_grid
  selected <- .gp_grade_exact_lambda_match(
    lambda_grid,
    selected_lambda,
    "selected_lambda"
  )$value
  matrix_opt <- .gp_grade_matrix_for_optimization(pairwise)
  solve_order <- .gp_grade_solve_order(lambda_grid)
  fits <- vector("list", length(lambda_grid))
  warm_start <- NULL
  warm_start_lambda <- NA_real_
  warm_start_status <- NA_character_
  use_warm_start <- identical(control$backend, "gurobi")

  for (index in solve_order) {
    lambda <- lambda_grid[[index]]
    result <- .gp_grade_solve_backend(
      pairwise = pairwise,
      lambda = lambda,
      control = control,
      matrix_opt = matrix_opt,
      warm_start = if (isTRUE(use_warm_start)) warm_start else NULL,
      acceptance_mode = acceptance_mode
    )
    fits[[index]] <- .gp_grade_fit_from_result(
      result = result,
      ids = pairwise$ids,
      lambda = lambda,
      control = control,
      warm_start_from = if (isTRUE(use_warm_start)) warm_start_lambda else NA_real_,
      warm_start_from_status = if (isTRUE(use_warm_start)) {
        warm_start_status
      } else {
        NA_character_
      }
    )
    if (isTRUE(use_warm_start)) {
      warm_start <- result$primal
      warm_start_lambda <- lambda
      warm_start_status <- result$status
    }
  }

  new_gp_grade_path(
    ids = pairwise$ids,
    lambda_grid = lambda_grid,
    fits = fits,
    summary = frontier_table(pairwise, fits),
    backend = .gp_grade_path_backend(control, fits),
    selection = list(
      selected_lambda = selected,
      selection_rule = selection_rule,
      endpoint_lambda = 1
    ),
    control = control,
    provenance = .gradepath_new_provenance(
      producer = "gp_grade_path",
      solve_order = lambda_grid[solve_order],
      warm_start_strategy = if (isTRUE(use_warm_start)) {
        "lambda_1_then_descending"
      } else {
        "none"
      },
      warm_start_nonoptimal_policy = if (isTRUE(use_warm_start)) {
        "allowed_labeled"
      } else {
        "not_applicable"
      },
      warm_start_sources = .gp_grade_path_warm_start_sources(lambda_grid, fits)
    ),
    warnings = .gp_grade_path_warnings(fits)
  )
}

#' Select a stored grade fit from a solved grade path
#'
#' `gp_select_grade()` returns the pre-solved `gp_grade_fit` at a chosen penalty
#' from a [gp_grade_path]. It never runs a solver: it indexes into the path's
#' stored fits, so it is the instant way to read off one assignment after
#' [gp_grade_path()] has done the work. If `lambda` is absent from the solved
#' grid, it errors and reports the available grid values.
#'
#' @param path A [gp_grade_path] object; the solved frontier to select from.
#' @param lambda Numeric scalar; the penalty to select, which must be present in
#'   `path$lambda_grid`. Defaults to the KRW baseline `0.25`.
#'
#' @return The stored `gp_grade_fit` at `lambda` (a list of class
#'   `c("gp_grade_fit", "list")`), exactly as solved on the path -- no re-solve.
#'   See [gp_grade()] for the slot-level description of a `gp_grade_fit`.
#'
#' @examples
#' # Solve-free: pull the baseline-penalty fit from the bundled tiny fit's path.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' gp_select_grade(fit$grade_path, lambda = 0.25)
#'
#' # The pipeline already stored this selection; it is the same object.
#' identical(gp_select_grade(fit$grade_path, lambda = 0.25), fit$selected_grade)
#'
#' @seealso [gp_grade_path()], [gp_grade()], [krw_report_card()], [gp_pairwise()]
#' @family gradepath-grade
#' @export
gp_select_grade <- function(path, lambda = 0.25) {
  path <- validate_gp_grade_path(path)
  match <- .gp_grade_exact_lambda_match(
    path$lambda_grid,
    lambda,
    "lambda"
  )
  path$fits[[match$index]]
}

#' Solve a single-penalty grade fit
#'
#' `gp_grade()` solves the grade integer program at exactly one penalty and
#' returns a `gp_grade_fit`. It is a targeted one-shot wrapper for diagnostics and
#' small workflows; for the whole frontier use [gp_grade_path()], and to read a
#' penalty already on a solved path use [gp_select_grade()] (which never solves).
#'
#' @param pairwise A [gp_pairwise] object: the `J x J` posterior outranking
#'   structure the integer program scores.
#' @param lambda Numeric scalar; the penalty to solve, which must be present in
#'   the active `control` lambda grid. Defaults to the KRW baseline `0.25`.
#' @param control Optional [gp_control]; the backend, solver tolerances, and
#'   active lambda grid. Defaults to `pairwise$control`.
#' @param acceptance_mode Logical; the solver solution-quality policy, honored by
#'   the gurobi backend only (ignored for open backends). `FALSE` (default) keeps
#'   the first Gurobi path that returns a normalized result, so a solve that stops
#'   at `gap_reached` / `time_limit` is reported honestly with that status (it
#'   lands `UNVERIFIED` downstream). `TRUE` additionally attempts later Gurobi
#'   paths to prove the optimum when an earlier path is not acceptance-ready; it
#'   stays honest and never relabels a non-optimal solve as `optimal`. The chosen
#'   value is recorded on the fit's `backend$acceptance_mode`.
#'
#' @return A validated `gp_grade_fit` object (a list of class
#'   `c("gp_grade_fit", "list")`) with the public slots: \describe{
#'   \item{`ids`}{Character vector; the canonical unit-id order of the assignment.}
#'   \item{`lambda`}{Numeric scalar; the penalty this fit was solved at.}
#'   \item{`assignment`}{Data frame with one row per id: `id` (character) and
#'     `grade` (contiguous integer label in `{1, ..., k}`).}
#'   \item{`summary`}{Named list with `grade_count`, `status`, and `n_units` for
#'     the realized assignment.}
#'   \item{`objective`}{Named list; the penalized and raw objective values at this
#'     penalty.}
#'   \item{`backend`}{Named list; the solver `name`, normalized `status`,
#'     `mipgap`, `runtime`, `acceptance_mode`, and the problem hash/signature.}
#'   \item{`control`}{The validated [gp_control] used for the solve.}
#'   \item{`warnings`, `schema_version`, `provenance`}{Internal audit slots.}
#' }
#'
#' @details Grade labels are contiguous integers in `{1, ..., k}` and carry no
#' ranking-superiority statement of any kind. Solver honesty matches
#' [gp_grade_path()]: a gurobi solve that stops at a gap or time limit keeps that
#' status, and `acceptance_mode = TRUE` only adds optimization attempts, never a
#' relabel. See KRW (2024, Section 4) and the grade-engine method vignette.
#'
#' @examples
#' # Instant: the bundled tiny fit already stores the baseline-penalty grade fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' fit$selected_grade                 # the gp_grade_fit at lambda = 0.25
#' gp_select_grade(fit$grade_path, lambda = 0.25)  # the same object, no solve
#'
#' # Live single-penalty solve from the stored pairwise; open HiGHS backend.
#' \donttest{
#' gp_grade(
#'   fit$pairwise,
#'   lambda = 0.25,
#'   control = gp_control(backend = "highs", lambda_grid = c(0.25, 1))
#' )
#' }
#'
#' @seealso [gp_grade_path()], [gp_select_grade()], [krw_report_card()], [gp_pairwise()]
#' @family gradepath-grade
#' @export
gp_grade <- function(pairwise, lambda = 0.25, control = NULL,
                     acceptance_mode = FALSE) {
  pairwise <- validate_gp_pairwise(pairwise)
  acceptance_mode <- .gradepath_validate_scalar_logical(
    acceptance_mode,
    "acceptance_mode"
  )
  control <- if (is.null(control)) pairwise$control else validate_gp_control(control)
  match <- .gp_grade_exact_lambda_match(
    control$lambda_grid,
    lambda,
    "lambda"
  )
  matrix_opt <- .gp_grade_matrix_for_optimization(pairwise)
  result <- .gp_grade_solve_backend(
    pairwise = pairwise,
    lambda = match$value,
    control = control,
    matrix_opt = matrix_opt,
    warm_start = NULL,
    acceptance_mode = acceptance_mode
  )

  .gp_grade_fit_from_result(
    result = result,
    ids = pairwise$ids,
    lambda = match$value,
    control = control,
    warm_start_from = NA_real_
  )
}
