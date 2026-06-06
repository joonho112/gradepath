# Decision-side typed objects for gradepath v2: gp_pairwise, gp_grade_fit,
# gp_grade_path (Chapter 9 §schemas; "reuse" objects per the Ch9 disposition
# table).
#
# These three classes are gradepath's OWN social-choice spine. Per Chapter 9
# they port "nearly verbatim" from v1's `class-constructors-decision.R`, with
# exactly these v2 changes:
#
#   1. Naming: gradepath_pairwise/_grade_fit/_grade_path  ->  gp_pairwise /
#      gp_grade_fit / gp_grade_path; new_gradepath_* / validate_gradepath_*  ->
#      new_gp_* / validate_gp_* (Chapter 9, GP-DEC-09-A; the constructors/
#      validators stay INTERNAL -- the Ch13 export surface is the four verbs +
#      monolith + gp_control + get_* accessors, not these new_/validate_ pairs).
#   2. `control` slot is re-typed to `gp_control`: validated by
#      `validate_gp_control()` (control.R), NOT v1's `validate_gradepath_control`.
#   3. The `replication_mode` branch v1 embedded in `validate_gradepath_grade_path`
#      (forcing selected_lambda = 0.25 / endpoint_lambda = 1.00 / a control-
#      coupled selection_rule) is DELETED. Per Chapter 9 §grade-path ("Selection
#      is first-class, not a control lock") the validator now only asserts grid
#      membership; the baseline value 0.25 is a default of gp_select_grade()
#      (Ch13), not a validator invariant. (gp_control has no `replication_mode`
#      field at all, so the branch could not run regardless.)
#
# Everything else is ported as-is: the slot field vectors, the canonical order
# enforced by `.gradepath_validate_named_fields()`, the schema_version tag
# (`.gradepath_schema_version`, "v2"), and the provenance-as-named-list contract.
#
# The W-seam shape assertion (Chapter 8 W1-W5; frozen invariant 5) lives at the
# gp_pairwise boundary. v1's `gp_pairwise` schema stores no `W` slot -- it stores
# the cleaned J x J `matrix` -- so the orientation hazard is caught HERE by the
# square + dimname + antisymmetry/diagonal/zero-floor checks: a mis-oriented or
# mis-normalized W fed through the pi_ij builder produces a matrix that fails one
# of those invariants with a gp_pairwise-named error before it can reach the
# solver (Chapter 8 callout "The W-seam shape assertion lives here"; Chapter 9
# §pairwise callout). See the note at validate_gp_pairwise() and SUMMARY.md for
# the full W-orientation discussion.
#
# Reuses (do NOT redefine -- built in Steps 1.2/1.3):
#   utils-validate.R : .gradepath_validate_* family, .gradepath_validate_named_fields,
#                      .gradepath_new_provenance, .gradepath_abort*,
#                      .gradepath_schema_version, gp_diag_positions
#   control.R        : gp_control, validate_gp_control, .gp_control_validate_lambda_grid

# ---------------------------------------------------------------------------
# gp_pairwise  -- the posterior pairwise probability matrix (Ch9 §pairwise)
# ---------------------------------------------------------------------------
#
# Ch9 slot schema (9 slots, canonical order; ends schema_version/provenance/
# warnings):
#   ids            chr[J]      canonical unit ids, unique
#   matrix         dbl[JxJ]    pi_ij in [0,1]; square; row/colnames = ids;
#                              diag = 0.5; antisymmetric (pi_ij + pi_ji = 1) to
#                              tolerance; no exact off-diagonal zeros
#   power          int         matrix power applied (public path: 0L)
#   cleanup        named list  antisymmetry = TRUE, diagonal = 0.5, zero_floor = 1e-7
#   source         named list  stage = "posterior", rule in {outer_product, ...},
#                              assumption in {one_level_independence,
#                                             grouped_industry_dependence}
#   control        gp_control  the run settings
#   schema_version chr         spine schema tag
#   provenance     named list  build stamp
#   warnings       chr[]       accumulated non-fatal notes

.gp_pairwise_fields <- c(
  "ids",
  "matrix",
  "power",
  "cleanup",
  "source",
  "control",
  "schema_version",
  "provenance",
  "warnings"
)

.gp_pairwise_cleanup_fields <- c("antisymmetry", "diagonal", "zero_floor")
.gp_pairwise_source_fields <- c("stage", "rule", "assumption")

# Cleanup invariant anchors (Chapter 4 §cleanup; frozen invariant 6). The 0.5
# diagonal and the 1e-7 zero-floor are `.frozen` bit-exact parity anchors -- do
# not relax.
.gp_pairwise_diagonal <- 0.5
.gp_pairwise_zero_floor <- 1e-7

# Accepted `source$rule` / `source$assumption` vocabularies (ported from v1;
# Ch9 lists rule in {outer_product, ...}, assumption in {one_level_independence,
# grouped_industry_dependence}). v1's second rule name is kept as an accepted
# alternative for the two-level (group_fx==1) archive matrix path.
.gp_pairwise_source_rules <- c(
  "outer_product",
  "groupfx1_archive_matrix"
)
.gp_pairwise_source_assumptions <- c(
  "one_level_independence",
  "grouped_industry_dependence"
)

# Low-level constructor: assemble in canonical order with the v2 schema tag, then
# validate. Internal; the public pi_ij verb (gp_pairwise(), Ch4/Ch13) calls this.
new_gp_pairwise <- function(ids,
                            matrix,
                            power,
                            cleanup,
                            source,
                            control,
                            schema_version = .gradepath_schema_version,
                            provenance = list(),
                            warnings = character()) {
  x <- structure(
    list(
      ids = ids,
      matrix = matrix,
      power = power,
      cleanup = cleanup,
      source = source,
      control = control,
      schema_version = schema_version,
      provenance = provenance,
      warnings = warnings
    ),
    class = c("gp_pairwise", "list")
  )

  validate_gp_pairwise(x)
}

# Validator. Asserts class + exact canonical field order, then every slot's
# type/contract. Crucially pins the Chapter 4 cleanup contract / frozen invariant
# 6 on the matrix itself:
#   * square J x J, entries in [0, 1], row/colnames = ids when present;
#   * diagonal == 0.5 (invariant 6) to tol 1e-8;
#   * antisymmetry pi_ij + pi_ji == 1 off-diagonal, to tolerance;
#   * zero-floor: no exact off-diagonal zeros (the 1e-7 one-sided floor; exact
#     ones are deliberately ALLOWED -- Ch4 §cleanup "Why the cap is one-sided").
#
# THE W-SEAM ASSERTION (Chapter 8, frozen invariant 5). The gp_pairwise schema
# does not carry W (it carries the cleaned J x J pi matrix), so there is no
# rowSums(W) check here. Instead the orientation hazard is caught structurally:
# the square + dimname-consistency + antisymmetry/diagonal/zero-floor invariants
# below are exactly the properties a correctly J x M row-stochastic W produces
# after the pi_ij build, and that a transposed (M x J, column-stochastic) W would
# violate -- so a mis-oriented W fails HERE with a gp_pairwise-named error before
# the solver. If a future schema revision adds an explicit W slot (or a
# W-orientation metadata flag), assert W is J x M with all(abs(rowSums(W)-1) <
# 1e-12) and reject the v1 transpose at that point (see the helper
# `.gp_assert_w_row_stochastic()` provided below for that eventuality).
validate_gp_pairwise <- function(x) {
  .gradepath_validate_list_class(x, "gp_pairwise")
  .gradepath_validate_named_fields(x, .gp_pairwise_fields, "gp_pairwise")

  ids <- .gradepath_validate_character_vector(
    x$ids,
    "gp_pairwise$ids",
    unique = TRUE
  )
  matrix <- .gradepath_validate_numeric_matrix(
    x$matrix,
    "gp_pairwise$matrix"
  )

  # --- square J x J + id/dimname consistency (the orientation guard) ---------
  if (!identical(dim(matrix), c(length(ids), length(ids)))) {
    .gradepath_abort(
      "`gp_pairwise$matrix` must be square with one row and column per `gp_pairwise$ids` value."
    )
  }

  if (any(matrix < 0) || any(matrix > 1)) {
    .gradepath_abort("`gp_pairwise$matrix` entries must lie in [0, 1].")
  }

  row_ids <- rownames(matrix)
  col_ids <- colnames(matrix)

  if ((!is.null(row_ids) && !identical(row_ids, ids)) ||
      (!is.null(col_ids) && !identical(col_ids, ids))) {
    .gradepath_abort(
      "When present, the row and column names of `gp_pairwise$matrix` must exactly match `gp_pairwise$ids`."
    )
  }

  power <- .gradepath_validate_integerish(
    x$power,
    "gp_pairwise$power",
    min = 0L
  )

  if (!identical(power, 0L)) {
    .gradepath_abort("`gp_pairwise$power` must be 0 for the public path.")
  }

  # --- cleanup invariant block (frozen: antisymmetry/diag 0.5/zero-floor) ----
  cleanup <- .gradepath_validate_named_list(
    x$cleanup,
    "gp_pairwise$cleanup"
  )
  .gradepath_validate_named_fields(
    cleanup,
    .gp_pairwise_cleanup_fields,
    "gp_pairwise$cleanup"
  )
  cleanup$antisymmetry <- .gradepath_validate_scalar_logical(
    cleanup$antisymmetry,
    "gp_pairwise$cleanup$antisymmetry"
  )
  cleanup$diagonal <- .gradepath_validate_scalar_numeric(
    cleanup$diagonal,
    "gp_pairwise$cleanup$diagonal",
    lower = 0,
    upper = 1,
    include_lower = TRUE,
    include_upper = TRUE
  )
  cleanup$zero_floor <- .gradepath_validate_scalar_numeric(
    cleanup$zero_floor,
    "gp_pairwise$cleanup$zero_floor",
    lower = 0,
    upper = 1,
    include_lower = TRUE,
    include_upper = FALSE
  )

  if (!isTRUE(cleanup$antisymmetry)) {
    .gradepath_abort("`gp_pairwise$cleanup$antisymmetry` must be TRUE.")
  }

  # Frozen anchor: diagonal flag must be exactly 0.5 (invariant 6).
  if (abs(cleanup$diagonal - .gp_pairwise_diagonal) > 1e-8) {
    .gradepath_abort("`gp_pairwise$cleanup$diagonal` must be 0.5.")
  }

  # Frozen anchor: zero-floor flag must be exactly 1e-7.
  if (abs(cleanup$zero_floor - .gp_pairwise_zero_floor) > 1e-12) {
    .gradepath_abort("`gp_pairwise$cleanup$zero_floor` must be 1e-7.")
  }

  # Matrix diagonal must equal the declared 0.5 diagonal (invariant 6 on the
  # data, not just the flag).
  if (any(abs(diag(matrix) - cleanup$diagonal) > 1e-8)) {
    .gradepath_abort(
      "The diagonal of `gp_pairwise$matrix` must equal `gp_pairwise$cleanup$diagonal` (0.5)."
    )
  }

  off_diag <- row(matrix) != col(matrix)

  # Zero-floor: no EXACT off-diagonal zeros after cleanup (the 1e-7 floor removed
  # them). One-sided by design -- exact ones remain allowed (Ch4 §cleanup).
  if (any(matrix[off_diag] == 0)) {
    .gradepath_abort(
      "`gp_pairwise$matrix` must not contain exact off-diagonal zeros after cleanup."
    )
  }

  # Antisymmetry: pi_ij + pi_ji == 1 off-diagonal, within the cleanup tolerance
  # (zero_floor + 1e-8, ~1e-7; the floor can perturb a partner by up to the floor
  # value, so the band is widened by it -- ported from v1).
  antisymmetry_gap <- abs(matrix + t(matrix) - 1)
  tolerance <- cleanup$zero_floor + 1e-8

  if (any(antisymmetry_gap[off_diag] > tolerance)) {
    .gradepath_abort(
      "`gp_pairwise$matrix` must satisfy the declared antisymmetry rule (pi_ij + pi_ji = 1) within the cleanup tolerance."
    )
  }

  # --- source provenance block ----------------------------------------------
  source <- .gradepath_validate_named_list(
    x$source,
    "gp_pairwise$source"
  )
  .gradepath_validate_named_fields(
    source,
    .gp_pairwise_source_fields,
    "gp_pairwise$source"
  )
  source$stage <- .gradepath_validate_scalar_character(
    source$stage,
    "gp_pairwise$source$stage",
    allowed = "posterior"
  )
  source$rule <- .gradepath_validate_scalar_character(
    source$rule,
    "gp_pairwise$source$rule",
    allowed = .gp_pairwise_source_rules
  )
  source$assumption <- .gradepath_validate_scalar_character(
    source$assumption,
    "gp_pairwise$source$assumption",
    allowed = .gp_pairwise_source_assumptions
  )

  # --- control (re-typed to gp_control) + audit slots -----------------------
  control <- validate_gp_control(x$control)
  schema_version <- .gradepath_validate_scalar_character(
    x$schema_version,
    "gp_pairwise$schema_version",
    allowed = .gradepath_schema_version
  )
  provenance <- .gradepath_validate_named_list(
    x$provenance,
    "gp_pairwise$provenance"
  )
  warnings <- .gradepath_validate_warning_vector(
    x$warnings,
    "gp_pairwise$warnings"
  )

  structure(
    list(
      ids = ids,
      matrix = matrix,
      power = power,
      cleanup = cleanup,
      source = source,
      control = control,
      schema_version = schema_version,
      provenance = provenance,
      warnings = warnings
    ),
    class = c("gp_pairwise", "list")
  )
}

# Standalone W-orientation guard (frozen invariant 5; Chapter 8 W1/W2). The
# gp_pairwise SCHEMA does not store W, so this is NOT called by
# validate_gp_pairwise above; it is the assertion the W-recompute seam
# (gp_posterior_weights()) and any future W-bearing schema must
# use to reject the v1 transpose (M x J, column-stochastic). Provided here so the
# orientation contract has a single named home at the decision boundary.
#
#   J = length(ids) units (rows), M = length(prior$support) grid points (cols).
#   W1 orientation : dim(W) == c(J, M)
#   W2 row-stochastic: all(abs(rowSums(W) - 1) < tol)   tol = 1e-12
# A transposed W (M x J) trips W1; a column-normalized W trips W2.
.gp_assert_w_row_stochastic <- function(W, n_units, n_support, name = "W", tol = 1e-12) {
  if (!is.matrix(W) || !is.numeric(W)) {
    .gradepath_abort("`%s` must be a numeric matrix.", name)
  }
  if (!identical(dim(W), c(as.integer(n_units), as.integer(n_support)))) {
    .gradepath_abort(
      "`%s` must be J x M (units x support): expected c(%d, %d), got c(%d, %d). (Reject the v1 M x J transpose -- frozen invariant 5.)",
      name, as.integer(n_units), as.integer(n_support), nrow(W), ncol(W)
    )
  }
  if (any(!is.finite(W)) || any(W < 0)) {
    .gradepath_abort("`%s` must be non-negative and finite (Chapter 8 W5).", name)
  }
  if (any(abs(rowSums(W) - 1) > tol)) {
    .gradepath_abort(
      "`%s` must be row-stochastic: each row (a unit's posterior mass) must sum to 1 (Chapter 8 W2). A column-normalized / transposed W is rejected.",
      name
    )
  }
  invisible(W)
}

# ---------------------------------------------------------------------------
# gp_grade_fit  -- one solved lambda (Ch9 §grade-fit)
# ---------------------------------------------------------------------------
#
# Ch9 slot schema (10 slots; ends schema_version/provenance/warnings):
#   ids            chr[J]            unique unit ids
#   lambda         dbl               loss-tradeoff in [0, 1]
#   assignment     data.frame[J]     columns id (= ids) + grade (contiguous int
#                                    1,2,..., normalized)
#   summary        named list        at least grade_count (int) = #distinct grades
#   objective      named list        at least value (dbl) = IP objective
#   backend        named list        at least name (chr) = the solver
#   control        gp_control        run settings
#   schema_version chr  / provenance named list / warnings chr[]   audit slots

.gp_grade_fit_fields <- c(
  "ids",
  "lambda",
  "assignment",
  "summary",
  "objective",
  "backend",
  "control",
  "schema_version",
  "provenance",
  "warnings"
)

new_gp_grade_fit <- function(ids,
                             lambda,
                             assignment,
                             summary,
                             objective,
                             backend,
                             control,
                             schema_version = .gradepath_schema_version,
                             provenance = list(),
                             warnings = character()) {
  x <- structure(
    list(
      ids = ids,
      lambda = lambda,
      assignment = assignment,
      summary = summary,
      objective = objective,
      backend = backend,
      control = control,
      schema_version = schema_version,
      provenance = provenance,
      warnings = warnings
    ),
    class = c("gp_grade_fit", "list")
  )

  validate_gp_grade_fit(x)
}

# Validator. Pins: integer grade labels (normalized to contiguous 1..k), ids
# length-match (assignment one row per id, assignment$id == ids), lambda in
# [0, 1], control = gp_control, and the audit slots. Structurally guarantees
# summary$grade_count == realized distinct-grade count (so the headline grade
# count K08-K11 cannot disagree with the assignment it summarizes -- Ch9).
validate_gp_grade_fit <- function(x) {
  .gradepath_validate_list_class(x, "gp_grade_fit")
  .gradepath_validate_named_fields(x, .gp_grade_fit_fields, "gp_grade_fit")

  ids <- .gradepath_validate_character_vector(
    x$ids,
    "gp_grade_fit$ids",
    unique = TRUE
  )
  lambda <- .gradepath_validate_scalar_numeric(
    x$lambda,
    "gp_grade_fit$lambda",
    lower = 0,
    upper = 1,
    include_lower = TRUE,
    include_upper = TRUE
  )
  assignment <- .gradepath_validate_data_frame(
    x$assignment,
    "gp_grade_fit$assignment"
  )

  if (nrow(assignment) != length(ids)) {
    .gradepath_abort(
      "`gp_grade_fit$assignment` must have one row per id."
    )
  }

  required_assignment_columns <- c("id", "grade")
  missing_assignment_columns <- setdiff(required_assignment_columns, names(assignment))

  if (length(missing_assignment_columns) > 0L) {
    .gradepath_abort(
      "`gp_grade_fit$assignment` must include required columns: %s.",
      paste(required_assignment_columns, collapse = ", ")
    )
  }

  assignment$id <- .gradepath_validate_character_vector(
    assignment$id,
    "gp_grade_fit$assignment$id",
    unique = TRUE
  )

  if (!identical(assignment$id, ids)) {
    .gradepath_abort(
      "`gp_grade_fit$assignment$id` must exactly match `gp_grade_fit$ids`."
    )
  }

  assignment$grade <- .gradepath_validate_numeric_vector(
    assignment$grade,
    "gp_grade_fit$assignment$grade"
  )

  # Integer grade labels (Ch9: contiguous integers 1,2,...). Reject non-integer
  # or < 1 before normalizing.
  if (any(abs(assignment$grade - round(assignment$grade)) > 1e-8) ||
      any(assignment$grade < 1)) {
    .gradepath_abort(
      "`gp_grade_fit$assignment$grade` must be positive contiguous integers."
    )
  }

  assignment$grade <- as.integer(round(assignment$grade))
  unique_grades <- sort(unique(assignment$grade))

  # Contiguity: the realized label set must be exactly {1, ..., k}.
  if (!identical(unique_grades, seq_len(length(unique_grades)))) {
    .gradepath_abort(
      "`gp_grade_fit$assignment$grade` must be normalized to contiguous integers starting at 1."
    )
  }

  summary <- .gradepath_validate_named_list(
    x$summary,
    "gp_grade_fit$summary"
  )
  .gradepath_validate_required_keys(
    summary,
    "grade_count",
    "gp_grade_fit$summary"
  )
  summary$grade_count <- .gradepath_validate_integerish(
    summary$grade_count,
    "gp_grade_fit$summary$grade_count",
    min = 1L
  )

  # The structural grade-count guarantee.
  if (!identical(summary$grade_count, as.integer(length(unique_grades)))) {
    .gradepath_abort(
      "`gp_grade_fit$summary$grade_count` must equal the number of contiguous grades in `assignment`."
    )
  }

  objective <- .gradepath_validate_named_list(
    x$objective,
    "gp_grade_fit$objective"
  )
  .gradepath_validate_required_keys(
    objective,
    "value",
    "gp_grade_fit$objective"
  )
  objective$value <- .gradepath_validate_scalar_numeric(
    objective$value,
    "gp_grade_fit$objective$value"
  )

  backend <- .gradepath_validate_named_list(
    x$backend,
    "gp_grade_fit$backend"
  )
  .gradepath_validate_required_keys(
    backend,
    "name",
    "gp_grade_fit$backend"
  )
  backend$name <- .gradepath_validate_scalar_character(
    backend$name,
    "gp_grade_fit$backend$name"
  )
  if ("warm_start_from_lambda" %in% names(backend)) {
    if (!is.numeric(backend$warm_start_from_lambda) ||
        length(backend$warm_start_from_lambda) != 1L) {
      .gradepath_abort(
        "`gp_grade_fit$backend$warm_start_from_lambda` must be a length-1 numeric value or NA."
      )
    }
    if (!is.na(backend$warm_start_from_lambda)) {
      backend$warm_start_from_lambda <- .gradepath_validate_scalar_numeric(
        backend$warm_start_from_lambda,
        "gp_grade_fit$backend$warm_start_from_lambda",
        lower = 0,
        upper = 1,
        include_lower = TRUE,
        include_upper = TRUE
      )
    }
  }
  if ("warm_start_from_status" %in% names(backend)) {
    if (!is.character(backend$warm_start_from_status) ||
        length(backend$warm_start_from_status) != 1L ||
        (!is.na(backend$warm_start_from_status) &&
          !nzchar(backend$warm_start_from_status))) {
      .gradepath_abort(
        "`gp_grade_fit$backend$warm_start_from_status` must be a length-1 character value or NA."
      )
    }
  }
  if ("warm_start_from_acceptance_ready" %in% names(backend)) {
    if (!is.logical(backend$warm_start_from_acceptance_ready) ||
        length(backend$warm_start_from_acceptance_ready) != 1L) {
      .gradepath_abort(
        "`gp_grade_fit$backend$warm_start_from_acceptance_ready` must be a length-1 logical value or NA."
      )
    }
  }
  if ("warm_start_used" %in% names(backend)) {
    if (!is.logical(backend$warm_start_used) ||
        length(backend$warm_start_used) != 1L ||
        is.na(backend$warm_start_used)) {
      .gradepath_abort(
        "`gp_grade_fit$backend$warm_start_used` must be TRUE or FALSE."
      )
    }
  }

  control <- validate_gp_control(x$control)
  schema_version <- .gradepath_validate_scalar_character(
    x$schema_version,
    "gp_grade_fit$schema_version",
    allowed = .gradepath_schema_version
  )
  provenance <- .gradepath_validate_named_list(
    x$provenance,
    "gp_grade_fit$provenance"
  )
  warnings <- .gradepath_validate_warning_vector(
    x$warnings,
    "gp_grade_fit$warnings"
  )

  structure(
    list(
      ids = ids,
      lambda = lambda,
      assignment = assignment,
      summary = summary,
      objective = objective,
      backend = backend,
      control = control,
      schema_version = schema_version,
      provenance = provenance,
      warnings = warnings
    ),
    class = c("gp_grade_fit", "list")
  )
}

# ---------------------------------------------------------------------------
# gp_grade_path  -- the lambda sweep + selection record (Ch9 §grade-path)
# ---------------------------------------------------------------------------
#
# Ch9 slot schema (10 slots; ends schema_version/provenance/warnings):
#   ids            chr[J]                  unique unit ids
#   lambda_grid    dbl[L]                  strictly increasing, unique, in [0, 1]
#   fits           list[L] of gp_grade_fit one validated fit per lambda; each
#                                          fit's ids == this ids; each fit's
#                                          lambda aligns to lambda_grid[i]
#   summary        data.frame[L]           columns lambda, grade_count; aligned
#                                          to lambda_grid and to stored fits'
#                                          counts
#   backend        named list              name (chr)
#   selection      named list              selected_lambda, selection_rule,
#                                          endpoint_lambda; selected & endpoint
#                                          must be members of lambda_grid
#   control        gp_control / schema_version chr / provenance named list /
#   warnings       chr[]                   run + audit slots

.gp_grade_path_fields <- c(
  "ids",
  "lambda_grid",
  "fits",
  "summary",
  "backend",
  "selection",
  "control",
  "schema_version",
  "provenance",
  "warnings"
)

.gp_grade_path_selection_fields <- c(
  "selected_lambda",
  "selection_rule",
  "endpoint_lambda"
)

new_gp_grade_path <- function(ids,
                              lambda_grid,
                              fits,
                              summary,
                              backend,
                              selection,
                              control,
                              schema_version = .gradepath_schema_version,
                              provenance = list(),
                              warnings = character()) {
  x <- structure(
    list(
      ids = ids,
      lambda_grid = lambda_grid,
      fits = fits,
      summary = summary,
      backend = backend,
      selection = selection,
      control = control,
      schema_version = schema_version,
      provenance = provenance,
      warnings = warnings
    ),
    class = c("gp_grade_path", "list")
  )

  validate_gp_grade_path(x)
}

# Validator. Pins:
#   * lambda_grid strictly increasing / unique / in [0, 1] (via the shared
#     control helper .gp_control_validate_lambda_grid -- which ADDITIONALLY
#     requires the parity anchors {0.25, 1.00}; see SUMMARY.md "Ch9 ambiguities");
#   * fits is a length-L list of valid gp_grade_fit, each fit$ids == path ids,
#     each fit$lambda aligned to lambda_grid[i];
#   * summary is an L-row df with columns lambda, grade_count, aligned to the
#     grid AND to the stored fits' counts;
#   * selection (selected_lambda / selection_rule / endpoint_lambda) with
#     selected_lambda and endpoint_lambda BOTH members of lambda_grid.
#
# v2 CHANGE: the v1 `replication_mode` branch (forcing selected_lambda = 0.25 /
# endpoint_lambda = 1.00 / selection_rule = control$selection_rule) is DELETED.
# Selection is first-class (Ch9 §grade-path); the validator asserts only grid
# membership.
validate_gp_grade_path <- function(x) {
  .gradepath_validate_list_class(x, "gp_grade_path")
  .gradepath_validate_named_fields(x, .gp_grade_path_fields, "gp_grade_path")

  ids <- .gradepath_validate_character_vector(
    x$ids,
    "gp_grade_path$ids",
    unique = TRUE
  )
  # Strictly-increasing / unique / in [0,1] -- the anchor-FREE core check
  # (Chapter 9 §grade-path: the solved grid need not contain the parity anchors;
  # a user may solve a sub-grid). The {0.25, 1.00} anchor requirement lives only
  # on gp_control's reference grid, not on a solved path's grid.
  lambda_grid <- .gp_validate_lambda_grid_core(
    x$lambda_grid,
    "gp_grade_path$lambda_grid"
  )

  if (!is.list(x$fits) || length(x$fits) != length(lambda_grid)) {
    .gradepath_abort(
      "`gp_grade_path$fits` must be a list with one fit per lambda value."
    )
  }

  fits <- lapply(seq_along(x$fits), function(i) {
    fit_i <- validate_gp_grade_fit(x$fits[[i]])

    if (!identical(fit_i$ids, ids)) {
      .gradepath_abort(
        "Each stored fit in `gp_grade_path$fits` must share the same ids as `gp_grade_path$ids`."
      )
    }

    if (abs(fit_i$lambda - lambda_grid[[i]]) > 1e-8) {
      .gradepath_abort(
        "Each stored fit in `gp_grade_path$fits` must align to `gp_grade_path$lambda_grid` in order."
      )
    }

    fit_i
  })

  summary <- .gradepath_validate_data_frame(
    x$summary,
    "gp_grade_path$summary"
  )

  if (nrow(summary) != length(lambda_grid)) {
    .gradepath_abort(
      "`gp_grade_path$summary` must have one row per lambda value."
    )
  }

  required_summary_columns <- c("lambda", "grade_count")
  missing_summary_columns <- setdiff(required_summary_columns, names(summary))

  if (length(missing_summary_columns) > 0L) {
    .gradepath_abort(
      "`gp_grade_path$summary` must include required columns: %s.",
      paste(required_summary_columns, collapse = ", ")
    )
  }

  summary$lambda <- .gradepath_validate_numeric_vector(
    summary$lambda,
    "gp_grade_path$summary$lambda"
  )
  summary$grade_count <- .gradepath_validate_numeric_vector(
    summary$grade_count,
    "gp_grade_path$summary$grade_count"
  )

  if (any(abs(summary$grade_count - round(summary$grade_count)) > 1e-8) ||
      any(summary$grade_count < 1)) {
    .gradepath_abort(
      "`gp_grade_path$summary$grade_count` must be positive integers."
    )
  }

  summary$grade_count <- as.integer(round(summary$grade_count))

  # summary$lambda aligned exactly to the grid.
  if (!isTRUE(all.equal(summary$lambda, lambda_grid, tolerance = 1e-8))) {
    .gradepath_abort(
      "`gp_grade_path$summary$lambda` must align exactly to `gp_grade_path$lambda_grid`."
    )
  }

  # summary$grade_count aligned to the stored fits' counts.
  fit_grade_counts <- vapply(
    fits,
    function(fit_i) fit_i$summary$grade_count,
    integer(1)
  )

  if (!identical(summary$grade_count, fit_grade_counts)) {
    .gradepath_abort(
      "`gp_grade_path$summary$grade_count` must align to the stored fits."
    )
  }

  backend <- .gradepath_validate_named_list(
    x$backend,
    "gp_grade_path$backend"
  )
  .gradepath_validate_required_keys(
    backend,
    "name",
    "gp_grade_path$backend"
  )
  backend$name <- .gradepath_validate_scalar_character(
    backend$name,
    "gp_grade_path$backend$name"
  )

  selection <- .gradepath_validate_named_list(
    x$selection,
    "gp_grade_path$selection"
  )
  .gradepath_validate_named_fields(
    selection,
    .gp_grade_path_selection_fields,
    "gp_grade_path$selection"
  )
  selection$selected_lambda <- .gradepath_validate_scalar_numeric(
    selection$selected_lambda,
    "gp_grade_path$selection$selected_lambda",
    lower = 0,
    upper = 1,
    include_lower = TRUE,
    include_upper = TRUE
  )
  selection$selection_rule <- .gradepath_validate_scalar_character(
    selection$selection_rule,
    "gp_grade_path$selection$selection_rule"
  )
  selection$endpoint_lambda <- .gradepath_validate_scalar_numeric(
    selection$endpoint_lambda,
    "gp_grade_path$selection$endpoint_lambda",
    lower = 0,
    upper = 1,
    include_lower = TRUE,
    include_upper = TRUE
  )

  # Grid membership -- the ONLY selection assertion v2 makes. Require exactly
  # one tolerance match so near-duplicate reconstructed grids cannot make a
  # selection ambiguous.
  selected_hits <- which(abs(lambda_grid - selection$selected_lambda) < 1e-8)
  endpoint_hits <- which(abs(lambda_grid - selection$endpoint_lambda) < 1e-8)

  if (length(selected_hits) != 1L || length(endpoint_hits) != 1L) {
    .gradepath_abort(
      "Both `selected_lambda` and `endpoint_lambda` must match exactly one member of `gp_grade_path$lambda_grid`."
    )
  }

  # NOTE (v2): v1's replication_mode branch is intentionally absent here.

  control <- validate_gp_control(x$control)
  schema_version <- .gradepath_validate_scalar_character(
    x$schema_version,
    "gp_grade_path$schema_version",
    allowed = .gradepath_schema_version
  )
  provenance <- .gradepath_validate_named_list(
    x$provenance,
    "gp_grade_path$provenance"
  )
  if ("warm_start_sources" %in% names(provenance)) {
    warm_start_sources <- .gradepath_validate_data_frame(
      provenance$warm_start_sources,
      "gp_grade_path$provenance$warm_start_sources"
    )
    required_warm_start_columns <- c(
      "lambda",
      "warm_start_from_lambda",
      "warm_start_from_status",
      "warm_start_from_acceptance_ready",
      "warm_start_used"
    )
    missing_warm_start_columns <- setdiff(
      required_warm_start_columns,
      names(warm_start_sources)
    )
    if (length(missing_warm_start_columns) > 0L) {
      .gradepath_abort(
        "`gp_grade_path$provenance$warm_start_sources` must include required columns: %s.",
        paste(required_warm_start_columns, collapse = ", ")
      )
    }
    if (nrow(warm_start_sources) != length(lambda_grid)) {
      .gradepath_abort(
        "`gp_grade_path$provenance$warm_start_sources` must have one row per lambda value."
      )
    }
    if (!is.numeric(warm_start_sources$lambda) ||
        !isTRUE(all.equal(warm_start_sources$lambda, lambda_grid, tolerance = 1e-8))) {
      .gradepath_abort(
        "`gp_grade_path$provenance$warm_start_sources$lambda` must align exactly to `gp_grade_path$lambda_grid`."
      )
    }
    if (!is.numeric(warm_start_sources$warm_start_from_lambda) ||
        any(!is.na(warm_start_sources$warm_start_from_lambda) &
          (warm_start_sources$warm_start_from_lambda < 0 |
            warm_start_sources$warm_start_from_lambda > 1))) {
      .gradepath_abort(
        "`gp_grade_path$provenance$warm_start_sources$warm_start_from_lambda` must contain lambda values or NA."
      )
    }
    if (!is.character(warm_start_sources$warm_start_from_status) ||
        any(!is.na(warm_start_sources$warm_start_from_status) &
          !nzchar(warm_start_sources$warm_start_from_status))) {
      .gradepath_abort(
        "`gp_grade_path$provenance$warm_start_sources$warm_start_from_status` must contain status labels or NA."
      )
    }
    if (!is.logical(warm_start_sources$warm_start_from_acceptance_ready)) {
      .gradepath_abort(
        "`gp_grade_path$provenance$warm_start_sources$warm_start_from_acceptance_ready` must be logical."
      )
    }
    if (!is.logical(warm_start_sources$warm_start_used) ||
        anyNA(warm_start_sources$warm_start_used)) {
      .gradepath_abort(
        "`gp_grade_path$provenance$warm_start_sources$warm_start_used` must contain TRUE/FALSE values."
      )
    }
    provenance$warm_start_sources <- warm_start_sources
  }
  warnings <- .gradepath_validate_warning_vector(
    x$warnings,
    "gp_grade_path$warnings"
  )

  structure(
    list(
      ids = ids,
      lambda_grid = lambda_grid,
      fits = fits,
      summary = summary,
      backend = backend,
      selection = selection,
      control = control,
      schema_version = schema_version,
      provenance = provenance,
      warnings = warnings
    ),
    class = c("gp_grade_path", "list")
  )
}
