# Validation, provenance, and small structural utilities for gradepath v2.
#
# This file ports the v1 `utils-validate.R` validator family (the design marks
# `.gradepath_validate_named_*` and `.gradepath_new_provenance` as **reuse**),
# retypes the schema tag to v2, and adds the structural helpers the
# downstream grade-IP and harness stages consume: `gp_diag_positions()`,
# `monotone_relabel()` (harness; also the contiguous-1..k
# normalizer the typed objects assert), and the stable row
# log-sum-exp `.gp_row_logsumexp()`.
#
# The internal namespace stays `.gradepath_*` (NOT `.gp_*`) for the validator
# family, matching the v1 surface the typed-object Steps validate against
# (Chapter 9 §schemas: "v1's `.gradepath_validate_named_fields` is reused
# verbatim"). The two leaner numeric helpers reused widely (`.gp_row_logsumexp`,
# `.gp_safe_normalize`) carry the shorter `.gp_*` prefix to match their math
# siblings in utils-math.R.

# ---------------------------------------------------------------------------
# Schema tag
# ---------------------------------------------------------------------------

# Single source of truth for the typed-object schema version (Chapter 9). Bumped
# to "v2" for the gradepath v2 spine; every `gp_*` object stamps this in its
# `schema_version` slot.
.gradepath_schema_version <- "v2"

# ---------------------------------------------------------------------------
# Condition / abort helpers  (Chapter 13 error-class contract)
# ---------------------------------------------------------------------------

# Base abort: builds and signals a classed condition. `message` is an sprintf
# template; `...` are its arguments. `class` is prepended to the standard
# c("error", "condition") tail so callers can `tryCatch()` on a gradepath class.
.gradepath_abort <- function(message, ..., class = "gradepath_error", call = FALSE) {
  if (!is.character(message) || length(message) != 1L || is.na(message) || !nzchar(message)) {
    stop("`message` must be a length-1 non-empty character string.", call. = FALSE)
  }

  if (!is.character(class) || length(class) < 1L || anyNA(class) || any(!nzchar(class))) {
    stop("`class` must be a non-empty character vector.", call. = FALSE)
  }

  condition <- structure(
    list(
      message = sprintf(message, ...),
      call = if (isTRUE(call)) sys.call(-1L) else NULL
    ),
    class = c(class, "error", "condition")
  )

  stop(condition)
}

# Internal-invariant violations (bugs, not user error).
.gradepath_abort_internal <- function(message, ..., call = FALSE) {
  .gradepath_abort(
    message,
    ...,
    class = c("gradepath_internal_error", "gradepath_error"),
    call = call
  )
}

# Not-yet-implemented surface.
.gradepath_not_implemented <- function(feature) {
  if (!is.character(feature) || length(feature) != 1L || is.na(feature) || !nzchar(feature)) {
    .gradepath_abort("`feature` must be a length-1 non-empty character string.")
  }

  .gradepath_abort(
    "%s is not yet implemented.",
    feature,
    class = c("gradepath_not_implemented", "gradepath_error")
  )
}

# Control-surface validation failures. Subclass of gradepath_error so generic
# handlers still catch it; carries the v2-named `gp_control_error` class so the
# control object's own validation can be discriminated from runtime status
# codes (Chapter 13 §status separates validation errors from solver/run status).
.gradepath_abort_control <- function(message, ..., call = FALSE) {
  .gradepath_abort(
    message,
    ...,
    class = c("gp_control_error", "gradepath_control_error", "gradepath_error"),
    call = call
  )
}

# Input-validation failures (the Chapter 13 `INPUT_ERROR` family): malformed
# columns, non-finite estimates, duplicated ids, invalid units. Provided here
# so the data/seam Steps can raise a single, greppable class.
.gradepath_abort_input <- function(message, ..., call = FALSE) {
  .gradepath_abort(
    message,
    ...,
    class = c("gp_input_error", "gradepath_input_error", "gradepath_error"),
    call = call
  )
}

# ---------------------------------------------------------------------------
# Provenance  (reused verbatim from v1)
# ---------------------------------------------------------------------------

# Build a provenance named list from `...`. Every entry must be named; NULL
# entries are dropped (so optional fields can be passed unconditionally). The
# canonical run-level audit stamp every `gp_*` object carries (Chapter 9
# §provenance): validated only as a named list, contents advisory/extensible.
#
# Callers append run facts -- producer, timestamp, R version, package version,
# seed. The standard stamp is produced by `.gradepath_new_provenance()` with the
# common keys filled, e.g.:
#   .gradepath_new_provenance(
#     producer        = "gp_control",
#     built_at        = Sys.time(),
#     r_version       = R.version.string,
#     package_version = .gradepath_package_version(),
#     seed            = seed
#   )
.gradepath_new_provenance <- function(...) {
  provenance <- list(...)

  if (length(provenance) == 0L) {
    return(list())
  }

  names_provenance <- names(provenance)

  if (is.null(names_provenance) || anyNA(names_provenance) || any(!nzchar(names_provenance))) {
    .gradepath_abort("All provenance entries must be named.")
  }

  provenance[!vapply(provenance, is.null, logical(1))]
}

# Resolve the installed package version as a plain string, tolerant of being
# sourced outside an installed package (returns NA_character_ so provenance
# stamping never errors during dev `load_all()` or a throwaway source()).
.gradepath_package_version <- function() {
  v <- tryCatch(
    as.character(utils::packageVersion("gradepath")),
    error = function(e) NA_character_
  )
  v
}

# Convenience: the standard provenance stamp for a producing function. Keeps the
# five common run facts in one place so every constructor stamps identically.
.gradepath_provenance_stamp <- function(producer, ..., seed = NULL) {
  producer <- .gradepath_validate_scalar_character(producer, "producer")
  .gradepath_new_provenance(
    producer        = producer,
    built_at        = Sys.time(),
    r_version       = R.version.string,
    package_version = .gradepath_package_version(),
    seed            = seed,
    ...
  )
}

# ---------------------------------------------------------------------------
# Structural validators  (named-fields + class -- reused from v1)
# ---------------------------------------------------------------------------

.gradepath_validate_list_class <- function(x, class_name) {
  if (!is.list(x) || !inherits(x, class_name)) {
    .gradepath_abort("Expected an object of class `%s`.", class_name)
  }

  invisible(x)
}

# Exact top-level field set in canonical order. The frozen-schema enforcer the
# typed objects rely on: an object with the right slots in the wrong order is
# rejected (Chapter 9 §schemas).
.gradepath_validate_named_fields <- function(x, required, class_name) {
  if (!is.character(required) || anyNA(required) || any(!nzchar(required))) {
    .gradepath_abort_internal("`required` must be a non-empty character vector.")
  }

  fields <- names(x)

  if (is.null(fields)) {
    fields <- character()
  }

  if (!identical(fields, required)) {
    .gradepath_abort(
      "%s must have exact top-level fields in canonical order: %s.",
      class_name,
      paste(required, collapse = ", ")
    )
  }

  invisible(x)
}

# ---------------------------------------------------------------------------
# Scalar / vector validators  (character / numeric / logical / probability /
# integerish, with allowed= and bounds -- reused from v1)
# ---------------------------------------------------------------------------

.gradepath_validate_scalar_character <- function(x, name, allowed = NULL, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .gradepath_abort("`%s` must be a length-1 non-empty character value.", name)
  }

  if (!is.null(allowed) && !x %in% allowed) {
    .gradepath_abort(
      "`%s` must be one of: %s.",
      name,
      paste(allowed, collapse = ", ")
    )
  }

  x
}

.gradepath_validate_scalar_logical <- function(x, name, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .gradepath_abort("`%s` must be a length-1 logical value.", name)
  }

  x
}

.gradepath_validate_scalar_numeric <- function(x,
                                               name,
                                               allow_null = FALSE,
                                               finite = TRUE,
                                               lower = -Inf,
                                               upper = Inf,
                                               include_lower = TRUE,
                                               include_upper = TRUE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    .gradepath_abort("`%s` must be a length-1 numeric value.", name)
  }

  if (isTRUE(finite) && !is.finite(x)) {
    .gradepath_abort("`%s` must be finite.", name)
  }

  lower_ok <- if (is.infinite(lower)) {
    TRUE
  } else if (isTRUE(include_lower)) {
    x >= lower
  } else {
    x > lower
  }

  upper_ok <- if (is.infinite(upper)) {
    TRUE
  } else if (isTRUE(include_upper)) {
    x <= upper
  } else {
    x < upper
  }

  if (!lower_ok || !upper_ok) {
    lower_bracket <- if (isTRUE(include_lower)) "[" else "("
    upper_bracket <- if (isTRUE(include_upper)) "]" else ")"

    .gradepath_abort(
      "`%s` must lie in %s%s, %s%s.",
      name,
      lower_bracket,
      lower,
      upper,
      upper_bracket
    )
  }

  as.numeric(x)
}

.gradepath_validate_numeric_vector <- function(x, name, allow_null = FALSE, finite = TRUE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.numeric(x) || length(x) < 1L || anyNA(x)) {
    .gradepath_abort("`%s` must be a non-empty numeric vector.", name)
  }

  if (isTRUE(finite) && any(!is.finite(x))) {
    .gradepath_abort("`%s` must be finite.", name)
  }

  as.numeric(x)
}

.gradepath_validate_named_list <- function(x, name, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.list(x)) {
    .gradepath_abort("`%s` must be a list.", name)
  }

  if (length(x) == 0L) {
    return(list())
  }

  names_x <- names(x)

  if (is.null(names_x) || anyNA(names_x) || any(!nzchar(names_x))) {
    .gradepath_abort("`%s` must be a named list.", name)
  }

  x
}

# Non-negative integerish scalar with a configurable lower bound. Used by the
# control object's `seed` field.
.gradepath_validate_integerish <- function(x, name, min = 0L, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  x <- .gradepath_validate_scalar_numeric(
    x,
    name = name,
    finite = TRUE,
    lower = min,
    include_lower = TRUE
  )

  if (abs(x - round(x)) > sqrt(.Machine$double.eps)) {
    .gradepath_abort("`%s` must be an integer >= %s.", name, min)
  }

  as.integer(round(x))
}

# Probability scalar -- thin wrapper over the bounded numeric validator. Default
# bracket is the OPEN unit interval (0, 1); callers widen via include_*.
.gradepath_validate_probability <- function(x,
                                            name,
                                            lower = 0,
                                            upper = 1,
                                            include_lower = FALSE,
                                            include_upper = FALSE,
                                            allow_null = FALSE) {
  .gradepath_validate_scalar_numeric(
    x,
    name = name,
    allow_null = allow_null,
    finite = TRUE,
    lower = lower,
    upper = upper,
    include_lower = include_lower,
    include_upper = include_upper
  )
}

.gradepath_validate_data_frame <- function(x, name, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.data.frame(x)) {
    .gradepath_abort("`%s` must be a data.frame.", name)
  }

  x
}

.gradepath_validate_character_vector <- function(x,
                                                 name,
                                                 allow_null = FALSE,
                                                 unique = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.character(x) || length(x) < 1L || anyNA(x) || any(!nzchar(x))) {
    .gradepath_abort("`%s` must be a non-empty character vector.", name)
  }

  if (isTRUE(unique) && anyDuplicated(x)) {
    .gradepath_abort("`%s` must contain unique values.", name)
  }

  x
}

.gradepath_validate_warning_vector <- function(x, name = "warnings", allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
    .gradepath_abort("`%s` must be a character vector of non-empty entries.", name)
  }

  x
}

.gradepath_validate_numeric_matrix <- function(x, name, allow_null = FALSE, finite = TRUE) {
  if (is.null(x)) {
    if (allow_null) {
      return(NULL)
    }

    .gradepath_abort("`%s` must not be NULL.", name)
  }

  if (!is.matrix(x) || !is.numeric(x) || length(x) < 1L || anyNA(x)) {
    .gradepath_abort("`%s` must be a non-empty numeric matrix.", name)
  }

  if (isTRUE(finite) && any(!is.finite(x))) {
    .gradepath_abort("`%s` must be finite.", name)
  }

  x
}

.gradepath_validate_required_keys <- function(x, required, name) {
  if (!is.character(required) || anyNA(required) || any(!nzchar(required))) {
    .gradepath_abort_internal("`required` must be a non-empty character vector.")
  }

  missing <- setdiff(required, names(x))

  if (length(missing) > 0L) {
    .gradepath_abort(
      "`%s` must include required keys: %s.",
      name,
      paste(required, collapse = ", ")
    )
  }

  invisible(x)
}

# ---------------------------------------------------------------------------
# Structural / combinatorial helpers
# ---------------------------------------------------------------------------

# Linear (column-major) indices of the diagonal of an n x n matrix:
#   seq.int(1, n*n, by = n + 1) == c(1, n+2, 2n+3, ...).
# The grade IP uses this to project d_ii = 0 and to read the
# diagonal of the pairwise matrix without forming `diag()`. Internal (Chapter 13
# exports only the four verbs + monolith + gp_control + accessors; the design
# lists `gp_diag_positions` as a grade-IP internal).
#
# `n` must be a single non-negative integer-valued number. n == 0 returns an
# integer(0); n == 1 returns 1L.
gp_diag_positions <- function(n) {
  n <- .gradepath_validate_integerish(n, "n", min = 0L)

  if (n == 0L) {
    return(integer(0))
  }

  seq.int(1L, n * n, by = n + 1L)
}

# Canonical relabel of an integer/numeric grade (or rank) vector to contiguous
# labels 1..k by SORTED-UNIQUE value: the smallest distinct value -> 1, the next
# -> 2, and so on, preserving ties and order.
#
#   monotone_relabel(c(5, 5, 2, 9)) == c(2, 2, 1, 3)
#   monotone_relabel(c(3, 1, 2))    == c(3, 1, 2)
#
# This is the contiguous-1..k normalizer the typed-object validators assert on
# `grade` columns (Chapter 9: "normalizes `grade` to contiguous integers
# starting at 1"), and the harness uses it to align label sets before comparison.
# It is order-preserving and monotone in the input values, so a
# grading's structure is unchanged -- only the labels are made canonical.
#
# Returns an integer vector the same length as `g`. Internal (not in the
# Chapter 13 export surface).
monotone_relabel <- function(g) {
  g <- .gradepath_validate_numeric_vector(g, "g")

  # match() against the sorted-unique values yields the rank-of-value, which is
  # exactly the contiguous 1..k label by ascending value.
  as.integer(match(g, sort(unique(g))))
}

# ---------------------------------------------------------------------------
# Stable row log-sum-exp
# ---------------------------------------------------------------------------

# Numerically stable row-wise log-sum-exp of a numeric matrix M:
#   out[i] = log(sum_j exp(M[i, j]))
# computed as m_i + log(sum_j exp(M[i, j] - m_i)) with m_i = max_j M[i, j], so
# it never overflows for large entries and returns -Inf for an all -Inf row.
#
# This is the W-recompute / posterior-weight workhorse (the W seam): grid
# log-likelihoods are turned into normalized log-posteriors via row
# log-sum-exp before the safe-normalize. Mirrors ebrecipe's `.eb_row_log_sum_exp`
# primitive (the binding map) so the native path and the cross-check agree.
#
# `M` must be a numeric matrix (a single row is fine). Returns a numeric vector
# of length nrow(M). Entries may be -Inf (empty-mass row); +Inf/NA in the input
# are rejected by the matrix validator unless finite = FALSE is wired upstream --
# here we allow -Inf explicitly because it is a legitimate "zero mass" code.
.gp_row_logsumexp <- function(M) {
  if (is.null(dim(M))) {
    M <- matrix(M, nrow = 1L)
  }
  if (!is.matrix(M) || !is.numeric(M)) {
    .gradepath_abort("`M` must be a numeric matrix.")
  }
  if (anyNA(M)) {
    .gradepath_abort("`M` must not contain NA.")
  }

  # Row maxima (na.rm not needed; NA already rejected). An all -Inf row gives
  # m = -Inf; the shifted exponentials are exp(-Inf - -Inf) = exp(NaN) which we
  # special-case to a clean -Inf result.
  m <- apply(M, 1L, max)

  finite_max <- is.finite(m)
  out <- rep(-Inf, nrow(M))

  if (any(finite_max)) {
    Mf <- M[finite_max, , drop = FALSE]
    mf <- m[finite_max]
    out[finite_max] <- mf + log(rowSums(exp(Mf - mf)))
  }

  out
}
