# The pruned gradepath v2 control surface: gp_control().
#
# v1's `gradepath_control()` carried twelve fields, five of which encoded the
# UX weaknesses. v2 keeps exactly SIX user-facing computational
# fields and drops the rest (Chapter 9 §control, Chapter 13 §uxfixes):
#
#   KEEP : lambda_grid, backend, precision_rule, interval_level,
#          solver_options (+ ergonomic time_limit / mip_gap aliases), seed
#   DROP : workflow, replication_mode, target_bundle, selection_rule,
#          backend_fallback, and the entire claim-tier vocabulary
#          (claim_bearing / claim_scope / backend_claim_tier).
#
# `backend` defaults to "gurobi" as the planned parity backend (the same solver
# KRW used; NOT v1's environment-dependent "auto"), read through
# getOption("gradepath.backend", "gurobi") so R/zzz.R remains consistent with
# this surface (Chapter 10; zzz.R mirrors -- gp_control() is the single source
# of truth). The accepted set is gurobi + explicit open backends (Chapter 10
# §backends, Chapter 13 §grade-path): "highs" uses the highs package directly;
# "glpk"/"symphony" route through ROI plugins. v1's "roi_*" spellings are
# accepted as aliases for back-compatibility with cached fixtures.
#
# Selection lambda is NOT here -- later selection APIs keep it as a first-class
# argument. The internal schema_version + provenance stamp are retained.

# ---------------------------------------------------------------------------
# Field schema + accepted vocabularies
# ---------------------------------------------------------------------------

# Canonical top-level field order for the stored gp_control object. Six user
# fields + the two internal audit slots. `.gradepath_validate_named_fields()`
# asserts this exact order.
.gp_control_fields <- c(
  "lambda_grid",
  "backend",
  "precision_rule",
  "interval_level",
  "solver_options",
  "seed",
  "schema_version",
  "provenance"
)

# Accepted grade-IP backends. The canonical v2 names are the short forms
# (Chapter 13 §grade-path); the "roi_*" forms are accepted aliases (v1 surface /
# cached fixtures) and normalized to the canonical short name on construction.
.gp_control_backends_canonical <- c("gurobi", "highs", "glpk", "symphony")

.gp_control_backend_aliases <- c(
  gurobi       = "gurobi",
  highs        = "highs",
  glpk         = "glpk",
  symphony     = "symphony",
  roi_highs    = "highs",
  roi_glpk     = "glpk",
  roi_symphony = "symphony"
)

# Accepted precision rules. "none" keeps the raw r-scale estimates; "krw_gmm"
# selects the gradepath-native beta-GMM standardization (Chapter 7;
# Chapter 9 lists precision_rule values as "none" | "krw_gmm"). v1's
# "log_linear" name is NOT carried -- it named ebrecipe's NLLS estimand, which
# the native core deliberately bypasses.
.gp_control_precision_rules <- c("none", "krw_gmm")

# Default lambda grid. The spec fixes the v1-surface default
# seq(0, 1, by = 0.01); the warmstarted *operational* default grid (the small
# ~11-point grid Chapter 13 uxfixes describes) belongs to later grade-path
# execution when lambda_grid = NULL, not here -- so the control object's stored grid
# stays the full reference grid unless the user narrows it. It MUST contain the
# two parity anchors 0.25 (KRW's selection) and 1.00 (the endpoint warm-start).
.gp_control_default_lambda_grid <- seq(0, 1, by = 0.01)

# The two lambda values every grid must contain (parity anchors): the published
# selection 0.25 and the endpoint 1.00 used for warm-starting the sweep.
.gp_control_required_lambda <- c(0.25, 1.00)

# ---------------------------------------------------------------------------
# Solver-option normalization  (ports v1, keeps compatibility option names)
# ---------------------------------------------------------------------------

# Re-raise any validation failure from `expr` as a control-surface error, so a
# bad solver-option *value* (whether supplied via an ergonomic alias or inside
# `solver_options`) carries the uniform `gp_control_error` class rather than the
# generic numeric-validator class. The original message is preserved.
.gp_control_as_control_error <- function(expr) {
  withCallingHandlers(
    expr,
    gradepath_error = function(cnd) {
      # Already a control error -> let it propagate unchanged.
      if (inherits(cnd, "gp_control_error")) {
        return(invisible(NULL))
      }
      .gradepath_abort_control("%s", conditionMessage(cnd))
    }
  )
}

# Validate a single time-limit value (seconds): finite, strictly positive.
# Range/shape failures are surfaced as control errors.
.gp_control_time_limit <- function(x, name) {
  x <- .gp_control_as_control_error(
    .gradepath_validate_scalar_numeric(
      x,
      name = name,
      finite = TRUE,
      lower = 0,
      include_lower = FALSE
    )
  )
  as.numeric(x)
}

# Validate a single MIP-gap value: finite, in the closed [0, 1]. Range/shape
# failures are surfaced as control errors.
.gp_control_mip_gap <- function(x, name) {
  x <- .gp_control_as_control_error(
    .gradepath_validate_probability(
      x,
      name = name,
      lower = 0,
      upper = 1,
      include_lower = TRUE,
      include_upper = TRUE
    )
  )
  as.numeric(x)
}

# Validate a user-supplied solver_options list. NULL or empty -> list(). Must be
# a fully-named list. Any portable keys present inside it (the canonical
# `time_limit` / `mip_gap`, OR the compatibility spellings `max_time` /
# `mip_rel_gap`) are range-checked in place. Backend-specific keys we don't
# recognize pass through untouched -- gp_control is portable and does not
# enumerate every backend's option surface.
.gp_control_validate_solver_options <- function(x, name = "solver_options") {
  if (is.null(x)) {
    return(list())
  }
  if (!is.list(x)) {
    .gradepath_abort_control("`%s` must be a named list.", name)
  }
  if (length(x) == 0L) {
    return(list())
  }

  names_x <- names(x)
  if (is.null(names_x) || anyNA(names_x) || any(!nzchar(names_x))) {
    .gradepath_abort_control("`%s` must be a named list.", name)
  }

  # Portable / known time-limit-like keys.
  for (k in intersect(c("time_limit", "max_time"), names_x)) {
    x[[k]] <- .gp_control_time_limit(x[[k]], paste0(name, "$", k))
  }
  # Portable / known gap-like keys.
  for (k in intersect(c("mip_gap", "mip_rel_gap"), names_x)) {
    x[[k]] <- .gp_control_mip_gap(x[[k]], paste0(name, "$", k))
  }

  x
}

# Merge the ergonomic top-level aliases (`time_limit`, `mip_gap`) into the
# solver_options list, normalizing values and refusing double-specification.
#
# Double-specification is detected against BOTH the canonical key and that
# compatibility spelling: passing `time_limit =` while also putting
# `max_time` in solver_options is an error (they mean the same time-limit knob),
# as is `time_limit` + `solver_options$time_limit`. Same for mip_gap /
# mip_rel_gap. This is the v1 "do not supply in both" guard, widened to the
# compatibility aliases.
.gp_control_normalize_solver_options <- function(solver_options,
                                                 time_limit = NULL,
                                                 mip_gap = NULL) {
  solver_options <- .gp_control_validate_solver_options(solver_options)
  solver_options <- solver_options[!vapply(solver_options, is.null, logical(1))]

  # Build the portable additions from the ergonomic aliases.
  portable <- list()
  if (!is.null(time_limit)) {
    portable$time_limit <- .gp_control_time_limit(time_limit, "time_limit")
  }
  if (!is.null(mip_gap)) {
    portable$mip_gap <- .gp_control_mip_gap(mip_gap, "mip_gap")
  }

  # Reject double-specification, treating each portable key and its backend
  # alias as the same knob.
  knob_aliases <- list(
    time_limit = c("time_limit", "max_time"),
    mip_gap    = c("mip_gap", "mip_rel_gap")
  )
  for (knob in names(portable)) {
    clash <- intersect(knob_aliases[[knob]], names(solver_options))
    if (length(clash) > 0L) {
      .gradepath_abort_control(
        "Do not supply `%s` both as an ergonomic argument and as %s in `solver_options`.",
        knob,
        paste(sprintf("`%s`", clash), collapse = " / ")
      )
    }
  }

  c(portable, solver_options)
}

# ---------------------------------------------------------------------------
# Backend / lambda-grid validators
# ---------------------------------------------------------------------------

# Validate and canonicalize the backend string. Accepts the canonical short
# names and the "roi_*" aliases; returns the canonical short name. No
# availability check here (that is the solver verb's job at call time, Chapter
# 10 / 13 SOLVER_BACKEND_UNAVAILABLE) and -- crucially -- no claim-tiering.
.gp_control_validate_backend <- function(x, name = "backend") {
  x <- .gradepath_validate_scalar_character(x, name)
  if (!x %in% names(.gp_control_backend_aliases)) {
    .gradepath_abort_control(
      "`%s` must be one of: %s.",
      name,
      paste(.gp_control_backends_canonical, collapse = ", ")
    )
  }
  unname(.gp_control_backend_aliases[[x]])
}

# Validate the lambda grid: non-empty numeric, all in [0, 1], unique, strictly
# increasing, and containing both parity anchors {0.25, 1.00}.
# Core lambda-grid shape check: a non-empty finite numeric vector, all in
# [0, 1], unique, strictly increasing. This is the structural contract a lambda
# grid must satisfy ANYWHERE (gp_control's grid AND a gp_grade_path's solved
# grid, Chapter 9 §grade-path). It does NOT require the parity anchors -- a
# user may legitimately solve a sub-grid such as c(0.1, 0.5, 0.9). Failures are
# generic `gradepath_error` so the helper is reusable outside the control
# surface; the caller may wrap in a control error where appropriate.
.gp_validate_lambda_grid_core <- function(x, name = "lambda_grid") {
  x <- .gradepath_validate_numeric_vector(x, name)

  if (any(x < 0 | x > 1)) {
    .gradepath_abort("`%s` must lie inside [0, 1].", name)
  }
  if (anyDuplicated(x)) {
    .gradepath_abort("`%s` must contain unique values.", name)
  }
  if (is.unsorted(x, strictly = TRUE)) {
    .gradepath_abort("`%s` must be strictly increasing.", name)
  }

  as.numeric(x)
}

# The gp_control lambda-grid validator: the core shape check PLUS the parity
# anchors {0.25, 1.00}. gp_control's grid is the *reference* sweep, so it must
# carry KRW's selection (0.25) and the warm-start endpoint (1.00). A solved
# gp_grade_path grid uses the anchor-free core above. Failures are control
# errors.
.gp_control_validate_lambda_grid <- function(x, name = "lambda_grid") {
  x <- .gp_control_as_control_error(.gp_validate_lambda_grid_core(x, name))

  # Parity anchors must be present (within floating tolerance).
  for (anchor in .gp_control_required_lambda) {
    if (!any(abs(x - anchor) < 1e-8)) {
      .gradepath_abort_control(
        "`%s` must contain the parity anchor %s.",
        name,
        format(anchor)
      )
    }
  }

  as.numeric(x)
}

# ---------------------------------------------------------------------------
# Constructor + validator
# ---------------------------------------------------------------------------

# Low-level constructor: assembles the structure with the canonical field order
# and the v2 schema tag + provenance stamp, then validates. Internal; user code
# calls gp_control().
new_gp_control <- function(lambda_grid,
                           backend,
                           precision_rule,
                           interval_level,
                           solver_options,
                           seed) {
  x <- structure(
    list(
      lambda_grid    = lambda_grid,
      backend        = backend,
      precision_rule = precision_rule,
      interval_level = interval_level,
      solver_options = solver_options,
      seed           = seed,
      schema_version = .gradepath_schema_version,
      provenance     = .gradepath_provenance_stamp("gp_control", seed = seed)
    ),
    class = c("gp_control", "list")
  )
  validate_gp_control(x)
}

# Validator: asserts class, exact field order, and each field's type/bounds.
# Idempotent on a well-formed object.
validate_gp_control <- function(x) {
  .gradepath_validate_list_class(x, "gp_control")
  .gradepath_validate_named_fields(x, .gp_control_fields, "gp_control")

  .gp_control_validate_lambda_grid(x$lambda_grid)
  .gp_control_validate_backend(x$backend)
  .gradepath_validate_scalar_character(
    x$precision_rule,
    "precision_rule",
    allowed = .gp_control_precision_rules
  )
  .gradepath_validate_probability(
    x$interval_level,
    "interval_level",
    lower = 0,
    upper = 1,
    include_lower = FALSE,
    include_upper = FALSE
  )
  .gradepath_validate_named_list(x$solver_options, "solver_options")
  if (!is.null(x$seed)) {
    .gradepath_validate_integerish(x$seed, "seed", min = 0L)
  }
  .gradepath_validate_scalar_character(x$schema_version, "schema_version")
  .gradepath_validate_named_list(x$provenance, "provenance", allow_null = FALSE)

  invisible(x)
}

# ---------------------------------------------------------------------------
# gp_control()  -- the user constructor
# ---------------------------------------------------------------------------

#' Create a gradepath run-control object
#'
#' `gp_control()` bundles the run settings the gradepath verbs read -- the
#' frontier penalty grid, the grade-IP solver backend, the precision-handling
#' rule, the credible-interval level, and solver runtime knobs -- into one
#' validated `gp_control` object. Every argument has a sane default, so
#' `gp_control()` with no arguments is a complete, valid configuration. Build one,
#' tweak the fields you care about, and pass it as the `control =` argument to
#' [krw_report_card()], [gp_grade_path()], or [gp_preview()].
#'
#' @param lambda_grid Numeric vector; the frontier penalty values the grade path
#'   is solved at. Must be unique, strictly increasing, and lie in `[0, 1]`, and
#'   must contain the parity anchors `0.25` (KRW's published selection) and `1.00`
#'   (the endpoint warm-start). Defaults to the reference grid `seq(0, 1, by =
#'   0.01)`; passing `NULL` restores that reference grid, and later grade-path
#'   execution may narrow it.
#' @param backend Character scalar; the planned grade-IP solver backend, one of
#'   `"gurobi"`, `"highs"`, `"glpk"`, or `"symphony"`. Defaults to the value of
#'   `getOption("gradepath.backend", "gurobi")` -- the Gurobi parity backend
#'   (the same solver KRW used), with `"highs"`, `"glpk"`, and `"symphony"` as
#'   license-free open alternatives (`"highs"` calls the `highs` package
#'   directly; `"glpk"` and `"symphony"` route through ROI plugins). The v1
#'   `"roi_highs"` / `"roi_glpk"` / `"roi_symphony"` spellings are accepted and
#'   normalized to the short names. There is no `"auto"`; backend availability is
#'   checked by the solver at call time, not here.
#' @param precision_rule Character scalar; the precision-handling rule. `"none"`
#'   keeps the raw estimate (r) scale; `"krw_gmm"` selects the gradepath-native
#'   KRW beta-GMM standardization. Defaults to `"none"`.
#' @param interval_level Numeric scalar; the credible-interval level stored on
#'   the downstream posterior and report-card objects. A probability in the open
#'   interval `(0, 1)`. Defaults to `0.90`.
#' @param solver_options Named list; backend-specific solver options. The
#'   portable keys `time_limit` and `mip_gap` (and their compatibility spellings
#'   `max_time` and `mip_rel_gap`) are range-checked in place; any other key
#'   passes through untouched. `mip_gap` is honored by Gurobi, HiGHS, and
#'   SYMPHONY but not by the current GLPK route. Defaults to an empty list.
#' @param time_limit Numeric scalar or `NULL`; an ergonomic top-level alias for
#'   the solver time limit in seconds (finite, strictly positive), merged into
#'   `solver_options$time_limit`. Supplying it alongside `time_limit` (or its
#'   `max_time` spelling) inside `solver_options` is an error. Defaults to `NULL`
#'   (no time limit set here).
#' @param mip_gap Numeric scalar or `NULL`; an ergonomic top-level alias for the
#'   relative MIP optimality gap in `[0, 1]`, merged into
#'   `solver_options$mip_gap`. Supplying it alongside `mip_gap` (or its
#'   `mip_rel_gap` spelling) inside `solver_options` is an error. Defaults to
#'   `NULL` (no gap target set here).
#' @param seed Integer scalar or `NULL`; a non-negative seed for the stochastic
#'   (two-level Monte-Carlo) paths. Defaults to `NULL` (unseeded).
#'
#' @return A validated S3 object of class `c("gp_control", "list")`: a named list
#'   carrying the six user-facing computational fields plus two internal audit
#'   slots. \describe{
#'   \item{`lambda_grid`}{Numeric vector; the validated frontier penalty grid.}
#'   \item{`backend`}{Character scalar; the canonical short backend name.}
#'   \item{`precision_rule`}{Character scalar; `"none"` or `"krw_gmm"`.}
#'   \item{`interval_level`}{Numeric scalar in `(0, 1)`; the credible-interval
#'     level.}
#'   \item{`solver_options`}{Named list; the normalized solver options, with any
#'     ergonomic `time_limit` / `mip_gap` aliases merged in.}
#'   \item{`seed`}{Non-negative integer or `NULL`; the stochastic-path seed.}
#'   \item{`schema_version`, `provenance`}{Internal audit slots.}
#' }
#'
#' @details
#' The surface is the pruned v2 control of six computational fields; v1's
#' `workflow`, `replication_mode`, `target_bundle`, `selection_rule`,
#' `backend_fallback`, and the legacy ranking vocabulary are gone. Selection of
#' the frontier penalty lives with the grading verbs, not here. Portable solver
#' runtime knobs may be passed either inside `solver_options` or via the
#' ergonomic top-level aliases `time_limit` and `mip_gap`; supplying the same
#' knob twice -- including a compatibility spelling alongside its ergonomic alias
#' -- is an error.
#'
#' @examples
#' # Runnable with no arguments and no solve: the friendly defaults give backend
#' # "gurobi", interval_level 0.90, and the reference lambda grid.
#' ctl <- gp_control()
#'
#' # A customized control on an open backend: a narrow grid (still carrying the
#' # 0.25 and 1.00 anchors) with an ergonomic time limit merged into
#' # solver_options.
#' ctl2 <- gp_control(
#'   lambda_grid = c(0, 0.25, 0.5, 1),
#'   backend     = "highs",
#'   time_limit  = 30
#' )
#' ctl2$solver_options$time_limit
#'
#' @seealso [krw_report_card()], [gp_grade_path()], [gp_preview()]
#' @family gradepath-pipeline
#' @export
gp_control <- function(lambda_grid    = seq(0, 1, by = 0.01),
                       backend        = getOption("gradepath.backend", "gurobi"),
                       precision_rule = c("none", "krw_gmm"),
                       interval_level = 0.90,
                       solver_options = NULL,
                       time_limit     = NULL,
                       mip_gap        = NULL,
                       seed           = NULL) {
  # Route a bad `precision_rule` through the package's classed control error
  # rather than match.arg()'s base simpleError, so EVERY gp_control validation
  # path raises a `gp_control_error` (a subclass of `gradepath_error`). The
  # no-arg default still resolves to the first formal ("none").
  precision_rule <- tryCatch(
    match.arg(precision_rule),
    error = function(e) {
      .gradepath_abort_control(
        "`precision_rule` must be one of: %s.",
        paste(.gp_control_precision_rules, collapse = ", ")
      )
    }
  )

  if (is.null(lambda_grid)) {
    lambda_grid <- .gp_control_default_lambda_grid
  }
  lambda_grid <- .gp_control_validate_lambda_grid(lambda_grid, "lambda_grid")

  backend <- .gp_control_validate_backend(backend, "backend")

  interval_level <- .gradepath_validate_probability(
    interval_level,
    "interval_level",
    lower = 0,
    upper = 1,
    include_lower = FALSE,
    include_upper = FALSE
  )

  seed <- .gradepath_validate_integerish(seed, "seed", min = 0L, allow_null = TRUE)

  solver_options <- .gp_control_normalize_solver_options(
    solver_options = solver_options,
    time_limit     = time_limit,
    mip_gap        = mip_gap
  )

  new_gp_control(
    lambda_grid    = lambda_grid,
    backend        = backend,
    precision_rule = precision_rule,
    interval_level = interval_level,
    solver_options = solver_options,
    seed           = seed
  )
}
