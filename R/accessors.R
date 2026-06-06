# The accessor family and base generics over `gp_fit`.
#
# This is gradepath v2's FIRST exported behavioural surface beyond gp_control().
# It implements the `ashr`-style accessor pattern adopted by GP-DEC-13-A
# (Chapter 13 §sec-ch13-patterns: "gp_fit is read through get_grades(),
# get_report_card(), get_pairwise(), get_prior(), get_posterior(),
# get_control() — *never* by reaching into `$` slots"; "the base-generic
# companions are also provided: coef(fit) returns the grade vector ... and
# as.data.frame(fit) returns the report-card data frame"). The full accessor set
# is specified by the .spec block in Chapter 14 §ux (the accessor + base-generic
# surface over gp_fit).
#
# Eight exported functions, ALL read-only, ALL S3 on `gp_fit`:
#   get_grades(fit)        -> integer vector, NAMED by ids
#   get_report_card(fit)   -> gp_report_card        (fit$report_card)
#   get_pairwise(fit)      -> gp_pairwise           (fit$pairwise)
#   get_prior(fit)         -> gp_prior              (fit$prior)
#   get_posterior(fit)     -> gp_posterior          (fit$posterior)
#   get_control(fit)       -> gp_control            (fit$control)
#   coef.gp_fit(object)    -> grade vector (== get_grades)
#   as.data.frame.gp_fit(x)-> report-card table (gp_report_card$table)
#
# DESIGN DECISION — generics, not plain functions. Chapter 13 §patterns models
# the family on `ashr`'s get_*() accessors, which are S3 generics so a downstream
# object (or a future gradepath object) can add a method. The six get_*() are
# therefore declared as S3 *generics* with a `.gp_fit` method each (and a
# `.default` that errors cleanly on non-gp_fit input). coef() and as.data.frame()
# are base generics, so here we only register `.gp_fit` methods. This matches the
# idiom in this exact lineage and keeps the surface extensible at zero present
# cost. See SUMMARY.md "generic-vs-plain decision".
#
# Reuses (do NOT redefine; built Steps 1.2-1.5):
#   utils-validate.R          : .gradepath_abort
#   class-objects-output.R    : gp_fit slots (ids, report_card, pairwise, prior,
#                               posterior, control, selected_grade)
#   class-objects-decision.R  : gp_grade_fit$assignment = data.frame(id, grade)

# ---------------------------------------------------------------------------
# internal guard — clean error on non-gp_fit input
# ---------------------------------------------------------------------------

# Shared off-contract guard. Every accessor's `.default` method routes here so a
# non-gp_fit input fails with one consistent, classed, greppable message rather
# than a cryptic `$` failure deep in the body.
#' @keywords internal
.gp_assert_fit <- function(x, what) {
  .gradepath_abort(
    "`%s()` requires a `gp_fit` object; got <%s>.",
    what, paste(class(x), collapse = "/")
  )
}

# ===========================================================================
# get_grades — the headline accessor: an integer vector NAMED by ids
# ===========================================================================

#' Extract per-unit grades from a gradepath fit
#'
#' Returns the selected grading as an integer vector, one entry per unit, named by
#' the unit ids. The grades are read from the selected `gp_grade_fit`
#' (`fit$selected_grade$assignment`) and aligned to the fit's canonical `ids` order,
#' so `get_grades(fit)[id]` is the grade of unit `id`. Use it as the headline
#' read of a `gp_fit`; [coef.gp_fit()] is the base-generic synonym and
#' [get_report_card()] returns the same grades alongside the posterior summaries.
#'
#' Per `GP-DEC-13-A` (Chapter 13 patterns) a `gp_fit` is read through the accessor
#' family -- `get_grades()`, [get_report_card()], [get_pairwise()], [get_prior()],
#' [get_posterior()], [get_control()] -- never by reaching into `$` slots. The
#' accessor never recomputes; it only reads materialized slots. Grade labels are
#' integers in `{1, ..., n}` and carry no ranking-superiority statement of any kind.
#'
#' @param fit A materialized `gp_fit` object.
#' @param ... Unused; accepted for S3 method extensibility. No arguments are
#'   forwarded; supplying any has no effect.
#'
#' @return A named `integer` vector of length `J` (the number of units): the
#'   selected grade per unit, with `names` equal to `fit$ids`.
#'
#' @examples
#' # Instant: read grades off the bundled tiny fit (no solver, no Gurobi).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' head(get_grades(fit))   # named int vector: names = ids, values = grades
#'
#' @seealso [coef.gp_fit()], [get_report_card()], [get_pairwise()],
#'   [krw_report_card()], [gp_report_card()]
#' @family gradepath-accessors
#' @export
get_grades <- function(fit, ...) {
  UseMethod("get_grades")
}

#' @rdname get_grades
#' @export
get_grades.default <- function(fit, ...) {
  .gp_assert_fit(fit, "get_grades")
}

#' @rdname get_grades
#' @export
get_grades.gp_fit <- function(fit, ...) {
  assignment <- fit$selected_grade$assignment
  ids <- fit$ids

  # Align the grade column to the canonical `ids` order by id. In a validated
  # gp_fit, assignment$id is already identical to ids (validate_gp_grade_fit
  # asserts assignment$id == fit$ids, and validate_gp_fit asserts the shared
  # canonical order), so this match is the identity — but doing it by name makes
  # the accessor robust to any future relaxation of that invariant and makes the
  # id->grade alignment explicit rather than positional.
  grades <- assignment$grade[match(ids, as.character(assignment$id))]

  grades <- as.integer(grades)
  names(grades) <- ids
  grades
}

# ===========================================================================
# get_report_card / get_pairwise / get_prior / get_posterior / get_control
#   — thin, validated slot readers
# ===========================================================================

#' Extract the report card from a gradepath fit
#'
#' Returns the `gp_report_card` slot of a `gp_fit` unchanged -- the endpoint-sorted
#' per-unit table of grade label, posterior mean, and credible interval, produced by
#' [gp_report_card()] inside the pipeline. Reach for it when you want the full
#' presentation table; [get_grades()] returns just the grade vector and the
#' base-generic [as.data.frame.gp_fit()] returns this card's underlying data frame.
#'
#' Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and only reads
#' the materialized slot. Grade labels are integers in `{1, ..., n}` and carry no
#' ranking-superiority statement of any kind.
#'
#' @param fit A `gp_fit` object.
#' @param ... Unused; accepted for S3 method extensibility. No arguments are
#'   forwarded; supplying any has no effect.
#'
#' @return The `gp_report_card` object stored on `fit` (returned unchanged).
#'
#' @examples
#' # Instant: read the report-card object off the bundled tiny fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' get_report_card(fit)   # <gp_report_card>
#'
#' @seealso [as.data.frame.gp_fit()], [get_grades()], [gp_report_card()],
#'   [krw_report_card()]
#' @family gradepath-accessors
#' @export
get_report_card <- function(fit, ...) {
  UseMethod("get_report_card")
}

#' @rdname get_report_card
#' @export
get_report_card.default <- function(fit, ...) {
  .gp_assert_fit(fit, "get_report_card")
}

#' @rdname get_report_card
#' @export
get_report_card.gp_fit <- function(fit, ...) {
  fit$report_card
}

#' Extract the pairwise outranking object from a gradepath fit
#'
#' Returns the `gp_pairwise` slot of a `gp_fit` unchanged: the \eqn{J \times J}
#' posterior outranking matrix \eqn{\pi_{ij}} (the posterior probability that unit
#' `i` outranks unit `j`) together with its cleanup and source metadata. Reach for it
#' to inspect the social-choice input that the grade path is solved over.
#'
#' Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and only reads
#' the materialized slot. The matrix encodes pairwise posterior order only and
#' carries no ranking-superiority statement of any kind.
#'
#' @param fit A `gp_fit` object.
#' @param ... Unused; accepted for S3 method extensibility. No arguments are
#'   forwarded; supplying any has no effect.
#'
#' @return The `gp_pairwise` object stored on `fit` (returned unchanged).
#'
#' @examples
#' # Instant: read the pairwise object off the bundled tiny fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' get_pairwise(fit)      # <gp_pairwise>
#'
#' @seealso [get_grades()], [get_report_card()], [gp_pairwise()]
#' @family gradepath-accessors
#' @export
get_pairwise <- function(fit, ...) {
  UseMethod("get_pairwise")
}

#' @rdname get_pairwise
#' @export
get_pairwise.default <- function(fit, ...) {
  .gp_assert_fit(fit, "get_pairwise")
}

#' @rdname get_pairwise
#' @export
get_pairwise.gp_fit <- function(fit, ...) {
  fit$pairwise
}

#' Extract the native prior from a gradepath fit
#'
#' Returns the `gp_prior` slot of a `gp_fit` unchanged. The native `gp_prior` is
#' `eb_prior`-shaped (it carries support, density, mean, scale, diagnostics, and
#' metadata) but remains a gradepath object rather than inheriting from `eb_prior`;
#' it is the deconvolved prior on the r-scale that the pipeline estimates. This is
#' how the empirical-Bayes front half -- absorbed into gradepath's native KRW core
#' -- is surfaced to the user.
#'
#' Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and only reads
#' the materialized slot.
#'
#' @param fit A `gp_fit` object.
#' @param ... Unused; accepted for S3 method extensibility. No arguments are
#'   forwarded; supplying any has no effect.
#'
#' @return The `gp_prior` object stored on `fit` (returned unchanged).
#'
#' @examples
#' # Instant: read the native prior off the bundled tiny fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' get_prior(fit)         # native gp_prior (eb_prior-shaped)
#'
#' @seealso [get_posterior()], [get_report_card()]
#' @family gradepath-accessors
#' @export
get_prior <- function(fit, ...) {
  UseMethod("get_prior")
}

#' @rdname get_prior
#' @export
get_prior.default <- function(fit, ...) {
  .gp_assert_fit(fit, "get_prior")
}

#' @rdname get_prior
#' @export
get_prior.gp_fit <- function(fit, ...) {
  fit$prior
}

#' Extract the native posterior from a gradepath fit
#'
#' Returns the `gp_posterior` slot of a `gp_fit` unchanged. The native
#' `gp_posterior` is `eb_posterior`-shaped (per-unit posterior means, standard
#' deviations, and credible intervals) but remains a gradepath object rather than
#' inheriting from `eb_posterior`. It holds the per-unit posterior summaries that the
#' report card presents.
#'
#' Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and only reads
#' the materialized slot.
#'
#' @param fit A `gp_fit` object.
#' @param ... Unused; accepted for S3 method extensibility. No arguments are
#'   forwarded; supplying any has no effect.
#'
#' @return The `gp_posterior` object stored on `fit` (returned unchanged).
#'
#' @examples
#' # Instant: read the native posterior off the bundled tiny fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' get_posterior(fit)     # native gp_posterior (eb_posterior-shaped)
#'
#' @seealso [get_prior()], [get_report_card()]
#' @family gradepath-accessors
#' @export
get_posterior <- function(fit, ...) {
  UseMethod("get_posterior")
}

#' @rdname get_posterior
#' @export
get_posterior.default <- function(fit, ...) {
  .gp_assert_fit(fit, "get_posterior")
}

#' @rdname get_posterior
#' @export
get_posterior.gp_fit <- function(fit, ...) {
  fit$posterior
}

#' Extract the run-control object from a gradepath fit
#'
#' Returns the `gp_control` slot of a `gp_fit` unchanged: the pruned run settings
#' actually used for the fit (`lambda_grid`, `backend`, `precision_rule`,
#' `interval_level`, `solver_options`, `seed`). Reach for it to see exactly which
#' [gp_control()] configuration produced the result, or to reuse those settings in a
#' follow-up run.
#'
#' Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and only reads
#' the materialized slot.
#'
#' @param fit A `gp_fit` object.
#' @param ... Unused; accepted for S3 method extensibility. No arguments are
#'   forwarded; supplying any has no effect.
#'
#' @return The `gp_control` object stored on `fit` (returned unchanged).
#'
#' @examples
#' # Instant: read the run-control object off the bundled tiny fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' get_control(fit)       # <gp_control>
#'
#' @seealso [gp_control()], [get_grades()]
#' @family gradepath-accessors
#' @export
get_control <- function(fit, ...) {
  UseMethod("get_control")
}

#' @rdname get_control
#' @export
get_control.default <- function(fit, ...) {
  .gp_assert_fit(fit, "get_control")
}

#' @rdname get_control
#' @export
get_control.gp_fit <- function(fit, ...) {
  fit$control
}

# ===========================================================================
# base-generic companions: coef() and as.data.frame()
# ===========================================================================

#' Grade vector from a gradepath fit (base-generic synonym for get_grades)
#'
#' `coef()` method for `gp_fit`: the "fitted summary" of a gradepath fit is its
#' per-unit grading, so `coef(fit)` returns the same named integer grade vector as
#' [get_grades()]. Provided because `coef()` is the base generic an R user reaches
#' for first (`GP-DEC-13-A`, Chapter 13 patterns; Chapter 14 ux). Grade labels are
#' integers in `{1, ..., n}` and carry no ranking-superiority statement of any kind.
#'
#' @param object A `gp_fit` object.
#' @param ... Unused; accepted for S3 method extensibility. No arguments are
#'   forwarded; supplying any has no effect.
#'
#' @return A named `integer` vector of grades, identical to `get_grades(object)`.
#'
#' @examples
#' # Instant: coef() is the grade vector of the bundled tiny fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' head(coef(fit))        # == head(get_grades(fit))
#'
#' @seealso [get_grades()], [as.data.frame.gp_fit()]
#' @family gradepath-accessors
#' @export
coef.gp_fit <- function(object, ...) {
  get_grades(object, ...)
}

#' Report-card table from a gradepath fit as a data frame
#'
#' `as.data.frame()` method for `gp_fit`: returns the report-card rows -- the
#' underlying data frame of the `gp_report_card` slot
#' (`get_report_card(fit)$table`) -- with row names reset to the default sequence.
#' Provided because `as.data.frame()` is the base generic an R user reaches for first
#' (`GP-DEC-13-A`, Chapter 13 patterns; Chapter 14 ux). Read-only; recomputes
#' nothing. Grade labels are integers in `{1, ..., n}` and carry no
#' ranking-superiority statement of any kind.
#'
#' @param x A `gp_fit` object.
#' @param ... Unused; accepted for S3 method extensibility. No arguments are
#'   forwarded; supplying any has no effect.
#'
#' @return A `data.frame`: the report-card table (one row per unit), row names reset
#'   to the default integer sequence.
#'
#' @examples
#' # Instant: the report-card rows of the bundled tiny fit as a data frame.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' head(as.data.frame(fit))   # the report-card rows
#'
#' @seealso [get_report_card()], [coef.gp_fit()]
#' @family gradepath-accessors
#' @export
as.data.frame.gp_fit <- function(x, ...) {
  tab <- get_report_card(x)$table
  rownames(tab) <- NULL
  tab
}
