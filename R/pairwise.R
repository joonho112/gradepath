# =============================================================================
# pairwise.R  --  pairwise outranking probabilities (gp_pairwise)
# -----------------------------------------------------------------------------
# Chapter 4. The pairwise outranking probability is the posterior probability
# that unit i's latent exceeds unit j's:
#
#   pi_ij = Pr(theta_i > theta_j | data) = sum_m W[i, m] * F_j(theta_i[m]^-)
#
# where theta_i = reporting_support[, i] is unit i's r-grid mapped onto ITS OWN
# reporting (theta) scale (sec-ch4-backtransform: per-unit, because the precision
# adjustment depends on the unit's own original_s), W[i, ] is unit i's posterior
# mass on the r-scale (the Chapter-8 W, never back-transformed -- only the axis
# is), and F_j(t^-) is unit j's posterior CDF STRICTLY below t. The CDF form
# (O(M) per pair, eq-ch4-cdf) is the v1 keeper:
#
#   cumulative_j <- c(0, cumsum(W[j, ]))
#   below_index  <- findInterval(theta_i, theta_j, left.open = TRUE)
#   pi_ij        <- sum(W[i, ] * cumulative_j[below_index + 1L])
#
# `left.open = TRUE` + the leading 0 implement the STRICT inequality
# theta_i > theta_j (ties contribute 0).
#
# Cleanup contract (sec-ch4-cleanup, frozen invariant 6), applied in order
# antisymmetry -> diagonal -> zero_floor:
#   1. Antisymmetry: compute the upper triangle, fill the lower by
#      pi_ji = 1 - pi_ij (not recomputed).
#   2. Diagonal = 0.5 exactly (a neutral self-comparison).
#   3. Zero-floor = 1e-7: ALL off-diagonal values below 1e-7 are raised to 1e-7
#      (a one-sided floor, `pmax(pmin(., 1), 1e-7)`), and values are capped at 1.
#      This catches more than exact 0 -- any tiny positive pi below 1e-7 is also
#      floored. Exact 1 remains exact 1 (no symmetric cap at 1 - 1e-7). This
#      floor-all rule follows the design (frozen invariant 6)
#      and is a deliberate v1 keeper. (KRW generate_grades.py:55 floors only
#      EXACT zeros; the package deliberately follows the design's floor-all
#      rule -- they coincide except for tiny-positive entries below 1e-7.)
#
# Orientation (sec-ch4-orient / Chapter 8): W is J x M (units = rows),
# reporting_support is M x J (column per unit), so pi is J x J with row i = how
# unit i outranks every j.
# =============================================================================


#' Pairwise outranking probability matrix
#'
#' `gp_pairwise()` returns the `J x J` pairwise outranking matrix whose entry
#' `pi[i, j]` is the posterior probability `Pr(theta_i > theta_j)`. Given a
#' materialized `gp_fit` it is a read-only accessor for the stored matrix; given
#' the output of `gp_posterior_weights()` (or a bare `J x M` weight matrix plus a
#' reporting axis) it computes the matrix by exact posterior-CDF integration and
#' applies the frozen cleanup contract. The result is the outranking structure the
#' grade engine ([gp_grade_path()], [gp_grade()]) scores.
#'
#' @param object A `gp_fit` (read-only accessor path), the list returned by
#'   `gp_posterior_weights()` carrying `$W` (`J x M`) and `$reporting_support`
#'   (`M x J`), or a bare `J x M` `W` matrix. A bare matrix is accepted only when
#'   `reporting_support` is supplied or a shared axis is otherwise intended.
#'   Default `NULL` requires the legacy `weights =` alias (see `...`).
#' @param reporting_support Optional `M x J` per-unit reporting axis; overrides
#'   `object$reporting_support` on the weight-list path. `NULL` (default) takes the
#'   axis from `object`. Ignored when `object` is a `gp_fit`.
#' @param ids Optional length-`J` identifiers. `NULL` (default) takes them from
#'   `object$ids`, the row names of `W`, or `seq_len(J)`. Ignored when `object` is
#'   a `gp_fit`.
#' @param control Optional [gp_control] recorded on the result. `NULL` (default)
#'   uses [gp_control()] on the weight-list path. Ignored when `object` is a
#'   `gp_fit`.
#' @param assumption Character scalar; the pairwise independence assumption
#'   recorded in `source`. Defaults to `"one_level_independence"` (the one-level
#'   analytic path). Ignored when `object` is a `gp_fit`.
#' @param ... Unused, except for the legacy named `weights =` alias for `object`;
#'   any other argument raises an error.
#' @param power Integer scalar; the matrix-power option. Only the KRW / public
#'   path `0L` (default) is implemented on the M1 surface; any other value errors.
#'
#' @return A validated [gp_pairwise] object (a list of class
#'   `c("gp_pairwise", "list")`) with the public slots: \describe{
#'   \item{`ids`}{Character vector of length `J`; the canonical unit-id order, and
#'     the row/column names of `matrix`.}
#'   \item{`matrix`}{Numeric `J x J` matrix of outranking probabilities in
#'     `[0, 1]`; `matrix[i, j] = Pr(theta_i > theta_j)`, diagonal `0.5`.}
#'   \item{`power`}{Integer; the matrix-power option (`0L`).}
#'   \item{`cleanup`}{Named list (`antisymmetry`, `diagonal`, `zero_floor`); the
#'     applied cleanup contract -- antisymmetry on, diagonal `0.5`, one-sided
#'     `1e-7` floor.}
#'   \item{`source`}{Named list (`stage`, `rule`, `assumption`); how the matrix
#'     was produced and under which independence assumption.}
#'   \item{`control`}{The [gp_control] recorded on the result.}
#'   \item{`warnings`, `schema_version`, `provenance`}{Internal audit slots.}
#' }
#'
#' @details
#' On the compute path the entry is the per-unit CDF integral
#' `pi[i, j] = sum_m W[i, m] * F_j(theta_i[m]^-)`, where `F_j` is unit `j`'s
#' posterior CDF strictly below the argument (ties contribute `0`). Cleanup is
#' applied in the frozen order antisymmetry (`pi[j, i] = 1 - pi[i, j]`) ->
#' diagonal set to `0.5` -> zero-floor: every off-diagonal value below `1e-7` is
#' raised to `1e-7` (a one-sided floor; exact `1` stays exact `1`), so final
#' antisymmetry is exact up to the floor tolerance. See KRW (2024, Section 4) and
#' the pairwise method vignette.
#'
#' @examples
#' # Read-only accessor: pull the stored pairwise matrix from the bundled fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' pw <- gp_pairwise(fit)
#' pw                              # a 24 x 24 gp_pairwise
#' pw$matrix[1:4, 1:4]             # the outranking probabilities, diagonal 0.5
#'
#' @seealso [gp_grade_path()], [gp_grade()], [krw_report_card()], [gp_report_card()]
#' @family gradepath-pairwise
#' @export
gp_pairwise <- function(object = NULL,
                        reporting_support = NULL,
                        ids = NULL,
                        control = NULL,
                        assumption = "one_level_independence",
                        ...,
                        power = 0L) {
  power <- .gradepath_validate_integerish(power, "power", min = 0L)
  if (!identical(power, 0L)) {
    .gradepath_abort("`gp_pairwise()` currently supports only `power = 0L`.")
  }

  dots <- list(...)
  if (is.null(object) && "weights" %in% names(dots)) {
    object <- dots$weights
    dots$weights <- NULL
  }
  if (is.null(object)) {
    .gradepath_abort("`object` is required.")
  }
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    if (is.null(dot_names)) {
      dot_names <- rep("<unnamed>", length(dots))
    } else {
      dot_names[is.na(dot_names) | !nzchar(dot_names)] <- "<unnamed>"
    }
    .gradepath_abort(
      "`gp_pairwise()` does not accept unused argument(s): %s.",
      paste(dot_names, collapse = ", ")
    )
  }

  if (inherits(object, "gp_fit")) {
    if (!is.null(reporting_support) || !is.null(ids) || !is.null(control) ||
        !identical(assumption, "one_level_independence")) {
      .gradepath_abort(
        "`gp_pairwise()` with a `gp_fit` input does not accept stage-wise arguments."
      )
    }
    return(validate_gp_pairwise(get_pairwise(validate_gp_fit(object))))
  }

  if (inherits(object, "gp_posterior")) {
    .gradepath_abort(
      paste0(
        "`gp_pairwise()` cannot rebuild pairwise probabilities from a ",
        "`gp_posterior` alone in this build because posterior weights `W` are ",
        "not stored on the posterior object. Pass a `gp_fit` or the list ",
        "returned by `gp_posterior_weights()`."
      )
    )
  }

  .gp_pairwise_from_weights(
    weights = object,
    reporting_support = reporting_support,
    ids = ids,
    control = control,
    assumption = assumption
  )
}

#' Pairwise outranking probability matrix from posterior weights
#' @keywords internal
#' @noRd
.gp_pairwise_from_weights <- function(weights, reporting_support = NULL,
                                      ids = NULL, control = NULL,
                                      assumption = "one_level_independence") {
  ## ---- extract W + per-unit reporting axis ---------------------------------
  W <- NULL; rs <- reporting_support
  if (is.matrix(weights)) {
    W <- weights
  } else if (is.list(weights) && !is.null(weights$W)) {
    W <- weights$W
    if (is.null(rs)) rs <- weights$reporting_support
    if (is.null(ids) && !is.null(weights$ids)) ids <- weights$ids
    if (is.null(rs) && !is.null(weights$support)) {
      rs <- matrix(as.numeric(weights$support), nrow = length(weights$support),
                   ncol = nrow(W))                 # shared-axis fallback
    }
  }
  if (is.null(W) || !is.matrix(W) || !is.numeric(W)) {
    .gradepath_abort(
      "`weights` must carry a numeric J x M posterior weight matrix `W`.")
  }
  J <- nrow(W); M <- ncol(W)
  .gp_assert_w_row_stochastic(
    W,
    n_units = J,
    n_support = M,
    name = if (is.matrix(weights)) "weights" else "weights$W"
  )
  if (is.null(rs)) {
    .gradepath_abort(
      "`reporting_support` (M x J) is required (from gp_posterior_weights()).")
  }
  rs <- as.matrix(rs)
  if (!is.numeric(rs) || any(!is.finite(rs))) {
    .gradepath_abort("`reporting_support` must be a finite numeric matrix.")
  }
  if (nrow(rs) != M || ncol(rs) != J) {
    .gradepath_abort("`reporting_support` must be M x J = %d x %d (got %d x %d).",
                     M, J, nrow(rs), ncol(rs))
  }

  ## ---- ids -----------------------------------------------------------------
  if (is.null(ids)) {
    rn <- rownames(W)
    ids <- if (!is.null(rn)) rn else as.character(seq_len(J))
  }
  ids <- as.character(ids)
  if (length(ids) != J) {
    .gradepath_abort("`ids` (%d) must have length J = %d.", length(ids), J)
  }

  ## ---- raw upper-triangle pi via the per-unit CDF integral -----------------
  raw <- .gp_pairwise_raw_matrix(W, rs)

  ## ---- cleanup: antisymmetry -> diagonal 0.5 -> zero-floor 1e-7 ------------
  pi <- .gp_pairwise_cleanup_matrix(raw)
  dimnames(pi) <- list(ids, ids)

  ## ---- assemble + validate -------------------------------------------------
  ctrl <- if (is.null(control)) gp_control() else control
  obj <- new_gp_pairwise(
    ids = ids,
    matrix = pi,
    power = 0L,
    cleanup = list(antisymmetry = TRUE,
                   diagonal = .gp_pairwise_diagonal,
                   zero_floor = .gp_pairwise_zero_floor),
    source = list(stage = "posterior",
                  rule = "outer_product",
                  assumption = assumption),
    control = ctrl,
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      step = "pairwise-pi", n_units = J, n_support = M,
      rule = "outer_product",
      cleanup_order = "antisymmetry_diagonal_zero_floor"),
    warnings = character(0))
  validate_gp_pairwise(obj)
}


# ---------------------------------------------------------------------------
# .gp_pairwise_raw_probability  --  single ordered pair via the CDF integral
# ---------------------------------------------------------------------------
# pi(i > j) = sum_m W[i, m] * F_j(theta_i[m]^-), with F_j the strictly-below
# posterior CDF on unit j's axis (v1 keeper, eq-ch4-cdf). theta_i / theta_j are
# the two units' reporting axes; weights_i / weights_j their r-scale posterior
# masses.
#' @keywords internal
#' @noRd
.gp_pairwise_raw_probability <- function(theta_i, weights_i, theta_j, weights_j) {
  cumulative_j <- c(0, cumsum(weights_j))                       # F_j step function
  below_index  <- findInterval(theta_i, theta_j, left.open = TRUE)
  sum(weights_i * cumulative_j[below_index + 1L])              # eq-ch4-cdf
}


# ---------------------------------------------------------------------------
# .gp_pairwise_raw_matrix  --  upper-triangle pi (no cleanups)
# ---------------------------------------------------------------------------
# Fills only the strict upper triangle (j > i) via the per-unit CDF integral;
# the diagonal and lower triangle are left as 0 for the cleanup step to set
# (diagonal -> 0.5, lower -> antisymmetry 1 - pi_ij). theta axes are the columns
# of `reporting_support` (M x J); weights are the rows of `W` (J x M).
#' @keywords internal
#' @noRd
.gp_pairwise_raw_matrix <- function(W, reporting_support) {
  J <- nrow(W)
  raw <- matrix(0, nrow = J, ncol = J)
  if (J >= 2L) {
    for (i in seq_len(J - 1L)) {
      theta_i <- reporting_support[, i]
      w_i <- W[i, ]
      for (j in seq.int(i + 1L, J)) {
        raw[i, j] <- .gp_pairwise_raw_probability(
          theta_i, w_i, reporting_support[, j], W[j, ])
      }
    }
  }
  raw
}


# ---------------------------------------------------------------------------
# .gp_pairwise_cleanup_matrix  --  antisymmetry -> diagonal 0.5 -> zero-floor
# ---------------------------------------------------------------------------
# Input: a J x J matrix with only the strict upper triangle populated (raw
# pi_ij for j > i). Output: the cleaned pairwise matrix. Order matters and is
# recorded in provenance as antisymmetry_diagonal_zero_floor (sec-ch4-cleanup):
#   1. antisymmetry: pi_ji = 1 - pi_ij for the lower triangle;
#   2. diagonal: 0.5;
#   3. zero-floor: ALL off-diagonal values below 1e-7 are raised to 1e-7 and
#      values are capped at 1 (`pmax(pmin(pi[off], 1), 1e-7)`) -- this floors
#      tiny positives, not only exact 0; exact 1 is kept exact rather than
#      capped at 1 - 1e-7.
#' @keywords internal
#' @noRd
.gp_pairwise_cleanup_matrix <- function(raw) {
  J <- nrow(raw)
  pi <- matrix(0, nrow = J, ncol = J)
  ## 1. antisymmetry: upper from raw, lower = 1 - upper.
  for (i in seq_len(J)) {
    for (j in seq_len(J)) {
      if (i < j) {
        pi[i, j] <- raw[i, j]
      } else if (i > j) {
        pi[i, j] <- 1 - raw[j, i]
      }
    }
  }
  ## 2. diagonal.
  diag(pi) <- .gp_pairwise_diagonal
  ## 3. one-sided zero-floor (off-diagonal only); exact ones remain exact.
  off <- row(pi) != col(pi)
  z <- .gp_pairwise_zero_floor
  pi[off] <- pmax(pmin(pi[off], 1), z)
  pi
}
