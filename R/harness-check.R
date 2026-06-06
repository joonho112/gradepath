# =============================================================================
# harness-check.R -- registry comparator and small harness adapters
# -----------------------------------------------------------------------------
# The harness never hard-codes paper values or tolerances. It reads gp_registry,
# expects the caller to pass replicated values already on the registry unit, and
# applies the single Chapter-16/17 absolute band comparator.
# =============================================================================

.gp_check_eps <- 1e-9

.gp_default_registry <- function() {
  ns <- asNamespace("gradepath")
  if (exists("gp_registry", envir = ns, inherits = FALSE)) {
    return(get("gp_registry", envir = ns, inherits = FALSE))
  }
  env <- new.env(parent = emptyenv())
  utils::data("gp_registry", package = "gradepath", envir = env)
  get("gp_registry", envir = env, inherits = FALSE)
}

.gp_registry_resolve <- function(registry = NULL) {
  if (is.null(registry)) {
    registry <- .gp_default_registry()
  }
  registry <- as.data.frame(registry, stringsAsFactors = FALSE)
  required <- c(
    "id",
    "paper_value",
    "unit",
    "tolerance",
    "class",
    "milestone",
    "quantity"
  )
  miss <- setdiff(required, names(registry))
  if (length(miss) > 0L) {
    .gradepath_abort("`registry` is missing column(s): %s.", paste(miss, collapse = ", "))
  }
  if (anyDuplicated(as.character(registry$id))) {
    .gradepath_abort("`registry$id` must be unique.")
  }
  registry
}

.gp_registry_row <- function(id, registry = NULL) {
  id <- .gradepath_validate_scalar_character(id, "id")
  registry <- .gp_registry_resolve(registry)
  row <- registry[as.character(registry$id) == id, , drop = FALSE]
  if (nrow(row) != 1L) {
    .gradepath_abort("Registry id not found: %s.", id)
  }
  row
}

.gp_check_scalar <- function(x, name = "replicated") {
  if (length(x) != 1L) {
    .gradepath_abort("`%s` must be a scalar.", name)
  }
  suppressWarnings(as.numeric(x))
}

.gp_check_status <- function(producer_status) {
  producer_status <- .gradepath_validate_scalar_character(
    producer_status,
    "producer_status"
  )
  .gp_status_normalize(producer_status)
}

#' Compare a replicated value against the gradepath registry
#'
#' `gp_check()` is the single registry comparator used by the replication harness.
#' It looks up one row of [gp_registry] by `id`, reads that row's published paper
#' value and absolute tolerance, compares your already-on-unit `replicated` value to
#' it, and returns a one-row verdict data frame. The harness embeds no tolerances of
#' its own -- every band comes from the registry -- and [gp_run_all()] is the vector
#' driver that calls this once per target.
#'
#' @param id Character scalar; a registry key present in `gp_registry$id`. An id not
#'   in the registry raises an error.
#' @param replicated Scalar replicated value, already converted to the registry unit
#'   (e.g. pass a percent for rows whose `unit` is `"percent"`). Defaults to
#'   `NA_real_`, which yields an `n/a` status (nothing to compare).
#' @param producer_status Character scalar; the runtime status of the producing
#'   stage, normalized internally. Defaults to `"OK"`. Any non-`OK` status short-
#'   circuits the numeric comparison and returns status `UNVERIFIED` (the producing
#'   stage was not acceptance-ready), regardless of `replicated`.
#' @param registry Registry data frame to look `id` up in. Defaults to `NULL`, which
#'   uses the bundled package data [gp_registry].
#' @param label Character scalar quantity label overriding `registry$quantity` in the
#'   returned row. Defaults to `NULL` (use the registry's own `quantity` text).
#'
#' @return A one-row `data.frame` describing the comparison. \describe{
#'   \item{`id`}{The registry id checked.}
#'   \item{`quantity`}{Human-readable target description (`label` if supplied, else
#'     the registry `quantity`).}
#'   \item{`paper`}{Numeric published value from the registry.}
#'   \item{`replicated`}{Numeric replicated value compared (`NA` when the producer
#'     status is non-`OK`).}
#'   \item{`delta`}{`replicated - paper` (`NA` when not compared).}
#'   \item{`tol`}{Absolute PASS tolerance from the registry.}
#'   \item{`unit`, `class`, `milestone`}{The registry row's unit, tolerance class,
#'     and milestone, passed through.}
#'   \item{`status`}{One of `PASS`, `FAIL`, `n/a` (missing value), `no-tol` (missing
#'     tolerance), or `UNVERIFIED` (non-`OK` producer status).}
#'   \item{`reason`}{Short reason string when not `PASS` (e.g. `outside_tolerance`,
#'     `missing_value`, `missing_tolerance`); `NA` on `PASS`.}
#'   \item{`producer_status`}{The normalized producer status used.}
#' }
#'
#' @details
#' The comparison is the single Chapter 16/17 absolute-band rule: `PASS` iff
#' `abs(replicated - paper) <= tol` (with a tiny `1e-9` numerical slack), `FAIL`
#' otherwise. This is a verification check against published values; it makes no
#' ranking-superiority statement of any kind.
#'
#' @examples
#' # Instant: compare a replicated count against a known registry target. The
#' # firm-count target (id "scale_n_firms_graded") is published as 97, tolerance 0.
#' gp_check("scale_n_firms_graded", replicated = 97)   # status PASS
#' gp_check("scale_n_firms_graded", replicated = 96)   # status FAIL (outside tol 0)
#'
#' @seealso [gp_run_all()], [gp_validate_targets()], [gp_registry]
#' @family gradepath-harness
#' @export
gp_check <- function(id,
                     replicated = NA_real_,
                     producer_status = "OK",
                     registry = NULL,
                     label = NULL) {
  row <- .gp_registry_row(id, registry = registry)
  status <- .gp_check_status(producer_status)
  paper <- suppressWarnings(as.numeric(row$paper_value[[1L]]))
  tol <- suppressWarnings(as.numeric(row$tolerance[[1L]]))
  class <- as.character(row$class[[1L]])
  unit <- as.character(row$unit[[1L]])
  milestone <- as.character(row$milestone[[1L]])
  quantity <- if (is.null(label)) as.character(row$quantity[[1L]]) else {
    .gradepath_validate_scalar_character(label, "label")
  }

  if (!identical(status, "OK")) {
    return(data.frame(
      id = as.character(row$id[[1L]]),
      quantity = quantity,
      paper = paper,
      replicated = NA_real_,
      delta = NA_real_,
      tol = tol,
      unit = unit,
      class = class,
      milestone = milestone,
      status = "UNVERIFIED",
      reason = status,
      producer_status = status,
      stringsAsFactors = FALSE
    ))
  }

  replicated <- .gp_check_scalar(replicated, "replicated")
  delta <- replicated - paper
  verdict <- if (is.na(paper) || is.na(replicated)) {
    "n/a"
  } else if (is.na(tol)) {
    "no-tol"
  } else if (abs(delta) <= tol + .gp_check_eps) {
    "PASS"
  } else {
    "FAIL"
  }
  reason <- switch(
    verdict,
    PASS = NA_character_,
    FAIL = "outside_tolerance",
    `n/a` = "missing_value",
    `no-tol` = "missing_tolerance",
    NA_character_
  )

  data.frame(
    id = as.character(row$id[[1L]]),
    quantity = quantity,
    paper = paper,
    replicated = replicated,
    delta = delta,
    tol = tol,
    unit = unit,
    class = class,
    milestone = milestone,
    status = verdict,
    reason = reason,
    producer_status = status,
    stringsAsFactors = FALSE
  )
}

#' Grade-level share of bias-corrected estimate variance
#'
#' `gp_r2()` reports the percent (or proportion) of the bias-corrected estimate
#' variance that is explained by the selected grade means -- a quick
#' between-grade variance-share diagnostic for any grade fit. It is the
#' lightweight generic helper, not the paper-faithful KRW recipe; for the
#' Kline-Rose-Walters race/gender R-squared target definition use [gp_krw_r2()].
#'
#' @param x Optional source of the three vectors below. A `gp_fit` (its
#'   per-unit `estimate`, `se`, and `grade` are read), or a data frame carrying
#'   `estimate`, `se`, and `grade` columns. `NULL` (default) means supply the
#'   vectors explicitly via `estimate`/`se`/`grade`.
#' @param estimate,se,grade Numeric vectors used when `x` is `NULL`: the
#'   bias-corrected point estimates, their strictly positive standard errors, and
#'   the integer grade labels (relabeled to be contiguous internally). All three
#'   must be the same length. Ignored when `x` is supplied.
#' @param scale Character scalar; `"percent"` (default) returns the share on a
#'   0-100 scale, `"proportion"` returns it on a 0-1 scale.
#'
#' @return Numeric scalar: the between-grade variance share on the requested
#'   scale, or `NA_real_` when it is undefined (fewer than two units, or a
#'   non-positive bias-corrected total variance).
#'
#' @examples
#' # Variance share explained by the selected grades of the bundled fit.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' gp_r2(fit)                       # percent scale
#' gp_r2(fit, scale = "proportion") # 0-1 scale
#'
#' @seealso [gp_krw_r2()], [gp_frontier()], [gp_report_card()]
#' @family gradepath-frontier
#' @export
gp_r2 <- function(x = NULL,
                  estimate = NULL,
                  se = NULL,
                  grade = NULL,
                  scale = c("percent", "proportion")) {
  scale <- match.arg(scale)
  if (inherits(x, "gp_fit")) {
    tab <- as.data.frame(validate_gp_fit(x))
    estimate <- tab$estimate
    se <- tab$se
    grade <- tab$grade
  } else if (is.data.frame(x)) {
    miss <- setdiff(c("estimate", "se", "grade"), names(x))
    if (length(miss) > 0L) {
      .gradepath_abort("`x` is missing column(s): %s.", paste(miss, collapse = ", "))
    }
    estimate <- x$estimate
    se <- x$se
    grade <- x$grade
  } else if (!is.null(x)) {
    .gradepath_abort("`x` must be a `gp_fit`, data frame, or NULL.")
  }

  estimate <- .gradepath_validate_numeric_vector(estimate, "estimate")
  se <- .gradepath_validate_numeric_vector(se, "se")
  if (length(se) != length(estimate)) {
    .gradepath_abort("`se` must have the same length as `estimate`.")
  }
  if (any(se <= 0)) {
    .gradepath_abort("`se` must be strictly positive.")
  }
  grade <- monotone_relabel(grade)
  if (length(grade) != length(estimate)) {
    .gradepath_abort("`grade` must have the same length as `estimate`.")
  }
  if (length(estimate) < 2L) {
    return(NA_real_)
  }

  total_var <- stats::var(estimate) - mean(se^2)
  if (!is.finite(total_var) || total_var <= 0) {
    return(NA_real_)
  }
  group_means <- as.numeric(tapply(estimate, grade, mean))
  group_n <- as.numeric(tabulate(grade, nbins = max(grade)))
  overall <- stats::weighted.mean(group_means, group_n)
  between_var <- sum(group_n * (group_means - overall)^2) / (sum(group_n) - 1)
  r2 <- between_var / total_var
  if (identical(scale, "percent")) {
    100 * r2
  } else {
    r2
  }
}

.gp_read_headerless_csv <- function(path) {
  path <- .gradepath_validate_scalar_character(path, "path")
  if (!file.exists(path)) {
    .gradepath_abort("File not found: %s.", path)
  }
  utils::read.csv(path, header = FALSE, stringsAsFactors = FALSE)
}

.gp_asin2p <- function(x) {
  x <- .gradepath_validate_numeric_vector(x, "x")
  sin(x)^2
}
