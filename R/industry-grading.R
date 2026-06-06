# =============================================================================
# industry-grading.R -- two-level industry grading/report cards
# -----------------------------------------------------------------------------
# Keeps the M2 industry grading path out of the M1 monolith.  The same solver
# and report-card assemblers are reused, but two-level Pi matrices enter through
# an explicit new-file dispatch:
#   industry_rfe : cleaned Pi_theta
#   btwn         : cleaned Pi_bar (= Pi_sbar_psi at industry representatives)
# =============================================================================

.gp_tlg_abort <- function(msg, ..., class = "gradepath_error") {
  args <- list(...)
  if (length(args) > 0L) msg <- do.call(sprintf, c(list(msg), args))
  .gradepath_abort(msg, class = class)
}

.gp_tlg_check_empty_dots <- function(caller, ...) {
  dots <- list(...)
  if (length(dots) == 0L) return(invisible(NULL))
  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("<unnamed>", length(dots))
  } else {
    dot_names[is.na(dot_names) | !nzchar(dot_names)] <- "<unnamed>"
  }
  .gp_tlg_abort(
    "`%s()` does not accept unused arguments: %s.",
    caller,
    paste(dot_names, collapse = ", ")
  )
}

.gp_tlg_validate_raw_pi <- function(x, name) {
  mat <- .gradepath_validate_numeric_matrix(x, name)
  if (!identical(nrow(mat), ncol(mat))) {
    .gp_tlg_abort("`%s` must be square.", name,
                  class = "gradepath_validation_error")
  }
  if (nrow(mat) < 2L) {
    .gp_tlg_abort("`%s` must contain at least two units.", name,
                  class = "gradepath_validation_error")
  }
  if (any(!is.finite(mat)) || any(mat < 0 | mat > 1)) {
    .gp_tlg_abort("`%s` must be finite probabilities in [0, 1].", name,
                  class = "gradepath_validation_error")
  }
  diag(mat) <- 0
  mat
}

.gp_tlg_ids <- function(ids, dimnames, n, prefix, name) {
  if (is.null(ids)) {
    ids <- dimnames[[1L]] %gp_or% dimnames[[2L]]
    if (is.null(ids)) ids <- paste0(prefix, seq_len(n))
  }
  ids <- .gradepath_validate_character_vector(
    as.character(ids),
    name,
    unique = TRUE
  )
  if (length(ids) != n) {
    .gp_tlg_abort("`%s` must have length %d.", name, n,
                  class = "gradepath_validation_error")
  }
  ids
}

.gp_tlg_lambda_grid <- function(lambda_grid, selected_lambda) {
  selected_lambda <- .gp_grade_validate_lambda(selected_lambda, "lambda")
  if (is.null(lambda_grid)) {
    return(sort(unique(c(.gp_control_required_lambda, selected_lambda, 1))))
  }
  lambda_grid
}

.gp_tlg_selected_status <- function(selected_grade) {
  .gp_producer_status_from_selected_grade(selected_grade)
}

.gp_tlg_combined_status <- function(statuses) {
  statuses <- as.character(statuses)
  if (length(statuses) == 0L || all(vapply(statuses, .gp_status_acceptance_ready,
                                          logical(1)))) {
    return("OK")
  }
  statuses[which(!vapply(statuses, .gp_status_acceptance_ready,
                         logical(1)))[1L]]
}

#' Wrap two-level Pi matrices as solver-ready pairwise objects
#'
#' `gp_twolevel_pairwise()` is the M2 pairwise dispatch for two-level industry
#' grading. It applies the same cleanup policy used by the one-level core
#' (antisymmetry, diagonal `0.5`, and the package zero floor) to raw two-level
#' outranking matrices, returning a bundle that carries both the firm-level
#' `industry_rfe` input (`Pi_theta`) and the between-industry `btwn` input
#' (`Pi_bar`) as validated `gp_pairwise` objects ready for [gp_twolevel_grade()].
#' It is solve-free.
#'
#' This is a structural M2 adapter, not an acceptance decision: it does not run the
#' fixture-promotion gate, set `PROMOTED` / `banded` status, compare industry DR,
#' tau, or R2 against paper targets, or update any M1 cache. M2 acceptance is the
#' job of [gp_m2_acceptance()] / [gp_m2_status()].
#'
#' @param Pi_theta Numeric square matrix; the firm-level two-level
#'   `Pr(theta_i > theta_j)` outranking probabilities. Must have at least two units
#'   and finite entries in `[0, 1]`; the diagonal is forced to zero on entry.
#' @param Pi_bar Numeric square matrix; the between-industry `Pi_bar` /
#'   `Pi_sbar_psi` outranking probabilities at the industry representatives. Same
#'   shape and value constraints as `Pi_theta`.
#' @param ids Optional character vector of firm ids labelling `Pi_theta`. When
#'   `NULL` (default), the matrix dimnames are used, or `unit_1`, `unit_2`, ... are
#'   generated. Must be unique and match `nrow(Pi_theta)`.
#' @param industry_levels Optional character vector of industry ids labelling
#'   `Pi_bar`. When `NULL` (default), the matrix dimnames are used, or `industry_1`,
#'   `industry_2`, ... are generated. Must be unique and match `nrow(Pi_bar)`.
#' @param control Optional [gp_control] object threaded onto both `gp_pairwise`
#'   wrappers. When `NULL` (default) a default `gp_control()` is used.
#'
#' @return A validated `gp_twolevel_pairwise_bundle` object (a list of class
#'   `c("gp_twolevel_pairwise_bundle", "list")`) with the public slots: \describe{
#'   \item{`ids`}{Character vector; the firm ids for `Pi_theta`.}
#'   \item{`industry_levels`}{Character vector; the industry ids for `Pi_bar`.}
#'   \item{`raw`}{Named list `Pi_theta` / `Pi_bar`: the inputs with only the
#'     diagonal zeroed (before the cleanup policy).}
#'   \item{`Pi_theta`, `Pi_bar`}{The cleaned firm- and industry-level matrices
#'     (antisymmetry, diagonal `0.5`, zero floor applied).}
#'   \item{`pairwise_theta`, `pairwise_bar`}{The matching `gp_pairwise` objects fed
#'     to [gp_twolevel_grade()].}
#'   \item{`cleanup`}{Named list recording the applied policy (`antisymmetry`,
#'     `diagonal`, `zero_floor`).}
#'   \item{`control`}{The validated [gp_control] used.}
#'   \item{`provenance`, `schema_version`, `warnings`}{Internal audit slots.}
#' }
#'
#' @examples
#' # A small synthetic strict-preference pair (the package's own test pattern):
#' # in strict(n, p) every off-diagonal entry is 0.95 above the diagonal and 0.05
#' # below, so the units are cleanly ordered. This step is solve-free.
#' strict <- function(n, p) {
#'   id <- paste0(p, seq_len(n)); m <- matrix(0, n, n, dimnames = list(id, id))
#'   for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) m[i, j] <- if (i < j) 0.95 else 0.05
#'   m
#' }
#' pw <- gp_twolevel_pairwise(strict(4, "f"), strict(3, "i"),
#'                            control = gp_control(backend = "highs"))
#' pw
#'
#' @seealso [gp_twolevel_grade()], [gp_twolevel_report_card()],
#'   [gp_m2_status()], [gp_control()]
#' @family gradepath-twolevel
#' @export
gp_twolevel_pairwise <- function(Pi_theta,
                                 Pi_bar,
                                 ids = NULL,
                                 industry_levels = NULL,
                                 control = NULL) {
  Pi_theta <- .gp_tlg_validate_raw_pi(Pi_theta, "Pi_theta")
  Pi_bar <- .gp_tlg_validate_raw_pi(Pi_bar, "Pi_bar")
  ids <- .gp_tlg_ids(ids, dimnames(Pi_theta), nrow(Pi_theta),
                     "unit_", "ids")
  industry_levels <- .gp_tlg_ids(
    industry_levels,
    dimnames(Pi_bar),
    nrow(Pi_bar),
    "industry_",
    "industry_levels"
  )
  control <- if (is.null(control)) gp_control() else validate_gp_control(control)

  dimnames(Pi_theta) <- list(ids, ids)
  dimnames(Pi_bar) <- list(industry_levels, industry_levels)
  Pi_theta_clean <- .gp_tl_clean_probability(Pi_theta)
  Pi_bar_clean <- .gp_tl_clean_probability(Pi_bar)
  dimnames(Pi_theta_clean) <- list(ids, ids)
  dimnames(Pi_bar_clean) <- list(industry_levels, industry_levels)

  pairwise_theta <- .gp_tl_new_pairwise(
    ids = ids,
    matrix = Pi_theta_clean,
    control = control,
    n_industries = length(industry_levels),
    provenance_step = "two-level-pairwise-theta-dispatch"
  )
  pairwise_bar <- .gp_tl_new_pairwise(
    ids = industry_levels,
    matrix = Pi_bar_clean,
    control = control,
    n_industries = length(industry_levels),
    provenance_step = "two-level-pairwise-bar-dispatch"
  )

  out <- list(
    ids = ids,
    industry_levels = industry_levels,
    raw = list(Pi_theta = Pi_theta, Pi_bar = Pi_bar),
    Pi_theta = Pi_theta_clean,
    Pi_bar = Pi_bar_clean,
    pairwise_theta = pairwise_theta,
    pairwise_bar = pairwise_bar,
    cleanup = list(
      antisymmetry = TRUE,
      diagonal = .gp_pairwise_diagonal,
      zero_floor = .gp_pairwise_zero_floor
    ),
    control = control,
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      producer = "gp_twolevel_pairwise",
      n_units = length(ids),
      n_industries = length(industry_levels),
      cleanup_order = "antisymmetry_diagonal_zero_floor"
    ),
    warnings = character(0)
  )
  validate_gp_twolevel_pairwise_bundle(
    structure(out, class = c("gp_twolevel_pairwise_bundle", "list"))
  )
}

validate_gp_twolevel_pairwise_bundle <- function(x) {
  if (!inherits(x, "gp_twolevel_pairwise_bundle")) {
    .gp_tlg_abort("Expected a gp_twolevel_pairwise_bundle object.",
                  class = "gradepath_validation_error")
  }
  req <- c("ids", "industry_levels", "raw", "Pi_theta", "Pi_bar",
           "pairwise_theta", "pairwise_bar", "cleanup", "control",
           "schema_version", "provenance", "warnings")
  if (any(!req %in% names(x))) {
    .gp_tlg_abort("`gp_twolevel_pairwise_bundle` is missing required fields.",
                  class = "gradepath_validation_error")
  }
  ids <- .gradepath_validate_character_vector(x$ids, "ids", unique = TRUE)
  industry_levels <- .gradepath_validate_character_vector(
    x$industry_levels,
    "industry_levels",
    unique = TRUE
  )
  for (nm in c("Pi_theta", "Pi_bar")) {
    mat <- .gradepath_validate_numeric_matrix(x[[nm]], nm)
    ids_nm <- if (identical(nm, "Pi_theta")) ids else industry_levels
    if (!identical(dim(mat), c(length(ids_nm), length(ids_nm)))) {
      .gp_tlg_abort("`%s` has inconsistent dimensions.", nm,
                    class = "gradepath_validation_error")
    }
  }
  if (!is.list(x$raw) || any(!c("Pi_theta", "Pi_bar") %in% names(x$raw))) {
    .gp_tlg_abort("`raw` must contain `Pi_theta` and `Pi_bar`.",
                  class = "gradepath_validation_error")
  }
  validate_gp_pairwise(x$pairwise_theta)
  validate_gp_pairwise(x$pairwise_bar)
  if (!isTRUE(all.equal(x$Pi_theta, x$pairwise_theta$matrix, tolerance = 0)) ||
      !isTRUE(all.equal(x$Pi_bar, x$pairwise_bar$matrix, tolerance = 0))) {
    .gp_tlg_abort("Cleaned matrices must equal their gp_pairwise wrappers.",
                  class = "gradepath_validation_error")
  }
  validate_gp_control(x$control)
  .gradepath_validate_scalar_character(
    x$schema_version,
    "schema_version",
    allowed = .gradepath_schema_version
  )
  .gradepath_validate_named_list(x$provenance, "provenance")
  .gradepath_validate_warning_vector(x$warnings, "warnings")
  x
}

.gp_tlg_resolve <- function(twolevel, pairwise_theta, pairwise_bar, posterior,
                            control) {
  pi <- NULL
  bundle <- NULL
  method <- "pairwise"

  if (!is.null(twolevel)) {
    if (inherits(twolevel, "gp_twolevel_grade")) {
      return(list(grade = validate_gp_twolevel_grade(twolevel)))
    }
    if (inherits(twolevel, "gp_twolevel_quadrature")) {
      twolevel <- validate_gp_twolevel_quadrature(twolevel)
      pi <- twolevel$pi
      posterior <- posterior %gp_or% twolevel$posterior
      pairwise_theta <- pairwise_theta %gp_or% twolevel$pairwise_theta
      pairwise_bar <- pairwise_bar %gp_or% twolevel$pairwise_bar
      method <- twolevel$method
    } else if (inherits(twolevel, "gp_twolevel_pi")) {
      pi <- validate_gp_twolevel_pi(twolevel)
      pairwise_theta <- pairwise_theta %gp_or% pi$pairwise_theta
      pairwise_bar <- pairwise_bar %gp_or% pi$pairwise_bar
      method <- "pi"
    } else if (inherits(twolevel, "gp_twolevel_pairwise_bundle")) {
      bundle <- validate_gp_twolevel_pairwise_bundle(twolevel)
      pairwise_theta <- pairwise_theta %gp_or% bundle$pairwise_theta
      pairwise_bar <- pairwise_bar %gp_or% bundle$pairwise_bar
      method <- "pairwise"
    } else if (is.list(twolevel) &&
               all(c("pairwise_theta", "pairwise_bar") %in% names(twolevel))) {
      pairwise_theta <- pairwise_theta %gp_or% twolevel$pairwise_theta
      pairwise_bar <- pairwise_bar %gp_or% twolevel$pairwise_bar
      posterior <- posterior %gp_or% twolevel$posterior
      method <- twolevel$method %gp_or% "pairwise"
    } else {
      .gp_tlg_abort(
        "`twolevel` must be a gp_twolevel_quadrature, gp_twolevel_pi, gp_twolevel_pairwise_bundle, or pairwise list.",
        class = "gradepath_validation_error"
      )
    }
  }

  if (is.null(pairwise_theta) || is.null(pairwise_bar)) {
    .gp_tlg_abort("`pairwise_theta` and `pairwise_bar` are required.",
                  class = "gradepath_validation_error")
  }
  pairwise_theta <- validate_gp_pairwise(pairwise_theta)
  pairwise_bar <- validate_gp_pairwise(pairwise_bar)
  if (!is.null(posterior)) posterior <- validate_gp_posterior(posterior)

  control <- if (is.null(control)) {
    pairwise_theta$control
  } else {
    validate_gp_control(control)
  }

  list(
    grade = NULL,
    pi = pi,
    bundle = bundle,
    posterior = posterior,
    pairwise_theta = pairwise_theta,
    pairwise_bar = pairwise_bar,
    control = control,
    method = method
  )
}

.gp_tlg_industry_posterior <- function(posterior, pi) {
  if (is.null(posterior) || is.null(pi)) return(NULL)
  posterior <- validate_gp_posterior(posterior)
  pi <- validate_gp_twolevel_pi(pi)
  tl <- posterior$metadata$two_level
  if (!is.list(tl) || is.null(tl$eta_posterior) ||
      is.null(tl$support_eta) || is.null(tl$sbar)) {
    return(NULL)
  }

  level <- as.numeric(posterior$metadata$interval_level %gp_or%
                        posterior$metadata$level %gp_or% 0.90)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    level <- 0.90
  }
  pct <- c((1 - level) / 2, 1 - ((1 - level) / 2))

  reps <- as.integer(pi$industry_representatives)
  ids <- as.character(pi$industry_levels)
  K <- length(ids)
  supp_eta <- as.numeric(tl$support_eta)
  eta_posterior <- as.matrix(tl$eta_posterior)
  sbar <- as.numeric(tl$sbar[reps])

  mean <- sd <- lower <- upper <- numeric(K)
  for (k in seq_len(K)) {
    values <- sbar[k] * supp_eta
    weights <- .gp_tl_normalize(eta_posterior[k, ], "eta_posterior")
    mean[k] <- sum(values * weights)
    second <- sum(values^2 * weights)
    sd[k] <- sqrt(max(0, second - mean[k]^2))
    lower[k] <- .gp_weighted_percentile(values, pct[1L], weights)
    upper[k] <- .gp_weighted_percentile(values, pct[2L], weights)
    lower[k] <- min(lower[k], mean[k])
    upper[k] <- max(upper[k], mean[k])
  }
  se <- pmax(sd, .Machine$double.eps)

  validate_gp_posterior(new_gp_posterior(
    estimate = mean,
    se = se,
    id = ids,
    label = ids,
    posterior_mean = mean,
    posterior_sd = sd,
    lower = lower,
    upper = upper,
    scale = "r",
    metadata = list(
      level = level,
      interval_level = level,
      component = "sbar_eta",
      two_level = list(
        characteristic = tl$characteristic %gp_or% NA_character_,
        industry_levels = ids,
        support_eta = supp_eta,
        eta_posterior = eta_posterior,
        sbar = sbar,
        percentile_convention = tl$percentile_convention %gp_or% NA_character_
      )
    )
  ))
}

.gp_tlg_solve_model <- function(pairwise, model, lambda_grid, lambda, control,
                                acceptance_mode) {
  path <- gp_grade_path(
    pairwise,
    lambda_grid = lambda_grid,
    control = control,
    selected_lambda = lambda,
    selection_rule = paste0(model, "_lambda_", format(lambda, trim = TRUE)),
    acceptance_mode = acceptance_mode
  )
  selected <- gp_select_grade(path, lambda = lambda)
  list(
    model = model,
    pairwise = pairwise,
    grade_path = path,
    selected_grade = selected,
    report_card = NULL,
    posterior = NULL,
    grade_count = selected$summary$grade_count,
    producer_status = .gp_tlg_selected_status(selected)
  )
}

#' Grade two-level industry matrices
#'
#' `gp_twolevel_grade()` routes the M2 two-level pairwise outputs into the package
#' grade solver without touching the one-level monolith: the `industry_rfe` route
#' grades the firm-level `Pi_theta`, and the `btwn` route grades the industry-level
#' `Pi_bar`. Feed it a bundle from [gp_twolevel_pairwise()] (or the lower-level
#' two-level carriers below). When a two-level posterior is supplied it also
#' assembles a firm `gp_report_card` for `industry_rfe` and an industry-level
#' `gp_report_card` for `btwn`, retrievable via [gp_twolevel_report_card()].
#'
#' The returned object is an M2 grading surface, not the M2 acceptance scorecard and
#' not a paper industry DR / tau / R2 reproduction. Use [gp_m2_acceptance()] to
#' classify recorded or supplied fixture-gate evidence into `PROMOTED`,
#' `APPROXIMATE_OK`, and the overall M2 status; see [gp_m2_status()].
#'
#' @param twolevel A `gp_twolevel_quadrature`, `gp_twolevel_pi`,
#'   `gp_twolevel_pairwise_bundle`, or a plain list carrying `pairwise_theta` and
#'   `pairwise_bar`. When `NULL` (default) the explicit `pairwise_theta` /
#'   `pairwise_bar` arguments must be supplied instead. Passing an already-graded
#'   `gp_twolevel_grade` returns it unchanged (no re-solve).
#' @param pairwise_theta,pairwise_bar Optional explicit `gp_pairwise` inputs for the
#'   firm-level and between-industry routes. When `NULL` (default) they are taken
#'   from `twolevel`; at least one source for each route must resolve.
#' @param posterior Optional two-level `gp_posterior` (e.g. from a two-level
#'   posterior builder). When `NULL` (default) the routes are graded structurally
#'   with no report cards. When supplied (and `build_report_cards = TRUE`) it backs
#'   the assembled report cards.
#' @param lambda Numeric in `[0, 1]`; the selected reporting penalty at which each
#'   route's grade is read off the path. Default `0.25`.
#' @param lambda_grid Optional numeric vector; the penalty grid solved for each
#'   route. When `NULL` (default) it is the package parity anchors plus `lambda`
#'   (`0.25`, `1`, and `lambda`), sorted and de-duplicated.
#' @param control Optional [gp_control] threaded to the solves. When `NULL`
#'   (default) the resolved `pairwise_theta`'s own control is used.
#' @param acceptance_mode Logical; the solver solution-quality / fallback policy
#'   forwarded to the internal [gp_grade_path()] (Gurobi backend only). `FALSE`
#'   (default) reports a gap/time-limit solve honestly with that status; `TRUE`
#'   adds optimization attempts and never relabels a non-optimal solve as optimal.
#' @param build_report_cards Logical; when `TRUE` (default) and a `posterior` is
#'   available, assemble the per-route `gp_report_card`s. When `FALSE` the routes
#'   are graded structurally only.
#' @param ... Reserved; an error is raised if any argument is passed.
#'
#' @return A validated `gp_twolevel_grade` object (a list of class
#'   `c("gp_twolevel_grade", "list")`) with the public slots: \describe{
#'   \item{`ids`, `industry_levels`}{Character vectors; the firm ids (`Pi_theta`
#'     route) and industry ids (`Pi_bar` route).}
#'   \item{`pairwise_theta`, `pairwise_bar`}{The `gp_pairwise` inputs actually
#'     graded for the two routes.}
#'   \item{`industry_rfe`}{Named list for the firm-level route: `grade_path`,
#'     `selected_grade` (at `lambda`), `grade_count`, `report_card` (or `NULL`),
#'     `posterior`, and `producer_status`.}
#'   \item{`btwn`}{The same payload for the between-industry route.}
#'   \item{`posterior`}{The supplied two-level `gp_posterior`, or `NULL`.}
#'   \item{`pi`, `pairwise`}{The two-level carriers (`gp_twolevel_pi`,
#'     `gp_twolevel_pairwise_bundle`) when supplied, else `NULL`.}
#'   \item{`selected_lambda`, `lambda_grid`}{The selected penalty and the solved
#'     grid.}
#'   \item{`method`}{Character tag for how the inputs were resolved (e.g.
#'     `"pairwise"`, `"pi"`).}
#'   \item{`producer_status`}{Combined honest solver status across the two routes
#'     (`"OK"` only when both routes are acceptance-ready).}
#'   \item{`control`}{The validated [gp_control] used.}
#'   \item{`provenance`, `schema_version`, `warnings`}{Internal audit slots.}
#' }
#'
#' @examples
#' # A small synthetic strict-preference pair (the package's own test pattern).
#' strict <- function(n, p) {
#'   id <- paste0(p, seq_len(n)); m <- matrix(0, n, n, dimnames = list(id, id))
#'   for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) m[i, j] <- if (i < j) 0.95 else 0.05
#'   m
#' }
#' pw <- gp_twolevel_pairwise(strict(4, "f"), strict(3, "i"),
#'                            control = gp_control(backend = "highs"))
#'
#' # The grade step solves a tiny integer program, so it needs a backend.
#' \donttest{
#' tlg <- gp_twolevel_grade(pw, control = gp_control(backend = "highs"),
#'                          lambda_grid = c(0.25, 1))
#' tlg
#' }
#'
#' @seealso [gp_twolevel_pairwise()], [gp_twolevel_report_card()],
#'   [gp_m2_acceptance()], [gp_grade_path()]
#' @family gradepath-twolevel
#' @export
gp_twolevel_grade <- function(twolevel = NULL,
                              pairwise_theta = NULL,
                              pairwise_bar = NULL,
                              posterior = NULL,
                              lambda = 0.25,
                              lambda_grid = NULL,
                              control = NULL,
                              acceptance_mode = FALSE,
                              build_report_cards = TRUE,
                              ...) {
  .gp_tlg_check_empty_dots("gp_twolevel_grade", ...)
  lambda <- .gp_grade_validate_lambda(lambda, "lambda")
  lambda_grid <- .gp_tlg_lambda_grid(lambda_grid, lambda)
  acceptance_mode <- .gradepath_validate_scalar_logical(
    acceptance_mode,
    "acceptance_mode"
  )
  build_report_cards <- .gradepath_validate_scalar_logical(
    build_report_cards,
    "build_report_cards"
  )
  resolved <- .gp_tlg_resolve(
    twolevel = twolevel,
    pairwise_theta = pairwise_theta,
    pairwise_bar = pairwise_bar,
    posterior = posterior,
    control = control
  )
  if (!is.null(resolved$grade)) return(resolved$grade)

  industry_rfe <- .gp_tlg_solve_model(
    resolved$pairwise_theta,
    model = "industry_rfe",
    lambda_grid = lambda_grid,
    lambda = lambda,
    control = resolved$control,
    acceptance_mode = acceptance_mode
  )
  btwn <- .gp_tlg_solve_model(
    resolved$pairwise_bar,
    model = "btwn",
    lambda_grid = lambda_grid,
    lambda = lambda,
    control = resolved$control,
    acceptance_mode = acceptance_mode
  )

  if (isTRUE(build_report_cards) && !is.null(resolved$posterior)) {
    industry_rfe$posterior <- resolved$posterior
    industry_rfe$report_card <- gp_report_card(
      posterior = resolved$posterior,
      selected_grade = industry_rfe$selected_grade,
      grade_path = industry_rfe$grade_path
    )

    btwn$posterior <- .gp_tlg_industry_posterior(
      resolved$posterior,
      resolved$pi
    )
    if (!is.null(btwn$posterior)) {
      btwn$report_card <- gp_report_card(
        posterior = btwn$posterior,
        selected_grade = btwn$selected_grade,
        grade_path = btwn$grade_path
      )
    }
  }

  statuses <- c(industry_rfe$producer_status, btwn$producer_status)
  out <- list(
    ids = resolved$pairwise_theta$ids,
    industry_levels = resolved$pairwise_bar$ids,
    pi = resolved$pi,
    pairwise = resolved$bundle,
    pairwise_theta = resolved$pairwise_theta,
    pairwise_bar = resolved$pairwise_bar,
    posterior = resolved$posterior,
    industry_rfe = industry_rfe,
    btwn = btwn,
    selected_lambda = lambda,
    lambda_grid = lambda_grid,
    control = resolved$control,
    method = resolved$method,
    producer_status = .gp_tlg_combined_status(statuses),
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      producer = "gp_twolevel_grade",
      n_units = length(resolved$pairwise_theta$ids),
      n_industries = length(resolved$pairwise_bar$ids),
      selected_lambda = lambda,
      route = "industry_rfe_Pi_theta__btwn_Pi_bar",
      m1_safe = TRUE
    ),
    warnings = unique(c(industry_rfe$grade_path$warnings,
                        btwn$grade_path$warnings))
  )
  validate_gp_twolevel_grade(
    structure(out, class = c("gp_twolevel_grade", "list"))
  )
}

validate_gp_twolevel_grade <- function(x) {
  if (!inherits(x, "gp_twolevel_grade")) {
    .gp_tlg_abort("Expected a gp_twolevel_grade object.",
                  class = "gradepath_validation_error")
  }
  req <- c("ids", "industry_levels", "pi", "pairwise", "pairwise_theta",
           "pairwise_bar", "posterior", "industry_rfe", "btwn",
           "selected_lambda", "lambda_grid", "control", "method",
           "producer_status", "schema_version", "provenance", "warnings")
  if (any(!req %in% names(x))) {
    .gp_tlg_abort("`gp_twolevel_grade` is missing required fields.",
                  class = "gradepath_validation_error")
  }
  ids <- .gradepath_validate_character_vector(x$ids, "ids", unique = TRUE)
  industry_levels <- .gradepath_validate_character_vector(
    x$industry_levels,
    "industry_levels",
    unique = TRUE
  )
  validate_gp_pairwise(x$pairwise_theta)
  validate_gp_pairwise(x$pairwise_bar)
  if (!identical(ids, x$pairwise_theta$ids) ||
      !identical(industry_levels, x$pairwise_bar$ids)) {
    .gp_tlg_abort("Two-level grade ids must match pairwise ids.",
                  class = "gradepath_validation_error")
  }
  if (!is.null(x$pi)) validate_gp_twolevel_pi(x$pi)
  if (!is.null(x$pairwise)) validate_gp_twolevel_pairwise_bundle(x$pairwise)
  if (!is.null(x$posterior)) validate_gp_posterior(x$posterior)
  lambda <- .gp_grade_validate_lambda(x$selected_lambda, "selected_lambda")
  lambda_grid <- .gp_validate_lambda_grid_core(x$lambda_grid, "lambda_grid")
  .gp_grade_exact_lambda_match(lambda_grid, lambda, "selected_lambda")
  validate_gp_control(x$control)

  validate_model <- function(model, expected_model, ids_expected) {
    if (!is.list(model) ||
        any(!c("model", "pairwise", "grade_path", "selected_grade",
               "report_card", "posterior", "grade_count",
               "producer_status") %in% names(model))) {
      .gp_tlg_abort("`%s` model payload is incomplete.", expected_model,
                    class = "gradepath_validation_error")
    }
    if (!identical(model$model, expected_model)) {
      .gp_tlg_abort("Model payload name mismatch.",
                    class = "gradepath_validation_error")
    }
    validate_gp_pairwise(model$pairwise)
    path <- validate_gp_grade_path(model$grade_path)
    selected <- validate_gp_grade_fit(model$selected_grade)
    if (!identical(path$ids, ids_expected) ||
        !identical(selected$ids, ids_expected)) {
      .gp_tlg_abort("`%s` ids do not match their pairwise route.",
                    expected_model, class = "gradepath_validation_error")
    }
    stored <- gp_select_grade(path, lambda = lambda)
    if (!identical(selected, stored)) {
      .gp_tlg_abort("`%s$selected_grade` must be the stored path member.",
                    expected_model, class = "gradepath_validation_error")
    }
    if (!identical(as.integer(model$grade_count),
                   as.integer(selected$summary$grade_count))) {
      .gp_tlg_abort("`%s$grade_count` must equal selected summary.",
                    expected_model, class = "gradepath_validation_error")
    }
    if (!is.null(model$posterior)) validate_gp_posterior(model$posterior)
    if (!is.null(model$report_card)) {
      card <- validate_gp_report_card(model$report_card)
      if (!setequal(card$ids, ids_expected)) {
        .gp_tlg_abort("`%s$report_card` ids do not match route ids.",
                      expected_model, class = "gradepath_validation_error")
      }
    }
    .gradepath_validate_scalar_character(
      model$producer_status,
      paste0(expected_model, "$producer_status")
    )
    model
  }

  industry_rfe <- validate_model(x$industry_rfe, "industry_rfe", ids)
  btwn <- validate_model(x$btwn, "btwn", industry_levels)
  .gradepath_validate_scalar_character(x$method, "method")
  .gradepath_validate_scalar_character(x$producer_status, "producer_status")
  .gradepath_validate_scalar_character(
    x$schema_version,
    "schema_version",
    allowed = .gradepath_schema_version
  )
  .gradepath_validate_named_list(x$provenance, "provenance")
  .gradepath_validate_warning_vector(x$warnings, "warnings")

  x$industry_rfe <- industry_rfe
  x$btwn <- btwn
  x
}

#' Return a two-level industry report card
#'
#' `gp_twolevel_report_card()` is the convenience selector for the report cards
#' assembled by [gp_twolevel_grade()]. Passing an existing `gp_twolevel_grade`
#' never re-solves; any other input is forwarded to [gp_twolevel_grade()] first.
#' Choose the firm-level (`industry_rfe`) or between-industry (`btwn`) card with
#' `model`.
#'
#' A report card exists only when the underlying two-level grade was built WITH
#' posterior information. This helper does not create new acceptance evidence, set
#' `PROMOTED` / `banded` status, or certify paper industry DR / tau / R2
#' reproduction; it simply returns the posterior-backed `gp_report_card` already
#' carried by the chosen route.
#'
#' @param twolevel A `gp_twolevel_grade` object, or any input accepted by
#'   [gp_twolevel_grade()] (which is then graded first).
#' @param model `"industry_rfe"` (default) for the firm-level two-level report card
#'   or `"btwn"` for the between-industry report card.
#' @param ... Arguments forwarded to [gp_twolevel_grade()] when `twolevel` is not
#'   already a graded `gp_twolevel_grade` (e.g. `posterior`, `control`, `lambda`).
#'
#' @return A validated `gp_report_card` for the chosen route: the per-unit
#'   grade-label table with posterior summaries.
#'
#' @note While the M2 industry surface is `PARTIAL_ACCEPTED`, no public
#'   input-to-industry-card builder is exported, so there is no public way to
#'   construct a posterior-backed two-level grade. Consequently a grade built from
#'   bare Pi matrices ([gp_twolevel_pairwise()] -> [gp_twolevel_grade()]) is
#'   structural only and carries no report card: calling this helper on it raises an
#'   informative error directing you to [gp_m2_status()]. The success path is shown
#'   in the package's two-level vignettes; see also [gp_m2_acceptance()].
#'
#' @examples
#' # The M2 industry surface is PARTIAL_ACCEPTED: there is no public posterior-backed
#' # two-level grade to build here, so the runnable example inspects the M2 contract.
#' gp_m2_status()
#' attr(gp_m2_status(), "m2_formal_status")   # "PARTIAL_ACCEPTED"
#'
#' \dontrun{
#' # By design this errors: a bare-Pi two-level grade is structural only and has no
#' # report card to return (no public builder is exported yet). See gp_m2_status().
#' strict <- function(n, p) {
#'   id <- paste0(p, seq_len(n)); m <- matrix(0, n, n, dimnames = list(id, id))
#'   for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) m[i, j] <- if (i < j) 0.95 else 0.05
#'   m
#' }
#' pw <- gp_twolevel_pairwise(strict(4, "f"), strict(3, "i"),
#'                            control = gp_control(backend = "highs"))
#' gp_twolevel_report_card(pw, control = gp_control(backend = "highs"))
#' }
#'
#' @seealso [gp_twolevel_grade()], [gp_twolevel_pairwise()],
#'   [gp_m2_status()], [gp_report_card()]
#' @family gradepath-twolevel
#' @export
gp_twolevel_report_card <- function(twolevel,
                                    model = c("industry_rfe", "btwn"),
                                    ...) {
  model <- match.arg(model)
  graded <- if (inherits(twolevel, "gp_twolevel_grade")) {
    validate_gp_twolevel_grade(twolevel)
  } else {
    gp_twolevel_grade(twolevel = twolevel, ...)
  }
  card <- graded[[model]]$report_card
  if (is.null(card)) {
    .gp_tlg_abort(
      paste0(
        "`gp_twolevel_report_card()` needs a posterior-backed two-level grade ",
        "(one carrying report cards). A grade built from bare Pi matrices ",
        "(`gp_twolevel_pairwise()` -> `gp_twolevel_grade()`) is structural only ",
        "and has no report card to return. The public input -> industry-card path ",
        "is not yet exported while the M2 industry surface is PARTIAL_ACCEPTED; ",
        "see `gp_m2_status()` and `?gp_twolevel_report_card`."
      ),
      class = "gradepath_validation_error"
    )
  }
  validate_gp_report_card(card)
}
