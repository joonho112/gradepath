# gradepath v2 output and composite typed objects.
#
# Output-layer + native EB-adjacent + composite S3 objects:
#   gp_report_card       (ADAPT of v1 gradepath_report_card)
#   gp_prior             (native; eb_prior-shaped     — mirrors ebrecipe::eb_prior)
#   gp_posterior         (native; eb_posterior-shaped — mirrors ebrecipe::eb_posterior)
#   gp_precision_fit     (native KRW beta-GMM standardization fit shell)
#   gp_fit               (composite; ADAPT->COMPOSE of v1 gradepath_fit)
#
# Schema source: the design's typed-objects spec + the
# real INSTALLED ebrecipe objects (v0.5.0) whose slot surface the native
# gp_prior / gp_posterior MIRROR so get_prior() / get_posterior() round-trip to
# ebrecipe's eb_* accessors.
#
# Reuses (do NOT redefine) — REAL names confirmed against the
# built R/ (utils-validate.R, control.R, class-objects-decision.R):
#   .gradepath_schema_version                                (utils-validate.R)
#   .gradepath_validate_list_class, .gradepath_validate_named_fields,
#   .gradepath_validate_named_list, .gradepath_validate_character_vector,
#   .gradepath_validate_numeric_vector, .gradepath_validate_scalar_numeric,
#   .gradepath_validate_scalar_character, .gradepath_validate_data_frame,
#   .gradepath_new_provenance, .gradepath_abort                (utils-validate.R)
#   validate_gp_pairwise, validate_gp_grade_fit,
#   validate_gp_grade_path                                   (class-objects-decision.R)
#   validate_gp_control                                      (control.R)
#
# NOTE on helper names: the design brief listed aspirational short helper
# names (.gradepath_validate_chr_vec / _int_vec / _df / _positive_dbl_vec /
# _scalar_in, and a .gradepath_match_lambda). Those do NOT exist in the built
# utils-validate.R; the real surface is the longer .gradepath_validate_*
# family above. Integer-vector and scalar-range checks are composed from those
# primitives here (no new shared helper is introduced — this is output
# objects, not the validator library).
#
# === ebrecipe slot surfaces MIRRORED (confirmed against installed v0.5.0,
#     built via eb_input(theta_hat=, s=, unit_id=); source:
#     ebrecipe-R-package-v2/R/class-constructors.R) ===
#   eb_estimates : theta_hat, s, unit_id, n, covariates, source, description,
#                  standardized, original_theta_hat, original_s,
#                  standardization_model, hyperparameters              (12 fields)
#   eb_prior     : support, density, mean, scale, diagnostics, metadata (6)
#   eb_posterior : estimate, se, id, label, posterior_mean, posterior_sd,
#                  lower, upper, scale, metadata                    (10 fields)
#                  (NB: eb_shrink leaves lower/upper = NA; gradepath's native
#                   posterior populates them.)
#   gp_prior   mirrors eb_prior  EXACTLY (same 6 fields, same order, scale="r").
#   gp_posterior mirrors eb_posterior EXACTLY (same 10 fields, same order).
#   -> a gp_prior/gp_posterior passes ebrecipe's validate_eb_prior /
#      validate_eb_posterior and is read by ebrecipe accessors unchanged.

# ---- field-name vectors (canonical slot order) -----------------------------
# (utils-validate.R deliberately leaves these output-object vectors here,
# co-located here next to their constructors.)

#' @keywords internal
.gp_report_card_fields <- c(
  "ids", "table", "selected_lambda", "grades", "control",
  "schema_version", "provenance", "warnings"
)

#' @keywords internal
.gp_report_card_table_columns <- c(
  "id",
  "label",
  "grade",
  "sort_rank",
  "selected_lambda",
  "posterior_mean",
  "lower",
  "upper",
  "estimate",
  "se"
)

#' @keywords internal
.gp_prior_fields <- c(
  "support", "density", "mean", "scale", "diagnostics", "metadata"
)

#' @keywords internal
.gp_posterior_fields <- c(
  "estimate", "se", "id", "label", "posterior_mean", "posterior_sd",
  "lower", "upper", "scale", "metadata"
)

#' @keywords internal
.gp_precision_fit_fields <- c(
  "parameters", "moments", "diagnostics", "scale",
  "schema_version", "provenance", "warnings"
)

#' @keywords internal
.gp_fit_fields <- c(
  "ids", "estimates", "prior", "posterior", "precision_fit",
  "pairwise", "grade_path", "selected_grade", "report_card", "control",
  "schema_version", "provenance", "warnings"
)

# ---- small internal helpers (local) ----------------------------------------

#' NULL-coalescing helper (base R has none until 4.4)
#' @keywords internal
`%gp_or%` <- function(a, b) if (is.null(a)) b else a

#' Integer-vector validator composed from the built numeric-vector primitive
#'
#' The brief named a `.gradepath_validate_int_vec`; it does not exist. We compose
#' the contract (numeric, finite, length n, integer-valued, >= min) from the real
#' `.gradepath_validate_numeric_vector`, returning an `integer` vector. Internal
#' to the output-object layer.
#' @keywords internal
.gp_validate_int_vec <- function(x, n = NULL, min = 1L, what = "value") {
  x <- .gradepath_validate_numeric_vector(x, what)
  if (!is.null(n) && length(x) != n) {
    .gradepath_abort("`%s` must have length %d.", what, as.integer(n))
  }
  if (any(abs(x - round(x)) > 1e-8) || any(x < min)) {
    .gradepath_abort("`%s` must be integers >= %d.", what, as.integer(min))
  }
  as.integer(round(x))
}

#' Read the per-unit point-estimate vector (length J) from an eb_estimates
#'
#' The installed ebrecipe eb_estimates stores point estimates as `theta_hat`.
#' The `estimate` fallback is retained only for older stand-ins and cached tests.
#' @keywords internal
.gp_estimates_theta <- function(estimates) {
  estimates$theta_hat %gp_or% estimates$estimate
}

#' Read unit ids from an eb_estimates-like container
#' @keywords internal
.gp_estimates_id <- function(estimates) {
  estimates$unit_id %gp_or% estimates$id %gp_or% estimates$ids
}

#' Read unit standard errors from an eb_estimates-like container
#' @keywords internal
.gp_estimates_se <- function(estimates) {
  estimates$s %gp_or% estimates$se
}

#' Read `original_s` (reporting-scale SEs) from an eb_estimates container
#'
#' The installed ebrecipe eb_estimates has a top-level `original_s` field. The
#' `metadata$original_s` fallback is retained for older cached stand-ins.
#' @keywords internal
.gp_estimates_original_s <- function(estimates) {
  estimates$original_s %gp_or% estimates$metadata$original_s
}

#' Read report-card posterior columns, preferring reporting-scale summaries
#' @keywords internal
.gp_posterior_report_card_values <- function(posterior) {
  reporting <- posterior$metadata$reporting
  if (is.list(reporting) &&
      all(c("posterior_mean", "lower", "upper") %in% names(reporting))) {
    return(list(
      posterior_mean = as.numeric(reporting$posterior_mean),
      lower = as.numeric(reporting$lower),
      upper = as.numeric(reporting$upper)
    ))
  }

  list(
    posterior_mean = as.numeric(posterior$posterior_mean),
    lower = as.numeric(posterior$lower),
    upper = as.numeric(posterior$upper)
  )
}

# ===========================================================================
# gp_report_card  (ADAPT of v1 gradepath_report_card)
# ===========================================================================
#
# Slot list per the design brief (the v1-ported surface):
#   ids, table (df[J]), selected_lambda, grades (int[J]), control,
#   schema_version, provenance, warnings.
# NOTE: the design's report-card *table* lists a richer `data/selection_rule/
# path_context/backend` surface; the brief pins the leaner ids/table/grades
# surface, which is what validate_gp_fit's report<->grade / report<->posterior
# cross-checks consume. See SUMMARY.md "report_card schema divergence".

#' @keywords internal
new_gp_report_card <- function(ids, table, selected_lambda, grades, control,
                               schema_version = .gradepath_schema_version,
                               provenance = list(),
                               warnings = character()) {
  structure(
    list(ids = ids, table = table, selected_lambda = selected_lambda,
         grades = grades, control = control,
         schema_version = schema_version, provenance = provenance,
         warnings = warnings),
    class = c("gp_report_card", "list")
  )
}

#' @keywords internal
validate_gp_report_card <- function(x) {
  .gradepath_validate_list_class(x, "gp_report_card")
  .gradepath_validate_named_fields(x, .gp_report_card_fields, "gp_report_card")

  ids <- .gradepath_validate_character_vector(x$ids, "gp_report_card$ids",
                                              unique = TRUE)
  J <- length(ids)

  tab <- .gradepath_validate_data_frame(x$table, "gp_report_card$table")
  if (nrow(tab) != J) {
    .gradepath_abort(
      "`gp_report_card$table` must be a J-row data.frame (one row per id).")
  }

  missing_columns <- setdiff(.gp_report_card_table_columns, names(tab))
  if (length(missing_columns) > 0L) {
    .gradepath_abort(
      "`gp_report_card$table` must include columns: %s.",
      paste(.gp_report_card_table_columns, collapse = ", ")
    )
  }

  grades <- .gp_validate_int_vec(x$grades, n = J, min = 1L,
                                 what = "gp_report_card$grades")
  unique_grades <- sort(unique(grades))
  if (!identical(unique_grades, seq_len(length(unique_grades)))) {
    .gradepath_abort(
      "`gp_report_card$grades` must be normalized to contiguous integers starting at 1."
    )
  }

  selected_lambda <- .gradepath_validate_scalar_numeric(
    x$selected_lambda, "gp_report_card$selected_lambda",
    lower = 0, upper = 1, include_lower = TRUE, include_upper = TRUE
  )

  # ids consistency: report-card ids are in the table's endpoint-sorted row
  # order, so the slot and table must match in that order.
  if (!identical(as.character(tab$id), ids)) {
    .gradepath_abort(
      "`gp_report_card$table$id` must equal `gp_report_card$ids`.")
  }
  if (!isTRUE(all.equal(as.numeric(tab$grade), as.numeric(grades)))) {
    .gradepath_abort(
      "`gp_report_card$table$grade` must equal `gp_report_card$grades`.")
  }
  if (!isTRUE(all.equal(as.numeric(tab$selected_lambda), rep(selected_lambda, J)))) {
    .gradepath_abort(
      "`gp_report_card$table$selected_lambda` must equal `gp_report_card$selected_lambda`."
    )
  }

  tab$sort_rank <- .gp_validate_int_vec(
    tab$sort_rank,
    n = J,
    min = 1L,
    what = "gp_report_card$table$sort_rank"
  )
  if (!identical(tab$sort_rank, seq_len(J))) {
    .gradepath_abort(
      "`gp_report_card$table$sort_rank` must be the endpoint row order 1, ..., J."
    )
  }

  for (name in c("posterior_mean", "estimate", "se")) {
    tab[[name]] <- .gradepath_validate_numeric_vector(
      tab[[name]],
      paste0("gp_report_card$table$", name)
    )
  }
  if (any(tab$se <= 0)) {
    .gradepath_abort("`gp_report_card$table$se` must be positive.")
  }

  # lower <= mean <= upper guard, only when the interval columns are present.
  lo <- tab[["lower"]]
  up <- tab[["upper"]]
  if (!is.numeric(lo) || !is.numeric(up) || length(lo) != J || length(up) != J) {
    .gradepath_abort("`gp_report_card$table$lower`/`upper` must be numeric of length J.")
  }
  if (any(!is.na(lo) & !is.finite(lo)) || any(!is.na(up) & !is.finite(up))) {
    .gradepath_abort("`gp_report_card$table$lower`/`upper` must be finite or NA.")
  }
  band <- is.finite(lo) & is.finite(tab$posterior_mean) & is.finite(up)
  if (any(band & (lo > tab$posterior_mean + 1e-8)) ||
      any(band & (tab$posterior_mean > up + 1e-8))) {
    .gradepath_abort(
      "`gp_report_card$table` must satisfy lower <= posterior_mean <= upper.")
  }

  # control slot must be a valid gp_control (the sibling decision-object
  # validators all do this; QA-1.5 found this validator was the lone omission).
  validate_gp_control(x$control)

  x
}

# ===========================================================================
# gp_prior  (NATIVE; eb_prior-shaped — MIRRORS ebrecipe::eb_prior exactly)
# ===========================================================================
#
# Full internals of the native two-level deconvolution are owned by the
# estimation layer. Here we build the validated eb_prior-shaped SHELL that gp_fit
# composes and that get_prior() hands to ebrecipe accessors. The two-level
# xi (x) eta fields (group_fx==1) ride in `diagnostics`/`metadata` (advisory),
# so adding them never breaks this schema. See SUMMARY.md "full-vs-shell".
#
# CONTRACT NOTE -- `density` is probability MASSES, not a PDF. Despite the field
# NAME (inherited from the eb_prior shape), `gp_prior$density` holds probability
# MASSES on the support grid: it is non-negative and SUMS TO 1 (validated below),
# and `mean == sum(support * density)` (a discrete expectation over masses, not a
# trapezoidal integral of a density). This follows the Matlab `g_xi/sum(g_xi)`
# convention. This is DISTINCT from `ebrecipe::eb_prior$density`, which is a
# per-point PDF (density values, integrated against grid widths). gradepath's W
# recipe (posterior-weights.R) reads these masses directly via
# `.gp_prior_grid_mass()`; on a uniform grid the constant spacing cancels under
# normalization so masses == density / sum(density).

#' @keywords internal
new_gp_prior <- function(support, density, mean, scale = "r",
                         diagnostics = list(), metadata = list()) {
  structure(
    list(support = support, density = density, mean = mean, scale = scale,
         diagnostics = diagnostics, metadata = metadata),
    class = c("gp_prior", "list")
  )
}

#' @keywords internal
validate_gp_prior <- function(x) {
  .gradepath_validate_list_class(x, "gp_prior")
  .gradepath_validate_named_fields(x, .gp_prior_fields, "gp_prior")

  sup <- .gradepath_validate_numeric_vector(x$support, "gp_prior$support")
  M <- length(sup)
  if (M < 2L || any(diff(sup) <= 0)) {
    .gradepath_abort("`gp_prior$support` must be strictly increasing.")
  }
  dens <- .gradepath_validate_numeric_vector(x$density, "gp_prior$density")
  if (length(dens) != M) {
    .gradepath_abort(
      "`gp_prior$density` must match the length of `gp_prior$support`.")
  }
  if (any(dens < 0)) {
    .gradepath_abort("`gp_prior$density` must be non-negative.")
  }
  dens_sum <- sum(dens)
  if (!is.finite(dens_sum) || abs(dens_sum - 1) > 1e-6) {
    .gradepath_abort("`gp_prior$density` must sum to 1 within tolerance 1e-6.")
  }
  prior_mean <- sum(sup * dens)
  mean <- .gradepath_validate_scalar_numeric(x$mean, "gp_prior$mean")
  if (abs(mean - prior_mean) > 1e-6) {
    .gradepath_abort("`gp_prior$mean` must equal sum(support * density).")
  }
  .gradepath_validate_scalar_character(x$scale, "gp_prior$scale",
                                       allowed = "r")
  .gradepath_validate_named_list(x$diagnostics, "gp_prior$diagnostics")
  .gradepath_validate_named_list(x$metadata, "gp_prior$metadata")
  x
}

# ===========================================================================
# gp_posterior  (NATIVE; eb_posterior-shaped — MIRRORS ebrecipe::eb_posterior)
# ===========================================================================
#
# Full internals of the native two-level Monte-Carlo posterior are owned by
# the estimation layer. Here we build the validated eb_posterior-shaped
# SHELL gp_fit composes; the five Pi matrices (group_fx==1) ride in `metadata`.
# Unlike ebrecipe's eb_shrink (CIs left NA), gradepath's native posterior DOES
# populate lower/upper, so this validator enforces lower<=pm<=upper on the
# finite band while still TOLERATING NA CIs (eb_posterior parity).

#' @keywords internal
new_gp_posterior <- function(estimate, se, id, label,
                             posterior_mean, posterior_sd, lower, upper,
                             scale = "r", metadata = list()) {
  structure(
    list(estimate = estimate, se = se, id = id, label = label,
         posterior_mean = posterior_mean, posterior_sd = posterior_sd,
         lower = lower, upper = upper, scale = scale, metadata = metadata),
    class = c("gp_posterior", "list")
  )
}

#' @keywords internal
validate_gp_posterior <- function(x) {
  .gradepath_validate_list_class(x, "gp_posterior")
  .gradepath_validate_named_fields(x, .gp_posterior_fields, "gp_posterior")

  id <- .gradepath_validate_character_vector(x$id, "gp_posterior$id")
  J <- length(id)

  pm <- .gradepath_validate_numeric_vector(x$posterior_mean,
                                           "gp_posterior$posterior_mean")
  if (length(pm) != J) {
    .gradepath_abort("`gp_posterior$posterior_mean` must have length J.")
  }
  est <- .gradepath_validate_numeric_vector(x$estimate, "gp_posterior$estimate")
  se <- .gradepath_validate_numeric_vector(x$se, "gp_posterior$se")
  psd <- .gradepath_validate_numeric_vector(x$posterior_sd,
                                            "gp_posterior$posterior_sd")
  for (nm in c("estimate", "se", "posterior_sd")) {
    if (length(x[[nm]]) != J) {
      .gradepath_abort("`gp_posterior$%s` must be numeric of length J.", nm)
    }
  }
  if (any(se <= 0)) {
    .gradepath_abort("`gp_posterior$se` must be positive.")
  }
  if (any(psd < 0)) {
    .gradepath_abort("`gp_posterior$posterior_sd` must be non-negative.")
  }
  lo <- x$lower; up <- x$upper
  if (!is.numeric(lo) || !is.numeric(up) || length(lo) != J || length(up) != J) {
    .gradepath_abort("`gp_posterior$lower`/`upper` must be numeric of length J.")
  }
  if (any(!is.na(lo) & !is.finite(lo)) || any(!is.na(up) & !is.finite(up))) {
    .gradepath_abort("`gp_posterior$lower`/`upper` must be finite or NA.")
  }
  # CIs MAY be NA (mirrors eb_posterior's permissible NA CIs); guard only the
  # finite band so the ordering invariant holds where it is populated.
  band <- is.finite(lo) & is.finite(up) & is.finite(pm)
  if (any(band & (lo > pm + 1e-8))) {
    .gradepath_abort("`gp_posterior`: lower must be <= posterior_mean.")
  }
  if (any(band & (pm > up + 1e-8))) {
    .gradepath_abort("`gp_posterior`: posterior_mean must be <= upper.")
  }
  .gradepath_validate_scalar_character(x$scale, "gp_posterior$scale",
                                       allowed = "r")
  .gradepath_validate_named_list(x$metadata, "gp_posterior$metadata")
  x
}

# ===========================================================================
# gp_precision_fit  (NATIVE KRW beta-GMM standardization fit)
# ===========================================================================
#
# Re-owned because ebrecipe's eb_standardize is a different
# estimand: conditional-mean NLLS, not KRW's 4-moment 2-step GMM. Full GMM
# internals (sandwich, V_m, penalty grid) are owned by the estimation layer.
# Here we build the validated SHELL gp_fit's krw_gmm cross-check
# consumes: `parameters` MUST carry `model_form`, `beta`, `mu`.

#' @keywords internal
new_gp_precision_fit <- function(parameters, moments = list(),
                                 diagnostics = list(), scale = "r",
                                 schema_version = .gradepath_schema_version,
                                 provenance = list(),
                                 warnings = character()) {
  structure(
    list(parameters = parameters, moments = moments,
         diagnostics = diagnostics, scale = scale,
         schema_version = schema_version, provenance = provenance,
         warnings = warnings),
    class = c("gp_precision_fit", "list")
  )
}

#' @keywords internal
validate_gp_precision_fit <- function(x) {
  .gradepath_validate_list_class(x, "gp_precision_fit")
  .gradepath_validate_named_fields(x, .gp_precision_fit_fields,
                                   "gp_precision_fit")
  .gradepath_validate_named_list(x$parameters, "gp_precision_fit$parameters")
  # The three load-bearing parameters consumed by the W/reporting-support seam.
  req <- c("model_form", "beta", "mu")
  miss <- setdiff(req, names(x$parameters))
  if (length(miss) > 0L) {
    .gradepath_abort("`gp_precision_fit$parameters` is missing: %s.",
                     paste(miss, collapse = ", "))
  }
  .gradepath_validate_scalar_character(
    x$parameters$model_form, "gp_precision_fit$parameters$model_form")
  .gradepath_validate_scalar_numeric(
    x$parameters$beta, "gp_precision_fit$parameters$beta")
  .gradepath_validate_scalar_numeric(
    x$parameters$mu, "gp_precision_fit$parameters$mu")
  .gradepath_validate_named_list(x$moments, "gp_precision_fit$moments")
  .gradepath_validate_named_list(x$diagnostics, "gp_precision_fit$diagnostics")
  x
}

# ===========================================================================
# gp_fit  (composite; ADAPT->COMPOSE of v1 gradepath_fit)
# ===========================================================================
#
# Slot list (13 slots, EXACT order):
#   ids, estimates (eb_estimates), prior (gp_prior), posterior (gp_posterior),
#   precision_fit (gp_precision_fit|NULL), pairwise, grade_path,
#   selected_grade, report_card, control, schema_version, provenance, warnings.

#' @keywords internal
new_gp_fit <- function(ids, estimates, prior, posterior, precision_fit,
                       pairwise, grade_path, selected_grade, report_card,
                       control, schema_version = .gradepath_schema_version,
                       provenance = list(), warnings = character()) {
  structure(
    list(ids = ids, estimates = estimates, prior = prior,
         posterior = posterior, precision_fit = precision_fit,
         pairwise = pairwise, grade_path = grade_path,
         selected_grade = selected_grade, report_card = report_card,
         control = control, schema_version = schema_version,
         provenance = provenance, warnings = warnings),
    class = c("gp_fit", "list")
  )
}

#' @keywords internal
validate_gp_fit <- function(x) {
  .gradepath_validate_list_class(x, "gp_fit")
  .gradepath_validate_named_fields(x, .gp_fit_fields, "gp_fit")

  # --- composition: stage-1 input validated through the ebrecipe seam;
  #     prior/posterior native (validated by gradepath) ----------------------
  if (!inherits(x$estimates, "eb_estimates")) {
    .gradepath_abort("`gp_fit$estimates` must be an `ebrecipe::eb_estimates`.")
  }
  # Delegate stage-1 validation to ebrecipe through the seam. The validator is
  # an internal ebrecipe symbol, so namespace introspection stays confined to
  # R/seam-ebrecipe.R.
  if (requireNamespace("ebrecipe", quietly = TRUE)) {
    .gp_eb_validate_estimates(x$estimates)
  }
  if (!inherits(x$prior, "gp_prior")) {
    .gradepath_abort("`gp_fit$prior` must be a gradepath `gp_prior`.")
  }
  if (!inherits(x$posterior, "gp_posterior")) {
    .gradepath_abort("`gp_fit$posterior` must be a gradepath `gp_posterior`.")
  }
  validate_gp_prior(x$prior)
  posterior <- validate_gp_posterior(x$posterior)

  # --- social-choice slots: gradepath's own validators ---------------------
  pairwise       <- validate_gp_pairwise(x$pairwise)
  grade_path     <- validate_gp_grade_path(x$grade_path)
  selected_grade <- validate_gp_grade_fit(x$selected_grade)
  report_card    <- validate_gp_report_card(x$report_card)

  # ======================================================================
  # CROSS-SLOT CHECK (a): one shared canonical `ids` order across
  #   gp_fit$ids / pairwise / grade_path / selected_grade / report_card /
  #   native gp_posterior.
  # ======================================================================
  ids <- .gradepath_validate_character_vector(x$ids, "gp_fit$ids",
                                              unique = TRUE)
  if (!identical(pairwise$ids, ids)) {
    .gradepath_abort("`gp_fit`: pairwise ids must match canonical `ids`.")
  }
  if (!identical(grade_path$ids, ids)) {
    .gradepath_abort("`gp_fit`: grade_path ids must match canonical `ids`.")
  }
  if (!identical(selected_grade$ids, ids)) {
    .gradepath_abort("`gp_fit`: selected_grade ids must match canonical `ids`.")
  }
  if (!setequal(report_card$ids, ids)) {
    .gradepath_abort("`gp_fit`: report_card ids must contain the canonical `ids`.")
  }
  if (!identical(as.character(posterior$id), ids)) {
    .gradepath_abort("`gp_fit`: posterior id must match canonical `ids`.")
  }
  estimate_ids <- .gp_estimates_id(x$estimates)
  if (is.null(estimate_ids) || !identical(as.character(estimate_ids), ids)) {
    .gradepath_abort("`gp_fit`: estimates ids must match canonical `ids`.")
  }

  # ======================================================================
  # CROSS-SLOT CHECK (b): selected_grade IS the path member at
  #   grade_path$selection$selected_lambda. (Grid membership is asserted by
  #   validate_gp_grade_path; here we locate the member by floating match —
  #   there is no .gradepath_match_lambda helper in the built library.)
  # ======================================================================
  sel_lambda <- grade_path$selection$selected_lambda
  idx        <- which(abs(grade_path$lambda_grid - sel_lambda) < 1e-8)
  if (length(idx) != 1L) {
    .gradepath_abort(
      "`gp_fit`: selected_lambda must match exactly one grade_path$lambda_grid point.")
  }
  member          <- grade_path$fits[[idx]]
  member_grades   <- member$assignment$grade
  selected_grades <- selected_grade$assignment$grade
  if (!identical(selected_grade, member)) {
    .gradepath_abort(
      "`gp_fit`: selected_grade must be the stored path member at selected_lambda.")
  }
  if (!isTRUE(all.equal(as.numeric(member_grades),
                        as.numeric(selected_grades)))) {
    .gradepath_abort(
      "`gp_fit`: selected_grade must equal the path member at selected_lambda.")
  }

  # ======================================================================
  # CROSS-SLOT CHECK (c): report_card agrees with selected_grade on GRADES,
  #   and with the native gp_posterior on posterior_mean / lower / upper.
  # ======================================================================
  report_selected <- selected_grades[match(report_card$ids, selected_grade$ids)]
  if (anyNA(report_selected) ||
      !isTRUE(all.equal(as.numeric(report_card$grades),
                        as.numeric(report_selected)))) {
    .gradepath_abort(
      "`gp_fit`: report_card grades must agree with selected_grade.")
  }
  # report_card$selected_lambda must agree with the path selection.
  if (!isTRUE(all.equal(report_card$selected_lambda, sel_lambda))) {
    .gradepath_abort(
      "`gp_fit`: report_card selected_lambda must agree with the grade_path selection.")
  }
  # report card table columns reconciled against the NATIVE gp_posterior.
  tab <- report_card$table
  posterior_index <- match(report_card$ids, as.character(posterior$id))
  if (anyNA(posterior_index)) {
    .gradepath_abort(
      "`gp_fit`: report_card ids must align to the native gp_posterior.")
  }
  posterior_values <- .gp_posterior_report_card_values(posterior)
  if (!isTRUE(all.equal(as.numeric(tab$posterior_mean),
                        as.numeric(posterior_values$posterior_mean[posterior_index])))) {
    .gradepath_abort(
      "`gp_fit`: report_card posterior_mean must agree with the native gp_posterior.")
  }
  if (!isTRUE(all.equal(as.numeric(tab$lower),
                        as.numeric(posterior_values$lower[posterior_index])))) {
    .gradepath_abort(
      "`gp_fit`: report_card lower must agree with the native gp_posterior.")
  }
  if (!isTRUE(all.equal(as.numeric(tab$upper),
                        as.numeric(posterior_values$upper[posterior_index])))) {
    .gradepath_abort(
      "`gp_fit`: report_card upper must agree with the native gp_posterior.")
  }
  estimate_index <- match(report_card$ids, ids)
  estimate_values <- .gp_estimates_theta(x$estimates)
  se_values <- .gp_estimates_se(x$estimates)
  if (length(estimate_values) != length(ids) || length(se_values) != length(ids)) {
    .gradepath_abort(
      "`gp_fit`: estimates theta/se must match canonical `ids` length.")
  }
  if (!isTRUE(all.equal(as.numeric(tab$estimate),
                        as.numeric(estimate_values[estimate_index])))) {
    .gradepath_abort(
      "`gp_fit`: report_card estimate must agree with estimates.")
  }
  if (!isTRUE(all.equal(as.numeric(tab$se),
                        as.numeric(se_values[estimate_index])))) {
    .gradepath_abort(
      "`gp_fit`: report_card se must agree with estimates.")
  }

  # ======================================================================
  # CROSS-SLOT CHECK (d): precision_fit presence matches the standardization
  #   rule — NULL iff precision_rule == "none".
  # ======================================================================
  rule <- x$control$precision_rule
  if (identical(rule, "none")) {
    if (!is.null(x$precision_fit)) {
      .gradepath_abort(
        "`gp_fit`: precision_fit must be NULL when precision_rule == \"none\".")
    }
  } else {
    if (is.null(x$precision_fit)) {
      .gradepath_abort(
        "`gp_fit`: precision_fit must be non-NULL when precision_rule != \"none\".")
    }
    validate_gp_precision_fit(x$precision_fit)
  }

  # ======================================================================
  # CROSS-SLOT CHECK (e): under "krw_gmm", estimates$original_s is non-null
  #   and length J, and precision_fit$parameters has model_form/beta/mu.
  #   (NB: the design's pseudocode names estimates$theta_hat / estimates$original_s
  #    at top level; the INSTALLED eb_estimates uses `estimate` and has no
  #    `original_s` slot — original_s lives in metadata to stay ebrecipe-valid.
  #    Accessors resolve either home. See SUMMARY.md.)
  # ======================================================================
  if (identical(rule, "krw_gmm")) {
    orig_s <- .gp_estimates_original_s(x$estimates)
    if (is.null(orig_s)) {
      .gradepath_abort(
        "`gp_fit`: estimates$original_s must be non-null under precision_rule == \"krw_gmm\".")
    }
    J <- length(.gp_estimates_theta(x$estimates))
    if (length(orig_s) != J) {
      .gradepath_abort(
        "`gp_fit`: estimates$original_s must have length J (= number of units).")
    }
    pf_par <- x$precision_fit$parameters
    miss <- setdiff(c("model_form", "beta", "mu"), names(pf_par))
    if (length(miss) > 0L) {
      .gradepath_abort(
        paste0("`gp_fit`: precision_fit$parameters must contain model_form, ",
               "beta, mu under precision_rule == \"krw_gmm\" (missing: %s)."),
        paste(miss, collapse = ", "))
    }
  }

  x
}
