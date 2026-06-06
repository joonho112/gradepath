# =============================================================================
# krw-r2.R -- KRW rsquared.do port
# -----------------------------------------------------------------------------
# This file ports the second-moment recipe in the KRW companion script
# `code/rsquared.do`. It is intentionally separate from `gp_r2()`, which remains
# a lightweight diagnostic variance-share helper for generic grade fits.
# =============================================================================

.gp_krw_r2_grade_col <- function(ranking, grade_col = NULL) {
  if (!is.data.frame(ranking)) {
    .gradepath_abort("`ranking` must be a data frame.")
  }
  if (!is.null(grade_col)) {
    grade_col <- .gradepath_validate_scalar_character(grade_col, "grade_col")
    if (!grade_col %in% names(ranking)) {
      .gradepath_abort("`ranking` is missing grade column `%s`.", grade_col)
    }
    return(grade_col)
  }

  candidates <- c("grades_lamb0.25", "grades_lamb025")
  hit <- candidates[candidates %in% names(ranking)]
  if (length(hit) == 0L) {
    .gradepath_abort(
      "`ranking` must contain `grades_lamb0.25`, `grades_lamb025`, or an explicit `grade_col`."
    )
  }
  hit[[1L]]
}

.gp_krw_r2_prepare <- function(theta, ranking, grade_col = NULL) {
  if (!is.data.frame(theta)) {
    .gradepath_abort("`theta` must be a data frame.")
  }
  if (!is.data.frame(ranking)) {
    .gradepath_abort("`ranking` must be a data frame.")
  }
  theta_required <- c("log_dif", "log_dif_se")
  ranking_required <- c("obs_idx", "post_mean_beta", "pmean_sq")
  theta_miss <- setdiff(theta_required, names(theta))
  ranking_miss <- setdiff(ranking_required, names(ranking))
  if (length(theta_miss) > 0L) {
    .gradepath_abort("`theta` is missing column(s): %s.", paste(theta_miss, collapse = ", "))
  }
  if (length(ranking_miss) > 0L) {
    .gradepath_abort("`ranking` is missing column(s): %s.", paste(ranking_miss, collapse = ", "))
  }

  if (!"obs_idx" %in% names(theta)) {
    theta$obs_idx <- seq_len(nrow(theta))
  }
  grade_col <- .gp_krw_r2_grade_col(ranking, grade_col)
  if (anyDuplicated(theta$obs_idx)) {
    .gradepath_abort("`theta$obs_idx` must be unique.")
  }
  if (anyDuplicated(ranking$obs_idx)) {
    .gradepath_abort("`ranking$obs_idx` must be unique.")
  }

  theta$obs_idx <- .gradepath_validate_numeric_vector(theta$obs_idx, "theta$obs_idx")
  ranking$obs_idx <- .gradepath_validate_numeric_vector(ranking$obs_idx, "ranking$obs_idx")
  idx <- match(theta$obs_idx, ranking$obs_idx)
  if (anyNA(idx)) {
    .gradepath_abort("Every `theta$obs_idx` value must be present in `ranking$obs_idx`.")
  }

  out <- data.frame(
    obs_idx = theta$obs_idx,
    log_dif = .gradepath_validate_numeric_vector(theta$log_dif, "theta$log_dif"),
    log_dif_se = .gradepath_validate_numeric_vector(theta$log_dif_se, "theta$log_dif_se"),
    grade = .gradepath_validate_numeric_vector(ranking[[grade_col]][idx], paste0("ranking$", grade_col)),
    post_mean_beta = .gradepath_validate_numeric_vector(
      ranking$post_mean_beta[idx],
      "ranking$post_mean_beta"
    ),
    pmean_sq = .gradepath_validate_numeric_vector(ranking$pmean_sq[idx], "ranking$pmean_sq"),
    stringsAsFactors = FALSE
  )
  if (any(out$log_dif_se <= 0)) {
    .gradepath_abort("`theta$log_dif_se` must be strictly positive.")
  }
  attr(out, "grade_col") <- grade_col
  out
}

.gp_krw_r2_components <- function(theta, ranking, grade_col = NULL) {
  dat <- .gp_krw_r2_prepare(theta, ranking, grade_col = grade_col)
  grade_col <- attr(dat, "grade_col", exact = TRUE)
  F <- nrow(dat)
  if (F < 2L) {
    .gradepath_abort("KRW R2 requires at least two rows.")
  }

  overall_var <- stats::var(dat$log_dif) - mean(dat$log_dif_se^2)
  if (!is.finite(overall_var) || overall_var <= 0) {
    .gradepath_abort("KRW R2 denominator is not positive.")
  }

  groups <- split(dat, dat$grade)
  grade_values <- suppressWarnings(as.numeric(names(groups)))
  ord <- order(grade_values)
  if (anyNA(grade_values)) {
    ord <- order(names(groups))
  }
  groups <- groups[ord]
  nfirms <- vapply(groups, nrow, integer(1))
  post_mean <- vapply(groups, function(x) mean(x$post_mean_beta), numeric(1))
  post_sum <- vapply(groups, function(x) sum(x$post_mean_beta), numeric(1))
  post_sq_sum <- vapply(groups, function(x) sum(x$post_mean_beta^2), numeric(1))
  pmean_sq_sum <- vapply(groups, function(x) sum(x$pmean_sq), numeric(1))

  tmp <- nfirms * post_mean
  summer <- 0
  ng <- length(groups)
  if (ng > 1L) {
    for (i in seq_len(ng - 1L)) {
      for (j in seq.int(i + 1L, ng)) {
        summer <- summer + tmp[[i]] * tmp[[j]]
      }
    }
  }

  e_tbar_sq <- (pmean_sq_sum + post_sum^2 - post_sq_sum) / (nfirms^2)
  one_firm <- nfirms == 1L
  if (any(one_firm)) {
    e_tbar_sq[one_firm] <- pmean_sq_sum[one_firm] / nfirms[one_firm]
  }
  summer2 <- sum(nfirms * (F - nfirms) * e_tbar_sq)

  between_var <- (summer2 / (F * (F - 1))) - (2 * summer / (F * (F - 1)))
  between_var <- between_var * ((F - 1) / F)
  if (!is.finite(between_var) || between_var < -sqrt(.Machine$double.eps)) {
    .gradepath_abort("KRW R2 between-grade variance is invalid.")
  }
  between_var <- max(0, between_var)

  data.frame(
    n = F,
    grade_count = ng,
    grade_col = grade_col,
    overall_variance = overall_var,
    overall_sd = sqrt(overall_var),
    between_variance = between_var,
    between_sd = sqrt(between_var),
    r2_proportion = between_var / overall_var,
    stringsAsFactors = FALSE
  )
}

#' KRW grade R-squared from the rsquared.do recipe
#'
#' `gp_krw_r2()` ports the Kline-Rose-Walters companion script
#' `code/rsquared.do`. It builds the bias-corrected denominator from the source
#' `theta_estimates_*` columns `log_dif` and `log_dif_se`, and the between-grade
#' numerator from the ranking-output posterior second moments
#' (`post_mean_beta`, `pmean_sq`). This is the paper-faithful producer for the
#' firm-level KRW race/gender R-squared rows; for the quick generic
#' variance-share diagnostic use [gp_r2()].
#'
#' @param theta Data frame shaped like `theta_estimates_race.csv` /
#'   `theta_estimates_gender.csv`: it must carry `log_dif` and (strictly positive)
#'   `log_dif_se`, and optionally `obs_idx`. When `obs_idx` is absent, row order
#'   defines it, matching `rsquared.do`.
#' @param ranking Data frame shaped like `ranking_results_*_*.csv`: it must carry
#'   `obs_idx`, `post_mean_beta`, `pmean_sq`, and a grade column. Joined to `theta`
#'   on `obs_idx`; every `theta` id must be present here.
#' @param grade_col Character scalar naming the grade column in `ranking`, or
#'   `NULL` (default) to auto-detect -- it accepts `grades_lamb0.25` and the
#'   Stata-style `grades_lamb025`.
#' @param scale Character scalar; `"percent"` (default) puts the returned `r2`
#'   column on a 0-100 scale, `"proportion"` on a 0-1 scale.
#'
#' @return A one-row data frame: `n` (units) and `grade_count`, the resolved
#'   `grade_col`, the variance components `overall_variance`/`overall_sd` and
#'   `between_variance`/`between_sd`, the `r2_proportion`, and `r2` on the
#'   requested `scale` (recorded in the `scale` column).
#'
#' @details The denominator subtracts the mean sampling variance from the raw
#'   estimate variance (the bias correction); the between-grade numerator applies
#'   the posterior second-moment correction from the ranking columns. The names
#'   pipeline uses different probability-scale quantities and is deferred to the
#'   `M1_NAMES` milestone, so this producer is firm-level race/gender only.
#'
#' @examples
#' # The KRW race/gender inputs are CSV products of the published solve, not part
#' # of the tiny example fixture; build the frame columns this recipe expects.
#' theta <- data.frame(
#'   log_dif    = c(-0.4, -0.1, 0.2, 0.5),
#'   log_dif_se = c(0.10, 0.12, 0.09, 0.11)
#' )
#' ranking <- data.frame(
#'   obs_idx         = 1:4,
#'   post_mean_beta  = c(-0.35, -0.08, 0.18, 0.46),
#'   pmean_sq        = c(0.130, 0.012, 0.038, 0.220),
#'   grades_lamb0.25 = c(1L, 2L, 2L, 3L)
#' )
#' gp_krw_r2(theta, ranking)
#'
#' @seealso [gp_r2()], [gp_frontier()], [gp_report_card()]
#' @family gradepath-frontier
#' @export
gp_krw_r2 <- function(theta,
                      ranking,
                      grade_col = NULL,
                      scale = c("percent", "proportion")) {
  scale <- match.arg(scale)
  out <- .gp_krw_r2_components(theta, ranking, grade_col = grade_col)
  out$r2 <- if (identical(scale, "percent")) {
    100 * out$r2_proportion
  } else {
    out$r2_proportion
  }
  out$scale <- scale
  out
}
