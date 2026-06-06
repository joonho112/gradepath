#' Load KRW's native beta-GMM input estimates (the Matlab lsqnonlin `data` arg)
#'
#' @description
#' The beta-GMM core (`gp_estimation_core`) reproduces KRW's published precision
#' exponent only when fed KRW's *actual* GMM input: the per-firm
#' `(theta_hat, s)` estimates KRW's Matlab `estimate_lsqnonlin.m` consumed,
#' shipped here as `inst/extdata/krw-gmm-input/theta_estimates_matlab_<char>.csv`
#' (a headerless `[industry_d, theta_hat, s]` matrix; 97 firms; the published
#' 19-SIC-group industry coding).
#'
#' NOTE (data provenance -- open question, flagged for the maintainer): this is
#' NOT the same numeric series as `ebrecipe::krw_firms` (the bundled
#' example). On `ebrecipe::krw_firms` the GMM degenerates to a spurious
#' large-beta optimum; on KRW's real input the faithful two-step GMM reproduces
#' KRW Table 3 (race beta = 0.5095 vs published 0.510; gender beta = 1.2554 vs
#' 1.255). The two series have only ~0.94 sorted correlation, so
#' `ebrecipe::krw_firms` is a *different* (related but not identical) dataset, not
#' a rescaling. Until that is reconciled, the beta-GMM uses KRW's real input as
#' the parity ground truth (Decision B: the Matlab is the authority).
#'
#' @param characteristic `"race"` or `"gender"`.
#' @return A list with `theta_hat`, `s`, unique row-level `unit_id`, `industry`
#'   (the industry-group code from column 1), and `metadata` documenting that the
#'   IDs are row IDs because the headerless Matlab GMM matrix carries industry but
#'   no firm identifier. Shaped so it can be passed directly to
#'   `gp_estimation_core()` (which reads `$theta_hat` / `$s`).
#' @keywords internal
#' @noRd
gp_krw_gmm_input <- function(characteristic = c("race", "gender")) {
  characteristic <- match.arg(characteristic)
  f <- sprintf("extdata/krw-gmm-input/theta_estimates_matlab_%s.csv", characteristic)
  path <- system.file(f, package = "gradepath")
  if (path == "" || !file.exists(path)) {
    # dev fallback: read straight from inst/ when not yet installed (load_all).
    path <- file.path("inst", f)
  }
  if (!file.exists(path)) {
    .gradepath_abort(
      sprintf("KRW GMM input for '%s' not found (expected %s).", characteristic, f),
      class = "gradepath_validation_error"
    )
  }
  # Headerless [industry_d, theta_hat, s]; header = FALSE so row 1 is data.
  m <- utils::read.csv(path, header = FALSE)
  n <- nrow(m)
  list(
    theta_hat = as.numeric(m[[2]]),
    s         = as.numeric(m[[3]]),
    unit_id   = sprintf("krw_%s_%03d", characteristic, seq_len(n)),
    industry  = as.integer(m[[1]]),
    metadata  = list(
      source_file = f,
      unit_id_source = "row_id",
      industry_source = "matlab_column_1"
    )
  )
}
