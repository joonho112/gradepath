#' KRW (2024) firm-level callback-discrimination estimates
#'
#' Per-firm empirical-Bayes input estimates from Kline, Rose and Walters
#' (2024), "A Discrimination Report Card": the race (White minus Black) and
#' gender (Male minus Female) callback-rate contact gaps for 97 large U.S.
#' employers, with their standard errors, plus the firm's industry grouping.
#' This is the worked example the `gradepath` pipeline runs on; the race
#' contrast is the M1 headline (97 firms, 3 grades 2/81/14 at lambda = 0.25).
#'
#' @format A data frame with 97 rows and 7 columns:
#' \describe{
#'   \item{firm_id}{integer. Firm identifier; non-contiguous (matches the KRW
#'     replication-archive firm numbering, which has gaps).}
#'   \item{theta_hat_race}{numeric. Estimated White minus Black callback-rate
#'     gap (the race discrimination contact gap).}
#'   \item{se_race}{numeric. Standard error of `theta_hat_race`.}
#'   \item{theta_hat_gender}{numeric. Estimated Male minus Female callback-rate
#'     gap (the gender contact gap).}
#'   \item{se_gender}{numeric. Standard error of `theta_hat_gender`.}
#'   \item{industry}{character. SIC-grouped industry code (used by the
#'     two-level / industry workflow).}
#'   \item{label}{character. Human-readable firm label.}
#' }
#'
#' @details
#' The estimates are taken verbatim from `ebrecipe::krw_firms`; the `industry`
#' and `label` columns are joined from the KRW replication metadata on
#' `firm_id`. A `sample_stats` attribute records the sampling-frame counts
#' (108 firms / 83,643 observations before filtering; 97 firms / 78,910 after).
#'
#' @section Provenance and parity (read before replicating KRW):
#' WARNING -- `krw_firms` is a PUBLIC EXAMPLE dataset, NOT the KRW parity input.
#' It does not reproduce the published KRW (2024) results and must not be used
#' for replication. The bundled series is on a different numeric scale from the
#' Matlab GMM input KRW actually consumed: running the beta-GMM on `krw_firms`
#' yields a spurious large-beta optimum (race beta = 2.1 and gender beta = 3.0,
#' versus the published 0.51 and 1.26), and the gender path errors out in the
#' deconvolution with `DECONV_BOUNDARY_ERROR`. It therefore does NOT reproduce
#' KRW Table 3.
#'
#' The parity-faithful input is KRW's real Matlab GMM series, shipped at
#' `inst/extdata/krw-gmm-input/` and read directly with, for example,
#' `read.csv(system.file("extdata/krw-gmm-input/theta_estimates_matlab_race.csv",
#' package = "gradepath"), header = FALSE)` (column 2 = `theta_hat`, column 3 =
#' `s`). Use that series (not `krw_firms`) whenever you
#' need to reproduce KRW. `krw_firms` is retained only as a small, public,
#' end-to-end worked example of the pipeline's mechanics.
#'
#' @examples
#' data(krw_firms)
#' str(krw_firms)
#'
#' @source Kline, P., Rose, E. K., and Walters, C. R. (2024). A Discrimination
#'   Report Card. *American Economic Review*, 114(8), 2472--2525.
#'   \doi{10.1257/aer.20230700}. Estimates via `ebrecipe::krw_firms`; industry
#'   grouping from the replication archive.
#' @keywords datasets
"krw_firms"

#' gradepath replication registry (curated KRW paper-value targets)
#'
#' A curated subset of the companion `paper_values.csv` (the published
#' Kline-Rose-Walters 2024 ground-truth registry), restricted to the
#' quantities `gradepath` owns and verifies. Each row carries the published
#' value, its unit, an absolute per-id tolerance, a tolerance class, and the
#' milestone at which it is asserted. It is read by the replication harness
#' (`gp_check()` / `gp_run_all()`); the harness embeds no tolerances of its
#' own. Pairwise outranking probabilities (K05--K07) deliberately have **no**
#' registry row -- they are verified by a golden-master on the posterior
#' weight matrix `W` instead.
#'
#' @format A data frame with 69 rows and 7 columns:
#' \describe{
#'   \item{id}{character. Stable target key: the `paper_values.csv` key
#'     verbatim where one exists; otherwise a synthetic key (the single such
#'     row is `n9_strict4_ngrades`, the N9 strict-order grade-label guard).}
#'   \item{paper_value}{character. Published value (stored character,
#'     CSV-faithful; coerce with `as.numeric`).}
#'   \item{unit}{character. One of `count`, `percent`, `proportion`,
#'     `correlation`, `sd`, `other`.}
#'   \item{tolerance}{character. Absolute per-id PASS tolerance (coerce to
#'     numeric); `0` for exact counts.}
#'   \item{class}{character. `exact`, `banded`, or `approximate` (tolerance
#'     provenance, per the Chapter 17 ledger).}
#'   \item{milestone}{character. `M1` (current one-level firm gate),
#'     `M1_NAMES` (deferred names pipeline), or `M2` (two-level / industry).}
#'   \item{quantity}{character. Human-readable description of the target.}
#' }
#'
#' @details
#' Values are sourced verbatim from the companion `paper_values.csv` (the
#' single source of truth); the `class`, `tolerance`, and `milestone` columns
#' are attached from the package's tolerance ledger. The
#' NPMLE name atoms use the verbatim CSV keys `names_npmle_mass1/2/3`. The
#' registry is intentionally an M1-and-M2 core subset; names rows are retained
#' under `M1_NAMES`, and Table F5/F6 industry cells are assigned to `M2` until
#' the two-level / industry report-card pipeline is live.
#'
#' @examples
#' data(gp_registry)
#' str(gp_registry)
#' # Count the rows per milestone gate.
#' table(gp_registry$milestone)
#'
#' @source Curated from the KRW (2024) companion `paper_values.csv`; tolerance
#'   classes from the package's tolerance ledger. Names rows are
#'   retained in the registry for paper-value traceability but assigned to the
#'   deferred `M1_NAMES` milestone until the names-specific NPMLE/ranking
#'   pipeline is bundled and tested.
#' @keywords datasets
"gp_registry"
