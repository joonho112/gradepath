#' gradepath: Discrimination Report Cards via a Native Kline-Rose-Walters Core
#'
#' @description
#' `gradepath` turns firm-level (or name-level) empirical-Bayes estimates into a
#' discrimination *report card*: it ranks units by a posterior outranking
#' probability, solves a grade integer-program along an information-reliability
#' frontier, and assembles a grade-plus-interval table you can print, plot, and
#' export. It is an R replication companion for the Kline-Rose-Walters (2024) "A
#' Discrimination Report Card" pipeline, built on a native one-level KRW core
#' (beta-GMM precision standardization -> deconvolution -> posterior weights ->
#' pairwise outranking -> grading -> frontier -> report card). The one-call entry
#' point is [krw_report_card()] (alias [gradepath()]); start with its bundled
#' tiny-input example.
#'
#' @details
#' `gradepath` is a development build for replicating the Kline-Rose-Walters
#' (2024) "A Discrimination Report Card" pipeline in R. The package implements
#' the scoped M1 one-level race/gender core: typed objects, KRW beta-GMM
#' precision standardization, one-level native deconvolution, posterior weight
#' recomputation, pairwise outranking, grade-IP solving, frontier metrics,
#' report-card assembly, public one-level APIs, and the M1 acceptance harness.
#' The bundled M1 gate records `ACCEPTED` under the predeclared scoped policy;
#' names and the two-level within-industry share remain outside that M1 gate.
#'
#' The package also includes the M2 two-level / industry surface: native
#' two-level GMM and deconvolution, theta pushforward, posterior summaries,
#' same-industry override and five Pi matrices, deterministic quadrature,
#' seeded simulation fallback, industry grading, N10 support guards, and an
#' explicit M2 acceptance ledger. M2 is `PARTIAL_ACCEPTED`: exact L01 industry
#' grade counts pass; gender continuous L02 rows are promoted to `banded` by
#' intermediate fixture parity; and race continuous L02 rows remain
#' `approximate` because the `Pi_theta` fixture gap is `0.0121756787 > 0.01`.
#' N10 support rows are synthetic evidence guards, not registered paper-value
#' passes. M2 `PROMOTED` / `banded` rows are fixture-parity evidence, not direct
#' reproduction of the paper's industry DR, tau, or R2 values.
#'
#' `gradepath` imports the `ebrecipe` package as a stage-1 input container, a
#' controlled source of low-level primitives behind the seam in
#' `R/seam-ebrecipe.R`, and a one-level cross-check oracle.
#'
#' @section Public surface:
#' The reference index (see the pkgdown site, or `help(package = "gradepath")`)
#' is grouped by `@family`. The headline entry points are:
#' \itemize{
#'   \item Pipeline -- [krw_report_card()] / [gradepath()] and [gp_control()].
#'   \item Grade engine -- [gp_grade()], [gp_grade_path()], [gp_select_grade()],
#'     [gp_preview()].
#'   \item Social choice -- [gp_pairwise()]; frontier metrics [gp_frontier()],
#'     [gp_r2()], [gp_krw_r2()]; and [gp_report_card()].
#'   \item Accessors -- [get_grades()], [get_report_card()], [get_pairwise()],
#'     [get_posterior()], [get_prior()], [get_control()], plus [coef()] and
#'     [as.data.frame()] methods for a `gp_fit`.
#'   \item Plots -- [gp_plot_frontier()], [gp_plot_posterior_contrast()],
#'     [gp_plot_report_card()], [gp_plot_discordance()], the [theme_gradepath()]
#'     and `scale_color_gradepath()` / `scale_fill_gradepath()` family, and the
#'     `autoplot()` methods for the gradepath objects.
#'   \item Calibration -- [gp_calibrate()], a seeded Monte-Carlo check of the
#'     grading method's frequentist guarantees.
#'   \item Console -- the `format_gp_*_cli()` family (`format_gp_fit_cli()`,
#'     `format_gp_report_card_cli()`, `format_gp_frontier_cli()`,
#'     `format_gp_grade_path_cli()`, `format_gp_grade_fit_cli()`,
#'     `format_gp_pairwise_cli()`), the result-first ASCII renderers behind
#'     `print()`.
#'   \item Harness -- [gp_check()], [gp_run_all()], [gp_validate_targets()].
#'   \item Two-level / industry (M2) -- [gp_twolevel_pairwise()],
#'     [gp_twolevel_grade()], [gp_twolevel_report_card()], and the M2 ledger
#'     [gp_m2_acceptance()], [gp_m2_promoted_registry()], [gp_m2_status()].
#' }
#'
#' @section Solver prerequisite:
#' Solver-backed grading uses Gurobi by default. Open backends are explicit:
#' HiGHS uses the `highs` package directly, while GLPK and SYMPHONY route
#' through ROI plugins. The package installs, loads, and most structural tests
#' run without a Gurobi license.
#'
#' @section Status:
#' This is a pre-release development build. The bundled M1 gate records
#' `ACCEPTED` for the scoped one-level race/gender surface. The bundled M2 gate
#' records `PARTIAL_ACCEPTED`: L01 industry grade counts pass exactly, gender L02
#' continuous rows are promoted to `banded` by fixture parity, race L02
#' continuous rows remain `approximate` because the `Pi_theta` fixture gap is
#' above the 0.01 promotion boundary, and N10 support rows are synthetic guard
#' evidence rather than registered paper-value passes. The M2 simulation
#' fallback is retained as approximate evidence; deterministic quadrature is the
#' primary M2 path.
#'
#' @section Getting started:
#' Read `vignette("a1-getting-started", package = "gradepath")` for a one-call
#' walkthrough, then reach for [krw_report_card()]. The applied vignettes
#' (`a1`--`a5`) cover the workflow, reading and exporting results, the figure
#' cookbook, and the solvers; the methodological vignettes (`m1`--`m5`) carry the
#' estimand-by-estimand derivations and the calibration scope.
#'
#' @references
#' Kline, P., Rose, E., & Walters, C. (2024). A Discrimination Report Card.
#' *American Economic Review*, 114(8), 2472-2525. \doi{10.1257/aer.20230700}
#'
#' @seealso [krw_report_card()], [gp_report_card()], [gp_control()]
#'
#' @keywords internal
"_PACKAGE"
