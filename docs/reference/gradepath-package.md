# gradepath: Discrimination Report Cards via a Native Kline-Rose-Walters Core

`gradepath` turns firm-level (or name-level) empirical-Bayes estimates
into a discrimination *report card*: it ranks units by a posterior
outranking probability, solves a grade integer-program along an
information-reliability frontier, and assembles a grade-plus-interval
table you can print, plot, and export. It is an R replication companion
for the Kline-Rose-Walters (2024) "A Discrimination Report Card"
pipeline, built on a native one-level KRW core (beta-GMM precision
standardization -\> deconvolution -\> posterior weights -\> pairwise
outranking -\> grading -\> frontier -\> report card). The one-call entry
point is
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
(alias
[`gradepath()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md));
start with its bundled tiny-input example.

## Details

`gradepath` is a development build for replicating the
Kline-Rose-Walters (2024) "A Discrimination Report Card" pipeline in R.
The package implements the scoped M1 one-level race/gender core: typed
objects, KRW beta-GMM precision standardization, one-level native
deconvolution, posterior weight recomputation, pairwise outranking,
grade-IP solving, frontier metrics, report-card assembly, public
one-level APIs, and the M1 acceptance harness. The bundled M1 gate
records `ACCEPTED` under the predeclared scoped policy; names and the
two-level within-industry share remain outside that M1 gate.

The package also includes the M2 two-level / industry surface: native
two-level GMM and deconvolution, theta pushforward, posterior summaries,
same-industry override and five Pi matrices, deterministic quadrature,
seeded simulation fallback, industry grading, N10 support guards, and an
explicit M2 acceptance ledger. M2 is `PARTIAL_ACCEPTED`: exact L01
industry grade counts pass; gender continuous L02 rows are promoted to
`banded` by intermediate fixture parity; and race continuous L02 rows
remain `approximate` because the `Pi_theta` fixture gap is
`0.0121756787 > 0.01`. N10 support rows are synthetic evidence guards,
not registered paper-value passes. M2 `PROMOTED` / `banded` rows are
fixture-parity evidence, not direct reproduction of the paper's industry
DR, tau, or R2 values.

`gradepath` imports the `ebrecipe` package as a stage-1 input container,
a controlled source of low-level primitives behind the seam in
`R/seam-ebrecipe.R`, and a one-level cross-check oracle.

## Public surface

The reference index (see the pkgdown site, or
[`help(package = "gradepath")`](https://rdrr.io/pkg/gradepath/man)) is
grouped by `@family`. The headline entry points are:

- Pipeline –
  [`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
  /
  [`gradepath()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
  and
  [`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md).

- Grade engine –
  [`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md),
  [`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
  [`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md),
  [`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md).

- Social choice –
  [`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md);
  frontier metrics
  [`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md),
  [`gp_r2()`](https://joonho112.github.io/gradepath/reference/gp_r2.md),
  [`gp_krw_r2()`](https://joonho112.github.io/gradepath/reference/gp_krw_r2.md);
  and
  [`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md).

- Accessors –
  [`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
  [`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md),
  [`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
  [`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md),
  [`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
  [`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md),
  plus [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  for a `gp_fit`.

- Plots –
  [`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
  [`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
  [`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
  [`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
  the
  [`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)
  and
  [`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md)
  /
  [`scale_fill_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md)
  family, and the
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods for the gradepath objects.

- Calibration –
  [`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md),
  a seeded Monte-Carlo check of the grading method's frequentist
  guarantees.

- Console – the `format_gp_*_cli()` family
  ([`format_gp_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md),
  [`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md),
  [`format_gp_frontier_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md),
  [`format_gp_grade_path_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md),
  [`format_gp_grade_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md),
  [`format_gp_pairwise_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md)),
  the result-first ASCII renderers behind
  [`print()`](https://rdrr.io/r/base/print.html).

- Harness –
  [`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md),
  [`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md),
  [`gp_validate_targets()`](https://joonho112.github.io/gradepath/reference/gp_validate_targets.md).

- Two-level / industry (M2) –
  [`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md),
  [`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
  [`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md),
  and the M2 ledger
  [`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
  [`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md),
  [`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md).

## Solver prerequisite

Solver-backed grading uses Gurobi by default. Open backends are
explicit: HiGHS uses the `highs` package directly, while GLPK and
SYMPHONY route through ROI plugins. The package installs, loads, and
most structural tests run without a Gurobi license.

## Status

This is a pre-release development build. The bundled M1 gate records
`ACCEPTED` for the scoped one-level race/gender surface. The bundled M2
gate records `PARTIAL_ACCEPTED`: L01 industry grade counts pass exactly,
gender L02 continuous rows are promoted to `banded` by fixture parity,
race L02 continuous rows remain `approximate` because the `Pi_theta`
fixture gap is above the 0.01 promotion boundary, and N10 support rows
are synthetic guard evidence rather than registered paper-value passes.
The M2 simulation fallback is retained as approximate evidence;
deterministic quadrature is the primary M2 path.

## Getting started

Read
[`vignette("a1-getting-started", package = "gradepath")`](https://joonho112.github.io/gradepath/articles/a1-getting-started.md)
for a one-call walkthrough, then reach for
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md).
The applied vignettes (`a1`–`a5`) cover the workflow, reading and
exporting results, the figure cookbook, and the solvers; the
methodological vignettes (`m1`–`m5`) carry the estimand-by-estimand
derivations and the calibration scope.

## References

Kline, P., Rose, E., & Walters, C. (2024). A Discrimination Report Card.
*American Economic Review*, 114(8), 2472-2525.
[doi:10.1257/aer.20230700](https://doi.org/10.1257/aer.20230700)

## See also

[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md),
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]

Authors:

- JoonHo Lee <jlee296@ua.edu>
  ([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]
