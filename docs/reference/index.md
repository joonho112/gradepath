# Package index

## Complete pipeline

One-call entry points that run the full Kline-Rose-Walters recipe, plus
the run-settings constructor.
([`gradepath()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
is an alias of
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md).)

- [`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
  [`gradepath()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
  : Run the one-level KRW report-card pipeline in a single call
- [`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)
  : Create a gradepath run-control object

## The grade engine

Solve the grade integer program, trace the grade path, select a grade,
and preview a run.

- [`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md)
  : Solve a single-penalty grade fit
- [`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
  : Solve the grade path over a lambda grid
- [`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)
  : Select a stored grade fit from a solved grade path
- [`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md)
  : Preview a gradepath run without solving

## Pairwise outranking

Posterior pairwise outranking probabilities between units.

- [`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)
  : Pairwise outranking probability matrix

## Frontier & reliability

The information-reliability frontier and the R-squared reliability
measures.

- [`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md)
  : Summarize the frontier of a solved grade path
- [`gp_r2()`](https://joonho112.github.io/gradepath/reference/gp_r2.md)
  : Grade-level share of bias-corrected estimate variance
- [`gp_krw_r2()`](https://joonho112.github.io/gradepath/reference/gp_krw_r2.md)
  : KRW grade R-squared from the rsquared.do recipe

## Report cards

Assemble the discrimination report card from a fit.

- [`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)
  : Assemble a unit-level report card

## Plotting

The four signature ggplot2 figures, the theme, the grade palette, and
autoplot methods.

- [`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md)
  : Plot the KRW information-reliability frontier
- [`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md)
  : Plot the grade-sorted posterior-contrast heatmap
- [`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md)
  : Plot the KRW report card (horizontal point-and-interval)
- [`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md)
  : Plot the KRW conditional-discordance heatmap (lower triangle)
- [`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)
  : gradepath ggplot2 theme (KRW house style)
- [`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md)
  [`scale_colour_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md)
  [`scale_fill_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md)
  : gradepath ordinal grade colour and fill scales
- [`autoplot(`*`<gp_frontier>`*`)`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md)
  [`autoplot(`*`<gp_fit>`*`)`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md)
  [`autoplot(`*`<gp_grade_path>`*`)`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md)
  : Autoplot methods for gradepath frontier objects
- [`autoplot(`*`<gp_report_card>`*`)`](https://joonho112.github.io/gradepath/reference/autoplot.gp_report_card.md)
  : Autoplot method for a gradepath report card

## Calibration

Seeded Monte Carlo harness for calibration and regret diagnostics.

- [`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md)
  : Seeded Monte-Carlo calibration of the one-level KRW grading method

## Two-level / industry (M2)

The two-level industry layer – pairwise, grading, acceptance, and the
promotion registry.

- [`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md)
  : Wrap two-level Pi matrices as solver-ready pairwise objects
- [`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md)
  : Grade two-level industry matrices
- [`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)
  : Return a two-level industry report card
- [`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md)
  : Build the M2 acceptance scorecard
- [`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md)
  : Effective M2 registry after fixture-gate promotion
- [`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md)
  : Summarize the M2 acceptance status

## Accessors

Typed accessors that pull components out of a fitted object.

- [`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md)
  : Extract per-unit grades from a gradepath fit
- [`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)
  : Extract the report card from a gradepath fit
- [`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md)
  : Extract the pairwise outranking object from a gradepath fit
- [`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md)
  : Extract the native posterior from a gradepath fit
- [`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md)
  : Extract the native prior from a gradepath fit
- [`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md)
  : Extract the run-control object from a gradepath fit

## Replication harness

Run the bundled replication, check it, and validate published targets.

- [`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md)
  : Compare a replicated value against the gradepath registry
- [`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md)
  : Run a registry-driven gradepath scorecard
- [`gp_validate_targets()`](https://joonho112.github.io/gradepath/reference/gp_validate_targets.md)
  : Validate replicated target values against the gradepath registry

## Console formatters (CLI)

Opt-in cli-decorated formatters that render rich console output.

- [`format_gp_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md)
  [`summary(`*`<gp_fit>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md)
  [`print(`*`<gp_fit>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md)
  :

  Format a `gp_fit` for the console

- [`format_gp_frontier_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md)
  [`summary(`*`<gp_frontier>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md)
  [`print(`*`<gp_frontier>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md)
  :

  Format a `gp_frontier` for the console

- [`format_gp_grade_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md)
  [`summary(`*`<gp_grade_fit>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md)
  [`print(`*`<gp_grade_fit>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md)
  :

  Format a `gp_grade_fit` for the console

- [`format_gp_grade_path_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md)
  [`summary(`*`<gp_grade_path>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md)
  [`print(`*`<gp_grade_path>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md)
  :

  Format a `gp_grade_path` for the console

- [`format_gp_pairwise_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md)
  [`summary(`*`<gp_pairwise>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md)
  [`print(`*`<gp_pairwise>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md)
  :

  Format a `gp_pairwise` for the console

- [`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)
  [`summary(`*`<gp_report_card>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)
  [`print(`*`<gp_report_card>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)
  [`as.data.frame(`*`<gp_report_card>`*`)`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)
  :

  Format a `gp_report_card` for the console

## S3 methods

Print, summary, coef, and data-frame coercion for gradepath classes.

- [`print(`*`<gp_calibration>`*`)`](https://joonho112.github.io/gradepath/reference/print.gp_calibration.md)
  :

  Print a `gp_calibration` – a compact, result-first one-liner

- [`print(`*`<gp_summary>`*`)`](https://joonho112.github.io/gradepath/reference/print.gp_summary.md)
  :

  Print a `gp_summary`

- [`summary(`*`<gp_calibration>`*`)`](https://joonho112.github.io/gradepath/reference/summary.gp_calibration.md)
  :

  Summarize a `gp_calibration` (typed one-row data frame)

- [`coef(`*`<gp_fit>`*`)`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md)
  : Grade vector from a gradepath fit (base-generic synonym for
  get_grades)

- [`as.data.frame(`*`<gp_fit>`*`)`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md)
  : Report-card table from a gradepath fit as a data frame

## Datasets

Bundled datasets used throughout the replication.

- [`krw_firms`](https://joonho112.github.io/gradepath/reference/krw_firms.md)
  : KRW (2024) firm-level callback-discrimination estimates
- [`gp_registry`](https://joonho112.github.io/gradepath/reference/gp_registry.md)
  : gradepath replication registry (curated KRW paper-value targets)
