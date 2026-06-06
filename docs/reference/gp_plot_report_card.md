# Plot the KRW report card (horizontal point-and-interval)

`gp_plot_report_card()` draws the Kline-Rose-Walters unit-level report
card (their Figures 7, 8, 14, 15) as a horizontal point-and-interval
("caterpillar") plot: each unit's posterior mean (point) with its stored
credible interval (`lower`, `upper`), one unit per row, ordered so the
most-extreme unit is on top, and coloured by the unit's selected grade
using the shared ordinal palette
([`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md)).
The interval already on the card is displayed as-is and never
recomputed. The figure is neutral: it carries no welfare-verdict text,
because grade 1's welfare meaning flips between firms and names. Pass it
a
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)
(or the `gp_fit` that carries one).

## Usage

``` r
gp_plot_report_card(
  x,
  order = "rank",
  ci_level = NULL,
  max_rows = NULL,
  zero_line = "auto",
  ...
)
```

## Arguments

- x:

  A `gp_report_card` (primary). A `gp_fit` is also accepted, in which
  case its `$report_card` is used.

- order:

  Row order, most-extreme unit on top. `"rank"` (default) uses the
  endpoint/Condorcet `sort_rank` (the canonical endpoint-sorted figure);
  `"mean"` orders by `posterior_mean`; `"grade"` groups by grade then
  `sort_rank`; `"label"` is alphabetical.

- ci_level:

  Optional number in `(0, 1)` used ONLY to label the caption, e.g. `0.9`
  -\> "90% credible interval". `NULL` (default) reads the level from
  `x$control$interval_level` (default `0.90`); the interval itself is
  never recomputed.

- max_rows:

  Optional positive integer keeping only the top `max_rows` units by the
  chosen order, with a "Showing N of M" subtitle (never silent). `NULL`
  (default) shows all units.

- zero_line:

  Dashed reference line at `x = 0`: `"auto"` (default) draws it only
  when the intervals straddle zero (`min(lower) < 0 < max(upper)`);
  `"on"` always; `"off"` never.

- ...:

  Reserved for future arguments; currently unused.

## Value

A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
printed.

## See also

[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md),
[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)

Other gradepath-plots:
[`autoplot.gp_frontier()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md),
[`autoplot.gp_report_card()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_report_card.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)

## Examples

``` r
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
card <- get_report_card(fit)

gp_plot_report_card(card)                  # endpoint-sorted, coloured by grade

gp_plot_report_card(card, order = "grade") # grouped by grade

gp_plot_report_card(card, max_rows = 10)   # top 10 only

gp_plot_report_card(fit)                   # convenience: a gp_fit

```
