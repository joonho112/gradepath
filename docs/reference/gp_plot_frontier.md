# Plot the KRW information-reliability frontier

`gp_plot_frontier()` draws the Kline-Rose-Walters
information-reliability frontier (their Figure 6b / 13b): the
optimal-grade path traced in reliability `1 - DR` (x-axis) against the
posterior expected Kendall agreement \\\bar{\tau}\\ (y-axis). Each
frontier-table row is one point and the points are connected as the
Pareto path. The selected-lambda slice is marked with a highlighted
point, and the five naive benchmark rules are overlaid as labelled
shapes when the frontier carries them. Feed it a
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md)
(or, for convenience, the `gp_fit` or `gp_grade_path` it was built
from).

## Usage

``` r
gp_plot_frontier(
  x,
  highlight = "selected",
  benchmarks = TRUE,
  annotate = TRUE,
  pairwise = NULL,
  ...
)
```

## Arguments

- x:

  A `gp_frontier` (primary). For convenience a `gp_fit` or
  `gp_grade_path` is also accepted and converted with
  [`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md);
  benchmarks are shown only if already present (they are never
  force-derived from a bare fit or path).

- highlight:

  Which slice to emphasise: `"selected"` (default) marks the selected
  lambda; `"all"` marks every path point; a numeric value marks the
  nearest path lambda.

- benchmarks:

  Logical; overlay the naive benchmark rules as labelled shapes when
  they are present on the frontier. Default `TRUE`.

- annotate:

  Logical; report the highlighted `(1 - DR, tau_bar)` point in the
  subtitle. Default `TRUE`.

- pairwise:

  Optional `gp_pairwise`, used only when `x` is a `gp_grade_path` (or to
  override a fit's stored pairwise) so the frontier can be rebuilt.

- ...:

  Reserved for future arguments; currently unused.

## Value

A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
printed and stays fully overridable (add `+ theme_*()` / `+ scale_*`).

## See also

[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md),
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md)

Other gradepath-plots:
[`autoplot.gp_frontier()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md),
[`autoplot.gp_report_card()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_report_card.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)

## Examples

``` r
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))

# Build a frontier with the three shrinkage benchmark scores on the fit's ids.
post <- get_posterior(fit)
ids <- get_pairwise(fit)$ids
fr <- gp_frontier(
  fit$grade_path,
  pairwise = get_pairwise(fit),
  benchmark_scores = list(
    raw_estimate     = setNames(post$estimate, ids),
    posterior_mean   = setNames(post$posterior_mean, ids),
    linear_shrinkage = setNames(0.5 * post$estimate, ids)
  )
)

gp_plot_frontier(fr)                    # selected slice + benchmarks

gp_plot_frontier(fr, highlight = "all") # mark every path point

gp_plot_frontier(fit)                   # convenience: a gp_fit

```
