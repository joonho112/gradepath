# Plot the KRW conditional-discordance heatmap (lower triangle)

`gp_plot_discordance()` draws the Kline-Rose-Walters conditional
discordance matrix (their Figure 10a / 17a) as a lower-triangular
heatmap. The cell for a grade pair is the conditional discordance rate
\\DR\_{g,g'}\\ between the more-extreme grade `g` (columns, grade 1 on
the left) and the less-extreme grade `g'` (rows, grade 2 at the top):
the posterior share of disagreeing rankings among cross-grade pairs.
Higher DR is darker; with `annotate = TRUE` the value is written in each
cell with auto-contrasting text. Pass it a
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md)
(or, for convenience, the `gp_fit` it was built from, or a bare DR
matrix).

## Usage

``` r
gp_plot_discordance(x, annotate = TRUE, ...)
```

## Arguments

- x:

  A `gp_frontier` (reads `x$dr_matrix`), a `gp_fit` (its frontier is
  rebuilt with
  [`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md)
  from `x$grade_path` and `x$pairwise`), or a numeric conditional-DR
  matrix (square; lower triangle populated, upper triangle and diagonal
  `NA`).

- annotate:

  Logical; write each cell's DR value (`"%.3f"`) with an
  auto-contrasting colour. Default `TRUE`.

- ...:

  Reserved for future arguments; currently unused.

## Value

A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
printed.

## See also

[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md),
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md)

Other gradepath-plots:
[`autoplot.gp_frontier()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md),
[`autoplot.gp_report_card()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_report_card.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)

## Examples

``` r
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))

# Build a frontier from the fit's grade path and pairwise.
fr <- gp_frontier(fit$grade_path, pairwise = get_pairwise(fit))

gp_plot_discordance(fr)                   # annotated lower-triangular heatmap

gp_plot_discordance(fr, annotate = FALSE) # heatmap without cell values

gp_plot_discordance(fit)                  # convenience: a gp_fit

gp_plot_discordance(fr$dr_matrix)         # a bare dr_matrix

```
