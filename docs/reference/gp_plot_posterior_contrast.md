# Plot the grade-sorted posterior-contrast heatmap

`gp_plot_posterior_contrast()` draws the Kline-Rose-Walters grade-sorted
posterior-contrast heatmap (their Figure 5a / 12a): the J x J pairwise
outranking matrix \\\Pi\\ (\\\pi\_{ij}\\, the posterior probability that
unit `i` is more extreme than unit `j`), with rows and columns reordered
by `(selected grade, endpoint grade, id)` so the grade blocks are
contiguous. Black segments separate the grade blocks; when
`threshold = TRUE`, cells whose contrast reaches the different-grade
threshold \\1 / (1 + \lambda)\\ are outlined in red. Pass it a `gp_fit`
(it reads the pairwise and grades from the fit).

## Usage

``` r
gp_plot_posterior_contrast(
  x,
  threshold = TRUE,
  selected_grades = NULL,
  endpoint_grades = NULL,
  selected_lambda = NULL,
  ...
)
```

## Arguments

- x:

  A `gp_fit` (primary): uses `x$pairwise$matrix`, `x$selected_grade`,
  and `x$grade_path`. A `gp_pairwise` is also accepted together with
  explicit `selected_grades` (and optionally `endpoint_grades`) vectors
  plus `selected_lambda`.

- threshold:

  Logical; outline cells reaching the different-grade threshold
  `1 / (1 + selected_lambda)`. Default `TRUE`.

- selected_grades:

  Integer grade vector (one per unit, in `ids` order), required only
  when `x` is a bare `gp_pairwise`.

- endpoint_grades:

  Integer grade vector used as the secondary sort key; defaults to
  `selected_grades` when omitted. Only used with a `gp_pairwise`.

- selected_lambda:

  Numeric selected lambda used for the threshold, required only when `x`
  is a bare `gp_pairwise` (read from the fit otherwise).

- ...:

  Reserved for future arguments; currently unused.

## Value

A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
printed.

## See also

[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md),
[`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)

Other gradepath-plots:
[`autoplot.gp_frontier()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md),
[`autoplot.gp_report_card()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_report_card.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)

## Examples

``` r
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))

gp_plot_posterior_contrast(fit)                    # grade-sorted Pi heatmap

gp_plot_posterior_contrast(fit, threshold = FALSE) # no different-grade outline

```
