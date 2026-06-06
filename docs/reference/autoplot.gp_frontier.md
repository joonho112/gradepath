# Autoplot methods for gradepath frontier objects

ggplot2
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods so `autoplot(fr)` (a `gp_frontier`), `autoplot(fit)` (a
`gp_fit`), or `autoplot(path)` (a `gp_grade_path`) draw the
information-reliability frontier via
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md)
without attaching ggplot2.

## Usage

``` r
# S3 method for class 'gp_frontier'
autoplot(object, ...)

# S3 method for class 'gp_fit'
autoplot(object, ...)

# S3 method for class 'gp_grade_path'
autoplot(object, ...)
```

## Arguments

- object:

  A `gp_frontier`, `gp_fit`, or `gp_grade_path`.

- ...:

  Further arguments forwarded to
  [`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md).

## Value

A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
printed.

## See also

[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md)

Other gradepath-plots:
[`autoplot.gp_report_card()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_report_card.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)

## Examples

``` r
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
autoplot(fit)   # the information-reliability frontier for a gp_fit

```
