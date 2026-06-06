# Autoplot method for a gradepath report card

A thin
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
method so `autoplot(card)` (a `gp_report_card`) draws the report-card
caterpillar via
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md)
without attaching ggplot2. (There is no discordance autoplot method:
`autoplot.gp_frontier` already maps to the frontier figure.)

## Usage

``` r
# S3 method for class 'gp_report_card'
autoplot(object, ...)
```

## Arguments

- object:

  A `gp_report_card`.

- ...:

  Further arguments forwarded to
  [`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md).

## Value

A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
printed.

## See also

[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)

Other gradepath-plots:
[`autoplot.gp_frontier()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)

## Examples

``` r
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
autoplot(get_report_card(fit))   # the report-card caterpillar for a gp_fit

```
