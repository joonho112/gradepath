# gradepath ggplot2 theme (KRW house style)

A clean ggplot2 theme matching the Kline-Rose-Walters figure house
style: built on
[`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html),
with the top and right spines removed, a light horizontal grid, no
vertical grid, and the legend at the bottom. It is the default theme
inside every gradepath plot verb, and is exported so you can re-apply or
extend it on any ggplot with `p + theme_gradepath()`.

## Usage

``` r
theme_gradepath(base_size = 11, base_family = "", ...)
```

## Arguments

- base_size:

  Base font size in points. Default `11`.

- base_family:

  Base font family. Default `""` (the device default).

- ...:

  Further arguments forwarded to
  [`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

## Value

A ggplot2 theme object (class `c("theme", "gg")`) to add to a plot with
`+`.

## See also

[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md)

Other gradepath-plots:
[`autoplot.gp_frontier()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md),
[`autoplot.gp_report_card()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_report_card.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`scale_color_gradepath()`](https://joonho112.github.io/gradepath/reference/scale_color_gradepath.md)

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  theme_gradepath()

```
