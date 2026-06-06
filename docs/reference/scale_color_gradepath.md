# gradepath ordinal grade colour and fill scales

Discrete ggplot2 scales that map grade levels to the shared gradepath
ordinal palette, so your own plots and the gradepath report-card figures
use identical grade colours. The palette is a perceptually uniform,
colour-blind-safe viridis ramp; grade 1 (the most-extreme theta block)
is the dark, high-contrast end. `scale_colour_gradepath()` is a spelling
alias of `scale_color_gradepath()`, and `scale_fill_gradepath()` is the
fill-aesthetic counterpart. Add one to a ggplot whose `colour`/`fill`
aesthetic is a grade factor.

## Usage

``` r
scale_color_gradepath(n_grades = NULL, ...)

scale_colour_gradepath(n_grades = NULL, ...)

scale_fill_gradepath(n_grades = NULL, ...)
```

## Arguments

- n_grades:

  Optional positive integer fixing the palette length. When `NULL` (the
  default) the palette adapts to the number of grade levels in the data.

- ...:

  Further arguments forwarded to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html)
  (for example `name`, `labels`, `guide`).

## Value

A ggplot2 discrete scale (a `Scale`/`ggproto` object) to add to a plot
with `+`.

## See also

[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md)

Other gradepath-plots:
[`autoplot.gp_frontier()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_frontier.md),
[`autoplot.gp_report_card()`](https://joonho112.github.io/gradepath/reference/autoplot.gp_report_card.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_posterior_contrast()`](https://joonho112.github.io/gradepath/reference/gp_plot_posterior_contrast.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`theme_gradepath()`](https://joonho112.github.io/gradepath/reference/theme_gradepath.md)

## Examples

``` r
library(ggplot2)
df <- data.frame(x = 1:3, y = 1:3, grade = factor(1:3))
ggplot(df, aes(x, y, colour = grade)) +
  geom_point(size = 4) +
  scale_color_gradepath() +
  theme_gradepath()


# The fill counterpart, with a fixed palette length:
ggplot(df, aes(x, y, fill = grade)) +
  geom_tile() +
  scale_fill_gradepath(n_grades = 3)

```
