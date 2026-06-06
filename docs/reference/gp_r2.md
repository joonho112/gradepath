# Grade-level share of bias-corrected estimate variance

`gp_r2()` reports the percent (or proportion) of the bias-corrected
estimate variance that is explained by the selected grade means – a
quick between-grade variance-share diagnostic for any grade fit. It is
the lightweight generic helper, not the paper-faithful KRW recipe; for
the Kline-Rose-Walters race/gender R-squared target definition use
[`gp_krw_r2()`](https://joonho112.github.io/gradepath/reference/gp_krw_r2.md).

## Usage

``` r
gp_r2(
  x = NULL,
  estimate = NULL,
  se = NULL,
  grade = NULL,
  scale = c("percent", "proportion")
)
```

## Arguments

- x:

  Optional source of the three vectors below. A `gp_fit` (its per-unit
  `estimate`, `se`, and `grade` are read), or a data frame carrying
  `estimate`, `se`, and `grade` columns. `NULL` (default) means supply
  the vectors explicitly via `estimate`/`se`/`grade`.

- estimate, se, grade:

  Numeric vectors used when `x` is `NULL`: the bias-corrected point
  estimates, their strictly positive standard errors, and the integer
  grade labels (relabeled to be contiguous internally). All three must
  be the same length. Ignored when `x` is supplied.

- scale:

  Character scalar; `"percent"` (default) returns the share on a 0-100
  scale, `"proportion"` returns it on a 0-1 scale.

## Value

Numeric scalar: the between-grade variance share on the requested scale,
or `NA_real_` when it is undefined (fewer than two units, or a
non-positive bias-corrected total variance).

## See also

[`gp_krw_r2()`](https://joonho112.github.io/gradepath/reference/gp_krw_r2.md),
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)

Other gradepath-frontier:
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md),
[`gp_krw_r2()`](https://joonho112.github.io/gradepath/reference/gp_krw_r2.md)

## Examples

``` r
# Variance share explained by the selected grades of the bundled fit.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
gp_r2(fit)                       # percent scale
#> [1] 11.72085
gp_r2(fit, scale = "proportion") # 0-1 scale
#> [1] 0.1172085
```
