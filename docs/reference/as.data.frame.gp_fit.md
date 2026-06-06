# Report-card table from a gradepath fit as a data frame

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
for `gp_fit`: returns the report-card rows – the underlying data frame
of the `gp_report_card` slot (`get_report_card(fit)$table`) – with row
names reset to the default sequence. Provided because
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) is the
base generic an R user reaches for first (`GP-DEC-13-A`, Chapter 13
patterns; Chapter 14 ux). Read-only; recomputes nothing. Grade labels
are integers in `{1, ..., n}` and carry no ranking-superiority statement
of any kind.

## Usage

``` r
# S3 method for class 'gp_fit'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `gp_fit` object.

- ...:

  Unused; accepted for S3 method extensibility. No arguments are
  forwarded; supplying any has no effect.

## Value

A `data.frame`: the report-card table (one row per unit), row names
reset to the default integer sequence.

## See also

[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md),
[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md)

Other gradepath-accessors:
[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md),
[`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md),
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
[`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md),
[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)

## Examples

``` r
# Instant: the report-card rows of the bundled tiny fit as a data frame.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
head(as.data.frame(fit))   # the report-card rows
#>       id  label grade sort_rank selected_lambda posterior_mean       lower
#> 1 firm03 firm03     1         1            0.25      0.1330361 0.004466502
#> 2 firm05 firm05     1         2            0.25      0.1547879 0.094098498
#> 3 firm10 firm10     1         3            0.25      0.1372338 0.024811477
#> 4 firm18 firm18     1         4            0.25      0.1166670 0.003301362
#> 5 firm24 firm24     1         5            0.25      0.1206680 0.004841016
#> 6 firm01 firm01     1         6            0.25      0.1323136 0.080808125
#>       upper estimate     se
#> 1 0.2014888   0.3431 0.2386
#> 2 0.1890643   0.4332 0.1309
#> 3 0.1749419   0.2895 0.1142
#> 4 0.1895925   0.1811 0.1902
#> 5 0.1782374   0.2017 0.1398
#> 6 0.1631853   0.2355 0.0838
```
