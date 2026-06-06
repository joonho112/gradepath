# Grade vector from a gradepath fit (base-generic synonym for get_grades)

[`coef()`](https://rdrr.io/r/stats/coef.html) method for `gp_fit`: the
"fitted summary" of a gradepath fit is its per-unit grading, so
`coef(fit)` returns the same named integer grade vector as
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md).
Provided because [`coef()`](https://rdrr.io/r/stats/coef.html) is the
base generic an R user reaches for first (`GP-DEC-13-A`, Chapter 13
patterns; Chapter 14 ux). Grade labels are integers in `{1, ..., n}` and
carry no ranking-superiority statement of any kind.

## Usage

``` r
# S3 method for class 'gp_fit'
coef(object, ...)
```

## Arguments

- object:

  A `gp_fit` object.

- ...:

  Unused; accepted for S3 method extensibility. No arguments are
  forwarded; supplying any has no effect.

## Value

A named `integer` vector of grades, identical to `get_grades(object)`.

## See also

[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md)

Other gradepath-accessors:
[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md),
[`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md),
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
[`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md),
[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)

## Examples

``` r
# Instant: coef() is the grade vector of the bundled tiny fit.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
head(coef(fit))        # == head(get_grades(fit))
#> firm01 firm02 firm03 firm04 firm05 firm06 
#>      1      2      1      1      1      2 
```
