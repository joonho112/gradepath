# Extract per-unit grades from a gradepath fit

Returns the selected grading as an integer vector, one entry per unit,
named by the unit ids. The grades are read from the selected
`gp_grade_fit` (`fit$selected_grade$assignment`) and aligned to the
fit's canonical `ids` order, so `get_grades(fit)[id]` is the grade of
unit `id`. Use it as the headline read of a `gp_fit`;
[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md)
is the base-generic synonym and
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)
returns the same grades alongside the posterior summaries.

## Usage

``` r
get_grades(fit, ...)

# Default S3 method
get_grades(fit, ...)

# S3 method for class 'gp_fit'
get_grades(fit, ...)
```

## Arguments

- fit:

  A materialized `gp_fit` object.

- ...:

  Unused; accepted for S3 method extensibility. No arguments are
  forwarded; supplying any has no effect.

## Value

A named `integer` vector of length `J` (the number of units): the
selected grade per unit, with `names` equal to `fit$ids`.

## Details

Per `GP-DEC-13-A` (Chapter 13 patterns) a `gp_fit` is read through the
accessor family – `get_grades()`,
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
[`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md),
[`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md)
– never by reaching into `$` slots. The accessor never recomputes; it
only reads materialized slots. Grade labels are integers in
`{1, ..., n}` and carry no ranking-superiority statement of any kind.

## See also

[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)

Other gradepath-accessors:
[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md),
[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md),
[`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
[`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md),
[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)

## Examples

``` r
# Instant: read grades off the bundled tiny fit (no solver, no Gurobi).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
head(get_grades(fit))   # named int vector: names = ids, values = grades
#> firm01 firm02 firm03 firm04 firm05 firm06 
#>      1      2      1      1      1      2 
```
