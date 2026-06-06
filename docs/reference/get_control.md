# Extract the run-control object from a gradepath fit

Returns the `gp_control` slot of a `gp_fit` unchanged: the pruned run
settings actually used for the fit (`lambda_grid`, `backend`,
`precision_rule`, `interval_level`, `solver_options`, `seed`). Reach for
it to see exactly which
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)
configuration produced the result, or to reuse those settings in a
follow-up run.

## Usage

``` r
get_control(fit, ...)

# Default S3 method
get_control(fit, ...)

# S3 method for class 'gp_fit'
get_control(fit, ...)
```

## Arguments

- fit:

  A `gp_fit` object.

- ...:

  Unused; accepted for S3 method extensibility. No arguments are
  forwarded; supplying any has no effect.

## Value

The `gp_control` object stored on `fit` (returned unchanged).

## Details

Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and
only reads the materialized slot.

## See also

[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md),
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md)

Other gradepath-accessors:
[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md),
[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md),
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
[`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md),
[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)

## Examples

``` r
# Instant: read the run-control object off the bundled tiny fit.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
get_control(fit)       # <gp_control>
#> $lambda_grid
#> [1] 0.25 1.00
#> 
#> $backend
#> [1] "highs"
#> 
#> $precision_rule
#> [1] "krw_gmm"
#> 
#> $interval_level
#> [1] 0.9
#> 
#> $solver_options
#> list()
#> 
#> $seed
#> NULL
#> 
#> $schema_version
#> [1] "v2"
#> 
#> $provenance
#> $provenance$producer
#> [1] "gp_control"
#> 
#> $provenance$built_at
#> [1] "2026-06-04 18:43:24 CDT"
#> 
#> $provenance$r_version
#> [1] "R version 4.6.0 (2026-04-24)"
#> 
#> $provenance$package_version
#> [1] "0.5.0"
#> 
#> 
#> attr(,"class")
#> [1] "gp_control" "list"      
```
