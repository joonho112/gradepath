# Extract the pairwise outranking object from a gradepath fit

Returns the `gp_pairwise` slot of a `gp_fit` unchanged: the \\J \times
J\\ posterior outranking matrix \\\pi\_{ij}\\ (the posterior probability
that unit `i` outranks unit `j`) together with its cleanup and source
metadata. Reach for it to inspect the social-choice input that the grade
path is solved over.

## Usage

``` r
get_pairwise(fit, ...)

# Default S3 method
get_pairwise(fit, ...)

# S3 method for class 'gp_fit'
get_pairwise(fit, ...)
```

## Arguments

- fit:

  A `gp_fit` object.

- ...:

  Unused; accepted for S3 method extensibility. No arguments are
  forwarded; supplying any has no effect.

## Value

The `gp_pairwise` object stored on `fit` (returned unchanged).

## Details

Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and
only reads the materialized slot. The matrix encodes pairwise posterior
order only and carries no ranking-superiority statement of any kind.

## See also

[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md),
[`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)

Other gradepath-accessors:
[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md),
[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md),
[`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md),
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md),
[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)

## Examples

``` r
# Instant: read the pairwise object off the bundled tiny fit.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
get_pairwise(fit)      # <gp_pairwise>
#> +---------------------------------------------------+
#> | gp_pairwise  .  24 x 24  .  552 ordered pairs     |
#> +---------------------------------------------------+
#> | pi range (off-diag): [0.009, 0.991]   diag = 0.50 |
#> | power = 0   rule = outer_product                  |
#> +---------------------------------------------------+
#> i  matrix: as.matrix(pw)   ids: pw$ids
```
