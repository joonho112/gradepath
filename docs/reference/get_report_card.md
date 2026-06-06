# Extract the report card from a gradepath fit

Returns the `gp_report_card` slot of a `gp_fit` unchanged – the
endpoint-sorted per-unit table of grade label, posterior mean, and
credible interval, produced by
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)
inside the pipeline. Reach for it when you want the full presentation
table;
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md)
returns just the grade vector and the base-generic
[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md)
returns this card's underlying data frame.

## Usage

``` r
get_report_card(fit, ...)

# Default S3 method
get_report_card(fit, ...)

# S3 method for class 'gp_fit'
get_report_card(fit, ...)
```

## Arguments

- fit:

  A `gp_fit` object.

- ...:

  Unused; accepted for S3 method extensibility. No arguments are
  forwarded; supplying any has no effect.

## Value

The `gp_report_card` object stored on `fit` (returned unchanged).

## Details

Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and
only reads the materialized slot. Grade labels are integers in
`{1, ..., n}` and carry no ranking-superiority statement of any kind.

## See also

[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md),
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md),
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)

Other gradepath-accessors:
[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md),
[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md),
[`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md),
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
[`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md),
[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md)

## Examples

``` r
# Instant: read the report-card object off the bundled tiny fit.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
get_report_card(fit)   # <gp_report_card>
#> <gp_report_card>  unit: unit | 24 rows | grades: 2 (21/3) | selected lambda: 0.25
#> sorted by Condorcet rank (endpoint lambda = 1; secondary key: id)
#> 
#>   sort_rank  id      label   grade  posterior_mean  CI (90%)        estimate
#>    1         firm03  firm03  1      0.133           [0.004, 0.201]  0.343   
#>    2         firm05  firm05  1      0.155           [0.094, 0.189]  0.433   
#>    3         firm10  firm10  1      0.137           [0.025, 0.175]  0.289   
#>    4         firm18  firm18  1      0.117           [0.003, 0.190]  0.181   
#>    5         firm24  firm24  1      0.121           [0.005, 0.178]  0.202   
#>    6         firm01  firm01  1      0.132           [0.081, 0.163]  0.235   
#>    7         firm15  firm15  1      0.107           [0.003, 0.171]  0.131   
#>   ... 16 more rows ...
#>   24         firm02  firm02  2      0.007           [0.000, 0.022]  -0.032  
#> 
#> i  grade 1 = most-extreme theta (firms: most discriminatory; names: best-treated). See ?gp_report_card and Appendix A.
#> i  full table: as.data.frame(card)   formatted CLI: format_gp_report_card_cli(card)
```
