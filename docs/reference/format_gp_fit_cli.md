# Format a `gp_fit` for the console

Returns the plain-ASCII lines that `print.gp_fit()` emits: a composite,
one-line-per-stage view leading with the result – unit count, the
selected grade count and composition, the information-reliability
metrics `(1 - DR, tau-bar)` (read from the enriched grade-path summary
when present), and the backend. Then pointers to
[`summary()`](https://rdrr.io/r/base/summary.html) /
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md).
Recomputes nothing and never calls the solver: the metrics are read from
materialized fields, and shown as `NA` when the path summary is not
enriched.

## Usage

``` r
format_gp_fit_cli(x, ..., width = NULL)

# S3 method for class 'gp_fit'
summary(object, ...)

# S3 method for class 'gp_fit'
print(x, ..., width = getOption("width"))
```

## Arguments

- x:

  A `gp_fit`.

- ...:

  Unused.

- width:

  Optional integer rendering width.

- object:

  A `gp_fit`.

## Value

A `character` vector, one element per line.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a typed
`gp_fit_summary`.

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## See also

`print.gp_fit()`,
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)

Other gradepath-cli:
[`format_gp_frontier_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md),
[`format_gp_grade_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md),
[`format_gp_grade_path_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md),
[`format_gp_pairwise_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md),
[`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)

## Examples

``` r
# The pre-solved tiny fit bundled with the package (instant; no solve).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
writeLines(format_gp_fit_cli(fit))
#> +---------------------------------------------------------------------+
#> | gp_fit  .  24 units  .  grades: 2 (21/3)  .  selected lambda = 0.25 |
#> +---------------------------------------------------------------------+
#> | units      : 24                                                     |
#> | grades     : 2 (21/3)                                               |
#> | reliability: (1 - DR) = 0.97   tau-bar = 0.17                       |
#> | backend    : highs                                                  |
#> +---------------------------------------------------------------------+
#> i  summary(fit) for backend/selection details; get_report_card(fit) for the ranked table.
```
