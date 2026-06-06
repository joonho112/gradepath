# Format a `gp_grade_path` for the console

Returns the plain-ASCII lines that `print.gp_grade_path()` emits. Leads
with the result: number of stored lambda values, the selected lambda,
the grade-count range across the path, and the endpoint lambda.
Recomputes nothing.

## Usage

``` r
format_gp_grade_path_cli(x, ..., width = NULL)

# S3 method for class 'gp_grade_path'
summary(object, ...)

# S3 method for class 'gp_grade_path'
print(x, ..., width = getOption("width"))
```

## Arguments

- x:

  A `gp_grade_path`.

- ...:

  Unused.

- width:

  Optional integer rendering width.

- object:

  A `gp_grade_path`.

## Value

A `character` vector, one element per line.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a typed
`gp_grade_path_summary`.

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## See also

`print.gp_grade_path()`,
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)

Other gradepath-cli:
[`format_gp_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md),
[`format_gp_frontier_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md),
[`format_gp_grade_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md),
[`format_gp_pairwise_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md),
[`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)

## Examples

``` r
# The bundled tiny fit carries a solved grade path (instant read; no solve).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
writeLines(format_gp_grade_path_cli(fit$grade_path))
#> +--------------------------------------------------------------+
#> | gp_grade_path  .  2 lambda values  .  selected lambda = 0.25 |
#> +--------------------------------------------------------------+
#> | grades across path in [2, 24]                                |
#> | endpoint lambda = 1.00                                       |
#> +--------------------------------------------------------------+
#> i  path table: x$summary   per-lambda fits: x$fits
```
