# Format a `gp_grade_fit` for the console

Returns the plain-ASCII lines that `print.gp_grade_fit()` emits. Leads
with the result: grade count, the grade-\>size table, solver status, and
lambda. Backend/channel detail is surfaced via
[`summary()`](https://rdrr.io/r/base/summary.html), not here. Recomputes
nothing.

## Usage

``` r
format_gp_grade_fit_cli(x, ..., width = NULL)

# S3 method for class 'gp_grade_fit'
summary(object, ...)

# S3 method for class 'gp_grade_fit'
print(x, ..., width = getOption("width"))
```

## Arguments

- x:

  A `gp_grade_fit`.

- ...:

  Unused.

- width:

  Optional integer rendering width.

- object:

  A `gp_grade_fit`.

## Value

A `character` vector, one element per line.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a typed
`gp_grade_fit_summary`.

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## See also

`print.gp_grade_fit()`,
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)

Other gradepath-cli:
[`format_gp_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md),
[`format_gp_frontier_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md),
[`format_gp_grade_path_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md),
[`format_gp_pairwise_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md),
[`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)

## Examples

``` r
# The bundled tiny fit; its selected slice is a gp_grade_fit (instant read).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
writeLines(format_gp_grade_fit_cli(fit$selected_grade))
#> +---------------------------------------------+
#> | gp_grade_fit  .  2 grades  .  lambda = 0.25 |
#> +---------------------------------------------+
#> | units = 24   status = optimal               |
#> | grade 1: 21                                 |
#> | grade 2: 3                                  |
#> +---------------------------------------------+
#> i  assignment: x$assignment   summary(x) for backend/channel details
```
