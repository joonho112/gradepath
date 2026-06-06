# Format a `gp_pairwise` for the console

Returns the plain-ASCII lines that `print.gp_pairwise()` emits. Leads
with the result: matrix dimensions (`J x J`), the off-diagonal pi range,
the number of ordered pairs, and the power. Recomputes nothing.

## Usage

``` r
format_gp_pairwise_cli(x, ..., width = NULL)

# S3 method for class 'gp_pairwise'
summary(object, ...)

# S3 method for class 'gp_pairwise'
print(x, ..., width = getOption("width"))
```

## Arguments

- x:

  A `gp_pairwise`.

- ...:

  Unused.

- width:

  Optional integer rendering width (defaults to `getOption("width")`).

- object:

  A `gp_pairwise`.

## Value

A `character` vector, one element per line.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a typed
`gp_pairwise_summary`.

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## See also

`print.gp_pairwise()`,
[`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md)

Other gradepath-cli:
[`format_gp_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md),
[`format_gp_frontier_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md),
[`format_gp_grade_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md),
[`format_gp_grade_path_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md),
[`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)

## Examples

``` r
# The bundled tiny fit; the formatter only reads materialized fields (instant).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
pw <- get_pairwise(fit)
writeLines(format_gp_pairwise_cli(pw))
#> +---------------------------------------------------+
#> | gp_pairwise  .  24 x 24  .  552 ordered pairs     |
#> +---------------------------------------------------+
#> | pi range (off-diag): [0.009, 0.991]   diag = 0.50 |
#> | power = 0   rule = outer_product                  |
#> +---------------------------------------------------+
#> i  matrix: as.matrix(pw)   ids: pw$ids
```
