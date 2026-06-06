# Format a `gp_frontier` for the console

Returns the plain-ASCII lines that `print.gp_frontier()` emits,
reproducing the Chapter 14 worked-plot box. Leads with the result: the
selected lambda and its `(1 - DR, tau-bar, grade_count)`, the `(1 - DR)`
and `tau-bar` ranges across the path, and how many naive benchmarks are
available. Recomputes nothing.

## Usage

``` r
format_gp_frontier_cli(x, ..., width = NULL)

# S3 method for class 'gp_frontier'
summary(object, ...)

# S3 method for class 'gp_frontier'
print(x, ..., width = getOption("width"))
```

## Arguments

- x:

  A `gp_frontier`.

- ...:

  Unused.

- width:

  Optional integer rendering width.

- object:

  A `gp_frontier`.

## Value

A `character` vector, one element per line.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a typed
`gp_frontier_summary`.

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## See also

`print.gp_frontier()`,
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md)

Other gradepath-cli:
[`format_gp_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md),
[`format_gp_grade_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md),
[`format_gp_grade_path_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md),
[`format_gp_pairwise_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md),
[`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)

## Examples

``` r
# Build the frontier from the bundled fit's grade path + pairwise (instant).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
fr <- gp_frontier(fit$grade_path, pairwise = get_pairwise(fit))
writeLines(format_gp_frontier_cli(fr))
#> +-------------------------------------------------------------+
#> | gp_frontier  .  2 lambda values  .  selected lambda = 0.25  |
#> +-------------------------------------------------------------+
#> | selected:  (1 - DR) = 0.97   tau-bar = 0.17   grades = 2    |
#> | range:     1 - DR in [0.72, 0.97]   tau-bar in [0.17, 0.43] |
#> | benchmarks: 0 naive overlays available (see gp_plot_*)      |
#> +-------------------------------------------------------------+
#> i  gp_plot_frontier(fr) draws the (1 - DR, tau-bar) frontier.
```
