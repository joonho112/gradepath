# Format a `gp_report_card` for the console

Returns the plain-ASCII lines that `print.gp_report_card()` emits,
reproducing the Chapter 12 target box. Leads with the result: unit
count, grade count and composition (e.g. `3 (2/81/14)`), selected
lambda, a HEAD of the ranked table (and the final row when truncated),
then the frozen INVARIANT 2 welfare footnote and pointers. Recomputes
nothing.

## Usage

``` r
format_gp_report_card_cli(x, ..., width = NULL, max_rows = 8L)

# S3 method for class 'gp_report_card'
summary(object, ...)

# S3 method for class 'gp_report_card'
print(x, ..., width = getOption("width"), max_rows = 8L)

# S3 method for class 'gp_report_card'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `gp_report_card`.

- ...:

  Unused.

- width:

  Optional integer rendering width.

- max_rows:

  Integer; how many ranked rows to show before an "... N more rows ..."
  line (default 8). Controls the HEAD length for large tables.

- object:

  A `gp_report_card`.

- row.names, optional:

  Passed through from
  [`base::as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html);
  `row.names` resets the payload's row names.

## Value

A `character` vector, one element per line.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a typed
`gp_report_card_summary`.

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the 10-column payload (`x$table`), the object a researcher writes to
CSV; this is the working target the print pointer advertises.

## Details

The welfare note (invariant 2) describes grade 1 as the most-extreme
theta and is NOT hard-coded to one application: it flips between firms
(most discriminatory) and names (best-treated). Rows are described as
"more discriminatory" / "higher posterior contact rate", not as a
contest.

## See also

`print.gp_report_card()`,
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md)

Other gradepath-cli:
[`format_gp_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_fit_cli.md),
[`format_gp_frontier_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_frontier_cli.md),
[`format_gp_grade_fit_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_fit_cli.md),
[`format_gp_grade_path_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_grade_path_cli.md),
[`format_gp_pairwise_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_pairwise_cli.md)

## Examples

``` r
# The bundled fit's report card; the formatter only renders materialized rows.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
card <- get_report_card(fit)
writeLines(format_gp_report_card_cli(card))
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
