# Assemble a unit-level report card

`gp_report_card()` assembles the per-unit report-card table: each unit
with its grade label, posterior summary, and original estimate, sorted
from most-extreme to least-extreme. It does not solve a grade problem –
the row order comes from the endpoint fit cached in `grade_path`, while
the displayed `grade` column comes from `selected_grade`. Pass a
finished `gp_fit` to read back its stored card, or the stage objects to
assemble one yourself; the result prints, plots with
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
and is what
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)
returns.

## Usage

``` r
gp_report_card(
  estimates = NULL,
  posterior = NULL,
  selected_grade = NULL,
  grade_path = NULL,
  ...
)
```

## Arguments

- estimates:

  One of: an `ebrecipe::eb_estimates`-like input container, a `gp_fit`
  (its stored report card is returned as-is), a `gp_grade_path`, or a
  `gp_grade_fit`. For stage-wise assembly it supplies the report-card
  `label`, `estimate`, and `se`; if omitted (`NULL`, the default) those
  columns are read from `posterior`.

- posterior:

  A `gp_posterior` (from
  [`get_posterior()`](https://joonho112.github.io/gradepath/reference/get_posterior.md))
  giving each unit's `posterior_mean`, `lower`, and `upper`. Required
  for stage-wise assembly; ignored when `estimates` is a `gp_fit`.

- selected_grade:

  A `gp_grade_fit`, usually from
  [`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md);
  supplies the displayed `grade` column and the selected penalty.
  Required for stage-wise assembly.

- grade_path:

  A solved `gp_grade_path`; its endpoint fit fixes the row sort order
  (most-extreme grade first, ties broken by id). Required for stage-wise
  assembly.

- ...:

  Unused; reserved for future arguments. Passing any value raises an
  error.

## Value

A validated `gp_report_card` object (a list of class
`c("gp_report_card", "list")`) with the public slots:

- `ids`:

  Character vector; the unit-id order of the table rows (most-extreme
  grade first).

- `table`:

  Data frame, one row per unit, with `id`, `label`, the integer `grade`,
  `sort_rank`, `selected_lambda`, the posterior summary
  `posterior_mean`/`lower`/`upper`, and the original `estimate`/`se`.

- `selected_lambda`:

  Numeric; the penalty the displayed grades were taken from.

- `grades`:

  Integer vector of the per-unit grade labels, aligned to `ids`.

- `control`:

  The
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  carried from `grade_path`.

- `provenance`, `warnings`, `schema_version`:

  Producer metadata and internal audit slots.

## Details

Grade labels are integers in `{1, ..., n}` and carry no
ranking-superiority statement of any kind: grade 1 is simply the
most-extreme theta block.

## Note

Welfare orientation (frozen INVARIANT 2). Grade 1 is the most-extreme
theta, and its welfare reading FLIPS with the application: for firms the
most-extreme block is the most discriminatory, while for names it is the
best-treated. The label is therefore not hard-coded to one application;
the console print and
[`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)
surface this same note. See Appendix A.

## See also

[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)

## Examples

``` r
# Read the stored card from the bundled fit (no Gurobi, no solve).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
gp_report_card(fit)              # returns fit$report_card
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

# Re-assemble the same card from the fit's stage slots (solve-free).
rc <- gp_report_card(
  estimates      = fit$estimates,
  posterior      = fit$posterior,
  selected_grade = fit$selected_grade,
  grade_path     = fit$grade_path
)
head(rc$table)                   # per-unit grade-label table
#>       id  label grade sort_rank selected_lambda posterior_mean       lower
#> 1 firm03 firm03     1         1            0.25      0.1330361 0.004466502
#> 2 firm05 firm05     1         2            0.25      0.1547879 0.094098498
#> 3 firm10 firm10     1         3            0.25      0.1372338 0.024811477
#> 4 firm18 firm18     1         4            0.25      0.1166670 0.003301362
#> 5 firm24 firm24     1         5            0.25      0.1206680 0.004841016
#> 6 firm01 firm01     1         6            0.25      0.1323136 0.080808125
#>       upper estimate     se
#> 1 0.2014888   0.3431 0.2386
#> 2 0.1890643   0.4332 0.1309
#> 3 0.1749419   0.2895 0.1142
#> 4 0.1895925   0.1811 0.1902
#> 5 0.1782374   0.2017 0.1398
#> 6 0.1631853   0.2355 0.0838
```
