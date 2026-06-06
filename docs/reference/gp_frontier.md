# Summarize the frontier of a solved grade path

`gp_frontier()` collects the frontier metrics for a solved grade path:
the per-penalty table of grade count, discordance rate, reliability, and
average rank distance, and – when you also supply the pairwise
outranking profile – the five naive benchmark scorings and the
selected-penalty conditional discordance matrix. Reach for it after
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
when you want the frontier as a tidy object to print, plot with
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
or compare grades against simple baselines.

## Usage

``` r
gp_frontier(
  grade_path,
  pairwise = NULL,
  selected_lambda = NULL,
  benchmark_scores = NULL
)
```

## Arguments

- grade_path:

  A solved `gp_grade_path` (from
  [`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
  or the `grade_path` slot of a `gp_fit`). Supplies the penalty grid and
  the per-penalty grade assignments the frontier is built over.

- pairwise:

  Optional `gp_pairwise` (from
  [`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)
  or
  [`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md)).
  When `NULL` (default) `gp_frontier()` returns only the frontier table
  cached by
  [`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md);
  the benchmarks and conditional discordance matrix are left empty.
  Supplying it must match `grade_path$ids` and unlocks both. Its `ids`
  must equal `grade_path$ids`.

- selected_lambda:

  Optional numeric penalty at which the conditional discordance matrix
  is evaluated. `NULL` (default) uses
  `grade_path$selection$selected_lambda`; any value supplied must match
  a penalty on the grid exactly.

- benchmark_scores:

  Optional named list (or one-row-per-unit data frame) of the per-unit
  scores the naive baselines are built from. Requires three numeric
  vectors named `raw_estimate`, `posterior_mean`, and
  `linear_shrinkage`, each one value per unit and named by (or
  column-aligned to) the unit ids in `pairwise$ids`; they are reordered
  to that id order. From these, five baseline grade scorings are formed
  – a rank cut on each of the three vectors plus decile and quartile
  bins of `posterior_mean` – and scored against the pairwise profile.
  `NULL` (default) leaves the benchmark table empty. Requires
  `pairwise`; passing scores without it raises an error.

## Value

A validated `gp_frontier` object (a list of class
`c("gp_frontier", "list")`) with the public slots:

- `ids`:

  Character vector; the canonical unit-id order, equal to
  `grade_path$ids`.

- `table`:

  Data frame, one row per penalty on the grid, with `lambda`,
  `grade_count`, `discordance_rate`, `reliability`, and `tau_bar` (the
  average between-grade rank distance), plus the solve
  `objective`/`status` when a `pairwise` profile is supplied.

- `benchmarks`:

  Data frame of the five naive baseline scorings (`raw_estimate_rank`,
  `posterior_mean_rank`, `linear_shrinkage_rank`,
  `posterior_mean_decile`, `posterior_mean_quartile`), each with its
  `grade_count`, `discordance_rate`, `reliability`, and `tau_bar`; empty
  when `benchmark_scores` (and `pairwise`) were not supplied.

- `dr_matrix`:

  Lower-triangular matrix of conditional discordance rates between grade
  pairs at the selected penalty (row = less-extreme grade, column =
  more-extreme grade; diagonal and upper triangle are `NA`); a 0 x 0
  matrix when no `pairwise` profile was supplied.

- `selection`:

  Named list recording `selected_lambda`, the penalty the conditional
  discordance matrix was evaluated at.

- `control`:

  The
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  carried from `grade_path`.

- `provenance`, `warnings`, `schema_version`:

  Producer metadata and internal audit slots; `warnings` notes when
  benchmarks and the conditional discordance matrix were skipped for
  lack of a `pairwise` profile.

## Details

The frontier metrics are computed over between-grade blocks of the
pairwise outranking profile and normalized by `2 * choose(n, 2)`. The
`dr_matrix` reports, for each pair of grades, the share of cross-grade
comparisons that run against the grade ordering (discordant),
conditional on that pair. Grade labels are integers and carry no
ranking-superiority statement; the benchmarks are descriptive baselines,
not a contest.

## See also

[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md),
[`gp_plot_frontier()`](https://joonho112.github.io/gradepath/reference/gp_plot_frontier.md),
[`gp_plot_discordance()`](https://joonho112.github.io/gradepath/reference/gp_plot_discordance.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md)

Other gradepath-frontier:
[`gp_krw_r2()`](https://joonho112.github.io/gradepath/reference/gp_krw_r2.md),
[`gp_r2()`](https://joonho112.github.io/gradepath/reference/gp_r2.md)

## Examples

``` r
# Inspect the frontier of the bundled pre-solved fit (no Gurobi, no solve).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
post <- get_posterior(fit)
ids <- get_pairwise(fit)$ids

fr <- gp_frontier(
  fit$grade_path,
  pairwise = get_pairwise(fit),
  benchmark_scores = list(
    raw_estimate     = setNames(post$estimate, ids),
    posterior_mean   = setNames(post$posterior_mean, ids),
    linear_shrinkage = setNames(0.5 * post$estimate, ids)
  )
)
fr$table        # per-penalty frontier metrics
#>   lambda grade_count discordance_rate reliability   tau_bar  objective
#> 1   0.25           2       0.02813416   0.9718658 0.1719926  -12.08743
#> 2   1.00          24       0.28414803   0.7158520 0.4317039 -238.30057
#>   raw_objective  status
#> 1     -6.043715 optimal
#> 2   -119.150287 optimal
fr$benchmarks   # the five naive baseline scorings
#>                 benchmark grade_count discordance_rate reliability   tau_bar
#> 1       raw_estimate_rank          24        0.3177080   0.6822920 0.3645841
#> 2     posterior_mean_rank          24        0.3204890   0.6795110 0.3590219
#> 3   linear_shrinkage_rank          24        0.3177080   0.6822920 0.3645841
#> 4   posterior_mean_decile          10        0.2903434   0.7096566 0.3540959
#> 5 posterior_mean_quartile           4        0.2357119   0.7642881 0.3111849
```
