# Run the one-level KRW report-card pipeline in a single call

`krw_report_card()` is the one-call entry point for the M1 one-level
Kline-Rose-Walters grading workflow. Give it firm-level estimates and a
demographic; it builds the stage-1 input through the ebrecipe seam, runs
gradepath's native one-level KRW core (beta-GMM precision seam -\>
deconvolution -\> posterior weights), then chains the public
social-choice verbs
([`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)
-\>
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
-\>
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)
-\>
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md))
into one validated `gp_fit` you can print, plot, and pass to the
accessors. `gradepath()` is an alias.

## Usage

``` r
krw_report_card(
  data,
  demographic = c("race", "gender"),
  groups = NULL,
  control = gp_control(precision_rule = "krw_gmm"),
  lambda = 0.25,
  acceptance_mode = FALSE,
  ...
)

gradepath(
  data,
  demographic = c("race", "gender"),
  groups = NULL,
  control = gp_control(precision_rule = "krw_gmm"),
  lambda = 0.25,
  acceptance_mode = FALSE,
  ...
)
```

## Arguments

- data:

  Data frame, list, or `ebrecipe::eb_estimates` carrying `theta_hat`/`s`
  (or the demographic-specific `theta_hat_race`/`se_race`,
  `theta_hat_gender`/`se_gender`) plus an optional
  `unit_id`/`firm_id`/`label`. Do not supply both a generic and a
  demographic-specific column for one quantity.

- demographic:

  Character scalar, `"race"` or `"gender"`; selects the estimate columns
  and the precision-dependence form (multiplicative for race, additive
  for gender). Defaults to `"race"`.

- groups:

  Optional grouping vector. Grouped/two-level execution is not
  implemented in this monolith; pass `NULL` (the default). A non-`NULL`
  value raises a `GROUPS_ERROR`. See
  [`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)
  for the two-level path.

- control:

  A
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  object. The monolith requires `precision_rule = "krw_gmm"` and
  defaults to that rule.

- lambda:

  Numeric scalar in `[0, 1]`; the frontier penalty at which the final
  grade is selected from the path. Defaults to `0.25` (KRW's published
  selection).

- acceptance_mode:

  Logical; solver solution-quality policy threaded to the internal
  [`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
  (gurobi backend only; ignored for open backends). `FALSE` (default)
  reports a `gap_reached`/`time_limit` solve honestly with that status
  (it lands `UNVERIFIED` downstream); `TRUE` additionally attempts later
  Gurobi paths to prove the optimum and never relabels a non-optimal
  solve as `optimal`.

- ...:

  Unused; an error is raised if any argument is passed.

## Value

A validated `gp_fit` object (a list of class `c("gp_fit", "list")`) with
the public slots:

- `ids`:

  Character vector; the one canonical unit-id order shared by every
  downstream slot (an enforced invariant).

- `estimates`:

  The stage-1 `ebrecipe::eb_estimates` input.

- `prior`:

  A `gp_prior`; the native deconvolved prior on the r-scale.

- `posterior`:

  A `gp_posterior`; per-unit posterior summaries (`posterior_mean`,
  `lower`, `upper`, `estimate`, `se`).

- `precision_fit`:

  A `gp_precision_fit`; the beta-GMM precision fit and its J-statistic
  moments.

- `pairwise`:

  A `gp_pairwise`; the J x J posterior outranking structure.

- `grade_path`:

  A `gp_grade_path`; the solved frontier of grade assignments across the
  penalty grid, with per-solve solver status.

- `selected_grade`:

  A `gp_grade_fit`; the assignment at the selected `lambda` (exactly the
  `grade_path` member at that penalty).

- `report_card`:

  A `gp_report_card`; the per-unit grade-label table with posterior
  summaries.

- `control`:

  The validated
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  used for the run.

- `provenance`:

  Named list: producer, demographic, workflow, and the honest
  `producer_status` (e.g. `ACCEPTED` / `UNVERIFIED`).

- `warnings`, `schema_version`:

  Internal audit slots.

## Details

The chain runs in order: stage-1 input (ebrecipe seam) -\> beta-GMM
precision seam -\> native deconvolution -\> posterior weights -\>
one-level posterior -\>
[`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)
-\>
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
-\>
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)
-\>
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md).
Grade labels are integers in `{1, ..., n}` and carry no
ranking-superiority statement of any kind. Solver honesty: a gurobi
solve that stops at a gap or time limit keeps that status (surfaced as
an `UNVERIFIED` producer status); `acceptance_mode = TRUE` only adds
optimization attempts, never a relabel.

## Note

Replicating KRW: do NOT pass the bundled
[krw_firms](https://joonho112.github.io/gradepath/reference/krw_firms.md)
dataset here. It is a public example on a different numeric scale and
will NOT reproduce the published KRW (2024) results – the beta-GMM lands
on a spurious large-beta optimum (race beta ~ 2.1 / gender beta ~ 3.0
versus the published 0.51 / 1.26) and the gender path errors in
deconvolution with `DECONV_BOUNDARY_ERROR`. For replication read KRW's
real Matlab GMM series (shipped under `inst/extdata/krw-gmm-input/`)
directly, e.g.
`read.csv(system.file("extdata/krw-gmm-input/theta_estimates_matlab_race.csv", package = "gradepath"), header = FALSE)`
(column 2 = `theta_hat`, column 3 = `s`); see
[krw_firms](https://joonho112.github.io/gradepath/reference/krw_firms.md)
and the applied vignettes.

## See also

[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md),
[`gp_plot_report_card()`](https://joonho112.github.io/gradepath/reference/gp_plot_report_card.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md),
[krw_firms](https://joonho112.github.io/gradepath/reference/krw_firms.md)

Other gradepath-pipeline:
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)

## Examples

``` r
# A small example input bundled with the package (24 firms; an example subset,
# not the parity input). The runnable solve uses the open HiGHS backend.
inp <- readRDS(system.file("extdata/examples/tiny_input.rds", package = "gradepath"))

# \donttest{
fit <- krw_report_card(inp, demographic = "race",
                       control = gp_control(backend = "highs",
                                            precision_rule = "krw_gmm"))
fit                       # result-first console print
#> +---------------------------------------------------------------------+
#> | gp_fit  .  24 units  .  grades: 2 (21/3)  .  selected lambda = 0.25 |
#> +---------------------------------------------------------------------+
#> | units      : 24                                                     |
#> | grades     : 2 (21/3)                                               |
#> | reliability: (1 - DR) = 0.97   tau-bar = 0.17                       |
#> | backend    : highs                                                  |
#> +---------------------------------------------------------------------+
#> i  summary(fit) for backend/selection details; get_report_card(fit) for the ranked table.
get_report_card(fit)      # the per-unit grade-label report card
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
# }
```
