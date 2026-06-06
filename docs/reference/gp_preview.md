# Preview a gradepath run without solving

`gp_preview()` is the solve-free dry run for the one-level KRW workflow.
Give it the same firm-level input and demographic you would hand to
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md);
it validates the input grammar and reports the inferred plan – the unit
ids, the inferred workflow, the chosen backend, the operational frontier
penalty grid, and the number of grade-IP solves the run would perform –
then returns a `gp_preview` object WITHOUT running the beta-GMM,
deconvolution, posterior, pairwise, or any solver stage. Reach for it to
size and sanity-check a run before committing to a full solve.

## Usage

``` r
gp_preview(
  data,
  demographic = c("race", "gender"),
  groups = NULL,
  control = gp_control(),
  ...
)
```

## Arguments

- data:

  Data frame, list, or `ebrecipe::eb_estimates` object carrying
  `theta_hat`/`s` (or the demographic-specific
  `theta_hat_race`/`se_race`, `theta_hat_gender`/`se_gender`) plus an
  optional `unit_id`/`firm_id`/`label`. Do not supply both a generic and
  a demographic-specific column for the same quantity in one input.

- demographic:

  Character scalar, `"race"` or `"gender"`; selects the estimate columns
  the input grammar reads. Defaults to `"race"`.

- groups:

  Optional grouping vector, or `NULL`. Grouped/two-level execution is
  not implemented in the M1 monolith; a non-`NULL` value is previewed as
  a pending grouped workflow (with an explanatory warning and a
  `GROUPS_ERROR` status), not run. Defaults to `NULL` (the one-level
  workflow).

- control:

  A
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  object, or `NULL` to use the defaults. Supplies the backend and the
  lambda grid the preview reports. Defaults to
  [`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md).

- ...:

  Unused; an error is raised if any argument is passed.

## Value

A validated S3 object of class `c("gp_preview", "list")`: the planned
run with no results.

- `ids`:

  Character vector; the canonical unit-id order parsed from the input.

- `n_units`:

  Integer scalar; the number of units (equals `length(ids)`).

- `demographic`:

  Character scalar; `"race"` or `"gender"`.

- `groups`:

  Character vector of length `n_units`, or `NULL`; the parsed grouping
  when one was supplied.

- `workflow`:

  Character scalar; the inferred workflow, either
  `"one_level_independence"` or `"grouped_pending"`.

- `lambda_grid`:

  Numeric vector; the operational frontier penalty grid the run would
  sweep.

- `estimated_solves`:

  Integer scalar; the number of grade-IP solves the run would perform
  (the length of `lambda_grid`).

- `backend`:

  Character scalar; the grade-IP backend the run would use.

- `control`:

  The validated
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  the preview was built from.

- `status`:

  A `gp_status`; `OK` for a one-level plan, or `GROUPS_ERROR` when a
  grouping was supplied.

- `schema_version`, `provenance`, `warnings`:

  Internal audit slots; the `provenance` stamp records that no solve was
  performed.

## See also

[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)

Other gradepath-grade:
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)

## Examples

``` r
# Solve-free: load the bundled tiny example input and preview the planned run
# (no Gurobi, no solve, instant).
inp <- readRDS(system.file("extdata/examples/tiny_input.rds", package = "gradepath"))
pv <- gp_preview(inp, demographic = "race")
pv                  # result-first console print of the planned run
#> $ids
#>  [1] "firm01" "firm02" "firm03" "firm04" "firm05" "firm06" "firm07" "firm08"
#>  [9] "firm09" "firm10" "firm11" "firm12" "firm13" "firm14" "firm15" "firm16"
#> [17] "firm17" "firm18" "firm19" "firm20" "firm21" "firm22" "firm23" "firm24"
#> 
#> $n_units
#> [1] 24
#> 
#> $demographic
#> [1] "race"
#> 
#> $groups
#> NULL
#> 
#> $workflow
#> [1] "one_level_independence"
#> 
#> $lambda_grid
#>  [1] 0.00 0.10 0.20 0.25 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00
#> 
#> $estimated_solves
#> [1] 12
#> 
#> $backend
#> [1] "gurobi"
#> 
#> $control
#> $lambda_grid
#>   [1] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14
#>  [16] 0.15 0.16 0.17 0.18 0.19 0.20 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29
#>  [31] 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44
#>  [46] 0.45 0.46 0.47 0.48 0.49 0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59
#>  [61] 0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74
#>  [76] 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
#>  [91] 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.00
#> 
#> $backend
#> [1] "gurobi"
#> 
#> $precision_rule
#> [1] "none"
#> 
#> $interval_level
#> [1] 0.9
#> 
#> $solver_options
#> list()
#> 
#> $seed
#> NULL
#> 
#> $schema_version
#> [1] "v2"
#> 
#> $provenance
#> $provenance$producer
#> [1] "gp_control"
#> 
#> $provenance$built_at
#> [1] "2026-06-06 08:37:51 CDT"
#> 
#> $provenance$r_version
#> [1] "R version 4.6.0 (2026-04-24)"
#> 
#> $provenance$package_version
#> [1] "0.5.0"
#> 
#> 
#> attr(,"class")
#> [1] "gp_control" "list"      
#> 
#> $status
#> $status
#> [1] "OK"
#> 
#> $message
#> [1] "Preview OK"
#> 
#> $detail
#> list()
#> 
#> attr(,"class")
#> [1] "gp_status" "list"     
#> 
#> $schema_version
#> [1] "v2"
#> 
#> $provenance
#> $provenance$producer
#> [1] "gp_preview"
#> 
#> $provenance$no_solve
#> [1] TRUE
#> 
#> 
#> $warnings
#> character(0)
#> 
#> attr(,"class")
#> [1] "gp_preview" "list"      
pv$n_units          # number of units the run would grade
#> [1] 24
pv$estimated_solves # number of grade-IP solves the run would perform
#> [1] 12
```
