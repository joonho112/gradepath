# Solve the grade path over a lambda grid

`gp_grade_path()` solves one grade integer program per penalty value on
a lambda grid, tracing the frontier between fewer grades and lower
discordance. It is the engine the one-level pipeline calls; reach for it
directly to obtain the whole solved path, then pull a single assignment
with
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)
(solve-free) or solve one penalty in isolation with
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md).
The Condorcet endpoint (`lambda = 1`) is solved first as a warm-start
anchor, the remaining penalties are solved in descending order, and fits
are stored back in ascending grid order.

## Usage

``` r
gp_grade_path(
  pairwise,
  lambda_grid = NULL,
  control = NULL,
  selected_lambda = 0.25,
  selection_rule = "baseline_lambda_0.25",
  acceptance_mode = FALSE
)
```

## Arguments

- pairwise:

  A
  [gp_pairwise](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)
  object: the `J x J` posterior outranking structure the integer program
  scores.

- lambda_grid:

  Numeric vector of penalties in `[0, 1]`, or `NULL`. `NULL` (default)
  uses the operational small grid unless `pairwise` or `control` already
  carries a custom grid. An explicit grid must include the parity
  anchors `0.25` and `1`.

- control:

  Optional
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md);
  the backend, solver tolerances, and active lambda grid. Defaults to
  `pairwise$control`.

- selected_lambda:

  Numeric scalar in the grid; the penalty recorded as the default
  reportable selection on the returned path. Defaults to `0.25` (KRW's
  published selection). Recording it never triggers a solve.

- selection_rule:

  Character scalar; a label stored alongside `selected_lambda` in the
  path's selection metadata. Defaults to `"baseline_lambda_0.25"`.

- acceptance_mode:

  Logical; the solver solution-quality policy, honored by the gurobi
  backend only (open backends solve a single invocation, so it is
  ignored there). `FALSE` (default) keeps the first Gurobi path that
  returns a normalized result, so a selected solve that stops at
  `gap_reached` / `time_limit` is reported honestly with that status (it
  lands `UNVERIFIED` downstream). `TRUE` additionally attempts later
  Gurobi paths to prove the optimum when an earlier path is not
  acceptance-ready; it stays honest and never relabels a non-optimal
  solve as `optimal`. The chosen value is recorded on each fit's
  `backend$acceptance_mode`.

## Value

A validated gp_grade_path object (a list of class
`c("gp_grade_path", "list")`) with the public slots:

- `ids`:

  Character vector; the canonical unit-id order shared by every fit on
  the path.

- `lambda_grid`:

  Numeric vector; the strictly increasing penalty grid that was solved,
  including the `0.25` and `1` anchors.

- `fits`:

  List of `gp_grade_fit`, one per `lambda_grid` value in ascending
  order; the per-penalty grade assignments and solver status.

- `summary`:

  Data frame, one row per penalty, with `lambda`, `grade_count`,
  `discordance_rate`, `reliability`, `tau_bar`, `objective`,
  `raw_objective`, and `status` – the frontier table.

- `backend`:

  Named list; the resolved solver name, paths, and encodings used across
  the solves.

- `selection`:

  Named list (`selected_lambda`, `selection_rule`, `endpoint_lambda`);
  the recorded default selection, not a re-solve.

- `control`:

  The validated
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  used for the run.

- `warnings`, `schema_version`, `provenance`:

  Internal audit slots (solve order, warm-start strategy, accumulated
  warnings).

## Details

Grade labels are contiguous integers in `{1, ..., k}` and carry no
ranking-superiority statement of any kind. Solving descends from the
`lambda = 1` endpoint so each gurobi solve can warm-start from the
adjacent penalty; open backends solve each penalty independently.
Selecting a fit later (with
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md))
only reads the stored list and never re-solves. See KRW (2024, Section
4) and the grade-engine method vignette for the integer program and the
frontier penalty.

## See also

[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md),
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md),
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)

Other gradepath-grade:
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md),
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md),
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)

## Examples

``` r
# Instant: read the solved path from the bundled tiny fit (24 firms).
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
fit$grade_path                 # the pre-solved gp_grade_path
#> +--------------------------------------------------------------+
#> | gp_grade_path  .  2 lambda values  .  selected lambda = 0.25 |
#> +--------------------------------------------------------------+
#> | grades across path in [2, 24]                                |
#> | endpoint lambda = 1.00                                       |
#> +--------------------------------------------------------------+
#> i  path table: x$summary   per-lambda fits: x$fits
fit$grade_path$summary         # the frontier table over the penalty grid
#>   lambda grade_count discordance_rate reliability   tau_bar  objective
#> 1   0.25           2       0.02813416   0.9718658 0.1719926  -12.08743
#> 2   1.00          24       0.28414803   0.7158520 0.4317039 -238.30057
#>   raw_objective  status
#> 1     -6.043715 optimal
#> 2   -119.150287 optimal

# Live solve from the stored pairwise structure; open HiGHS backend, small grid.
# \donttest{
path <- gp_grade_path(
  fit$pairwise,
  control = gp_control(backend = "highs", lambda_grid = c(0.25, 1))
)
gp_select_grade(path, lambda = 0.25)
#> +---------------------------------------------+
#> | gp_grade_fit  .  2 grades  .  lambda = 0.25 |
#> +---------------------------------------------+
#> | units = 24   status = optimal               |
#> | grade 1: 21                                 |
#> | grade 2: 3                                  |
#> +---------------------------------------------+
#> i  assignment: x$assignment   summary(x) for backend/channel details
# }
```
