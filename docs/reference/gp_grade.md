# Solve a single-penalty grade fit

`gp_grade()` solves the grade integer program at exactly one penalty and
returns a `gp_grade_fit`. It is a targeted one-shot wrapper for
diagnostics and small workflows; for the whole frontier use
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
and to read a penalty already on a solved path use
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)
(which never solves).

## Usage

``` r
gp_grade(pairwise, lambda = 0.25, control = NULL, acceptance_mode = FALSE)
```

## Arguments

- pairwise:

  A
  [gp_pairwise](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)
  object: the `J x J` posterior outranking structure the integer program
  scores.

- lambda:

  Numeric scalar; the penalty to solve, which must be present in the
  active `control` lambda grid. Defaults to the KRW baseline `0.25`.

- control:

  Optional
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md);
  the backend, solver tolerances, and active lambda grid. Defaults to
  `pairwise$control`.

- acceptance_mode:

  Logical; the solver solution-quality policy, honored by the gurobi
  backend only (ignored for open backends). `FALSE` (default) keeps the
  first Gurobi path that returns a normalized result, so a solve that
  stops at `gap_reached` / `time_limit` is reported honestly with that
  status (it lands `UNVERIFIED` downstream). `TRUE` additionally
  attempts later Gurobi paths to prove the optimum when an earlier path
  is not acceptance-ready; it stays honest and never relabels a
  non-optimal solve as `optimal`. The chosen value is recorded on the
  fit's `backend$acceptance_mode`.

## Value

A validated `gp_grade_fit` object (a list of class
`c("gp_grade_fit", "list")`) with the public slots:

- `ids`:

  Character vector; the canonical unit-id order of the assignment.

- `lambda`:

  Numeric scalar; the penalty this fit was solved at.

- `assignment`:

  Data frame with one row per id: `id` (character) and `grade`
  (contiguous integer label in `{1, ..., k}`).

- `summary`:

  Named list with `grade_count`, `status`, and `n_units` for the
  realized assignment.

- `objective`:

  Named list; the penalized and raw objective values at this penalty.

- `backend`:

  Named list; the solver `name`, normalized `status`, `mipgap`,
  `runtime`, `acceptance_mode`, and the problem hash/signature.

- `control`:

  The validated
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  used for the solve.

- `warnings`, `schema_version`, `provenance`:

  Internal audit slots.

## Details

Grade labels are contiguous integers in `{1, ..., k}` and carry no
ranking-superiority statement of any kind. Solver honesty matches
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md):
a gurobi solve that stops at a gap or time limit keeps that status, and
`acceptance_mode = TRUE` only adds optimization attempts, never a
relabel. See KRW (2024, Section 4) and the grade-engine method vignette.

## See also

[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md),
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)

Other gradepath-grade:
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md),
[`gp_select_grade()`](https://joonho112.github.io/gradepath/reference/gp_select_grade.md)

## Examples

``` r
# Instant: the bundled tiny fit already stores the baseline-penalty grade fit.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
fit$selected_grade                 # the gp_grade_fit at lambda = 0.25
#> +---------------------------------------------+
#> | gp_grade_fit  .  2 grades  .  lambda = 0.25 |
#> +---------------------------------------------+
#> | units = 24   status = optimal               |
#> | grade 1: 21                                 |
#> | grade 2: 3                                  |
#> +---------------------------------------------+
#> i  assignment: x$assignment   summary(x) for backend/channel details
gp_select_grade(fit$grade_path, lambda = 0.25)  # the same object, no solve
#> +---------------------------------------------+
#> | gp_grade_fit  .  2 grades  .  lambda = 0.25 |
#> +---------------------------------------------+
#> | units = 24   status = optimal               |
#> | grade 1: 21                                 |
#> | grade 2: 3                                  |
#> +---------------------------------------------+
#> i  assignment: x$assignment   summary(x) for backend/channel details

# Live single-penalty solve from the stored pairwise; open HiGHS backend.
# \donttest{
gp_grade(
  fit$pairwise,
  lambda = 0.25,
  control = gp_control(backend = "highs", lambda_grid = c(0.25, 1))
)
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
