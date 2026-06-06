# Create a gradepath run-control object

`gp_control()` bundles the run settings the gradepath verbs read – the
frontier penalty grid, the grade-IP solver backend, the
precision-handling rule, the credible-interval level, and solver runtime
knobs – into one validated `gp_control` object. Every argument has a
sane default, so `gp_control()` with no arguments is a complete, valid
configuration. Build one, tweak the fields you care about, and pass it
as the `control =` argument to
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
or
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md).

## Usage

``` r
gp_control(
  lambda_grid = seq(0, 1, by = 0.01),
  backend = getOption("gradepath.backend", "gurobi"),
  precision_rule = c("none", "krw_gmm"),
  interval_level = 0.9,
  solver_options = NULL,
  time_limit = NULL,
  mip_gap = NULL,
  seed = NULL
)
```

## Arguments

- lambda_grid:

  Numeric vector; the frontier penalty values the grade path is solved
  at. Must be unique, strictly increasing, and lie in `[0, 1]`, and must
  contain the parity anchors `0.25` (KRW's published selection) and
  `1.00` (the endpoint warm-start). Defaults to the reference grid
  `seq(0, 1, by = 0.01)`; passing `NULL` restores that reference grid,
  and later grade-path execution may narrow it.

- backend:

  Character scalar; the planned grade-IP solver backend, one of
  `"gurobi"`, `"highs"`, `"glpk"`, or `"symphony"`. Defaults to the
  value of `getOption("gradepath.backend", "gurobi")` – the Gurobi
  parity backend (the same solver KRW used), with `"highs"`, `"glpk"`,
  and `"symphony"` as license-free open alternatives (`"highs"` calls
  the `highs` package directly; `"glpk"` and `"symphony"` route through
  ROI plugins). The v1 `"roi_highs"` / `"roi_glpk"` / `"roi_symphony"`
  spellings are accepted and normalized to the short names. There is no
  `"auto"`; backend availability is checked by the solver at call time,
  not here.

- precision_rule:

  Character scalar; the precision-handling rule. `"none"` keeps the raw
  estimate (r) scale; `"krw_gmm"` selects the gradepath-native KRW
  beta-GMM standardization. Defaults to `"none"`.

- interval_level:

  Numeric scalar; the credible-interval level stored on the downstream
  posterior and report-card objects. A probability in the open interval
  `(0, 1)`. Defaults to `0.90`.

- solver_options:

  Named list; backend-specific solver options. The portable keys
  `time_limit` and `mip_gap` (and their compatibility spellings
  `max_time` and `mip_rel_gap`) are range-checked in place; any other
  key passes through untouched. `mip_gap` is honored by Gurobi, HiGHS,
  and SYMPHONY but not by the current GLPK route. Defaults to an empty
  list.

- time_limit:

  Numeric scalar or `NULL`; an ergonomic top-level alias for the solver
  time limit in seconds (finite, strictly positive), merged into
  `solver_options$time_limit`. Supplying it alongside `time_limit` (or
  its `max_time` spelling) inside `solver_options` is an error. Defaults
  to `NULL` (no time limit set here).

- mip_gap:

  Numeric scalar or `NULL`; an ergonomic top-level alias for the
  relative MIP optimality gap in `[0, 1]`, merged into
  `solver_options$mip_gap`. Supplying it alongside `mip_gap` (or its
  `mip_rel_gap` spelling) inside `solver_options` is an error. Defaults
  to `NULL` (no gap target set here).

- seed:

  Integer scalar or `NULL`; a non-negative seed for the stochastic
  (two-level Monte-Carlo) paths. Defaults to `NULL` (unseeded).

## Value

A validated S3 object of class `c("gp_control", "list")`: a named list
carrying the six user-facing computational fields plus two internal
audit slots.

- `lambda_grid`:

  Numeric vector; the validated frontier penalty grid.

- `backend`:

  Character scalar; the canonical short backend name.

- `precision_rule`:

  Character scalar; `"none"` or `"krw_gmm"`.

- `interval_level`:

  Numeric scalar in `(0, 1)`; the credible-interval level.

- `solver_options`:

  Named list; the normalized solver options, with any ergonomic
  `time_limit` / `mip_gap` aliases merged in.

- `seed`:

  Non-negative integer or `NULL`; the stochastic-path seed.

- `schema_version`, `provenance`:

  Internal audit slots.

## Details

The surface is the pruned v2 control of six computational fields; v1's
`workflow`, `replication_mode`, `target_bundle`, `selection_rule`,
`backend_fallback`, and the legacy ranking vocabulary are gone.
Selection of the frontier penalty lives with the grading verbs, not
here. Portable solver runtime knobs may be passed either inside
`solver_options` or via the ergonomic top-level aliases `time_limit` and
`mip_gap`; supplying the same knob twice – including a compatibility
spelling alongside its ergonomic alias – is an error.

## See also

[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md)

Other gradepath-pipeline:
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)

## Examples

``` r
# Runnable with no arguments and no solve: the friendly defaults give backend
# "gurobi", interval_level 0.90, and the reference lambda grid.
ctl <- gp_control()

# A customized control on an open backend: a narrow grid (still carrying the
# 0.25 and 1.00 anchors) with an ergonomic time limit merged into
# solver_options.
ctl2 <- gp_control(
  lambda_grid = c(0, 0.25, 0.5, 1),
  backend     = "highs",
  time_limit  = 30
)
ctl2$solver_options$time_limit
#> [1] 30
```
