# A5: Solvers and calibration

Abstract

This is the practical-operations vignette — the last of the applied
track, for running gradepath day to day rather than reading its output.
It covers four things. First, choosing a grade integer-program backend:
Gurobi is the recommended default and by far the fastest here, while the
license-free HiGHS is the open fallback. Second, sizing a run before you
commit to a minutes-long solve, with the solve-free gp_preview(). Third,
the seeded Monte-Carlo gp_calibrate() harness — what it checks, the
budget heads-up it prints, and its honest scope as an in-package
diagnostic rather than a Table E1 reproduction. Fourth, reading the
package’s sharpened error messages so a mistake names exactly what went
wrong. The package installs, loads, and runs its tests with no Gurobi
license; only an actual grade solve needs one.

## 1. Running gradepath in practice

The earlier applied vignettes were about *output*: load a fit, read its
grades, draw its figures. This one is about *operations* — the four
things you actually manage when you run gradepath on your own input
rather than off a bundled fit:

- **The solver.** The grade step is an integer program, and it needs a
  backend. Which one, and what does each cost you?
- **Sizing.** A real 97-firm solve takes minutes. How do you confirm a
  run is well-formed and learn how big it is *before* paying for it?
- **Calibration.** Does the grading method’s frequentist behaviour hold
  up under simulation, and how far can you trust the check?
- **Errors.** When you assemble something by hand and get it wrong, does
  the package tell you exactly what is missing?

Nothing in this vignette solves a grade problem on render. The
solve-free pieces —
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md),
building a control, the sharpened-error demos — run live; the solving
pieces —
[`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md),
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
— are shown but not executed, exactly as the package documents. For the
live demos we use a real bundled fit and the bundled tiny raw input:

``` r

fit <- gp_parity_fit("race")
inp <- readRDS(system.file("extdata/examples/tiny_input.rds", package = "gradepath"))
```

`fit` is the proven 97-firm RACE report card — grades `2 / 81 / 14` at
the KRW baseline `lambda = 0.25` (Kline et al. 2024); `inp` is a tiny
raw input, the per-unit `theta_hat` and `s` series, that we will preview
without solving. Throughout, the grades are a *sorting by posterior
evidence*, never a contest: grade 1 is simply the most-extreme block,
and the one-level race and gender core is accepted within its declared
scope.

## 2. Choosing a solver

gradepath grades units by solving an integer program, and that program
needs a backend solver. The package follows Kline, Rose & Walters (Kline
et al. 2024) on the default, and offers an open fallback.

- **Gurobi — the recommended default.** It is the solver KRW used, and
  on this problem it is by far the fastest backend; a real 97-firm solve
  lands in a few minutes. The one cost is that Gurobi is commercial and
  needs a license (free academic licenses are available).
- **HiGHS — the open, license-free fallback** (`backend = "highs"`, via
  the `highs` package). It is several times slower on a real 97-firm
  solve, but it installs with no license and reaches the same grading.
  (Two more license-free routes, `"glpk"` and `"symphony"`, go through
  ROI plugins; HiGHS is the one to reach for first.)

You select the backend two ways. Per run, set it on the control object —
this just *builds a configuration* and does not solve, so it runs
instantly:

``` r

ctrl <- gp_control(backend = "highs")   # license-free; builds instantly (no solve)
ctrl$backend
#> [1] "highs"
```

Or set it once as the session-wide default.
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)
with no backend argument reads
`getOption("gradepath.backend", "gurobi")`, so this option becomes the
default every later
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md),
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md),
and
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
will use:

``` r

options(gradepath.backend = "gurobi")   # global default for the session
```

A control carries more than the backend — the `lambda_grid`, the
`precision_rule`, the credible-interval level, and solver knobs like
`mip_gap` and `time_limit` — but for solver choice the `backend` field
is the one that matters. You then hand the control to
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
or
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md)
via `control =`.

Two honest points, stated plainly. First, the backend is checked by the
solver *at call time*, not when you build the control —
`gp_control(backend = "highs")` succeeds whether or not the `highs`
package is installed; you find out when you actually solve. Second, and
more reassuring: **the package installs, loads, and runs its full
structural test suite with no Gurobi license at all.** Only an actual
grade solve needs a backend — Gurobi or HiGHS. So you can explore the
bundled fits, preview runs, read tables, and draw figures license-free,
and you reach for a solver only when you grade fresh data. On attach,
gradepath prints a one-line note reporting the backend it resolved, so
you always know which solver a fresh solve *would* use. A HiGHS solve
reaches the same grading Gurobi does; it is simply slower, and a
non-optimal solve is never relabelled optimal — the solver status is
reported honestly.

## 3. Sizing a run with `gp_preview()`

A real solve costs minutes, so before you commit you want two things:
confirmation that your input is well-formed, and an estimate of how much
work the run is.
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md)
gives you both *without solving anything*.

It takes the same input and `demographic` you would hand to
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
validates the input grammar, and reports the **plan** — the unit ids,
the inferred workflow, the backend, the frontier penalty grid, and the
number of grade-IP solves the run would perform. It runs **no**
beta-GMM, **no** deconvolution, **no** posterior, **no** pairwise, and
**no** solver stage. So it is instant, needs no license, and we run it
live:

``` r

pv <- gp_preview(inp, demographic = "race")
pv                   # result-first print of the planned run
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
#> [1] "2026-06-06 08:39:31 CDT"
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
```

The print is result-first: it shows the planned run at a glance. The
object is a `gp_preview` you can read field by field — the operational
ones are:

``` r

pv$n_units           # units the run would grade
#> [1] 24
pv$estimated_solves  # number of grade-IP solves (= length of the lambda grid)
#> [1] 12
pv$backend           # the backend the run would use
#> [1] "gurobi"
pv$workflow          # the inferred workflow ("one_level_independence" here)
#> [1] "one_level_independence"
```

The full slot set is all thirteen of `ids` (the canonical unit-id order
parsed from the input), `n_units` (the unit count, `length(ids)`),
`demographic`, `groups` (the parsed grouping, or `NULL`), `workflow`,
`lambda_grid` (the penalty grid the run would sweep),
`estimated_solves`, `backend`, the `control` it was built from, the
`status` — a `gp_status` object that reads `OK` for a clean one-level
plan; and three bookkeeping slots: `schema_version` (the object’s schema
tag, here `"v2"`), `provenance` (a short record of how the preview was
produced — `producer = "gp_preview"`, `no_solve = TRUE`), and `warnings`
(any diagnostics raised while building the plan — empty for a clean
preview).

`estimated_solves` is the headline number: each entry in the lambda grid
is one integer-program solve, so on a real fit this is how many
minutes-long solves the run will perform. The default grid is dense
(`seq(0, 1, by = 0.01)`, 101 points), so a production run narrows it —
for the parity result the grid is `c(0.25, 1)`, two solves. This is
exactly how you size a 97-firm run before committing: preview it, read
`n_units` and `estimated_solves`, multiply by the per-solve cost your
backend gives you (Gurobi minutes, HiGHS more), and decide whether to
run it now or overnight. Preview the tiny input first to learn the
shape; preview the real input to learn the bill.

## 4. Calibration with `gp_calibrate()`

[`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md)
is the package’s seeded Monte-Carlo calibration harness (GP-DEC-18-A).
It asks a different question from everything above: *does the grading
method’s own frequentist behaviour hold up under simulation?* It treats
a fitted prior as ground truth, simulates `n_sim` synthetic datasets
from it, re-runs the full one-level pipeline
([`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md))
on each, and checks three guarantees *qualitatively*:

1.  the empirical discordance rate (DR) tracks `dr_target`;
2.  the posterior credible intervals cover the truth near `ci_level`;
3.  the empirical grading’s regret against an oracle that *knows* the
    true prior is small and **non-negative by construction** (both
    gradings are scored on the same truth, so regret cannot go below
    zero).

The signature is

``` r

gp_calibrate(fit, n_sim = 200L, seed = 15238L,
             dr_target = 0.05, ci_level = 0.9, min_ok = 1L)
```

and it returns a compact `gp_calibration` carrying `n_sim`, `seed`,
`characteristic` (`"race"` or `"gender"`), `n_ok` and `n_failed`
(successful and skipped draw counts, which sum to `n_sim`), `dr_mean`
and `dr_target`, `coverage` and `ci_level`, `regret_mean`, the two
boolean verdicts `dr_ok` (DR within `0.01` of target) and `coverage_ok`
(coverage within `0.02` of `ci_level`), plus provenance. The harness is
robust to refit non-convergence: a draw whose refit fails or degenerates
is skipped and counted in `n_failed`, the means are taken over the
successful draws only, and if fewer than `min_ok` succeed it errors
rather than averaging noise.

**The call — shown, not run.**
[`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md)
re-solves the grade integer program *once per draw*, and the package
states it **never runs on package load or on vignette render**. So we
show it only, with `eval = FALSE`; keep `n_sim` small because each draw
costs a solve:

``` r

cal <- gp_calibrate(fit, n_sim = 5, seed = 15238)   # each draw re-solves the IP
cal                                                  # compact, result-first one-liner
summary(cal)                                         # one-row typed data frame
```

**The budget heads-up.** Because a run can be long,
[`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md)
prints a one-line, suppressible
[`message()`](https://rdrr.io/r/base/message.html) *before* the draw
loop whenever `n_sim * n_units` exceeds 600 — informational, never
changing the result. On this 97-unit fit, `n_sim = 20` gives
`20 * 97 = 1940 > 600`, so you would see:

> `gp_calibrate(): about to run 20 refit+grade solves on 97 units. The grade step is an integer program (often minutes per solve on a large fit), so this can take a long time. For a quick check use a smaller n_sim or a smaller, well-identified fit; gp_preview() gauges one solve's cost.`

Silence it with
[`suppressMessages()`](https://rdrr.io/r/base/message.html) once you
have sized the run.

**The honest scope.** State this plainly:

- [`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md)
  is a **qualitative, in-package diagnostic — not a reproduction of KRW
  Table E1.** There is no parseable E1 artifact and the Matlab
  Monte-Carlo stream cannot be reproduced, so an exact E1 match is out
  of reach. The `15238` default seed is a **cosmetic homage** to the
  reference seed, not Matlab parity (the provenance records
  `matlab_rng_parity = FALSE`); it only fixes in-package reproducibility
  — the same seed yields identical output, and your global RNG stream is
  restored afterward.
- It is **not a headline result and not part of any default scorecard.**
- On **tiny or weakly identified fits the verdicts are frequently
  `FALSE`** — the small-N beta-GMM is noisy and its r-scale can be
  ill-identified — so a meaningful calibration needs a **well-identified
  fit** (KRW’s real multi-firm input). The harness reports honestly
  either way; it never tunes the data to force a pass.

For the full calibration scope, including how this sits relative to the
deferred E1 Monte-Carlo, see
[M5](https://joonho112.github.io/gradepath/articles/m5-two-level-and-calibration.md).

## 5. Troubleshooting — the sharpened errors

gradepath’s errors are written to name *exactly* what went wrong and how
to fix it. Each demo below errors **immediately**, with no solve, so we
run them live through [`try()`](https://rdrr.io/r/base/try.html) (which
lets the chunk continue past the error).

**Stage-wise assembly with a missing argument.** If you assemble a
report card from individual stage objects and forget one, the error
names every missing argument and hands you a copy-paste recipe built
from a fit’s slots:

``` r

try(gp_report_card(grade_path = fit$grade_path))
#> Error : `gp_report_card()` stage-wise assembly is missing required argument(s): posterior, selected_grade. Supply `posterior`, `selected_grade`, and `grade_path` together -- e.g. from a fit's slots: gp_report_card(estimates = fit$estimates, posterior = get_posterior(fit), selected_grade = fit$selected_grade, grade_path = fit$grade_path).
```

The message names the missing `posterior` and `selected_grade` and gives
the full
`gp_report_card(estimates = fit$estimates, posterior = get_posterior(fit), selected_grade = fit$selected_grade, grade_path = fit$grade_path)`
form — no guessing which slot you dropped.

**A bad control value.**
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)
validates each field as you build it, so a typo errors clearly at
construction time rather than surfacing later inside a solve:

``` r

try(gp_control(precision_rule = "nope"))
#> Error : `precision_rule` must be one of: none, krw_gmm.
```

**A grouped input is previewed, not run.** Grouped, two-level execution
is the M2 surface, not the one-level M1 monolith. If you hand
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md)
a non-`NULL` `groups` vector, it does not error outright — it returns a
`gp_preview` whose `status` is `GROUPS_ERROR`, with an explanatory
warning, so you can see the planned (but unrun) grouped workflow:

``` r

pv_grouped <- suppressWarnings(
  gp_preview(inp, demographic = "race",
             groups = rep(c("A", "B"), length.out = pv$n_units))
)
pv_grouped$status
#> $status
#> [1] "GROUPS_ERROR"
#> 
#> $message
#> [1] "Grouped/two-level execution is not implemented in the M1 preview."
#> 
#> $detail
#> list()
#> 
#> attr(,"class")
#> [1] "gp_status" "list"
pv_grouped$workflow
#> [1] "grouped_pending"
```

The `GROUPS_ERROR` status and the `"grouped_pending"` workflow are the
signal: the one-level pipeline will not run a grouped input, and the
two-level surface lives in the m-track — see
[M5](https://joonho112.github.io/gradepath/articles/m5-two-level-and-calibration.md).
Treat the preview as a routing check, not a failure.

## 6. Where to next

You now have the operational toolkit: pick a backend, size a run with
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md),
calibrate within scope with
[`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md),
and read the sharpened errors. The applied track ends here; the method
track derives what these operations sit on top of.

- **The applied track** —
  [A1](https://joonho112.github.io/gradepath/articles/a1-getting-started.md)
  is the gentle first contact,
  [A2](https://joonho112.github.io/gradepath/articles/a2-the-grading-workflow.md)
  walks the seven stages behind a fit,
  [A3](https://joonho112.github.io/gradepath/articles/a3-reading-and-exporting-results.md)
  covers the accessors and export, and
  [A4](https://joonho112.github.io/gradepath/articles/a4-figures.md)
  builds the four signature figures.
- **Calibration in full** —
  [M5](https://joonho112.github.io/gradepath/articles/m5-two-level-and-calibration.md)
  gives the complete calibration scope (and the two-level surface that
  `GROUPS_ERROR` points to), including how
  [`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md)
  relates to the deferred Table E1 Monte-Carlo.
- **The grading the solver produces** —
  [M4](https://joonho112.github.io/gradepath/articles/m4-grading-frontier-and-report-cards.md)
  derives the grade integer program, the frontier, and the report card;
  [M1](https://joonho112.github.io/gradepath/articles/m1-foundations-and-notation.md)
  fixes the notation.

The honest bottom line, restated: Gurobi is the recommended default and
HiGHS the open fallback; a non-optimal solve is never described as
optimal;
[`gp_calibrate()`](https://joonho112.github.io/gradepath/reference/gp_calibrate.md)
is a diagnostic, not an E1 reproduction; and the one-level race and
gender core is accepted within its declared scope.

### Provenance

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US/en_US/en_US/C/en_US/en_US
#> 
#> time zone: America/Chicago
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.3   gradepath_0.5.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5       ebrecipe_0.5.0     gtable_0.3.6       jsonlite_2.0.0    
#>  [5] dplyr_1.2.1        compiler_4.6.0     tidyselect_1.2.1   slam_0.1-55       
#>  [9] dichromat_2.0-0.1  jquerylib_0.1.4    splines_4.6.0      systemfonts_1.3.2 
#> [13] scales_1.4.0       textshaping_1.0.5  yaml_2.3.12        fastmap_1.2.0     
#> [17] lattice_0.22-9     R6_2.6.1           generics_0.1.4     knitr_1.51        
#> [21] htmlwidgets_1.6.4  gurobi_13.0-1      tibble_3.3.1       desc_1.4.3        
#> [25] bslib_0.11.0       pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.2.0       
#> [29] cachem_1.1.0       xfun_0.57          fs_2.1.0           sass_0.4.10       
#> [33] S7_0.2.2           otel_0.2.0         cli_3.6.6          withr_3.0.2       
#> [37] pkgdown_2.2.0      magrittr_2.0.5     digest_0.6.39      grid_4.6.0        
#> [41] lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5     glue_1.8.1        
#> [45] farver_2.1.2       ragg_1.5.2         rmarkdown_2.31     tools_4.6.0       
#> [49] pkgconfig_2.0.3    htmltools_0.5.9
```

## References

Kline, Patrick, Evan K. Rose, and Christopher R. Walters. 2024. “A
Discrimination Report Card.” *American Economic Review* 114 (8):
2472–525. <https://doi.org/10.1257/aer.20230700>.
