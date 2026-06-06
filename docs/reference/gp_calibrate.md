# Seeded Monte-Carlo calibration of the one-level KRW grading method

`gp_calibrate()` is a standalone, seeded, qualitative Monte-Carlo
calibration harness. It treats a fitted prior as ground truth, simulates
`n_sim` synthetic datasets from it, re-runs the full one-level pipeline
([`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md))
on each, and checks that the grading method's own frequentist guarantees
hold qualitatively:

- the empirical discordance rate (DR) tracks `dr_target`;

- posterior credible intervals cover the truth near `ci_level`;

- the empirical grading's regret against an oracle that KNOWS the true
  prior is small and non-negative.

It is NOT a replication of KRW Table E1 (there is no parseable E1
artifact and the Matlab Monte-Carlo stream cannot be reproduced; see the
calibration scope in the `m5-two-level-and-calibration` vignette), is
NOT seeded to match Matlab (the `15238` default is a cosmetic homage;
the seed only fixes in-package reproducibility), and is NOT a headline
result or part of any default scorecard. It never runs on package load
or on vignette render. The result is a compact `gp_calibration` object:
scalar summaries plus the two boolean verdicts `dr_ok` / `coverage_ok`.

## Usage

``` r
gp_calibrate(
  fit,
  n_sim = 200L,
  seed = 15238L,
  dr_target = 0.05,
  ci_level = 0.9,
  min_ok = 1L
)
```

## Arguments

- fit:

  A `gp_fit` from
  [`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
  /
  [`gradepath()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
  whose fitted prior, precision exponent `beta` (and additive `mu` for
  gender), and observed standard errors `s_i` define the data-generating
  truth.

- n_sim:

  Integer; number of Monte-Carlo synthetic datasets to draw and re-fit.
  Default `200L`. Use a small value (e.g. `5L`) for a fast check; each
  draw triggers one grade-IP solve.

- seed:

  Integer; RNG seed for in-package reproducibility. Default `15238L` (a
  cosmetic homage to the reference seed – it does NOT reproduce the
  Matlab stream). The same `seed` yields identical output; the caller's
  RNG stream is restored afterward.

- dr_target:

  Numeric in `[0, 1]`; the target discordance rate the mean empirical DR
  is compared against for the `dr_ok` verdict (passes when within
  `0.01`). Default `0.05`.

- ci_level:

  Numeric in the open interval `(0, 1)`; the nominal credible-interval
  level, also used as each refit's reporting-interval level. The
  `coverage_ok` verdict passes when mean coverage is within `0.02` of
  this. Default `0.90`.

- min_ok:

  Integer `>= 1`; the minimum number of successful (non-skipped) draws
  required before the harness reports. If fewer draws succeed it errors
  rather than averaging noise. Default `1L`.

## Value

A validated `gp_calibration` object (a list of class
`c("gp_calibration", "list")`) summarizing the run:

- `n_sim`, `seed`:

  Integer; the requested number of draws and the RNG seed used.

- `characteristic`:

  `"race"` or `"gender"`; the truth's demographic.

- `n_ok`, `n_failed`:

  Integer; successful and skipped draw counts
  (`n_ok + n_failed == n_sim`).

- `dr_mean`, `dr_target`:

  Numeric in `[0, 1]`; the mean empirical discordance rate over
  successful draws and the target it is compared with.

- `coverage`, `ci_level`:

  Numeric; mean theta-scale interval coverage over successful draws and
  the nominal level it is compared with.

- `regret_mean`:

  Numeric; mean regret of the empirical grading versus the true-prior
  oracle (non-negative by construction).

- `dr_ok`, `coverage_ok`:

  Logical; the two qualitative verdicts (DR within `0.01` of target;
  coverage within `0.02` of `ci_level`).

- `provenance`:

  Named list: producer, selected `lambda`, `beta_truth`, `n_units`,
  oracle and coverage-scale tags, and `matlab_rng_parity = FALSE`.

## Details

The harness is robust to refit non-convergence. Each draw's refit is
wrapped in `tryCatch`; a failed refit (e.g. a singular step-1 GMM moment
covariance) or a degenerate refit (an implausibly large `|beta|` that
corrupts the r-scale) is skipped and counted in `n_failed`. The reported
means are over the successful draws only. If fewer than `min_ok` draws
succeed, `gp_calibrate()` errors rather than averaging noise.

Two scale subtleties are handled internally: the ORACLE grading is built
from the known true prior WITHOUT re-deconvolving (its pairwise comes
from `gp_posterior_weights(prior = truth, ...)`), and both gradings are
scored on that same truth pairwise, so regret is non-negative by
construction. COVERAGE is assessed on the beta-invariant theta scale
(the refit's theta-scale reporting interval against the true
`theta_i = s_i^beta * v_i`), which stays coherent even when a refit's
estimated beta differs from the truth's.

Reproducibility: the whole draw loop runs under a single
state-preserving seeded context, so the same `seed` yields identical
results and the caller's global RNG stream is left untouched.

Cost: each draw re-solves a grade integer program once, so a run costs
about `n_sim` refit-and-grade solves. The function prints a one-line
budget heads-up (an ordinary, suppressible
[`message()`](https://rdrr.io/r/base/message.html)) before the loop
whenever `n_sim * <number of units>` exceeds 600; the message is
informational and never changes the result. For a quick check use a
small `n_sim`.

On tiny or weakly identified fits the verdicts are frequently `FALSE` –
the small-N KRW beta-GMM is noisy and its r-scale can be ill-identified
– so meaningful calibration needs a well-identified fit (KRW's real
multi-firm input). The harness reports honestly either way; it never
tunes the data to force a pass.

## Note

Not a headline result: the published Table E1 Monte-Carlo is deferred –
gp_calibrate is a qualitative in-package diagnostic, not an E1
reproduction. See the `m5-two-level-and-calibration` vignette.

## See also

[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)

## Examples

``` r
# A pre-solved tiny one-level fit bundled with the package.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))

# A live run re-solves the grade IP once per draw, so it needs a backend and is
# slow; keep n_sim small. (The harness also prints a budget heads-up before the
# loop when n_sim * n_units > 600.)
# \donttest{
cal <- gp_calibrate(fit, n_sim = 5, seed = 15238)
cal                       # compact, result-first one-liner
#> <gp_calibration> n_sim=5 seed=15238 race | DR mean=0.045 (target 0.050) [OK] | coverage=0.717 (90%) [--] | regret mean=0.0172
summary(cal)              # one-row typed data frame
#>   n_sim n_ok n_failed characteristic    dr_mean dr_target dr_ok  coverage
#> 1     5    5        0           race 0.04453321      0.05  TRUE 0.7166667
#>   ci_level coverage_ok regret_mean
#> 1      0.9       FALSE  0.01715955
# }
```
