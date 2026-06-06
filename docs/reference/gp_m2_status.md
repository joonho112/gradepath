# Summarize the M2 acceptance status

`gp_m2_status()` returns a compact public summary of the current M2
acceptance ledger – one row per status component – without touching the
M1 status contract in `R/status.R`. It is the one-screen read of
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md)
and the verb
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)
points users to while the M2 surface is `PARTIAL_ACCEPTED`.

## Usage

``` r
gp_m2_status(acceptance = NULL)
```

## Arguments

- acceptance:

  Optional `gp_m2_acceptance` object to summarize. When `NULL` (default)
  a fresh
  [`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md)
  (the recorded-snapshot scorecard) is built and summarized.

## Value

A `gp_m2_status` object: a data frame of class
`c("gp_m2_status", "data.frame")` with one row per component (`formal`,
`L01_industry_grade_counts`, `L02_race_continuous`,
`L02_gender_continuous`, `N10_support_guards`, `M1_status_boundary`) and
the columns `component`, `status`, `detail`, `evidence`,
`producer_status`, `source`. It carries the attributes `m2_acceptance`
(the backing `gp_m2_acceptance`), `m2_formal_status`
(`"PARTIAL_ACCEPTED"` by default), and `m1_status_recalculated` (always
`FALSE`).

## Details

The default summary is derived from
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md).
It reports the formal M2 status, L01 exact industry grade-count
evidence, L02 race/gender fixture-promotion status, and N10 synthetic
support-guard evidence. The default M2 state remains `PARTIAL_ACCEPTED`:
race L02 continuous rows remain approximate because the recorded
`Pi_theta` fixture gap is above `0.01`, while gender L02 rows are
promoted to `banded` by fixture parity.

This helper is M2-only. It does not use `new_gp_status()`, does not edit
`R/status.R`, and does not recalculate the M1 acceptance state; the
returned object's `m1_status_recalculated` attribute is always `FALSE`.

## See also

[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
[`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)

Other gradepath-twolevel:
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
[`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md),
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
[`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)

## Examples

``` r
# Instant; no solve. The per-component M2 ledger summary.
gp_m2_status()
#>                   component           status
#> 1                    formal PARTIAL_ACCEPTED
#> 2 L01_industry_grade_counts             PASS
#> 3       L02_race_continuous   APPROXIMATE_OK
#> 4     L02_gender_continuous         PROMOTED
#> 5        N10_support_guards      EVIDENCE_OK
#> 6        M1_status_boundary NOT_RECALCULATED
#>                                                              detail
#> 1                           ONE_OR_MORE_L02_ROWS_REMAIN_APPROXIMATE
#> 2                                        2/2 exact count rows pass.
#> 3   Race L02 remains approximate: Pi_theta gap 0.0121756787 > 0.01.
#> 4 Gender L02 promoted to banded: Pi_theta gap 0.0066335820 <= 0.01.
#> 5                     4/4 N10 rows carry synthetic evidence status.
#> 6       M1 acceptance status is not recalculated by gp_m2_status().
#>                                                                                                               evidence
#> 1                          Exact L01/N10 gates pass; L02 continuous rows promote panel-by-panel from fixture evidence.
#> 2                                          Race and gender industry grade-count targets are registered L01 exact rows.
#> 3                  Race industry DR/tau/R2 rows are classified by fixture parity, not direct paper-value reproduction.
#> 4                Gender industry DR/tau/R2 rows are classified by fixture parity, not direct paper-value reproduction.
#> 5                                    N10 support rows are synthetic guard evidence, not registered paper-value passes.
#> 6 This helper summarizes M2 only and does not edit or recompute R/status.R, the M1 cache, or M1 acceptance scorecards.
#>   producer_status                    source
#> 1              OK          gp_m2_acceptance
#> 2              OK          gp_m2_acceptance
#> 3  APPROXIMATE_OK  gp_twolevel_fixture_gate
#> 4              OK  gp_twolevel_fixture_gate
#> 5              OK test-n10-support-parity.R
#> 6  not_applicable              gp_m2_status

# The headline formal status is also an attribute.
attr(gp_m2_status(), "m2_formal_status")   # "PARTIAL_ACCEPTED"
#> [1] "PARTIAL_ACCEPTED"
```
