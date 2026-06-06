# Return a two-level industry report card

`gp_twolevel_report_card()` is the convenience selector for the report
cards assembled by
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md).
Passing an existing `gp_twolevel_grade` never re-solves; any other input
is forwarded to
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md)
first. Choose the firm-level (`industry_rfe`) or between-industry
(`btwn`) card with `model`.

## Usage

``` r
gp_twolevel_report_card(twolevel, model = c("industry_rfe", "btwn"), ...)
```

## Arguments

- twolevel:

  A `gp_twolevel_grade` object, or any input accepted by
  [`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md)
  (which is then graded first).

- model:

  `"industry_rfe"` (default) for the firm-level two-level report card or
  `"btwn"` for the between-industry report card.

- ...:

  Arguments forwarded to
  [`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md)
  when `twolevel` is not already a graded `gp_twolevel_grade` (e.g.
  `posterior`, `control`, `lambda`).

## Value

A validated `gp_report_card` for the chosen route: the per-unit
grade-label table with posterior summaries.

## Details

A report card exists only when the underlying two-level grade was built
WITH posterior information. This helper does not create new acceptance
evidence, set `PROMOTED` / `banded` status, or certify paper industry DR
/ tau / R2 reproduction; it simply returns the posterior-backed
`gp_report_card` already carried by the chosen route.

## Note

While the M2 industry surface is `PARTIAL_ACCEPTED`, no public
input-to-industry-card builder is exported, so there is no public way to
construct a posterior-backed two-level grade. Consequently a grade built
from bare Pi matrices
([`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md)
-\>
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md))
is structural only and carries no report card: calling this helper on it
raises an informative error directing you to
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md).
The success path is shown in the package's two-level vignettes; see also
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md).

## See also

[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
[`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md),
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)

Other gradepath-twolevel:
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
[`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md),
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
[`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md)

## Examples

``` r
# The M2 industry surface is PARTIAL_ACCEPTED: there is no public posterior-backed
# two-level grade to build here, so the runnable example inspects the M2 contract.
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
attr(gp_m2_status(), "m2_formal_status")   # "PARTIAL_ACCEPTED"
#> [1] "PARTIAL_ACCEPTED"

if (FALSE) { # \dontrun{
# By design this errors: a bare-Pi two-level grade is structural only and has no
# report card to return (no public builder is exported yet). See gp_m2_status().
strict <- function(n, p) {
  id <- paste0(p, seq_len(n)); m <- matrix(0, n, n, dimnames = list(id, id))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) m[i, j] <- if (i < j) 0.95 else 0.05
  m
}
pw <- gp_twolevel_pairwise(strict(4, "f"), strict(3, "i"),
                           control = gp_control(backend = "highs"))
gp_twolevel_report_card(pw, control = gp_control(backend = "highs"))
} # }
```
