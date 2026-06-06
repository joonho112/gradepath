# Effective M2 registry after fixture-gate promotion

`gp_m2_promoted_registry()` returns a copy of the package registry whose
M2 L02 continuous-quantity `class` column is promoted from `approximate`
to `banded` only for panels whose fixture gate passes. The registry
keeps these rows `approximate` by default; this helper applies the
panel-by-panel promotion that
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md)
consumes.

## Usage

``` r
gp_m2_promoted_registry(fixture_gates = NULL, registry = NULL)
```

## Arguments

- fixture_gates:

  Optional fixture-gate evidence. When `NULL` (default) the frozen
  race/gender fixture snapshot is used instead of rerunning the fixture
  gate. May also be a live `gp_twolevel_fixture_gate` object, a list of
  such objects, or an equivalent data frame.

- registry:

  Optional registry data frame to promote. When `NULL` (default) the
  package
  [gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md)
  is used.

## Value

A registry data frame (a copy of the input) with the eligible L02
`class` entries promoted to `banded`, and two attributes attached:
`m2_fixture_gates` (the race/gender gate data frame used) and
`m2_l02_ids` (the L02 row ids considered for promotion).

## Details

With the default `fixture_gates = NULL` it uses the recorded race/gender
fixture snapshot: that snapshot leaves the race L02 rows `approximate`
(the race `Pi_theta` gap is `0.0121756787 > 0.01`) and promotes the
gender L02 rows to `banded` (the gender fixture gate passed). Supply
live fixture-gate evidence to recompute the promotion from fresh
fixtures.

`PROMOTED` / `banded` means only that a panel's fixture artifacts met
the intermediate parity band used for M2 promotion. It does NOT mean the
paper's industry DR, tau, or R2 targets were directly reproduced or
accepted.

## See also

[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md)

Other gradepath-twolevel:
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
[`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)

## Examples

``` r
# Instant; no solve. The effective registry after fixture promotion.
reg <- gp_m2_promoted_registry()
attr(reg, "m2_l02_ids")
#> [1] "race_industry_dr"    "race_industry_tau"   "race_industry_r2"   
#> [4] "gender_industry_dr"  "gender_industry_tau" "gender_industry_r2" 
```
