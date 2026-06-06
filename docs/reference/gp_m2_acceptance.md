# Build the M2 acceptance scorecard

`gp_m2_acceptance()` assembles the M2 L01 / L02 / N10 ledger and the
overall formal M2 status as a `gp_m2_acceptance` object. With the
default `fixture_gates = NULL` it uses the recorded fixture snapshot
rather than rerunning the slow fixture gate; pass live fixture-gate
evidence to rebuild the scorecard. For a one-screen summary use
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md).

## Usage

``` r
gp_m2_acceptance(fixture_gates = NULL, registry = NULL)
```

## Arguments

- fixture_gates:

  Optional fixture-gate evidence. When `NULL` (default) the frozen
  race/gender fixture snapshot is used instead of rerunning the slow
  fixture gate. May also be a live `gp_twolevel_fixture_gate` object, a
  list of such objects, or an equivalent data frame, to recompute the
  promotion decision from fresh evidence.

- registry:

  Optional registry data frame supplying the L01 / L02 / N10 target
  rows. When `NULL` (default) the package
  [gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md)
  is used.

## Value

A validated `gp_m2_acceptance` object (a list of class
`c("gp_m2_acceptance", "list")`) with the public slots:

- `table`:

  Data frame; the full M2 ledger, one row per gate (`M2_ACCEPTANCE`,
  `L01`, `L02`, `N10`, `OVERRIDE`) with `status`, `reason`,
  `producer_status`, registry/effective class, paper/replicated values,
  and notes. The `M2_ACCEPTANCE` row carries the overall formal status
  (`PARTIAL_ACCEPTED` by default).

- `checks`:

  An alias of `table`.

- `effective_registry`:

  Registry data frame after panel-by-panel L02 promotion (see
  [`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md)).

- `fixture_gates`:

  The race/gender fixture-gate data frame used.

- `summary`:

  One-row data frame counting rows by status (`pass`, `promoted`,
  `approximate_ok`, `evidence_ok`, `fail`).

- `provenance`, `warnings`:

  Internal audit slots.

## Details

The default recorded scorecard is `PARTIAL_ACCEPTED`: L01 industry
grade-count rows pass exactly; N10 rows are synthetic support-guard
evidence rather than registered paper-value passes; gender L02
continuous rows are `PROMOTED` to `banded` by fixture parity; and race
L02 continuous rows remain `APPROXIMATE_OK` because the recorded race
`Pi_theta` fixture gap is `0.0121756787 > 0.01`.

`PROMOTED` / `banded` means only that a panel's fixture artifacts met
the intermediate parity band used for M2 promotion. It does NOT mean the
paper's industry DR, tau, or R2 targets were directly reproduced or
accepted; some continuous L02 rows are fixture-parity evidence, not
direct reproduction.

## See also

[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md),
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
[`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md)

Other gradepath-twolevel:
[`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md),
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
[`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)

## Examples

``` r
# Instant; no solve. The default recorded M2 scorecard.
acc <- gp_m2_acceptance()
acc$summary
#>   rows pass promoted approximate_ok evidence_ok fail
#> 1   14    2        3              3           5    0
```
