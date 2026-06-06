# gradepath replication registry (curated KRW paper-value targets)

A curated subset of the companion `paper_values.csv` (the published
Kline-Rose-Walters 2024 ground-truth registry), restricted to the
quantities `gradepath` owns and verifies. Each row carries the published
value, its unit, an absolute per-id tolerance, a tolerance class, and
the milestone at which it is asserted. It is read by the replication
harness
([`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md)
/
[`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md));
the harness embeds no tolerances of its own. Pairwise outranking
probabilities (K05–K07) deliberately have **no** registry row – they are
verified by a golden-master on the posterior weight matrix `W` instead.

## Usage

``` r
gp_registry
```

## Format

A data frame with 69 rows and 7 columns:

- id:

  character. Stable target key: the `paper_values.csv` key verbatim
  where one exists; otherwise a synthetic key (the single such row is
  `n9_strict4_ngrades`, the N9 strict-order grade-label guard).

- paper_value:

  character. Published value (stored character, CSV-faithful; coerce
  with `as.numeric`).

- unit:

  character. One of `count`, `percent`, `proportion`, `correlation`,
  `sd`, `other`.

- tolerance:

  character. Absolute per-id PASS tolerance (coerce to numeric); `0` for
  exact counts.

- class:

  character. `exact`, `banded`, or `approximate` (tolerance provenance,
  per the Chapter 17 ledger).

- milestone:

  character. `M1` (current one-level firm gate), `M1_NAMES` (deferred
  names pipeline), or `M2` (two-level / industry).

- quantity:

  character. Human-readable description of the target.

## Source

Curated from the KRW (2024) companion `paper_values.csv`; tolerance
classes from the package's tolerance ledger. Names rows are retained in
the registry for paper-value traceability but assigned to the deferred
`M1_NAMES` milestone until the names-specific NPMLE/ranking pipeline is
bundled and tested.

## Details

Values are sourced verbatim from the companion `paper_values.csv` (the
single source of truth); the `class`, `tolerance`, and `milestone`
columns are attached from the package's tolerance ledger. The NPMLE name
atoms use the verbatim CSV keys `names_npmle_mass1/2/3`. The registry is
intentionally an M1-and-M2 core subset; names rows are retained under
`M1_NAMES`, and Table F5/F6 industry cells are assigned to `M2` until
the two-level / industry report-card pipeline is live.

## Examples

``` r
data(gp_registry)
str(gp_registry)
#> 'data.frame':    69 obs. of  7 variables:
#>  $ id         : chr  "scale_n_firms_graded" "scale_n_names" "race_n_industries" "race_baseline_ngrades" ...
#>  $ paper_value: chr  "97" "76" "19" "3" ...
#>  $ unit       : chr  "count" "count" "count" "count" ...
#>  $ tolerance  : chr  "0" "0" "0" "0" ...
#>  $ class      : chr  "exact" "exact" "exact" "exact" ...
#>  $ milestone  : chr  "M1" "M1_NAMES" "M1" "M1" ...
#>  $ quantity   : chr  "design: firms graded" "design: names ranked" "design: industry groupings" "race baseline grade count" ...
#>  - attr(*, "source")=List of 5
#>   ..$ name        : chr "KRW companion paper_values.csv"
#>   ..$ path        : chr "inst/extdata/registry/paper_values.csv"
#>   ..$ sha256      : chr "711a0afed4d6a30dda8385717138c46ef855c6c5307abe3a187a4d699db23573"
#>   ..$ rows        : int 297
#>   ..$ override_env: chr NA
# Count the rows per milestone gate.
table(gp_registry$milestone)
#> 
#>       M1 M1_NAMES       M2 
#>       31       17       21 
```
