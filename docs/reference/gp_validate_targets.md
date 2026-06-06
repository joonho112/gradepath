# Validate replicated target values against the gradepath registry

`gp_validate_targets()` is the explicit verification step that replaces
v1's replication-mode control flag: call it with the values you
replicated to confirm they match the published KRW (2024) registry
targets within tolerance. It is a thin convenience wrapper over
[`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md)
and returns the same `gp_scorecard`.

## Usage

``` r
gp_validate_targets(
  replicated = NULL,
  targets = NULL,
  registry = NULL,
  producer_status = NULL
)
```

## Arguments

- replicated:

  What to score. A `gp_fit` (its registry targets – firm count, grade
  count, discordance rate, tau, and for race the worst/best-grade counts
  – are derived and scored); a data frame with `id` and `replicated`
  columns (optionally `producer_status`, `group`); a named atomic vector
  or named list of `id = value`; or `NULL` (default) for an `UNVERIFIED`
  skeleton over `targets`.

- targets:

  Optional character vector of registry ids to score (and the row
  order). Defaults to `NULL`, which means all `M1` registry rows when
  `replicated` is `NULL`, otherwise the ids present in `replicated`.
  Requested ids absent from `replicated` are added as `UNVERIFIED` rows.

- registry:

  Registry data frame. Defaults to `NULL`, which uses the bundled
  package data
  [gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md).

- producer_status:

  Optional character producer-status override. Defaults to `NULL`: a
  `gp_fit` input then uses `fit$provenance$producer_status` when present
  and the fail-safe `"UNVERIFIED"` when absent; other replicated inputs
  (named vectors / lists / data frames) default to `"OK"`, the
  caller-asserts-OK contract for explicit replicated values.

## Value

A `gp_scorecard` (see
[`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md)
for the slot layout): `table`, `summary`, `pass_rate`, and internal
`provenance` / `warnings`.

## Details

A verification step against published values; it makes no
ranking-superiority statement of any kind.

## See also

[`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md),
[`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md),
[gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md)

Other gradepath-harness:
[`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md),
[`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md)

## Examples

``` r
# Instant, no solve: validate explicit replicated values against the registry.
sc <- gp_validate_targets(c(scale_n_firms_graded = 97))
sc$table[, c("id", "paper", "replicated", "status")]
#>                     id paper replicated status
#> 1 scale_n_firms_graded    97         97   PASS
```
