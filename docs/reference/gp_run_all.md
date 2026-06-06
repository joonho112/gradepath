# Run a registry-driven gradepath scorecard

`gp_run_all()` is the vector driver over
[`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md):
it applies the registry comparator to a set of replicated values (or to
a `gp_fit`, whose own targets it extracts) and returns a `gp_scorecard`
with the per-target table plus a grouped pass/fail summary. Passing
`replicated = NULL` returns the requested targets as `UNVERIFIED` (a
coverage skeleton with nothing solved).
[`gp_validate_targets()`](https://joonho112.github.io/gradepath/reference/gp_validate_targets.md)
is the thin verification-step alias.

## Usage

``` r
gp_run_all(
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

A list of class `c("gp_scorecard", "list")`.

- `table`:

  Per-target verdict data frame (the bound
  [`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md)
  rows with a leading `group` column); `checks` is a synonym slot with
  the same content.

- `summary`:

  Per-`group` and `TOTAL` tally of `pass` / `fail` / `unverified` /
  `n_a` / `no_tol` counts and a `pass_rate`.

- `pass_rate`:

  Numeric overall PASS fraction among scored (`PASS` + `FAIL`) rows, or
  `NA` when none are scored.

- `provenance`, `warnings`:

  Internal audit slots.

## Details

This is a verification scorecard against published KRW (2024) values; it
makes no ranking-superiority statement of any kind. Bundled artifacts
record the current M1 `NOT_ACCEPTED` gate: unresolved or non-`OK`
producer rows stay `UNVERIFIED` rather than counting as `PASS`.

## See also

[`gp_validate_targets()`](https://joonho112.github.io/gradepath/reference/gp_validate_targets.md),
[`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md),
[gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md)

Other gradepath-harness:
[`gp_check()`](https://joonho112.github.io/gradepath/reference/gp_check.md),
[`gp_validate_targets()`](https://joonho112.github.io/gradepath/reference/gp_validate_targets.md)

## Examples

``` r
# Instant, no solve: score explicit replicated values against the registry.
# Named vector of id = value; the firm-count target is published as 97.
sc <- gp_run_all(c(scale_n_firms_graded = 97))
sc$table[, c("id", "paper", "replicated", "status")]
#>                     id paper replicated status
#> 1 scale_n_firms_graded    97         97   PASS
sc$pass_rate
#> [1] 1

# NULL replicated returns the requested target as an UNVERIFIED skeleton.
gp_run_all(NULL, targets = "scale_n_firms_graded")$table$status
#> [1] "UNVERIFIED"

# \donttest{
# A gp_fit input derives and scores the fit's own registry targets.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
gp_run_all(fit)$summary
#>     group checks pass fail unverified n_a no_tol pass_rate
#> 1   TOTAL      6    0    6          0   0      0         0
#> 2 Targets      6    0    6          0   0      0         0
# }
```
