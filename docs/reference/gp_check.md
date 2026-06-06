# Compare a replicated value against the gradepath registry

`gp_check()` is the single registry comparator used by the replication
harness. It looks up one row of
[gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md)
by `id`, reads that row's published paper value and absolute tolerance,
compares your already-on-unit `replicated` value to it, and returns a
one-row verdict data frame. The harness embeds no tolerances of its own
– every band comes from the registry – and
[`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md)
is the vector driver that calls this once per target.

## Usage

``` r
gp_check(
  id,
  replicated = NA_real_,
  producer_status = "OK",
  registry = NULL,
  label = NULL
)
```

## Arguments

- id:

  Character scalar; a registry key present in `gp_registry$id`. An id
  not in the registry raises an error.

- replicated:

  Scalar replicated value, already converted to the registry unit (e.g.
  pass a percent for rows whose `unit` is `"percent"`). Defaults to
  `NA_real_`, which yields an `n/a` status (nothing to compare).

- producer_status:

  Character scalar; the runtime status of the producing stage,
  normalized internally. Defaults to `"OK"`. Any non-`OK` status short-
  circuits the numeric comparison and returns status `UNVERIFIED` (the
  producing stage was not acceptance-ready), regardless of `replicated`.

- registry:

  Registry data frame to look `id` up in. Defaults to `NULL`, which uses
  the bundled package data
  [gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md).

- label:

  Character scalar quantity label overriding `registry$quantity` in the
  returned row. Defaults to `NULL` (use the registry's own `quantity`
  text).

## Value

A one-row `data.frame` describing the comparison.

- `id`:

  The registry id checked.

- `quantity`:

  Human-readable target description (`label` if supplied, else the
  registry `quantity`).

- `paper`:

  Numeric published value from the registry.

- `replicated`:

  Numeric replicated value compared (`NA` when the producer status is
  non-`OK`).

- `delta`:

  `replicated - paper` (`NA` when not compared).

- `tol`:

  Absolute PASS tolerance from the registry.

- `unit`, `class`, `milestone`:

  The registry row's unit, tolerance class, and milestone, passed
  through.

- `status`:

  One of `PASS`, `FAIL`, `n/a` (missing value), `no-tol` (missing
  tolerance), or `UNVERIFIED` (non-`OK` producer status).

- `reason`:

  Short reason string when not `PASS` (e.g. `outside_tolerance`,
  `missing_value`, `missing_tolerance`); `NA` on `PASS`.

- `producer_status`:

  The normalized producer status used.

## Details

The comparison is the single Chapter 16/17 absolute-band rule: `PASS`
iff `abs(replicated - paper) <= tol` (with a tiny `1e-9` numerical
slack), `FAIL` otherwise. This is a verification check against published
values; it makes no ranking-superiority statement of any kind.

## See also

[`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md),
[`gp_validate_targets()`](https://joonho112.github.io/gradepath/reference/gp_validate_targets.md),
[gp_registry](https://joonho112.github.io/gradepath/reference/gp_registry.md)

Other gradepath-harness:
[`gp_run_all()`](https://joonho112.github.io/gradepath/reference/gp_run_all.md),
[`gp_validate_targets()`](https://joonho112.github.io/gradepath/reference/gp_validate_targets.md)

## Examples

``` r
# Instant: compare a replicated count against a known registry target. The
# firm-count target (id "scale_n_firms_graded") is published as 97, tolerance 0.
gp_check("scale_n_firms_graded", replicated = 97)   # status PASS
#>                     id             quantity paper replicated delta tol  unit
#> 1 scale_n_firms_graded design: firms graded    97         97     0   0 count
#>   class milestone status reason producer_status
#> 1 exact        M1   PASS   <NA>              OK
gp_check("scale_n_firms_graded", replicated = 96)   # status FAIL (outside tol 0)
#>                     id             quantity paper replicated delta tol  unit
#> 1 scale_n_firms_graded design: firms graded    97         96    -1   0 count
#>   class milestone status            reason producer_status
#> 1 exact        M1   FAIL outside_tolerance              OK
```
