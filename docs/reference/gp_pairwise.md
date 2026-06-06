# Pairwise outranking probability matrix

`gp_pairwise()` returns the `J x J` pairwise outranking matrix whose
entry `pi[i, j]` is the posterior probability `Pr(theta_i > theta_j)`.
Given a materialized `gp_fit` it is a read-only accessor for the stored
matrix; given the output of `gp_posterior_weights()` (or a bare `J x M`
weight matrix plus a reporting axis) it computes the matrix by exact
posterior-CDF integration and applies the frozen cleanup contract. The
result is the outranking structure the grade engine
([`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md))
scores.

## Usage

``` r
gp_pairwise(
  object = NULL,
  reporting_support = NULL,
  ids = NULL,
  control = NULL,
  assumption = "one_level_independence",
  ...,
  power = 0L
)
```

## Arguments

- object:

  A `gp_fit` (read-only accessor path), the list returned by
  `gp_posterior_weights()` carrying `$W` (`J x M`) and
  `$reporting_support` (`M x J`), or a bare `J x M` `W` matrix. A bare
  matrix is accepted only when `reporting_support` is supplied or a
  shared axis is otherwise intended. Default `NULL` requires the legacy
  `weights =` alias (see `...`).

- reporting_support:

  Optional `M x J` per-unit reporting axis; overrides
  `object$reporting_support` on the weight-list path. `NULL` (default)
  takes the axis from `object`. Ignored when `object` is a `gp_fit`.

- ids:

  Optional length-`J` identifiers. `NULL` (default) takes them from
  `object$ids`, the row names of `W`, or `seq_len(J)`. Ignored when
  `object` is a `gp_fit`.

- control:

  Optional
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  recorded on the result. `NULL` (default) uses
  [`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)
  on the weight-list path. Ignored when `object` is a `gp_fit`.

- assumption:

  Character scalar; the pairwise independence assumption recorded in
  `source`. Defaults to `"one_level_independence"` (the one-level
  analytic path). Ignored when `object` is a `gp_fit`.

- ...:

  Unused, except for the legacy named `weights =` alias for `object`;
  any other argument raises an error.

- power:

  Integer scalar; the matrix-power option. Only the KRW / public path
  `0L` (default) is implemented on the M1 surface; any other value
  errors.

## Value

A validated gp_pairwise object (a list of class
`c("gp_pairwise", "list")`) with the public slots:

- `ids`:

  Character vector of length `J`; the canonical unit-id order, and the
  row/column names of `matrix`.

- `matrix`:

  Numeric `J x J` matrix of outranking probabilities in `[0, 1]`;
  `matrix[i, j] = Pr(theta_i > theta_j)`, diagonal `0.5`.

- `power`:

  Integer; the matrix-power option (`0L`).

- `cleanup`:

  Named list (`antisymmetry`, `diagonal`, `zero_floor`); the applied
  cleanup contract – antisymmetry on, diagonal `0.5`, one-sided `1e-7`
  floor.

- `source`:

  Named list (`stage`, `rule`, `assumption`); how the matrix was
  produced and under which independence assumption.

- `control`:

  The
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  recorded on the result.

- `warnings`, `schema_version`, `provenance`:

  Internal audit slots.

## Details

On the compute path the entry is the per-unit CDF integral
`pi[i, j] = sum_m W[i, m] * F_j(theta_i[m]^-)`, where `F_j` is unit
`j`'s posterior CDF strictly below the argument (ties contribute `0`).
Cleanup is applied in the frozen order antisymmetry
(`pi[j, i] = 1 - pi[i, j]`) -\> diagonal set to `0.5` -\> zero-floor:
every off-diagonal value below `1e-7` is raised to `1e-7` (a one-sided
floor; exact `1` stays exact `1`), so final antisymmetry is exact up to
the floor tolerance. See KRW (2024, Section 4) and the pairwise method
vignette.

## See also

[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md),
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)

## Examples

``` r
# Read-only accessor: pull the stored pairwise matrix from the bundled fit.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
pw <- gp_pairwise(fit)
pw                              # a 24 x 24 gp_pairwise
#> +---------------------------------------------------+
#> | gp_pairwise  .  24 x 24  .  552 ordered pairs     |
#> +---------------------------------------------------+
#> | pi range (off-diag): [0.009, 0.991]   diag = 0.50 |
#> | power = 0   rule = outer_product                  |
#> +---------------------------------------------------+
#> i  matrix: as.matrix(pw)   ids: pw$ids
pw$matrix[1:4, 1:4]             # the outranking probabilities, diagonal 0.5
#>             firm01    firm02     firm03    firm04
#> firm01 0.500000000 0.9904685 0.34493294 0.8709196
#> firm02 0.009531549 0.5000000 0.06438113 0.1364485
#> firm03 0.655067061 0.9356189 0.50000000 0.7970593
#> firm04 0.129080374 0.8635515 0.20294068 0.5000000
```
