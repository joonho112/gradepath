# Wrap two-level Pi matrices as solver-ready pairwise objects

`gp_twolevel_pairwise()` is the M2 pairwise dispatch for two-level
industry grading. It applies the same cleanup policy used by the
one-level core (antisymmetry, diagonal `0.5`, and the package zero
floor) to raw two-level outranking matrices, returning a bundle that
carries both the firm-level `industry_rfe` input (`Pi_theta`) and the
between-industry `btwn` input (`Pi_bar`) as validated `gp_pairwise`
objects ready for
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md).
It is solve-free.

## Usage

``` r
gp_twolevel_pairwise(
  Pi_theta,
  Pi_bar,
  ids = NULL,
  industry_levels = NULL,
  control = NULL
)
```

## Arguments

- Pi_theta:

  Numeric square matrix; the firm-level two-level
  `Pr(theta_i > theta_j)` outranking probabilities. Must have at least
  two units and finite entries in `[0, 1]`; the diagonal is forced to
  zero on entry.

- Pi_bar:

  Numeric square matrix; the between-industry `Pi_bar` / `Pi_sbar_psi`
  outranking probabilities at the industry representatives. Same shape
  and value constraints as `Pi_theta`.

- ids:

  Optional character vector of firm ids labelling `Pi_theta`. When
  `NULL` (default), the matrix dimnames are used, or `unit_1`, `unit_2`,
  ... are generated. Must be unique and match `nrow(Pi_theta)`.

- industry_levels:

  Optional character vector of industry ids labelling `Pi_bar`. When
  `NULL` (default), the matrix dimnames are used, or `industry_1`,
  `industry_2`, ... are generated. Must be unique and match
  `nrow(Pi_bar)`.

- control:

  Optional
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  object threaded onto both `gp_pairwise` wrappers. When `NULL`
  (default) a default
  [`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)
  is used.

## Value

A validated `gp_twolevel_pairwise_bundle` object (a list of class
`c("gp_twolevel_pairwise_bundle", "list")`) with the public slots:

- `ids`:

  Character vector; the firm ids for `Pi_theta`.

- `industry_levels`:

  Character vector; the industry ids for `Pi_bar`.

- `raw`:

  Named list `Pi_theta` / `Pi_bar`: the inputs with only the diagonal
  zeroed (before the cleanup policy).

- `Pi_theta`, `Pi_bar`:

  The cleaned firm- and industry-level matrices (antisymmetry, diagonal
  `0.5`, zero floor applied).

- `pairwise_theta`, `pairwise_bar`:

  The matching `gp_pairwise` objects fed to
  [`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md).

- `cleanup`:

  Named list recording the applied policy (`antisymmetry`, `diagonal`,
  `zero_floor`).

- `control`:

  The validated
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  used.

- `provenance`, `schema_version`, `warnings`:

  Internal audit slots.

## Details

This is a structural M2 adapter, not an acceptance decision: it does not
run the fixture-promotion gate, set `PROMOTED` / `banded` status,
compare industry DR, tau, or R2 against paper targets, or update any M1
cache. M2 acceptance is the job of
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md)
/
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md).

## See also

[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md),
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)

Other gradepath-twolevel:
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
[`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md),
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[`gp_twolevel_grade()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_grade.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)

## Examples

``` r
# A small synthetic strict-preference pair (the package's own test pattern):
# in strict(n, p) every off-diagonal entry is 0.95 above the diagonal and 0.05
# below, so the units are cleanly ordered. This step is solve-free.
strict <- function(n, p) {
  id <- paste0(p, seq_len(n)); m <- matrix(0, n, n, dimnames = list(id, id))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) m[i, j] <- if (i < j) 0.95 else 0.05
  m
}
pw <- gp_twolevel_pairwise(strict(4, "f"), strict(3, "i"),
                           control = gp_control(backend = "highs"))
pw
#> $ids
#> [1] "f1" "f2" "f3" "f4"
#> 
#> $industry_levels
#> [1] "i1" "i2" "i3"
#> 
#> $raw
#> $raw$Pi_theta
#>      f1   f2   f3   f4
#> f1 0.00 0.95 0.95 0.95
#> f2 0.05 0.00 0.95 0.95
#> f3 0.05 0.05 0.00 0.95
#> f4 0.05 0.05 0.05 0.00
#> 
#> $raw$Pi_bar
#>      i1   i2   i3
#> i1 0.00 0.95 0.95
#> i2 0.05 0.00 0.95
#> i3 0.05 0.05 0.00
#> 
#> 
#> $Pi_theta
#>      f1   f2   f3   f4
#> f1 0.50 0.95 0.95 0.95
#> f2 0.05 0.50 0.95 0.95
#> f3 0.05 0.05 0.50 0.95
#> f4 0.05 0.05 0.05 0.50
#> 
#> $Pi_bar
#>      i1   i2   i3
#> i1 0.50 0.95 0.95
#> i2 0.05 0.50 0.95
#> i3 0.05 0.05 0.50
#> 
#> $pairwise_theta
#> +---------------------------------------------------+
#> | gp_pairwise  .  4 x 4  .  12 ordered pairs        |
#> +---------------------------------------------------+
#> | pi range (off-diag): [0.050, 0.950]   diag = 0.50 |
#> | power = 0   rule = groupfx1_archive_matrix        |
#> +---------------------------------------------------+
#> i  matrix: as.matrix(pw)   ids: pw$ids
#> 
#> $pairwise_bar
#> +---------------------------------------------------+
#> | gp_pairwise  .  3 x 3  .  6 ordered pairs         |
#> +---------------------------------------------------+
#> | pi range (off-diag): [0.050, 0.950]   diag = 0.50 |
#> | power = 0   rule = groupfx1_archive_matrix        |
#> +---------------------------------------------------+
#> i  matrix: as.matrix(pw)   ids: pw$ids
#> 
#> $cleanup
#> $cleanup$antisymmetry
#> [1] TRUE
#> 
#> $cleanup$diagonal
#> [1] 0.5
#> 
#> $cleanup$zero_floor
#> [1] 1e-07
#> 
#> 
#> $control
#> $lambda_grid
#>   [1] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14
#>  [16] 0.15 0.16 0.17 0.18 0.19 0.20 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29
#>  [31] 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44
#>  [46] 0.45 0.46 0.47 0.48 0.49 0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59
#>  [61] 0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74
#>  [76] 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
#>  [91] 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.00
#> 
#> $backend
#> [1] "highs"
#> 
#> $precision_rule
#> [1] "none"
#> 
#> $interval_level
#> [1] 0.9
#> 
#> $solver_options
#> list()
#> 
#> $seed
#> NULL
#> 
#> $schema_version
#> [1] "v2"
#> 
#> $provenance
#> $provenance$producer
#> [1] "gp_control"
#> 
#> $provenance$built_at
#> [1] "2026-06-06 08:37:52 CDT"
#> 
#> $provenance$r_version
#> [1] "R version 4.6.0 (2026-04-24)"
#> 
#> $provenance$package_version
#> [1] "0.5.0"
#> 
#> 
#> attr(,"class")
#> [1] "gp_control" "list"      
#> 
#> $schema_version
#> [1] "v2"
#> 
#> $provenance
#> $provenance$producer
#> [1] "gp_twolevel_pairwise"
#> 
#> $provenance$n_units
#> [1] 4
#> 
#> $provenance$n_industries
#> [1] 3
#> 
#> $provenance$cleanup_order
#> [1] "antisymmetry_diagonal_zero_floor"
#> 
#> 
#> $warnings
#> character(0)
#> 
#> attr(,"class")
#> [1] "gp_twolevel_pairwise_bundle" "list"                       
```
