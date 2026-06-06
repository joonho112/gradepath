# Grade two-level industry matrices

`gp_twolevel_grade()` routes the M2 two-level pairwise outputs into the
package grade solver without touching the one-level monolith: the
`industry_rfe` route grades the firm-level `Pi_theta`, and the `btwn`
route grades the industry-level `Pi_bar`. Feed it a bundle from
[`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md)
(or the lower-level two-level carriers below). When a two-level
posterior is supplied it also assembles a firm `gp_report_card` for
`industry_rfe` and an industry-level `gp_report_card` for `btwn`,
retrievable via
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md).

## Usage

``` r
gp_twolevel_grade(
  twolevel = NULL,
  pairwise_theta = NULL,
  pairwise_bar = NULL,
  posterior = NULL,
  lambda = 0.25,
  lambda_grid = NULL,
  control = NULL,
  acceptance_mode = FALSE,
  build_report_cards = TRUE,
  ...
)
```

## Arguments

- twolevel:

  A `gp_twolevel_quadrature`, `gp_twolevel_pi`,
  `gp_twolevel_pairwise_bundle`, or a plain list carrying
  `pairwise_theta` and `pairwise_bar`. When `NULL` (default) the
  explicit `pairwise_theta` / `pairwise_bar` arguments must be supplied
  instead. Passing an already-graded `gp_twolevel_grade` returns it
  unchanged (no re-solve).

- pairwise_theta, pairwise_bar:

  Optional explicit `gp_pairwise` inputs for the firm-level and
  between-industry routes. When `NULL` (default) they are taken from
  `twolevel`; at least one source for each route must resolve.

- posterior:

  Optional two-level `gp_posterior` (e.g. from a two-level posterior
  builder). When `NULL` (default) the routes are graded structurally
  with no report cards. When supplied (and `build_report_cards = TRUE`)
  it backs the assembled report cards.

- lambda:

  Numeric in `[0, 1]`; the selected reporting penalty at which each
  route's grade is read off the path. Default `0.25`.

- lambda_grid:

  Optional numeric vector; the penalty grid solved for each route. When
  `NULL` (default) it is the package parity anchors plus `lambda`
  (`0.25`, `1`, and `lambda`), sorted and de-duplicated.

- control:

  Optional
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  threaded to the solves. When `NULL` (default) the resolved
  `pairwise_theta`'s own control is used.

- acceptance_mode:

  Logical; the solver solution-quality / fallback policy forwarded to
  the internal
  [`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
  (Gurobi backend only). `FALSE` (default) reports a gap/time-limit
  solve honestly with that status; `TRUE` adds optimization attempts and
  never relabels a non-optimal solve as optimal.

- build_report_cards:

  Logical; when `TRUE` (default) and a `posterior` is available,
  assemble the per-route `gp_report_card`s. When `FALSE` the routes are
  graded structurally only.

- ...:

  Reserved; an error is raised if any argument is passed.

## Value

A validated `gp_twolevel_grade` object (a list of class
`c("gp_twolevel_grade", "list")`) with the public slots:

- `ids`, `industry_levels`:

  Character vectors; the firm ids (`Pi_theta` route) and industry ids
  (`Pi_bar` route).

- `pairwise_theta`, `pairwise_bar`:

  The `gp_pairwise` inputs actually graded for the two routes.

- `industry_rfe`:

  Named list for the firm-level route: `grade_path`, `selected_grade`
  (at `lambda`), `grade_count`, `report_card` (or `NULL`), `posterior`,
  and `producer_status`.

- `btwn`:

  The same payload for the between-industry route.

- `posterior`:

  The supplied two-level `gp_posterior`, or `NULL`.

- `pi`, `pairwise`:

  The two-level carriers (`gp_twolevel_pi`,
  `gp_twolevel_pairwise_bundle`) when supplied, else `NULL`.

- `selected_lambda`, `lambda_grid`:

  The selected penalty and the solved grid.

- `method`:

  Character tag for how the inputs were resolved (e.g. `"pairwise"`,
  `"pi"`).

- `producer_status`:

  Combined honest solver status across the two routes (`"OK"` only when
  both routes are acceptance-ready).

- `control`:

  The validated
  [gp_control](https://joonho112.github.io/gradepath/reference/gp_control.md)
  used.

- `provenance`, `schema_version`, `warnings`:

  Internal audit slots.

## Details

The returned object is an M2 grading surface, not the M2 acceptance
scorecard and not a paper industry DR / tau / R2 reproduction. Use
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md)
to classify recorded or supplied fixture-gate evidence into `PROMOTED`,
`APPROXIMATE_OK`, and the overall M2 status; see
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md).

## See also

[`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md),
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)

Other gradepath-twolevel:
[`gp_m2_acceptance()`](https://joonho112.github.io/gradepath/reference/gp_m2_acceptance.md),
[`gp_m2_promoted_registry()`](https://joonho112.github.io/gradepath/reference/gp_m2_promoted_registry.md),
[`gp_m2_status()`](https://joonho112.github.io/gradepath/reference/gp_m2_status.md),
[`gp_twolevel_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_pairwise.md),
[`gp_twolevel_report_card()`](https://joonho112.github.io/gradepath/reference/gp_twolevel_report_card.md)

## Examples

``` r
# A small synthetic strict-preference pair (the package's own test pattern).
strict <- function(n, p) {
  id <- paste0(p, seq_len(n)); m <- matrix(0, n, n, dimnames = list(id, id))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) m[i, j] <- if (i < j) 0.95 else 0.05
  m
}
pw <- gp_twolevel_pairwise(strict(4, "f"), strict(3, "i"),
                           control = gp_control(backend = "highs"))

# The grade step solves a tiny integer program, so it needs a backend.
# \donttest{
tlg <- gp_twolevel_grade(pw, control = gp_control(backend = "highs"),
                         lambda_grid = c(0.25, 1))
tlg
#> $ids
#> [1] "f1" "f2" "f3" "f4"
#> 
#> $industry_levels
#> [1] "i1" "i2" "i3"
#> 
#> $pi
#> NULL
#> 
#> $pairwise
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
#> $posterior
#> NULL
#> 
#> $industry_rfe
#> $industry_rfe$model
#> [1] "industry_rfe"
#> 
#> $industry_rfe$pairwise
#> +---------------------------------------------------+
#> | gp_pairwise  .  4 x 4  .  12 ordered pairs        |
#> +---------------------------------------------------+
#> | pi range (off-diag): [0.050, 0.950]   diag = 0.50 |
#> | power = 0   rule = groupfx1_archive_matrix        |
#> +---------------------------------------------------+
#> i  matrix: as.matrix(pw)   ids: pw$ids
#> 
#> $industry_rfe$grade_path
#> +--------------------------------------------------------------+
#> | gp_grade_path  .  2 lambda values  .  selected lambda = 0.25 |
#> +--------------------------------------------------------------+
#> | grades across path in [4, 4]                                 |
#> | endpoint lambda = 1.00                                       |
#> +--------------------------------------------------------------+
#> i  path table: x$summary   per-lambda fits: x$fits
#> 
#> $industry_rfe$selected_grade
#> +---------------------------------------------+
#> | gp_grade_fit  .  4 grades  .  lambda = 0.25 |
#> +---------------------------------------------+
#> | units = 4   status = optimal                |
#> | grade 1: 1                                  |
#> | grade 2: 1                                  |
#> | grade 3: 1                                  |
#> | grade 4: 1                                  |
#> +---------------------------------------------+
#> i  assignment: x$assignment   summary(x) for backend/channel details
#> 
#> $industry_rfe$report_card
#> NULL
#> 
#> $industry_rfe$posterior
#> NULL
#> 
#> $industry_rfe$grade_count
#> [1] 4
#> 
#> $industry_rfe$producer_status
#> [1] "OK"
#> 
#> 
#> $btwn
#> $btwn$model
#> [1] "btwn"
#> 
#> $btwn$pairwise
#> +---------------------------------------------------+
#> | gp_pairwise  .  3 x 3  .  6 ordered pairs         |
#> +---------------------------------------------------+
#> | pi range (off-diag): [0.050, 0.950]   diag = 0.50 |
#> | power = 0   rule = groupfx1_archive_matrix        |
#> +---------------------------------------------------+
#> i  matrix: as.matrix(pw)   ids: pw$ids
#> 
#> $btwn$grade_path
#> +--------------------------------------------------------------+
#> | gp_grade_path  .  2 lambda values  .  selected lambda = 0.25 |
#> +--------------------------------------------------------------+
#> | grades across path in [3, 3]                                 |
#> | endpoint lambda = 1.00                                       |
#> +--------------------------------------------------------------+
#> i  path table: x$summary   per-lambda fits: x$fits
#> 
#> $btwn$selected_grade
#> +---------------------------------------------+
#> | gp_grade_fit  .  3 grades  .  lambda = 0.25 |
#> +---------------------------------------------+
#> | units = 3   status = optimal                |
#> | grade 1: 1                                  |
#> | grade 2: 1                                  |
#> | grade 3: 1                                  |
#> +---------------------------------------------+
#> i  assignment: x$assignment   summary(x) for backend/channel details
#> 
#> $btwn$report_card
#> NULL
#> 
#> $btwn$posterior
#> NULL
#> 
#> $btwn$grade_count
#> [1] 3
#> 
#> $btwn$producer_status
#> [1] "OK"
#> 
#> 
#> $selected_lambda
#> [1] 0.25
#> 
#> $lambda_grid
#> [1] 0.25 1.00
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
#> $method
#> [1] "pairwise"
#> 
#> $producer_status
#> [1] "OK"
#> 
#> $schema_version
#> [1] "v2"
#> 
#> $provenance
#> $provenance$producer
#> [1] "gp_twolevel_grade"
#> 
#> $provenance$n_units
#> [1] 4
#> 
#> $provenance$n_industries
#> [1] 3
#> 
#> $provenance$selected_lambda
#> [1] 0.25
#> 
#> $provenance$route
#> [1] "industry_rfe_Pi_theta__btwn_Pi_bar"
#> 
#> $provenance$m1_safe
#> [1] TRUE
#> 
#> 
#> $warnings
#> character(0)
#> 
#> attr(,"class")
#> [1] "gp_twolevel_grade" "list"             
# }
```
