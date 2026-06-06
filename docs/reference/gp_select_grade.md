# Select a stored grade fit from a solved grade path

`gp_select_grade()` returns the pre-solved `gp_grade_fit` at a chosen
penalty from a
[gp_grade_path](https://joonho112.github.io/gradepath/reference/gp_grade_path.md).
It never runs a solver: it indexes into the path's stored fits, so it is
the instant way to read off one assignment after
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
has done the work. If `lambda` is absent from the solved grid, it errors
and reports the available grid values.

## Usage

``` r
gp_select_grade(path, lambda = 0.25)
```

## Arguments

- path:

  A
  [gp_grade_path](https://joonho112.github.io/gradepath/reference/gp_grade_path.md)
  object; the solved frontier to select from.

- lambda:

  Numeric scalar; the penalty to select, which must be present in
  `path$lambda_grid`. Defaults to the KRW baseline `0.25`.

## Value

The stored `gp_grade_fit` at `lambda` (a list of class
`c("gp_grade_fit", "list")`), exactly as solved on the path – no
re-solve. See
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md)
for the slot-level description of a `gp_grade_fit`.

## See also

[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md),
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md),
[`gp_pairwise()`](https://joonho112.github.io/gradepath/reference/gp_pairwise.md)

Other gradepath-grade:
[`gp_grade()`](https://joonho112.github.io/gradepath/reference/gp_grade.md),
[`gp_grade_path()`](https://joonho112.github.io/gradepath/reference/gp_grade_path.md),
[`gp_preview()`](https://joonho112.github.io/gradepath/reference/gp_preview.md)

## Examples

``` r
# Solve-free: pull the baseline-penalty fit from the bundled tiny fit's path.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
gp_select_grade(fit$grade_path, lambda = 0.25)
#> +---------------------------------------------+
#> | gp_grade_fit  .  2 grades  .  lambda = 0.25 |
#> +---------------------------------------------+
#> | units = 24   status = optimal               |
#> | grade 1: 21                                 |
#> | grade 2: 3                                  |
#> +---------------------------------------------+
#> i  assignment: x$assignment   summary(x) for backend/channel details

# The pipeline already stored this selection; it is the same object.
identical(gp_select_grade(fit$grade_path, lambda = 0.25), fit$selected_grade)
#> [1] TRUE
```
