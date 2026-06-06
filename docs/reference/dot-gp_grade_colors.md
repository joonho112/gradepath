# Ordinal colour values for grade levels

Internal helper returning `n_grades` colours along a perceptually
uniform, colour-blind-safe viridis ramp (via
[`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html))
so that grade 1 (the most-extreme theta block) is the dark,
high-contrast end of the scale. Shared by every gradepath figure so
grade colours stay identical across plots.

## Usage

``` r
.gp_grade_colors(n_grades)
```

## Arguments

- n_grades:

  Positive integer: the number of grade levels.

## Value

A character vector of `n_grades` hex colours, named `"1"` ..
`"n_grades"`, ordered from grade 1 to grade `n_grades`.
