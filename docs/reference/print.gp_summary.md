# Print a `gp_summary`

Shared printer for the typed summary objects returned by
[`summary()`](https://rdrr.io/r/base/summary.html) on a gradepath
result. Renders the headline scalars first, then a `provenance:` block
(backend, channel, selected-lambda rule) – provenance appears ONLY in
the summary, never in [`print()`](https://rdrr.io/r/base/print.html) of
the object itself.

## Usage

``` r
# S3 method for class 'gp_summary'
print(x, ...)
```

## Arguments

- x:

  A `gp_summary` (e.g. from `summary(gp_frontier_object)`).

- ...:

  Unused.

## Value

`x`, invisibly.
