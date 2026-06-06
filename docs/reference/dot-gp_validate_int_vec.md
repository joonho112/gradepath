# Integer-vector validator composed from the built numeric-vector primitive

The brief named a `.gradepath_validate_int_vec`; it does not exist. We
compose the contract (numeric, finite, length n, integer-valued, \>=
min) from the real `.gradepath_validate_numeric_vector`, returning an
`integer` vector. Internal to the output-object layer.

## Usage

``` r
.gp_validate_int_vec(x, n = NULL, min = 1L, what = "value")
```
