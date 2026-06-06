# KRW (2024) firm-level callback-discrimination estimates

Per-firm empirical-Bayes input estimates from Kline, Rose and Walters
(2024), "A Discrimination Report Card": the race (White minus Black) and
gender (Male minus Female) callback-rate contact gaps for 97 large U.S.
employers, with their standard errors, plus the firm's industry
grouping. This is the worked example the `gradepath` pipeline runs on;
the race contrast is the M1 headline (97 firms, 3 grades 2/81/14 at
lambda = 0.25).

## Usage

``` r
krw_firms
```

## Format

A data frame with 97 rows and 7 columns:

- firm_id:

  integer. Firm identifier; non-contiguous (matches the KRW
  replication-archive firm numbering, which has gaps).

- theta_hat_race:

  numeric. Estimated White minus Black callback-rate gap (the race
  discrimination contact gap).

- se_race:

  numeric. Standard error of `theta_hat_race`.

- theta_hat_gender:

  numeric. Estimated Male minus Female callback-rate gap (the gender
  contact gap).

- se_gender:

  numeric. Standard error of `theta_hat_gender`.

- industry:

  character. SIC-grouped industry code (used by the two-level / industry
  workflow).

- label:

  character. Human-readable firm label.

## Source

Kline, P., Rose, E. K., and Walters, C. R. (2024). A Discrimination
Report Card. *American Economic Review*, 114(8), 2472–2525.
[doi:10.1257/aer.20230700](https://doi.org/10.1257/aer.20230700) .
Estimates via
[`ebrecipe::krw_firms`](https://rdrr.io/pkg/ebrecipe/man/krw_firms.html);
industry grouping from the replication archive.

## Details

The estimates are taken verbatim from
[`ebrecipe::krw_firms`](https://rdrr.io/pkg/ebrecipe/man/krw_firms.html);
the `industry` and `label` columns are joined from the KRW replication
metadata on `firm_id`. A `sample_stats` attribute records the
sampling-frame counts (108 firms / 83,643 observations before filtering;
97 firms / 78,910 after).

## Provenance and parity (read before replicating KRW)

WARNING – `krw_firms` is a PUBLIC EXAMPLE dataset, NOT the KRW parity
input. It does not reproduce the published KRW (2024) results and must
not be used for replication. The bundled series is on a different
numeric scale from the Matlab GMM input KRW actually consumed: running
the beta-GMM on `krw_firms` yields a spurious large-beta optimum (race
beta = 2.1 and gender beta = 3.0, versus the published 0.51 and 1.26),
and the gender path errors out in the deconvolution with
`DECONV_BOUNDARY_ERROR`. It therefore does NOT reproduce KRW Table 3.

The parity-faithful input is KRW's real Matlab GMM series, shipped at
`inst/extdata/krw-gmm-input/` and read directly with, for example,
`read.csv(system.file("extdata/krw-gmm-input/theta_estimates_matlab_race.csv", package = "gradepath"), header = FALSE)`
(column 2 = `theta_hat`, column 3 = `s`). Use that series (not
`krw_firms`) whenever you need to reproduce KRW. `krw_firms` is retained
only as a small, public, end-to-end worked example of the pipeline's
mechanics.

## Examples

``` r
data(krw_firms)
str(krw_firms)
#> 'data.frame':    97 obs. of  7 variables:
#>  $ firm_id         : int  1 2 3 4 5 7 8 9 10 11 ...
#>  $ theta_hat_race  : num  0.04696 0.022 0.04216 0.00571 0.03408 ...
#>  $ se_race         : num  0.0162 0.0153 0.023 0.015 0.0215 ...
#>  $ theta_hat_gender: num  -0.02287 0.058 -0.09103 0.01499 -0.00699 ...
#>  $ se_gender       : num  0.0251 0.032 0.0355 0.0249 0.0252 ...
#>  $ industry        : chr  "59" "58" "56" "72-73" ...
#>  $ label           : chr  "Rite Aid" "KFC" "Victoria's Secret" "United Rentals" ...
#>  - attr(*, "sample_stats")=List of 5
#>   ..$ full_observations    : int 83643
#>   ..$ full_firms           : int 108
#>   ..$ dropped_observations : int 4733
#>   ..$ filtered_firms       : int 97
#>   ..$ filtered_observations: int 78910
```
