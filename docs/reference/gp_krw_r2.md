# KRW grade R-squared from the rsquared.do recipe

`gp_krw_r2()` ports the Kline-Rose-Walters companion script
`code/rsquared.do`. It builds the bias-corrected denominator from the
source `theta_estimates_*` columns `log_dif` and `log_dif_se`, and the
between-grade numerator from the ranking-output posterior second moments
(`post_mean_beta`, `pmean_sq`). This is the paper-faithful producer for
the firm-level KRW race/gender R-squared rows; for the quick generic
variance-share diagnostic use
[`gp_r2()`](https://joonho112.github.io/gradepath/reference/gp_r2.md).

## Usage

``` r
gp_krw_r2(theta, ranking, grade_col = NULL, scale = c("percent", "proportion"))
```

## Arguments

- theta:

  Data frame shaped like `theta_estimates_race.csv` /
  `theta_estimates_gender.csv`: it must carry `log_dif` and (strictly
  positive) `log_dif_se`, and optionally `obs_idx`. When `obs_idx` is
  absent, row order defines it, matching `rsquared.do`.

- ranking:

  Data frame shaped like `ranking_results_*_*.csv`: it must carry
  `obs_idx`, `post_mean_beta`, `pmean_sq`, and a grade column. Joined to
  `theta` on `obs_idx`; every `theta` id must be present here.

- grade_col:

  Character scalar naming the grade column in `ranking`, or `NULL`
  (default) to auto-detect – it accepts `grades_lamb0.25` and the
  Stata-style `grades_lamb025`.

- scale:

  Character scalar; `"percent"` (default) puts the returned `r2` column
  on a 0-100 scale, `"proportion"` on a 0-1 scale.

## Value

A one-row data frame: `n` (units) and `grade_count`, the resolved
`grade_col`, the variance components `overall_variance`/`overall_sd` and
`between_variance`/`between_sd`, the `r2_proportion`, and `r2` on the
requested `scale` (recorded in the `scale` column).

## Details

The denominator subtracts the mean sampling variance from the raw
estimate variance (the bias correction); the between-grade numerator
applies the posterior second-moment correction from the ranking columns.
The names pipeline uses different probability-scale quantities and is
deferred to the `M1_NAMES` milestone, so this producer is firm-level
race/gender only.

## See also

[`gp_r2()`](https://joonho112.github.io/gradepath/reference/gp_r2.md),
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md),
[`gp_report_card()`](https://joonho112.github.io/gradepath/reference/gp_report_card.md)

Other gradepath-frontier:
[`gp_frontier()`](https://joonho112.github.io/gradepath/reference/gp_frontier.md),
[`gp_r2()`](https://joonho112.github.io/gradepath/reference/gp_r2.md)

## Examples

``` r
# The KRW race/gender inputs are CSV products of the published solve, not part
# of the tiny example fixture; build the frame columns this recipe expects.
theta <- data.frame(
  log_dif    = c(-0.4, -0.1, 0.2, 0.5),
  log_dif_se = c(0.10, 0.12, 0.09, 0.11)
)
ranking <- data.frame(
  obs_idx         = 1:4,
  post_mean_beta  = c(-0.35, -0.08, 0.18, 0.46),
  pmean_sq        = c(0.130, 0.012, 0.038, 0.220),
  grades_lamb0.25 = c(1L, 2L, 2L, 3L)
)
gp_krw_r2(theta, ranking)
#>   n grade_count       grade_col overall_variance overall_sd between_variance
#> 1 4           3 grades_lamb0.25          0.13885  0.3726258           0.0857
#>   between_sd r2_proportion       r2   scale
#> 1  0.2927456     0.6172128 61.72128 percent
```
