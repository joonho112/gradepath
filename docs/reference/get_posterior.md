# Extract the native posterior from a gradepath fit

Returns the `gp_posterior` slot of a `gp_fit` unchanged. The native
`gp_posterior` is `eb_posterior`-shaped (per-unit posterior means,
standard deviations, and credible intervals) but remains a gradepath
object rather than inheriting from `eb_posterior`. It holds the per-unit
posterior summaries that the report card presents.

## Usage

``` r
get_posterior(fit, ...)

# Default S3 method
get_posterior(fit, ...)

# S3 method for class 'gp_fit'
get_posterior(fit, ...)
```

## Arguments

- fit:

  A `gp_fit` object.

- ...:

  Unused; accepted for S3 method extensibility. No arguments are
  forwarded; supplying any has no effect.

## Value

The `gp_posterior` object stored on `fit` (returned unchanged).

## Details

Read-only accessor adopted in `GP-DEC-13-A`; it never recomputes and
only reads the materialized slot.

## See also

[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)

Other gradepath-accessors:
[`as.data.frame.gp_fit()`](https://joonho112.github.io/gradepath/reference/as.data.frame.gp_fit.md),
[`coef.gp_fit()`](https://joonho112.github.io/gradepath/reference/coef.gp_fit.md),
[`get_control()`](https://joonho112.github.io/gradepath/reference/get_control.md),
[`get_grades()`](https://joonho112.github.io/gradepath/reference/get_grades.md),
[`get_pairwise()`](https://joonho112.github.io/gradepath/reference/get_pairwise.md),
[`get_prior()`](https://joonho112.github.io/gradepath/reference/get_prior.md),
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)

## Examples

``` r
# Instant: read the native posterior off the bundled tiny fit.
fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
get_posterior(fit)     # native gp_posterior (eb_posterior-shaped)
#> $estimate
#>  [1]  0.411152808 -0.071090509  0.473473744  0.112719140  0.684172551
#>  [6]  0.003436935  0.026503544  0.173350995  0.075307252  0.471463554
#> [11] -0.652549945  0.285023870 -0.020468156 -0.163023397  0.209276813
#> [16]  0.088788468  0.213391253  0.262980500  0.227317019  0.303749869
#> [21]  0.160705668  0.163141700 -0.033699669  0.313878987
#> 
#> $se
#>  [1] 0.14630406 0.06580824 0.32926504 0.10659518 0.20673635 0.05728224
#>  [7] 0.16072021 0.07435398 0.15302709 0.18597975 0.53288978 0.11539352
#> [13] 0.11379544 0.33182972 0.19612685 0.12072278 0.12143768 0.27619487
#> [19] 0.17835897 0.13480130 0.12527737 0.16999499 0.14886915 0.21755222
#> 
#> $id
#>  [1] "firm01" "firm02" "firm03" "firm04" "firm05" "firm06" "firm07" "firm08"
#>  [9] "firm09" "firm10" "firm11" "firm12" "firm13" "firm14" "firm15" "firm16"
#> [17] "firm17" "firm18" "firm19" "firm20" "firm21" "firm22" "firm23" "firm24"
#> 
#> $label
#>  [1] "firm01" "firm02" "firm03" "firm04" "firm05" "firm06" "firm07" "firm08"
#>  [9] "firm09" "firm10" "firm11" "firm12" "firm13" "firm14" "firm15" "firm16"
#> [17] "firm17" "firm18" "firm19" "firm20" "firm21" "firm22" "firm23" "firm24"
#> 
#> $posterior_mean
#>  [1] 0.23100264 0.01577917 0.18358822 0.13517484 0.24446367 0.02124816
#>  [7] 0.10892519 0.19846637 0.12635147 0.22349138 0.12156350 0.22278278
#> [13] 0.05236847 0.12256963 0.17134316 0.12183574 0.19616020 0.16941551
#> [19] 0.17984562 0.21753547 0.16835229 0.16216836 0.07789707 0.18777959
#> 
#> $posterior_sd
#>  [1] 0.05465029 0.01676989 0.09604267 0.09870348 0.08518254 0.02297699
#>  [7] 0.10089756 0.06648436 0.10231160 0.06876965 0.10532778 0.05589827
#> [13] 0.07183101 0.10481905 0.09644554 0.10022565 0.07954877 0.09911810
#> [19] 0.09242694 0.06391997 0.09387572 0.09843229 0.09051652 0.09050439
#> 
#> $lower
#>  [1] 0.1410806262 0.0006848574 0.0061637167 0.0034242870 0.1486140577
#>  [6] 0.0006848574 0.0020545722 0.0301337260 0.0027394296 0.0404065871
#> [11] 0.0020545722 0.0821828891 0.0013697148 0.0020545722 0.0047940019
#> [16] 0.0027394296 0.0123274334 0.0047940019 0.0061637167 0.0362974427
#> [21] 0.0054788593 0.0041091445 0.0013697148 0.0075334315
#> 
#> $upper
#>  [1] 0.28490068 0.04794002 0.27805211 0.26230039 0.29859783 0.06369174
#>  [7] 0.26024582 0.26777925 0.26435496 0.28490068 0.26709439 0.28079154
#> [13] 0.22737266 0.26640953 0.27394296 0.26093067 0.27531268 0.27531268
#> [19] 0.27531268 0.28079154 0.27051868 0.27188839 0.24928810 0.27736725
#> 
#> $scale
#> [1] "r"
#> 
#> $metadata
#> $metadata$level
#> [1] 0.9
#> 
#> $metadata$interval_level
#> [1] 0.9
#> 
#> $metadata$reporting
#> $metadata$reporting$posterior_mean
#>  [1] 0.132313631 0.007169272 0.133036141 0.070633949 0.154787945 0.009273449
#>  [7] 0.064113427 0.093422339 0.073320421 0.137233841 0.101285845 0.119120184
#> [13] 0.027888021 0.089019395 0.106845485 0.066002933 0.106449307 0.116667009
#> [19] 0.109101868 0.121676684 0.092187176 0.097017697 0.044843263 0.120667981
#> 
#> $metadata$reporting$posterior_sd
#>  [1] 0.031302581 0.007619408 0.069596764 0.051576289 0.053935337 0.010027972
#>  [7] 0.059388360 0.031295604 0.059370339 0.042227684 0.087758360 0.029888363
#> [13] 0.038252492 0.076127573 0.060141128 0.054295946 0.043168345 0.068257106
#> [19] 0.056070040 0.035753109 0.051404927 0.058887403 0.052107944 0.058158515
#> 
#> $metadata$reporting$lower
#>  [1] 0.0808081249 0.0003111652 0.0044665015 0.0017893191 0.0940984985
#>  [6] 0.0002988960 0.0012093223 0.0141845857 0.0015896620 0.0248114767
#> [11] 0.0017118551 0.0439425383 0.0007294205 0.0014921867 0.0029894246
#> [16] 0.0014840504 0.0066896687 0.0033013616 0.0037391680 0.0203026771
#> [21] 0.0030001407 0.0024583077 0.0007885083 0.0048410158
#> 
#> $metadata$reporting$upper
#>  [1] 0.16318534 0.02178157 0.20148885 0.13706184 0.18906426 0.02779733
#>  [7] 0.15318083 0.12604939 0.15340238 0.17494194 0.22254116 0.15013701
#> [13] 0.12108380 0.19348688 0.17082426 0.14135580 0.14940260 0.18959248
#> [19] 0.16701617 0.15705845 0.14813195 0.16265803 0.14350851 0.17823740
#> 
#> $metadata$reporting$scale
#> [1] "theta"
#> 
#> $metadata$reporting$level
#> [1] 0.9
#> 
#> 
#> $metadata$has_reporting
#> [1] TRUE
#> 
#> 
#> attr(,"class")
#> [1] "gp_posterior" "list"        
```
