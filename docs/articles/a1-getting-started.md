# A1: Getting started — your first discrimination report card

Abstract

gradepath turns firm-level empirical-Bayes estimates into a
discrimination *report card*: it ranks units by a posterior outranking
probability, solves a grade integer program along an
information-reliability frontier, and assembles a grade-plus-interval
table you can print, plot, and export. This vignette is a gentle first
contact. In under 60 lines of code you install the package, see the one
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
call you would write, load the bundled RACE report card (97 firms,
grades 2 / 81 / 14), and read it three ways: as a console box, as a tidy
data frame, and as one figure. You leave knowing what input a report
card needs, how to read its grades and reliability honestly, and where
the full workflow lives next.

## 1. Why you are here

You have a vector of unit-level estimates — say, 97 discrimination point
estimates across large U.S. firms — each with its own standard error.
The noisiest estimates look the most extreme purely by chance, so a raw
ordering mistakes imprecision for evidence. How do you sort units by
*what the data actually support*, and report each one with an honest
interval, without overstating a difference the noise cannot bear?

> “Estimates of firm-level discrimination are inherently noisy … A
> report card summarizes this evidence by sorting firms into a small
> number of grades.” — Kline, Rose & Walters (Kline et al. 2024)

gradepath is a native-R replication companion for that paper. This
vignette runs the full report-card workflow in one conceptual call, then
unpacks what just happened.

## 2. What you will leave with

After this vignette, you will be able to:

- Install **gradepath** and load a bundled report card with one call.
- Read the package’s console box and locate the grade composition, the
  reliability `(1 - DR)`, and `tau-bar`.
- Pull the report card into your own workflow as a tidy data frame and
  export it to CSV.
- Read the one signature figure — each firm’s posterior gap with its
  interval, ranked and coloured by grade.
- Understand, in plain steps, the three-part KRW pipeline (estimate →
  rank → grade) without touching a derivation.

Reading + reflection ≈ 25–35 min. It runs in seconds, because it loads a
pre-computed fit rather than solving anything.

## 3. Install and load

``` r

# install.packages("remotes")
remotes::install_github("joonho112/gradepath")
```

The setup chunk above already attached **gradepath** and **ggplot2** and
set the package theme, so you do not need to load anything else. For
your own scripts the one line to remember is:

``` r

library(gradepath)
```

## 4. The input — what a report card needs

A report card needs very little: a vector of per-firm gaps and a vector
of their standard errors. Each gap is one firm’s estimated
discrimination (here, the White-minus-Black callback gap); each standard
error says how precisely that gap was measured.

The package ships a public example, `krw_firms`, that carries exactly
this shape (97 firms with race and gender gaps and SEs). It is handy for
trying the API, but it is an example on a *different numeric scale* and
does **not** reproduce the paper’s grades — so we do not drive the
headline off it. The figures in this vignette come from the bundled
parity input: KRW’s real per-firm beta-GMM series, shipped under
`inst/extdata/krw-gmm-input/`, which is what reproduces the published
RACE result.

## 5. The one call

The whole workflow is one function,
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md).
You assemble the input from the shipped beta-GMM series and hand it over
with a control:

``` r

# Read the real per-firm beta-GMM series shipped with the package (97 firms,
# headerless: column 2 = theta_hat gap, column 3 = s standard error).
gmm  <- read.csv(system.file("extdata/krw-gmm-input/theta_estimates_matlab_race.csv",
                             package = "gradepath"), header = FALSE)
race <- list(theta_hat = gmm[[2]], s = gmm[[3]])

fit <- krw_report_card(
  race, demographic = "race",
  control = gp_control(backend = "gurobi", precision_rule = "krw_gmm",
                       lambda_grid = c(0.25, 0.5, 1)))   # the grid must include 1
```

That single call runs all seven stages in order: beta-GMM precision
standardization → deconvolution of the prior → posterior shrinkage → the
pairwise outranking matrix → the grade integer program → the
information-reliability frontier → the assembled report card.

The catch is the fifth stage. The grade step **solves an integer program
over 97 firms**, which takes roughly three minutes with Gurobi. So we do
not run it in this vignette — gradepath ships the proven result, and we
load it here:

``` r

fit <- gp_parity_fit("race")
```

This bundled object is exactly what the call above produces: the RACE
report card with grades 2 / 81 / 14 at the KRW baseline `lambda = 0.25`.
Everything below reads off it; nothing below re-solves.

## 6. Read it: the console

The fastest read is to print the fit. The console box is result-first:
how many units, the grade composition, the selected `lambda`, the
reliability, and the backend.

``` r

fit
#> +------------------------------------------------------------------------+
#> | gp_fit  .  97 units  .  grades: 3 (2/81/14)  .  selected lambda = 0.25 |
#> +------------------------------------------------------------------------+
#> | units      : 97                                                        |
#> | grades     : 3 (2/81/14)                                               |
#> | reliability: (1 - DR) = 0.96   tau-bar = 0.21                          |
#> | backend    : gurobi                                                    |
#> +------------------------------------------------------------------------+
#> i  summary(fit) for backend/selection details; get_report_card(fit) for the ranked table.
```

What to look at:

- **`grades`** — `3 (2/81/14)`: 97 firms sorted into three grades, with
  2 in the most-extreme grade, 81 in the middle, and 14 in the
  least-extreme grade.
- **`reliability: (1 - DR)`** — `0.96`: the share of comparable firm
  pairs on which the grades agree with the underlying posterior ranking.
  Higher is cleaner.
- **`tau-bar`** — `0.21`: the posterior expected Kendall rank agreement
  the grading carries across grade boundaries (between-grade blocks), in
  \[-1, 1\]; higher means the grades preserve more rank information.

For the plain-English read, call
[`summary()`](https://rdrr.io/r/base/summary.html). On top of the same
numbers it adds an **Interpretation** block (what the grades,
reliability, and `tau-bar` mean in words), a **Glossary** (each term
defined), and a provenance block (backend and selection rule):

``` r

summary(fit)
#> <summary: gp_fit>
#>   units           : 97
#>   grade_count     : 3
#>   composition     : 3 (2/81/14)
#>   selected_lambda : 0.250
#>   reliability     : 0.961
#>   tau_bar         : 0.207
#>   Interpretation:
#>     97 units sorted into 3 grades (2/81/14) at lambda = 0.25.
#>     Grades agree with the pairwise (posterior) ranking on 96% of comparable
#>     pairs (reliability = 1 - discordance rate = 0.961); higher is cleaner.
#>     The grades carry tau-bar = 0.207 of posterior Kendall rank agreement across
#>     grade boundaries (in [-1, 1]); higher means more rank information.
#>   Glossary:
#>     reliability     : 1 - discordance rate: share of comparable pairs the gra...
#>     tau_bar         : posterior expected Kendall rank agreement over between-...
#>     selected_lambda : the grade-count selection knob (KRW baseline 0.25); lar...
#> 
#>   provenance:
#>     backend         : gurobi
#>     selection_rule  : baseline_lambda_0.25
#>     frontier_source : grade_path$summary (enriched)
```

The Interpretation block states the headline directly: the grades agree
with the pairwise (posterior) ranking on 96% of comparable pairs. Grades
are a sorting by posterior evidence, never a contest, and grade 1 is not
a verdict — it is simply the most-extreme grade under the data.

## 7. Read it: the table

For joins, regressions, or plots in your own style, pull the report card
into a data frame. `as.data.frame(fit)` returns a tidy 97 × 10 table,
ranked with the most-extreme unit (grade 1) first:

``` r

head(as.data.frame(fit), 8)
#>             id        label grade sort_rank selected_lambda posterior_mean
#> 1 krw_race_062 krw_race_062     1         1            0.25      0.2498718
#> 2 krw_race_015 krw_race_015     1         2            0.25      0.2328297
#> 3 krw_race_017 krw_race_017     2         3            0.25      0.1907685
#> 4 krw_race_037 krw_race_037     2         4            0.25      0.1877637
#> 5 krw_race_095 krw_race_095     2         5            0.25      0.1870007
#> 6 krw_race_010 krw_race_010     2         6            0.25      0.1775977
#> 7 krw_race_085 krw_race_085     2         7            0.25      0.1698248
#> 8 krw_race_088 krw_race_088     2         8            0.25      0.1672324
#>        lower     upper  estimate        se
#> 1 0.11742058 0.3424767 0.3303573 0.0713513
#> 2 0.08092269 0.4474549 0.4332195 0.1308592
#> 3 0.03242811 0.3827927 0.3772942 0.2828115
#> 4 0.03993012 0.3950010 0.4022631 0.2158661
#> 5 0.03117536 0.3726872 0.3499199 0.2856306
#> 6 0.03232195 0.3523093 0.3431389 0.2385807
#> 7 0.02297781 0.3417023 0.1686227 0.3120664
#> 8 0.02915253 0.3276491 0.2950459 0.2294773
```

The columns are `id`, `label`, `grade`, `sort_rank`, `selected_lambda`,
`posterior_mean`, `lower`, `upper`, `estimate`, and `se`. Compare
`estimate` (the raw gap) with `posterior_mean` (the shrunken gap) to
watch empirical Bayes work: the noisier a firm’s estimate, the more its
posterior mean is pulled toward the center. The interval
`[lower, upper]` is the credible interval carried on the card; it is
reported as-is and never recomputed downstream.

Two one-line reads confirm the headline numbers:

``` r

round(get_prior(fit)$metadata$beta, 4)   # the precision exponent
#> [1] 0.5095
table(get_grades(fit))                    # the grade composition
#> 
#>  1  2  3 
#>  2 81 14
```

The precision exponent is $`\beta`$ in $`\theta = s^{\beta}\, v`$ —
gradepath lands on `0.5095`, matching the `0.510` KRW report. The grade
table is the `2 / 81 / 14` composition again.

Exporting is a one-liner:

``` r

write.csv(as.data.frame(fit), "report_card_race.csv", row.names = FALSE)
```

## 8. Read it: one figure

The single most informative figure is the report card itself (KRW Figure
7): each firm’s posterior gap shown as a point with its interval, one
firm per row, ranked with the most-extreme firm on top and coloured by
grade.

``` r

gp_plot_report_card(get_report_card(fit))
```

![A horizontal caterpillar plot of 97 firms, one per row, ranked from
most extreme at the top to least extreme at the bottom. Each row shows a
point for the firm's posterior gap and a horizontal line for its
credible interval. Points are coloured by grade into three groups: a
small most-extreme grade at the top, a large middle grade, and a
least-extreme grade at the bottom, in the 2 / 81 / 14
split.](a1-getting-started_files/figure-html/figure-report-card-1.png)

RACE report card (KRW Figure 7): each of the 97 firms shown as its
posterior White-minus-Black callback gap (point) with its credible
interval (horizontal line), one firm per row, ranked with the
most-extreme firm on top and coloured by its grade (2 / 81 / 14).

Reading the figure: the vertical position is the firm’s rank; the
horizontal position is the posterior gap; the line length is the firm’s
uncertainty. The colour blocks are the grades. Because the figure shows
each interval, you can see at a glance which firms are separated by the
data and which overlap — the grades respect that separation rather than
the raw point order.

## 9. What you just did, formally

The KRW pipeline is three plain steps, no equations required.

1.  **Estimate.** Take each firm’s raw gap and standard error and put
    them on a common precision scale via the beta-GMM standardization
    (the precision exponent $`\beta`$ above).
2.  **Rank.** Deconvolve the prior distribution of true gaps, shrink
    each firm’s estimate toward it to get a posterior, and compare every
    pair of firms by the posterior probability that one outranks the
    other.
3.  **Grade.** Solve an integer program that sorts firms into a small
    number of grades along the information-reliability frontier,
    selecting the grade count at the KRW baseline $`\lambda = 0.25`$.

Grades are a sorting by posterior evidence, not a ranking contest: a
grade label is an integer in `1 … n` that carries no superiority
statement. The one-level race and gender core is accepted within its
declared scope, and the numbers here match the paper ($`\beta = 0.5095`$
vs the reported `0.510`; $`\lambda = 0.25`$ is the KRW baseline).

## 10. Where to next

- **Full workflow** — [the grading
  workflow](https://joonho112.github.io/gradepath/articles/a2-the-grading-workflow.md)
  walks the seven stages individually and shows race and gender side by
  side.
- **Reading and exporting** —
  [`vignette("a3-reading-and-exporting-results")`](https://joonho112.github.io/gradepath/articles/a3-reading-and-exporting-results.md)
  goes deeper on the accessors, the table, and CSV/round-trip export.
- **Figures** —
  [`vignette("a4-figures")`](https://joonho112.github.io/gradepath/articles/a4-figures.md)
  builds all four KRW signature figures.
- **Solvers and calibration** —
  [`vignette("a5-solvers-and-calibration")`](https://joonho112.github.io/gradepath/articles/a5-solvers-and-calibration.md)
  covers the open-source backends and the calibration harness.
- **Foundations** —
  [`vignette("m1-foundations-and-notation")`](https://joonho112.github.io/gradepath/articles/m1-foundations-and-notation.md)
  formalizes the estimate → rank → grade pipeline with the math this
  vignette deliberately skipped.

### Provenance

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US/en_US/en_US/C/en_US/en_US
#> 
#> time zone: America/Chicago
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.3   gradepath_0.5.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5       ebrecipe_0.5.0     gtable_0.3.6       jsonlite_2.0.0    
#>  [5] dplyr_1.2.1        compiler_4.6.0     tidyselect_1.2.1   slam_0.1-55       
#>  [9] dichromat_2.0-0.1  jquerylib_0.1.4    splines_4.6.0      systemfonts_1.3.2 
#> [13] scales_1.4.0       textshaping_1.0.5  yaml_2.3.12        fastmap_1.2.0     
#> [17] lattice_0.22-9     R6_2.6.1           labeling_0.4.3     generics_0.1.4    
#> [21] knitr_1.51         htmlwidgets_1.6.4  gurobi_13.0-1      tibble_3.3.1      
#> [25] desc_1.4.3         bslib_0.11.0       pillar_1.11.1      RColorBrewer_1.1-3
#> [29] rlang_1.2.0        cachem_1.1.0       xfun_0.57          fs_2.1.0          
#> [33] sass_0.4.10        S7_0.2.2           otel_0.2.0         cli_3.6.6         
#> [37] withr_3.0.2        pkgdown_2.2.0      magrittr_2.0.5     digest_0.6.39     
#> [41] grid_4.6.0         lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5    
#> [45] glue_1.8.1         farver_2.1.2       ragg_1.5.2         rmarkdown_2.31    
#> [49] tools_4.6.0        pkgconfig_2.0.3    htmltools_0.5.9
```

## References

Kline, Patrick, Evan K. Rose, and Christopher R. Walters. 2024. “A
Discrimination Report Card.” *American Economic Review* 114 (8):
2472–525. <https://doi.org/10.1257/aer.20230700>.
