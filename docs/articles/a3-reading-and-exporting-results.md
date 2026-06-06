# A3: Reading and exporting results

Abstract

A fitted gp_fit is a rich object: it carries the grades, the deconvolved
prior, the per-firm posterior, the pairwise outranking matrix, the
assembled report card, and the run-control that produced them. This
vignette is the reference for getting all of that out and into your own
workflow. We read the fit several ways for the console and the page (the
result-first print box, the fuller summary with its Interpretation and
Glossary, the tidy ten-column data frame, the six get\_\* accessors, and
coef for a quick firm-to-grade lookup), then export it three ways (a CSV
you drop into a paper table, an RDS that round-trips a gp_fit exactly,
and an ASCII card you embed in a log via the format\_\*\_cli renderers).
Everything runs on the bundled RACE fit, and nothing here re-solves:
every read is a fast lookup off a loaded object.

## 1. What a fit holds

By the end of
[A2](https://joonho112.github.io/gradepath/articles/a2-the-grading-workflow.md)
you have a `gp_fit`: a single object that bundles everything the KRW
pipeline produced (Kline et al. 2024) — the deconvolved prior, the
per-firm posterior, the $`J \times J`$ pairwise outranking matrix, the
selected grading, the assembled report card, and the run-control that
drove the solve. This vignette is the reference for reading that object
and exporting it, with each read pointed at a concrete downstream use:
joining the card to firm covariates, dropping a block into a paper
table, handing off a reproducible artifact, or pasting a result-first
card into a build log.

We load the proven RACE fit **once** and read everything below off it.
Loading is instant because the fit ships pre-solved; **nothing in this
vignette re-solves**. Every read is a fast lookup off materialized slots
— the accessors never call the solver and recompute nothing.

``` r

fit <- gp_parity_fit("race")
```

For the record, the one call that *produces* a fit in the first place is
[`krw_report_card()`](https://joonho112.github.io/gradepath/reference/krw_report_card.md)
— it solves the integer program and takes minutes, so it is shown here
only as a pointer and never run:

``` r

fit <- krw_report_card(
  race, demographic = "race",
  control = gp_control(backend = "gurobi", precision_rule = "krw_gmm",
                       lambda_grid = c(0.25, 0.5, 1)))
```

## 2. The quick read: `print()` / the console box

The fastest read is to print the fit.
[`print()`](https://rdrr.io/r/base/print.html) (which you get by typing
`fit` alone) emits a result-first ASCII box — it leads with the answer,
not the inputs:

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

Line by line, this is the headline you would quote in a slide or a
commit message:

- The header —
  `gp_fit . 97 units . grades: 3 (2/81/14) . selected lambda = 0.25`: 97
  firms sorted into three grades, with 2 in the most-extreme grade, 81
  in the middle, and 14 in the least-extreme grade, at the KRW baseline
  $`\lambda = 0.25`$.
- `reliability: (1 - DR)` — `0.96`: the share of comparable firm pairs
  on which the grades agree with the underlying posterior ordering.
  Higher is cleaner.
- `tau-bar` — `0.21`: the posterior expected Kendall rank agreement the
  grading carries across grade boundaries (between-grade blocks), in
  \[-1, 1\]; higher means the grades preserve more rank information.
- `backend` — `gurobi`: the solver that produced the grading.

A footer line points you onward —
[`summary()`](https://rdrr.io/r/base/summary.html) for the full read,
[`get_report_card()`](https://joonho112.github.io/gradepath/reference/get_report_card.md)
for the ranked table. Grade 1 is simply the most-extreme $`\theta`$
block under the data, not a verdict.

## 3. The full read: `summary()`

When you are writing up the result rather than glancing at it,
[`summary()`](https://rdrr.io/r/base/summary.html) is the read you want.
It keeps every number from the box and adds three plain-English blocks;
it returns a typed `gp_fit_summary`, and printing it shows the full
read:

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

The three additions:

- **Interpretation** — the headline in words you can lift into a methods
  paragraph: the grades agree with the pairwise (posterior) ordering on
  about 96% of comparable pairs, and grade 1 is the most-extreme grade,
  never a contest.
- **Glossary** — each term defined in place: `reliability`
  $`(1 - \text{DR})`$, `tau_bar`, and `selected_lambda`, so a reader who
  has never opened the package can follow the output.
- **Provenance** — the backend (`gurobi`) and the selection rule
  (`baseline_lambda_0.25`), exactly what you cite to make the run
  reproducible.

Use [`print()`](https://rdrr.io/r/base/print.html) for a glance and
[`summary()`](https://rdrr.io/r/base/summary.html) when you need the
reading and the provenance on the page.

## 4. The tidy table: `as.data.frame(fit)`

For anything analytical — a join to firm covariates, a regression of
grade on firm size, a plot in your own house style — coerce the fit to a
data frame. `as.data.frame(fit)` returns the report-card rows as a tidy
`97 × 10` table, ranked with the most-extreme firm (grade 1) on top and
row names reset:

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

The ten columns are the working surface for downstream use:

- `id`, `label` — the firm identifier and its display label, your **join
  keys** against any external firm table.
- `grade` — the integer grade in $`\{1, \dots, n\}`$ (grade 1 is most
  extreme).
- `sort_rank` — the rank used to order the rows; order on it, facet on
  `grade`.
- `selected_lambda` — the penalty that produced this grading (`0.25`),
  carried on every row so the table is self-documenting.
- `posterior_mean`, `lower`, `upper` — the shrunken gap and its credible
  interval, reported as-is and never recomputed downstream.
- `estimate`, `se` — the raw gap and its standard error. Compare
  `estimate` with `posterior_mean` to watch empirical Bayes work: the
  noisier a firm’s estimate, the harder its posterior is pulled toward
  the centre.

Because it is an ordinary data frame, this is the object you
[`merge()`](https://rdrr.io/r/base/merge.html),
[`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html),
or hand to [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html)
for a paper table. The report-card object coerces to the identical ten
columns, so you can start from either end:

``` r

identical(as.data.frame(get_report_card(fit)), as.data.frame(fit))
#> [1] TRUE
```

## 5. The accessors: the six `get_*` functions

A `gp_fit` is read through its accessor family, never by reaching into
`$` slots. Each accessor is a read-only lookup that recomputes nothing —
reach for the one whose output your downstream step needs. There are
six.

**`get_grades(fit)`** — the selected grading as a named integer vector
of length 97, names equal to the firm ids. This is what you
[`cbind()`](https://rdrr.io/r/base/cbind.html) onto a firm table or pass
to [`table()`](https://rdrr.io/r/base/table.html) for the composition:

``` r

head(get_grades(fit))
#> krw_race_001 krw_race_002 krw_race_003 krw_race_004 krw_race_005 krw_race_006 
#>            2            2            2            2            2            3
table(get_grades(fit))      # the grade composition: 2 / 81 / 14
#> 
#>  1  2  3 
#>  2 81 14
```

**`get_report_card(fit)`** — the assembled `gp_report_card` object, the
source of the tidy table in section 4. Reach for it when you want the
object that prints and plots, not just the rows:

``` r

get_report_card(fit)
#> <gp_report_card>  unit: unit | 97 rows | grades: 3 (2/81/14) | selected lambda: 0.25
#> sorted by Condorcet rank (endpoint lambda = 1; secondary key: id)
#> 
#>   sort_rank  id            label         grade  posterior_mean  CI (90%)        estimate
#>    1         krw_race_062  krw_race_062  1      0.250           [0.117, 0.342]  0.330   
#>    2         krw_race_015  krw_race_015  1      0.233           [0.081, 0.447]  0.433   
#>    3         krw_race_017  krw_race_017  2      0.191           [0.032, 0.383]  0.377   
#>    4         krw_race_037  krw_race_037  2      0.188           [0.040, 0.395]  0.402   
#>    5         krw_race_095  krw_race_095  2      0.187           [0.031, 0.373]  0.350   
#>    6         krw_race_010  krw_race_010  2      0.178           [0.032, 0.352]  0.343   
#>    7         krw_race_085  krw_race_085  2      0.170           [0.023, 0.342]  0.169   
#>   ... 89 more rows ...
#>   97         krw_race_006  krw_race_006  3      0.017           [0.001, 0.044]  -0.032  
#> 
#> i  grade 1 = most-extreme theta (firms: most discriminatory; names: best-treated). See ?gp_report_card and Appendix A.
#> i  full table: as.data.frame(card)   formatted CLI: format_gp_report_card_cli(card)
```

**`get_posterior(fit)`** — the native `gp_posterior` (per-firm means,
standard deviations, and intervals). Its data frame carries `estimate`,
`se`, `id`, `label`, `posterior_mean`, `posterior_sd`, `lower`, `upper`,
`scale`, and more; `estimate` versus `posterior_mean` is the shrinkage:

``` r

head(as.data.frame(get_posterior(fit))[, c("id", "estimate", "posterior_mean", "posterior_sd")])
#>             id    estimate posterior_mean posterior_sd
#> 1 krw_race_001  0.83306863      0.4930476    0.2219116
#> 2 krw_race_002  0.31162778      0.2957276    0.1429102
#> 3 krw_race_003  0.47406122      0.3543188    0.1605673
#> 4 krw_race_004  0.09661802      0.2536574    0.1528156
#> 5 krw_race_005  0.43743703      0.3302823    0.1713923
#> 6 krw_race_006 -0.19290522      0.1021265    0.0815843
```

**`get_prior(fit)`** — the deconvolved native `gp_prior`, a list
carrying `support`, `density`, `mean`, `scale`, `diagnostics`, and
`metadata`. The precision exponent $`\beta`$ — the number you cite to
show the replication lands on the paper — lives in the metadata:

``` r

names(get_prior(fit))
#> [1] "support"     "density"     "mean"        "scale"       "diagnostics"
#> [6] "metadata"
round(get_prior(fit)$metadata$beta, 4)   # the precision exponent beta
#> [1] 0.5095
```

gradepath lands on $`\beta = 0.5095`$, matching the `0.510` KRW report.

**`get_pairwise(fit)`** — the `gp_pairwise` social-choice input the
grading is solved over: the $`J \times J`$ matrix
$`\Pi_{ij} = \Pr(\theta_i > \theta_j)`$ in `$matrix`, with the firm
order in `$ids`:

``` r

pw <- get_pairwise(fit)
pw$matrix[1:5, 1:5]         # a corner of the J x J outranking matrix
#>              krw_race_001 krw_race_002 krw_race_003 krw_race_004 krw_race_005
#> krw_race_001    0.5000000    0.8912488    0.7623452    0.7589266    0.6559360
#> krw_race_002    0.1087512    0.5000000    0.2996195    0.3715744    0.2459148
#> krw_race_003    0.2376548    0.7003805    0.5000000    0.5365897    0.4039615
#> krw_race_004    0.2410734    0.6284256    0.4634103    0.5000000    0.3848517
#> krw_race_005    0.3440640    0.7540852    0.5960385    0.6151483    0.5000000
head(pw$ids)
#> [1] "krw_race_001" "krw_race_002" "krw_race_003" "krw_race_004" "krw_race_005"
#> [6] "krw_race_006"
```

An entry near $`1`$ says the posterior strongly orders that row’s firm
above the column’s; an entry near $`0.5`$ says the data cannot separate
them. The matrix encodes pairwise posterior order only and carries no
superiority statement.

**`get_control(fit)`** — the exact
[`gp_control()`](https://joonho112.github.io/gradepath/reference/gp_control.md)
that produced the fit, your reproducibility receipt: the grid, the
backend, the precision rule, the interval level, the solver options, the
seed, the schema version, and the provenance:

``` r

ctrl <- get_control(fit)
names(ctrl)
#> [1] "lambda_grid"    "backend"        "precision_rule" "interval_level"
#> [5] "solver_options" "seed"           "schema_version" "provenance"
c(backend = ctrl$backend, precision_rule = ctrl$precision_rule)
#>        backend precision_rule 
#>       "gurobi"      "krw_gmm"
```

## 6. `coef()`: a firm → grade lookup

[`coef()`](https://rdrr.io/r/stats/coef.html) is the base generic an R
user reaches for first, so a `gp_fit` answers it with the grading:
`coef(fit)` returns the same named integer vector as `get_grades(fit)` —
names = firm ids, values = the assigned grade — surfaced under the name
your muscle memory already knows:

``` r

identical(coef(fit), get_grades(fit))
#> [1] TRUE
head(coef(fit))                 # e.g. krw_race_001 = 2
#> krw_race_001 krw_race_002 krw_race_003 krw_race_004 krw_race_005 krw_race_006 
#>            2            2            2            2            2            3
coef(fit)["krw_race_001"]       # the grade of one firm, by id
#> krw_race_001 
#>            2
```

Because it is named, it doubles as a quick firm-to-grade lookup table:
index it by id to read any one firm’s grade, or
`data.frame(id = names(coef(fit)), grade = coef(fit))` to spread the
whole vector into two columns for a join. When you simply need “which
firm got which grade” without the posterior summaries, this is the
lightest read in the package; a grade label carries no
ranking-superiority statement.

## 7. The CLI renderers: embed a card in a log

[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) write to the console,
but sometimes you want the same result-first card as a *character
vector* you can splice into a build log, a plain-text report, or an HTML
`<pre>` block. The `format_*_cli` renderers are exactly that: the
formatter behind the console output, returned as one string per line.

`format_gp_fit_cli(fit)` returns the box of section 2 as lines you
[`cat()`](https://rdrr.io/r/base/cat.html) or
[`writeLines()`](https://rdrr.io/r/base/writeLines.html):

``` r

writeLines(format_gp_fit_cli(fit))
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

[`format_gp_report_card_cli()`](https://joonho112.github.io/gradepath/reference/format_gp_report_card_cli.md)
does the same for the ranked card: it leads with the unit count, grade
composition, and selected $`\lambda`$, shows a HEAD of the ranked table,
and closes with the welfare footnote:

``` r

writeLines(format_gp_report_card_cli(get_report_card(fit)))
#> <gp_report_card>  unit: unit | 97 rows | grades: 3 (2/81/14) | selected lambda: 0.25
#> sorted by Condorcet rank (endpoint lambda = 1; secondary key: id)
#> 
#>   sort_rank  id            label         grade  posterior_mean  CI (90%)        estimate
#>    1         krw_race_062  krw_race_062  1      0.250           [0.117, 0.342]  0.330   
#>    2         krw_race_015  krw_race_015  1      0.233           [0.081, 0.447]  0.433   
#>    3         krw_race_017  krw_race_017  2      0.191           [0.032, 0.383]  0.377   
#>    4         krw_race_037  krw_race_037  2      0.188           [0.040, 0.395]  0.402   
#>    5         krw_race_095  krw_race_095  2      0.187           [0.031, 0.373]  0.350   
#>    6         krw_race_010  krw_race_010  2      0.178           [0.032, 0.352]  0.343   
#>    7         krw_race_085  krw_race_085  2      0.170           [0.023, 0.342]  0.169   
#>   ... 89 more rows ...
#>   97         krw_race_006  krw_race_006  3      0.017           [0.001, 0.044]  -0.032  
#> 
#> i  grade 1 = most-extreme theta (firms: most discriminatory; names: best-treated). See ?gp_report_card and Appendix A.
#> i  full table: as.data.frame(card)   formatted CLI: format_gp_report_card_cli(card)
```

Because each returns plain ASCII as a character vector, both compose
with everything that takes text — dropping a result-first card into a
pipeline log is a one-liner,
`message(paste(format_gp_fit_cli(fit), collapse = "\n"))`, so a teammate
reading the log sees the grading without opening R. The renderer’s
welfare note describes grade 1 as the most-extreme $`\theta`$ and is not
hard-coded to one reading: it adapts between the firms application (most
discriminatory) and the names application (best-treated), so the same
code is honest for either.

## 8. Exporting

Three exports cover almost every hand-off.

**CSV.** The first is the analytical artifact: write the tidy table to
CSV and it lands in a spreadsheet, a paper table, or a co-author’s
pipeline. This is the command you would run on your own machine (we do
not write into the package working directory here):

``` r

write.csv(as.data.frame(fit), "report_card_race.csv", row.names = FALSE)
```

**RDS round-trip.** The second is the reproducible artifact: the whole
`gp_fit`, saved with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html),
round-trips **exactly** — every slot, the prior, the posterior, the
pairwise matrix, and the control come back identical. We demonstrate it
live through a [`tempfile()`](https://rdrr.io/r/base/tempfile.html), so
the chunk both runs and proves the integrity:

``` r

path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
identical(fit, readRDS(path))   # a gp_fit round-trips through RDS exactly
#> [1] TRUE
```

A `TRUE` here is the guarantee you rely on when you ship a solved fit
alongside a paper: a collaborator who
[`readRDS()`](https://rdrr.io/r/base/readRDS.html)s the file recovers
the identical object you solved, with no re-solve and no drift. (This is
precisely how the bundled parity fits travel inside the package.)

**ASCII card to a file.** The third is the human-readable artifact: the
result-first card written to a text file, ready to attach to an issue or
paste into a README. Like the CSV, this is the command you would run
yourself:

``` r

writeLines(format_gp_report_card_cli(get_report_card(fit)), "card.txt")
```

## 9. The same reads for gender

Every read and export above is demographic-agnostic — the accessors,
[`coef()`](https://rdrr.io/r/stats/coef.html), the renderers, and the
export idioms do not know or care which gap they describe. The bundled
gender fit answers every call of sections 2–8 identically, only with its
own numbers:

``` r

gfit <- gp_parity_fit("gender")
table(get_grades(gfit))                     # the gender composition: 1 / 3 / 89 / 4
#> 
#>  1  2  3  4 
#>  1  3 89  4
round(get_prior(gfit)$metadata$beta, 4)     # the gender precision exponent (1.2554)
#> [1] 1.2554
head(coef(gfit))
#> krw_gender_001 krw_gender_002 krw_gender_003 krw_gender_004 krw_gender_005 
#>              3              3              3              3              3 
#> krw_gender_006 
#>              3
```

The gender fit sorts the 97 firms into four grades (`1 / 3 / 89 / 4`),
and its precision exponent is $`\beta = 1.2554`$ — a different number
from the race fit’s $`0.5095`$, as it should be, since the two
demographics carry different gap series. Its
[`summary()`](https://rdrr.io/r/base/summary.html),
`as.data.frame(gfit)`, the six `get_*` accessors, the `format_*_cli`
renderers, and the CSV / RDS export all behave exactly as above:
`as.data.frame(gfit)` is the same ten columns, `coef(gfit)` is the same
named lookup, and `saveRDS(gfit, ...)` round-trips with the same
[`identical()`](https://rdrr.io/r/base/identical.html) guarantee. The
welfare reading of grade 1 follows the application and is never
hard-coded to one interpretation.

## 10. Where to next

- **Getting started** —
  [A1](https://joonho112.github.io/gradepath/articles/a1-getting-started.md)
  is the gentle first contact, if you arrived here without a fit in
  hand;
  [A2](https://joonho112.github.io/gradepath/articles/a2-the-grading-workflow.md)
  walks the seven stages that fill the slots this vignette reads.
- **Figures** —
  [A4](https://joonho112.github.io/gradepath/articles/a4-figures.md) is
  the cookbook for all four KRW signature figures (frontier, posterior
  contrast, report card, discordance), built from these same accessors.
- **Solvers and calibration** —
  [A5](https://joonho112.github.io/gradepath/articles/a5-solvers-and-calibration.md)
  covers the open-source backends and the calibration harness.
- **Foundations** —
  [`vignette("m1-foundations-and-notation")`](https://joonho112.github.io/gradepath/articles/m1-foundations-and-notation.md)
  formalizes the estimate → rank → grade pipeline, and the method track
  ([`vignette("m2-precision-and-standardization")`](https://joonho112.github.io/gradepath/articles/m2-precision-and-standardization.md),
  [`vignette("m3-deconvolution-and-posterior")`](https://joonho112.github.io/gradepath/articles/m3-deconvolution-and-posterior.md),
  [`vignette("m4-grading-frontier-and-report-cards")`](https://joonho112.github.io/gradepath/articles/m4-grading-frontier-and-report-cards.md),
  [`vignette("m5-two-level-and-calibration")`](https://joonho112.github.io/gradepath/articles/m5-two-level-and-calibration.md))
  derives each stage whose output you read here.

A closing reminder on reading the output honestly: a grade is a sorting
of firms by posterior evidence, never a contest, and grade 1 is the
most-extreme grade rather than a judgment. The one-level race and gender
core is accepted within its declared scope, and the numbers here match
the paper ($`\beta = 0.5095`$ versus the reported `0.510`;
$`\lambda = 0.25`$ is the KRW baseline).

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
#> [17] lattice_0.22-9     R6_2.6.1           generics_0.1.4     knitr_1.51        
#> [21] htmlwidgets_1.6.4  gurobi_13.0-1      tibble_3.3.1       desc_1.4.3        
#> [25] bslib_0.11.0       pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.2.0       
#> [29] cachem_1.1.0       xfun_0.57          fs_2.1.0           sass_0.4.10       
#> [33] S7_0.2.2           otel_0.2.0         cli_3.6.6          withr_3.0.2       
#> [37] pkgdown_2.2.0      magrittr_2.0.5     digest_0.6.39      grid_4.6.0        
#> [41] lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5     glue_1.8.1        
#> [45] farver_2.1.2       ragg_1.5.2         rmarkdown_2.31     tools_4.6.0       
#> [49] pkgconfig_2.0.3    htmltools_0.5.9
```

## References

Kline, Patrick, Evan K. Rose, and Christopher R. Walters. 2024. “A
Discrimination Report Card.” *American Economic Review* 114 (8):
2472–525. <https://doi.org/10.1257/aer.20230700>.
