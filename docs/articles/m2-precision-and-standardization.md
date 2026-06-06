# M2: Precision and standardization

Abstract

This is the beta-GMM precision core of the KRW pipeline, made rigorous.
Heteroskedastic stage-1 estimates make the noisiest units look the most
extreme; KRW standardize this away with the precision-dependence model
theta = s^beta v, estimating the precision exponent beta by two-step
optimally-weighted GMM rather than assuming it. This vignette derives
the moment conditions and the GMM, explains the multiplicative (race)
versus additive (gender) split, and verifies beta = 0.5095 (race) and
1.2554 (gender) live against the bundled parity fits, with the GMM
overidentification J-test. The reader leaves with the standardized
r-scale the deconvolution of m3 operates on.

This is the second method-track vignette, and the first to derive a
stage. It takes the one paragraph of
[m1](https://joonho112.github.io/gradepath/articles/m1-foundations-and-notation.md)’s
§3 — the line that said “the precision dependence is standardized away,
the derivation is m2” — and makes it rigorous. We summarize the theory
from Kline, Rose & Walters (Kline et al. 2024) (their §5.1, Equations
8–10) rather than re-derive every algebra step, and we follow every
displayed equation with a verification panel that confirms it live
against **both** bundled parity fits.

We load both fits once, at the top. Every panel below is a fast read off
an already-fitted `precision_fit` slot — there is no GMM solve anywhere
in this document. The fits are pre-solved, so the beta-GMM stage we are
deriving has already run; here it is the ground truth we check the
precision core against.

``` r

fit_race   <- gp_parity_fit("race")     # 97-firm RACE card, grades 2 / 81 / 14
fit_gender <- gp_parity_fit("gender")   # 97-name GENDER card, grades 1 / 3 / 89 / 4
```

These are the 97-firm RACE and gender parity fits at the KRW baseline
`lambda = 0.25`. The applied track
([a2](https://joonho112.github.io/gradepath/articles/a2-the-grading-workflow.md))
shows how they are produced. The precision fit lives in one slot,
`fit$precision_fit`, with the same layout on both demographics — the
crosswalk we verify against throughout:

``` r

names(fit_race$precision_fit)
#> [1] "parameters"     "moments"        "diagnostics"    "scale"         
#> [5] "schema_version" "provenance"     "warnings"
```

`parameters` carries the estimated model (`model_form`, `beta`, `mu`,
`characteristic`); `moments` carries the GMM moment vector `m_hat`, its
covariance `V_m`, and the overidentification test (`J`, `df`,
`p_value`); and `scale` records that the prior lives on the standardized
r-scale.

## 1. Why standardize

Start with the problem, because the whole stage exists to solve it. You
hold $`n`$ units — 97 large U.S. firms — each with a noisy stage-1
estimate $`\hat\theta_i`$ of a latent effect $`\theta_i`$ and a standard
error $`s_i`$ that says how precisely $`\hat\theta_i`$ was measured. The
estimates are **heteroskedastic**: some firms were measured with far
more correspondence applications than others, so their $`s_i`$ varies by
an order of magnitude across the panel.

Heteroskedasticity has a sharp consequence for any procedure that
compares units by how extreme they look. A firm measured with few
applications has a large $`s_i`$, and a large $`s_i`$ means
$`\hat\theta_i`$ can land far from $`\theta_i`$ by chance alone. So the
noisiest units are *exactly* the ones most likely to produce an
eye-catching $`\hat\theta_i`$ — not because their latent effect is
extreme, but because their estimate is loose. A grading that took
$`\hat\theta_i`$ at face value would hand the most-extreme grades to the
least-precisely-measured firms, which is precisely backwards.

> “Estimates of firm-level discrimination are inherently noisy … A
> report card summarizes this evidence by sorting firms into a small
> number of grades.” — Kline, Rose & Walters (Kline et al. 2024)

There is a second, subtler reason, and it is the one KRW build the model
around. Empirical Bayes shrinkage (Walters 2024) borrows strength across
units by shrinking each estimate toward a common prior, with the amount
of shrinkage set by each unit’s noise. That machinery assumes the latent
effect $`\theta_i`$ is drawn *independently* of its standard error
$`s_i`$. In the firm setting that assumption is implausible: firms with
many postings generate more applications (small $`s_i`$) and may also
differ systematically in their discrimination. If $`\theta_i`$ and
$`s_i`$ are correlated, a single prior fit to the whole sample
over-shrinks some units and under-shrinks others, and the bias is
invisible — it hides inside the very precision the method relies on.

KRW’s fix is to put every unit on a **common precision scale** before
any comparison or shrinkage happens. The deconvolution of
[m3](https://joonho112.github.io/gradepath/articles/m3-deconvolution-and-posterior.md),
the pairwise posteriors, and the grading frontier of
[m4](https://joonho112.github.io/gradepath/articles/m4-grading-frontier-and-report-cards.md)
all operate on that scale — never on the raw $`\hat\theta_i`$. The rest
of this vignette derives the scale and checks the one number that
defines it.

## 2. The precision-dependence model and the r-scale

The model is a power law linking the latent effect to the precision of
its measurement (Kline et al. 2024, Eq. 8):

``` math
\theta_i \;=\; s_i^{\beta}\, v_i,
\qquad
v_i \mid s_i \;\overset{\text{iid}}{\sim}\; G_v .
```

Read it as a decomposition. The factor $`s_i^{\beta}`$ carries
everything about $`\theta_i`$ that travels with precision; the residual
$`v_i = \theta_i/s_i^{\beta}`$ is what is left once that dependence is
divided out, and *by construction* $`v_i`$ is conditionally independent
of $`s_i`$. The exponent $`\beta`$ is the one number that controls how
much the latent effect rides on precision:

- $`\beta = 0`$ gives $`\theta_i = v_i`$ — the homoskedastic textbook
  case, no precision dependence, latent effect independent of standard
  error.
- $`\beta \approx 0.5`$ is the square-root-of-precision regime: the
  latent scale grows like the standard deviation of the noise, halfway
  between “no dependence” and “fully proportional.”
- $`\beta = 1`$ makes the latent scale move one-for-one with $`s_i`$.

The crucial methodological point — the reason this stage is GMM and not
a guess — is that $`\beta`$ is **estimated, not assumed**. Picking
$`\beta = 0`$ (assume away the problem) or $`\beta = 1`$ (assume it
fully) would each bias every downstream residual, prior, posterior, and
grade. The data are allowed to say where between those poles the truth
sits.

Once $`\beta`$ is in hand, the object the rest of the pipeline lives on
is the **standardized r-scale**: the residual you get by dividing the
*estimate* by its own $`s_i^{\beta}`$,

``` math
r_i \;=\; \hat\theta_i \,/\, s_i^{\beta} .
```

This is the scale on which the prior $`G`$ is deconvolved (m3) and the
comparison weights are built. It is also where the noise rescales: under
the model the estimate is $`\hat\theta_i = \theta_i + \varepsilon_i`$
with $`\varepsilon_i \sim N(0, s_i^2)`$, so dividing by $`s_i^{\beta}`$
leaves the residual with noise scale $`s_{v,i} = s_i^{1-\beta}`$. That
$`s_{v,i} = s_i^{1-\beta}`$ is the *likelihood* scale m3 will use; keep
it in view, because the moments in §3 are built from it.

**Verification panel.** Read $`\beta`$ off the race fit, then build the
r-scale by hand from $`\hat\theta`$ and $`s`$. Note the standard-error
column is `est$s` (not `est$se`), and there is no stored
`precision_fit$r` slot — the r-scale is *computed*, never stored:

``` r

beta_race <- fit_race$precision_fit$parameters$beta
round(beta_race, 4)                           # the race precision exponent
#> [1] 0.5095

est <- fit_race$estimates
r   <- est$theta_hat / est$s^beta_race        # the standardized r-scale
head(r, 3)
#> [1] 0.8330686 0.3116278 0.4740612
```

gradepath lands on $`\beta = 0.5095`$ for race, a faithful match to the
$`\sim 0.510`$ KRW report. The same value is reachable through the prior
accessor, `get_prior(fit_race)$metadata$beta`, because the prior is the
object that lives on the r-scale. The first three r-scale values are
0.8331, 0.3116, 0.4741.

## 3. The moment conditions

Why is $`\beta`$ identifiable at all? Because the model makes a
*testable* promise: after standardizing, the residual should look the
same regardless of how precisely it was measured. KRW turn that promise
into moment conditions (Kline et al. 2024, Eq. 9).

Form the **studentized residual** — the r-scale point estimate
recentered and scaled to be mean-zero, unit-variance under the model.
With $`\hat v_i = \hat\theta_i/s_i^{\beta}`$, its noise scale
$`s_{v,i} = s_i^{1-\beta}`$, and the implied latent SD $`\sigma`$,

``` math
r_i \;=\; \frac{\hat v_i - \mu}{\sqrt{\sigma^2 + s_{v,i}^2}}
\quad\text{(race)},
\qquad
r_i \;=\; \frac{\hat v_i}{\sqrt{\sigma^2 + s_{v,i}^2}}
\quad\text{(gender, with }\hat v_i = (\hat\theta_i - \mu)/s_i^{\beta}\text{)} ,
```

so that under correct specification $`r_i`$ is approximately i.i.d. with
mean zero and unit variance — *regardless of $`s_i`$*. (The two forms —
subtracting $`\mu`$ in the numerator versus inside $`\hat v_i`$ — are
the multiplicative-versus-additive split of §6; the four moments below
are identical for both.) The four base moment conditions ask exactly the
mean-zero, unit-variance, $`s`$-free property, in two pairs:

``` math
\underbrace{E[r_i] = 0}_{\text{mean is zero}},
\quad
\underbrace{E[r_i\, s_i] = 0}_{\text{...even conditional on }s_i},
\quad
\underbrace{E[r_i^2 - 1] = 0}_{\text{variance is one}},
\quad
\underbrace{E[(r_i^2 - 1)\, s_i] = 0}_{\text{...even conditional on }s_i},
```

stacked into the sample vector
$`g_i = [\,r_i,\; r_i s_i,\; r_i^2 - 1,\; (r_i^2 - 1) s_i\,]`$.

The structure *is* the identification story, and it is worth slowing
down on. The two **level** moments ($`r`$ and $`r^2-1`$) pin the
location $`\mu`$ and scale $`\sigma`$ of the standardized residual. The
two **interactions with $`s_i`$** ($`r\,s`$ and $`(r^2-1)\,s`$) are the
ones that earn their keep: they test whether the residual’s mean or
variance still trends with precision *after* standardizing. If $`\beta`$
were wrong, some precision dependence would survive into the residual,
those interaction moments would not vanish, and the GMM would move
$`\beta`$ to make them vanish. That residual-versus-$`s`$ dependence is
precisely what identifies $`\beta`$*separately* from $`\sigma`$ — which
is why this is GMM and not a conditional-mean least-squares fit. A mean
regression has no analogue of the variance moments and so cannot recover
$`\beta`$ this way.

These four moments identify the one-level parameter vector
$`\delta = [\mu,\ \log\sigma_\xi,\ \beta]`$ — **three parameters, four
moments**, so **one over-identifying restriction** is left over and
becomes the specification test of §5. (The intercept enters differently
under the additive gender form, but the count still leaves the same
single over-identifying degree of freedom — which is what the J-test
reads for both fits in §5.) The moment **covariance** is estimated
cluster-robust on industry — one-level it is the ordinary outer-product
average $`\hat\Omega = \tfrac1N\sum_i g_i g_i'`$; with industry
structure it sums cross-products only within shared industries — so the
standard errors respect the within-industry coupling. The two-level
extension that appends two between-industry moments is
[m5](https://joonho112.github.io/gradepath/articles/m5-two-level-and-calibration.md);
here we use the four base moments. This stored moment vector and its
covariance are not an abstraction — they are slots on the fit, which is
what lets us verify the GMM rather than assert it.

## 4. Two-step optimally-weighted GMM

GMM chooses $`\delta = [\mu,\ \log\sigma_\xi,\ \beta]`$ to make the
sample moment average $`g(\delta)`$ as close to zero as a weight matrix
$`\mathsf W`$ allows, minimizing the quadratic form

``` math
J(\delta;\,\mathsf W) \;=\; N\, g(\delta)'\,\mathsf W\, g(\delta) .
```

The choice of $`\mathsf W`$ is the whole game with an over-identified
system, because with four moments and three parameters you cannot drive
all four to zero at once — you must decide which to favor. The
**efficient** answer, standard GMM theory, is to weight each moment by
the inverse of its variance: set $`\mathsf W = \hat\Omega^{-1}`$, so
noisy moments count for less. But $`\hat\Omega`$ itself depends on
$`\delta`$, which is the chicken-and-egg the **two-step** procedure
resolves:

1.  **First step ($`\mathsf W = I`$).** Minimize $`J(\delta; I)`$ with
    the identity weight to get a consistent first-pass
    $`\hat\delta^{(1)}`$. KRW seed the search at their published
    estimates (race one-level $`\approx [-1.18,\ -1.58,\ 0.51]`$) and
    otherwise multistart.
2.  **Reweight.** Evaluate
    $`\hat\Omega_1 = \hat\Omega(\hat\delta^{(1)})`$ and set the optimal
    weight to its inverse, $`\mathsf W = \hat\Omega_1^{-1}`$
    (implemented as the Cholesky factor of $`\hat\Omega_1^{-1}`$, so the
    objective is the squared norm $`\lVert \mathsf W^{1/2} g\rVert^2`$
    the optimizer minimizes directly).
3.  **Second step.** Re-minimize $`J(\delta; \mathsf W)`$ from
    $`\hat\delta^{(1)}`$ to get the reported $`\hat\delta`$. Standard
    errors come from the optimal-GMM sandwich
    $`C = (G'\hat\Omega_2^{-1} G)^{-1}/N`$ with
    $`G = \partial g/\partial\delta'`$; because $`\hat\Omega`$ was
    clustered on industry, the SEs are cluster-robust.

The intuition for “optimal weight = inverse covariance” is the same one
behind weighted least squares: a moment you can measure tightly deserves
a tight constraint, and a moment drowning in sampling noise should not
be allowed to drag the estimate around. Inverse-variance weighting is
just that instinct made precise — and it is also what makes the leftover
quadratic at the optimum into a clean $`\chi^2`$ statistic (§5).

From this *single* fit the core emits two carriers the deconvolution of
m3 consumes: the matched **moment vector** $`\hat m`$ (the GMM-implied
prior location and scale) and its delta-method **covariance** $`V_m`$.
These are what m3 scores its penalty grid and builds its support caps
against, which is why the precision and deconvolution steps are one
unified core rather than two independent fits.

**Verification panel.** Read $`\hat m`$ and $`V_m`$ off the race fit —
what GMM minimized is encoded right here:

``` r

mom <- fit_race$precision_fit$moments
mom$m_hat                                      # GMM-implied [mu, sigma_xi]
#> [1] 0.3073472 0.2066007
round(mom$V_m, 6)                              # their (delta-method) covariance
#>          [,1]     [,2]
#> [1,] 0.021440 0.013937
#> [2,] 0.013937 0.011146
isSymmetric(mom$V_m)                           # symmetric: TRUE
#> [1] TRUE
eigen(mom$V_m, symmetric = TRUE, only.values = TRUE)$values  # both > 0 => positive-definite
#> [1] 0.031149336 0.001435902
```

The implied prior mean is $`\hat m_1 = 0.3073`$ — exactly the $`\mu`$ in
`parameters` (the race intercept) — and the implied scale is
$`\hat m_2 = 0.2066`$. The covariance $`V_m`$ is symmetric (so
[`isSymmetric()`](https://rdrr.io/r/base/isSymmetric.html) returns
`TRUE`) with both eigenvalues positive — symmetric positive-definite,
exactly as a delta-method covariance must be. Together they are the
coupling the next stage’s penalty grid and support caps are built from;
here they confirm the GMM emits both from a single fit.

## 5. The overidentification J-test

With four moments and three parameters, one restriction is
over-identifying — and that surplus is a free specification test, no
extra data required. At the optimum, the weighted quadratic does not
quite reach zero (it cannot, with a surplus moment), and **how far it
falls short** is the test statistic. Hansen’s J-statistic is the
second-step objective at the optimum,

``` math
J \;=\; N\, g(\hat\delta)'\,\hat\Omega^{-1} g(\hat\delta)
\;\overset{a}{\sim}\; \chi^2(\#\text{moments} - \#\text{params}),
```

with one degree of freedom here ($`4 - 3 = 1`$; the $`\chi^2(1)`$
critical value at the 5% level is $`3.84`$). Read it plainly: if the
precision model is right, the standardized residual really is mean-zero,
unit-variance, and free of $`s`$-dependence, so all four moments are
near zero, $`J`$ is small, and the p-value is large. A large $`J`$
(small p-value) would say the moments *cannot* be jointly satisfied —
the $`\theta = s^{\beta} v`$ form is too rigid for the data. A $`J`$
that grew with the sample size would be the signature of a misspecified
model, not sampling noise.

**Verification panel.** Read the J-statistic, its degrees of freedom,
and its p-value off both fits:

``` r

jt <- function(fit) {
  m <- fit$precision_fit$moments
  c(J = round(m$J, 4), df = m$df, p = round(m$p_value, 4))
}
rbind(race = jt(fit_race), gender = jt(fit_gender))
#>             J df      p
#> race   0.1021  1 0.7493
#> gender 0.0112  1 0.9158
```

For race, $`J = 0.1021`$ on 1 degree of freedom, p-value 0.7493; for
gender, $`J = 0.0112`$, p-value 0.9158. Both p-values are far above any
conventional threshold, so **the moment conditions are not rejected**:
the precision-dependence model is consistent with the data for both
demographics. That is the honest reading — the J-test does not establish
the model, it fails to find evidence against it, which is exactly the
assurance you want before deconvolving on the r-scale.

## 6. The multiplicative-versus-additive split

The model in §2 was written multiplicatively:
$`\theta_i = s_i^{\beta} v_i`$, with the residual formed as
$`\hat v_i = \hat\theta_i / s_i^{\beta}`$. That is the race form. Gender
uses an **additive** variant that subtracts a location *before*
dividing:

``` math
\text{multiplicative (race):}\quad
\hat v_i = \frac{\hat\theta_i}{s_i^{\beta}},
\qquad
\text{additive (gender):}\quad
\hat v_i = \frac{\hat\theta_i - \mu}{s_i^{\beta}} .
```

Why two forms? Because the location of the latent effect means different
things for the two demographics, and the model’s intercept $`\mu`$ has
to sit where the data centers. For **race**, the average firm
discriminates against minority applicants, so the latent location is
bounded away from zero and *positive* — KRW parameterize it
multiplicatively (internally $`\mu = e^{\delta_1}`$, which keeps it
positive), and the standardization scales around that positive center.
For **gender**, the average gap can fall on either side of zero, so the
location is modeled as a free additive intercept $`\mu`$ that is
subtracted off before standardizing. The split is not a tuning knob the
analyst turns: gradepath **selects the form per demographic** to match
the centering KRW use, and the choice propagates through the moment
conditions and the variance formulas.

**Verification panel.** Read the model form, the precision exponent, and
the intercept for both fits side by side:

``` r

split_row <- function(fit) {
  p <- fit$precision_fit$parameters
  c(model_form = p$model_form,
    beta = round(p$beta, 4),
    mu   = round(p$mu, 4))
}
rbind(race = split_row(fit_race), gender = split_row(fit_gender))
#>        model_form       beta     mu       
#> race   "multiplicative" "0.5095" "0.3073" 
#> gender "additive"       "1.2554" "-0.0088"
```

Race comes back **multiplicative** with $`\beta = 0.5095`$ and a
positive intercept $`\mu = 0.3073`$; gender comes back **additive** with
$`\beta = 1.2554`$ and a near-zero intercept $`\mu = -0.0088`$. The
positive race intercept and the near-zero gender intercept are exactly
the centering story above, read straight off the fits. Note too that the
two $`\beta`$ values differ substantially ($`\approx 0.51`$ vs
$`\approx 1.26`$): precision dependence is genuinely stronger in the
gender panel, which is another reason $`\beta`$ must be estimated
separately rather than fixed at a shared value.

## 7. What the next stage receives

The precision core hands one scale and two carriers forward. The
**scale** is the r-scale: the prior $`G`$ and the comparison weights are
deconvolved and built on $`r_i = \hat\theta_i / s_i^{\beta}`$, with the
per-unit likelihood noise $`s_{v,i} = s_i^{1-\beta}`$. The **carriers**
are the matched moment vector $`\hat m`$ and its covariance $`V_m`$ from
§4, which m3’s deconvolution uses to score its penalty grid and set its
support caps. Everything downstream — the deconvolved prior, the
posteriors, the pairwise outranking matrix, the grades — lives on this
r-scale and never sees the raw $`\hat\theta_i`$ again.

**Verification panel.** Confirm the recorded scale, and rebuild the
gender r-scale in its additive form — `(theta_hat - mu)/s^beta`, using
the gender fit’s own $`\beta`$ and intercept $`\mu`$ (the additive
variant of §6, not the multiplicative race transform):

``` r

fit_race$precision_fit$scale                  # the scale the prior lives on: "r"
#> [1] "r"

est_g  <- fit_gender$estimates
beta_g <- fit_gender$precision_fit$parameters$beta
mu_g   <- fit_gender$precision_fit$parameters$mu
head((est_g$theta_hat - mu_g) / est_g$s^beta_g, 3)   # gender r-scale: additive form, (theta_hat - mu)/s^beta
#> [1] -1.752534  3.427092 -4.098203
```

The scale slot reads “r”, and the gender r-scale is computed by the
additive transform of §6 — subtracting the gender intercept $`\mu`$
before dividing by $`s^\beta`$ — using the gender fit’s own $`\beta`$
and $`\mu`$. Because the gender intercept is near zero, this lands close
to the raw quotient, but the form is the additive one the deconvolution
actually receives. That single column is what
[m3](https://joonho112.github.io/gradepath/articles/m3-deconvolution-and-posterior.md)
deconvolves a prior over.

## 8. Where to next

With $`\beta`$ estimated, the r-scale built, and the precision model
checked against both fits and its own J-test, the pipeline continues:

- **Foundations and notation** —
  [m1](https://joonho112.github.io/gradepath/articles/m1-foundations-and-notation.md)
  is the symbol table and the two frozen invariants this vignette
  assumed; its §3 is the one-paragraph version of everything above.
- **Deconvolution and posterior** —
  [`vignette("m3-deconvolution-and-posterior")`](https://joonho112.github.io/gradepath/articles/m3-deconvolution-and-posterior.md)
  takes the r-scale and the carriers $`\hat m`$, $`V_m`$ from here,
  deconvolves the prior $`G`$ over $`v`$, and forms the posteriors and
  the pairwise matrix $`\pi_{ij}`$.
- **Grading frontier and report cards** —
  [`vignette("m4-grading-frontier-and-report-cards")`](https://joonho112.github.io/gradepath/articles/m4-grading-frontier-and-report-cards.md)
  turns the posteriors into grades — integer labels in
  $`\{1, \dots, n\}`$, never letter grades — via the grade integer
  program and the $`\lambda`$ frontier.
- **Two-level and calibration** —
  [`vignette("m5-two-level-and-calibration")`](https://joonho112.github.io/gradepath/articles/m5-two-level-and-calibration.md)
  adds the two between-industry moments that extend §3 to the six-moment
  industry model, and validates the whole pipeline by Monte Carlo.

For the workflow side,
[a1](https://joonho112.github.io/gradepath/articles/a1-getting-started.md)
is the first contact,
[a2](https://joonho112.github.io/gradepath/articles/a2-the-grading-workflow.md)
runs this same precision stage hands-on with race and gender side by
side, and
[a5](https://joonho112.github.io/gradepath/articles/a5-solvers-and-calibration.md)
covers the calibration harness that the J-test logic above feeds into.
Throughout, gradepath’s printed output is result-first and
quantity-anchored (the package decision GP-DEC-14-A): it names the
quantity rather than a verdict, so the same code path serves both
demographics without the printed text inverting on one of them.

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

Walters, Christopher R. 2024. “Empirical Bayes Methods in Labor
Economics.” In *Handbook of Labor Economics*, vol. 5. Elsevier.
<https://doi.org/10.1016/bs.heslab.2024.11.001>.
