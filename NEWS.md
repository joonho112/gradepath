# gradepath 0.5.0

This is the first public pre-release of **gradepath**, a native-R replication of the discrimination-report-card analysis of Kline, Rose & Walters (2024, AER). The one-level race and gender core (M1) is accepted under a predeclared scoped policy, the two-level industry surface (M2) is partially accepted, and the release ships a full documentation layer: a reader-first reference, a two-track vignette set, a pkgdown site, a hex logo, and this README/NEWS.

## Major changes

- **Native KRW estimation core, entirely in R.** The full pipeline — beta-GMM precision standardization; empirical-Bayes prior deconvolution and posterior shrinkage; the pairwise posterior outranking matrix; the Bayes-risk grade integer program; and the information-reliability frontier — now runs in R, replacing the original Stata + R + Python + MATLAB + Gurobi pipeline.
- **License-free to install, load, and test.** The grade integer program solves with **Gurobi by default**; open backends **HiGHS** (direct) and **ROI** plugins (GLPK / SYMPHONY) are supported, so the workflow runs without a commercial license. The package installs, loads, and tests without any solver present.
- **M1 (one-level race + gender) is accepted under a predeclared scoped policy.** It recovers the published grade distributions — race **2 / 81 / 14** and gender **1 / 3 / 89 / 4** — as proven-optimal solves (`mip_gap = 0`) at the KRW baseline `lambda = 0.25`, with `beta = 0.5095` (race) and `beta = 1.2554` (gender) matching the paper.
- **M2 (two-level / industry) is partially accepted.** Exact industry grade counts are reproduced; some continuous rows rest on fixture-parity evidence rather than direct reproduction (notably the race two-level L02 row is an approximation).
- **ebrecipe seam.** Low-level empirical-Bayes primitives are supplied by the **ebrecipe** package behind a controlled seam — an input container plus a one-level cross-check.

## New features

- **The four signature KRW figures as composable ggplot2 verbs:** `gp_plot_frontier()`, `gp_plot_posterior_contrast()`, `gp_plot_report_card()`, and `gp_plot_discordance()`, together with `theme_gradepath()`, gradepath colour/fill scales, and `autoplot` methods.
- **Result-first console output.** `print` and `summary` methods for every gradepath object, each with an interpretation block and a glossary. Grades are integer labels {1..n}, a sorting by posterior evidence.
- **A seeded Monte-Carlo calibration harness,** `gp_calibrate()` — a method check on the estimation pipeline (not a reproduction of the paper's Table E1).
- **43 exported functions** spanning the estimation core, the grade integer program and its solver backends, the figure verbs, and the print/summary layer.

## Documentation

- **Reader-first reference on all exports:** roxygen documentation organized by an 11-family `@family` taxonomy, with runnable examples.
- **A pkgdown site:** <https://joonho112.github.io/gradepath/>.
- **A two-track vignette set:** applied **a1–a5** for users who want to run the analysis, and method **m1–m5** for readers who want the estimation theory.
- **A hex logo and a README** to orient new users.

## Known limitations

- **Pre-release status (v0.5.0).** The M2 two-level surface is only partially accepted: the race two-level **L02** continuous row is an approximation that rests on fixture-parity evidence rather than direct reproduction.
- **`gp_krw_gmm_input()` is not yet exported** in this release. To assemble the bundled per-firm beta-GMM input today, read the shipped CSV directly, e.g.

    ```r
    read.csv(
      system.file("extdata/krw-gmm-input/theta_estimates_matlab_race.csv",
                  package = "gradepath"),
      header = FALSE
    )
    ```

- **Solver performance.** **Gurobi** is the recommended default backend (a real 97-firm grade solve takes about 3 minutes); the license-free **HiGHS** backend works but is several times slower on the same solve.
- **`gp_calibrate()` is a calibration method check,** not a reproduction of KRW's Table E1 (there is no parseable artifact, and the MATLAB RNG is not reproducible).
