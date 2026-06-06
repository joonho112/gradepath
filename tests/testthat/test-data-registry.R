# Bundled data, registry, and archive fixtures.
# krw_firms and gp_registry are lazy-loaded package data; fixtures ship under
# inst/extdata/fixtures/ (headerless CSVs — read with header = FALSE).

# ---- krw_firms -------------------------------------------------------------
test_that("krw_firms has the expected shape and id structure", {
  expect_s3_class(krw_firms, "data.frame")
  expect_identical(nrow(krw_firms), 97L)
  expect_true(all(c("firm_id", "theta_hat_race", "se_race",
                    "theta_hat_gender", "se_gender",
                    "industry", "label") %in% names(krw_firms)))
  # firm_id is the gappy KRW numbering, NOT 1:97
  expect_false(all(krw_firms$firm_id == seq_len(nrow(krw_firms))))
  expect_identical(length(unique(krw_firms$firm_id)), 97L)
})

test_that("krw_firms estimate/SE columns are complete and SEs positive", {
  expect_false(anyNA(krw_firms$theta_hat_race))
  expect_false(anyNA(krw_firms$se_race))
  expect_false(anyNA(krw_firms$theta_hat_gender))
  expect_false(anyNA(krw_firms$se_gender))
  expect_true(all(krw_firms$se_race > 0))
  expect_true(all(krw_firms$se_gender > 0))
})

test_that("krw_firms industry/label are complete, ASCII, and 19 industries", {
  expect_false(anyNA(krw_firms$industry))
  expect_false(anyNA(krw_firms$label))
  expect_identical(length(unique(krw_firms$industry)), 19L)
  # CRAN-clean: data must be ASCII (no typographic punctuation leaked through)
  expect_false(any(grepl("[^ -~]", krw_firms$label)))
  expect_false(any(grepl("[^ -~]", krw_firms$industry)))
})

test_that("krw_firms carries sample_stats provenance", {
  ss <- attr(krw_firms, "sample_stats")
  expect_false(is.null(ss))
  expect_equal(as.integer(ss$filtered_firms), 97L)
})

# ---- gp_registry: schema + domains ----------------------------------------
test_that("gp_registry has the C16 contract schema", {
  expect_s3_class(gp_registry, "data.frame")
  expect_identical(nrow(gp_registry), 69L)               # row-count regression guard
  expect_identical(sum(gp_registry$milestone == "M1"), 31L)
  expect_identical(sum(gp_registry$milestone == "M1_NAMES"), 17L)
  expect_identical(sum(gp_registry$milestone == "M2"), 21L)
  req <- c("id", "paper_value", "unit", "tolerance", "class", "milestone")
  expect_true(all(req %in% names(gp_registry)))
  expect_true(all(gp_registry$class %in% c("exact", "banded", "approximate")))
  expect_true(all(gp_registry$milestone %in% c("M1", "M1_NAMES", "M2")))
  expect_false(anyNA(gp_registry$id))
  expect_identical(anyDuplicated(gp_registry$id), 0L)
})

test_that("gp_registry tolerances parse where present; exact rows are tol 0", {
  tol <- suppressWarnings(as.numeric(gp_registry$tolerance))
  # every NON-NA tolerance string parses to a number
  expect_length(which(is.na(tol) & !is.na(gp_registry$tolerance)), 0L)
  # exact rows that carry a tolerance must be exactly 0
  ex <- gp_registry$class == "exact" & !is.na(gp_registry$tolerance)
  expect_true(all(as.numeric(gp_registry$tolerance[ex]) == 0))
  # banded rows that carry a tolerance must be strictly positive
  bd <- gp_registry$class == "banded" & !is.na(gp_registry$tolerance)
  expect_true(all(as.numeric(gp_registry$tolerance[bd]) > 0))
})

# ---- gp_registry: headline values (verbatim CSV ground truth) -------------
test_that("gp_registry headline values match KRW published numbers", {
  hv <- function(id) as.numeric(gp_registry$paper_value[gp_registry$id == id])
  expect_equal(hv("race_baseline_ngrades"), 3)
  expect_equal(hv("gender_baseline_ngrades"), 4)
  expect_equal(hv("names_ngrades_l025"), 2)
  expect_equal(hv("scale_n_firms_graded"), 97)
  expect_equal(hv("scale_n_names"), 76)
  expect_equal(hv("race_n_industries"), 19)
  expect_equal(hv("t3_race_ni_beta"), 0.510)
  expect_equal(hv("t3b_race_withinshare"), 0.366)
  expect_equal(hv("race_industry_ngrades"), 4)
  expect_equal(hv("gender_industry_ngrades"), 5)
  # banded-row value regression guards (verbatim paper_values.csv)
  expect_equal(hv("race_baseline_dr"), 3.9)
  expect_equal(hv("gender_baseline_dr"), 1.8)
  expect_equal(hv("race_baseline_r2"), 25)
  expect_equal(hv("race_baseline_betweengrade_sd"), 0.034)
  expect_equal(hv("race_baseline_worst_n"), 2)
  expect_equal(hv("race_baseline_best_n"), 14)
  expect_equal(hv("f5_genuineparts_theta"), 0.33)
  expect_equal(hv("f5_charterspectrum_condrank_baseline"), 97)
  expect_equal(hv("f6_buildersfirstsource_postmean_baseline"), 0.90)
  expect_equal(hv("f6_ascena_condrank_baseline"), 96)
})

test_that("F5/F6 report-card cells use verbatim ids and milestone scope", {
  report_rows <- gp_registry[grepl("^f[56]_", gp_registry$id), , drop = FALSE]
  expect_identical(nrow(report_rows), 16L)
  expect_equal(sum(report_rows$milestone == "M1"), 13L)
  expect_equal(sum(report_rows$milestone == "M2"), 3L)
  expect_true(all(report_rows$class %in% c("exact", "banded")))
  expect_true(all(c(
    "f5_genuineparts_theta",
    "f5_genuineparts_postmean_baseline",
    "f5_charterspectrum_condrank_baseline",
    "f6_buildersfirstsource_theta",
    "f6_ascena_condrank_baseline"
  ) %in% report_rows$id))
  expect_setequal(
    report_rows$id[report_rows$milestone == "M2"],
    c(
      "f5_genuineparts_postmean_industry",
      "f5_autonation_condrank_industry",
      "f6_buildersfirstsource_postmean_industry"
    )
  )
})

test_that("gp_registry NPMLE atoms use verbatim CSV ids with per-atom tols", {
  hv  <- function(id) as.numeric(gp_registry$paper_value[gp_registry$id == id])
  tol <- function(id) gp_registry$tolerance[gp_registry$id == id]
  expect_equal(hv("names_npmle_mass1"), 0.226)
  expect_equal(hv("names_npmle_mass2"), 0.244)
  expect_equal(hv("names_npmle_mass3"), 0.260)
  expect_equal(tol("names_npmle_mass1"), "0.001")
  expect_equal(tol("names_npmle_mass2"), "0.002")
  expect_equal(tol("names_npmle_mass3"), "0.005")
})

test_that("names paper targets are deferred to the M1_NAMES milestone", {
  names_rows <- gp_registry[
    grepl("^names_", gp_registry$id) | gp_registry$id == "scale_n_names",
    ,
    drop = FALSE
  ]
  expect_identical(nrow(names_rows), 17L)
  expect_true(all(names_rows$milestone == "M1_NAMES"))
  expect_false(any(names_rows$id %in% gp_registry$id[gp_registry$milestone == "M1"]))
  expect_true(all(c(
    "scale_n_names",
    "names_ngrades_l025",
    "names_dr_l025",
    "names_tau_l025",
    "names_grade_r2",
    "names_npmle_mass1",
    "names_npmle_mass2",
    "names_npmle_mass3"
  ) %in% names_rows$id))
})

test_that("gp_registry uses the MC-aware band for approximate rows", {
  # The design widens approximate M2 rows beyond the CSV's published precision.
  appr <- gp_registry[gp_registry$class == "approximate", ]
  expect_true(nrow(appr) >= 1)
  expect_equal(gp_registry$tolerance[gp_registry$id == "race_industry_dr"], "1")
})

test_that("gp_registry excludes pairwise-pi rows (K05-K07 are golden-master)", {
  expect_length(grep("pi_golden|pairwise", gp_registry$id), 0L)
})

test_that("synthetic guard rows (N9, N10) are present", {
  n9 <- gp_registry[gp_registry$id == "n9_strict4_ngrades", , drop = FALSE]
  expect_identical(nrow(n9), 1L)
  expect_equal(as.numeric(n9$paper_value), 4)
  expect_identical(n9$class, "exact")
  for (id in c("n10_race_mult_ngrades", "n10_race_mult_support",
               "n10_gender_add_ngrades", "n10_gender_add_support")) {
    expect_true(id %in% gp_registry$id)
  }
  # the N10 support guards carry TBD (NA) targets at this stage
  expect_true(is.na(gp_registry$paper_value[gp_registry$id == "n10_race_mult_support"]))
})

# ---- archive fixtures (headerless) ----------------------------------------
test_that("archive fixtures load to documented dims (headerless)", {
  fixdir <- system.file("extdata/fixtures", package = "gradepath")
  skip_if(fixdir == "" || !file.exists(file.path(fixdir, "Pi_groupfx1_race.csv")),
          "fixtures not installed")
  pi_race <- as.matrix(read.csv(file.path(fixdir, "Pi_groupfx1_race.csv"),
                                header = FALSE))
  expect_identical(dim(pi_race), c(97L, 97L))
  expect_false(anyNA(pi_race))
  post <- read.csv(file.path(fixdir, "posteriors_groupfx1_race.csv"), header = FALSE)
  expect_identical(dim(post), c(97L, 4L))
  gth <- read.csv(file.path(fixdir, "g_theta_groupfx1_race.csv"), header = FALSE)
  expect_identical(dim(gth), c(250L, 2L))
})

test_that("one-level (groupfx0) and gender fixtures are present and numeric", {
  fixdir <- system.file("extdata/fixtures", package = "gradepath")
  skip_if(fixdir == "" || !file.exists(file.path(fixdir, "Pi_groupfx0_race.csv")),
          "fixtures not installed")
  for (f in c("Pi_groupfx0_race.csv", "posteriors_groupfx0_race.csv",
              "Pi_groupfx1_gender.csv", "posteriors_groupfx1_gender.csv")) {
    m <- as.matrix(read.csv(file.path(fixdir, f), header = FALSE))
    expect_true(is.numeric(m) && !anyNA(m), info = f)
  }
})
