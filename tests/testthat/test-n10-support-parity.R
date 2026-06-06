# =============================================================================
# test-n10-support-parity.R -- N10 two-level support guards
# =============================================================================

.n10_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}

.n10_pairwise <- function(...) .n10_get("gp_twolevel_pairwise")(...)
.n10_grade <- function(...) .n10_get("gp_twolevel_grade")(...)
.n10_push <- function(...) .n10_get("gp_pushforward_theta")(...)

.n10_control <- function() {
  gp_control(backend = "highs")
}

.n10_strict_matrix <- function(n, prefix) {
  ids <- paste0(prefix, seq_len(n))
  m <- matrix(0, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      m[i, j] <- if (i < j) 0.97 else 0.03
    }
  }
  m
}

.n10_case <- function(characteristic = c("race", "gender")) {
  characteristic <- match.arg(characteristic)
  if (identical(characteristic, "race")) {
    list(
      characteristic = "race",
      n = 4L,
      expected_grade_count = 4L,
      supp_xi = c(0.50, 1.40),
      g_xi = c(0.35, 0.65),
      supp_eta = c(0.75, 1.25),
      g_eta = c(0.45, 0.55),
      s = c(0.80, 1.10, 1.40, 0.95),
      mu = 0,
      beta = 0.50
    )
  } else {
    list(
      characteristic = "gender",
      n = 5L,
      expected_grade_count = 5L,
      supp_xi = c(-0.60, 0.20),
      g_xi = c(0.40, 0.60),
      supp_eta = c(-0.30, 0.70),
      g_eta = c(0.50, 0.50),
      s = c(0.80, 1.05, 1.20, 1.45, 0.90),
      mu = 0.12,
      beta = 0.65
    )
  }
}

.n10_expected_support <- function(case) {
  if (identical(case$characteristic, "race")) {
    vals <- numeric(0)
    xi_eta <- as.vector(outer(case$supp_xi, case$supp_eta, `*`))
    for (s in case$s) vals <- c(vals, (s^case$beta) * xi_eta)
  } else {
    vals <- numeric(0)
    xi_eta <- as.vector(outer(case$supp_xi, case$supp_eta, `+`))
    for (s in case$s) vals <- c(vals, case$mu + (s^case$beta) * xi_eta)
  }
  range(vals)
}

.n10_run <- function(case) {
  bundle <- .n10_pairwise(
    Pi_theta = .n10_strict_matrix(case$n, paste0(case$characteristic, "_")),
    Pi_bar = .n10_strict_matrix(2L, paste0(case$characteristic, "_industry_")),
    control = .n10_control()
  )
  graded <- .n10_grade(
    bundle,
    control = .n10_control(),
    lambda_grid = c(0.25, 1),
    build_report_cards = FALSE
  )
  push <- .n10_push(
    supp_xi = case$supp_xi,
    g_xi = case$g_xi,
    supp_eta = case$supp_eta,
    g_eta = case$g_eta,
    s = case$s,
    mu = case$mu,
    beta = case$beta,
    characteristic = case$characteristic,
    supp_pts_theta = 41L
  )
  list(bundle = bundle, graded = graded, push = push)
}

test_that("N10 registry rows carry the synthetic exact/banded M2 contract", {
  rows <- gp_registry[gp_registry$id %in% c(
    "n10_race_mult_ngrades",
    "n10_race_mult_support",
    "n10_gender_add_ngrades",
    "n10_gender_add_support"
  ), ]
  rows <- rows[match(c("n10_race_mult_ngrades",
                       "n10_race_mult_support",
                       "n10_gender_add_ngrades",
                       "n10_gender_add_support"),
                     rows$id), ]

  expect_equal(nrow(rows), 4L)
  expect_identical(as.character(rows$milestone), rep("M2", 4L))
  expect_identical(as.character(rows$unit), c("count", "other", "count", "other"))
  expect_identical(as.character(rows$class), c("exact", "banded", "exact", "banded"))
  expect_equal(as.numeric(rows$tolerance[c(1L, 3L)]), c(0, 0), tolerance = 0)
  expect_true(all(is.na(rows$paper_value)))
})

test_that("N10 race multiplicative guard has exact grade count and banded support endpoints", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  case <- .n10_case("race")
  out <- .n10_run(case)
  expected <- .n10_expected_support(case)
  observed <- range(out$push$support)
  band <- 1e-12

  expect_identical(out$graded$industry_rfe$grade_count,
                   case$expected_grade_count)
  expect_equal(out$graded$industry_rfe$selected_grade$assignment$grade,
               seq_len(case$expected_grade_count))
  expect_lte(max(abs(observed - expected)), band)
  expect_identical(out$push$diagnostics$characteristic, "race")
  expect_identical(out$push$diagnostics$method,
                   "native-two-level-get_g_theta")
})

test_that("N10 gender additive guard has exact grade count and banded support endpoints", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  case <- .n10_case("gender")
  out <- .n10_run(case)
  expected <- .n10_expected_support(case)
  observed <- range(out$push$support)
  band <- 1e-12

  expect_identical(out$graded$industry_rfe$grade_count,
                   case$expected_grade_count)
  expect_equal(out$graded$industry_rfe$selected_grade$assignment$grade,
               seq_len(case$expected_grade_count))
  expect_lte(max(abs(observed - expected)), band)
  expect_identical(out$push$diagnostics$characteristic, "gender")
  expect_identical(out$push$diagnostics$method,
                   "native-two-level-get_g_theta")
})
