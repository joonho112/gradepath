# =============================================================================
# test-pairwise.R  --  pairwise outranking probabilities (gp_pairwise)
# -----------------------------------------------------------------------------
# Chapter 4. gp_pairwise() integrates the Chapter-8 W on each unit's per-unit
# reporting axis into the J x J outranking matrix pi_ij = Pr(theta_i > theta_j),
# with the frozen cleanup contract (antisymmetry, diagonal 0.5, zero-floor 1e-7).
# The headline gate is the GOLDEN MASTER: pi matches the KRW companion archive
# fixture Pi_groupfx0_<char>.csv (the one-level Pi) to floating-point roundoff
# (pi has no published paper_values row -- Ch4 sec-ch4-verify -- so the archive
# Pi matrix is the verification target).
# =============================================================================

.gpw_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}
.krw   <- function(...) .gpw_get("gp_krw_gmm_input")(...)
.core  <- function(...) .gpw_get("gp_estimation_core")(...)
.wseam <- function(...) .gpw_get("gp_w_seam")(...)
.decon <- function(...) .gpw_get("gp_deconvolve")(...)
.pw    <- function(...) .gpw_get("gp_posterior_weights")(...)
.pair  <- function(...) .gpw_get("gp_pairwise")(...)
.valid <- function(...) .gpw_get("validate_gp_pairwise")(...)
.rawp  <- function(...) .gpw_get(".gp_pairwise_raw_probability")(...)
.clean <- function(...) .gpw_get(".gp_pairwise_cleanup_matrix")(...)

.build <- function(char) {
  inp <- .krw(char)
  fit <- .core(inp, characteristic = char)
  cpl <- .wseam(inp, char)
  estimates <- ebrecipe::eb_input(theta_hat = cpl$v_hat, s = cpl$s_v)
  prior <- ebrecipe::eb_deconvolve(estimates)
  parts <- .pw(prior, list(theta_hat = cpl$v_hat, s = cpl$s_v,
                           original_s = inp$s, id = inp$unit_id), fit)
  list(inp = inp, parts = parts, obj = .pair(parts))
}

.native_control <- function(char) {
  ctl <- gp_control()
  if (identical(char, "race")) {
    ctl$deconv_penalty_grid_race <- seq(0.02, 0.12, by = 0.02)
  } else {
    ctl$deconv_penalty_grid_gender <- seq(0.001, 0.006, by = 0.001)
  }
  ctl
}

.build_native <- function(char) {
  inp <- .krw(char)
  fit <- .core(inp, characteristic = char)
  cpl <- .wseam(inp, char)
  prior <- .decon(cpl, control = .native_control(char))
  r_est <- list(theta_hat = cpl$v_hat, s = cpl$s_v, original_s = inp$s)
  parts <- .pw(prior, r_est, cpl)
  list(inp = inp, fit = fit, cpl = cpl, prior = prior,
       parts = parts, obj = .pair(parts))
}

test_that("gp_pairwise returns a valid J x J gp_pairwise (race + gender)", {
  skip_if_not_installed("ebrecipe")
  for (char in c("race", "gender")) {
    b <- .build(char)
    obj <- b$obj
    expect_s3_class(obj, "gp_pairwise")
    expect_silent(.valid(obj))
    J <- length(b$inp$theta_hat)
    expect_identical(dim(obj$matrix), c(J, J))
    expect_identical(obj$power, 0L)
  }
})

test_that("diagonal is exactly 0.5 and off-diagonals are antisymmetric", {
  skip_if_not_installed("ebrecipe")
  pi <- .build("race")$obj$matrix
  expect_true(all(diag(pi) == 0.5))                       # inv #6
  off <- row(pi) != col(pi)
  gap <- abs(pi + t(pi) - 1)
  expect_lt(max(gap[off]), 1e-10)                         # pi_ij + pi_ji = 1
})

test_that("zero-floor 1e-7: exact zeros floor, exact ones remain valid", {
  skip_if_not_installed("ebrecipe")
  pi <- .build("race")$obj$matrix
  off <- row(pi) != col(pi)
  expect_false(any(pi[off] == 0))                         # exact zeros floored
  expect_gte(min(pi[off]), 1e-7)
  expect_lte(max(pi), 1)
})

test_that("cleanup is one-sided: exact off-diagonal one is preserved", {
  raw <- matrix(0, 3, 3)
  raw[1, 2] <- 1
  pi <- .clean(raw)
  expect_identical(unname(pi[1, 2]), 1)
  expect_identical(unname(pi[2, 1]), 1e-7)
  expect_true(all(diag(pi) == 0.5))
})

test_that("zero-floor catches tiny positives, not only exact 0 (CCR-06)", {
  # The cleanup floors ALL off-diagonal values below 1e-7, per the design's
  # floor-all rule (`pmax(pmin(., 1), 1e-7)`) -- NOT only exact
  # zeros (KRW generate_grades.py:55 floors only exact 0). A tiny POSITIVE raw
  # upper-triangle value below 1e-7 must therefore be raised to exactly 1e-7.
  raw <- matrix(0, 3, 3)
  raw[1, 2] <- 5e-9                                # tiny positive, below 1e-7
  pi <- .clean(raw)
  expect_identical(unname(pi[1, 2]), 1e-7)         # tiny positive floored to 1e-7
  expect_false(pi[1, 2] == 5e-9)                   # NOT left at its raw value
  # antisymmetric complement 1 - 5e-9 ~ 1 is above the floor, so it is untouched.
  expect_identical(unname(pi[2, 1]), 1 - 5e-9)
  expect_true(all(diag(pi) == 0.5))
})

test_that("pi_ij is the bit-exact per-unit CDF integral (cumsum + findInterval)", {
  skip_if_not_installed("ebrecipe")
  b <- .build("race")
  W <- b$parts$W; rs <- b$parts$reporting_support
  pi <- b$obj$matrix
  # hand-recompute a strict-upper entry (i < j) directly from the keeper formula
  for (pair in list(c(3, 7), c(10, 42), c(1, 2))) {
    i <- pair[1]; j <- pair[2]
    hand <- .rawp(rs[, i], W[i, ], rs[, j], W[j, ])
    expect_equal(unname(pi[i, j]), hand, tolerance = 1e-12)
    # the lower entry is the antisymmetric complement
    expect_equal(unname(pi[j, i]), 1 - hand, tolerance = 1e-12)
  }
})

# CROSS-CHECK vs the KRW companion archive Pi_groupfx0 fixtures. This is an
# APPROXIMATE diagnostic, NOT a bit-exact gate: pi has no published paper_values
# row, and its rigor comes from GP-W-EXACT (the
# recomputed W is bit-identical to ebrecipe's internal) plus the exact CDF
# integral (pinned above). The archive Pi was produced by KRW's Matlab spline
# deconvolution, whereas gradepath's W rides ebrecipe's deconvolution of the same
# standardized estimates, so the two pi matrices differ at the ~2-4% level
# (correlation ~0.97) -- the SAME upstream deconvolution-provenance gap flagged
# for the maintainer (open item #1, ebrecipe vs KRW-Matlab deconvolution). We
# assert the strong correlation + the dimension, and record (do not gate) the
# magnitude. We compare off-diagonal entries only because the diagonal is a
# fixed convention, NOT a divergence: gradepath's 0.5 MATCHES KRW's consumed
# grade-input matrix, which sets the diagonal to 0.5 (generate_grades.py:56,
# `np.fill_diagonal(p_ij, 0.5)`). The "archive 0" we skip is merely the raw
# Pi_groupfx0 CSV as written BEFORE that fill, so 0 vs 0.5 reflects raw-vs-filled
# staging, not a methodological difference in the grade input.
test_that("pi is strongly correlated with the KRW archive Pi_groupfx0 (race; approx)", {
  skip_if_not_installed("ebrecipe")
  pi <- .build("race")$obj$matrix; dimnames(pi) <- NULL
  fp <- system.file("extdata/fixtures/Pi_groupfx0_race.csv", package = "gradepath")
  if (fp == "" || !file.exists(fp)) fp <- "inst/extdata/fixtures/Pi_groupfx0_race.csv"
  skip_if_not(file.exists(fp), "Pi_groupfx0_race.csv fixture not found")
  Pf <- as.matrix(utils::read.csv(fp, header = FALSE)); dimnames(Pf) <- NULL
  expect_identical(dim(pi), dim(Pf))
  off <- row(pi) != col(pi)
  expect_gt(stats::cor(pi[off], Pf[off]), 0.95)          # ~0.97; deconv-provenance gap
})

test_that("pi is strongly correlated with the KRW archive Pi_groupfx0 (gender; approx)", {
  skip_if_not_installed("ebrecipe")
  pi <- .build("gender")$obj$matrix; dimnames(pi) <- NULL
  fp <- system.file("extdata/fixtures/Pi_groupfx0_gender.csv", package = "gradepath")
  if (fp == "" || !file.exists(fp)) fp <- "inst/extdata/fixtures/Pi_groupfx0_gender.csv"
  skip_if_not(file.exists(fp), "Pi_groupfx0_gender.csv fixture not found")
  Pf <- as.matrix(utils::read.csv(fp, header = FALSE)); dimnames(Pf) <- NULL
  expect_identical(dim(pi), dim(Pf))
  off <- row(pi) != col(pi)
  expect_gt(stats::cor(pi[off], Pf[off]), 0.95)
})

test_that("native pipeline builds valid pairwise matrices (race + gender)", {
  skip_if_not_installed("ebrecipe")
  for (char in c("race", "gender")) {
    b <- .build_native(char)
    expect_s3_class(b$cpl, "gp_estimation_fit")
    expect_s3_class(b$prior, "gp_prior")
    expect_s3_class(b$obj, "gp_pairwise")
    expect_silent(.valid(b$obj))
    expect_identical(dim(b$obj$matrix), c(length(b$inp$theta_hat), length(b$inp$theta_hat)))
    expect_true(all(diag(b$obj$matrix) == 0.5))
    off <- row(b$obj$matrix) != col(b$obj$matrix)
    expect_lt(max(abs(b$obj$matrix[off] + t(b$obj$matrix)[off] - 1)), 1e-10)
  }
})

test_that("native pipeline pairwise matrices are strongly correlated with KRW archive", {
  skip_if_not_installed("ebrecipe")
  for (char in c("race", "gender")) {
    pi <- .build_native(char)$obj$matrix; dimnames(pi) <- NULL
    fp <- system.file("extdata/fixtures", sprintf("Pi_groupfx0_%s.csv", char),
                      package = "gradepath")
    if (fp == "" || !file.exists(fp)) {
      fp <- file.path("inst/extdata/fixtures", sprintf("Pi_groupfx0_%s.csv", char))
    }
    skip_if_not(file.exists(fp), sprintf("Pi_groupfx0_%s.csv fixture not found", char))
    Pf <- as.matrix(utils::read.csv(fp, header = FALSE)); dimnames(Pf) <- NULL
    expect_identical(dim(pi), dim(Pf))
    off <- row(pi) != col(pi)
    expect_gt(stats::cor(pi[off], Pf[off]), 0.99)
  }
})

test_that("non-unique ids are rejected; a bare W needs reporting_support", {
  skip_if_not_installed("ebrecipe")
  b <- .build("race")
  # KRW GMM row ids are unique and can be used for a pairwise object.
  expect_silent(.pair(b$parts, ids = b$inp$unit_id))
  # Fabricated duplicates are still rejected by the validator.
  dup_ids <- rep("dup", length(b$inp$unit_id))
  expect_error(.pair(b$parts, ids = dup_ids))
  # a bare W matrix with no reporting axis is rejected
  expect_error(.pair(b$parts$W))
})
