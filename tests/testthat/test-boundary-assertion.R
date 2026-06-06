# gp_assert_ebrecipe_boundary(): the seam's loud-failure schema guard
# (the middle layer of the GP-DEC-07-A drift guard).
#
# A real ebrecipe eb_estimates + eb_prior (+ the posterior W) must PASS; each
# schema violation (missing field, length mismatch, non-'r' scale, non-row-
# stochastic W, Pi diag != 0.5) must fail LOUDLY with a gradepath_error.

# Build real ebrecipe objects via the seam (skip if ebrecipe absent).
make_boundary_inputs <- function() {
  est  <- .gp_eb_input(theta_hat = krw_firms$theta_hat_race,
                       s = krw_firms$se_race,
                       unit_id = krw_firms$firm_id)
  prio <- ebrecipe::eb_deconvolve(est)
  W    <- .gp_eb_posterior_weights(est, prio)
  list(est = est, prio = prio, W = W)
}

test_that("a real ebrecipe boundary passes the shape-assertion", {
  skip_if_not_installed("ebrecipe")
  io <- make_boundary_inputs()
  expect_s3_class(.gp_eb_validate_estimates(io$est), "eb_estimates")
  expect_true(gp_assert_ebrecipe_boundary(io$est))
  expect_true(gp_assert_ebrecipe_boundary(io$est, prior = io$prio))
  expect_true(gp_assert_ebrecipe_boundary(io$est, prior = io$prio, weights = io$W))
  # a valid pairwise Pi (diag 0.5) passes
  J  <- 4L
  Pi <- matrix(0.4, J, J); diag(Pi) <- 0.5
  expect_true(gp_assert_ebrecipe_boundary(io$est, pairwise = Pi))
})

test_that("eb_estimates missing a required field fails loudly", {
  skip_if_not_installed("ebrecipe")
  io <- make_boundary_inputs()
  bad <- io$est
  bad$theta_hat <- NULL
  expect_error(gp_assert_ebrecipe_boundary(bad),
               class = "gradepath_error", regexp = "eb_estimates missing")
})

test_that("eb_estimates with mismatched theta_hat/s lengths fails loudly", {
  skip_if_not_installed("ebrecipe")
  io <- make_boundary_inputs()
  bad <- io$est
  bad$s <- bad$s[-1]
  expect_error(gp_assert_ebrecipe_boundary(bad),
               class = "gradepath_error", regexp = "length-matched")
})

test_that("eb_prior missing a required field fails loudly", {
  skip_if_not_installed("ebrecipe")
  io <- make_boundary_inputs()
  bad <- io$prio
  bad$density <- NULL
  expect_error(gp_assert_ebrecipe_boundary(io$est, prior = bad),
               class = "gradepath_error", regexp = "eb_prior missing")
})

test_that("eb_prior with a non-'r' scale fails loudly", {
  skip_if_not_installed("ebrecipe")
  io <- make_boundary_inputs()
  bad <- io$prio
  bad$scale <- "theta"
  expect_error(gp_assert_ebrecipe_boundary(io$est, prior = bad),
               class = "gradepath_error", regexp = "scale must be 'r'")
})

test_that("a non-row-stochastic W fails loudly (orientation guard, inv. 5)", {
  skip_if_not_installed("ebrecipe")
  io <- make_boundary_inputs()
  # transpose -> M x J, columns (not rows) sum to 1 -> rejected
  expect_error(gp_assert_ebrecipe_boundary(io$est, weights = t(io$W)),
               class = "gradepath_error", regexp = "row-stochastic")
  # scaled rows -> rows no longer sum to 1
  bad <- io$W * 2
  expect_error(gp_assert_ebrecipe_boundary(io$est, weights = bad),
               class = "gradepath_error", regexp = "row-stochastic")
})

test_that("a pairwise Pi with wrong diagonal fails loudly (inv. 6)", {
  skip_if_not_installed("ebrecipe")
  io <- make_boundary_inputs()
  J  <- 4L
  Pi <- matrix(0.4, J, J); diag(Pi) <- 0.0   # solver-scale diag, not the stored 0.5
  expect_error(gp_assert_ebrecipe_boundary(io$est, pairwise = Pi),
               class = "gradepath_error", regexp = "diagonal must be 0.5")
  # non-square
  expect_error(gp_assert_ebrecipe_boundary(io$est, pairwise = matrix(0.5, 3, 4)),
               class = "gradepath_error", regexp = "square")
})
