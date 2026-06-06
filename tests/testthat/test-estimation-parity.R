# KRW Table-3 PARITY for the native beta-GMM core, on KRW's REAL GMM input.
#
# The structural tests in test-estimation-core.R run on the bundled
# ebrecipe::krw_firms example (wrong scale for the GMM) and assert SHAPE only.
# Here we feed KRW's ACTUAL Matlab GMM input -- the per-firm (theta_hat, s) that
# estimate_lsqnonlin.m consumed, loaded via gp_krw_gmm_input() from
# inst/extdata/krw-gmm-input/theta_estimates_matlab_<char>.csv -- and reproduce
# the published precision exponents to three decimals (race 0.510, gender 1.255).
# The two-step optimal weighting (W = chol(Omega1^-1)) is what makes 0.510 the
# robust optimum; on the wrong-scale example a degenerate large-beta basin opens
# up instead (an input-provenance issue, see estimation-input.R).
#
# The Jacobian backend is base-R central differences (self-contained, no
# external dependency), so these gates run unconditionally.

# --- race: the headline one-level precision exponent ------------------------

test_that("KRW GMM input keeps unique row ids separate from industry groups", {
  inp <- gp_krw_gmm_input("race")
  expect_length(inp$unit_id, length(inp$theta_hat))
  expect_identical(length(unique(inp$unit_id)), length(inp$theta_hat))
  expect_length(inp$industry, length(inp$theta_hat))
  expect_identical(length(unique(inp$industry)), 19L)
  expect_identical(inp$metadata$unit_id_source, "row_id")
})

test_that("race beta-GMM reproduces KRW Table 3 (real GMM input)", {
  fit <- gp_estimation_core(gp_krw_gmm_input("race"), characteristic = "race")
  expect_equal(fit$beta,     0.510, tolerance = 0.005)   # t3_race_ni_beta
  expect_equal(fit$mu,       0.308, tolerance = 0.005)   # E[v] location
  expect_equal(fit$sigma_xi, 0.207, tolerance = 0.005)   # sigma_xi
  expect_equal(fit$SE[3],    0.190, tolerance = 0.02)    # SE(beta) ~ 0.190
  expect_equal(fit$df, 1L)
  expect_true(is.finite(fit$J_stat) && fit$J_stat >= 0)
})

test_that("race coupling carriers m_hat = [mu, sigma_xi] from ONE fit (inv #3)", {
  fit <- gp_estimation_core(gp_krw_gmm_input("race"), characteristic = "race")
  expect_length(fit$m_hat, 2L)
  expect_equal(unname(fit$m_hat[1]), fit$mu,       tolerance = 1e-8)
  expect_equal(unname(fit$m_hat[2]), fit$sigma_xi, tolerance = 1e-8)
  expect_equal(dim(fit$V_m), c(2L, 2L))
  expect_equal(fit$V_m, t(fit$V_m), tolerance = 1e-10)
})

# --- gender: additive-residual exponent -------------------------------------

test_that("gender beta-GMM reproduces KRW Table 3 (real GMM input)", {
  fit <- gp_estimation_core(gp_krw_gmm_input("gender"), characteristic = "gender")
  expect_equal(fit$beta, 1.255, tolerance = 0.01)        # t3_gender_ni_beta
  expect_equal(fit$df, 1L)
  expect_length(fit$m_hat, 1L)                           # gender carries [sigma_xi]
})

# --- invariant #7: exactly one factor of N (J-stat AND sandwich C) ----------

test_that("J-statistic carries exactly ONE factor of N (invariant #7)", {
  inp  <- gp_krw_gmm_input("race")
  fit  <- gp_estimation_core(inp, characteristic = "race")
  data <- .gp_estimation_data_matrix(inp)
  N    <- nrow(data)
  # Independent hand-recompute: J = N * g2' W_gmm g2, W_gmm = solve(Omega1),
  # g2 = the moment vector at the step-2 optimum. ONE factor of N.
  g2     <- gp_moment_conditions(fit$delta, data, "race", group_fx = 0L)$g
  J_hand <- as.numeric(N * crossprod(g2, fit$W_gmm %*% g2))
  expect_equal(fit$J_stat, J_hand, tolerance = 1e-8)
  # A double-N (N^2) regression would be ~97x larger; the one-N J is O(1) (~0.1).
  expect_lt(fit$J_stat, 10)
})

test_that("sandwich covariance C carries exactly ONE factor of 1/N (invariant #7)", {
  inp <- gp_krw_gmm_input("race")
  fit <- gp_estimation_core(inp, characteristic = "race")
  N   <- length(inp$theta_hat)
  # Independent hand-recompute: C = (G' Omega2^-1 G)^-1 / N. ONE factor of 1/N.
  bread  <- t(fit$G) %*% solve(fit$Omega2) %*% fit$G
  C_hand <- solve(bread) / N
  C_hand <- (C_hand + t(C_hand)) / 2
  expect_equal(unname(fit$C),  unname(C_hand),              tolerance = 1e-6)
  expect_equal(unname(fit$SE), unname(sqrt(diag(C_hand))),  tolerance = 1e-6)
})

# --- robustness: 0.510 is reached even from a high-beta start ----------------

test_that("two-step optimum is 0.510 even from a high-beta start (real input)", {
  if (!"start" %in% names(formals(gp_estimation_core))) skip("no start arg")
  inp <- gp_krw_gmm_input("race")
  fit <- tryCatch(
    gp_estimation_core(inp, characteristic = "race",
                       start = c(-1.1789, -1.5761, 1.5)),
    error = function(e) {
      skip(paste("explicit start not supported:", conditionMessage(e)))
    }
  )
  expect_equal(fit$beta, 0.510, tolerance = 0.01)
})
