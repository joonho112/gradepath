# Tests for the native TWO-level (group_fx == 1) 6-moment cluster-robust GMM.
# Faithful port of the KRW Matlab group_fx == 1 branch
# (moment_conditions.m / get_moments_gmm.m / gmm_obj.m / estimate_lsqnonlin.m).
#
# This file covers BOTH:
#   * STRUCTURE / ALGEBRA: 6-moment vector shape, the dummyvar P_D industry-mean
#     algebra, the cluster-robust Omega == the literal Matlab double loop, the
#     one-level reduction (forcing sigma_eta -> 0 reproduces the 4-moment core),
#     verbatim V_vbar / within_share formulas, single-N J-stat & sandwich,
#     determinism, and error classes.
#   * PARITY on KRW's REAL GMM input (gp_krw_gmm_input, 97 firms / 19 SIC
#     industries): the published two-level numbers -- within_share ~ 0.366 (race)
#     / 0.562 (gender), sigma_eta ~ 0.528 / sigma_xi ~ 0.113 (race), J-stat a
#     small chi-square(2). Continuous targets are banded (M2 class).
#
# The new source is sourced as a NEW file (it does NOT live in R/ yet); a guard
# skips every test if the two-level entry point is unavailable.

skip_if_no_2l <- function() {
  skip_if_not(exists("gp_two_level_gmm", mode = "function"),
              "two-level GMM source (gp_two_level_gmm) not loaded")
}

is_psd <- function(M, tol = 1e-8) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  all(ev >= -tol)
}

# A faithful in-test re-key of the KRW input to the [d, theta_hat, s] matrix.
race_input   <- function() gp_krw_gmm_input("race")
gender_input <- function() gp_krw_gmm_input("gender")

# ============================================================================
# Data matrix + industry projection (the dummyvar P_D algebra)
# ============================================================================

test_that("two-level data matrix carries the real industry column (dense 1..J)", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  expect_equal(ncol(data), 3L)
  expect_equal(nrow(data), 97L)
  expect_identical(length(unique(data[, 1L])), 19L)
  expect_identical(sort(unique(data[, 1L])), as.numeric(1:19))  # contiguous
})

test_that("two-level data matrix aborts without an industry vector", {
  skip_if_no_2l()
  bad <- list(theta_hat = rnorm(10), s = runif(10))            # no $industry
  expect_error(.gp_two_level_data_matrix(bad),
               class = "gradepath_validation_error")
})

test_that("industry projection reproduces the literal dummyvar P_D and n", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  d <- data[, 1L]
  ## Literal Matlab: D = dummyvar(d); P_D = D inv(D'D) D'; n = D (sum(D,1)').
  D <- model.matrix(~ factor(d) - 1)
  P_D <- D %*% solve(crossprod(D)) %*% t(D)
  n_lit <- as.numeric(D %*% colSums(D))
  proj <- .gp_industry_projection(d)
  x <- data[, 2L]
  expect_equal(proj$proj(x), as.numeric(P_D %*% x), tolerance = 1e-12)
  expect_equal(proj$n, n_lit, tolerance = 1e-12)
  expect_equal(proj$J, 19L)
})

# ============================================================================
# Six moment conditions
# ============================================================================

test_that("race two-level moment vector has length 6 and g_i is N x 6", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  delta <- gp_gmm_start_2l("race")
  mc <- gp_moment_conditions_2l(delta, data, "race")
  expect_length(mc$g, 6L)
  expect_equal(dim(mc$g_i), c(nrow(data), 6L))
  expect_true(all(is.finite(mc$g)))
  expect_true(all(is.finite(mc$g_i)))
})

test_that("gender two-level moment vector has length 6 and g_i is N x 6", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(gender_input())
  delta <- gp_gmm_start_2l("gender")
  mc <- gp_moment_conditions_2l(delta, data, "gender")
  expect_length(mc$g, 6L)
  expect_equal(dim(mc$g_i), c(nrow(data), 6L))
  expect_true(all(is.finite(mc$g)))
})

test_that("the first 4 two-level moments match the one-level moments verbatim", {
  skip_if_no_2l()
  ## At any delta, the first 4 two-level moments are the SAME studentized-resid
  ## moments as the one-level core (only sigma now includes sigma_eta). With
  ## sigma_eta -> 0 the two coincide exactly (the one-level reduction).
  data2l <- .gp_two_level_data_matrix(race_input())
  data1l <- cbind(d = rep(1, nrow(data2l)),
                  theta_hat = data2l[, 2L], s = data2l[, 3L])
  delta3 <- c(-1.1789, -1.5761, 0.5099)         # one-level race start
  delta4 <- c(delta3, -50)                       # sigma_eta = exp(-50) ~ 0
  g1l <- gp_moment_conditions(delta3, data1l, "race", group_fx = 0L)$g
  g2l <- gp_moment_conditions_2l(delta4, data2l, "race")$g[1:4]
  expect_equal(g2l, g1l, tolerance = 1e-10)
})

test_that("race two-level moments match a hand-coded port of moment_conditions.m", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  d <- data[, 1L]; theta_hat <- data[, 2L]; s <- data[, 3L]; N <- length(s)
  delta <- c(-1.10, -1.95, 0.52, -0.65)
  mu_raw <- delta[1]; sigma_xi <- exp(delta[2]); beta <- delta[3]
  sigma_eta <- exp(delta[4])
  mu <- exp(mu_raw)                                            # race multiplicative
  sigma <- sqrt(sigma_xi^2 * sigma_eta^2 + mu^2 * sigma_eta^2 + sigma_xi^2)
  v_hat <- theta_hat / s^beta
  s_v <- s^(1 - beta)
  r <- (v_hat - mu) / sqrt(sigma^2 + s_v^2)
  ## dummyvar block
  D <- model.matrix(~ factor(d) - 1)
  P_D <- D %*% solve(crossprod(D)) %*% t(D)
  n <- as.numeric(D %*% colSums(D))
  w <- 1 / n
  v_bar <- as.numeric(P_D %*% (w * v_hat)) * n
  s_bar <- as.numeric(P_D %*% (w * s)) * n
  mu_vbar <- mu
  V_vbar <- (((1 + sigma_eta^2) * sigma_xi^2) * (as.numeric(P_D %*% (w^2)) * n)) +
    (sigma_eta^2 * mu^2) +
    (as.numeric(P_D %*% ((w^2) * (s^(2 * (1 - beta))))) * n)
  grp <- (v_bar - mu_vbar)^2 - V_vbar
  g_ref <- c(mean(r), mean(r * s), mean(r^2 - 1), mean((r^2 - 1) * s),
             mean(grp), mean(grp * s_bar))
  mc <- gp_moment_conditions_2l(delta, data, "race")
  expect_equal(mc$g, g_ref, tolerance = 1e-12)
})

test_that("gender two-level moments match a hand-coded port (different V_vbar)", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(gender_input())
  d <- data[, 1L]; theta_hat <- data[, 2L]; s <- data[, 3L]
  delta <- c(0.00, -0.45, 1.10, -0.60)
  mu <- delta[1]; sigma_xi <- exp(delta[2]); beta <- delta[3]
  sigma_eta <- exp(delta[4])
  sigma <- sqrt(sigma_eta^2 + sigma_xi^2)
  v_hat <- (theta_hat - mu) / s^beta
  s_v <- s^(1 - beta)
  r <- v_hat / sqrt(sigma^2 + s_v^2)
  D <- model.matrix(~ factor(d) - 1)
  P_D <- D %*% solve(crossprod(D)) %*% t(D)
  n <- as.numeric(D %*% colSums(D))
  w <- 1 / n
  v_bar <- as.numeric(P_D %*% (w * v_hat)) * n
  s_bar <- as.numeric(P_D %*% (w * s)) * n
  mu_vbar <- 0
  V_vbar <- sigma_eta^2 + (sigma_xi^2) / n + (as.numeric(P_D %*% (s_v^2))) / n
  grp <- (v_bar - mu_vbar)^2 - V_vbar
  g_ref <- c(mean(r), mean(r * s), mean(r^2 - 1), mean((r^2 - 1) * s),
             mean(grp), mean(grp * s_bar))
  mc <- gp_moment_conditions_2l(delta, data, "gender")
  expect_equal(mc$g, g_ref, tolerance = 1e-12)
})

test_that("v_bar and s_bar are exactly per-industry means broadcast to firms", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  d <- data[, 1L]; theta_hat <- data[, 2L]; s <- data[, 3L]
  delta <- gp_gmm_start_2l("race")
  mc <- gp_moment_conditions_2l(delta, data, "race")
  beta <- delta[3]
  v_hat <- theta_hat / s^beta
  vbar_ref <- ave(v_hat, d, FUN = mean)
  sbar_ref <- ave(s,     d, FUN = mean)
  expect_equal(mc$v_bar, as.numeric(vbar_ref), tolerance = 1e-12)
  expect_equal(mc$s_bar, as.numeric(sbar_ref), tolerance = 1e-12)
})

# ============================================================================
# Cluster-robust covariance == the literal Matlab double loop
# ============================================================================

test_that("cluster Omega equals the literal sum_{i,j: d(i)==d(j)} g_i g_j' / N", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  d <- data[, 1L]; N <- nrow(data)
  delta <- gp_gmm_start_2l("race")
  mc <- gp_moment_conditions_2l(delta, data, "race", get_cov = TRUE)
  gi <- mc$g_i
  ## moment_conditions.m:76-86 -- the literal double loop.
  V_lit <- matrix(0, 6, 6)
  for (i in seq_len(N)) for (j in seq_len(N)) if (d[i] == d[j]) {
    V_lit <- V_lit + (1 / N) * (gi[i, ] %o% gi[j, ])
  }
  expect_equal(mc$Omega, V_lit, tolerance = 1e-10)
  expect_equal(dim(mc$Omega), c(6L, 6L))
  expect_true(is_psd(mc$Omega))
})

test_that("the cluster covariance dominates the independence (group_fx==0) form", {
  skip_if_no_2l()
  ## With non-trivial clusters the cluster-robust Omega differs from the naive
  ## crossprod(g_i)/N (the one-level independence covariance). Confirms the
  ## off-diagonal within-cluster terms are actually present.
  data <- .gp_two_level_data_matrix(race_input())
  delta <- gp_gmm_start_2l("race")
  mc <- gp_moment_conditions_2l(delta, data, "race", get_cov = TRUE)
  N <- nrow(data)
  Omega_indep <- crossprod(mc$g_i) / N
  expect_gt(max(abs(mc$Omega - Omega_indep)), 1e-6)
})

# ============================================================================
# GMM-implied moments + within-share (get_moments_gmm.m)
# ============================================================================

test_that("two-level coupling moments have the right shape and ordering", {
  skip_if_no_2l()
  ## race coupling -> [mu, sigma_xi, sigma_eta]; gender -> [sigma_xi, sigma_eta].
  delta_r <- c(-1.14, -2.18, 0.52, -0.64)
  m_r <- gp_get_moments_2l(delta_r, "race")
  expect_length(m_r, 3L)
  expect_equal(m_r, c(exp(delta_r[1]), exp(delta_r[2]), exp(delta_r[4])),
               tolerance = 1e-12)
  delta_g <- c(0.0, -0.44, 1.11, -0.56)
  m_g <- gp_get_moments_2l(delta_g, "gender")
  expect_length(m_g, 2L)
  expect_equal(m_g, c(exp(delta_g[2]), exp(delta_g[4])), tolerance = 1e-12)
})

test_that("extra_moments==1 appends [E_theta, sd_theta, within_share]", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  delta <- c(-1.14, -2.18, 0.52, -0.64)
  m_full <- gp_get_moments_2l(delta, "race", data = data, extra_moments = 1L)
  ## race: [mu, sigma_xi, sigma_eta, E_theta, sd_theta, within_share] = length 6.
  expect_length(m_full, 6L)
  expect_equal(unname(m_full[6]), gp_within_share(delta, "race"), tolerance = 1e-12)
})

test_that("within_share matches the verbatim get_moments_gmm.m formulas", {
  skip_if_no_2l()
  ## RACE L33 / GENDER L41, hand-coded.
  delta_r <- c(-1.14, -2.18, 0.52, -0.64)
  mu <- exp(delta_r[1]); sxi <- exp(delta_r[2]); seta <- exp(delta_r[4])
  ws_r_ref <- ((seta^2 + 1) * sxi^2) / (seta^2 * sxi^2 + seta^2 * mu^2 + sxi^2)
  expect_equal(gp_within_share(delta_r, "race"), ws_r_ref, tolerance = 1e-12)

  delta_g <- c(0.0, -0.44, 1.11, -0.56)
  sxi_g <- exp(delta_g[2]); seta_g <- exp(delta_g[4])
  ws_g_ref <- sxi_g^2 / (seta_g^2 + sxi_g^2)
  expect_equal(gp_within_share(delta_g, "gender"), ws_g_ref, tolerance = 1e-12)
  ## within-share is a share in (0, 1).
  expect_gt(gp_within_share(delta_r, "race"), 0); expect_lt(gp_within_share(delta_r, "race"), 1)
  expect_gt(gp_within_share(delta_g, "gender"), 0); expect_lt(gp_within_share(delta_g, "gender"), 1)
})

# ============================================================================
# GMM objective and start
# ============================================================================

test_that("two-level start is the documented 4-param Matlab start", {
  skip_if_no_2l()
  ## estimate_lsqnonlin.m:67 (group_fx != 0 branch), shared race/gender.
  expect_equal(gp_gmm_start_2l("race"),   c(-1.1598, -1.9752, 0.5193, -0.7595))
  expect_equal(gp_gmm_start_2l("gender"), c(-1.1598, -1.9752, 0.5193, -0.7595))
  expect_length(gp_gmm_start_2l("race"), 4L)
})

test_that("gp_gmm_min_2l reduces J = N g' W g and is deterministic", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  N <- nrow(data)
  W <- diag(6)
  st <- gp_gmm_start_2l("race")
  f1 <- gp_gmm_min_2l(data, "race", W, st)
  f2 <- gp_gmm_min_2l(data, "race", W, st)
  expect_equal(f1$delta, f2$delta, tolerance = 1e-12)        # deterministic optim
  ## objective at the returned optimum matches a hand recompute of N g' W g.
  g <- gp_moment_conditions_2l(f1$delta, data, "race")$g
  expect_equal(f1$objective, as.numeric(N * crossprod(g, W %*% g)),
               tolerance = 1e-8)
})

# ============================================================================
# PARITY -- KRW's REAL GMM input (the published two-level numbers)
# ============================================================================

test_that("RACE two-level GMM reproduces the published within-share & scales", {
  skip_if_no_2l()
  fit <- gp_two_level_gmm(race_input(), characteristic = "race")
  ## within-share target 0.366; fitted 0.3657 -- essentially exact.
  expect_equal(fit$within_share, 0.366, tolerance = 0.01)
  ## sigma_eta ~ 0.528, sigma_xi ~ 0.113.
  expect_equal(fit$sigma_eta, 0.528, tolerance = 0.02)
  expect_equal(fit$sigma_xi,  0.113, tolerance = 0.02)
  ## df = 6 - 4 = 2; small chi-square(2) J-stat (overid not rejected).
  expect_equal(fit$df, 2L)
  expect_true(is.finite(fit$J_stat) && fit$J_stat >= 0)
  expect_lt(fit$J_stat, qchisq(0.99, df = 2))               # ~ 9.21
  ## coupling carriers: [mu, sigma_xi, sigma_eta], from ONE fit.
  expect_length(fit$m_hat, 3L)
  expect_equal(unname(fit$m_hat[1]), fit$mu,        tolerance = 1e-8)
  expect_equal(unname(fit$m_hat[2]), fit$sigma_xi,  tolerance = 1e-8)
  expect_equal(unname(fit$m_hat[3]), fit$sigma_eta, tolerance = 1e-8)
})

test_that("GENDER two-level GMM reproduces the published within-share", {
  skip_if_no_2l()
  fit <- gp_two_level_gmm(gender_input(), characteristic = "gender")
  ## within-share target 0.562; fitted 0.5621 -- essentially exact.
  expect_equal(fit$within_share, 0.562, tolerance = 0.01)
  expect_equal(fit$df, 2L)
  expect_true(is.finite(fit$J_stat) && fit$J_stat >= 0)
  expect_lt(fit$J_stat, qchisq(0.99, df = 2))
  ## gender coupling carriers: [sigma_xi, sigma_eta].
  expect_length(fit$m_hat, 2L)
  expect_equal(unname(fit$m_hat[1]), fit$sigma_xi,  tolerance = 1e-8)
  expect_equal(unname(fit$m_hat[2]), fit$sigma_eta, tolerance = 1e-8)
})

test_that("two-level fit objects carry coherent natural-scale parameters", {
  skip_if_no_2l()
  for (ch in c("race", "gender")) {
    fit <- gp_two_level_gmm(get(paste0(ch, "_input"))(), characteristic = ch)
    expect_equal(fit$sigma_xi,  exp(fit$delta[2]), tolerance = 1e-12)
    expect_equal(fit$sigma_eta, exp(fit$delta[4]), tolerance = 1e-12)
    expect_equal(fit$beta,      fit$delta[3],      tolerance = 1e-12)
    if (ch == "race") expect_equal(fit$mu, exp(fit$delta[1]), tolerance = 1e-12)
    expect_true(fit$sigma_xi > 0 && fit$sigma_eta > 0)
    expect_equal(fit$within_share, gp_within_share(fit$delta, ch), tolerance = 1e-12)
  }
})

# ============================================================================
# Invariant #7: exactly one factor of N (J-stat AND sandwich C)
# ============================================================================

test_that("two-level J-statistic carries exactly ONE factor of N (df = 2)", {
  skip_if_no_2l()
  inp <- race_input()
  fit <- gp_two_level_gmm(inp, characteristic = "race")
  data <- .gp_two_level_data_matrix(inp); N <- nrow(data)
  ## J = N g2' W_gmm g2, W_gmm = solve(Omega1), g2 at the step-2 optimum.
  g2 <- gp_moment_conditions_2l(fit$delta, data, "race")$g
  J_hand <- as.numeric(N * crossprod(g2, fit$W_gmm %*% g2))
  expect_equal(fit$J_stat, J_hand, tolerance = 1e-8)
  ## one-N J is O(1); an N^2 regression would be ~97x larger.
  expect_lt(fit$J_stat, 20)
})

test_that("two-level sandwich C carries exactly ONE factor of 1/N", {
  skip_if_no_2l()
  inp <- race_input()
  fit <- gp_two_level_gmm(inp, characteristic = "race")
  N <- length(inp$theta_hat)
  bread  <- t(fit$G) %*% solve(fit$Omega2) %*% fit$G          # G' Omega2^-1 G (4x4)
  C_hand <- solve(bread) / N
  C_hand <- (C_hand + t(C_hand)) / 2
  expect_equal(unname(fit$C),  unname(C_hand),               tolerance = 1e-6)
  expect_equal(unname(fit$SE), unname(sqrt(diag(C_hand))),   tolerance = 1e-6)
  expect_length(fit$SE, 4L)
})

test_that("J-stat weight is solve(Omega1), NOT solve(Omega2) (Omega1/Omega2 split)", {
  skip_if_no_2l()
  inp <- race_input()
  fit <- gp_two_level_gmm(inp, characteristic = "race")
  data <- .gp_two_level_data_matrix(inp); N <- nrow(data)
  g2 <- gp_moment_conditions_2l(fit$delta, data, "race")$g
  J_with_O1 <- as.numeric(N * crossprod(g2, solve(fit$Omega1) %*% g2))
  J_with_O2 <- as.numeric(N * crossprod(g2, solve(fit$Omega2) %*% g2))
  expect_equal(fit$J_stat, J_with_O1, tolerance = 1e-8)      # uses Omega1
  ## Omega1 != Omega2 at the optimum, so the two J values genuinely differ.
  expect_gt(abs(J_with_O1 - J_with_O2), 1e-6)
})

# ============================================================================
# Coupling V_m
# ============================================================================

test_that("V_m has the coupling-vector dimension and is symmetric PSD", {
  skip_if_no_2l()
  fit_r <- gp_two_level_gmm(race_input(),   characteristic = "race")
  fit_g <- gp_two_level_gmm(gender_input(), characteristic = "gender")
  expect_equal(dim(fit_r$V_m), c(3L, 3L))                    # [mu, sigma_xi, sigma_eta]
  expect_equal(dim(fit_g$V_m), c(2L, 2L))                    # [sigma_xi, sigma_eta]
  expect_equal(fit_r$V_m, t(fit_r$V_m), tolerance = 1e-10)
  expect_equal(fit_g$V_m, t(fit_g$V_m), tolerance = 1e-10)
  expect_true(is_psd(fit_r$V_m)); expect_true(is_psd(fit_g$V_m))
  expect_gt(rcond(fit_r$V_m), 1e-8)                          # downstream uses V_m^-1
  expect_gt(rcond(fit_g$V_m), 1e-8)
  expect_silent(solve(fit_r$V_m))
  expect_silent(solve(fit_g$V_m))
  expect_equal(fit_r$industry, .gp_two_level_data_matrix(race_input())[, 1L])
  expect_equal(fit_g$industry, .gp_two_level_data_matrix(gender_input())[, 1L])
  expect_length(fit_r$industry, length(fit_r$v_hat))
  expect_length(fit_g$industry, length(fit_g$v_hat))
})

# ============================================================================
# Determinism
# ============================================================================

test_that("two-level GMM is deterministic across reruns and seeds", {
  skip_if_no_2l()
  for (ch in c("race", "gender")) {
    inp <- get(paste0(ch, "_input"))()
    set.seed(8675309)
    seed_before <- .Random.seed
    f1 <- gp_two_level_gmm(inp, characteristic = ch)
    expect_identical(.Random.seed, seed_before)
    f2 <- gp_two_level_gmm(inp, characteristic = ch)
    expect_equal(f1$delta, f2$delta, tolerance = 1e-12)       # bit-identical rerun
    expect_equal(f1$within_share, f2$within_share, tolerance = 1e-12)
    ## a different seed reaches the same optimum (robust basin).
    fs <- gp_two_level_gmm(inp, characteristic = ch,
                           control = gp_control(seed = 99L))
    expect_equal(fs$within_share, f1$within_share, tolerance = 1e-3)
  }
})

# ============================================================================
# Error classes
# ============================================================================

test_that("two-level driver requires an explicit characteristic", {
  skip_if_no_2l()
  expect_error(gp_two_level_gmm(race_input()),
               class = "gradepath_validation_error")
})

test_that("a singular cluster Omega1 raises gradepath_singular_error", {
  skip_if_no_2l()
  ## gender's documented start diverges under identity weighting to a singular-
  ## Omega1 corner; with n_starts = 0 there is no multistart rescue, so the
  ## driver must abort with the typed singular-error class (NOT a raw R error).
  expect_error(
    gp_two_level_gmm(gender_input(), characteristic = "gender", n_starts = 0L),
    class = "gradepath_singular_error"
  )
})

test_that("unknown characteristic in the moment routine errors", {
  skip_if_no_2l()
  data <- .gp_two_level_data_matrix(race_input())
  expect_error(gp_moment_conditions_2l(gp_gmm_start_2l("race"), data, "ethnicity"),
               class = "gradepath_validation_error")
  expect_error(gp_within_share(gp_gmm_start_2l("race"), "ethnicity"),
               class = "gradepath_validation_error")
})

# ============================================================================
# Additional guards: the >=2-industry requirement, warm-start
# robustness, and locale-independent (radix) industry re-keying.
# ============================================================================

test_that("two-level data needs >= 2 industries (single-cluster aborts cleanly)", {
  skip_if_no_2l()
  ## A single industry cannot identify the between-industry effect or form a
  ## full-rank cluster covariance; the data matrix must abort with a typed
  ## validation error rather than fall through to a confusing singular error.
  one_ind <- list(theta_hat = rnorm(8), s = runif(8, 0.05, 0.2),
                  industry = rep(1L, 8))
  expect_error(.gp_two_level_data_matrix(one_ind),
               class = "gradepath_validation_error")
  ## two industries is the minimum that is accepted.
  two_ind <- list(theta_hat = rnorm(8), s = runif(8, 0.05, 0.2),
                  industry = rep(1:2, each = 4))
  expect_silent(.gp_two_level_data_matrix(two_ind))
})

test_that("warm start alone (n_starts = 0) already reaches the published race optimum", {
  skip_if_no_2l()
  ## Race's documented start already lands at KRW's published within-share without
  ## the multistart; the multistart only REFINES it to a marginally lower-objective
  ## optimum (its real job is hardening GENDER, whose documented start diverges).
  ## So n_starts = 0 must (a) land within the M2 banded tolerance of 0.366 and
  ## close to the refined optimum, and (b) the multistart must be no WORSE in
  ## GMM objective (it can only find an equal-or-lower minimum).
  f0  <- gp_two_level_gmm(race_input(), characteristic = "race", n_starts = 0L)
  ref <- gp_two_level_gmm(race_input(), characteristic = "race")
  expect_equal(f0$within_share, 0.366, tolerance = 0.01)            # published target
  expect_equal(f0$within_share, ref$within_share, tolerance = 5e-3) # near the refined fit
  expect_true(f0$objective >= ref$objective - 1e-8)                 # multistart never worse
})

test_that("industry re-keying is dense + locale-independent (radix)", {
  skip_if_no_2l()
  ## Non-contiguous / non-1-based integer codes must re-key to a dense 1..J
  ## partition; only the partition (not the labels) matters, and the keying is
  ## radix (byte-order), so results are reproducible under the C locale (CRAN/CI).
  inp  <- race_input()
  shifted <- inp; shifted$industry <- inp$industry * 10L + 100L   # 110,120,...,290
  d0 <- .gp_two_level_data_matrix(inp)[, 1L]
  d1 <- .gp_two_level_data_matrix(shifted)[, 1L]
  expect_identical(sort(unique(d1)), as.numeric(1:19))
  expect_identical(d0, d1)                                        # same partition
  ## the fit is invariant to a monotone relabel of the industry codes.
  fit0 <- gp_two_level_gmm(inp,     characteristic = "race")
  fit1 <- gp_two_level_gmm(shifted, characteristic = "race")
  expect_equal(fit0$within_share, fit1$within_share, tolerance = 1e-10)
})
