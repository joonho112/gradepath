## tests for the W-seam coupling object (gp_estimation_fit): object validity,
## the chi-square p-value, the report list, V_m PSD, the support caps, and the
## frozen invariant #3 (every slot from the ONE fit -- recompute and compare).
##
## The helper below resolves package internals namespace-first so the file runs
## against a devtools::load_all()-ed package.

## ---- resolve internals -----------------------------------------------------
.gpc_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = if (name == "gp_registry") "any" else "function")
}
.krw       <- function(...) .gpc_get("gp_krw_gmm_input")(...)
.core      <- function(...) .gpc_get("gp_estimation_core")(...)
.couple    <- function(...) .gpc_get("gp_estimation_coupling")(...)
.wseam     <- function(...) .gpc_get("gp_w_seam")(...)
.validate  <- function(...) .gpc_get("validate_gp_estimation_fit")(...)
.datamat   <- function(...) .gpc_get(".gp_estimation_data_matrix")(...)
.moments   <- function(...) .gpc_get("gp_get_moments")(...)

.matlab_xi_caps <- function(fit, char) {
  v <- fit$v_hat
  if (identical(char, "race")) {
    mu <- fit$m_hat[1L]
    sig <- fit$m_hat[2L]
    return(c(lo = 0, hi = min(max(max(v), mu + 5 * sig), mu + 7 * sig)))
  }
  sig <- fit$m_hat[1L]
  c(lo = max(min(min(v), -5 * sig), -7 * sig),
    hi = min(max(max(v),  5 * sig),  7 * sig))
}

## Build both characteristics' coupling objects once.
.build <- function(char) {
  inp <- .krw(char)
  fit <- .core(inp, characteristic = char)
  list(input = inp, fit = fit, obj = .couple(fit, inp))
}

test_that("coupling object is a valid gp_estimation_fit (race + gender)", {
  for (char in c("race", "gender")) {
    b <- .build(char)
    obj <- b$obj
    expect_s3_class(obj, "gp_estimation_fit")
    ## validator passes
    expect_silent(.validate(obj))
    ## exact contract slots, no extras like J_stat / m_report / reject
    expect_setequal(
      names(obj),
      c("beta", "m_hat", "V_m", "v_hat", "s_v", "J", "df", "p_value",
        "report", "caps", "characteristic", "provenance")
    )
    expect_false("J_stat"   %in% names(obj))
    expect_false("m_report" %in% names(obj))
    expect_false("V_report" %in% names(obj))
    expect_false("reject"   %in% names(obj))
    expect_identical(obj$characteristic, char)
    ## m_hat shape: race length 2, gender length 1
    expect_length(obj$m_hat, if (char == "race") 2L else 1L)
    expect_true(all(is.finite(obj$m_hat)))
    ## beta scalar finite
    expect_true(is.numeric(obj$beta) && length(obj$beta) == 1L &&
                  is.finite(obj$beta))
  }
})

test_that("J is finite, df == 1, and p_value == pchisq(J, df) in [0, 1]", {
  for (char in c("race", "gender")) {
    obj <- .build(char)$obj
    expect_true(is.finite(obj$J) && obj$J >= 0)
    ## df is 1 (carried verbatim from the core, where it is the integer 4L - 3L);
    ## assert value equality, not type identity -- the contract is "df == 1".
    expect_equal(obj$df, 1)
    expect_true(is.finite(obj$p_value))
    expect_gte(obj$p_value, 0)
    expect_lte(obj$p_value, 1)
    expect_equal(obj$p_value,
                 stats::pchisq(obj$J, obj$df, lower.tail = FALSE))
  }
})

test_that("invariant #7: J is the core's J_stat verbatim (never rescaled)", {
  for (char in c("race", "gender")) {
    b <- .build(char)
    expect_identical(b$obj$J, b$fit$J_stat)
  }
})

test_that("report is a named list(E_theta, SD_theta) with finite values", {
  for (char in c("race", "gender")) {
    obj <- .build(char)$obj
    expect_true(is.list(obj$report))
    expect_setequal(names(obj$report), c("E_theta", "SD_theta"))
    expect_true(is.finite(obj$report$E_theta))
    expect_true(is.finite(obj$report$SD_theta))
    expect_length(obj$report$E_theta, 1L)
    expect_length(obj$report$SD_theta, 1L)
    ## reporting carries NO covariance
    expect_false("V_report" %in% names(obj$report))
    expect_false("cov" %in% names(obj$report))
  }
})

test_that("V_m is square, symmetric, and positive semidefinite", {
  for (char in c("race", "gender")) {
    obj <- .build(char)$obj
    m <- if (char == "race") 2L else 1L
    expect_true(is.matrix(obj$V_m))
    expect_equal(dim(obj$V_m), c(m, m))
    expect_true(all(is.finite(obj$V_m)))
    expect_lte(max(abs(obj$V_m - t(obj$V_m))),
               1e-8 * max(1, max(abs(obj$V_m))))
    ev <- eigen(obj$V_m, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev >= -1e-8 * max(1, max(abs(obj$V_m)))))
  }
})

test_that("caps are finite with lo < hi", {
  for (char in c("race", "gender")) {
    obj <- .build(char)$obj
    expect_true(is.numeric(obj$caps) && length(obj$caps) == 2L)
    expect_setequal(names(obj$caps), c("lo", "hi"))
    expect_true(all(is.finite(obj$caps)))
    expect_lt(obj$caps[["lo"]], obj$caps[["hi"]])
  }
})

test_that("caps match Matlab one-level xi support formulas", {
  ## Caps are on the standardized xi/v_hat scale consumed by deconvolution, not
  ## raw theta_hat. Matlab estimate_lsqnonlin.m constructs xi_hat = v_hat before
  ## applying the 5/7 sigma cap rules.
  br <- .build("race")
  expect_equal(unname(br$obj$caps), unname(.matlab_xi_caps(br$fit, "race")),
               tolerance = 1e-12)
  expect_equal(br$obj$caps[["lo"]], 0)            # one-sided floor, by design

  bg <- .build("gender")
  expect_equal(unname(bg$obj$caps), unname(.matlab_xi_caps(bg$fit, "gender")),
               tolerance = 1e-12)
  expect_equal(bg$obj$caps[["lo"]], -bg$obj$caps[["hi"]], tolerance = 1e-10)
  expect_lt(bg$obj$caps[["lo"]], bg$obj$caps[["hi"]])
})

test_that("v_hat and s_v are length N (matching the input)", {
  for (char in c("race", "gender")) {
    b <- .build(char)
    N <- length(b$input$theta_hat)
    expect_length(b$obj$v_hat, N)
    expect_length(b$obj$s_v, N)
  }
})

test_that("invariant #3: m_hat / V_m / report recomputed at fit$delta are identical", {
  ## The coupling object must NOT refit: every moment slot is a pass-through or
  ## a deterministic function of the ONE fit's delta (+ the fit's C, baked into
  ## V_m). Recompute from fit$delta and the rebuilt data matrix; require
  ## bit-identical equality with the object's slots.
  for (char in c("race", "gender")) {
    b <- .build(char)
    obj <- b$obj; fit <- b$fit
    data <- .datamat(b$input)

    ## m_hat + V_m are pure pass-throughs of the fit.
    expect_identical(obj$m_hat, fit$m_hat)
    expect_identical(obj$V_m, fit$V_m)
    expect_identical(obj$v_hat, fit$v_hat)
    expect_identical(obj$s_v, fit$s_v)

    ## report is gp_get_moments(extra=1) at fit$delta -- recompute and compare.
    m_full <- .moments(fit$delta, char, data = data, 0L, 1L)
    k <- length(m_full)
    expect_identical(obj$report$E_theta, unname(m_full[k - 1L]))
    expect_identical(obj$report$SD_theta, unname(m_full[k]))

    ## p_value is a deterministic function of the carried-through J.
    expect_identical(obj$p_value,
                     stats::pchisq(fit$J_stat, fit$df, lower.tail = FALSE))
  }
})

test_that("gp_w_seam one-call wrapper equals the two-step path on stable slots", {
  for (char in c("race", "gender")) {
    b <- .build(char)
    w <- .wseam(b$input, char)
    expect_s3_class(w, "gp_estimation_fit")
    expect_equal(w$beta, b$obj$beta)
    expect_equal(w$J, b$obj$J)
    expect_equal(unname(w$caps), unname(b$obj$caps))
    expect_equal(w$report$E_theta, b$obj$report$E_theta)
    expect_equal(w$report$SD_theta, b$obj$report$SD_theta)
  }
})

test_that("validator rejects malformed objects", {
  good <- .build("race")$obj
  ctor <- .gpc_get("new_gp_estimation_fit")

  ## df != 1
  bad_df <- good; bad_df$df <- 2
  expect_error(.validate(bad_df), class = "gradepath_validation_error")

  ## p_value out of range
  bad_p <- good; bad_p$p_value <- 1.5
  expect_error(.validate(bad_p), class = "gradepath_validation_error")

  ## report not a named list
  bad_rep <- good; bad_rep$report <- c(E_theta = 1, SD_theta = 2)
  expect_error(.validate(bad_rep), class = "gradepath_validation_error")

  ## caps lo >= hi
  bad_caps <- good; bad_caps$caps <- c(lo = 1, hi = 1)
  expect_error(.validate(bad_caps), class = "gradepath_validation_error")

  ## missing slot
  bad_miss <- good; bad_miss$beta <- NULL
  expect_error(.validate(bad_miss), class = "gradepath_validation_error")
})
