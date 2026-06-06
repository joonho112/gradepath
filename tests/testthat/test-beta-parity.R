## KRW Table-3 parity for the W-seam coupling object `gp_estimation_fit`, on
## KRW's real GMM input.
##
## The built `gp_registry` contains only ONE non-interactive parity row for these
## quantities: t3_race_ni_beta = 0.510. There is NO gender no-industry beta row,
## and NO E_theta / SD_theta rows (verified live). So we:
##   * GATE race beta against the registry row t3_race_ni_beta;
##   * GATE gender beta against the published KRW (2024) Table 3 literal 1.255;
##   * PIN the realized report moments + caps as regression locks, recomputed
##     from the ONE raw core fit (invariant #3);
##   * FLAG the missing registry rows via skip() as a maintainer open item.
##
## NOTE: the coupling object `gp_estimation_fit` has NO `delta` slot, and the RAW
## core fit's `$provenance` has no `delta_hat`. Moment recomputes use the RAW core
## fit's numeric `$delta` (length-3), obtained from gp_estimation_core directly.

.gpc_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}
.krw     <- function(...) .gpc_get("gp_krw_gmm_input")(...)
.core    <- function(...) .gpc_get("gp_estimation_core")(...)
.couple  <- function(...) .gpc_get("gp_estimation_coupling")(...)
.datamat <- function(...) .gpc_get(".gp_estimation_data_matrix")(...)
.moments <- function(...) .gpc_get("gp_get_moments")(...)

## Raw core fit (has numeric $delta) + coupling object, once per characteristic.
.build <- function(char) {
  inp <- .krw(char)
  fit <- .core(inp, characteristic = char)
  list(input = inp, fit = fit, obj = .couple(fit, inp))
}

## Registry lookup by id (data.frame keyed by `id`, value in `paper_value`).
.reg_id <- function(id) {
  reg <- tryCatch(get("gp_registry", envir = asNamespace("gradepath")),
                  error = function(e) NULL)
  if (is.null(reg)) {
    reg <- tryCatch({
      e <- new.env()
      utils::data("gp_registry", package = "gradepath", envir = e)
      get("gp_registry", envir = e)
    }, error = function(e) NULL)
  }
  if (!is.data.frame(reg) || !all(c("id", "paper_value") %in% names(reg))) {
    return(NULL)
  }
  i <- which(reg$id == id)
  if (length(i) != 1L) return(NULL)
  list(value = as.numeric(reg$paper_value[i]),
       tol   = if ("tolerance" %in% names(reg))
                 suppressWarnings(as.numeric(reg$tolerance[i])) else NA_real_)
}

test_that("race beta matches the registry row t3_race_ni_beta (0.510)", {
  r <- .reg_id("t3_race_ni_beta")
  expect_false(is.null(r))
  expect_equal(r$value, 0.510)
  tol <- if (is.finite(r$tol) && r$tol > 0) r$tol else 1e-3
  expect_equal(.build("race")$obj$beta, r$value, tolerance = max(tol, 1e-3))
})

test_that("gender beta reproduces KRW (2024) Table 3 = 1.255 (no registry row)", {
  expect_null(.reg_id("gender_ni_beta"))
  expect_equal(.build("gender")$obj$beta, 1.255, tolerance = 0.01)
})

test_that("race report moments recompute from the ONE fit and lock to realized", {
  b <- .build("race")
  data <- .datamat(b$input)
  m_full <- .moments(b$fit$delta, "race", data, 0L, 1L)
  k <- length(m_full)
  expect_equal(b$obj$report$E_theta,  unname(m_full[k - 1L]))
  expect_equal(b$obj$report$SD_theta, unname(m_full[k]))
  expect_equal(b$obj$report$E_theta,  0.092358, tolerance = 1e-4)
  expect_equal(b$obj$report$SD_theta, 0.071515, tolerance = 1e-4)
})

test_that("gender report moments recompute from the ONE fit and lock to realized", {
  b <- .build("gender")
  data <- .datamat(b$input)
  m_full <- .moments(b$fit$delta, "gender", data, 0L, 1L)
  k <- length(m_full)
  expect_equal(b$obj$report$E_theta,  unname(m_full[k - 1L]))
  expect_equal(b$obj$report$SD_theta, unname(m_full[k]))
  expect_equal(b$obj$report$E_theta,  -0.008814, tolerance = 1e-4)
  expect_equal(b$obj$report$SD_theta,  0.180479, tolerance = 1e-4)
})

test_that("support caps are finite, ordered, and race is floored at 0", {
  br <- .build("race")
  expect_equal(br$obj$caps[["lo"]], 0)
  expect_gte(br$obj$caps[["hi"]], max(br$input$theta_hat))
  expect_lt(br$obj$caps[["lo"]], br$obj$caps[["hi"]])
  bg <- .build("gender")
  expect_true(all(is.finite(bg$obj$caps)))
  expect_lt(bg$obj$caps[["lo"]], bg$obj$caps[["hi"]])
})

test_that("missing registry rows (gender beta, E_theta/SD_theta) are an open item", {
  expect_null(.reg_id("race_ni_E_theta"))
  expect_null(.reg_id("race_ni_sd_theta"))
  expect_null(.reg_id("gender_ni_E_theta"))
  expect_null(.reg_id("gender_ni_sd_theta"))
  expect_null(.reg_id("gender_ni_beta"))
  skip(paste0(
    "OPEN ITEM (maintainer): gp_registry has only t3_race_ni_beta for these ",
    "quantities. Realized: race E_theta=0.0924 SD=0.0715, gender ",
    "E_theta=-0.0088 SD=0.1805. Add rows to data-raw/make-registry.R to ",
    "promote the locks above to registry gates."))
})
