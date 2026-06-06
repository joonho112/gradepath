#' W-seam coupling object: the typed `gp_estimation_fit`
#'
#' @description
#' The W-seam coupling stage of `gradepath`. The estimation core (`gp_estimation_core()`, in
#' `R/estimation-core.R`) produces the raw two-step efficient beta-GMM fit as a
#' plain named list. This file wraps that raw fit into the typed, validated S3
#' object `gp_estimation_fit` -- the W-seam's ONLY downstream contract. Every
#' field is read off the ONE fit; nothing is re-estimated (frozen invariant #3).
#' The over-identification statistic `J` is carried through verbatim from the
#' core (single-N already applied there) and is never rescaled (frozen invariant
#' #7); we only attach its chi-square p-value.
#'
#' The deconvolution reads the slots of this object by
#' NAME. The slot names are therefore the contract and must not drift:
#'
#'   - `beta`           precision exponent (scalar)
#'   - `m_hat`          GMM-implied coupling moments (extra_moments = 0):
#'                      race `c(mu, sigma_xi)`, gender `c(sigma_xi)`
#'   - `V_m`            delta-method covariance of `m_hat` (the fit's `V_m`)
#'   - `v_hat`, `s_v`   standardized-residual carriers (length N)
#'   - `J`              the over-id statistic (= `fit$J_stat`; named `J`)
#'   - `df`             degrees of freedom (1)
#'   - `p_value`        `pchisq(J, df, lower.tail = FALSE)`
#'   - `report`         named list `list(E_theta = ., SD_theta = .)`
#'                      (paper reporting only; NO covariance)
#'   - `caps`           support caps for the deconvolution prior search,
#'                      named numeric `c(lo = ., hi = .)`
#'   - `characteristic` "race" / "gender"
#'   - `provenance`     one-fit record plus reporting metadata such as `mu`
#'
#' There is deliberately NO `m_report`, NO `V_report`, NO `beta_se`, and NO
#' `reject` slot. Reporting moments carry no covariance.
#'
#' RELATIONSHIP TO `gp_precision_fit`: see NOTES.md. In brief,
#' `gp_estimation_fit` is the *concrete W-seam estimation output* (a populated,
#' GMM-specific contract consumed by the deconvolution), whereas the
#' `gp_precision_fit` is a generic precision-stage shell. They are distinct S3
#' classes; this file does not subclass or mutate `gp_precision_fit`.
#'
#' @section Invariants enforced here:
#' #3 (coupling): `m_hat`, `V_m`, `report`, and `caps` are ALL derived from the
#'   single `fit$delta` (and, for `V_m`, the single sandwich `fit$C` already
#'   baked into `fit$V_m`). Nothing is refit. #7 (single-N): `J <- fit$J_stat`
#'   verbatim; only `p_value <- pchisq(J, df, lower.tail = FALSE)` is added.
#'
#' @keywords internal
#' @noRd
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Resolve the package abort helper whether sourced or namespaced
#'
#' During development this file is `source()`-d into the global environment on
#' top of a `devtools::load_all()`-ed namespace. `.gradepath_abort` then lives in
#' the package namespace, not the global env. This shim finds it either way and
#' falls back to a plain `stop()` with a classed condition if neither is present.
#'
#' @keywords internal
#' @noRd
.gp_coupling_abort <- function(msg, class = "gradepath_error") {
  fn <- tryCatch(get(".gradepath_abort", envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (is.function(fn)) return(fn(msg, class = class))
  if (exists(".gradepath_abort", mode = "function")) {
    return(get(".gradepath_abort")(msg, class = class))
  }
  cnd <- structure(
    class = c(class, "gradepath_error", "error", "condition"),
    list(message = msg, call = NULL)
  )
  stop(cnd)
}

#' Resolve a package internal by name (namespace-first), erroring if absent
#' @keywords internal
#' @noRd
.gp_coupling_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  .gp_coupling_abort(
    sprintf("Required gradepath internal '%s' not found.", name),
    class = "gradepath_internal_error"
  )
}

#' Support caps for the deconvolution prior search
#'
#' @description
#' Port of the one-level xi-support block of `estimate_lsqnonlin.m`. Caps are
#' computed on the deconvolution carrier scale (`xi_hat == v_hat`), not on raw
#' `theta_hat`. Race is one-sided and floored at 0. Gender is centered at zero
#' on the xi scale and clipped by the 5/7 sigma rule. These are the BOUNDS only;
#' the grid resolution (`supp_pts`) is a deconvolution concern.
#'
#' @param m_hat One-level coupling moments: race `c(mu_xi, sigma_xi)`, gender
#'   `c(sigma_xi)`.
#' @param characteristic "race" or "gender".
#' @param v_hat Numeric vector of standardized carriers (`xi_hat` in Matlab).
#' @return A named numeric `c(lo = ., hi = .)`.
#' @keywords internal
#' @noRd
.gp_coupling_caps <- function(m_hat, characteristic, v_hat) {
  if (!is.numeric(v_hat) || length(v_hat) < 1L || any(!is.finite(v_hat))) {
    .gp_coupling_abort(
      "`v_hat` must be a non-empty finite numeric vector.",
      class = "gradepath_validation_error"
    )
  }
  if (identical(characteristic, "race")) {
    if (!is.numeric(m_hat) || length(m_hat) != 2L || any(!is.finite(m_hat))) {
      .gp_coupling_abort(
        "`m_hat` must be c(mu_xi, sigma_xi) for race.",
        class = "gradepath_validation_error"
      )
    }
    mu_xi <- m_hat[1L]
    sigma_xi <- m_hat[2L]
    lo <- 0
    hi <- min(max(max(v_hat), mu_xi + 5 * sigma_xi), mu_xi + 7 * sigma_xi)
  } else if (identical(characteristic, "gender")) {
    if (!is.numeric(m_hat) || length(m_hat) != 1L || any(!is.finite(m_hat))) {
      .gp_coupling_abort(
        "`m_hat` must be c(sigma_xi) for gender.",
        class = "gradepath_validation_error"
      )
    }
    sigma_xi <- m_hat[1L]
    lo <- max(min(min(v_hat), -5 * sigma_xi), -7 * sigma_xi)
    hi <- min(max(max(v_hat),  5 * sigma_xi),  7 * sigma_xi)
  } else {
    .gp_coupling_abort(
      sprintf("Unknown `characteristic`: '%s'.", characteristic),
      class = "gradepath_validation_error"
    )
  }
  c(lo = unname(lo), hi = unname(hi))
}

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Bare constructor for `gp_estimation_fit`
#'
#' @description
#' Low-level: assembles the list and stamps the class. Performs NO validation
#' (use `validate_gp_estimation_fit()`). All slots are taken as given.
#'
#' @param beta,m_hat,V_m,v_hat,s_v,J,df,p_value,report,caps,characteristic,provenance
#'   The contract slots; see the file header.
#' @return An object of class `"gp_estimation_fit"`.
#' @keywords internal
#' @noRd
new_gp_estimation_fit <- function(beta, m_hat, V_m, v_hat, s_v, J, df,
                                  p_value, report, caps, characteristic,
                                  provenance = NULL) {
  structure(
    list(
      beta           = beta,
      m_hat          = m_hat,
      V_m            = V_m,
      v_hat          = v_hat,
      s_v            = s_v,
      J              = J,
      df             = df,
      p_value        = p_value,
      report         = report,
      caps           = caps,
      characteristic = characteristic,
      provenance     = provenance
    ),
    class = "gp_estimation_fit"
  )
}

# ---------------------------------------------------------------------------
# Validator
# ---------------------------------------------------------------------------

#' Validate a `gp_estimation_fit`
#'
#' @description
#' Asserts the full contract: slot presence; shapes consistent with the
#' characteristic (race `m_hat` length 2, gender length 1); `V_m` square,
#' symmetric, and positive semidefinite (PSD); `df == 1`; `J` finite and `>= 0`;
#' `p_value` in `[0, 1]` and finite; `report` a named list with names exactly
#' `c("E_theta", "SD_theta")` and finite values; `caps` finite with
#' `caps["lo"] < caps["hi"]`; `v_hat` and `s_v` equal length. Returns `x`
#' invisibly on success; aborts (classed `gradepath_validation_error`) otherwise.
#'
#' @param x A candidate `gp_estimation_fit`.
#' @return `x`, invisibly.
#' @keywords internal
#' @noRd
validate_gp_estimation_fit <- function(x) {
  bad <- function(msg) {
    .gp_coupling_abort(msg, class = "gradepath_validation_error")
  }

  if (!inherits(x, "gp_estimation_fit")) {
    bad("Object is not of class 'gp_estimation_fit'.")
  }

  required <- c("beta", "m_hat", "V_m", "v_hat", "s_v", "J", "df",
                "p_value", "report", "caps", "characteristic", "provenance")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    bad(sprintf("gp_estimation_fit is missing slot(s): %s.",
                paste(missing, collapse = ", ")))
  }

  ## characteristic
  if (!is.character(x$characteristic) || length(x$characteristic) != 1L ||
      !x$characteristic %in% c("race", "gender")) {
    bad("`characteristic` must be a single string in c('race','gender').")
  }
  expected_m <- if (identical(x$characteristic, "race")) 2L else 1L

  ## beta
  if (!is.numeric(x$beta) || length(x$beta) != 1L || !is.finite(x$beta)) {
    bad("`beta` must be a finite numeric scalar.")
  }

  ## m_hat shape + finiteness
  if (!is.numeric(x$m_hat) || length(x$m_hat) != expected_m ||
      any(!is.finite(x$m_hat))) {
    bad(sprintf("`m_hat` must be a finite numeric vector of length %d for %s.",
                expected_m, x$characteristic))
  }

  ## V_m: square (expected_m x expected_m), symmetric, PSD, finite
  Vm <- x$V_m
  if (!is.matrix(Vm) || nrow(Vm) != expected_m || ncol(Vm) != expected_m ||
      any(!is.finite(Vm))) {
    bad(sprintf("`V_m` must be a finite %d x %d matrix for %s.",
                expected_m, expected_m, x$characteristic))
  }
  if (max(abs(Vm - t(Vm))) > 1e-8 * max(1, max(abs(Vm)))) {
    bad("`V_m` must be symmetric.")
  }
  ev <- tryCatch(eigen(Vm, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || any(ev < -1e-8 * max(1, max(abs(Vm))))) {
    bad("`V_m` must be positive semidefinite.")
  }

  ## v_hat / s_v equal length
  if (!is.numeric(x$v_hat) || !is.numeric(x$s_v)) {
    bad("`v_hat` and `s_v` must be numeric vectors.")
  }
  if (length(x$v_hat) != length(x$s_v)) {
    bad(sprintf("`v_hat` (%d) and `s_v` (%d) must be the same length.",
                length(x$v_hat), length(x$s_v)))
  }
  if (length(x$v_hat) < 1L) {
    bad("`v_hat`/`s_v` must be non-empty (length N).")
  }

  ## df == 1
  if (!is.numeric(x$df) || length(x$df) != 1L || x$df != 1) {
    bad("`df` must equal 1.")
  }

  ## J finite >= 0
  if (!is.numeric(x$J) || length(x$J) != 1L || !is.finite(x$J) || x$J < 0) {
    bad("`J` must be a finite numeric scalar >= 0.")
  }

  ## p_value in [0, 1], finite
  if (!is.numeric(x$p_value) || length(x$p_value) != 1L ||
      !is.finite(x$p_value) || x$p_value < 0 || x$p_value > 1) {
    bad("`p_value` must be a finite numeric scalar in [0, 1].")
  }

  ## report: named list c('E_theta','SD_theta'), finite
  rep <- x$report
  if (!is.list(rep) || is.null(names(rep)) ||
      !identical(sort(names(rep)), c("E_theta", "SD_theta"))) {
    bad("`report` must be a named list with names c('E_theta','SD_theta').")
  }
  if (!is.numeric(rep$E_theta) || length(rep$E_theta) != 1L ||
      !is.finite(rep$E_theta) ||
      !is.numeric(rep$SD_theta) || length(rep$SD_theta) != 1L ||
      !is.finite(rep$SD_theta)) {
    bad("`report$E_theta` and `report$SD_theta` must be finite numeric scalars.")
  }

  ## caps: finite, named lo/hi, lo < hi
  cp <- x$caps
  if (!is.numeric(cp) || length(cp) != 2L ||
      !all(c("lo", "hi") %in% names(cp)) || any(!is.finite(cp))) {
    bad("`caps` must be a finite named numeric c(lo = ., hi = .).")
  }
  if (!(cp[["lo"]] < cp[["hi"]])) {
    bad(sprintf("`caps` must satisfy lo < hi (got lo=%g, hi=%g).",
                cp[["lo"]], cp[["hi"]]))
  }

  invisible(x)
}

# ---------------------------------------------------------------------------
# The coupling builder: raw core fit + input -> typed gp_estimation_fit
# ---------------------------------------------------------------------------

#' Assemble the typed W-seam coupling object from a raw core fit
#'
#' @description
#' Takes the raw list returned by `gp_estimation_core()` plus the SAME `input`
#' that was passed to the core, and produces a validated `gp_estimation_fit`.
#' The `input` is needed only to rebuild the `c(d, theta_hat, s)` data matrix
#' (via `.gp_estimation_data_matrix()`), because the core fit does not carry
#' the data; the reporting moments (extra_moments = 1) and the support caps both
#' need the per-firm `(theta_hat, s)` values. `R/estimation-core.R` is left
#' untouched (see NOTES.md for the alternative of caching `data` in the core
#' return).
#'
#' Invariant #3: `m_hat`, `V_m`, `report`, and `caps` are all derived from the
#' single `fit$delta`; nothing is refit, and `V_m` is the fit's own delta-method
#' covariance. Invariant #7: `J <- fit$J_stat` verbatim; only the chi-square
#' p-value is attached.
#'
#' @param fit The raw list from `gp_estimation_core()`.
#' @param input The eb_estimates-like object originally passed to the core
#'   (carries `theta_hat`, `s`).
#' @param control Optional; accepted for signature symmetry with the core and
#'   the wrapper. Not read here (the coupling step has no tunables).
#' @return A validated object of class `"gp_estimation_fit"`.
#' @keywords internal
#' @noRd
gp_estimation_coupling <- function(fit, input, control = NULL) {
  ## --- argument sanity ----------------------------------------------------
  needed <- c("delta", "J_stat", "df", "m_hat", "V_m", "v_hat", "s_v",
              "characteristic", "C", "provenance")
  if (!is.list(fit) || any(!needed %in% names(fit))) {
    miss <- setdiff(needed, names(fit))
    .gp_coupling_abort(
      sprintf(paste0("`fit` does not look like a gp_estimation_core() result; ",
                     "missing: %s."), paste(miss, collapse = ", ")),
      class = "gradepath_validation_error"
    )
  }
  characteristic <- fit$characteristic
  if (!characteristic %in% c("race", "gender")) {
    .gp_coupling_abort(
      sprintf("`fit$characteristic` must be 'race' or 'gender' (got '%s').",
              characteristic),
      class = "gradepath_validation_error"
    )
  }

  ## --- rebuild the data matrix from the SAME input ------------------------
  dm_fn <- .gp_coupling_get(".gp_estimation_data_matrix")
  data <- dm_fn(input)
  theta_hat <- data[, 2L]

  ## Defensive: the data must be consistent with the fit's carriers (one fit).
  if (length(theta_hat) != length(fit$v_hat)) {
    .gp_coupling_abort(
      sprintf(paste0("`input` (N=%d) is inconsistent with the fit's carriers ",
                     "(N=%d): the same input passed to gp_estimation_core() ",
                     "must be passed here (invariant #3)."),
              length(theta_hat), length(fit$v_hat)),
      class = "gradepath_validation_error"
    )
  }

  ## --- report (paper moments) via extra_moments = 1, ONE fit's delta ------
  get_moments <- .gp_coupling_get("gp_get_moments")
  m_full <- get_moments(fit$delta, characteristic, data = data,
                        group_fx = 0L, extra_moments = 1L)
  k <- length(m_full)
  report <- list(E_theta = unname(m_full[k - 1L]),
                 SD_theta = unname(m_full[k]))

  ## --- support caps from the ONE fit's standardized carrier ----------------
  caps <- .gp_coupling_caps(fit$m_hat, characteristic, fit$v_hat)

  ## --- J carried through verbatim (invariant #7); only attach p-value -----
  J <- fit$J_stat
  df <- fit$df
  p_value <- stats::pchisq(J, df, lower.tail = FALSE)

  ## --- assemble + validate ------------------------------------------------
  provenance <- fit$provenance
  provenance$mu <- fit$mu
  provenance$mu_source <- if (identical(characteristic, "race")) {
    "exp(delta[1])"
  } else {
    "delta[1]"
  }

  obj <- new_gp_estimation_fit(
    beta           = fit$beta,
    m_hat          = fit$m_hat,
    V_m            = fit$V_m,
    v_hat          = fit$v_hat,
    s_v            = fit$s_v,
    J              = J,
    df             = df,
    p_value        = p_value,
    report         = report,
    caps           = caps,
    characteristic = characteristic,
    provenance     = provenance
  )
  validate_gp_estimation_fit(obj)
}

#' One-call W-seam: estimate then couple
#'
#' @description
#' Convenience wrapper that runs the estimation core and then builds the
#' typed coupling object in one call. Equivalent to
#' `gp_estimation_coupling(gp_estimation_core(input, characteristic = ...),
#' input, control)`. The same `input` is threaded to both stages so the data
#' matrix the coupling rebuilds is byte-for-byte the one the core estimated on.
#'
#' This is the natural user-facing entry point for the W-seam and may be
#' EXPORTED in a later step; for now it is internal (see NOTES.md).
#'
#' @param input The eb_estimates-like object (carries `theta_hat`, `s`).
#' @param characteristic "race" or "gender".
#' @param control Optional `gp_control`; passed through to the core.
#' @return A validated `gp_estimation_fit`.
#' @keywords internal
#' @noRd
gp_w_seam <- function(input, characteristic, control = NULL) {
  core_fn <- .gp_coupling_get("gp_estimation_core")
  fit <- core_fn(input, control = control, characteristic = characteristic)
  gp_estimation_coupling(fit, input, control = control)
}

# ---------------------------------------------------------------------------
# print method (optional convenience)
# ---------------------------------------------------------------------------

#' Print a `gp_estimation_fit`
#' @param x A `gp_estimation_fit`.
#' @param digits Number of significant digits.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
#' @noRd
print.gp_estimation_fit <- function(x, digits = 4L, ...) {
  cat("<gp_estimation_fit>  (W-seam coupling object)\n")
  cat(sprintf("  characteristic : %s\n", x$characteristic))
  cat(sprintf("  beta           : %s\n", format(x$beta, digits = digits)))
  cat(sprintf("  m_hat          : %s\n",
              paste(format(x$m_hat, digits = digits), collapse = ", ")))
  cat(sprintf("  J (df=%g)       : %s   p_value = %s\n",
              x$df, format(x$J, digits = digits),
              format(x$p_value, digits = digits)))
  cat(sprintf("  report         : E_theta = %s,  SD_theta = %s\n",
              format(x$report$E_theta, digits = digits),
              format(x$report$SD_theta, digits = digits)))
  cat(sprintf("  caps           : lo = %s,  hi = %s\n",
              format(x$caps[["lo"]], digits = digits),
              format(x$caps[["hi"]], digits = digits)))
  cat(sprintf("  carriers       : v_hat/s_v length %d\n", length(x$v_hat)))
  invisible(x)
}
