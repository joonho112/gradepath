# =============================================================================
# two-level-fixture-gate.R  --  quadrature-vs-fixture gate
# -----------------------------------------------------------------------------
# Compares quadrature artifacts against the archive group_fx == 1
# fixtures.  This is the M2 decision input: if deterministic quadrature reaches
# fixture parity for Pi_theta, posteriors, and g_theta within the L02 band, the
# continuous two-level targets are eligible for approximate -> banded promotion.
# =============================================================================

#' Resolve package abort helper
#' @keywords internal
#' @noRd
.gp_tlfg_abort <- function(msg, ..., class = "gradepath_error") {
  args <- list(...)
  if (length(args) > 0L) msg <- do.call(sprintf, c(list(msg), args))
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

#' Default L02 fixture bands
#' @keywords internal
#' @noRd
.gp_tlfg_default_tolerances <- function() {
  list(
    Pi_theta = 0.01,
    posteriors = 0.01,
    g_theta_support = 5e-4,
    g_theta_density = 5e-3
  )
}

#' Resolve fixture directory
#' @keywords internal
#' @noRd
.gp_tlfg_fixture_dir <- function(fixture_dir = NULL) {
  if (!is.null(fixture_dir)) return(fixture_dir)
  p <- system.file("extdata", "fixtures", package = "gradepath")
  if (nzchar(p) && dir.exists(p)) return(p)
  file.path("inst", "extdata", "fixtures")
}

#' Read a headerless numeric fixture if present
#' @keywords internal
#' @noRd
.gp_tlfg_read_fixture <- function(fixture_dir, file, expected_cols = NULL) {
  path <- file.path(fixture_dir, file)
  if (!file.exists(path)) {
    return(list(path = path, available = FALSE, value = NULL))
  }
  x <- utils::read.csv(path, header = FALSE)
  if (!is.null(expected_cols) && ncol(x) != expected_cols) {
    .gp_tlfg_abort("Fixture `%s` must have %d column(s).",
                   file, expected_cols, class = "gradepath_validation_error")
  }
  mat <- as.matrix(x)
  storage.mode(mat) <- "double"
  if (any(!is.finite(mat))) {
    .gp_tlfg_abort("Fixture `%s` must be finite.", file,
                   class = "gradepath_validation_error")
  }
  list(path = path, available = TRUE, value = mat)
}

#' One-row comparison result
#' @keywords internal
#' @noRd
.gp_tlfg_row <- function(artifact, status, max_abs = NA_real_,
                         tolerance = NA_real_, n = NA_integer_,
                         path = NA_character_, reason = "OK") {
  data.frame(
    artifact = artifact,
    status = status,
    max_abs = as.numeric(max_abs),
    tolerance = as.numeric(tolerance),
    n = as.integer(n),
    fixture_path = as.character(path),
    reason = as.character(reason),
    stringsAsFactors = FALSE
  )
}

#' Matrix comparison row
#' @keywords internal
#' @noRd
.gp_tlfg_compare_matrix <- function(artifact, actual, expected, tolerance, path) {
  actual <- as.matrix(actual)
  storage.mode(actual) <- "double"
  if (!identical(dim(actual), dim(expected))) {
    return(.gp_tlfg_row(
      artifact = artifact,
      status = "FAIL",
      tolerance = tolerance,
      n = length(actual),
      path = path,
      reason = sprintf("shape_mismatch_actual_%s_expected_%s",
                       paste(dim(actual), collapse = "x"),
                       paste(dim(expected), collapse = "x"))
    ))
  }
  if (any(!is.finite(actual))) {
    return(.gp_tlfg_row(artifact, "FAIL", tolerance = tolerance,
                        n = length(actual), path = path,
                        reason = "actual_nonfinite"))
  }
  gap <- max(abs(actual - expected))
  .gp_tlfg_row(
    artifact = artifact,
    status = if (gap <= tolerance + 1e-12) "PASS" else "FAIL",
    max_abs = gap,
    tolerance = tolerance,
    n = length(actual),
    path = path,
    reason = if (gap <= tolerance + 1e-12) "OK" else "outside_l02_band"
  )
}

#' Compare quadrature/artifact outputs against groupfx1 fixtures
#'
#' @description
#' The default path consumes a `gp_twolevel_quadrature` object.  Tests and cache
#' builders may also pass an `artifacts` list directly when they want to compare
#' `posteriors`/`g_theta` without materializing the full 97 x 97 Pi matrix.
#'
#' @return A validated `gp_twolevel_fixture_gate` object.
#' @keywords internal
#' @noRd
gp_twolevel_fixture_gate <- function(quadrature = NULL, artifacts = NULL,
                                     characteristic = NULL,
                                     fixture_dir = NULL,
                                     tolerances = NULL,
                                     require_pi = TRUE) {
  if (is.null(artifacts)) {
    if (is.null(quadrature)) {
      .gp_tlfg_abort("Either `quadrature` or `artifacts` is required.",
                     class = "gradepath_validation_error")
    }
    quadrature <- validate_gp_twolevel_quadrature(quadrature)
    artifacts <- quadrature$artifacts
    characteristic <- characteristic %gp_or%
      quadrature$posterior$metadata$two_level$characteristic
  }
  characteristic <- match.arg(as.character(characteristic), c("race", "gender"))
  fixture_dir <- .gp_tlfg_fixture_dir(fixture_dir)
  tol <- .gp_tlfg_default_tolerances()
  if (!is.null(tolerances)) {
    for (nm in names(tolerances)) tol[[nm]] <- tolerances[[nm]]
  }

  rows <- list()

  pi_file <- sprintf("Pi_groupfx1_%s.csv", characteristic)
  pi_fixture <- .gp_tlfg_read_fixture(fixture_dir, pi_file)
  if (isTRUE(require_pi)) {
    if (!isTRUE(pi_fixture$available)) {
      rows[[length(rows) + 1L]] <- .gp_tlfg_row(
        "Pi_theta", "SKIP", tolerance = tol$Pi_theta,
        path = pi_fixture$path, reason = "fixture_missing"
      )
    } else if (is.null(artifacts$Pi_theta)) {
      rows[[length(rows) + 1L]] <- .gp_tlfg_row(
        "Pi_theta", "FAIL", tolerance = tol$Pi_theta,
        path = pi_fixture$path, reason = "artifact_missing"
      )
    } else {
      rows[[length(rows) + 1L]] <- .gp_tlfg_compare_matrix(
        "Pi_theta", artifacts$Pi_theta, pi_fixture$value, tol$Pi_theta,
        pi_fixture$path
      )
    }
  } else {
    rows[[length(rows) + 1L]] <- .gp_tlfg_row(
      "Pi_theta", "SKIP", tolerance = tol$Pi_theta,
      path = pi_fixture$path, reason = "pi_check_disabled"
    )
  }

  post_file <- sprintf("posteriors_groupfx1_%s.csv", characteristic)
  post_fixture <- .gp_tlfg_read_fixture(fixture_dir, post_file, expected_cols = 4L)
  if (!isTRUE(post_fixture$available)) {
    rows[[length(rows) + 1L]] <- .gp_tlfg_row(
      "posteriors", "SKIP", tolerance = tol$posteriors,
      path = post_fixture$path, reason = "fixture_missing"
    )
  } else if (is.null(artifacts$posteriors)) {
    rows[[length(rows) + 1L]] <- .gp_tlfg_row(
      "posteriors", "FAIL", tolerance = tol$posteriors,
      path = post_fixture$path, reason = "artifact_missing"
    )
  } else {
    rows[[length(rows) + 1L]] <- .gp_tlfg_compare_matrix(
      "posteriors", artifacts$posteriors, post_fixture$value,
      tol$posteriors, post_fixture$path
    )
  }

  g_file <- sprintf("g_theta_groupfx1_%s.csv", characteristic)
  g_fixture <- .gp_tlfg_read_fixture(fixture_dir, g_file, expected_cols = 2L)
  if (!isTRUE(g_fixture$available)) {
    rows[[length(rows) + 1L]] <- .gp_tlfg_row(
      "g_theta_support", "SKIP", tolerance = tol$g_theta_support,
      path = g_fixture$path, reason = "fixture_missing"
    )
    rows[[length(rows) + 1L]] <- .gp_tlfg_row(
      "g_theta_density", "SKIP", tolerance = tol$g_theta_density,
      path = g_fixture$path, reason = "fixture_missing"
    )
  } else if (is.null(artifacts$g_theta)) {
    rows[[length(rows) + 1L]] <- .gp_tlfg_row(
      "g_theta_support", "FAIL", tolerance = tol$g_theta_support,
      path = g_fixture$path, reason = "artifact_missing"
    )
    rows[[length(rows) + 1L]] <- .gp_tlfg_row(
      "g_theta_density", "FAIL", tolerance = tol$g_theta_density,
      path = g_fixture$path, reason = "artifact_missing"
    )
  } else {
    rows[[length(rows) + 1L]] <- .gp_tlfg_compare_matrix(
      "g_theta_support", artifacts$g_theta$support,
      g_fixture$value[, 1L, drop = FALSE],
      tol$g_theta_support, g_fixture$path
    )
    rows[[length(rows) + 1L]] <- .gp_tlfg_compare_matrix(
      "g_theta_density", artifacts$g_theta$density,
      g_fixture$value[, 2L, drop = FALSE],
      tol$g_theta_density, g_fixture$path
    )
  }

  checks <- do.call(rbind, rows)
  all_pass <- all(checks$status == "PASS")
  any_fail <- any(checks$status == "FAIL")
  decision <- if (all_pass) "banded_candidate" else "approximate"
  producer_status <- if (all_pass) "OK" else "APPROXIMATE_OK"
  reason <- if (all_pass) {
    "all_fixture_checks_passed"
  } else if (any_fail) {
    "one_or_more_fixture_checks_failed"
  } else {
    "one_or_more_fixture_checks_skipped"
  }

  out <- list(
    characteristic = characteristic,
    checks = checks,
    pass = all_pass,
    class_decision = decision,
    producer_status = producer_status,
    reason = reason,
    tolerances = tol,
    fixture_dir = fixture_dir,
    require_pi = isTRUE(require_pi),
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(
      step = "two-level-fixture-gate",
      characteristic = characteristic,
      class_decision = decision,
      producer_status = producer_status,
      require_pi = isTRUE(require_pi)
    ),
    warnings = character(0)
  )
  validate_gp_twolevel_fixture_gate(
    structure(out, class = c("gp_twolevel_fixture_gate", "list"))
  )
}

#' Validate a fixture gate object
#' @keywords internal
#' @noRd
validate_gp_twolevel_fixture_gate <- function(x) {
  if (!inherits(x, "gp_twolevel_fixture_gate")) {
    .gp_tlfg_abort("Expected a gp_twolevel_fixture_gate object.",
                   class = "gradepath_validation_error")
  }
  req <- c("characteristic", "checks", "pass", "class_decision",
           "producer_status", "reason", "tolerances", "fixture_dir",
           "require_pi", "schema_version", "provenance", "warnings")
  if (any(!req %in% names(x))) {
    .gp_tlfg_abort("`gp_twolevel_fixture_gate` is missing required fields.",
                   class = "gradepath_validation_error")
  }
  if (!identical(x$characteristic, "race") &&
      !identical(x$characteristic, "gender")) {
    .gp_tlfg_abort("Fixture gate characteristic must be 'race' or 'gender'.",
                   class = "gradepath_validation_error")
  }
  checks <- x$checks
  needed_cols <- c("artifact", "status", "max_abs", "tolerance", "n",
                   "fixture_path", "reason")
  if (!is.data.frame(checks) || any(!needed_cols %in% names(checks)) ||
      nrow(checks) < 1L) {
    .gp_tlfg_abort("Fixture gate checks must be a non-empty data.frame.",
                   class = "gradepath_validation_error")
  }
  if (any(!checks$status %in% c("PASS", "FAIL", "SKIP"))) {
    .gp_tlfg_abort("Fixture gate status must be PASS, FAIL, or SKIP.",
                   class = "gradepath_validation_error")
  }
  expected_pass <- all(checks$status == "PASS")
  if (!identical(x$pass, expected_pass)) {
    .gp_tlfg_abort("Fixture gate pass flag is inconsistent with checks.",
                   class = "gradepath_validation_error")
  }
  if (expected_pass) {
    if (!identical(x$class_decision, "banded_candidate") ||
        !identical(x$producer_status, "OK")) {
      .gp_tlfg_abort("Passing fixture gate must be a banded OK candidate.",
                     class = "gradepath_validation_error")
    }
  } else if (!identical(x$class_decision, "approximate") ||
             !identical(x$producer_status, "APPROXIMATE_OK")) {
    .gp_tlfg_abort("Non-passing fixture gate must remain approximate.",
                   class = "gradepath_validation_error")
  }
  x
}
