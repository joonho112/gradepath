# =============================================================================
# m2-status.R -- Public M2-only status surface
# -----------------------------------------------------------------------------
# Keeps M2 status reporting out of R/status.R, which is part of the M1 cache
# manifest.
# =============================================================================

.gp_m2_status_row <- function(component, status, detail, evidence,
                              producer_status = "OK",
                              source = "gp_m2_acceptance") {
  data.frame(
    component = component,
    status = status,
    detail = detail,
    evidence = evidence,
    producer_status = producer_status,
    source = source,
    stringsAsFactors = FALSE
  )
}

.gp_m2_l02_status <- function(rows) {
  if (nrow(rows) > 0L && all(rows$status == "PROMOTED")) {
    "PROMOTED"
  } else if (nrow(rows) > 0L && all(rows$status == "APPROXIMATE_OK")) {
    "APPROXIMATE_OK"
  } else {
    "UNVERIFIED"
  }
}

.gp_m2_l02_detail <- function(panel, status, gate) {
  panel_title <- paste0(toupper(substr(panel, 1L, 1L)), substr(panel, 2L, nchar(panel)))
  if (identical(status, "PROMOTED")) {
    sprintf(
      "%s L02 promoted to banded: Pi_theta gap %.10f <= 0.01.",
      panel_title,
      gate$pi_theta_max_abs[[1L]]
    )
  } else if (identical(status, "APPROXIMATE_OK")) {
    sprintf(
      "%s L02 remains approximate: Pi_theta gap %.10f > 0.01.",
      panel_title,
      gate$pi_theta_max_abs[[1L]]
    )
  } else {
    sprintf("%s L02 status is not verified.", panel_title)
  }
}

#' Summarize the M2 acceptance status
#'
#' `gp_m2_status()` returns a compact public summary of the current M2 acceptance
#' ledger -- one row per status component -- without touching the M1 status
#' contract in `R/status.R`. It is the one-screen read of [gp_m2_acceptance()] and
#' the verb [gp_twolevel_report_card()] points users to while the M2 surface is
#' `PARTIAL_ACCEPTED`.
#'
#' The default summary is derived from [gp_m2_acceptance()]. It reports the formal
#' M2 status, L01 exact industry grade-count evidence, L02 race/gender
#' fixture-promotion status, and N10 synthetic support-guard evidence. The default
#' M2 state remains `PARTIAL_ACCEPTED`: race L02 continuous rows remain approximate
#' because the recorded `Pi_theta` fixture gap is above `0.01`, while gender L02
#' rows are promoted to `banded` by fixture parity.
#'
#' This helper is M2-only. It does not use `new_gp_status()`, does not edit
#' `R/status.R`, and does not recalculate the M1 acceptance state; the returned
#' object's `m1_status_recalculated` attribute is always `FALSE`.
#'
#' @param acceptance Optional `gp_m2_acceptance` object to summarize. When `NULL`
#'   (default) a fresh [gp_m2_acceptance()] (the recorded-snapshot scorecard) is
#'   built and summarized.
#'
#' @return A `gp_m2_status` object: a data frame of class
#'   `c("gp_m2_status", "data.frame")` with one row per component
#'   (`formal`, `L01_industry_grade_counts`, `L02_race_continuous`,
#'   `L02_gender_continuous`, `N10_support_guards`, `M1_status_boundary`) and the
#'   columns `component`, `status`, `detail`, `evidence`, `producer_status`,
#'   `source`. It carries the attributes `m2_acceptance` (the backing
#'   `gp_m2_acceptance`), `m2_formal_status` (`"PARTIAL_ACCEPTED"` by default), and
#'   `m1_status_recalculated` (always `FALSE`).
#'
#' @examples
#' # Instant; no solve. The per-component M2 ledger summary.
#' gp_m2_status()
#'
#' # The headline formal status is also an attribute.
#' attr(gp_m2_status(), "m2_formal_status")   # "PARTIAL_ACCEPTED"
#'
#' @seealso [gp_m2_acceptance()], [gp_m2_promoted_registry()],
#'   [gp_twolevel_report_card()]
#' @family gradepath-twolevel
#' @export
gp_m2_status <- function(acceptance = NULL) {
  acceptance <- if (is.null(acceptance)) {
    gp_m2_acceptance()
  } else {
    validate_gp_m2_acceptance(acceptance)
  }
  table <- acceptance$table
  gates <- acceptance$fixture_gates

  formal <- table[table$id == "m2_acceptance", , drop = FALSE]
  if (nrow(formal) != 1L) {
    .gp_m2_abort("M2 scorecard must have exactly one formal status row.",
                 class = "gradepath_validation_error")
  }

  l01 <- table[table$gate == "L01", , drop = FALSE]
  n10 <- table[table$gate == "N10", , drop = FALSE]
  race_l02 <- table[table$gate == "L02" & grepl("^race_", table$id), ,
                    drop = FALSE]
  gender_l02 <- table[table$gate == "L02" & grepl("^gender_", table$id), ,
                      drop = FALSE]
  race_gate <- gates[gates$characteristic == "race", , drop = FALSE]
  gender_gate <- gates[gates$characteristic == "gender", , drop = FALSE]
  if (nrow(race_gate) != 1L || nrow(gender_gate) != 1L) {
    .gp_m2_abort("M2 fixture gates must have one race row and one gender row.",
                 class = "gradepath_validation_error")
  }

  l01_status <- if (nrow(l01) > 0L && all(l01$status == "PASS")) {
    "PASS"
  } else {
    "FAIL"
  }
  n10_status <- if (nrow(n10) > 0L && all(n10$status == "EVIDENCE_OK")) {
    "EVIDENCE_OK"
  } else {
    "UNVERIFIED"
  }
  race_l02_status <- .gp_m2_l02_status(race_l02)
  gender_l02_status <- .gp_m2_l02_status(gender_l02)

  rows <- list(
    .gp_m2_status_row(
      "formal",
      formal$status[[1L]],
      formal$reason[[1L]],
      formal$notes[[1L]],
      formal$producer_status[[1L]]
    ),
    .gp_m2_status_row(
      "L01_industry_grade_counts",
      l01_status,
      sprintf("%d/%d exact count rows pass.",
              sum(l01$status == "PASS"), nrow(l01)),
      "Race and gender industry grade-count targets are registered L01 exact rows."
    ),
    .gp_m2_status_row(
      "L02_race_continuous",
      race_l02_status,
      .gp_m2_l02_detail("race", race_l02_status, race_gate),
      "Race industry DR/tau/R2 rows are classified by fixture parity, not direct paper-value reproduction.",
      race_gate$producer_status[[1L]],
      "gp_twolevel_fixture_gate"
    ),
    .gp_m2_status_row(
      "L02_gender_continuous",
      gender_l02_status,
      .gp_m2_l02_detail("gender", gender_l02_status, gender_gate),
      "Gender industry DR/tau/R2 rows are classified by fixture parity, not direct paper-value reproduction.",
      gender_gate$producer_status[[1L]],
      "gp_twolevel_fixture_gate"
    ),
    .gp_m2_status_row(
      "N10_support_guards",
      n10_status,
      sprintf("%d/%d N10 rows carry synthetic evidence status.",
              sum(n10$status == "EVIDENCE_OK"), nrow(n10)),
      "N10 support rows are synthetic guard evidence, not registered paper-value passes.",
      source = "test-n10-support-parity.R"
    ),
    .gp_m2_status_row(
      "M1_status_boundary",
      "NOT_RECALCULATED",
      "M1 acceptance status is not recalculated by gp_m2_status().",
      "This helper summarizes M2 only and does not edit or recompute R/status.R, the M1 cache, or M1 acceptance scorecards.",
      producer_status = "not_applicable",
      source = "gp_m2_status"
    )
  )

  out <- do.call(rbind, rows)
  class(out) <- c("gp_m2_status", "data.frame")
  attr(out, "m2_acceptance") <- acceptance
  attr(out, "m2_formal_status") <- formal$status[[1L]]
  attr(out, "m1_status_recalculated") <- FALSE
  out
}
