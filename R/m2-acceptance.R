# =============================================================================
# m2-acceptance.R -- M2 acceptance gate and promotion ledger
# -----------------------------------------------------------------------------
# M2 remains separate from the M1 acceptance cache.  This file records the
# effective two-level registry after applying the fixture-promotion
# decision and assembles an explicit M2 scorecard.
# =============================================================================

.gp_m2_abort <- function(msg, ..., class = "gradepath_error") {
  args <- list(...)
  if (length(args) > 0L) msg <- do.call(sprintf, c(list(msg), args))
  .gradepath_abort(msg, class = class)
}

.gp_m2_l02_ids <- function() {
  c(
    "race_industry_dr",
    "race_industry_tau",
    "race_industry_r2",
    "gender_industry_dr",
    "gender_industry_tau",
    "gender_industry_r2"
  )
}

.gp_m2_l01_ids <- function() {
  c("race_industry_ngrades", "gender_industry_ngrades")
}

.gp_m2_n10_ids <- function() {
  c(
    "n10_race_mult_ngrades",
    "n10_race_mult_support",
    "n10_gender_add_ngrades",
    "n10_gender_add_support"
  )
}

.gp_m2_panel_for_id <- function(id) {
  id <- as.character(id)
  out <- rep(NA_character_, length(id))
  out[grepl("^race_", id) | grepl("_race_", id)] <- "race"
  out[grepl("^gender_", id) | grepl("_gender_", id)] <- "gender"
  out
}

# Recorded fixture snapshot.  These constants are frozen
# measurements from the full fixture gate, not derived at runtime. This is the
# default for the public M2 scorecard helpers; pass live
# gp_twolevel_fixture_gate objects through `fixture_gates` when the slow fixture
# backstop is being recomputed. The numeric gaps are pinned by the M2 acceptance
# tests and slow live backstop.
.gp_m2_recorded_fixture_gates <- function() {
  data.frame(
    characteristic = c("race", "gender"),
    pass = c(FALSE, TRUE),
    class_decision = c("approximate", "banded_candidate"),
    producer_status = c("APPROXIMATE_OK", "OK"),
    reason = c("Pi_theta outside 1pp fixture band",
               "all fixture checks passed"),
    pi_theta_max_abs = c(0.0121756787, 0.0066335820),
    posteriors_max_abs = c(0.0021779164, 0.0020897030),
    g_theta_support_max_abs = c(0.0001578607, 0.0000888980),
    g_theta_density_max_abs = c(0.0015353847, 0.0007618801),
    stringsAsFactors = FALSE
  )
}

.gp_m2_gate_frame <- function(fixture_gates = NULL) {
  if (is.null(fixture_gates)) {
    return(.gp_m2_recorded_fixture_gates())
  }
  if (inherits(fixture_gates, "gp_twolevel_fixture_gate")) {
    fixture_gates <- list(fixture_gates)
  }
  if (is.list(fixture_gates) && !is.data.frame(fixture_gates)) {
    rows <- lapply(fixture_gates, function(gate) {
      gate <- validate_gp_twolevel_fixture_gate(gate)
      checks <- gate$checks
      get_gap <- function(artifact) {
        vals <- checks$max_abs[checks$artifact == artifact]
        if (length(vals) == 0L) NA_real_ else vals[[1L]]
      }
      data.frame(
        characteristic = gate$characteristic,
        pass = gate$pass,
        class_decision = gate$class_decision,
        producer_status = gate$producer_status,
        reason = gate$reason,
        pi_theta_max_abs = get_gap("Pi_theta"),
        posteriors_max_abs = get_gap("posteriors"),
        g_theta_support_max_abs = get_gap("g_theta_support"),
        g_theta_density_max_abs = get_gap("g_theta_density"),
        stringsAsFactors = FALSE
      )
    })
    fixture_gates <- do.call(rbind, rows)
  }
  fixture_gates <- as.data.frame(fixture_gates, stringsAsFactors = FALSE)
  required <- c("characteristic", "pass", "class_decision", "producer_status")
  miss <- setdiff(required, names(fixture_gates))
  if (length(miss) > 0L) {
    .gp_m2_abort("`fixture_gates` is missing column(s): %s.",
                 paste(miss, collapse = ", "),
                 class = "gradepath_validation_error")
  }
  fixture_gates$characteristic <- match.arg(
    as.character(fixture_gates$characteristic),
    c("race", "gender"),
    several.ok = TRUE
  )
  if (anyDuplicated(fixture_gates$characteristic)) {
    .gp_m2_abort("`fixture_gates$characteristic` must be unique.",
                 class = "gradepath_validation_error")
  }
  fixture_gates$pass <- as.logical(fixture_gates$pass)
  if (any(is.na(fixture_gates$pass))) {
    .gp_m2_abort("`fixture_gates$pass` must be logical.",
                 class = "gradepath_validation_error")
  }
  if (!"reason" %in% names(fixture_gates)) {
    fixture_gates$reason <- ifelse(fixture_gates$pass,
                                   "all_fixture_checks_passed",
                                   "not_promoted")
  }
  for (nm in c("pi_theta_max_abs", "posteriors_max_abs",
               "g_theta_support_max_abs", "g_theta_density_max_abs")) {
    if (!nm %in% names(fixture_gates)) fixture_gates[[nm]] <- NA_real_
    fixture_gates[[nm]] <- suppressWarnings(as.numeric(fixture_gates[[nm]]))
  }
  fixture_gates
}

.gp_m2_effective_class <- function(gate_row) {
  promoted <- isTRUE(gate_row$pass) &&
    identical(as.character(gate_row$class_decision), "banded_candidate") &&
    identical(.gp_status_normalize(gate_row$producer_status), "OK")
  if (promoted) "banded" else "approximate"
}

#' Effective M2 registry after fixture-gate promotion
#'
#' `gp_m2_promoted_registry()` returns a copy of the package registry whose M2 L02
#' continuous-quantity `class` column is promoted from `approximate` to `banded`
#' only for panels whose fixture gate passes. The registry keeps these
#' rows `approximate` by default; this helper applies the panel-by-panel promotion
#' that [gp_m2_acceptance()] consumes.
#'
#' With the default `fixture_gates = NULL` it uses the recorded
#' race/gender fixture snapshot: that snapshot leaves the race L02 rows
#' `approximate` (the race `Pi_theta` gap is `0.0121756787 > 0.01`) and promotes the
#' gender L02 rows to `banded` (the gender fixture gate passed). Supply live
#' fixture-gate evidence to recompute the promotion from fresh fixtures.
#'
#' `PROMOTED` / `banded` means only that a panel's fixture artifacts met
#' the intermediate parity band used for M2 promotion. It does NOT mean the paper's
#' industry DR, tau, or R2 targets were directly reproduced or accepted.
#'
#' @param fixture_gates Optional fixture-gate evidence. When `NULL`
#'   (default) the frozen race/gender fixture snapshot is used
#'   instead of rerunning the fixture gate. May also be a live
#'   `gp_twolevel_fixture_gate` object, a list of such objects, or an equivalent
#'   data frame.
#' @param registry Optional registry data frame to promote. When `NULL` (default)
#'   the package [gp_registry] is used.
#'
#' @return A registry data frame (a copy of the input) with the eligible L02
#'   `class` entries promoted to `banded`, and two attributes attached:
#'   `m2_fixture_gates` (the race/gender gate data frame used) and `m2_l02_ids`
#'   (the L02 row ids considered for promotion).
#'
#' @examples
#' # Instant; no solve. The effective registry after fixture promotion.
#' reg <- gp_m2_promoted_registry()
#' attr(reg, "m2_l02_ids")
#'
#' @seealso [gp_m2_acceptance()], [gp_m2_status()], [gp_registry]
#' @family gradepath-twolevel
#' @export
gp_m2_promoted_registry <- function(fixture_gates = NULL, registry = NULL) {
  registry <- .gp_registry_resolve(registry)
  gates <- .gp_m2_gate_frame(fixture_gates)
  out <- registry
  l02 <- .gp_m2_l02_ids()
  for (panel in c("race", "gender")) {
    row <- gates[gates$characteristic == panel, , drop = FALSE]
    if (nrow(row) == 0L) next
    effective <- .gp_m2_effective_class(row[1L, , drop = FALSE])
    ids <- l02[.gp_m2_panel_for_id(l02) == panel]
    out$class[out$id %in% ids] <- effective
  }
  attr(out, "m2_fixture_gates") <- gates
  attr(out, "m2_l02_ids") <- l02
  out
}

.gp_m2_score_row <- function(gate, layer, target, id, group, status,
                             reason = NA_character_, notes = "",
                             registry_class = NA_character_,
                             effective_class = NA_character_,
                             producer_status = "OK",
                             paper = NA_real_,
                             replicated = NA_real_,
                             delta = NA_real_,
                             tol = NA_real_,
                             unit = NA_character_,
                             fixture_decision = NA_character_,
                             source = "gp_m2_acceptance") {
  data.frame(
    gate = gate,
    layer = layer,
    target = target,
    id = id,
    group = group,
    status = status,
    reason = reason,
    producer_status = producer_status,
    registry_class = registry_class,
    effective_class = effective_class,
    paper = as.numeric(paper),
    replicated = as.numeric(replicated),
    delta = as.numeric(delta),
    tol = as.numeric(tol),
    unit = unit,
    fixture_decision = fixture_decision,
    source = source,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

.gp_m2_registered_row <- function(id, replicated, group, target,
                                  registry, notes = "") {
  check <- gp_check(id, replicated = replicated, registry = registry)
  .gp_m2_score_row(
    gate = "L01",
    layer = "registered_target",
    target = target,
    id = id,
    group = group,
    status = check$status,
    reason = check$reason,
    registry_class = check$class,
    effective_class = check$class,
    producer_status = check$producer_status,
    paper = check$paper,
    replicated = check$replicated,
    delta = check$delta,
    tol = check$tol,
    unit = check$unit,
    fixture_decision = "not_applicable",
    source = "gp_check",
    notes = notes
  )
}

.gp_m2_n10_row <- function(id, replicated, status, reason, notes, registry) {
  row <- .gp_registry_row(id, registry = registry)
  .gp_m2_score_row(
    gate = "N10",
    layer = "synthetic_guard",
    target = id,
    id = id,
    group = "N10 support parity",
    status = status,
    reason = reason,
    registry_class = as.character(row$class[[1L]]),
    effective_class = as.character(row$class[[1L]]),
    producer_status = "OK",
    paper = suppressWarnings(as.numeric(row$paper_value[[1L]])),
    replicated = replicated,
    delta = NA_real_,
    tol = suppressWarnings(as.numeric(row$tolerance[[1L]])),
    unit = as.character(row$unit[[1L]]),
    fixture_decision = "synthetic_guard",
    source = "test-n10-support-parity.R",
    notes = notes
  )
}

.gp_m2_l02_rows <- function(effective_registry, gates, registry) {
  rows <- list()
  for (id in .gp_m2_l02_ids()) {
    panel <- .gp_m2_panel_for_id(id)
    gate <- gates[gates$characteristic == panel, , drop = FALSE]
    row <- .gp_registry_row(id, registry = registry)
    eff <- .gp_registry_row(id, registry = effective_registry)
    promoted <- identical(as.character(eff$class[[1L]]), "banded")
    rows[[length(rows) + 1L]] <- .gp_m2_score_row(
      gate = "L02",
      layer = "promotion_gate",
      target = id,
      id = id,
      group = sprintf("%s industry continuous", panel),
      status = if (promoted) "PROMOTED" else "APPROXIMATE_OK",
      reason = if (promoted) "fixture_gate_passed" else "fixture_gate_not_promoted",
      registry_class = as.character(row$class[[1L]]),
      effective_class = as.character(eff$class[[1L]]),
      producer_status = if (nrow(gate) == 1L) {
        .gp_status_normalize(gate$producer_status[[1L]])
      } else {
        "UNVERIFIED"
      },
      paper = suppressWarnings(as.numeric(row$paper_value[[1L]])),
      replicated = NA_real_,
      delta = NA_real_,
      tol = suppressWarnings(as.numeric(row$tolerance[[1L]])),
      unit = as.character(row$unit[[1L]]),
      fixture_decision = if (nrow(gate) == 1L) {
        gate$class_decision[[1L]]
      } else {
        "fixture_gate_missing"
      },
      source = "gp_twolevel_fixture_gate",
      notes = if (nrow(gate) == 1L) {
        sprintf(
          "%s gate: class_decision=%s; reason=%s; Pi_theta gap=%s.",
          panel,
          gate$class_decision[[1L]],
          gate$reason[[1L]],
          format(gate$pi_theta_max_abs[[1L]], digits = 10)
        )
      } else {
        "No fixture gate supplied for this panel."
      }
    )
  }
  do.call(rbind, rows)
}

.gp_m2_score_summary <- function(table) {
  data.frame(
    rows = nrow(table),
    pass = sum(table$status == "PASS"),
    promoted = sum(table$status == "PROMOTED"),
    approximate_ok = sum(table$status == "APPROXIMATE_OK"),
    evidence_ok = sum(table$status == "EVIDENCE_OK"),
    fail = sum(table$status == "FAIL"),
    stringsAsFactors = FALSE
  )
}

#' Build the M2 acceptance scorecard
#'
#' `gp_m2_acceptance()` assembles the M2 L01 / L02 / N10 ledger and the overall
#' formal M2 status as a `gp_m2_acceptance` object. With the default
#' `fixture_gates = NULL` it uses the recorded fixture snapshot
#' rather than rerunning the slow fixture gate; pass live fixture-gate evidence to
#' rebuild the scorecard. For a one-screen summary use [gp_m2_status()].
#'
#' The default recorded scorecard is `PARTIAL_ACCEPTED`: L01 industry grade-count
#' rows pass exactly; N10 rows are synthetic support-guard evidence rather than
#' registered paper-value passes; gender L02 continuous rows are `PROMOTED` to
#' `banded` by fixture parity; and race L02 continuous rows remain `APPROXIMATE_OK`
#' because the recorded race `Pi_theta` fixture gap is `0.0121756787 > 0.01`.
#'
#' `PROMOTED` / `banded` means only that a panel's fixture artifacts met
#' the intermediate parity band used for M2 promotion. It does NOT mean the paper's
#' industry DR, tau, or R2 targets were directly reproduced or accepted; some
#' continuous L02 rows are fixture-parity evidence, not direct reproduction.
#'
#' @param fixture_gates Optional fixture-gate evidence. When `NULL`
#'   (default) the frozen race/gender fixture snapshot is used
#'   instead of rerunning the slow fixture gate. May also be a live
#'   `gp_twolevel_fixture_gate` object, a list of such objects, or an equivalent
#'   data frame, to recompute the promotion decision from fresh evidence.
#' @param registry Optional registry data frame supplying the L01 / L02 / N10
#'   target rows. When `NULL` (default) the package [gp_registry] is used.
#'
#' @return A validated `gp_m2_acceptance` object (a list of class
#'   `c("gp_m2_acceptance", "list")`) with the public slots: \describe{
#'   \item{`table`}{Data frame; the full M2 ledger, one row per gate
#'     (`M2_ACCEPTANCE`, `L01`, `L02`, `N10`, `OVERRIDE`) with `status`, `reason`,
#'     `producer_status`, registry/effective class, paper/replicated values, and
#'     notes. The `M2_ACCEPTANCE` row carries the overall formal status
#'     (`PARTIAL_ACCEPTED` by default).}
#'   \item{`checks`}{An alias of `table`.}
#'   \item{`effective_registry`}{Registry data frame after panel-by-panel L02
#'     promotion (see [gp_m2_promoted_registry()]).}
#'   \item{`fixture_gates`}{The race/gender fixture-gate data frame used.}
#'   \item{`summary`}{One-row data frame counting rows by status (`pass`,
#'     `promoted`, `approximate_ok`, `evidence_ok`, `fail`).}
#'   \item{`provenance`, `warnings`}{Internal audit slots.}
#' }
#'
#' @examples
#' # Instant; no solve. The default recorded M2 scorecard.
#' acc <- gp_m2_acceptance()
#' acc$summary
#'
#' @seealso [gp_m2_status()], [gp_m2_promoted_registry()], [gp_twolevel_grade()],
#'   [gp_check()]
#' @family gradepath-twolevel
#' @export
gp_m2_acceptance <- function(fixture_gates = NULL, registry = NULL) {
  registry <- .gp_registry_resolve(registry)
  gates <- .gp_m2_gate_frame(fixture_gates)
  effective_registry <- gp_m2_promoted_registry(gates, registry = registry)

  rows <- list(
    .gp_m2_registered_row(
      "race_industry_ngrades",
      replicated = 4,
      group = "Industry grade counts",
      target = "race_industry_grade_count",
      registry = registry,
      notes = "L01 exact race industry grade count: 4."
    ),
    .gp_m2_registered_row(
      "gender_industry_ngrades",
      replicated = 5,
      group = "Industry grade counts",
      target = "gender_industry_grade_count",
      registry = registry,
      notes = "L01 exact gender industry grade count: 5."
    ),
    .gp_m2_score_row(
      gate = "OVERRIDE",
      layer = "component_evidence",
      target = "same_industry_override",
      id = "same_industry_override",
      group = "Two-level pairwise",
      status = "EVIDENCE_OK",
      reason = NA_character_,
      registry_class = "n/a",
      effective_class = "n/a",
      producer_status = "OK",
      unit = "evidence",
      fixture_decision = "same_industry_override_test",
      source = "test-same-industry-override.R",
      notes = "Tests pin L_ij = L_i for same-industry pairs and route all five Pi matrices through the same denominator."
    )
  )
  rows <- c(
    rows,
    list(
      .gp_m2_n10_row(
        "n10_race_mult_ngrades",
        replicated = 4,
        "EVIDENCE_OK",
        "synthetic_exact_count",
        "N10 race multiplicative fixture returns 4 grades exactly.",
        registry
      ),
      .gp_m2_n10_row(
        "n10_race_mult_support",
        replicated = 0,
        "EVIDENCE_OK",
        "synthetic_banded_support",
        "N10 race support endpoint matches analytic range within 1e-12.",
        registry
      ),
      .gp_m2_n10_row(
        "n10_gender_add_ngrades",
        replicated = 5,
        "EVIDENCE_OK",
        "synthetic_exact_count",
        "N10 gender additive fixture returns 5 grades exactly.",
        registry
      ),
      .gp_m2_n10_row(
        "n10_gender_add_support",
        replicated = 0,
        "EVIDENCE_OK",
        "synthetic_banded_support",
        "N10 gender support endpoint matches analytic range within 1e-12.",
        registry
      )
    )
  )
  rows[[length(rows) + 1L]] <- .gp_m2_l02_rows(effective_registry, gates, registry)
  table <- do.call(rbind, rows)

  has_fail <- any(table$status == "FAIL")
  has_approx <- any(table$status == "APPROXIMATE_OK")
  formal_status <- if (has_fail) {
    "NOT_ACCEPTED"
  } else if (has_approx) {
    "PARTIAL_ACCEPTED"
  } else {
    "ACCEPTED"
  }
  formal <- .gp_m2_score_row(
    gate = "M2_ACCEPTANCE",
    layer = "formal_gate",
    target = "overall",
    id = "m2_acceptance",
    group = "M2 acceptance",
    status = formal_status,
    reason = if (has_fail) {
      "M2_FAIL"
    } else if (has_approx) {
      "ONE_OR_MORE_L02_ROWS_REMAIN_APPROXIMATE"
    } else {
      "ALL_M2_ROWS_PROMOTED_OR_PASS"
    },
    registry_class = "mixed",
    effective_class = "mixed",
    producer_status = if (has_fail) "UNVERIFIED" else "OK",
    unit = "gate",
    fixture_decision = "panel_specific",
    source = "gp_m2_acceptance",
    notes = "Exact L01/N10 gates pass; L02 continuous rows promote panel-by-panel from fixture evidence."
  )
  table <- rbind(formal, table)

  out <- structure(
    list(
      table = table,
      checks = table,
      effective_registry = effective_registry,
      fixture_gates = gates,
      summary = .gp_m2_score_summary(table),
      provenance = .gradepath_new_provenance(
        producer = "gp_m2_acceptance",
        step = "8.7",
        promotion_rule = "panel_l02_banded_iff_step8.5_fixture_gate_passes"
      ),
      warnings = character(0)
    ),
    class = c("gp_m2_acceptance", "list")
  )
  validate_gp_m2_acceptance(out)
}

validate_gp_m2_acceptance <- function(x) {
  if (!inherits(x, "gp_m2_acceptance")) {
    .gp_m2_abort("Expected a gp_m2_acceptance object.",
                 class = "gradepath_validation_error")
  }
  req <- c("table", "checks", "effective_registry", "fixture_gates",
           "summary", "provenance", "warnings")
  if (any(!req %in% names(x))) {
    .gp_m2_abort("`gp_m2_acceptance` is missing required fields.",
                 class = "gradepath_validation_error")
  }
  table <- as.data.frame(x$table, stringsAsFactors = FALSE)
  required_cols <- c("gate", "layer", "target", "id", "group", "status",
                     "reason", "producer_status", "registry_class",
                     "effective_class", "paper", "replicated", "delta",
                     "tol", "unit", "fixture_decision", "source", "notes")
  if (any(!required_cols %in% names(table))) {
    .gp_m2_abort("M2 scorecard table is missing required columns.",
                 class = "gradepath_validation_error")
  }
  allowed_status <- c("ACCEPTED", "PARTIAL_ACCEPTED", "NOT_ACCEPTED",
                      "PASS", "FAIL", "PROMOTED", "APPROXIMATE_OK",
                      "EVIDENCE_OK", "DEFERRED")
  if (any(!table$status %in% allowed_status)) {
    .gp_m2_abort("M2 scorecard has an unknown status.",
                 class = "gradepath_validation_error")
  }
  effective_registry <- .gp_registry_resolve(x$effective_registry)
  gates <- .gp_m2_gate_frame(x$fixture_gates)
  summary <- as.data.frame(x$summary, stringsAsFactors = FALSE)
  if (!identical(x$checks, x$table)) {
    .gp_m2_abort("`checks` must alias `table`.",
                 class = "gradepath_validation_error")
  }
  .gradepath_validate_named_list(x$provenance, "provenance")
  .gradepath_validate_warning_vector(x$warnings, "warnings")
  x$table <- table
  x$effective_registry <- effective_registry
  x$fixture_gates <- gates
  x$summary <- summary
  x
}
