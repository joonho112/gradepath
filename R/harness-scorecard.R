# =============================================================================
# harness-scorecard.R -- scorecard and target validation wrappers
# =============================================================================

.gp_default_producer_status <- function(producer_status, fallback = "OK") {
  if (is.null(producer_status)) {
    return(fallback)
  }
  .gp_status_normalize(producer_status)
}

.gp_fit_producer_status <- function(fit, producer_status = NULL) {
  if (!is.null(producer_status)) {
    return(.gp_status_normalize(producer_status))
  }
  # Fail-SAFE: a `gp_fit` SHOULD always carry `provenance$producer_status`. If it
  # is missing we cannot assert the producing stage was acceptance-ready, so we
  # default to "UNVERIFIED" (which routes targets away from PASS) rather than the
  # fail-open "OK". An explicit producer_status override (above) still wins.
  .gp_status_normalize(fit$provenance$producer_status %gp_or% "UNVERIFIED")
}

.gp_replicated_frame <- function(replicated, producer_status = NULL) {
  producer_status <- .gp_default_producer_status(producer_status)

  if (is.data.frame(replicated)) {
    miss <- setdiff(c("id", "replicated"), names(replicated))
    if (length(miss) > 0L) {
      .gradepath_abort("`replicated` is missing column(s): %s.", paste(miss, collapse = ", "))
    }
    out <- data.frame(
      id = as.character(replicated$id),
      replicated = replicated$replicated,
      producer_status = if ("producer_status" %in% names(replicated)) {
        as.character(replicated$producer_status)
      } else {
        rep(producer_status, nrow(replicated))
      },
      group = if ("group" %in% names(replicated)) {
        as.character(replicated$group)
      } else {
        rep("Targets", nrow(replicated))
      },
      stringsAsFactors = FALSE
    )
  } else if (is.atomic(replicated) || is.list(replicated)) {
    ids <- names(replicated)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids))) {
      .gradepath_abort("Named replicated values must have non-empty names.")
    }
    out <- data.frame(
      id = as.character(ids),
      replicated = unlist(replicated, use.names = FALSE),
      producer_status = rep(producer_status, length(ids)),
      group = rep("Targets", length(ids)),
      stringsAsFactors = FALSE
    )
  } else {
    .gradepath_abort(
      "`replicated` must be a data frame, named atomic vector, or named list."
    )
  }
  if (anyDuplicated(out$id)) {
    .gradepath_abort("`replicated$id` values must be unique.")
  }
  out
}

.gp_fit_targets <- function(fit, producer_status = NULL) {
  fit <- validate_gp_fit(fit)
  producer_status <- .gp_fit_producer_status(fit, producer_status)
  demographic <- fit$provenance$demographic %gp_or% NULL
  if (is.null(demographic) || !demographic %in% c("race", "gender")) {
    .gradepath_abort(
      "`gp_validate_targets()` with a `gp_fit` requires provenance$demographic in c('race', 'gender')."
    )
  }
  prefix <- if (identical(demographic, "race")) "race_baseline" else "gender_baseline"
  selected_summary <- fit$selected_grade$summary
  selected_frontier <- fit$grade_path$summary[
    abs(fit$grade_path$summary$lambda - fit$grade_path$selection$selected_lambda) < 1e-8,
    ,
    drop = FALSE
  ]
  if (nrow(selected_frontier) != 1L) {
    .gradepath_abort("`gp_fit$grade_path` must have exactly one selected summary row.")
  }
  grades <- fit$selected_grade$assignment$grade
  values <- list(
    stats::setNames(length(fit$ids), "scale_n_firms_graded"),
    stats::setNames(selected_summary$grade_count, paste0(prefix, "_ngrades")),
    stats::setNames(100 * selected_frontier$discordance_rate, paste0(prefix, "_dr")),
    stats::setNames(selected_frontier$tau_bar, paste0(prefix, "_tau"))
  )
  if (identical(demographic, "race")) {
    values <- c(
      values,
      list(
        race_baseline_worst_n = sum(grades == min(grades)),
        race_baseline_best_n = sum(grades == max(grades))
      )
    )
  }
  values <- unlist(values, use.names = TRUE)
  .gp_replicated_frame(values, producer_status = producer_status)
}

.gp_target_rows <- function(replicated = NULL,
                            targets = NULL,
                            registry = NULL,
                            producer_status = NULL) {
  registry <- .gp_registry_resolve(registry)
  if (inherits(replicated, "gp_fit")) {
    frame <- .gp_fit_targets(replicated, producer_status = producer_status)
  } else if (is.null(replicated)) {
    if (is.null(targets)) {
      targets <- registry$id[registry$milestone == "M1"]
    }
    targets <- as.character(targets)
    frame <- data.frame(
      id = targets,
      replicated = NA_real_,
      producer_status = rep("UNVERIFIED", length(targets)),
      group = rep("Targets", length(targets)),
      stringsAsFactors = FALSE
    )
  } else {
    frame <- .gp_replicated_frame(replicated, producer_status = producer_status)
  }

  if (!is.null(targets)) {
    targets <- as.character(targets)
    frame <- frame[frame$id %in% targets, , drop = FALSE]
    missing <- setdiff(targets, frame$id)
    if (length(missing) > 0L) {
      frame <- rbind(
        frame,
        data.frame(
          id = missing,
          replicated = NA_real_,
          producer_status = rep("UNVERIFIED", length(missing)),
          group = rep("Targets", length(missing)),
          stringsAsFactors = FALSE
        )
      )
    }
    frame <- frame[match(targets, frame$id), , drop = FALSE]
  }
  frame
}

.gp_scorecard_summary <- function(checks) {
  scored <- checks$status %in% c("PASS", "FAIL")
  overall <- data.frame(
    group = "TOTAL",
    checks = nrow(checks),
    pass = sum(checks$status == "PASS"),
    fail = sum(checks$status == "FAIL"),
    unverified = sum(checks$status == "UNVERIFIED"),
    n_a = sum(checks$status == "n/a"),
    no_tol = sum(checks$status == "no-tol"),
    pass_rate = if (any(scored)) {
      sum(checks$status == "PASS") / sum(scored)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  groups <- unique(checks$group)
  rows <- lapply(groups, function(group) {
    subset <- checks[checks$group == group, , drop = FALSE]
    scored <- subset$status %in% c("PASS", "FAIL")
    pass <- sum(subset$status == "PASS")
    fail <- sum(subset$status == "FAIL")
    data.frame(
      group = group,
      checks = nrow(subset),
      pass = pass,
      fail = fail,
      unverified = sum(subset$status == "UNVERIFIED"),
      n_a = sum(subset$status == "n/a"),
      no_tol = sum(subset$status == "no-tol"),
      pass_rate = if (any(scored)) pass / sum(scored) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  rbind(overall, do.call(rbind, rows))
}

.gp_scorecard_provenance <- function() {
  .gradepath_provenance_stamp("gp_run_all")
}

#' Run a registry-driven gradepath scorecard
#'
#' `gp_run_all()` is the vector driver over [gp_check()]: it applies the registry
#' comparator to a set of replicated values (or to a `gp_fit`, whose own targets it
#' extracts) and returns a `gp_scorecard` with the per-target table plus a grouped
#' pass/fail summary. Passing `replicated = NULL` returns the requested targets as
#' `UNVERIFIED` (a coverage skeleton with nothing solved). [gp_validate_targets()] is
#' the thin verification-step alias.
#'
#' @param replicated What to score. A `gp_fit` (its registry targets -- firm count,
#'   grade count, discordance rate, tau, and for race the worst/best-grade counts --
#'   are derived and scored); a data frame with `id` and `replicated` columns
#'   (optionally `producer_status`, `group`); a named atomic vector or named list of
#'   `id = value`; or `NULL` (default) for an `UNVERIFIED` skeleton over `targets`.
#' @param targets Optional character vector of registry ids to score (and the row
#'   order). Defaults to `NULL`, which means all `M1` registry rows when `replicated`
#'   is `NULL`, otherwise the ids present in `replicated`. Requested ids absent from
#'   `replicated` are added as `UNVERIFIED` rows.
#' @param registry Registry data frame. Defaults to `NULL`, which uses the bundled
#'   package data [gp_registry].
#' @param producer_status Optional character producer-status override. Defaults to
#'   `NULL`: a `gp_fit` input then uses `fit$provenance$producer_status` when present
#'   and the fail-safe `"UNVERIFIED"` when absent; other replicated inputs (named
#'   vectors / lists / data frames) default to `"OK"`, the caller-asserts-OK contract
#'   for explicit replicated values.
#'
#' @return A list of class `c("gp_scorecard", "list")`. \describe{
#'   \item{`table`}{Per-target verdict data frame (the bound [gp_check()] rows with a
#'     leading `group` column); `checks` is a synonym slot with the same content.}
#'   \item{`summary`}{Per-`group` and `TOTAL` tally of `pass` / `fail` /
#'     `unverified` / `n_a` / `no_tol` counts and a `pass_rate`.}
#'   \item{`pass_rate`}{Numeric overall PASS fraction among scored (`PASS` + `FAIL`)
#'     rows, or `NA` when none are scored.}
#'   \item{`provenance`, `warnings`}{Internal audit slots.}
#' }
#'
#' @details
#' This is a verification scorecard against published KRW (2024) values; it makes no
#' ranking-superiority statement of any kind. Bundled artifacts record the
#' current M1 `NOT_ACCEPTED` gate: unresolved or non-`OK` producer rows stay
#' `UNVERIFIED` rather than counting as `PASS`.
#'
#' @examples
#' # Instant, no solve: score explicit replicated values against the registry.
#' # Named vector of id = value; the firm-count target is published as 97.
#' sc <- gp_run_all(c(scale_n_firms_graded = 97))
#' sc$table[, c("id", "paper", "replicated", "status")]
#' sc$pass_rate
#'
#' # NULL replicated returns the requested target as an UNVERIFIED skeleton.
#' gp_run_all(NULL, targets = "scale_n_firms_graded")$table$status
#'
#' \donttest{
#' # A gp_fit input derives and scores the fit's own registry targets.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' gp_run_all(fit)$summary
#' }
#'
#' @seealso [gp_validate_targets()], [gp_check()], [gp_registry]
#' @family gradepath-harness
#' @export
gp_run_all <- function(replicated = NULL,
                       targets = NULL,
                       registry = NULL,
                       producer_status = NULL) {
  registry <- .gp_registry_resolve(registry)
  frame <- .gp_target_rows(
    replicated = replicated,
    targets = targets,
    registry = registry,
    producer_status = producer_status
  )
  rows <- vector("list", nrow(frame))
  for (i in seq_len(nrow(frame))) {
    rows[[i]] <- gp_check(
      frame$id[[i]],
      replicated = frame$replicated[[i]],
      producer_status = frame$producer_status[[i]],
      registry = registry
    )
    rows[[i]]$group <- frame$group[[i]]
  }
  checks <- do.call(rbind, rows)
  checks <- checks[, c(
    "group",
    "id",
    "quantity",
    "paper",
    "replicated",
    "delta",
    "tol",
    "unit",
    "class",
    "milestone",
    "status",
    "reason",
    "producer_status"
  )]
  scored <- checks$status %in% c("PASS", "FAIL")
  pass_rate <- if (any(scored)) {
    mean(checks$status[scored] == "PASS")
  } else {
    NA_real_
  }
  structure(
    list(
      checks = checks,
      table = checks,
      summary = .gp_scorecard_summary(checks),
      pass_rate = pass_rate,
      provenance = .gp_scorecard_provenance(),
      warnings = character(0)
    ),
    class = c("gp_scorecard", "list")
  )
}

#' Validate replicated target values against the gradepath registry
#'
#' `gp_validate_targets()` is the explicit verification step that replaces v1's
#' replication-mode control flag: call it with the values you replicated to confirm
#' they match the published KRW (2024) registry targets within tolerance. It is a
#' thin convenience wrapper over [gp_run_all()] and returns the same `gp_scorecard`.
#'
#' @inheritParams gp_run_all
#'
#' @return A `gp_scorecard` (see [gp_run_all()] for the slot layout): `table`,
#'   `summary`, `pass_rate`, and internal `provenance` / `warnings`.
#'
#' @details
#' A verification step against published values; it makes no ranking-superiority
#' statement of any kind.
#'
#' @examples
#' # Instant, no solve: validate explicit replicated values against the registry.
#' sc <- gp_validate_targets(c(scale_n_firms_graded = 97))
#' sc$table[, c("id", "paper", "replicated", "status")]
#'
#' @seealso [gp_run_all()], [gp_check()], [gp_registry]
#' @family gradepath-harness
#' @export
gp_validate_targets <- function(replicated = NULL,
                                targets = NULL,
                                registry = NULL,
                                producer_status = NULL) {
  gp_run_all(
    replicated = replicated,
    targets = targets,
    registry = registry,
    producer_status = producer_status
  )
}
