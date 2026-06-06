# test-print-summary.R
# Result-first console surface (print / summary / format_gp_*_cli)
# for gp_pairwise, gp_grade_fit, gp_grade_path, gp_frontier, gp_report_card,
# gp_fit. GP-DEC-14-A; frozen invariant 2.
#
# Self-contained: local solver-free fixtures (ported from test-frontier.R /
# test-report-card.R). testthat edition 3.

# ---------------------------------------------------------------------------
# Local fixtures (copies of the constructor patterns; no solver involved)
# ---------------------------------------------------------------------------

ps_pairwise <- function(pairwise_matrix = NULL, ids = NULL,
                        control = gp_control(backend = "highs")) {
  if (is.null(pairwise_matrix)) pairwise_matrix <- matrix(c(0.5, 0.9, 0.1, 0.5), 2, byrow = TRUE)
  if (is.null(ids)) ids <- paste0("unit_", seq_len(nrow(pairwise_matrix)))
  pairwise_matrix <- as.matrix(pairwise_matrix)
  storage.mode(pairwise_matrix) <- "double"
  dimnames(pairwise_matrix) <- list(ids, ids)
  new_gp_pairwise(
    ids = ids, matrix = pairwise_matrix, power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = .gp_pairwise_diagonal,
                   zero_floor = .gp_pairwise_zero_floor),
    source = list(stage = "posterior", rule = "outer_product",
                  assumption = "one_level_independence"),
    control = control
  )
}

ps_grade_fit <- function(ids, lambda, grades, control) {
  new_gp_grade_fit(
    ids = ids, lambda = lambda,
    assignment = data.frame(id = ids, grade = as.integer(grades), stringsAsFactors = FALSE),
    summary = list(grade_count = as.integer(length(unique(grades))),
                   status = "optimal", n_units = length(ids)),
    objective = list(value = 0, raw = 0, canonical = 0),
    backend = list(name = control$backend, path = "fixture", encoding = "fixture",
                   status = "optimal", objbound = NA_real_, mipgap = NA_real_,
                   runtime = NA_real_, warmstart = "none",
                   warm_start_from_lambda = NA_real_, warm_start_used = FALSE,
                   problem_hash = "fixture"),
    control = control
  )
}

ps_grade_path <- function(pairwise,
                          fit_specs = list(list(lambda = 0.25, grades = c(1L, 1L)),
                                           list(lambda = 1, grades = c(1L, 2L))),
                          enriched = TRUE) {
  control <- pairwise$control
  fits <- lapply(fit_specs, function(spec) ps_grade_fit(pairwise$ids, spec$lambda, spec$grades, control))
  lambda_grid <- vapply(fit_specs, function(spec) spec$lambda, numeric(1))
  summary <- if (isTRUE(enriched)) frontier_table(pairwise, fits) else
    data.frame(lambda = lambda_grid,
               grade_count = vapply(fits, function(f) f$summary$grade_count, integer(1)))
  new_gp_grade_path(
    ids = pairwise$ids, lambda_grid = lambda_grid, fits = fits, summary = summary,
    backend = list(name = control$backend),
    selection = list(selected_lambda = 0.25, selection_rule = "baseline_lambda_0.25",
                     endpoint_lambda = 1),
    control = control
  )
}

ps_ids <- c("u1", "u2", "u3", "u4")
ps_ctl <- function() gp_control(backend = "highs")
ps_rc_grade_fit <- function(lambda, grades, ids = ps_ids, control = ps_ctl()) {
  new_gp_grade_fit(
    ids = ids, lambda = lambda,
    assignment = data.frame(id = ids, grade = as.integer(grades), stringsAsFactors = FALSE),
    summary = list(grade_count = as.integer(length(unique(grades)))),
    objective = list(value = 0, raw = 0, canonical = 0),
    backend = list(name = control$backend, path = "fixture"), control = control
  )
}
ps_rc_grade_path <- function(ids = ps_ids, control = ps_ctl()) {
  selected <- ps_rc_grade_fit(0.25, c(2L, 1L, 2L, 3L), ids, control)
  endpoint <- ps_rc_grade_fit(1, c(2L, 1L, 1L, 2L), ids, control)
  fits <- list(selected, endpoint); lambda_grid <- c(0.25, 1)
  new_gp_grade_path(
    ids = ids, lambda_grid = lambda_grid, fits = fits,
    summary = data.frame(lambda = lambda_grid,
                         grade_count = vapply(fits, function(f) f$summary$grade_count, integer(1))),
    backend = list(name = control$backend),
    selection = list(selected_lambda = 0.25, selection_rule = "baseline_lambda_0.25",
                     endpoint_lambda = 1),
    control = control
  )
}
ps_rc_posterior <- function(ids = ps_ids, labels = paste("Posterior", ids)) {
  pm <- c(0.4, 0.1, 0.2, 0.3)
  new_gp_posterior(estimate = c(4, 1, 2, 3), se = rep(0.2, length(ids)), id = ids,
                   label = labels, posterior_mean = pm, posterior_sd = rep(0.05, length(ids)),
                   lower = pm - 0.01, upper = pm + 0.01, scale = "r")
}
ps_rc_estimates <- function(ids = ps_ids, labels = paste("Label", ids)) {
  list(unit_id = ids, label = labels, theta_hat = c(10, 20, 30, 40), s = rep(0.4, length(ids)))
}
ps_rc_pairwise <- function(ids = ps_ids, control = ps_ctl()) {
  m <- matrix(0.5, length(ids), length(ids), dimnames = list(ids, ids))
  for (i in seq_along(ids)) for (j in seq_along(ids)) if (i != j) m[i, j] <- if (i < j) 0.6 else 0.4
  new_gp_pairwise(ids = ids, matrix = m, power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = .gp_pairwise_diagonal,
                   zero_floor = .gp_pairwise_zero_floor),
    source = list(stage = "posterior", rule = "outer_product",
                  assumption = "one_level_independence"),
    control = control)
}
ps_rc_prior <- function() {
  support <- seq(-1, 1, length.out = 11); density <- rep(1 / length(support), length(support))
  new_gp_prior(support = support, density = density, mean = sum(support * density), scale = "r")
}
ps_report_card <- function() {
  path <- ps_rc_grade_path(); posterior <- ps_rc_posterior()
  gp_report_card(ps_rc_estimates(), posterior = posterior,
                 selected_grade = path$fits[[1L]], grade_path = path)
}
ps_fit <- function() {
  path <- ps_rc_grade_path(); posterior <- ps_rc_posterior()
  card <- gp_report_card(ps_rc_estimates(), posterior = posterior,
                         selected_grade = path$fits[[1L]], grade_path = path)
  new_gp_fit(
    ids = ps_ids,
    estimates = ebrecipe::eb_input(theta_hat = c(10, 20, 30, 40),
                                   s = rep(0.4, length(ps_ids)), unit_id = ps_ids),
    prior = ps_rc_prior(), posterior = posterior, precision_fit = NULL,
    pairwise = ps_rc_pairwise(), grade_path = path, selected_grade = path$fits[[1L]],
    report_card = card, control = ps_ctl()
  )
}

ps_objects <- function() {
  pw <- ps_pairwise()
  path <- ps_grade_path(pw)
  list(
    gp_pairwise    = pw,
    gp_grade_fit   = path$fits[[1L]],
    gp_grade_path  = path,
    gp_frontier    = gp_frontier(path, pairwise = pw),
    gp_report_card = ps_report_card(),
    gp_fit         = ps_fit()
  )
}

ps_formatters <- list(
  gp_pairwise    = format_gp_pairwise_cli,
  gp_grade_fit   = format_gp_grade_fit_cli,
  gp_grade_path  = format_gp_grade_path_cli,
  gp_frontier    = format_gp_frontier_cli,
  gp_report_card = format_gp_report_card_cli,
  gp_fit         = format_gp_fit_cli
)

# ---------------------------------------------------------------------------
# Per-object contract: print outputs + invisibility + summary + formatter type
# ---------------------------------------------------------------------------

test_that("print() produces output, returns x invisibly, and recomputes nothing", {
  objs <- ps_objects()
  for (type in names(objs)) {
    x <- objs[[type]]
    expect_output(print(x), regexp = ".", info = type)
    expect_invisible(print(x))
    ret <- withVisible(print(x))
    expect_identical(ret$value, x)
    expect_false(ret$visible)
  }
})

test_that("format_gp_*_cli() returns a non-empty character vector", {
  objs <- ps_objects()
  for (type in names(objs)) {
    out <- ps_formatters[[type]](objs[[type]])
    expect_type(out, "character")
    expect_gt(length(out), 0L)
    expect_false(anyNA(out))
  }
})

test_that("formatter output equals what print() emits (print delegates)", {
  objs <- ps_objects()
  for (type in names(objs)) {
    x <- objs[[type]]
    printed <- capture.output(print(x))
    formatted <- ps_formatters[[type]](x)
    expect_identical(printed, formatted, info = type)
  }
})

# ---------------------------------------------------------------------------
# LEADS WITH THE RESULT (object-specific headline assertions)
# ---------------------------------------------------------------------------

test_that("each view leads with the result", {
  objs <- ps_objects()

  pw <- format_gp_pairwise_cli(objs$gp_pairwise)
  expect_true(any(grepl("2 x 2", pw)))
  expect_true(any(grepl("ordered pairs", pw)))
  expect_true(any(grepl("pi range", pw)))

  gf <- format_gp_grade_fit_cli(objs$gp_grade_fit)
  expect_true(any(grepl("grades", gf[[2L]])))
  expect_true(any(grepl("lambda = 0.25", gf)))
  expect_true(any(grepl("status = optimal", gf)))

  gp <- format_gp_grade_path_cli(objs$gp_grade_path)
  expect_true(any(grepl("2 lambda values", gp)))
  expect_true(any(grepl("selected lambda = 0.25", gp)))
  expect_true(any(grepl("grades across path", gp)))

  fr <- format_gp_frontier_cli(objs$gp_frontier)
  expect_true(any(grepl("selected lambda = 0.25", fr)))
  expect_true(any(grepl("1 - DR", fr)))
  expect_true(any(grepl("tau-bar", fr)))
  expect_true(any(grepl("benchmarks", fr)))

  rc <- format_gp_report_card_cli(objs$gp_report_card)
  expect_true(any(grepl("<gp_report_card>", rc)))
  expect_true(any(grepl("selected lambda: 0.25", rc)))
  expect_true(any(grepl("grades: 3 \\(1/2/1\\)", rc)))
  expect_true(any(grepl("posterior_mean", rc)))

  ft <- format_gp_fit_cli(objs$gp_fit)
  expect_true(any(grepl("4 units", ft)))
  expect_true(any(grepl("selected lambda = 0.25", ft)))
  expect_true(any(grepl("grades:", ft)))
})

test_that("the result headline precedes any pointer/footnote (result-first)", {
  objs <- ps_objects()
  for (type in names(objs)) {
    out <- ps_formatters[[type]](objs[[type]])
    note_idx <- grep("^i  ", out)
    first_note <- if (length(note_idx)) min(note_idx) else length(out) + 1L
    expect_gt(first_note, 1L)
    expect_true(any(nzchar(out[seq_len(first_note - 1L)])))
  }
})

# ---------------------------------------------------------------------------
# summary() -> typed gp_<type>_summary with a working shared printer
# ---------------------------------------------------------------------------

test_that("summary() returns a typed gp_<type>_summary with a working print", {
  objs <- ps_objects()
  for (type in names(objs)) {
    s <- summary(objs[[type]])
    expect_s3_class(s, paste0(type, "_summary"))
    expect_s3_class(s, "gp_summary")
    expect_true("provenance" %in% names(s))
    expect_output(print(s), regexp = ".")
    expect_invisible(print(s))
    out <- capture.output(print(s))
    expect_true(any(grepl("provenance", out)))
  }
})

test_that("provenance (backend/channel/rule) is surfaced ONLY in summary, never print", {
  objs <- ps_objects()

  s <- summary(objs$gp_grade_fit)
  expect_true("backend" %in% names(s$provenance))
  expect_true("channel" %in% names(s$provenance))
  expect_identical(s$provenance$backend, "highs")

  for (type in names(objs)) {
    printed <- capture.output(print(objs[[type]]))
    expect_false(any(grepl("provenance", printed, ignore.case = TRUE)), info = type)
  }
  expect_false(any(grepl("fixture", capture.output(print(objs$gp_grade_fit)))))

  expect_true("selection_rule" %in% names(summary(objs$gp_grade_path)$provenance))
  expect_true("selected_lambda_rule" %in% names(summary(objs$gp_frontier)$provenance))
  expect_true("selection_rule" %in% names(summary(objs$gp_fit)$provenance))
  for (type in c("gp_grade_path", "gp_frontier", "gp_fit")) {
    printed <- capture.output(print(objs[[type]]))
    expect_false(any(grepl("selection_rule|baseline_lambda", printed, ignore.case = TRUE)),
                 info = type)
  }
})

# ---------------------------------------------------------------------------
# GP-DEC-14-A (hard): NO claim/tier vocabulary in ANY rendered output.
#
# NOTE on the grep. The bare substring "tier" appears inside the legitimate word
# "fronTIER" (the object is gp_frontier and the design's own target box prints
# it). GP-DEC-14-A forbids the v1 claim/tier VOCABULARY, not the substring inside
# "frontier". So we forbid the claim/tier tokens with WORD BOUNDARIES (\\btier\\b
# matches "claim tier" but NOT "frontier") plus the exact v1 disclaimer phrases as
# fixed strings -- a strong grep that catches real violations with zero false
# positives. We additionally run the verbatim substring regex after stripping the
# single allowed word "frontier", proving no OTHER token smuggles in.
# ---------------------------------------------------------------------------

test_that("GP-DEC-14-A: no claim/tier vocabulary in any print/summary/format output", {
  objs <- ps_objects()
  word_re <- "\\bclaim\\b|claim-bearing|claim tier|claim_bearing|research_only|backend_claim|\\bbeats\\b|\\bwins\\b|\\btier\\b"
  v1_phrases <- c("claim_bearing", "claim-bearing", "claim tier",
                  "research_only", "backend_claim")

  for (type in names(objs)) {
    x <- objs[[type]]
    rendered <- c(
      ps_formatters[[type]](x),
      capture.output(print(x)),
      capture.output(print(summary(x)))
    )
    expect_false(any(grepl(word_re, rendered, ignore.case = TRUE)), info = type)
    for (phrase in v1_phrases) {
      expect_false(any(grepl(phrase, rendered, fixed = TRUE)), info = paste(type, phrase))
    }
    stripped <- gsub("frontier", "FRONT", rendered, ignore.case = TRUE)
    expect_false(
      any(grepl("claim|tier|research_only|beats|wins|backend_claim",
                stripped, ignore.case = TRUE)),
      info = type
    )
  }
})

# ---------------------------------------------------------------------------
# INVARIANT 2 (frozen): the firms<->names welfare-flip note in the report card.
# ---------------------------------------------------------------------------

test_that("INVARIANT 2: report card carries the firms<->names welfare-flip note", {
  card <- ps_report_card()
  out <- format_gp_report_card_cli(card)

  expect_true(any(grepl("grade 1", out)))
  expect_true(any(grepl("most discriminatory", out)))   # firms application
  expect_true(any(grepl("best-treated", out)))          # names application (it FLIPS)
  expect_true(any(grepl("Appendix A", out)))
  expect_true(any(grepl("\\?gp_report_card", out)))
  expect_false(any(grepl("\\bbeats\\b|\\bwins\\b", out, ignore.case = TRUE)))

  printed <- capture.output(print(card))
  note <- printed[grepl("grade 1", printed)]
  expect_length(note, 1L)
  expect_true(grepl("most discriminatory", note) && grepl("best-treated", note))
  expect_true(grepl("Appendix A", note))
})

test_that("INVARIANT 2: the welfare label is NOT hard-coded to one application", {
  card <- ps_report_card()
  note <- grep("grade 1", format_gp_report_card_cli(card), value = TRUE)
  expect_length(note, 1L)
  expect_true(grepl("firms", note))
  expect_true(grepl("names", note))
})

# ---------------------------------------------------------------------------
# Report-card table head/overflow (head + ... + final row) + NA CI.
# ---------------------------------------------------------------------------

test_that("report-card table shows a HEAD with overflow AND the final row, NA-safe", {
  ids <- paste0("u", sprintf("%02d", 1:12))
  grades <- as.integer(c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3))
  pm <- seq(0.5, 0.05, length.out = 12)
  tab <- data.frame(
    id = ids, label = c(NA, paste("Firm", 2:12)), grade = grades, sort_rank = 1:12,
    selected_lambda = 0.25, posterior_mean = pm, lower = pm - 0.05, upper = pm + 0.05,
    estimate = pm + 0.001, se = rep(0.03, 12), stringsAsFactors = FALSE
  )
  tab$lower[3] <- NA_real_; tab$upper[3] <- NA_real_
  card <- new_gp_report_card(ids = ids, table = tab, selected_lambda = 0.25,
                             grades = grades, control = ps_ctl())

  out <- format_gp_report_card_cli(card)                 # default head = 8 rows
  expect_true(any(grepl("\\[NA, NA\\]", out)))           # NA CI rendered cleanly
  expect_false(any(grepl("NaN|Inf", out)))               # no non-finite leakage
  expect_true(any(grepl("\\.\\.\\. 4 more rows \\.\\.\\.", out)))  # 12 - 8 = 4

  # The final (least-extreme) row appears below the ellipsis; the
  # hidden middle rows do not.
  expect_true(any(grepl("\\bu12\\b", out)))              # final row shown
  expect_false(any(grepl("\\bu09\\b", out)))             # a hidden middle row absent

  out3 <- format_gp_report_card_cli(card, max_rows = 3L)
  expect_true(any(grepl("\\.\\.\\. 9 more rows \\.\\.\\.", out3)))
  expect_true(any(grepl("\\bu12\\b", out3)))             # final row always kept
})

test_that("report-card truncates long labels to a prefix + '...' (no NA-corruption)", {
  # Regression guard: a length-1 `keep` indexed by the logical `too_long` used to
  # inject NA into every over-long label past the first, printing firm names as
  # "NA...". Two labels > 16 chars must BOTH truncate to a real prefix.
  ids <- paste0("f", 1:3)
  longlabs <- c("International Business Machines Corp",
                "Berkshire Hathaway Incorporated", "Short")
  pm <- c(0.3, 0.2, 0.1)
  tab <- data.frame(
    id = ids, label = longlabs, grade = c(1L, 2L, 3L), sort_rank = 1:3,
    selected_lambda = 0.25, posterior_mean = pm, lower = pm - 0.05, upper = pm + 0.05,
    estimate = pm, se = rep(0.03, 3), stringsAsFactors = FALSE
  )
  card <- new_gp_report_card(ids = ids, table = tab, selected_lambda = 0.25,
                             grades = c(1L, 2L, 3L), control = ps_ctl())
  out <- format_gp_report_card_cli(card)
  expect_false(any(grepl("NA\\.\\.\\.", out)))                 # the bug symptom
  expect_true(any(grepl("International", out)))                # 1st long label prefix
  expect_true(any(grepl("Berkshire Hat", out, fixed = TRUE)))  # 2nd long label (the victim)
})

test_that("as.data.frame.gp_report_card returns the 10-column payload (the print pointer works)", {
  card <- ps_report_card()
  df <- as.data.frame(card)
  expect_s3_class(df, "data.frame")
  expect_identical(nrow(df), nrow(card$table))
  expect_true(all(c("id", "label", "grade", "sort_rank", "posterior_mean",
                    "lower", "upper", "estimate", "se") %in% names(df)))
})

# ---------------------------------------------------------------------------
# No-recompute / no-solver sanity on minimal fixtures
# ---------------------------------------------------------------------------

test_that("print/summary do not error and do not call the solver on minimal fixtures", {
  objs <- ps_objects()
  for (type in names(objs)) {
    expect_error(format(ps_formatters[[type]](objs[[type]])), NA)
    expect_error(summary(objs[[type]]), NA)
    expect_error(capture.output(print(objs[[type]])), NA)
  }

  # gp_fit whose grade_path summary is NOT enriched: the metrics are not
  # materialized, so the view must show NA (a pure read), never recompute/solve.
  ft <- ps_fit()
  out <- format_gp_fit_cli(ft)
  expect_true(any(grepl("\\(1 - DR\\) = NA", out)))
  s <- summary(ft)
  expect_true(is.na(s$reliability))
  expect_identical(s$provenance$frontier_source, "not materialized")
})

test_that("gp_fit metrics ARE shown when the path summary is enriched", {
  pw <- ps_pairwise()
  enriched <- ps_grade_path(pw, enriched = TRUE)
  ft <- structure(
    list(ids = pw$ids, selected_grade = enriched$fits[[1L]],
         grade_path = enriched, control = pw$control),
    class = c("gp_fit", "list")
  )
  out <- format_gp_fit_cli(ft)
  expect_false(any(grepl("= NA", out)))
  expect_true(any(grepl("\\(1 - DR\\) = 1.00", out)))
  expect_identical(summary(ft)$provenance$frontier_source, "grade_path$summary (enriched)")
})

# ---------------------------------------------------------------------------
# ASCII-only + width plumbing
# ---------------------------------------------------------------------------

test_that("all rendered output is ASCII-only (no ANSI, no Unicode box glyphs)", {
  objs <- ps_objects()
  for (type in names(objs)) {
    out <- c(ps_formatters[[type]](objs[[type]]), capture.output(print(summary(objs[[type]]))))
    expect_false(any(grepl("[^\x01-\x7F]", out)), info = type)
    expect_false(any(grepl("\033\\[", out)), info = type)
  }
})

test_that("box formatters honour the width argument", {
  objs <- ps_objects()
  for (type in c("gp_pairwise", "gp_grade_fit", "gp_grade_path", "gp_frontier", "gp_fit")) {
    out <- ps_formatters[[type]](objs[[type]], width = 48L)
    box <- out[grepl("^\\+", out)]
    expect_true(all(nchar(box) <= 48L), info = type)
  }
})
