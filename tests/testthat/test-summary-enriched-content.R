# The enriched summary.* methods carry a genuinely useful
# Interpretation + Glossary block on REAL objects. Cache-gated: skips cleanly
# (e.g. under R CMD check) when the precomputed workflow-test fit is unavailable.

.gp_find_cache <- function(name) {
  cands <- c(
    file.path("log/054_workflow-test-v1/outputs/cache", name),
    file.path("../../log/054_workflow-test-v1/outputs/cache", name),
    file.path("../../../log/054_workflow-test-v1/outputs/cache", name)
  )
  hit <- cands[file.exists(cands)]
  if (length(hit)) normalizePath(hit[[1]]) else NA_character_
}

test_that("summary.gp_fit carries a correct Interpretation + Glossary", {
  p <- .gp_find_cache("fit_race_parity.rds")
  skip_if(is.na(p), "cached real race fit unavailable")
  fit <- readRDS(p)
  out <- capture.output(print(summary(fit)))
  joined <- paste(out, collapse = " ")

  expect_true(any(grepl("^  Interpretation:$", out)))
  expect_true(any(grepl("^  Glossary:$", out)))
  # the grade distribution is NOT doubled (bug guard): "3 grades (2/81/14)"
  expect_match(joined, "3 grades (2/81/14)", fixed = TRUE)
  expect_false(grepl("(3 (2/81/14))", joined, fixed = TRUE))
  # reliability is glossed and interpreted (1 - discordance rate)
  expect_match(joined, "1 - discordance rate")
  expect_match(joined, "reliability")
  expect_match(joined, "tau-bar")
  # provenance still present and summary-only; GP-DEC-14-A clean
  expect_true(any(grepl("^  provenance:$", out)))
  expect_false(grepl("\\btier\\b", joined))
  expect_false(grepl("claim", joined, ignore.case = TRUE))
})

test_that("summary.gp_report_card carries the frozen INVARIANT 2 welfare note", {
  p <- .gp_find_cache("fit_race_parity.rds")
  skip_if(is.na(p), "cached real race fit unavailable")
  fit <- readRDS(p)
  card <- get_report_card(fit)
  out <- capture.output(print(summary(card)))
  joined <- paste(out, collapse = " ")

  expect_true(any(grepl("^  Interpretation:$", out)))
  # INVARIANT 2 (firms<->names welfare flip), verbatim from the frozen note --
  # joined so the assertion survives line-wrapping of the long note.
  expect_match(joined, "grade 1 = most-extreme theta")
  expect_match(joined, "firms: *most discriminatory; *names: *best-treated")
  expect_match(joined, "Appendix A")
  expect_match(joined, "3 grades (2/81/14)", fixed = TRUE)
  expect_false(grepl("\\btier\\b", joined))
  expect_false(grepl("\\bbeats\\b|\\bwins\\b", joined))
})

test_that("summary.gp_frontier + summary.gp_grade_path carry interpretation + glossary", {
  p <- .gp_find_cache("fit_race_parity.rds")
  skip_if(is.na(p), "cached real race fit unavailable")
  fit <- readRDS(p)

  # frontier built solve-free from the fit's existing grade_path + pairwise
  fr <- gp_frontier(fit$grade_path, pairwise = get_pairwise(fit))
  of <- paste(capture.output(print(summary(fr))), collapse = " ")
  expect_match(of, "information-reliability frontier")
  expect_match(of, "1 - discordance rate")
  expect_match(of, "tau_bar")
  expect_false(grepl("\\btier\\b", of))

  og <- paste(capture.output(print(summary(fit$grade_path))), collapse = " ")
  expect_match(og, "path of optimal gradings")
  expect_match(og, "selected lambda = 0.25")
  expect_match(og, "number of grades ranges from")
  expect_false(grepl("\\btier\\b", og))
})

test_that("summary.gp_pairwise + summary.gp_grade_fit carry interpretation + glossary", {
  p <- .gp_find_cache("fit_race_parity.rds")
  skip_if(is.na(p), "cached real race fit unavailable")
  fit <- readRDS(p)

  # normalize whitespace so assertions survive the printer's line-wrap + indent
  norm <- function(x) gsub("[[:space:]]+", " ", paste(x, collapse = " "))
  op <- norm(capture.output(print(summary(get_pairwise(fit)))))
  expect_match(op, "pairwise outranking matrix Pi")
  expect_match(op, "more extreme than unit j")
  expect_match(op, "diagonal is fixed at 0.5")
  expect_false(grepl("\\btier\\b", op))

  gf <- gp_select_grade(fit$grade_path, lambda = 0.25)
  ogf <- norm(capture.output(print(summary(gf))))
  expect_match(ogf, "single solved grading slice")
  expect_match(ogf, "97 units graded into 3 grades at lambda = 0.25")
  expect_match(ogf, "Solver status")
  expect_false(grepl("\\btier\\b", ogf))
})

test_that("all six summary types render consistently and stay GP-DEC-14-A clean", {
  p <- .gp_find_cache("fit_race_parity.rds")
  skip_if(is.na(p), "cached real race fit unavailable")
  fit <- readRDS(p)
  pw <- get_pairwise(fit)
  objs <- list(
    gp_fit         = fit,
    gp_report_card = get_report_card(fit),
    gp_pairwise    = pw,
    gp_grade_fit   = gp_select_grade(fit$grade_path, lambda = 0.25),
    gp_grade_path  = fit$grade_path,
    gp_frontier    = gp_frontier(fit$grade_path, pairwise = pw)
  )
  word_re <- paste0("\\bclaim\\b|claim-bearing|claim tier|\\btier\\b|\\bbeats\\b|",
                    "\\bwins\\b|research_only|backend_claim")
  for (nm in names(objs)) {
    out <- capture.output(print(summary(objs[[nm]])))
    # every enriched summary renders both blocks
    expect_true(any(grepl("^  Interpretation:$", out)), info = nm)
    expect_true(any(grepl("^  Glossary:$", out)), info = nm)
    # provenance present, summary-only, and after the blocks
    i_prov <- grep("^  provenance:$", out)
    expect_length(i_prov, 1L)
    expect_gt(i_prov, max(grep("^  (Interpretation|Glossary):$", out)))
    # ASCII-only; GP-DEC-14-A clean (the word "frontier" is allowed)
    stripped <- gsub("frontier", "FRONT", out, ignore.case = TRUE)
    expect_false(any(grepl(word_re, stripped, ignore.case = TRUE)), info = nm)
    expect_false(any(grepl("[^\x01-\x7F]", out)), info = nm)
    # the matching format_*_cli still returns without error
    fmt <- get(paste0("format_", nm, "_cli"), mode = "function")
    expect_silent(invisible(fmt(objs[[nm]])))
  }
})
