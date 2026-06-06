# shared Interpretation: / Glossary: scaffold on gp_summary.
# Extras are attribute-stored (gp_interpretation / gp_glossary); they render AFTER
# the scalars and BEFORE provenance; a no-extras summary is byte-identical to
# before; everything is ASCII-only and GP-DEC-14-A clean.

ss_bare <- function() {
  .gp_new_summary(
    "gp_fit",
    units = 97L, grade_count = 3L, composition = "3 (2/81/14)",
    selected_lambda = 0.25, reliability = 0.961, tau_bar = 0.207,
    provenance = list(backend = "gurobi")
  )
}
ss_extras <- function() {
  .gp_new_summary(
    "gp_fit",
    units = 97L, grade_count = 3L, composition = "3 (2/81/14)",
    selected_lambda = 0.25, reliability = 0.961, tau_bar = 0.207,
    interpretation = c(
      "97 units sorted into 3 grades (2/81/14) at lambda = 0.25.",
      "Grades agree with the pairwise ranking on 96.1% of comparable pairs.",
      "Within a grade, units are on average 0.207-reliable; higher is finer."
    ),
    glossary = c(
      reliability = "1 - discordance rate: share of comparable pairs ranked the same way.",
      tau_bar     = "average within-grade reliability across grades; 0 = coarse, 1 = fine."
    ),
    provenance = list(backend = "gurobi")
  )
}

test_that("a no-extras summary prints with no Interpretation/Glossary block (byte-exact)", {
  s <- ss_bare()
  out <- capture.output(print(s))
  expect_false(any(grepl("Interpretation:", out, fixed = TRUE)))
  expect_false(any(grepl("Glossary:", out, fixed = TRUE)))
  expect_identical(out, c(
    "<summary: gp_fit>",
    "  units           : 97",
    "  grade_count     : 3",
    "  composition     : 3 (2/81/14)",
    "  selected_lambda : 0.250",
    "  reliability     : 0.961",
    "  tau_bar         : 0.207",
    "  provenance:",
    "    backend : gurobi"
  ))
  expect_null(attr(s, "gp_interpretation", exact = TRUE))
  expect_null(attr(s, "gp_glossary", exact = TRUE))
})

test_that("an extras-bearing summary renders Interpretation and Glossary", {
  out <- capture.output(print(ss_extras()))
  expect_true(any(grepl("^  Interpretation:$", out)))
  expect_true(any(grepl("^  Glossary:$", out)))
  expect_true(any(grepl("97 units sorted into 3 grades", out)))
  expect_true(any(grepl("Grades agree with the pairwise ranking", out)))
  expect_true(any(grepl("Within a grade, units are on average", out)))
  expect_true(any(grepl("^    reliability : ", out)))
  expect_true(any(grepl("^    tau_bar     : ", out)))
})

test_that("Interpretation/Glossary sit after scalars and before provenance", {
  out <- capture.output(print(ss_extras()))
  i_scalar <- grep("^  selected_lambda : ", out)
  i_interp <- grep("^  Interpretation:$", out)
  i_gloss  <- grep("^  Glossary:$", out)
  i_prov   <- grep("^  provenance:$", out)
  expect_length(i_interp, 1L); expect_length(i_gloss, 1L); expect_length(i_prov, 1L)
  expect_true(i_scalar < i_interp)
  expect_true(i_interp < i_gloss)
  expect_true(i_gloss  < i_prov)
  expect_identical(max(c(i_interp, i_gloss, i_prov)), i_prov)
})

test_that("extras render is ASCII-only", {
  out <- capture.output(print(ss_extras()))
  expect_false(any(grepl("[^\x01-\x7F]", out)))
  expect_false(any(grepl("\033\\[", out)))
})

test_that("extras render carries no claim/tier vocabulary", {
  out <- capture.output(print(ss_extras()))
  word_re <- paste0("\\bclaim\\b|claim-bearing|claim tier|claim_bearing|",
                    "research_only|backend_claim|\\bbeats\\b|\\bwins\\b|\\btier\\b")
  expect_false(any(grepl(word_re, out, ignore.case = TRUE)))
  stripped <- gsub("frontier", "FRONT", out, ignore.case = TRUE)
  expect_false(any(grepl("claim|tier|research_only|beats|wins|backend_claim",
                         stripped, ignore.case = TRUE)))
})

test_that("extras honour width (wrap interpretation, truncate glossary)", {
  s <- .gp_new_summary(
    "gp_fit",
    units = 1L,
    interpretation = strrep("word ", 20L),
    glossary = c(metric = strrep("definition ", 20L)),
    provenance = list(backend = "highs")
  )
  old <- options(width = 50L); on.exit(options(old))
  out <- capture.output(print(s))
  expect_true(all(nchar(out) <= 50L))
  expect_length(grep("^  Interpretation:$", out), 1L)
  gloss_line <- grep("^    metric : ", out, value = TRUE)
  expect_length(gloss_line, 1L)
  expect_match(gloss_line, "\\.\\.\\.$")
})

test_that("empty/NULL/NA extras render no block; unnamed glossary entry dropped", {
  expect_identical(
    capture.output(print(.gp_new_summary(
      "gp_fit", units = 97L, grade_count = 3L, composition = "3 (2/81/14)",
      selected_lambda = 0.25, reliability = 0.961, tau_bar = 0.207,
      interpretation = NULL, glossary = NULL,
      provenance = list(backend = "gurobi")
    ))),
    capture.output(print(ss_bare()))
  )
  s_empty <- .gp_new_summary(
    "gp_fit", units = 1L,
    interpretation = character(0),
    glossary = stats::setNames(character(0), character(0)),
    provenance = list(backend = "highs")
  )
  expect_false(any(grepl("Interpretation:|Glossary:", capture.output(print(s_empty)))))

  s_blank <- .gp_new_summary(
    "gp_fit", units = 1L,
    interpretation = c(NA_character_, ""),
    glossary = c(NA_character_, ""),
    provenance = list(backend = "highs")
  )
  expect_false(any(grepl("Interpretation:|Glossary:", capture.output(print(s_blank)))))

  s_mix <- .gp_new_summary(
    "gp_fit", units = 1L,
    glossary = c(reliability = "1 - discordance rate.", "orphan with no term"),
    provenance = list(backend = "highs")
  )
  out_mix <- capture.output(print(s_mix))
  expect_true(any(grepl("^    reliability : 1 - discordance rate\\.$", out_mix)))
  expect_false(any(grepl("^    +: ", out_mix)))
  expect_false(any(grepl("orphan with no term", out_mix)))
})

test_that("a glossary term with no matching scalar still renders", {
  s <- .gp_new_summary(
    "gp_fit", units = 97L, reliability = 0.961,
    glossary = c(discordance = "share of comparable pairs the grades rank oppositely."),
    provenance = list(backend = "gurobi")
  )
  out <- capture.output(print(s))
  expect_false("discordance" %in% names(s))
  expect_true(any(grepl("^    discordance : ", out)))
})
