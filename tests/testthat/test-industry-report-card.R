# =============================================================================
# test-industry-report-card.R -- two-level INDUSTRY report cards
# -----------------------------------------------------------------------------
# The firm baseline, the INDUSTRY flavor (industry random effects -> a FIRM
# report card with sharper posteriors, ~4 grades), and the BETWEEN-industry
# flavor (`btwn` -> an INDUSTRY report card of SIC-coded units) are the SAME
# gp_report_card object on different inputs (per the design). gp_report_card
# (R/report-card.R) is generic over units; gp_twolevel_report_card
# (R/industry-grading.R) already returns a valid gp_report_card for both models.
# This test adds NO new R code -- it proves that contract AND that the Phase-10
# M3 surfaces built in Steps 10.1/10.2 work on INDUSTRY-flavored cards:
#   * gp_plot_report_card()              (the KRW Fig 9 / 16 figure)
#   * format_gp_report_card_cli() / print.gp_report_card  (the console box)
#   * as.data.frame.gp_report_card       (the 10-column payload)
#   * summary.gp_report_card             (the typed summary)
#
# Two complementary fixture strategies (BOTH used):
#   A. A fast SYNTHETIC industry-flavored gp_report_card built directly via the
#      stage-wise constructor (no solver) -- drives the M3-surface tests, which
#      therefore ALWAYS execute (no highs gate).
#   B. The REAL two-level path via the cheap deterministic .ircg_* helpers
#      (copied from test-industry-grading.R) -- gated on the highs backend.
#
# testthat edition 3, self-contained, ASCII-only.
# =============================================================================

# ---------------------------------------------------------------------------
# Fixture A -- synthetic industry-flavored gp_report_card (FAST, no solver)
# ---------------------------------------------------------------------------

.ircg_control <- function() {
  gp_control(backend = "highs")
}

# Build a gp_report_card directly over arbitrary INDUSTRY-flavored units/labels
# (same stage-wise pattern as tests/testthat/test-report-card.R). `grades_sel`
# is the SELECTED-lambda assignment (drives the displayed grade column / grade
# count); `grades_end` is the endpoint (lambda = 1) assignment (drives the row
# sort order). Both must be contiguous integers from 1. No solve happens here.
.ircg_card <- function(ids, labels, grades_sel, grades_end, pm,
                       control = .ircg_control()) {
  n <- length(ids)
  fit <- function(lambda, grades) {
    new_gp_grade_fit(
      ids = ids,
      lambda = lambda,
      assignment = data.frame(
        id = ids, grade = as.integer(grades), stringsAsFactors = FALSE
      ),
      summary = list(
        grade_count = length(unique(grades)), status = "optimal"
      ),
      objective = list(value = 0, raw = 0, canonical = 0),
      backend = list(name = "highs", path = "fixture"),
      control = control
    )
  }
  fits <- list(fit(0.25, grades_sel), fit(1, grades_end))
  path <- new_gp_grade_path(
    ids = ids,
    lambda_grid = c(0.25, 1),
    fits = fits,
    summary = data.frame(
      lambda = c(0.25, 1),
      grade_count = c(length(unique(grades_sel)), length(unique(grades_end)))
    ),
    backend = list(name = "highs"),
    selection = list(
      selected_lambda = 0.25,
      selection_rule = "synthetic_fixture",
      endpoint_lambda = 1
    ),
    control = control
  )
  posterior <- new_gp_posterior(
    estimate = pm, se = rep(0.04, n), id = ids, label = labels,
    posterior_mean = pm, posterior_sd = rep(0.03, n),
    lower = pm - 0.06, upper = pm + 0.06, scale = "r"
  )
  gp_report_card(
    list(unit_id = ids, label = labels, theta_hat = pm, s = rep(0.04, n)),
    posterior = posterior,
    selected_grade = path$fits[[1L]],
    grade_path = path
  )
}

# btwn INDUSTRY card: 6 SIC-coded industry units, 4 grades at the selected
# lambda. (The endpoint fit splits all 6 apart, fixing the row order.)
.ircg_btwn_ids <- paste0("sic_", c(20, 28, 35, 36, 48, 50))
.ircg_btwn_labels <- paste0(
  c("Manufacturing", "Chemicals", "Machinery", "Electronics", "Telecom",
    "Wholesale"),
  " (", c(20, 28, 35, 36, 48, 50), ")"
)
.ircg_btwn_card <- function() {
  .ircg_card(
    ids = .ircg_btwn_ids,
    labels = .ircg_btwn_labels,
    grades_sel = c(1, 1, 2, 3, 3, 4),
    grades_end = c(1, 2, 3, 4, 5, 6),
    pm = c(0.30, 0.22, 0.12, 0.05, 0.02, -0.03)
  )
}

# industry_rfe FIRM card: 12 firm units, 4 grades at the selected lambda -- one
# more than the firm baseline's 3, reflecting the sharper two-level posterior
# (per the design: industry flavor -> 4 grades vs baseline 3).
.ircg_firm_ids <- paste0("firm_", sprintf("%02d", 1:12))
.ircg_firm_labels <- paste("Firm", 1:12)
.ircg_firm_pm <- seq(0.40, -0.15, length.out = 12)
.ircg_rfe_firm_card <- function() {
  .ircg_card(
    ids = .ircg_firm_ids,
    labels = .ircg_firm_labels,
    grades_sel = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4), # 4 grades
    grades_end = 1:12,
    pm = .ircg_firm_pm
  )
}

# firm BASELINE card: the SAME 12 firms graded more coarsely into 3 grades, the
# one-level baseline against which the industry flavor is compared.
.ircg_baseline_firm_card <- function() {
  .ircg_card(
    ids = .ircg_firm_ids,
    labels = .ircg_firm_labels,
    grades_sel = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3), # 3 grades
    grades_end = 1:12,
    pm = .ircg_firm_pm
  )
}

# ---------------------------------------------------------------------------
# Fixture B -- REAL two-level path helpers (copied from test-industry-grading.R)
# ---------------------------------------------------------------------------

.ircg_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  get(name, mode = "function")
}
.ircg_quad <- function(...) .ircg_get("gp_twolevel_quadrature")(...)
.ircg_grade <- function(...) .ircg_get("gp_twolevel_grade")(...)
.ircg_tlcard <- function(...) .ircg_get("gp_twolevel_report_card")(...)

# Cheap deterministic race-characteristic two-level case (3 firms, 2 industries).
.ircg_tiny_prior <- function() {
  supp_xi <- c(0.45, 1.25)
  g_xi <- c(0.45, 0.55)
  structure(
    list(
      support = supp_xi,
      density = g_xi,
      mean = sum(supp_xi * g_xi),
      scale = "r",
      diagnostics = list(group_fx = 1L, support_eta = c(0.70, 1.45),
                         g_eta = c(0.40, 0.60)),
      metadata = list(characteristic = "race", beta = 0.55, mu = 0)
    ),
    class = c("gp_prior", "list")
  )
}

.ircg_tiny_case <- function() {
  s <- c(0.80, 1.35, 1.05)
  beta <- 0.55
  mu <- 0
  v_hat <- c(0.72, 1.34, 1.08)
  theta_hat <- v_hat * s^beta
  prior <- .ircg_tiny_prior()
  prior$metadata$beta <- beta
  prior$metadata$mu <- mu
  list(
    input = list(
      theta_hat = theta_hat,
      s = s,
      industry = c(1, 1, 2),
      unit_id = paste0("race_", 1:3),
      label = paste("race", 1:3)
    ),
    prior = prior,
    fit = list(
      characteristic = "race",
      beta = beta,
      mu = mu,
      v_hat = v_hat,
      s_v = c(0.30, 0.34, 0.28),
      industry = c(1, 1, 2)
    )
  )
}

# Solve the real cheap two-level grade once (requires the highs backend). Returns
# the graded object plus the quadrature carrier (for the pairwise-id mirror).
.ircg_real_graded <- function() {
  case <- .ircg_tiny_case()
  q <- .ircg_quad(
    input = case$input,
    prior = case$prior,
    fit = case$fit,
    include_g_theta = FALSE,
    control = .ircg_control()
  )
  graded <- .ircg_grade(q, control = .ircg_control(), lambda_grid = c(0.25, 1))
  list(quad = q, graded = graded)
}

# ===========================================================================
# 1. The synthetic btwn industry card is a valid gp_report_card
# ===========================================================================

test_that("synthetic btwn industry card is a valid gp_report_card with SIC labels", {
  card <- .ircg_btwn_card()

  expect_s3_class(card, "gp_report_card")
  expect_silent(validate_gp_report_card(card))

  # 10-column payload with the canonical column set, one row per industry.
  expect_identical(names(card$table), .gp_report_card_table_columns)
  expect_identical(ncol(card$table), 10L)
  expect_identical(nrow(card$table), 6L)

  # ids are the SIC-coded industry units; slot and table agree in row order.
  expect_identical(card$ids, .ircg_btwn_ids)
  expect_identical(as.character(card$table$id), .ircg_btwn_ids)

  # Grades contiguous 1..k (k = 4) and the SIC labels are carried verbatim.
  grade_levels <- sort(unique(as.integer(card$table$grade)))
  expect_identical(grade_levels, 1:4)
  expect_identical(length(unique(card$grades)), 4L)
  expect_true(all(.ircg_btwn_labels %in% card$table$label))

  # The selected lambda is materialized consistently on the card and the table.
  expect_equal(card$selected_lambda, 0.25)
  expect_equal(card$table$selected_lambda, rep(0.25, 6L))
})

# ===========================================================================
# 2. The M3 FIGURE works on the industry card (KRW Fig 9 / 16)
# ===========================================================================

# Locate the built data for the caterpillar's geom_pointrange layer by geom
# class -- robust to the optional leading geom_vline (the zero reference line is
# drawn first when the intervals straddle zero, so it cannot be assumed to be
# layer 1 or 2).
.ircg_pointrange_data <- function(built) {
  idx <- which(vapply(
    built$plot$layers,
    function(layer) inherits(layer$geom, "GeomPointrange"),
    logical(1)
  ))
  testthat::expect_length(idx, 1L)
  built$data[[idx]]
}

.ircg_quiet_build <- function(p) {
  testthat::expect_warning(out <- ggplot2::ggplot_build(p), NA)
  out
}

test_that("gp_plot_report_card() builds on the industry card with grade colours", {
  card <- .ircg_btwn_card()

  p <- gp_plot_report_card(card)
  expect_s3_class(p, "ggplot")
  built <- expect_no_error(.ircg_quiet_build(p))

  # The industry units appear on the y axis (the figure shows full SIC labels;
  # only the console box truncates). The endpoint-sorted figure draws the
  # most-extreme unit on top, so the labels are present in reverse row order.
  y_labels <- built$layout$panel_params[[1L]]$y$get_labels()
  expect_true(all(.ircg_btwn_labels %in% y_labels))

  # Grade colours come from the shared ordinal palette: one colour per grade
  # present on the FULL card (4 grades -> 4 distinct colours). One point-range
  # row per industry unit.
  pr <- .ircg_pointrange_data(built)
  expect_identical(nrow(pr), 6L)
  expect_identical(length(unique(pr$colour)), 4L)

  # order= and max_rows= work on the industry card.
  expect_no_error(.ircg_quiet_build(gp_plot_report_card(card, order = "grade")))
  expect_no_error(.ircg_quiet_build(gp_plot_report_card(card, order = "mean")))
  expect_no_error(.ircg_quiet_build(gp_plot_report_card(card, order = "label")))

  p_top <- gp_plot_report_card(card, max_rows = 3)
  built_top <- expect_no_error(.ircg_quiet_build(p_top))
  expect_identical(nrow(.ircg_pointrange_data(built_top)), 3L) # only top 3 rows

  # autoplot() dispatches to the same figure.
  expect_s3_class(ggplot2::autoplot(card), "ggplot")
})

# ===========================================================================
# 3. The M3 CONSOLE works on the industry card (result-first box + payload)
# ===========================================================================

test_that("console surface leads with the result and carries SIC labels + welfare note", {
  card <- .ircg_btwn_card()

  lines <- format_gp_report_card_cli(card)
  expect_type(lines, "character")

  # Leads with the result: row count, grade composition, selected lambda.
  expect_match(lines[[1L]], "<gp_report_card>", fixed = TRUE)
  expect_match(lines[[1L]], "6 rows", fixed = TRUE)
  expect_match(lines[[1L]], "selected lambda: 0.25", fixed = TRUE)
  expect_match(lines[[1L]], "grades: 4", fixed = TRUE)

  # A SIC label that survives the 16-char body truncation appears in the table.
  expect_true(any(grepl("Chemicals", lines, fixed = TRUE)))

  # Frozen INVARIANT 2 welfare note (firms <-> names application flip), and NO
  # "beats"/"wins"/"claim"/"tier" verdict prose.
  expect_true(any(grepl("most-extreme theta", lines, fixed = TRUE)))
  joined <- paste(lines, collapse = "\n")
  expect_false(grepl("beats|wins|\\bclaim\\b|\\btier\\b", joined))

  # print() emits the same lines and returns the object invisibly.
  expect_output(ret <- withVisible(print(card)), "<gp_report_card>")
  expect_false(ret$visible)
  expect_identical(ret$value, card)

  # as.data.frame() returns the 10-column payload with all SIC labels intact.
  df <- as.data.frame(card)
  expect_s3_class(df, "data.frame")
  expect_identical(dim(df), c(6L, 10L))
  expect_identical(names(df), .gp_report_card_table_columns)
  expect_true(all(.ircg_btwn_labels %in% df$label))

  # summary() is a typed gp_report_card summary leading with the headline scalars.
  s <- summary(card)
  expect_s3_class(s, "gp_report_card_summary")
  expect_s3_class(s, "gp_summary")
  expect_identical(s$units, 6L)
  expect_identical(s$grade_count, 4L)
  expect_equal(s$selected_lambda, 0.25)
})

# ===========================================================================
# 4. Industry FLAVOR vs firm BASELINE: the industry flavor supports MORE grades
# ===========================================================================

test_that("industry-flavor firm card supports more grades than the firm baseline", {
  industry_flavor <- .ircg_rfe_firm_card()
  baseline <- .ircg_baseline_firm_card()

  expect_s3_class(industry_flavor, "gp_report_card")
  expect_s3_class(baseline, "gp_report_card")
  expect_silent(validate_gp_report_card(industry_flavor))
  expect_silent(validate_gp_report_card(baseline))

  flavor_grades <- length(unique(industry_flavor$grades))
  baseline_grades <- length(unique(baseline$grades))

  # Per the design: the sharper two-level posterior supports 4 grades
  # against the one-level baseline's 3. Both grade sets are contiguous from 1.
  expect_identical(flavor_grades, 4L)
  expect_identical(baseline_grades, 3L)
  expect_gt(flavor_grades, baseline_grades)
  expect_identical(sort(unique(industry_flavor$table$grade)), 1:4)
  expect_identical(sort(unique(baseline$table$grade)), 1:3)
})

# ===========================================================================
# 5. btwn industry units vs industry_rfe firm units (same object, diff inputs)
# ===========================================================================

test_that("btwn card has industry units and rfe card has firm units, differing", {
  btwn <- .ircg_btwn_card()
  rfe_firm <- .ircg_rfe_firm_card()

  # One row per industry on the btwn card; one row per firm on the rfe card.
  expect_identical(nrow(btwn$table), 6L)
  expect_identical(nrow(rfe_firm$table), 12L)
  expect_false(nrow(btwn$table) == nrow(rfe_firm$table))

  # The unit identifiers are different kinds of object: SIC industry codes vs
  # firm ids, with disjoint id sets.
  expect_true(all(grepl("^sic_", btwn$ids)))
  expect_true(all(grepl("^firm_", rfe_firm$ids)))
  expect_length(intersect(btwn$ids, rfe_firm$ids), 0L)

  # Both nonetheless render through the shared M3 figure (same object, diff input).
  expect_no_error(.ircg_quiet_build(gp_plot_report_card(btwn)))
  expect_no_error(.ircg_quiet_build(gp_plot_report_card(rfe_firm)))
})

# ===========================================================================
# 6. REAL two-level path (gated on highs): firm + between-industry cards + M3
# ===========================================================================

test_that("real two-level firm and between-industry cards are valid gp_report_cards", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  real <- .ircg_real_graded()
  q <- real$quad
  graded <- real$graded

  firm_card <- .ircg_tlcard(graded, model = "industry_rfe")
  industry_card <- .ircg_tlcard(graded, model = "btwn")

  expect_s3_class(firm_card, "gp_report_card")
  expect_s3_class(industry_card, "gp_report_card")
  expect_silent(validate_gp_report_card(firm_card))
  expect_silent(validate_gp_report_card(industry_card))

  # The firm card's ids are the firm pairwise ids; the industry card's ids are
  # the industry pairwise ids (mirrors test-industry-grading.R:159-189). The two
  # id sets differ -- firm units vs industry units.
  expect_identical(sort(firm_card$ids), sort(q$pairwise_theta$ids))
  expect_identical(sort(industry_card$ids), sort(q$pairwise_bar$ids))
  expect_identical(nrow(firm_card$table), length(q$pairwise_theta$ids))
  expect_identical(nrow(industry_card$table), length(q$pairwise_bar$ids))
  expect_false(identical(sort(firm_card$ids), sort(industry_card$ids)))

  # The public wrapper forwards to gp_twolevel_grade() when handed the ungraded
  # quadrature carrier directly, returning the same firm card.
  wrapped <- gp_twolevel_report_card(
    q, model = "industry_rfe",
    control = .ircg_control(), lambda_grid = c(0.25, 1)
  )
  expect_s3_class(wrapped, "gp_report_card")
  expect_identical(sort(wrapped$ids), sort(firm_card$ids))
})

test_that("M3 surfaces work on the real two-level industry cards", {
  skip_if(!("highs" %in% .gp_available_open_backends()), "highs backend unavailable")
  real <- .ircg_real_graded()
  graded <- real$graded

  firm_card <- .ircg_tlcard(graded, model = "industry_rfe")
  industry_card <- .ircg_tlcard(graded, model = "btwn")

  for (card in list(firm_card, industry_card)) {
    # Figure builds.
    expect_s3_class(gp_plot_report_card(card), "ggplot")
    expect_no_error(.ircg_quiet_build(gp_plot_report_card(card)))

    # Console leads with the result and carries the frozen welfare note;
    # print() returns invisibly.
    lines <- format_gp_report_card_cli(card)
    expect_match(lines[[1L]], "<gp_report_card>", fixed = TRUE)
    expect_true(any(grepl("most-extreme theta", lines, fixed = TRUE)))
    expect_output(ret <- withVisible(print(card)), "<gp_report_card>")
    expect_false(ret$visible)

    # as.data.frame() is the 10-column payload; summary() is the typed summary.
    df <- as.data.frame(card)
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), .gp_report_card_table_columns)
    expect_identical(nrow(df), nrow(card$table))
    expect_s3_class(summary(card), "gp_report_card_summary")
  }
})
