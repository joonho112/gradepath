# =============================================================================
# test-plots-report-card.R -- report-card + discordance figures
# -----------------------------------------------------------------------------
# Structural tests always run (they only need ggplot2). The strongest cheap
# guarantee is ggplot2::ggplot_build(), which forces every layer/scale to
# compute. vdiffr snapshot tests are gated and run only when vdiffr is present.
# Fixture builders are defined locally so the file is self-contained (no solver).
#
# ggplot2 4.x notes: the caterpillar uses geom_pointrange(orientation = "y")
# with point size controlled by `size` (not deprecated `fatten`). coord_equal()
# yields a CoordCartesian with ratio 1.
# =============================================================================

rc_ids <- c("u1", "u2", "u3", "u4")
rc_control <- function() gp_control(backend = "highs")

# ---- report-card fixture builders (mirror test-report-card.R) --------------

rc_grade_fit <- function(lambda, grades, ids = rc_ids, control = rc_control()) {
  new_gp_grade_fit(
    ids = ids,
    lambda = lambda,
    assignment = data.frame(id = ids, grade = as.integer(grades), stringsAsFactors = FALSE),
    summary = list(grade_count = as.integer(length(unique(grades))), status = "optimal"),
    objective = list(value = 0, raw = 0, canonical = 0),
    backend = list(name = control$backend, path = "fixture"),
    control = control
  )
}

rc_grade_path <- function(ids = rc_ids, control = rc_control(),
                          selected = c(2L, 1L, 2L, 3L),
                          endpoint = c(2L, 1L, 1L, 2L)) {
  fits <- list(rc_grade_fit(0.25, selected, ids, control), rc_grade_fit(1, endpoint, ids, control))
  new_gp_grade_path(
    ids = ids,
    lambda_grid = c(0.25, 1),
    fits = fits,
    summary = data.frame(
      lambda = c(0.25, 1),
      grade_count = vapply(fits, function(f) f$summary$grade_count, integer(1))
    ),
    backend = list(name = control$backend),
    selection = list(selected_lambda = 0.25, selection_rule = "baseline_lambda_0.25", endpoint_lambda = 1),
    control = control
  )
}

rc_posterior <- function(ids = rc_ids, labels = paste("L", ids)) {
  pm <- c(0.4, -0.1, 0.2, 0.3) # u2 interval straddles 0
  new_gp_posterior(
    estimate = c(4, 1, 2, 3), se = rep(0.2, length(ids)), id = ids, label = labels,
    posterior_mean = pm, posterior_sd = rep(0.05, length(ids)),
    lower = pm - 0.2, upper = pm + 0.2, scale = "r"
  )
}

rc_posterior_positive <- function(ids = rc_ids, labels = paste("L", ids)) {
  pm <- c(0.9, 0.4, 0.6, 0.7)
  new_gp_posterior(
    estimate = c(4, 1, 2, 3), se = rep(0.2, length(ids)), id = ids, label = labels,
    posterior_mean = pm, posterior_sd = rep(0.05, length(ids)),
    lower = pm - 0.1, upper = pm + 0.1, scale = "r"
  )
}

make_card <- function(posterior = rc_posterior()) {
  gp_report_card(rc_grade_path(), posterior = posterior)
}

make_card_single_grade <- function() {
  path <- rc_grade_path(selected = c(1L, 1L, 1L, 1L), endpoint = c(1L, 1L, 1L, 1L))
  gp_report_card(path, posterior = rc_posterior())
}

# ---- frontier / dr_matrix fixture builders (mirror test-frontier.R) --------

strict_matrix <- function(n = 4, high = 0.9) {
  out <- matrix(0.5, n, n)
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      out[i, j] <- high
      out[j, i] <- 1 - high
    }
  }
  out
}

fr_pairwise <- function(m, ids, control = rc_control()) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  dimnames(m) <- list(ids, ids)
  new_gp_pairwise(
    ids = ids, matrix = m, power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = .gp_pairwise_diagonal, zero_floor = .gp_pairwise_zero_floor),
    source = list(stage = "posterior", rule = "outer_product", assumption = "one_level_independence"),
    control = control
  )
}

fr_fit <- function(ids, lambda, grades, control = rc_control()) {
  new_gp_grade_fit(
    ids = ids, lambda = lambda,
    assignment = data.frame(id = ids, grade = as.integer(grades), stringsAsFactors = FALSE),
    summary = list(grade_count = as.integer(length(unique(grades))), status = "optimal", n_units = length(ids)),
    objective = list(value = 0, raw = 0, canonical = 0),
    backend = list(name = control$backend, path = "fixture", encoding = "fixture", status = "optimal",
      objbound = NA_real_, mipgap = NA_real_, runtime = NA_real_, warmstart = "none",
      warm_start_from_lambda = NA_real_, warm_start_used = FALSE, problem_hash = "fixture"),
    control = control
  )
}

# A 4-grade SELECTED slice (lambda 0.25) -> fr$dr_matrix is 4x4 (6 lower cells).
make_frontier <- function(control = rc_control()) {
  ids <- paste0("unit_", 1:4)
  pw <- fr_pairwise(strict_matrix(4), ids, control)
  fits <- list(fr_fit(ids, 0.25, c(1L, 2L, 3L, 4L), control), fr_fit(ids, 1, c(1L, 2L, 3L, 4L), control))
  path <- new_gp_grade_path(
    ids = ids, lambda_grid = c(0.25, 1), fits = fits, summary = frontier_table(pw, fits),
    backend = list(name = control$backend),
    selection = list(selected_lambda = 0.25, selection_rule = "x", endpoint_lambda = 1), control = control
  )
  gp_frontier(path, pairwise = pw)
}

# A single-grade SELECTED slice -> 1x1 dr_matrix (no grade pairs).
make_frontier_single_grade <- function(control = rc_control()) {
  ids <- paste0("unit_", 1:4)
  pw <- fr_pairwise(strict_matrix(4), ids, control)
  fits <- list(fr_fit(ids, 0.25, c(1L, 1L, 1L, 1L), control), fr_fit(ids, 1, c(1L, 1L, 1L, 1L), control))
  path <- new_gp_grade_path(
    ids = ids, lambda_grid = c(0.25, 1), fits = fits, summary = frontier_table(pw, fits),
    backend = list(name = control$backend),
    selection = list(selected_lambda = 0.25, selection_rule = "x", endpoint_lambda = 1), control = control
  )
  gp_frontier(path, pairwise = pw)
}

make_fit <- function() {
  path <- rc_grade_path()
  posterior <- rc_posterior()
  card <- gp_report_card(
    list(unit_id = rc_ids, label = paste("Label", rc_ids), theta_hat = c(10, 20, 30, 40), s = rep(0.4, 4)),
    posterior = posterior, selected_grade = path$fits[[1L]], grade_path = path
  )
  new_gp_fit(
    ids = rc_ids,
    estimates = ebrecipe::eb_input(theta_hat = c(10, 20, 30, 40), s = rep(0.4, 4), unit_id = rc_ids),
    prior = new_gp_prior(support = seq(-1, 1, length.out = 11), density = rep(1 / 11, 11), mean = 0, scale = "r"),
    posterior = posterior, precision_fit = NULL,
    pairwise = fr_pairwise(strict_matrix(4), rc_ids),
    grade_path = path, selected_grade = path$fits[[1L]], report_card = card, control = rc_control()
  )
}

skip_no_ggplot2 <- function() {
  skip_if(!requireNamespace("ggplot2", quietly = TRUE), "ggplot2 not installed")
}
geom_classes <- function(p) vapply(p$layers, function(l) class(l$geom)[[1L]], character(1))
layer_index <- function(p, geom_class) which(geom_classes(p) == geom_class)
quiet_build <- function(p) {
  expect_warning(out <- ggplot2::ggplot_build(p), NA)
  out
}

# =============================================================================
# gp_plot_report_card
# =============================================================================

test_that("gp_plot_report_card() returns a buildable horizontal interval plot", {
  skip_no_ggplot2()
  card <- make_card()
  p <- gp_plot_report_card(card)

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "Posterior mean")
  expect_identical(p$labels$title, "Report card")
  expect_identical(p$labels$colour, "Grade")
  expect_no_error(quiet_build(p))

  expect_true("GeomPointrange" %in% geom_classes(p))
  built <- quiet_build(p)
  idx <- layer_index(p, "GeomPointrange")
  expect_length(idx, 1L)
  expect_equal(nrow(built$data[[idx]]), nrow(card$table))
})

test_that("gp_plot_report_card() builds without ggplot2 warnings", {
  skip_no_ggplot2()
  card <- make_card()
  fit <- make_fit()

  for (p in list(
    gp_plot_report_card(card),
    gp_plot_report_card(card, order = "grade"),
    gp_plot_report_card(card, order = "mean"),
    gp_plot_report_card(card, order = "label"),
    gp_plot_report_card(card, max_rows = 2),
    gp_plot_report_card(make_card_single_grade()),
    gp_plot_report_card(fit),
    ggplot2::autoplot(card)
  )) {
    quiet_build(p)
  }
})

test_that("gp_plot_report_card() puts the most-extreme unit on top", {
  skip_no_ggplot2()
  card <- make_card()
  p <- gp_plot_report_card(card, order = "rank")
  # sort_rank == 1 is id "u2"; the y factor's LAST level is drawn at the top.
  expect_identical(utils::tail(levels(p$data$id), 1L), "u2")
  expect_identical(utils::head(levels(p$data$id), 1L), "u4")
})

test_that("gp_plot_report_card() order= changes the y ordering", {
  skip_no_ggplot2()
  card <- make_card()
  p_rank <- gp_plot_report_card(card, order = "rank")
  p_mean <- gp_plot_report_card(card, order = "mean")
  p_grade <- gp_plot_report_card(card, order = "grade")
  p_label <- gp_plot_report_card(card, order = "label")
  for (p in list(p_rank, p_mean, p_grade, p_label)) expect_no_error(ggplot2::ggplot_build(p))

  # mean order: largest posterior_mean (u1 = 0.4) on top.
  expect_identical(utils::tail(levels(p_mean$data$id), 1L), "u1")
  expect_false(identical(levels(p_rank$data$id), levels(p_mean$data$id)))
})

test_that("gp_plot_report_card() zero_line toggles a dashed vline", {
  skip_no_ggplot2()
  card <- make_card() # u2 straddles 0
  expect_true("GeomVline" %in% geom_classes(gp_plot_report_card(card, zero_line = "auto")))
  expect_true("GeomVline" %in% geom_classes(gp_plot_report_card(card, zero_line = "on")))
  expect_false("GeomVline" %in% geom_classes(gp_plot_report_card(card, zero_line = "off")))

  card_pos <- make_card(posterior = rc_posterior_positive())
  expect_false("GeomVline" %in% geom_classes(gp_plot_report_card(card_pos, zero_line = "auto")))
  expect_true("GeomVline" %in% geom_classes(gp_plot_report_card(card_pos, zero_line = "on")))
})

test_that("gp_plot_report_card() max_rows truncates and annotates, never silently", {
  skip_no_ggplot2()
  card <- make_card()
  p_all <- gp_plot_report_card(card)
  p_top2 <- gp_plot_report_card(card, max_rows = 2)

  expect_equal(nrow(ggplot2::ggplot_build(p_all)$data[[layer_index(p_all, "GeomPointrange")]]), 4L)
  expect_equal(nrow(ggplot2::ggplot_build(p_top2)$data[[layer_index(p_top2, "GeomPointrange")]]), 2L)
  expect_match(p_top2$labels$subtitle, "Showing 2 of 4")
  expect_null(p_all$labels$subtitle)

  p_big <- gp_plot_report_card(card, max_rows = 99)
  expect_null(p_big$labels$subtitle)
  expect_error(gp_plot_report_card(card, max_rows = 0), "positive integer")
})

test_that("gp_plot_report_card() captions the CI level without recomputing it", {
  skip_no_ggplot2()
  card <- make_card()
  p_default <- gp_plot_report_card(card)
  expect_match(p_default$labels$caption, "90")
  expect_match(p_default$labels$caption, "credible interval")

  p_95 <- gp_plot_report_card(card, ci_level = 0.95)
  expect_match(p_95$labels$caption, "95")
  # The interval geometry is identical (read straight off the card).
  i_d <- layer_index(p_default, "GeomPointrange")
  i_9 <- layer_index(p_95, "GeomPointrange")
  expect_equal(
    ggplot2::ggplot_build(p_default)$data[[i_d]]$xmin,
    ggplot2::ggplot_build(p_95)$data[[i_9]]$xmin
  )
  expect_error(gp_plot_report_card(card, ci_level = 1.5), "ci_level")
})

test_that("gp_plot_report_card() carries the shared grade colour scale", {
  skip_no_ggplot2()
  p <- gp_plot_report_card(make_card())
  expect_true(any(vapply(p$scales$scales, function(s) "colour" %in% s$aesthetics, logical(1))))
})

test_that("gp_plot_report_card() accepts a gp_fit and rejects junk", {
  skip_no_ggplot2()
  fit <- make_fit()
  p <- gp_plot_report_card(fit)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  idx <- layer_index(p, "GeomPointrange")
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[idx]]), nrow(fit$report_card$table))
  expect_error(gp_plot_report_card(list(a = 1)), "gp_report_card")
})

test_that("gp_plot_report_card() builds for a degenerate single-grade card", {
  skip_no_ggplot2()
  card <- make_card_single_grade()
  expect_length(unique(card$table$grade), 1L)
  p <- gp_plot_report_card(card)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

# =============================================================================
# gp_plot_discordance
# =============================================================================

test_that("gp_plot_discordance() returns a buildable lower-triangular heatmap", {
  skip_no_ggplot2()
  fr <- make_frontier()
  p <- gp_plot_discordance(fr)

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "Conditional discordance")
  expect_no_error(ggplot2::ggplot_build(p))

  geoms <- geom_classes(p)
  expect_true("GeomTile" %in% geoms)
  expect_true("GeomText" %in% geoms) # annotate = TRUE default
  expect_s3_class(p$coordinates, "CoordCartesian")
  expect_equal(p$coordinates$ratio, 1)
})

test_that("gp_plot_discordance() draws exactly k*(k-1)/2 cells", {
  skip_no_ggplot2()
  fr <- make_frontier()
  k <- nrow(fr$dr_matrix)
  expected <- k * (k - 1L) / 2L
  expect_equal(sum(!is.na(fr$dr_matrix)), expected)

  p <- gp_plot_discordance(fr)
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[layer_index(p, "GeomTile")[[1L]]]]), expected)
  expect_equal(nrow(built$data[[layer_index(p, "GeomText")[[1L]]]]), expected)
})

test_that("gp_plot_discordance() annotate toggles the cell text", {
  skip_no_ggplot2()
  fr <- make_frontier()
  p_yes <- gp_plot_discordance(fr, annotate = TRUE)
  p_no <- gp_plot_discordance(fr, annotate = FALSE)
  expect_true("GeomText" %in% geom_classes(p_yes))
  expect_false("GeomText" %in% geom_classes(p_no))
  expect_no_error(ggplot2::ggplot_build(p_no))

  labs <- ggplot2::ggplot_build(p_yes)$data[[layer_index(p_yes, "GeomText")[[1L]]]]$label
  expect_true(all(grepl("^[0-9]+\\.[0-9]{3}$", labs)))
})

test_that("gp_plot_discordance() accepts gp_frontier, gp_fit, and a raw matrix", {
  skip_no_ggplot2()
  fr <- make_frontier()
  p_fr <- gp_plot_discordance(fr)
  p_fit <- gp_plot_discordance(make_fit())
  p_mat <- gp_plot_discordance(fr$dr_matrix)
  for (p in list(p_fr, p_fit, p_mat)) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
    expect_true("GeomTile" %in% geom_classes(p))
  }
  expect_error(gp_plot_discordance(list(a = 1)), "gp_frontier")
  expect_error(gp_plot_discordance(matrix(0, 2, 3)), "square")
})

test_that("gp_plot_discordance() uses plain ASCII grade labels", {
  skip_no_ggplot2()
  p <- gp_plot_discordance(make_frontier())
  x_labels <- ggplot2::ggplot_build(p)$layout$panel_params[[1L]]$x$get_labels()
  x_labels <- x_labels[!is.na(x_labels)]
  expect_true(any(grepl("Grade", x_labels)))
  expect_false(any(grepl("grade_", x_labels)))
})

test_that("gp_plot_discordance() builds a placeholder for a single-grade frontier", {
  skip_no_ggplot2()
  fr <- make_frontier_single_grade()
  expect_equal(sum(!is.na(fr$dr_matrix)), 0L)
  p <- gp_plot_discordance(fr)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_false("GeomTile" %in% geom_classes(p))
})

test_that(".gp_dr_text_color switches white/dark by luminance", {
  cols <- .gp_dr_text_color(c(0.01, 0.05, 0.90, 0.99), vmax = 1)
  expect_length(cols, 4L)
  expect_identical(cols[[1L]], "grey10") # light cell -> dark text
  expect_identical(cols[[4L]], "white") # dark cell -> white text
  expect_length(.gp_dr_text_color(numeric(0), vmax = 1), 0L)
})

test_that("gp_plot_discordance() auto-contrasts a spread of cell fills", {
  skip_no_ggplot2()
  # A raw lower-triangular dr_matrix spanning light -> dark fills.
  m <- matrix(NA_real_, 3, 3, dimnames = list(paste0("grade_", 1:3), paste0("grade_", 1:3)))
  m[2, 1] <- 0.02
  m[3, 1] <- 0.05
  m[3, 2] <- 0.45
  p <- gp_plot_discordance(m)
  built <- ggplot2::ggplot_build(p)
  txt <- built$data[[layer_index(p, "GeomText")[[1L]]]]
  # A varied fill range must produce both light- and dark-text decisions.
  expect_gt(length(unique(txt$colour)), 1L)
})

# =============================================================================
# autoplot
# =============================================================================

test_that("autoplot.gp_report_card dispatches to the report-card plot", {
  skip_no_ggplot2()
  card <- make_card()
  p <- ggplot2::autoplot(card)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Report card")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_match(ggplot2::autoplot(card, max_rows = 2)$labels$subtitle, "Showing 2 of 4")
})

# =============================================================================
# vdiffr snapshot tests (gated)
# =============================================================================

test_that("report-card figure matches its snapshot", {
  skip_if_not_installed("vdiffr")
  skip_no_ggplot2()
  vdiffr::expect_doppelganger("report-card-rank", gp_plot_report_card(make_card()))
  vdiffr::expect_doppelganger("report-card-grade-order", gp_plot_report_card(make_card(), order = "grade"))
})

test_that("discordance figure matches its snapshot", {
  skip_if_not_installed("vdiffr")
  skip_no_ggplot2()
  vdiffr::expect_doppelganger("discordance-4grade", gp_plot_discordance(make_frontier()))
  vdiffr::expect_doppelganger("discordance-no-annotate", gp_plot_discordance(make_frontier(), annotate = FALSE))
})
