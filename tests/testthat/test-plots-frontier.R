# =============================================================================
# test-plots-frontier.R -- frontier + posterior-contrast figures
# -----------------------------------------------------------------------------
# Structural tests always run (they only need ggplot2). The strongest cheap
# guarantee is ggplot2::ggplot_build(), which forces every layer/scale to
# compute. vdiffr snapshot tests are gated and run only when vdiffr is present.
# Fixture builders are defined locally so the file is self-contained (no solver).
# =============================================================================

# ---- frontier fixture builders (mirror test-frontier.R) --------------------

fr_pairwise <- function(m = NULL, ids = NULL, control = gp_control(backend = "highs")) {
  if (is.null(m)) m <- matrix(c(0.5, 0.9, 0.1, 0.5), 2, byrow = TRUE)
  if (is.null(ids)) ids <- paste0("unit_", seq_len(nrow(m)))
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  dimnames(m) <- list(ids, ids)
  new_gp_pairwise(
    ids = ids,
    matrix = m,
    power = 0L,
    cleanup = list(
      antisymmetry = TRUE,
      diagonal = .gp_pairwise_diagonal,
      zero_floor = .gp_pairwise_zero_floor
    ),
    source = list(
      stage = "posterior",
      rule = "outer_product",
      assumption = "one_level_independence"
    ),
    control = control
  )
}

fr_fit <- function(ids, lambda, grades, control) {
  new_gp_grade_fit(
    ids = ids,
    lambda = lambda,
    assignment = data.frame(id = ids, grade = as.integer(grades), stringsAsFactors = FALSE),
    summary = list(
      grade_count = as.integer(length(unique(grades))),
      status = "optimal",
      n_units = length(ids)
    ),
    objective = list(value = 0, raw = 0, canonical = 0),
    backend = list(
      name = control$backend, path = "fixture", encoding = "fixture",
      status = "optimal", objbound = NA_real_, mipgap = NA_real_,
      runtime = NA_real_, warmstart = "none", warm_start_from_lambda = NA_real_,
      warm_start_used = FALSE, problem_hash = "fixture"
    ),
    control = control
  )
}

fr_path <- function(pairwise, fit_specs, enriched = TRUE) {
  control <- pairwise$control
  fits <- lapply(fit_specs, function(s) fr_fit(pairwise$ids, s$lambda, s$grades, control))
  lambda_grid <- vapply(fit_specs, function(s) s$lambda, numeric(1))
  summary <- if (enriched) {
    frontier_table(pairwise, fits)
  } else {
    data.frame(
      lambda = lambda_grid,
      grade_count = vapply(fits, function(f) f$summary$grade_count, integer(1))
    )
  }
  new_gp_grade_path(
    ids = pairwise$ids,
    lambda_grid = lambda_grid,
    fits = fits,
    summary = summary,
    backend = list(name = control$backend),
    selection = list(
      selected_lambda = 0.25,
      selection_rule = "baseline_lambda_0.25",
      endpoint_lambda = 1
    ),
    control = control
  )
}

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

# A 4-unit frontier with a non-trivial Pareto shape AND benchmarks.
make_frontier <- function(benchmarks = TRUE) {
  pw <- fr_pairwise(strict_matrix(4), ids = paste0("unit_", 1:4))
  path <- fr_path(pw, list(
    list(lambda = 0.25, grades = c(1L, 1L, 2L, 2L)),
    list(lambda = 0.5, grades = c(1L, 2L, 3L, 3L)),
    list(lambda = 1, grades = c(1L, 2L, 3L, 4L))
  ))
  if (benchmarks) {
    gp_frontier(path, pairwise = pw, benchmark_scores = list(
      raw_estimate = stats::setNames(c(1, 2, 3, 4), rev(pw$ids)),
      posterior_mean = c(4, 3, 2, 1),
      linear_shrinkage = c(4, 3, 1, 2)
    ))
  } else {
    gp_frontier(path, pairwise = pw)
  }
}

# A degenerate frontier whose selected slice has a single grade -> point (1, 0).
make_frontier_degenerate <- function() {
  pw <- fr_pairwise(strict_matrix(4), ids = paste0("unit_", 1:4))
  path <- fr_path(pw, list(
    list(lambda = 0.25, grades = c(1L, 1L, 1L, 1L)),
    list(lambda = 1, grades = c(1L, 2L, 3L, 4L))
  ))
  gp_frontier(path, pairwise = pw)
}

# ---- gp_fit fixture (report_fit recipe, adapted) ---------------------------

rc_ids <- c("u1", "u2", "u3", "u4")
rc_control <- function() gp_control(backend = "highs")

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

rc_grade_path <- function(ids = rc_ids, control = rc_control()) {
  selected <- rc_grade_fit(0.25, c(2L, 1L, 2L, 3L), ids, control)
  endpoint <- rc_grade_fit(1, c(2L, 1L, 1L, 2L), ids, control)
  fits <- list(selected, endpoint)
  lambda_grid <- c(0.25, 1)
  new_gp_grade_path(
    ids = ids,
    lambda_grid = lambda_grid,
    fits = fits,
    summary = data.frame(
      lambda = lambda_grid,
      grade_count = vapply(fits, function(f) f$summary$grade_count, integer(1))
    ),
    backend = list(name = control$backend),
    selection = list(
      selected_lambda = 0.25,
      selection_rule = "baseline_lambda_0.25",
      endpoint_lambda = 1
    ),
    control = control
  )
}

rc_posterior <- function(ids = rc_ids, labels = paste("Posterior", ids)) {
  pm <- c(0.4, 0.1, 0.2, 0.3)
  new_gp_posterior(
    estimate = c(4, 1, 2, 3),
    se = rep(0.2, length(ids)),
    id = ids,
    label = labels,
    posterior_mean = pm,
    posterior_sd = rep(0.05, length(ids)),
    lower = pm - 0.01,
    upper = pm + 0.01,
    scale = "r"
  )
}

rc_estimates <- function(ids = rc_ids, labels = paste("Label", ids)) {
  list(unit_id = ids, label = labels, theta_hat = c(10, 20, 30, 40), s = rep(0.4, length(ids)))
}

rc_pairwise <- function(ids = rc_ids, control = rc_control()) {
  m <- matrix(0.5, length(ids), length(ids), dimnames = list(ids, ids))
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      if (i != j) m[i, j] <- if (i < j) 0.6 else 0.4
    }
  }
  new_gp_pairwise(
    ids = ids,
    matrix = m,
    power = 0L,
    cleanup = list(
      antisymmetry = TRUE,
      diagonal = .gp_pairwise_diagonal,
      zero_floor = .gp_pairwise_zero_floor
    ),
    source = list(
      stage = "posterior",
      rule = "outer_product",
      assumption = "one_level_independence"
    ),
    control = control
  )
}

rc_prior <- function() {
  support <- seq(-1, 1, length.out = 11)
  density <- rep(1 / length(support), length(support))
  new_gp_prior(support = support, density = density, mean = sum(support * density), scale = "r")
}

# A gp_fit whose grade fits carry status, so both the posterior-contrast plot
# AND the frontier-from-fit convenience path work.
make_fit <- function() {
  path <- rc_grade_path()
  posterior <- rc_posterior()
  card <- gp_report_card(
    rc_estimates(),
    posterior = posterior,
    selected_grade = path$fits[[1L]],
    grade_path = path
  )
  new_gp_fit(
    ids = rc_ids,
    estimates = ebrecipe::eb_input(
      theta_hat = c(10, 20, 30, 40),
      s = rep(0.4, length(rc_ids)),
      unit_id = rc_ids
    ),
    prior = rc_prior(),
    posterior = posterior,
    precision_fit = NULL,
    pairwise = rc_pairwise(),
    grade_path = path,
    selected_grade = path$fits[[1L]],
    report_card = card,
    control = rc_control()
  )
}

# A gp_pairwise whose strict structure puts cells above the 0.8 threshold.
make_strict_pairwise <- function(n = 4) {
  ids <- paste0("u", seq_len(n))
  fr_pairwise(strict_matrix(n), ids = ids)
}

skip_no_ggplot2 <- function() {
  skip_if(!requireNamespace("ggplot2", quietly = TRUE), "ggplot2 not installed")
}

geom_classes <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[[1L]], character(1))
}

# =============================================================================
# theme + palette
# =============================================================================

test_that("theme_gradepath() returns a ggplot2 theme with a bottom legend", {
  skip_no_ggplot2()
  th <- theme_gradepath()
  expect_s3_class(th, "theme")
  expect_identical(th$legend.position, "bottom")
})

test_that(".gp_grade_colors() returns a named, ordinal, colour-blind palette", {
  expect_silent(cols <- .gp_grade_colors(4))
  expect_length(cols, 4L)
  expect_identical(names(cols), as.character(1:4))
  expect_true(all(grepl("^#", cols)))
  expect_false(identical(cols[[1L]], cols[[4L]]))
  expect_length(.gp_grade_colors(1), 1L)
  expect_error(.gp_grade_colors(0))
})

test_that("scale_*_gradepath() build discrete ggplot2 scales", {
  skip_no_ggplot2()
  expect_s3_class(scale_color_gradepath(), "Scale")
  expect_s3_class(scale_colour_gradepath(), "Scale")
  expect_s3_class(scale_fill_gradepath(), "Scale")
})

# =============================================================================
# gp_plot_frontier
# =============================================================================

test_that("gp_plot_frontier() returns a buildable ggplot with the right axes", {
  skip_no_ggplot2()
  fr <- make_frontier()
  p <- gp_plot_frontier(fr)

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "1 .* DR")
  expect_true(is.expression(p$labels$y) || grepl("tau", deparse(p$labels$y)))

  expect_no_error(ggplot2::ggplot_build(p))

  geoms <- geom_classes(p)
  expect_true("GeomLine" %in% geoms)
  expect_true("GeomPoint" %in% geoms)

  # The frontier line layer has one row per frontier table row.
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[1L]]), nrow(fr$table))
})

test_that("gp_plot_frontier() shows benchmarks as colour+shape series only when present and requested", {
  skip_no_ggplot2()
  fr <- make_frontier(benchmarks = TRUE)
  fr_nb <- make_frontier(benchmarks = FALSE)

  p_on <- gp_plot_frontier(fr, benchmarks = TRUE)
  p_off <- gp_plot_frontier(fr, benchmarks = FALSE)
  p_none <- gp_plot_frontier(fr_nb, benchmarks = TRUE)

  expect_no_error(ggplot2::ggplot_build(p_on))
  expect_no_error(ggplot2::ggplot_build(p_off))
  expect_no_error(ggplot2::ggplot_build(p_none))

  # Benchmarks are NOT a separate layer in the redesign: the plot is always two
  # layers (the grey Pareto line + the colour/shape-mapped points). No text.
  expect_identical(length(p_on$layers), 2L)
  expect_identical(length(p_off$layers), 2L)
  expect_identical(length(p_none$layers), 2L)
  expect_false(any(grepl("Text", geom_classes(p_on))))

  # With benchmarks present + requested, the mapped point layer carries exactly
  # nrow(benchmarks) MORE rows (the benchmark series) than without, and more
  # distinct colours.
  rows_on  <- nrow(ggplot2::ggplot_build(p_on)$data[[2L]])
  rows_off <- nrow(ggplot2::ggplot_build(p_off)$data[[2L]])
  expect_equal(rows_on - rows_off, nrow(fr$benchmarks))
  n_col_on  <- length(unique(ggplot2::ggplot_build(p_on)$data[[2L]]$colour))
  n_col_off <- length(unique(ggplot2::ggplot_build(p_off)$data[[2L]]$colour))
  expect_gt(n_col_on, n_col_off)
  # When the frontier carries no benchmarks, requesting them adds nothing.
  expect_identical(
    nrow(ggplot2::ggplot_build(p_none)$data[[2L]]),
    nrow(ggplot2::ggplot_build(p_off)$data[[2L]])
  )
})

test_that("gp_plot_frontier() handles highlight variants", {
  skip_no_ggplot2()
  fr <- make_frontier()

  expect_no_error(ggplot2::ggplot_build(gp_plot_frontier(fr, highlight = "selected")))
  expect_no_error(ggplot2::ggplot_build(gp_plot_frontier(fr, highlight = "all")))
  expect_no_error(ggplot2::ggplot_build(gp_plot_frontier(fr, highlight = 1)))

  # highlight = "all" puts every path point in one (Selected slice) series -> the
  # point layer shows a single colour; "selected" marks one point -> two colours.
  p_all <- gp_plot_frontier(fr, highlight = "all", benchmarks = FALSE, annotate = FALSE)
  p_sel <- gp_plot_frontier(fr, highlight = "selected", benchmarks = FALSE, annotate = FALSE)
  c_all <- ggplot2::ggplot_build(p_all)$data[[2L]]$colour
  c_sel <- ggplot2::ggplot_build(p_sel)$data[[2L]]$colour
  expect_identical(length(unique(c_all)), 1L)
  expect_identical(length(unique(c_sel)), 2L)
})

test_that("gp_plot_frontier() encodes series by colour AND shape with a bottom legend (no text labels)", {
  skip_no_ggplot2()
  p <- gp_plot_frontier(make_frontier())

  # exactly two layers: grey Pareto line + the colour/shape-mapped points; no text
  expect_identical(length(p$layers), 2L)
  expect_identical(sort(unique(geom_classes(p))), c("GeomLine", "GeomPoint"))
  expect_false(any(grepl("Text", geom_classes(p))))

  # the points map BOTH colour and shape (so each rule reads even in grayscale)
  mapping <- p$layers[[2L]]$mapping
  expect_false(is.null(mapping$colour))
  expect_false(is.null(mapping$shape))

  # a manual colour scale AND a manual shape scale are present (they merge into
  # one legend because they share name/breaks/labels)
  scale_aes <- vapply(p$scales$scales,
                      function(s) paste(s$aesthetics, collapse = ","), character(1))
  expect_true(any(grepl("colour|color", scale_aes)))
  expect_true(any(grepl("shape", scale_aes)))

  # legend at the bottom (2 rows is set via guide_legend(nrow = 2); position here)
  expect_identical(p$theme$legend.position, "bottom")
})

test_that("gp_plot_frontier() annotation toggles the selected-slice subtitle", {
  skip_no_ggplot2()
  fr <- make_frontier(benchmarks = FALSE)
  p_yes <- gp_plot_frontier(fr, annotate = TRUE)
  p_no <- gp_plot_frontier(fr, annotate = FALSE)
  # Annotation lives in the subtitle, not a floating layer (avoids overlap).
  expect_false(is.null(p_yes$labels$subtitle))
  expect_match(p_yes$labels$subtitle, "1 - DR")
  expect_null(p_no$labels$subtitle)
  expect_no_error(ggplot2::ggplot_build(p_yes))
})

test_that("gp_plot_frontier() accepts a gp_fit and a gp_grade_path", {
  skip_no_ggplot2()
  fit <- make_fit()
  p_fit <- gp_plot_frontier(fit)
  expect_s3_class(p_fit, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p_fit))

  pw <- fr_pairwise(strict_matrix(4), ids = paste0("unit_", 1:4))
  path <- fr_path(pw, list(
    list(lambda = 0.25, grades = c(1L, 1L, 2L, 2L)),
    list(lambda = 1, grades = c(1L, 2L, 3L, 4L))
  ))
  p_path <- gp_plot_frontier(path, pairwise = pw)
  expect_s3_class(p_path, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p_path))
})

test_that("gp_plot_frontier() builds for a degenerate single-grade selected slice", {
  skip_no_ggplot2()
  fr <- make_frontier_degenerate()
  sel_row <- fr$table[fr$table$lambda == 0.25, ]
  expect_equal(sel_row$reliability, 1, tolerance = 1e-8)
  expect_equal(sel_row$tau_bar, 0, tolerance = 1e-8)

  p <- gp_plot_frontier(fr)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gp_plot_frontier() rejects unsupported inputs", {
  skip_no_ggplot2()
  expect_error(gp_plot_frontier(list(a = 1)), "gp_frontier")
})

# =============================================================================
# gp_plot_posterior_contrast
# =============================================================================

test_that("gp_plot_posterior_contrast() returns a buildable square heatmap", {
  skip_no_ggplot2()
  fit <- make_fit()
  p <- gp_plot_posterior_contrast(fit)

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "grade-sorted")
  expect_match(p$labels$y, "grade-sorted")

  expect_no_error(ggplot2::ggplot_build(p))

  geoms <- geom_classes(p)
  expect_true("GeomRaster" %in% geoms)
  expect_true(any(c("GeomVline", "GeomHline") %in% geoms))

  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[1L]]), length(fit$ids)^2)

  # coord_equal yields a CoordCartesian with a fixed ratio of 1 (ggplot2 >= 4).
  expect_s3_class(p$coordinates, "CoordCartesian")
  expect_equal(p$coordinates$ratio, 1)
})

test_that("gp_plot_posterior_contrast() reorders by (selected, endpoint, id)", {
  skip_no_ggplot2()
  fit <- make_fit()
  parts <- .gp_contrast_inputs(fit)
  ord <- order(parts$selected, parts$endpoint, parts$ids, method = "radix")
  # selected = c(u1=2,u2=1,u3=2,u4=3); endpoint = c(u1=2,u2=1,u3=1,u4=2)
  # => order u2(1,1), u3(2,1), u1(2,2), u4(3,2)
  expect_identical(parts$ids[ord], c("u2", "u3", "u1", "u4"))
  expect_identical(which(diff(parts$selected[ord]) != 0L), c(1L, 3L))
})

test_that("gp_plot_posterior_contrast() threshold outline toggles and is correct", {
  skip_no_ggplot2()
  # strict pairwise: six upper-triangle 0.9 cells >= threshold 0.8 (lambda 0.25)
  pw <- make_strict_pairwise(4)
  p_thr <- gp_plot_posterior_contrast(
    pw, selected_grades = c(1L, 2L, 3L, 4L), selected_lambda = 0.25, threshold = TRUE
  )
  p_no <- gp_plot_posterior_contrast(
    pw, selected_grades = c(1L, 2L, 3L, 4L), selected_lambda = 0.25, threshold = FALSE
  )

  expect_no_error(ggplot2::ggplot_build(p_thr))
  expect_no_error(ggplot2::ggplot_build(p_no))

  expect_true("GeomTile" %in% geom_classes(p_thr))
  expect_gt(length(p_thr$layers), length(p_no$layers))

  built <- ggplot2::ggplot_build(p_thr)
  tile_idx <- which(vapply(
    p_thr$layers,
    function(l) inherits(l$geom, "GeomTile") && !inherits(l$geom, "GeomRaster"),
    logical(1)
  ))
  expect_length(tile_idx, 1L)
  expect_equal(nrow(built$data[[tile_idx]]), 6L)
})

test_that("gp_plot_posterior_contrast() accepts gp_pairwise + explicit grades", {
  skip_no_ggplot2()
  fit <- make_fit()
  p <- gp_plot_posterior_contrast(
    fit$pairwise,
    selected_grades = fit$selected_grade$assignment$grade,
    endpoint_grades = c(2L, 1L, 1L, 2L),
    selected_lambda = 0.25
  )
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  expect_error(gp_plot_posterior_contrast(fit$pairwise), "selected_grades")
  expect_error(gp_plot_posterior_contrast(list(a = 1)), "gp_fit")
})

test_that("gp_plot_posterior_contrast() suppresses dense ticks for large J", {
  skip_no_ggplot2()
  n <- 30L
  ids <- paste0("u", sprintf("%02d", seq_len(n)))
  m <- matrix(0.5, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) m[i, j] <- if (i < j) 0.7 else 0.3
    }
  }
  pw <- new_gp_pairwise(
    ids = ids, matrix = m, power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = .gp_pairwise_diagonal, zero_floor = .gp_pairwise_zero_floor),
    source = list(stage = "posterior", rule = "outer_product", assumption = "one_level_independence"),
    control = gp_control(backend = "highs")
  )
  grades <- rep(1:5, each = 6)
  p <- gp_plot_posterior_contrast(pw, selected_grades = grades, selected_lambda = 0.25)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[1L]]), n^2)
  x_breaks <- built$layout$panel_params[[1L]]$x$breaks
  x_breaks <- x_breaks[is.finite(x_breaks)]
  expect_lt(length(x_breaks), n)
})

# =============================================================================
# autoplot
# =============================================================================

test_that("autoplot dispatches to the frontier plot", {
  skip_no_ggplot2()
  fr <- make_frontier()
  fit <- make_fit()

  p_fr <- ggplot2::autoplot(fr)
  expect_s3_class(p_fr, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p_fr))

  p_fit <- ggplot2::autoplot(fit)
  expect_s3_class(p_fit, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p_fit))
})

# =============================================================================
# vdiffr snapshot tests (gated -- only run if vdiffr is installed)
# =============================================================================

test_that("frontier figure matches its snapshot", {
  skip_if_not_installed("vdiffr")
  skip_no_ggplot2()
  vdiffr::expect_doppelganger("frontier-4grade-benchmarks", gp_plot_frontier(make_frontier()))
  vdiffr::expect_doppelganger(
    "frontier-no-benchmarks",
    gp_plot_frontier(make_frontier(benchmarks = FALSE), benchmarks = FALSE, annotate = FALSE)
  )
})

test_that("posterior-contrast figure matches its snapshot", {
  skip_if_not_installed("vdiffr")
  skip_no_ggplot2()
  vdiffr::expect_doppelganger("posterior-contrast-4unit", gp_plot_posterior_contrast(make_fit()))
  vdiffr::expect_doppelganger(
    "posterior-contrast-threshold",
    gp_plot_posterior_contrast(
      make_strict_pairwise(4),
      selected_grades = c(1L, 2L, 3L, 4L),
      selected_lambda = 0.25
    )
  )
})
