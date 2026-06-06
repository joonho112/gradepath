# =============================================================================
# plots-report-card.R -- KRW report-card + discordance figures
# -----------------------------------------------------------------------------
# Two more reader-facing display figures from Kline, Rose & Walters (2024),
# "A Discrimination Report Card" (KRW), rendered as ggplot2 objects, reusing the
# house style + ordinal grade palette from plots-frontier.R:
#   * gp_plot_report_card()  -- horizontal point-and-interval ("caterpillar")
#                               plot of the unit-level report card (KRW Fig 7 /
#                               8 / 14 / 15): posterior mean + credible interval
#                               per unit, units ordered by the endpoint/Condorcet
#                               rank (default), coloured by the selected grade.
#   * gp_plot_discordance()  -- lower-triangular heatmap of the conditional
#                               discordance matrix DR_{g,g'} (KRW Fig 10a / 17a),
#                               annotated with the DR value in each cell.
# plus autoplot.gp_report_card() dispatching to gp_plot_report_card().
#
# Design notes
#   * Grade colours come from scale_color_gradepath() / .gp_grade_colors() so a
#     grade's colour is identical across every gradepath figure.
#     theme_gradepath() is reused; nothing here redefines those helpers.
#   * The caterpillar is a single geom_pointrange(orientation = "y") layer:
#     x = posterior_mean, interval from the lower/upper ALREADY ON THE CARD
#     (never recomputed -- the payload is computed once and only displayed).
#     orientation = "y" is the ggplot2 >= 3.4 / 4.x idiom that supersedes the
#     deprecated geom_errorbarh().
#   * "Most-extreme unit on TOP": ggplot draws factor level 1 at the bottom, so
#     the y factor (keyed on the unique id) is releveled with the most-extreme
#     unit as the LAST level; scale_y_discrete() maps ids back to labels, so
#     duplicate labels still get one row each.
#   * The figure is welfare-NEUTRAL: grade 1's meaning flips firms <-> names, so
#     there is no "beats"/"wins"/"claim"/verdict text.
#   * The discordance heatmap keeps only the lower triangle (the matrix is NA on
#     the diagonal and upper triangle by construction) and draws it on discrete
#     grade-label axes: the more-extreme grade across the columns (grade 1 on the
#     left), the less-extreme grade up the rows (grade 2 at the top). Annotated
#     cells auto-contrast (white text on dark fills) by relative luminance.
#   * Every gp_plot_*() RETURNS a ggplot and renders nothing; ASCII-only source;
#     aes() columns are referenced through the rlang .data pronoun (imported
#     package-wide in plots-frontier.R).
# =============================================================================

# ---------------------------------------------------------------------------
# gp_plot_report_card() -- KRW Fig 7 / 8 / 14 / 15
# ---------------------------------------------------------------------------

# Coerce any accepted input to a gp_report_card (a gp_fit carries one in
# $report_card; a gp_report_card is returned unchanged).
.gp_as_report_card <- function(x) {
  if (inherits(x, "gp_report_card")) {
    return(x)
  }
  if (inherits(x, "gp_fit")) {
    rc <- x$report_card
    if (!inherits(rc, "gp_report_card")) {
      .gradepath_abort("`gp_plot_report_card()` got a `gp_fit` without a `$report_card`.")
    }
    return(rc)
  }
  .gradepath_abort(
    "`gp_plot_report_card()` accepts a `gp_report_card` or a `gp_fit`; got <%s>.",
    paste(class(x), collapse = "/")
  )
}

# Row permutation (FIRST element is the most-extreme unit, drawn at the TOP).
# "rank" = endpoint/Condorcet sort_rank; "mean" = largest posterior_mean first;
# "grade" = grade block then sort_rank; "label" = alphabetical. Radix sort keeps
# the order locale-independent.
.gp_report_card_row_order <- function(tbl, order) {
  switch(order,
    rank  = order(as.integer(tbl$sort_rank), method = "radix"),
    mean  = order(-as.numeric(tbl$posterior_mean), as.integer(tbl$sort_rank), method = "radix"),
    grade = order(as.integer(tbl$grade), as.integer(tbl$sort_rank), method = "radix"),
    label = order(as.character(tbl$label), as.character(tbl$id), method = "radix")
  )
}

#' Plot the KRW report card (horizontal point-and-interval)
#'
#' `gp_plot_report_card()` draws the Kline-Rose-Walters unit-level report card
#' (their Figures 7, 8, 14, 15) as a horizontal point-and-interval
#' ("caterpillar") plot: each unit's posterior mean (point) with its stored
#' credible interval (`lower`, `upper`), one unit per row, ordered so the
#' most-extreme unit is on top, and coloured by the unit's selected grade using
#' the shared ordinal palette ([scale_color_gradepath()]). The interval already
#' on the card is displayed as-is and never recomputed. The figure is neutral: it
#' carries no welfare-verdict text, because grade 1's welfare meaning flips
#' between firms and names. Pass it a [gp_report_card()] (or the `gp_fit` that
#' carries one).
#'
#' @param x A `gp_report_card` (primary). A `gp_fit` is also accepted, in which
#'   case its `$report_card` is used.
#' @param order Row order, most-extreme unit on top. `"rank"` (default) uses the
#'   endpoint/Condorcet `sort_rank` (the canonical endpoint-sorted figure);
#'   `"mean"` orders by `posterior_mean`; `"grade"` groups by grade then
#'   `sort_rank`; `"label"` is alphabetical.
#' @param ci_level Optional number in `(0, 1)` used ONLY to label the caption,
#'   e.g. `0.9` -> "90% credible interval". `NULL` (default) reads the level from
#'   `x$control$interval_level` (default `0.90`); the interval itself is never
#'   recomputed.
#' @param max_rows Optional positive integer keeping only the top `max_rows`
#'   units by the chosen order, with a "Showing N of M" subtitle (never silent).
#'   `NULL` (default) shows all units.
#' @param zero_line Dashed reference line at `x = 0`: `"auto"` (default) draws it
#'   only when the intervals straddle zero (`min(lower) < 0 < max(upper)`);
#'   `"on"` always; `"off"` never.
#' @param ... Reserved for future arguments; currently unused.
#'
#' @return A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
#'   printed.
#'
#' @examples
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' card <- get_report_card(fit)
#'
#' gp_plot_report_card(card)                  # endpoint-sorted, coloured by grade
#' gp_plot_report_card(card, order = "grade") # grouped by grade
#' gp_plot_report_card(card, max_rows = 10)   # top 10 only
#' gp_plot_report_card(fit)                   # convenience: a gp_fit
#'
#' @seealso [gp_plot_discordance()], [theme_gradepath()],
#'   [scale_color_gradepath()], [gp_report_card()]
#' @family gradepath-plots
#' @export
gp_plot_report_card <- function(x,
                                order = "rank",
                                ci_level = NULL,
                                max_rows = NULL,
                                zero_line = "auto",
                                ...) {
  card <- .gp_as_report_card(x)
  order <- match.arg(order, c("rank", "mean", "grade", "label"))
  zero_line <- match.arg(as.character(zero_line), c("auto", "on", "off"))

  tbl <- card$table
  n_total <- nrow(tbl)
  if (n_total == 0L) {
    .gradepath_abort("`gp_plot_report_card()` requires a non-empty report card.")
  }

  # Order most-extreme-first, then optional top-N truncation.
  tbl <- tbl[.gp_report_card_row_order(tbl, order), , drop = FALSE]
  subtitle <- NULL
  if (!is.null(max_rows)) {
    if (length(max_rows) != 1L || !is.numeric(max_rows) || is.na(max_rows) ||
      max_rows < 1 || abs(max_rows - round(max_rows)) > 1e-9) {
      .gradepath_abort("`max_rows` must be a single positive integer or NULL.")
    }
    max_rows <- as.integer(round(max_rows))
    if (max_rows < n_total) {
      tbl <- tbl[seq_len(max_rows), , drop = FALSE]
      subtitle <- sprintf("Showing %d of %d units (top by %s order)", nrow(tbl), n_total, order)
    }
  }

  # Grade levels present on the FULL card fix the palette length, so colours are
  # stable even when max_rows hides some grades.
  grade_levels <- sort(unique(as.integer(card$table$grade)))

  df <- data.frame(
    id = as.character(tbl$id),
    label = as.character(tbl$label),
    posterior_mean = as.numeric(tbl$posterior_mean),
    lower = as.numeric(tbl$lower),
    upper = as.numeric(tbl$upper),
    grade = factor(as.integer(tbl$grade), levels = grade_levels),
    stringsAsFactors = FALSE
  )
  # y keyed on the unique id; level order reversed so the first row (most
  # extreme) is the LAST level -> drawn at the TOP. Labels shown via the lookup.
  df$id <- factor(df$id, levels = rev(df$id))
  axis_labels <- stats::setNames(df$label, as.character(df$id))

  # Zero-line decision uses the FULL card so truncation cannot flip it. The
  # report-card validator permits NA interval bounds, so na.rm guards them.
  draw_zero <- switch(zero_line,
    on = TRUE,
    off = FALSE,
    auto = isTRUE(min(card$table$lower, na.rm = TRUE) < 0 &&
      max(card$table$upper, na.rm = TRUE) > 0)
  )

  # Caption: the stored credible-interval level (display only; not recomputed).
  level <- if (is.null(ci_level)) {
    lv <- tryCatch(card$control$interval_level, error = function(e) NULL)
    if (is.null(lv) || !is.finite(lv)) 0.90 else as.numeric(lv)
  } else {
    lv <- suppressWarnings(as.numeric(ci_level))
    if (length(lv) != 1L || is.na(lv) || lv <= 0 || lv >= 1) {
      .gradepath_abort("`ci_level` must be a single number in (0, 1) or NULL.")
    }
    lv
  }
  caption <- sprintf("%g%% credible interval", round(100 * level, 6))

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["posterior_mean"]], y = .data[["id"]],
      xmin = .data[["lower"]], xmax = .data[["upper"]],
      colour = .data[["grade"]]
    )
  )

  # Zero line drawn first so the intervals sit on top of it.
  if (isTRUE(draw_zero)) {
    p <- p + ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4
    )
  }

  p +
    ggplot2::geom_pointrange(orientation = "y", size = 0.275, linewidth = 0.5) +
    scale_color_gradepath(n_grades = length(grade_levels)) +
    ggplot2::scale_y_discrete(labels = axis_labels) +
    ggplot2::labs(
      x = "Posterior mean",
      y = NULL,
      colour = "Grade",
      title = "Report card",
      subtitle = subtitle,
      caption = caption
    ) +
    theme_gradepath() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey92", linewidth = 0.3)
    )
}

# ---------------------------------------------------------------------------
# gp_plot_discordance() -- KRW Fig 10a / 17a
# ---------------------------------------------------------------------------

# Coerce any accepted input to a k x k conditional-DR matrix (lower triangle
# populated, upper triangle + diagonal NA, dimnames grade_1..grade_k -- or a
# validated bare numeric matrix).
.gp_as_dr_matrix <- function(x) {
  if (inherits(x, "gp_frontier")) {
    return(x$dr_matrix)
  }
  if (inherits(x, "gp_fit")) {
    return(gp_frontier(x$grade_path, pairwise = x$pairwise)$dr_matrix)
  }
  if (is.matrix(x) && is.numeric(x)) {
    if (nrow(x) != ncol(x)) {
      .gradepath_abort("`gp_plot_discordance()` matrix input must be square.")
    }
    return(x)
  }
  .gradepath_abort(
    "`gp_plot_discordance()` accepts a `gp_frontier`, `gp_fit`, or a numeric DR matrix; got <%s>.",
    paste(class(x), collapse = "/")
  )
}

# Plain ASCII grade labels from dr_matrix dimnames (strip "grade_"; fall back to
# the integer index when dimnames are absent).
.gp_dr_grade_labels <- function(mat) {
  nm <- rownames(mat)
  if (is.null(nm)) nm <- colnames(mat)
  if (is.null(nm)) {
    return(paste("Grade", seq_len(nrow(mat))))
  }
  paste("Grade", sub("^grade_", "", nm))
}

# White-or-dark text per cell so annotations stay legible at any DR scale.
# Evaluate the EXACT fill that scale_fill_distiller(palette = "Blues",
# direction = 1) renders for each cell's rescaled DR value -- the same 7-stop
# Brewer "Blues" ramp interpolated in Lab space that distiller builds -- then
# pick the higher-contrast of white / grey10 by relative luminance. Using the
# real fill (not an approximation) keeps the text->fill contrast decision
# correct across the whole DR range.
.gp_dr_text_color <- function(values, vmax) {
  if (length(values) == 0L) {
    return(character(0))
  }
  rng <- if (is.finite(vmax) && vmax > 0) vmax else 1
  t <- pmin(pmax(values / rng, 0), 1)
  fill_pal <- scales::gradient_n_pal(
    scales::brewer_pal(type = "seq", palette = "Blues", direction = 1)(7),
    space = "Lab"
  )
  rgb_mat <- grDevices::col2rgb(fill_pal(t)) / 255 # 3 x n
  lum <- 0.299 * rgb_mat[1, ] + 0.587 * rgb_mat[2, ] + 0.114 * rgb_mat[3, ]
  ifelse(lum < 0.5, "white", "grey10")
}

#' Plot the KRW conditional-discordance heatmap (lower triangle)
#'
#' `gp_plot_discordance()` draws the Kline-Rose-Walters conditional discordance
#' matrix (their Figure 10a / 17a) as a lower-triangular heatmap. The cell for a
#' grade pair is the conditional discordance rate \eqn{DR_{g,g'}} between the
#' more-extreme grade `g` (columns, grade 1 on the left) and the less-extreme
#' grade `g'` (rows, grade 2 at the top): the posterior share of disagreeing
#' rankings among cross-grade pairs. Higher DR is darker; with `annotate = TRUE`
#' the value is written in each cell with auto-contrasting text. Pass it a
#' [gp_frontier()] (or, for convenience, the `gp_fit` it was built from, or a
#' bare DR matrix).
#'
#' @param x A `gp_frontier` (reads `x$dr_matrix`), a `gp_fit` (its frontier is
#'   rebuilt with [gp_frontier()] from `x$grade_path` and `x$pairwise`), or a
#'   numeric conditional-DR matrix (square; lower triangle populated, upper
#'   triangle and diagonal `NA`).
#' @param annotate Logical; write each cell's DR value (`"%.3f"`) with an
#'   auto-contrasting colour. Default `TRUE`.
#' @param ... Reserved for future arguments; currently unused.
#'
#' @return A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
#'   printed.
#'
#' @examples
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#'
#' # Build a frontier from the fit's grade path and pairwise.
#' fr <- gp_frontier(fit$grade_path, pairwise = get_pairwise(fit))
#'
#' gp_plot_discordance(fr)                   # annotated lower-triangular heatmap
#' gp_plot_discordance(fr, annotate = FALSE) # heatmap without cell values
#' gp_plot_discordance(fit)                  # convenience: a gp_fit
#' gp_plot_discordance(fr$dr_matrix)         # a bare dr_matrix
#'
#' @seealso [gp_plot_report_card()], [theme_gradepath()], [gp_frontier()]
#' @family gradepath-plots
#' @export
gp_plot_discordance <- function(x, annotate = TRUE, ...) {
  mat <- .gp_as_dr_matrix(x)
  k <- nrow(mat)
  labels <- .gp_dr_grade_labels(mat)

  base_labs <- ggplot2::labs(
    x = "More-extreme grade",
    y = "Less-extreme grade",
    title = "Conditional discordance (DR g,g')"
  )

  # Degenerate frontier (< 2 grades): no grade pairs -> a valid placeholder.
  if (k < 2L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No grade boundaries", colour = "grey40") +
        base_labs +
        theme_gradepath()
    )
  }

  # Lower triangle only (row > col, finite): row = less-extreme grade,
  # col = more-extreme grade.
  rc <- which(lower.tri(mat) & is.finite(mat), arr.ind = TRUE)
  cells <- data.frame(
    row = as.integer(rc[, "row"]),
    col = as.integer(rc[, "col"]),
    dr = as.numeric(mat[lower.tri(mat) & is.finite(mat)]),
    stringsAsFactors = FALSE
  )

  # Discrete grade-label axes: more-extreme grade across columns (grade 1 left),
  # less-extreme grade up the rows (grade 2 at the TOP via reversed y levels).
  col_levels <- labels[sort(unique(cells$col))]
  row_levels <- rev(labels[sort(unique(cells$row))])
  cells$more <- factor(labels[cells$col], levels = col_levels)
  cells$less <- factor(labels[cells$row], levels = row_levels)

  dr_max <- if (nrow(cells) > 0L) max(cells$dr) else 1

  p <- ggplot2::ggplot(
    cells,
    ggplot2::aes(x = .data[["more"]], y = .data[["less"]], fill = .data[["dr"]])
  ) +
    ggplot2::geom_tile(colour = "grey90", linewidth = 0.3) +
    ggplot2::scale_fill_distiller(
      name = "DR",
      palette = "Blues",
      direction = 1,
      limits = c(0, dr_max)
    )

  if (isTRUE(annotate) && nrow(cells) > 0L) {
    cells$txt_colour <- .gp_dr_text_color(cells$dr, dr_max)
    cells$txt <- sprintf("%.3f", cells$dr)
    p <- p + ggplot2::geom_text(
      data = cells,
      ggplot2::aes(x = .data[["more"]], y = .data[["less"]], label = .data[["txt"]]),
      colour = cells$txt_colour, size = 3, inherit.aes = FALSE
    )
  }

  p +
    ggplot2::coord_equal() +
    base_labs +
    theme_gradepath() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(2.5, "lines"),
      legend.key.height = ggplot2::unit(0.6, "lines")
    )
}

# ---------------------------------------------------------------------------
# autoplot() method for the report card (dispatch to gp_plot_report_card())
# ---------------------------------------------------------------------------

#' Autoplot method for a gradepath report card
#'
#' A thin [ggplot2::autoplot()] method so `autoplot(card)` (a `gp_report_card`)
#' draws the report-card caterpillar via [gp_plot_report_card()] without
#' attaching ggplot2. (There is no discordance autoplot method:
#' `autoplot.gp_frontier` already maps to the frontier figure.)
#'
#' @param object A `gp_report_card`.
#' @param ... Further arguments forwarded to [gp_plot_report_card()].
#'
#' @return A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
#'   printed.
#'
#' @examples
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' autoplot(get_report_card(fit))   # the report-card caterpillar for a gp_fit
#'
#' @seealso [gp_plot_report_card()], [gp_report_card()]
#' @family gradepath-plots
#' @export
autoplot.gp_report_card <- function(object, ...) {
  gp_plot_report_card(object, ...)
}
