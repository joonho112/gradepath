# =============================================================================
# plots-frontier.R -- KRW frontier + posterior-contrast figures
# -----------------------------------------------------------------------------
# Two reader-facing display figures from Kline, Rose & Walters (2024),
# "A Discrimination Report Card" (KRW), rendered as ggplot2 objects:
#   * gp_plot_frontier()           -- information-reliability frontier
#                                     (KRW Fig 6b / 13b): the (1 - DR, tau_bar)
#                                     optimal-grade path with the selected-lambda
#                                     slice and the five dominated naive rules.
#   * gp_plot_posterior_contrast() -- grade-sorted posterior-contrast heatmap
#                                     (KRW Fig 5a / 12a): the J x J pairwise
#                                     matrix Pi reordered by grade blocks.
# plus theme_gradepath() (exported KRW house-style theme), a reusable ordinal
# grade palette (scale_color_gradepath() / scale_fill_gradepath()), and
# autoplot() methods.
#
# Design notes
#   * Palette: a perceptually uniform, colour-blind-safe viridis ramp via base
#     grDevices::hcl.colors() (no extra hard dependency); grade 1 (the
#     most-extreme theta block) is the dark/high-contrast end.
#   * Heatmap fill: sequential Brewer "Blues" on [0, 1] -- Pi is a probability,
#     so a single-hue sequential ramp keeps Pi = 0.5 visually neutral and the
#     >= threshold cells dark, matching the published figure.
#   * Reorder key: (selected grade, endpoint grade, id) with order(method =
#     "radix") -- locale-independent, so the figure is reproducible under any
#     LC_COLLATE (the same fix the cache hashing uses).
#   * Every gp_plot_*() RETURNS a ggplot and renders nothing itself; the
#     KRW-faithful theme/palette/labels are applied inside the function but stay
#     overridable (p + theme_*() / + scale_*).
#   * Source is ASCII-only (R CMD check cleanliness); math labels use plotmath
#     expressions, not Unicode glyphs.
# =============================================================================

#' gradepath plotting internals
#'
#' Import anchor for the gradepath figure layer. The plotting code uses many
#' ggplot2 verbs, so the whole namespace is imported; data-frame columns are
#' referenced through the `rlang` `.data` pronoun inside `aes()` to avoid
#' R CMD check "no visible binding for global variable" notes.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @keywords internal
#' @name gradepath-plots
NULL

# ---------------------------------------------------------------------------
# Reusable ordinal grade palette (shared by the frontier and contrast figures)
# ---------------------------------------------------------------------------

#' Ordinal colour values for grade levels
#'
#' Internal helper returning `n_grades` colours along a perceptually uniform,
#' colour-blind-safe viridis ramp (via [grDevices::hcl.colors()]) so that grade
#' 1 (the most-extreme theta block) is the dark, high-contrast end of the scale.
#' Shared by every gradepath figure so grade colours stay identical across
#' plots.
#'
#' @param n_grades Positive integer: the number of grade levels.
#'
#' @return A character vector of `n_grades` hex colours, named `"1"` ..
#'   `"n_grades"`, ordered from grade 1 to grade `n_grades`.
#'
#' @keywords internal
.gp_grade_colors <- function(n_grades) {
  n_grades <- as.integer(n_grades)
  if (length(n_grades) != 1L || is.na(n_grades) || n_grades < 1L) {
    .gradepath_abort("`n_grades` must be a single positive integer.")
  }
  cols <- grDevices::hcl.colors(max(n_grades, 2L), palette = "viridis")
  cols <- cols[seq_len(n_grades)]
  stats::setNames(cols, as.character(seq_len(n_grades)))
}

#' gradepath ordinal grade colour and fill scales
#'
#' Discrete ggplot2 scales that map grade levels to the shared gradepath ordinal
#' palette, so your own plots and the gradepath report-card figures use identical
#' grade colours. The palette is a perceptually uniform, colour-blind-safe
#' viridis ramp; grade 1 (the most-extreme theta block) is the dark,
#' high-contrast end. `scale_colour_gradepath()` is a spelling alias of
#' `scale_color_gradepath()`, and `scale_fill_gradepath()` is the fill-aesthetic
#' counterpart. Add one to a ggplot whose `colour`/`fill` aesthetic is a grade
#' factor.
#'
#' @param n_grades Optional positive integer fixing the palette length. When
#'   `NULL` (the default) the palette adapts to the number of grade levels in the
#'   data.
#' @param ... Further arguments forwarded to [ggplot2::discrete_scale()] (for
#'   example `name`, `labels`, `guide`).
#'
#' @return A ggplot2 discrete scale (a `Scale`/`ggproto` object) to add to a plot
#'   with `+`.
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x = 1:3, y = 1:3, grade = factor(1:3))
#' ggplot(df, aes(x, y, colour = grade)) +
#'   geom_point(size = 4) +
#'   scale_color_gradepath() +
#'   theme_gradepath()
#'
#' # The fill counterpart, with a fixed palette length:
#' ggplot(df, aes(x, y, fill = grade)) +
#'   geom_tile() +
#'   scale_fill_gradepath(n_grades = 3)
#'
#' @seealso [theme_gradepath()], [gp_plot_report_card()], [gp_plot_frontier()]
#' @family gradepath-plots
#' @export
scale_color_gradepath <- function(n_grades = NULL, ...) {
  pal <- function(n) unname(.gp_grade_colors(if (is.null(n_grades)) n else n_grades))
  ggplot2::discrete_scale(aesthetics = "colour", palette = pal, ...)
}

#' @rdname scale_color_gradepath
#' @export
scale_colour_gradepath <- scale_color_gradepath

#' @rdname scale_color_gradepath
#' @export
scale_fill_gradepath <- function(n_grades = NULL, ...) {
  pal <- function(n) unname(.gp_grade_colors(if (is.null(n_grades)) n else n_grades))
  ggplot2::discrete_scale(aesthetics = "fill", palette = pal, ...)
}

# ---------------------------------------------------------------------------
# KRW house-style theme
# ---------------------------------------------------------------------------

#' gradepath ggplot2 theme (KRW house style)
#'
#' A clean ggplot2 theme matching the Kline-Rose-Walters figure house style:
#' built on [ggplot2::theme_minimal()], with the top and right spines removed, a
#' light horizontal grid, no vertical grid, and the legend at the bottom. It is
#' the default theme inside every gradepath plot verb, and is exported so you can
#' re-apply or extend it on any ggplot with `p + theme_gradepath()`.
#'
#' @param base_size Base font size in points. Default `11`.
#' @param base_family Base font family. Default `""` (the device default).
#' @param ... Further arguments forwarded to [ggplot2::theme_minimal()].
#'
#' @return A ggplot2 theme object (class `c("theme", "gg")`) to add to a plot
#'   with `+`.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_gradepath()
#'
#' @seealso [scale_color_gradepath()], [gp_plot_frontier()],
#'   [gp_plot_posterior_contrast()]
#' @family gradepath-plots
#' @export
theme_gradepath <- function(base_size = 11, base_family = "", ...) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family, ...) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      axis.line.x = ggplot2::element_line(colour = "grey20", linewidth = 0.4),
      axis.line.y = ggplot2::element_line(colour = "grey20", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "grey20", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(colour = "grey30", hjust = 0),
      plot.caption = ggplot2::element_text(colour = "grey40", hjust = 1),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.9, "lines"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# ---------------------------------------------------------------------------
# Frontier coercion: accept gp_frontier (primary), gp_fit, gp_grade_path
# ---------------------------------------------------------------------------

# Resolve any accepted input to a gp_frontier. For a gp_fit / gp_grade_path the
# frontier is rebuilt from the stored pairwise via gp_frontier(); benchmarks are
# NOT force-derived (they appear only if the resulting object already carries
# them, e.g. the input was a gp_frontier built with benchmark_scores).
.gp_as_frontier <- function(x, pairwise = NULL) {
  if (inherits(x, "gp_frontier")) {
    return(x)
  }
  if (inherits(x, "gp_fit")) {
    pw <- if (is.null(pairwise)) x$pairwise else pairwise
    return(gp_frontier(x$grade_path, pairwise = pw))
  }
  if (inherits(x, "gp_grade_path")) {
    return(gp_frontier(x, pairwise = pairwise))
  }
  .gradepath_abort(
    "`gp_plot_frontier()` accepts a `gp_frontier`, `gp_fit`, or `gp_grade_path`; got <%s>.",
    paste(class(x), collapse = "/")
  )
}

# Human-readable labels for the five naive benchmark rules.
.gp_benchmark_labels <- function(benchmark) {
  lut <- c(
    raw_estimate_rank = "Raw rank",
    posterior_mean_rank = "Posterior-mean rank",
    linear_shrinkage_rank = "Linear-shrinkage rank",
    posterior_mean_decile = "Posterior-mean decile",
    posterior_mean_quartile = "Posterior-mean quartile"
  )
  out <- unname(lut[as.character(benchmark)])
  out[is.na(out)] <- as.character(benchmark)[is.na(out)]
  out
}

# ---------------------------------------------------------------------------
# gp_plot_frontier() -- KRW Fig 6b / 13b
# ---------------------------------------------------------------------------

#' Plot the KRW information-reliability frontier
#'
#' `gp_plot_frontier()` draws the Kline-Rose-Walters information-reliability
#' frontier (their Figure 6b / 13b): the optimal-grade path traced in reliability
#' `1 - DR` (x-axis) against the posterior expected Kendall agreement
#' \eqn{\bar{\tau}} (y-axis). Each frontier-table row is one point and the points
#' are connected as the Pareto path. The selected-lambda slice is marked with a
#' highlighted point, and the five naive benchmark rules are overlaid as labelled
#' shapes when the frontier carries them. Feed it a [gp_frontier()] (or, for
#' convenience, the `gp_fit` or `gp_grade_path` it was built from).
#'
#' @param x A `gp_frontier` (primary). For convenience a `gp_fit` or
#'   `gp_grade_path` is also accepted and converted with [gp_frontier()];
#'   benchmarks are shown only if already present (they are never force-derived
#'   from a bare fit or path).
#' @param highlight Which slice to emphasise: `"selected"` (default) marks the
#'   selected lambda; `"all"` marks every path point; a numeric value marks the
#'   nearest path lambda.
#' @param benchmarks Logical; overlay the naive benchmark rules as labelled
#'   shapes when they are present on the frontier. Default `TRUE`.
#' @param annotate Logical; report the highlighted `(1 - DR, tau_bar)` point in
#'   the subtitle. Default `TRUE`.
#' @param pairwise Optional `gp_pairwise`, used only when `x` is a `gp_grade_path`
#'   (or to override a fit's stored pairwise) so the frontier can be rebuilt.
#' @param ... Reserved for future arguments; currently unused.
#'
#' @return A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
#'   printed and stays fully overridable (add `+ theme_*()` / `+ scale_*`).
#'
#' @examples
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#'
#' # Build a frontier with the three shrinkage benchmark scores on the fit's ids.
#' post <- get_posterior(fit)
#' ids <- get_pairwise(fit)$ids
#' fr <- gp_frontier(
#'   fit$grade_path,
#'   pairwise = get_pairwise(fit),
#'   benchmark_scores = list(
#'     raw_estimate     = setNames(post$estimate, ids),
#'     posterior_mean   = setNames(post$posterior_mean, ids),
#'     linear_shrinkage = setNames(0.5 * post$estimate, ids)
#'   )
#' )
#'
#' gp_plot_frontier(fr)                    # selected slice + benchmarks
#' gp_plot_frontier(fr, highlight = "all") # mark every path point
#' gp_plot_frontier(fit)                   # convenience: a gp_fit
#'
#' @seealso [gp_plot_posterior_contrast()], [gp_plot_discordance()],
#'   [theme_gradepath()], [gp_frontier()]
#' @family gradepath-plots
#' @export
gp_plot_frontier <- function(x,
                             highlight = "selected",
                             benchmarks = TRUE,
                             annotate = TRUE,
                             pairwise = NULL,
                             ...) {
  frontier <- .gp_as_frontier(x, pairwise = pairwise)

  tbl <- frontier$table
  selected_lambda <- frontier$selection$selected_lambda

  path_df <- data.frame(
    lambda = as.numeric(tbl$lambda),
    reliability = as.numeric(tbl$reliability),
    tau_bar = as.numeric(tbl$tau_bar),
    stringsAsFactors = FALSE
  )
  # Order along the curve by reliability so the connecting line traces the
  # Pareto path monotonically.
  path_df <- path_df[order(path_df$reliability, path_df$tau_bar), , drop = FALSE]

  # Which path point(s) to highlight.
  if (length(highlight) != 1L || is.na(highlight)) {
    .gradepath_abort(
      "`highlight` must be a single non-missing value: \"selected\", \"all\", or a numeric lambda."
    )
  }
  if (is.numeric(highlight)) {
    hl_df <- path_df[which.min(abs(path_df$lambda - highlight)), , drop = FALSE]
  } else {
    highlight <- match.arg(as.character(highlight), c("selected", "all"))
    hl_df <- if (identical(highlight, "all")) {
      path_df
    } else {
      path_df[which.min(abs(path_df$lambda - selected_lambda)), , drop = FALSE]
    }
  }

  # One tidy point frame with a `series` factor (legend order: the optimal-grade
  # path, the selected slice, then the up-to-5 benchmark rules). The connecting
  # Pareto line stays a separate unmapped grey geom_line under the points, so it
  # never reaches the legend. Under highlight = "all" every path point is the
  # selected slice, so those rows move to the selected series (no double-marking).
  lvl_path <- "Optimal-grade path"
  lvl_selected <- "Selected slice"

  is_hl <- path_df$lambda %in% hl_df$lambda
  pts_path <- path_df[!is_hl, c("reliability", "tau_bar"), drop = FALSE]
  pts_sel <- path_df[is_hl, c("reliability", "tau_bar"), drop = FALSE]

  point_frames <- list()
  if (nrow(pts_path) > 0L) {
    point_frames[[length(point_frames) + 1L]] <- data.frame(
      reliability = pts_path$reliability, tau_bar = pts_path$tau_bar,
      series = lvl_path, stringsAsFactors = FALSE
    )
  }
  if (nrow(pts_sel) > 0L) {
    point_frames[[length(point_frames) + 1L]] <- data.frame(
      reliability = pts_sel$reliability, tau_bar = pts_sel$tau_bar,
      series = lvl_selected, stringsAsFactors = FALSE
    )
  }

  # Benchmark points (optional): each rule is its own series -> a legend swatch.
  bench <- frontier$benchmarks
  show_bench <- isTRUE(benchmarks) && is.data.frame(bench) && nrow(bench) > 0L
  bench_levels <- character(0)
  if (show_bench) {
    bench_labels <- .gp_benchmark_labels(bench$benchmark)
    bench_levels <- unique(bench_labels)
    point_frames[[length(point_frames) + 1L]] <- data.frame(
      reliability = as.numeric(bench$reliability),
      tau_bar = as.numeric(bench$tau_bar),
      series = bench_labels, stringsAsFactors = FALSE
    )
  }

  point_df <- do.call(rbind, point_frames)

  # Fixed, legend-ordered levels; keep only those actually drawn so the manual
  # scales (values keyed to these levels) and the legend agree when benchmarks
  # are absent.
  all_levels <- c(lvl_path, lvl_selected, bench_levels)
  present <- all_levels[all_levels %in% point_df$series]
  point_df$series <- factor(point_df$series, levels = present)

  # Colour + shape per series, keyed by name to the present levels. Path neutral
  # dark, selected a strong highlight, the benchmarks a colour-blind-aware
  # qualitative ramp (base hcl.colors, no extra dependency) + solid shapes that
  # also separate in grayscale.
  col_map <- stats::setNames(character(length(present)), present)
  shp_map <- stats::setNames(integer(length(present)), present)
  col_map[[lvl_path]] <- "grey30"
  shp_map[[lvl_path]] <- 16L
  if (lvl_selected %in% present) {
    col_map[[lvl_selected]] <- "firebrick"
    shp_map[[lvl_selected]] <- 18L
  }
  if (length(bench_levels) > 0L) {
    bench_cols <- grDevices::hcl.colors(max(length(bench_levels), 2L), palette = "Dark 3")
    bench_cols <- bench_cols[seq_along(bench_levels)]
    bench_shapes <- rep_len(c(17L, 15L, 8L, 7L, 9L, 10L), length(bench_levels))
    for (i in seq_along(bench_levels)) {
      col_map[[bench_levels[[i]]]] <- bench_cols[[i]]
      shp_map[[bench_levels[[i]]]] <- bench_shapes[[i]]
    }
  }
  col_values <- col_map[present]
  shp_values <- shp_map[present]
  legend_title <- "Grade rule"

  # Assemble: unmapped grey Pareto line, then the colour+shape mapped points. The
  # two manual scales share identical name/breaks/limits/labels so ggplot2 merges
  # them into ONE legend; guides() forces 2 rows at the bottom.
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = path_df,
      ggplot2::aes(x = .data[["reliability"]], y = .data[["tau_bar"]]),
      colour = "grey40", linewidth = 0.5, inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(
        x = .data[["reliability"]], y = .data[["tau_bar"]],
        colour = .data[["series"]], shape = .data[["series"]]
      ),
      size = 3, stroke = 0.9, inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(
      name = legend_title, values = col_values,
      breaks = present, limits = present, labels = present, drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = legend_title, values = shp_values,
      breaks = present, limits = present, labels = present, drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    )

  # Selected-slice coordinates in the subtitle (clean, no in-plot text labels).
  lab_subtitle <- if (isTRUE(annotate)) {
    sel_row <- path_df[which.min(abs(path_df$lambda - selected_lambda)), , drop = FALSE]
    sprintf(
      "Selected lambda = %s:  1 - DR = %.3f,  tau-bar = %.3f",
      format(selected_lambda, trim = TRUE),
      sel_row$reliability[[1L]], sel_row$tau_bar[[1L]]
    )
  } else {
    NULL
  }

  p +
    ggplot2::labs(
      x = "1 - DR (reliability)",
      y = expression(bar(tau)),
      title = "Information-reliability frontier",
      subtitle = lab_subtitle,
      colour = legend_title,
      shape = legend_title
    ) +
    theme_gradepath() +
    ggplot2::theme(legend.position = "bottom")
}

# ---------------------------------------------------------------------------
# gp_plot_posterior_contrast() -- KRW Fig 5a / 12a
# ---------------------------------------------------------------------------

#' Plot the grade-sorted posterior-contrast heatmap
#'
#' `gp_plot_posterior_contrast()` draws the Kline-Rose-Walters grade-sorted
#' posterior-contrast heatmap (their Figure 5a / 12a): the J x J pairwise
#' outranking matrix \eqn{\Pi} (\eqn{\pi_{ij}}, the posterior probability that
#' unit `i` is more extreme than unit `j`), with rows and columns reordered by
#' `(selected grade, endpoint grade, id)` so the grade blocks are contiguous.
#' Black segments separate the grade blocks; when `threshold = TRUE`, cells whose
#' contrast reaches the different-grade threshold \eqn{1 / (1 + \lambda)} are
#' outlined in red. Pass it a `gp_fit` (it reads the pairwise and grades from the
#' fit).
#'
#' @param x A `gp_fit` (primary): uses `x$pairwise$matrix`, `x$selected_grade`,
#'   and `x$grade_path`. A `gp_pairwise` is also accepted together with explicit
#'   `selected_grades` (and optionally `endpoint_grades`) vectors plus
#'   `selected_lambda`.
#' @param threshold Logical; outline cells reaching the different-grade threshold
#'   `1 / (1 + selected_lambda)`. Default `TRUE`.
#' @param selected_grades Integer grade vector (one per unit, in `ids` order),
#'   required only when `x` is a bare `gp_pairwise`.
#' @param endpoint_grades Integer grade vector used as the secondary sort key;
#'   defaults to `selected_grades` when omitted. Only used with a `gp_pairwise`.
#' @param selected_lambda Numeric selected lambda used for the threshold,
#'   required only when `x` is a bare `gp_pairwise` (read from the fit otherwise).
#' @param ... Reserved for future arguments; currently unused.
#'
#' @return A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
#'   printed.
#'
#' @examples
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#'
#' gp_plot_posterior_contrast(fit)                    # grade-sorted Pi heatmap
#' gp_plot_posterior_contrast(fit, threshold = FALSE) # no different-grade outline
#'
#' @seealso [gp_plot_frontier()], [theme_gradepath()], [gp_pairwise()]
#' @family gradepath-plots
#' @export
gp_plot_posterior_contrast <- function(x,
                                       threshold = TRUE,
                                       selected_grades = NULL,
                                       endpoint_grades = NULL,
                                       selected_lambda = NULL,
                                       ...) {
  parts <- .gp_contrast_inputs(
    x,
    selected_grades = selected_grades,
    endpoint_grades = endpoint_grades,
    selected_lambda = selected_lambda
  )
  ids <- parts$ids
  n <- length(ids)

  # Row/column order: (selected grade, endpoint grade, id). Radix is
  # locale-independent and so reproducible under any LC_COLLATE.
  ord <- order(parts$selected, parts$endpoint, ids, method = "radix")
  ord_ids <- ids[ord]
  ord_sel <- parts$selected[ord]
  mat <- parts$matrix[ord, ord, drop = FALSE]

  # Long-format raster data. Row index 1 is drawn at the TOP (y flipped) so the
  # matrix reads like the printed Pi with grade 1 in the top-left block.
  long <- data.frame(
    row = rep(seq_len(n), times = n),
    col = rep(seq_len(n), each = n),
    pi = as.numeric(mat),
    stringsAsFactors = FALSE
  )
  long$y <- n - long$row + 1L
  long$x <- long$col

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["pi"]])
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_distiller(
      name = expression(pi[ij]),
      palette = "Blues",
      direction = 1,
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    )

  # Grade-block boundaries between consecutive selected-grade groups.
  boundaries <- which(diff(ord_sel) != 0L)
  if (length(boundaries) > 0L) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = boundaries + 0.5, colour = "black", linewidth = 0.5
      ) +
      ggplot2::geom_hline(
        yintercept = n - boundaries + 0.5, colour = "black", linewidth = 0.5
      )
  }

  # Different-grade threshold outline.
  if (isTRUE(threshold)) {
    thr <- 1 / (1 + parts$lambda)
    hot <- long[is.finite(long$pi) & long$pi >= thr, , drop = FALSE]
    if (nrow(hot) > 0L) {
      p <- p + ggplot2::geom_tile(
        data = hot,
        ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
        fill = NA, colour = "red", linewidth = 0.3, inherit.aes = FALSE
      )
    }
  }

  # Per-unit tick labels only for small matrices; numeric index breaks (ggplot2
  # default) with blanked text otherwise.
  use_labels <- n <= 25L
  if (use_labels) {
    x_scale <- ggplot2::scale_x_continuous(
      breaks = seq_len(n), labels = ord_ids, expand = ggplot2::expansion(0, 0)
    )
    y_scale <- ggplot2::scale_y_continuous(
      breaks = seq_len(n), labels = rev(ord_ids), expand = ggplot2::expansion(0, 0)
    )
    axis_text_x <- ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
    axis_text_y <- ggplot2::element_text(size = 7)
  } else {
    x_scale <- ggplot2::scale_x_continuous(expand = ggplot2::expansion(0, 0))
    y_scale <- ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, 0))
    axis_text_x <- ggplot2::element_blank()
    axis_text_y <- ggplot2::element_blank()
  }

  p +
    x_scale +
    y_scale +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Unit j (grade-sorted)",
      y = "Unit i (grade-sorted)",
      title = "Posterior-contrast matrix"
    ) +
    theme_gradepath() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(7, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      axis.text.x = axis_text_x,
      axis.text.y = axis_text_y
    )
}

# Resolve posterior-contrast inputs from a gp_fit or a bare gp_pairwise. Returns
# a list with $matrix, $ids, $selected, $endpoint, $lambda.
.gp_contrast_inputs <- function(x,
                                selected_grades = NULL,
                                endpoint_grades = NULL,
                                selected_lambda = NULL) {
  if (inherits(x, "gp_fit")) {
    mat <- x$pairwise$matrix
    ids <- x$pairwise$ids
    if (is.null(ids)) {
      ids <- rownames(mat)
    }
    sel_fit <- x$selected_grade
    selected <- as.integer(sel_fit$assignment$grade[match(ids, sel_fit$assignment$id)])
    endpoint <- .gp_endpoint_grades(x$grade_path, ids)
    if (is.null(endpoint)) {
      endpoint <- selected
    }
    return(list(
      matrix = mat, ids = ids,
      selected = selected, endpoint = endpoint,
      lambda = sel_fit$lambda
    ))
  }

  if (inherits(x, "gp_pairwise")) {
    mat <- x$matrix
    ids <- x$ids
    if (is.null(ids)) {
      ids <- rownames(mat)
    }
    if (is.null(selected_grades)) {
      .gradepath_abort(
        "`gp_plot_posterior_contrast()` needs `selected_grades` when given a bare `gp_pairwise`."
      )
    }
    selected <- as.integer(selected_grades)
    if (length(selected) != length(ids)) {
      .gradepath_abort("`selected_grades` must have one grade per unit.")
    }
    endpoint <- if (is.null(endpoint_grades)) selected else as.integer(endpoint_grades)
    if (length(endpoint) != length(ids)) {
      .gradepath_abort("`endpoint_grades` must have one grade per unit.")
    }
    if (is.null(selected_lambda)) {
      .gradepath_abort(
        "`gp_plot_posterior_contrast()` needs `selected_lambda` when given a bare `gp_pairwise`."
      )
    }
    return(list(
      matrix = mat, ids = ids,
      selected = selected, endpoint = endpoint,
      lambda = as.numeric(selected_lambda)
    ))
  }

  .gradepath_abort(
    "`gp_plot_posterior_contrast()` accepts a `gp_fit` or a `gp_pairwise` with explicit grade vectors; got <%s>.",
    paste(class(x), collapse = "/")
  )
}

# Endpoint grade vector (in `ids` order) from a grade path. Prefers the fit at
# selection$endpoint_lambda; falls back to the lambda == 1 fit, then the last
# fit. Returns NULL when no grade path is available.
.gp_endpoint_grades <- function(grade_path, ids) {
  if (!inherits(grade_path, "gp_grade_path")) {
    return(NULL)
  }
  fits <- grade_path$fits
  lambdas <- vapply(fits, function(f) f$lambda, numeric(1))
  target <- grade_path$selection$endpoint_lambda
  idx <- NULL
  if (!is.null(target)) {
    hit <- which(abs(lambdas - target) < 1e-8)
    if (length(hit) > 0L) idx <- hit[[1L]]
  }
  if (is.null(idx)) {
    hit <- which(abs(lambdas - 1) < 1e-8)
    idx <- if (length(hit) > 0L) hit[[1L]] else length(fits)
  }
  fit <- fits[[idx]]
  as.integer(fit$assignment$grade[match(ids, fit$assignment$id)])
}

# ---------------------------------------------------------------------------
# autoplot() methods (thin wrappers; dispatch to gp_plot_frontier())
# ---------------------------------------------------------------------------

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Autoplot methods for gradepath frontier objects
#'
#' ggplot2 [ggplot2::autoplot()] methods so `autoplot(fr)` (a `gp_frontier`),
#' `autoplot(fit)` (a `gp_fit`), or `autoplot(path)` (a `gp_grade_path`) draw the
#' information-reliability frontier via [gp_plot_frontier()] without attaching
#' ggplot2.
#'
#' @param object A `gp_frontier`, `gp_fit`, or `gp_grade_path`.
#' @param ... Further arguments forwarded to [gp_plot_frontier()].
#'
#' @return A `ggplot` object (class `c("gg", "ggplot")`); it renders only when
#'   printed.
#'
#' @examples
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' autoplot(fit)   # the information-reliability frontier for a gp_fit
#'
#' @seealso [gp_plot_frontier()], [gp_frontier()]
#' @family gradepath-plots
#' @export
autoplot.gp_frontier <- function(object, ...) {
  gp_plot_frontier(object, ...)
}

#' @rdname autoplot.gp_frontier
#' @export
autoplot.gp_fit <- function(object, ...) {
  gp_plot_frontier(object, ...)
}

#' @rdname autoplot.gp_frontier
#' @export
autoplot.gp_grade_path <- function(object, ...) {
  gp_plot_frontier(object, ...)
}
