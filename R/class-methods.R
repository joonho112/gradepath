# ===========================================================================
# class-methods.R -- result-first console surface for the gradepath objects
# (GP-DEC-14-A).
#
# DESIGN. Every typed object gets three faces:
#   * format_gp_<type>_cli(x, ..., width)  EXPORTED; returns a character vector
#     (one element per line) drawing a PLAIN-ASCII box (`+ - |`) that LEADS WITH
#     THE RESULT and nothing else. It cat()s nothing and recomputes nothing -- it
#     reads only materialized fields. This is the canonical tier-2/3 surface
#     (SPEC-14.2/14.3): byte-stable, ASCII-only, no ANSI colour, no Unicode box
#     glyphs (an `ansi` tier may be layered later).
#   * print.gp_<type>(x, ..., width)  delegates to the formatter and cat()s the
#     lines, then returns invisible(x).
#   * summary.gp_<type>(object, ...)  returns a small typed gp_<type>_summary
#     (class c("gp_<type>_summary", "gp_summary")) carrying the headline scalars
#     plus a $provenance slot (backend / channel / selected-lambda rule). A
#     shared print.gp_summary renders it. PROVENANCE IS SURFACED ONLY HERE --
#     never in print() (GP-DEC-14-A).
#
# GP-DEC-14-A (hard): the strings "claim", "claim-bearing", "claim tier",
# "research_only", "backend_claim" and backend-disclaimer prose appear in NO
# print/summary/format output. v1 lectured "claim_bearing: FALSE" before the
# answer; v2 leads with the result. (See test-print-summary.R for the grep gate,
# which word-anchors "tier" so the legitimate object name "frontier" is allowed.)
#
# INVARIANT 2 (frozen, pinned here): print.gp_report_card /
# format_gp_report_card_cli carry the firms<->names welfare-flip footnote in
# HOUSE STYLE, citing Appendix A. The welfare label is NOT hard-coded to one
# application (it FLIPS): a row is described as "more discriminatory" or "higher
# posterior contact rate", never "beats"/"wins".
#
# Uses the package's internal `%gp_or%` (class-objects-output.R) for NULL
# coalescing; defines no operators of its own.
# ===========================================================================

# ---------------------------------------------------------------------------
# Shared low-level helpers (ASCII box drawing + number formatting). Internal.
# ---------------------------------------------------------------------------

# Resolve the rendering width. NULL -> the console width option -> a sane 80.
.gp_cli_width <- function(width = NULL) {
  if (is.null(width)) {
    width <- getOption("width", 80L)
  }
  width <- suppressWarnings(as.integer(width)[[1L]])
  if (is.na(width) || width < 24L) {
    width <- 24L
  }
  width
}

# Truncate a single line to `width` columns with an ASCII ellipsis, NA-safe.
# `width` is a single cap; `keep` is therefore scalar and recycles across the
# over-long elements (indexing it by the logical `too_long` would inject NA into
# every element past the first -- the classic length-1[logical] footgun).
.gp_cli_truncate <- function(line, width) {
  line <- as.character(line)
  line[is.na(line)] <- ""
  too_long <- nchar(line) > width
  if (any(too_long)) {
    keep <- max(as.integer(width)[[1L]] - 3L, 1L)
    line[too_long] <- paste0(substr(line[too_long], 1L, keep), "...")
  }
  line
}

# Draw a plain-ASCII box around a title line + body lines. The box auto-sizes to
# the widest content but never exceeds `width`. Returns a character vector
# (no trailing cat). Uses ONLY `+`, `-`, `|`, and space -- ASCII tier.
.gp_cli_box <- function(title, body = character(), width = NULL) {
  width <- .gp_cli_width(width)
  inner_cap <- width - 4L # 2 borders + 2 pad spaces

  title <- .gp_cli_truncate(title, inner_cap)
  body <- if (length(body)) .gp_cli_truncate(body, inner_cap) else character()

  content_width <- max(nchar(c(title, body)), 0L)
  inner <- min(max(content_width, 1L), inner_cap)

  bar <- paste0("+", strrep("-", inner + 2L), "+")
  pad_line <- function(s) {
    s <- as.character(s)
    s[is.na(s)] <- ""
    paste0("| ", formatC(s, width = -inner, flag = "-"), " |")
  }

  out <- c(bar, pad_line(title), bar)
  if (length(body)) {
    out <- c(out, pad_line(body), bar)
  }
  out
}

# Footnote / pointer lines that sit BELOW the box (ASCII `i ` info prefix).
# These are NEVER width-truncated: they carry load-bearing references (e.g. the
# frozen INVARIANT 2 welfare note, the `?gp_report_card` / Appendix A pointers).
# Truncating them could drop the invariant, so they always render in full. The
# `width` argument is accepted for a uniform call signature but ignored here.
.gp_cli_note <- function(text, width = NULL) {
  paste0("i  ", text)
}

# Fixed-precision number formatter. NA / non-finite -> "NA" (never "NaN"/"Inf"
# leakage in a clean box). Scalars and vectors both supported.
.gp_fmt_num <- function(x, digits = 3L) {
  out <- rep("NA", length(x))
  ok <- is.finite(x)
  if (any(ok)) {
    out[ok] <- formatC(x[ok], format = "f", digits = digits)
  }
  out
}

# Format a closed interval [lo, hi] with NA-safe endpoints.
.gp_fmt_range <- function(lo, hi, digits = 3L) {
  sprintf("[%s, %s]", .gp_fmt_num(lo, digits), .gp_fmt_num(hi, digits))
}

# Compose a grade-composition string like "3 (2/81/14)" from an integer grade
# vector: K distinct grades, then the per-grade counts in grade order 1..K.
.gp_grade_composition <- function(grades) {
  grades <- grades[!is.na(grades)]
  if (!length(grades)) {
    return("0 ()")
  }
  k <- length(unique(grades))
  counts <- tabulate(match(grades, sort(unique(grades))))
  sprintf("%d (%s)", k, paste(counts, collapse = "/"))
}

# Read the enriched frontier-metric columns from a grade-path summary IF present.
# Returns a one-row list (reliability/tau_bar at the selected lambda) or NULL.
# Pure read; never recomputes, never calls the solver.
.gp_path_selected_metrics <- function(path) {
  summary <- path$summary
  sel <- path$selection$selected_lambda
  needed <- c("lambda", "reliability", "discordance_rate", "tau_bar")
  if (!is.data.frame(summary) || !all(needed %in% names(summary)) ||
    is.null(sel) || !is.finite(sel)) {
    return(NULL)
  }
  idx <- which(abs(summary$lambda - sel) < 1e-8)
  if (length(idx) != 1L) {
    return(NULL)
  }
  list(
    reliability = summary$reliability[[idx]], # (1 - DR)
    discordance_rate = summary$discordance_rate[[idx]],
    tau_bar = summary$tau_bar[[idx]]
  )
}

# Render a scalar slot value as a short text token (NA-safe, no recompute).
# Whole numbers (counts) render without decimals; fractionals to 3 digits.
.gp_scalar_to_text <- function(v) {
  if (is.null(v) || length(v) == 0L) {
    return("NA")
  }
  if (is.numeric(v)) {
    out <- .gp_fmt_num(v, digits = 3L)
    whole <- is.finite(v) & (v == round(v))
    out[whole] <- formatC(v[whole], format = "d")
    return(paste(out, collapse = ", "))
  }
  v <- as.character(v)
  v[is.na(v)] <- "NA"
  paste(v, collapse = ", ")
}

# Typed summary constructor (shared). Builds a c("gp_<type>_summary",
# "gp_summary") object from named scalars; a `provenance` argument (named list)
# is stored last and rendered as its own block by print.gp_summary.
.gp_new_summary <- function(type, ..., provenance = list(),
                            interpretation = NULL, glossary = NULL) {
  fields <- list(...)
  fields$provenance <- provenance
  # `interpretation` / `glossary` are stored as ATTRIBUTES (not list fields), so
  # the scalar loop in print.gp_summary (names(x) minus "provenance") is unchanged
  # and a summary built without them is byte-identical to before. structure()
  # drops a NULL attribute, so the no-extras object carries no new attributes.
  structure(
    fields,
    class = c(paste0(type, "_summary"), "gp_summary", "list"),
    gp_type = type,
    gp_interpretation = interpretation,
    gp_glossary = glossary
  )
}

# Coerce a value to a clean character vector for rendering: NULL/length-0 ->
# empty; non-character coerced; NA -> ""; names preserved. Never errors -- the
# printer must be unbreakable on a hand-built or partially-populated summary.
.gp_summary_norm_chr <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(stats::setNames(character(0), character(0)))
  }
  nm <- names(x)
  x <- as.character(x)
  x[is.na(x)] <- ""
  if (!is.null(nm)) {
    nm[is.na(nm)] <- ""
    names(x) <- nm
  }
  x
}

# Render the optional Interpretation: and Glossary: blocks (attribute-stored
# extras) for print.gp_summary. Headers at 2-space indent (parallel to the scalar
# labels and `provenance:`); body at 4-space indent. Interpretation lines are
# width-wrapped with a hanging indent; glossary lines are aligned `term : def`
# and width-truncated. Returns character(0) when neither block has content, so a
# no-extras summary prints byte-identically to before.
.gp_render_summary_extras <- function(x, width = NULL) {
  width <- .gp_cli_width(width)
  head_indent <- "  "
  body_indent <- "    "
  out <- character(0)

  interp <- .gp_summary_norm_chr(attr(x, "gp_interpretation", exact = TRUE))
  interp <- interp[nzchar(interp)]
  if (length(interp)) {
    wrapped <- unlist(lapply(interp, function(line) {
      strwrap(line, width = width, prefix = body_indent, initial = body_indent)
    }), use.names = FALSE)
    out <- c(out, paste0(head_indent, "Interpretation:"), wrapped)
  }

  gloss <- .gp_summary_norm_chr(attr(x, "gp_glossary", exact = TRUE))
  terms <- names(gloss)
  if (is.null(terms)) terms <- rep("", length(gloss))
  keep <- nzchar(terms)
  gloss <- gloss[keep]
  terms <- terms[keep]
  if (length(gloss)) {
    pad <- max(nchar(terms), 0L)
    lines <- sprintf(
      "%s%s : %s",
      body_indent,
      formatC(terms, width = -pad, flag = "-"),
      gloss
    )
    lines <- .gp_cli_truncate(lines, width)
    out <- c(out, paste0(head_indent, "Glossary:"), lines)
  }

  out
}

# ---------------------------------------------------------------------------
# Shared summary printer. Every summary(x) returns an object of class
# c("gp_<type>_summary", "gp_summary"); this renders it. Provenance lives here.
# ---------------------------------------------------------------------------

#' Print a `gp_summary`
#'
#' Shared printer for the typed summary objects returned by `summary()` on a
#' gradepath result. Renders the headline scalars first, then a `provenance:`
#' block (backend, channel, selected-lambda rule) -- provenance appears ONLY in
#' the summary, never in `print()` of the object itself.
#'
#' @param x A `gp_summary` (e.g. from `summary(gp_frontier_object)`).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.gp_summary <- function(x, ...) {
  type <- attr(x, "gp_type") %gp_or% sub("_summary$", "", class(x)[[1L]])
  cat(sprintf("<summary: %s>\n", type))

  scalars <- x[setdiff(names(x), "provenance")]
  if (length(scalars)) {
    labels <- names(scalars)
    pad <- max(nchar(labels), 0L)
    for (i in seq_along(scalars)) {
      cat(sprintf(
        "  %s : %s\n",
        formatC(labels[[i]], width = -pad, flag = "-"),
        .gp_scalar_to_text(scalars[[i]])
      ))
    }
  }

  # Optional Interpretation: / Glossary: blocks (attribute-stored extras), AFTER
  # the scalars and BEFORE provenance. character(0) when absent, so a no-extras
  # summary prints byte-identically to before.
  extras <- .gp_render_summary_extras(x, width = getOption("width"))
  if (length(extras)) {
    cat(extras, sep = "\n")
    cat("\n")
  }

  prov <- x$provenance
  if (length(prov)) {
    cat("  provenance:\n")
    labels <- names(prov)
    pad <- max(nchar(labels), 0L)
    for (i in seq_along(prov)) {
      cat(sprintf(
        "    %s : %s\n",
        formatC(labels[[i]], width = -pad, flag = "-"),
        .gp_scalar_to_text(prov[[i]])
      ))
    }
  }
  invisible(x)
}

# ===========================================================================
# gp_pairwise
# ===========================================================================

#' Format a `gp_pairwise` for the console
#'
#' Returns the plain-ASCII lines that [print.gp_pairwise()] emits. Leads with the
#' result: matrix dimensions (`J x J`), the off-diagonal pi range, the number of
#' ordered pairs, and the power. Recomputes nothing.
#'
#' @param x A `gp_pairwise`.
#' @param ... Unused.
#' @param width Optional integer rendering width (defaults to `getOption("width")`).
#' @return A `character` vector, one element per line.
#'
#' @examples
#' # The bundled tiny fit; the formatter only reads materialized fields (instant).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' pw <- get_pairwise(fit)
#' writeLines(format_gp_pairwise_cli(pw))
#'
#' @seealso [print.gp_pairwise()], [gp_pairwise()], [get_pairwise()]
#' @family gradepath-cli
#' @export
format_gp_pairwise_cli <- function(x, ..., width = NULL) {
  J <- length(x$ids)
  m <- x$matrix
  off <- row(m) != col(m)
  vals <- m[off]
  rng <- if (length(vals)) range(vals) else c(NA_real_, NA_real_)
  n_pairs <- J * (J - 1L)

  title <- sprintf("gp_pairwise  .  %d x %d  .  %d ordered pairs", J, J, n_pairs)
  body <- c(
    sprintf(
      "pi range (off-diag): %s   diag = %s",
      .gp_fmt_range(rng[[1L]], rng[[2L]]),
      .gp_fmt_num(0.5, 2L)
    ),
    sprintf(
      "power = %d   rule = %s",
      as.integer(x$power %gp_or% 0L),
      as.character(x$source$rule %gp_or% "NA")
    )
  )
  c(
    .gp_cli_box(title, body, width = width),
    .gp_cli_note("matrix: as.matrix(pw)   ids: pw$ids", width = width)
  )
}

#' @rdname format_gp_pairwise_cli
#' @param object A `gp_pairwise`.
#' @return `summary()` returns a typed `gp_pairwise_summary`.
#' @export
summary.gp_pairwise <- function(object, ...) {
  J <- length(object$ids)
  m <- object$matrix
  off <- row(m) != col(m)
  vals <- m[off]
  rng <- if (length(vals)) range(vals) else c(NA_real_, NA_real_)
  pi_min <- rng[[1L]]
  pi_max <- rng[[2L]]
  ordered_pairs <- J * (J - 1L)
  power <- as.integer(object$power %gp_or% 0L)

  interp <- sprintf(
    "The %d x %d pairwise outranking matrix Pi over %d units (%d ordered pairs).",
    J, J, J, ordered_pairs)
  if (is.finite(pi_min) && is.finite(pi_max)) {
    interp <- c(interp, sprintf(paste0(
      "Off-diagonal pi ranges from %s to %s; pi_ij is the posterior probability ",
      "that unit i is more extreme than unit j (the diagonal is fixed at 0.5)."),
      .gp_fmt_num(pi_min, 2L), .gp_fmt_num(pi_max, 2L)))
  }

  .gp_new_summary(
    "gp_pairwise",
    units = J,
    ordered_pairs = ordered_pairs,
    pi_min = pi_min,
    pi_max = pi_max,
    power = power,
    interpretation = interp,
    glossary = c(
      pi_ij =
        "posterior probability that unit i is more extreme than unit j; the diagonal is 0.5 by convention.",
      ordered_pairs =
        "units x (units - 1): the count of off-diagonal cells (each ordered pair i != j).",
      power =
        "the outer-product power applied to Pi; 0 = the raw posterior outranking matrix."
    ),
    provenance = list(
      stage = object$source$stage %gp_or% NA_character_,
      rule = object$source$rule %gp_or% NA_character_,
      assumption = object$source$assumption %gp_or% NA_character_,
      backend = object$control$backend %gp_or% NA_character_
    )
  )
}

#' @rdname format_gp_pairwise_cli
#' @return `print()` returns `x` invisibly.
#' @export
print.gp_pairwise <- function(x, ..., width = getOption("width")) {
  cat(format_gp_pairwise_cli(x, width = width), sep = "\n")
  invisible(x)
}

# ===========================================================================
# gp_grade_fit
# ===========================================================================

# Build the grade -> size body lines (HEAD with overflow), from an assignment df.
.gp_grade_size_lines <- function(assignment, max_rows = 8L) {
  grades <- sort(unique(assignment$grade))
  sizes <- tabulate(match(assignment$grade, grades))
  lines <- sprintf("grade %d: %d", grades, sizes)
  if (length(lines) > max_rows) {
    head_lines <- lines[seq_len(max_rows)]
    lines <- c(
      head_lines,
      sprintf("... %d more grades ...", length(lines) - max_rows)
    )
  }
  lines
}

#' Format a `gp_grade_fit` for the console
#'
#' Returns the plain-ASCII lines that [print.gp_grade_fit()] emits. Leads with
#' the result: grade count, the grade->size table, solver status, and lambda.
#' Backend/channel detail is surfaced via `summary()`, not here. Recomputes
#' nothing.
#'
#' @param x A `gp_grade_fit`.
#' @param ... Unused.
#' @param width Optional integer rendering width.
#' @return A `character` vector, one element per line.
#'
#' @examples
#' # The bundled tiny fit; its selected slice is a gp_grade_fit (instant read).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' writeLines(format_gp_grade_fit_cli(fit$selected_grade))
#'
#' @seealso [print.gp_grade_fit()], [gp_select_grade()], [gp_grade_path()]
#' @family gradepath-cli
#' @export
format_gp_grade_fit_cli <- function(x, ..., width = NULL) {
  grade_count <- as.integer(x$summary$grade_count %gp_or% length(unique(x$assignment$grade)))
  status <- as.character(x$summary$status %gp_or% x$backend$status %gp_or% "NA")
  n_units <- as.integer(x$summary$n_units %gp_or% length(x$ids))

  title <- sprintf(
    "gp_grade_fit  .  %d grades  .  lambda = %s",
    grade_count, .gp_fmt_num(x$lambda %gp_or% NA_real_, 2L)
  )
  body <- c(
    sprintf("units = %d   status = %s", n_units, status),
    .gp_grade_size_lines(x$assignment)
  )
  c(
    .gp_cli_box(title, body, width = width),
    .gp_cli_note("assignment: x$assignment   summary(x) for backend/channel details", width = width)
  )
}

#' @rdname format_gp_grade_fit_cli
#' @param object A `gp_grade_fit`.
#' @return `summary()` returns a typed `gp_grade_fit_summary`.
#' @export
summary.gp_grade_fit <- function(object, ...) {
  grade_count <- as.integer(object$summary$grade_count %gp_or% length(unique(object$assignment$grade)))
  lambda <- object$lambda
  n_units <- as.integer(object$summary$n_units %gp_or% length(object$ids))
  obj <- object$objective$value %gp_or% NA_real_
  status <- as.character(object$summary$status %gp_or% object$backend$status %gp_or% NA_character_)

  interp <- sprintf(
    "A single solved grading slice: %d units graded into %d grades at lambda = %s.",
    n_units, grade_count, .gp_fmt_num(lambda %gp_or% NA_real_, 2L))
  have_status <- length(status) == 1L && !is.na(status) && nzchar(trimws(status))
  have_obj <- is.finite(obj)
  if (have_status && have_obj) {
    interp <- c(interp, sprintf(
      "Solver status: %s; objective value = %s at the solved optimum.",
      status, .gp_fmt_num(obj, 3L)))
  } else if (have_status) {
    interp <- c(interp, sprintf("Solver status: %s.", status))
  } else if (have_obj) {
    interp <- c(interp, sprintf(
      "Objective value = %s at the solved optimum.", .gp_fmt_num(obj, 3L)))
  }

  .gp_new_summary(
    "gp_grade_fit",
    grade_count = grade_count,
    lambda = lambda,
    n_units = n_units,
    objective = obj,
    status = status,
    interpretation = interp,
    glossary = c(
      lambda =
        "the single grade-count selection knob for this slice (KRW baseline 0.25); larger lambda yields more, finer grades (smaller lambda pools units into fewer, coarser grades).",
      objective =
        "the grade integer-program objective value at the solved optimum.",
      status =
        "the solver's termination status (e.g. optimal / time_limit / gap)."
    ),
    provenance = list(
      backend = object$backend$name %gp_or% object$control$backend %gp_or% NA_character_,
      channel = object$backend$path %gp_or% NA_character_,
      mipgap = object$backend$mipgap %gp_or% NA_real_,
      runtime = object$backend$runtime %gp_or% NA_real_
    )
  )
}

#' @rdname format_gp_grade_fit_cli
#' @return `print()` returns `x` invisibly.
#' @export
print.gp_grade_fit <- function(x, ..., width = getOption("width")) {
  cat(format_gp_grade_fit_cli(x, width = width), sep = "\n")
  invisible(x)
}

# ===========================================================================
# gp_grade_path
# ===========================================================================

#' Format a `gp_grade_path` for the console
#'
#' Returns the plain-ASCII lines that [print.gp_grade_path()] emits. Leads with
#' the result: number of stored lambda values, the selected lambda, the
#' grade-count range across the path, and the endpoint lambda. Recomputes
#' nothing.
#'
#' @param x A `gp_grade_path`.
#' @param ... Unused.
#' @param width Optional integer rendering width.
#' @return A `character` vector, one element per line.
#'
#' @examples
#' # The bundled tiny fit carries a solved grade path (instant read; no solve).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' writeLines(format_gp_grade_path_cli(fit$grade_path))
#'
#' @seealso [print.gp_grade_path()], [gp_grade_path()], [gp_select_grade()]
#' @family gradepath-cli
#' @export
format_gp_grade_path_cli <- function(x, ..., width = NULL) {
  n_lambda <- length(x$lambda_grid)
  sel <- x$selection$selected_lambda %gp_or% NA_real_
  endpoint <- x$selection$endpoint_lambda %gp_or% NA_real_
  gc <- x$summary$grade_count
  gc_rng <- if (length(gc)) range(gc) else c(NA_real_, NA_real_)

  title <- sprintf(
    "gp_grade_path  .  %d lambda values  .  selected lambda = %s",
    n_lambda, .gp_fmt_num(sel, 2L)
  )
  body <- c(
    sprintf(
      "grades across path in [%d, %d]",
      as.integer(gc_rng[[1L]]), as.integer(gc_rng[[2L]])
    ),
    sprintf("endpoint lambda = %s", .gp_fmt_num(endpoint, 2L))
  )
  c(
    .gp_cli_box(title, body, width = width),
    .gp_cli_note("path table: x$summary   per-lambda fits: x$fits", width = width)
  )
}

#' @rdname format_gp_grade_path_cli
#' @param object A `gp_grade_path`.
#' @return `summary()` returns a typed `gp_grade_path_summary`.
#' @export
summary.gp_grade_path <- function(object, ...) {
  gc <- object$summary$grade_count
  gc_rng <- if (length(gc)) range(gc) else c(NA_real_, NA_real_)

  n_lambda <- length(object$lambda_grid)
  sel <- object$selection$selected_lambda %gp_or% NA_real_
  endpoint <- object$selection$endpoint_lambda %gp_or% NA_real_
  gc_min <- as.integer(gc_rng[[1L]])
  gc_max <- as.integer(gc_rng[[2L]])

  head <- sprintf("A path of optimal gradings over %d lambda values.", n_lambda)
  if (is.finite(sel)) {
    sel_clause <- if (is.finite(endpoint)) {
      sprintf("The selected lambda = %s (endpoint lambda = %s).",
              .gp_fmt_num(sel, 2L), .gp_fmt_num(endpoint, 2L))
    } else {
      sprintf("The selected lambda = %s.", .gp_fmt_num(sel, 2L))
    }
    head <- paste(head, sel_clause)
  } else if (is.finite(endpoint)) {
    head <- paste(head, sprintf("The endpoint lambda = %s.",
                                .gp_fmt_num(endpoint, 2L)))
  }
  interp <- head
  if (is.finite(gc_min) && is.finite(gc_max)) {
    interp <- c(interp, sprintf(paste0(
      "Across the sweep the number of grades ranges from %d to %d ",
      "(smaller lambda yields more, finer grades)."),
      gc_min, gc_max))
  }

  .gp_new_summary(
    "gp_grade_path",
    n_lambda = n_lambda,
    selected_lambda = sel,
    endpoint_lambda = endpoint,
    grade_count_min = gc_min,
    grade_count_max = gc_max,
    interpretation = interp,
    glossary = c(
      selected_lambda =
        "the grade-count selection knob (KRW baseline 0.25); larger lambda yields more, finer grades (smaller lambda pools units into fewer, coarser grades).",
      endpoint_lambda =
        "the lambda = 1 anchor used for the Condorcet endpoint ranking (the finest, fully ordered end of the path).",
      grade_count_min =
        "fewest grades seen along the path (the coarsest grading, at the smallest lambda).",
      grade_count_max =
        "most grades seen along the path (the finest grading, at the largest lambda)."
    ),
    provenance = list(
      backend = object$backend$name %gp_or% object$control$backend %gp_or% NA_character_,
      selection_rule = object$selection$selection_rule %gp_or% NA_character_
    )
  )
}

#' @rdname format_gp_grade_path_cli
#' @return `print()` returns `x` invisibly.
#' @export
print.gp_grade_path <- function(x, ..., width = getOption("width")) {
  cat(format_gp_grade_path_cli(x, width = width), sep = "\n")
  invisible(x)
}

# ===========================================================================
# gp_frontier  (Chapter 14 worked-plot target box)
# ===========================================================================

# Locate the selected frontier row (by lambda match). Pure read.
.gp_frontier_selected_row <- function(x) {
  sel <- x$selection$selected_lambda
  tab <- x$table
  if (is.null(sel) || !is.finite(sel) || !nrow(tab)) {
    return(NULL)
  }
  idx <- which(abs(tab$lambda - sel) < 1e-8)
  if (length(idx) != 1L) {
    return(NULL)
  }
  tab[idx, , drop = FALSE]
}

#' Format a `gp_frontier` for the console
#'
#' Returns the plain-ASCII lines that [print.gp_frontier()] emits, reproducing
#' the Chapter 14 worked-plot box. Leads with the result: the selected lambda and
#' its `(1 - DR, tau-bar, grade_count)`, the `(1 - DR)` and `tau-bar` ranges
#' across the path, and how many naive benchmarks are available. Recomputes
#' nothing.
#'
#' @param x A `gp_frontier`.
#' @param ... Unused.
#' @param width Optional integer rendering width.
#' @return A `character` vector, one element per line.
#'
#' @examples
#' # Build the frontier from the bundled fit's grade path + pairwise (instant).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' fr <- gp_frontier(fit$grade_path, pairwise = get_pairwise(fit))
#' writeLines(format_gp_frontier_cli(fr))
#'
#' @seealso [print.gp_frontier()], [gp_frontier()], [gp_plot_frontier()]
#' @family gradepath-cli
#' @export
format_gp_frontier_cli <- function(x, ..., width = NULL) {
  tab <- x$table
  n_lambda <- nrow(tab)
  sel <- x$selection$selected_lambda %gp_or% NA_real_
  row <- .gp_frontier_selected_row(x)

  dr_rng <- if (n_lambda) range(tab$reliability) else c(NA_real_, NA_real_)
  tau_rng <- if (n_lambda) range(tab$tau_bar) else c(NA_real_, NA_real_)
  n_bench <- nrow(x$benchmarks)

  title <- sprintf(
    "gp_frontier  .  %d lambda values  .  selected lambda = %s",
    n_lambda, .gp_fmt_num(sel, 2L)
  )

  sel_line <- if (is.null(row)) {
    "selected:  (not on the stored path)"
  } else {
    sprintf(
      "selected:  (1 - DR) = %s   tau-bar = %s   grades = %d",
      .gp_fmt_num(row$reliability, 2L),
      .gp_fmt_num(row$tau_bar, 2L),
      as.integer(row$grade_count)
    )
  }
  body <- c(
    sel_line,
    sprintf(
      "range:     1 - DR in %s   tau-bar in %s",
      .gp_fmt_range(dr_rng[[1L]], dr_rng[[2L]], 2L),
      .gp_fmt_range(tau_rng[[1L]], tau_rng[[2L]], 2L)
    ),
    sprintf("benchmarks: %d naive overlays available (see gp_plot_*)", n_bench)
  )
  c(
    .gp_cli_box(title, body, width = width),
    .gp_cli_note("gp_plot_frontier(fr) draws the (1 - DR, tau-bar) frontier.", width = width)
  )
}

#' @rdname format_gp_frontier_cli
#' @param object A `gp_frontier`.
#' @return `summary()` returns a typed `gp_frontier_summary`.
#' @export
summary.gp_frontier <- function(object, ...) {
  tab <- object$table
  row <- .gp_frontier_selected_row(object)
  dr_rng <- if (nrow(tab)) range(tab$reliability) else c(NA_real_, NA_real_)
  tau_rng <- if (nrow(tab)) range(tab$tau_bar) else c(NA_real_, NA_real_)
  n_lambda <- nrow(tab)
  sel <- object$selection$selected_lambda %gp_or% NA_real_
  sel_reliability <- if (is.null(row)) NA_real_ else row$reliability # (1 - DR)
  sel_tau_bar <- if (is.null(row)) NA_real_ else row$tau_bar
  sel_grades <- if (is.null(row)) NA_integer_ else as.integer(row$grade_count)
  n_bench <- nrow(object$benchmarks)

  interp <- sprintf(
    "The information-reliability frontier across %d lambda values.", n_lambda)
  if (!is.null(row) && is.finite(sel_reliability) && is.finite(sel_tau_bar) &&
      is.finite(sel_grades)) {
    interp <- c(interp, sprintf(paste0(
      "At the selected lambda = %s the grading reaches reliability ",
      "(1 - DR) = %s and tau-bar = %s with %d grades."),
      .gp_fmt_num(sel, 2L), .gp_fmt_num(sel_reliability, 3L),
      .gp_fmt_num(sel_tau_bar, 3L), sel_grades))
  }
  if (is.finite(n_bench) && n_bench > 0L) {
    interp <- c(interp, sprintf(paste0(
      "%d naive benchmark rules are shown for comparison ",
      "(the frontier path is at least as reliable at equal coarseness)."),
      n_bench))
  }

  .gp_new_summary(
    "gp_frontier",
    n_lambda = n_lambda,
    selected_lambda = sel,
    selected_reliability = sel_reliability, # (1 - DR)
    selected_tau_bar = sel_tau_bar,
    selected_grades = sel_grades,
    reliability_min = dr_rng[[1L]],
    reliability_max = dr_rng[[2L]],
    tau_bar_min = tau_rng[[1L]],
    tau_bar_max = tau_rng[[2L]],
    benchmarks = n_bench,
    interpretation = interp,
    glossary = c(
      reliability =
        "1 - discordance rate: share of comparable pairs the grades rank the same way as the posterior.",
      tau_bar =
        "posterior expected Kendall rank agreement over between-grade blocks, in [-1, 1]; higher = more rank information carried by the grading.",
      benchmarks =
        "count of naive grading rules overlaid on the frontier for comparison (see gp_plot_frontier)."
    ),
    provenance = list(
      backend = object$control$backend %gp_or% NA_character_,
      selected_lambda_rule = "table row matched on selection$selected_lambda"
    )
  )
}

#' @rdname format_gp_frontier_cli
#' @return `print()` returns `x` invisibly.
#' @export
print.gp_frontier <- function(x, ..., width = getOption("width")) {
  cat(format_gp_frontier_cli(x, width = width), sep = "\n")
  invisible(x)
}

# ===========================================================================
# gp_report_card  (Chapter 12 sec.6 target box; carries frozen INVARIANT 2)
# ===========================================================================

# House-style welfare footnote (frozen invariant 2). The label is NOT hard-coded
# to one application -- it FLIPS between firms and names. Cites Appendix A.
.gp_report_card_welfare_note <- function() {
  paste0(
    "grade 1 = most-extreme theta ",
    "(firms: most discriminatory; names: best-treated). ",
    "See ?gp_report_card and Appendix A."
  )
}

# Render the head of the report-card table as fixed-width ASCII rows. Matching
# the Chapter 12 sec.6 box, a long table shows the first (head_rows - 1) ranked
# rows, then a "... N more rows ..." line, then the FINAL row -- so the
# most-extreme and least-extreme units both appear.
.gp_report_card_table_lines <- function(tab, level, head_rows = 8L) {
  J <- nrow(tab)
  head_rows <- max(as.integer(head_rows), 2L)
  ci_label <- sprintf("CI (%d%%)", as.integer(round(100 * level)))

  fmt_rows <- function(idx) {
    data.frame(
      sort_rank = format(tab$sort_rank[idx]),
      id = format(as.character(tab$id[idx])),
      label = format(.gp_cli_truncate(as.character(tab$label[idx]), 16L)),
      grade = format(tab$grade[idx]),
      posterior_mean = .gp_fmt_num(tab$posterior_mean[idx], 3L),
      ci = .gp_fmt_range(tab$lower[idx], tab$upper[idx], 3L),
      estimate = .gp_fmt_num(tab$estimate[idx], 3L),
      stringsAsFactors = FALSE
    )
  }
  header <- data.frame(
    sort_rank = "sort_rank", id = "id", label = "label", grade = "grade",
    posterior_mean = "posterior_mean", ci = ci_label, estimate = "estimate",
    stringsAsFactors = FALSE
  )

  truncated <- J > head_rows
  shown_idx <- if (truncated) c(seq_len(head_rows - 1L), J) else seq_len(J)
  block <- rbind(header, fmt_rows(shown_idx))

  # Per-column width = max over header+body; right-justify numerics, left labels.
  widths <- vapply(block, function(col) max(nchar(col)), integer(1))
  just <- c(
    sort_rank = "right", id = "left", label = "left", grade = "left",
    posterior_mean = "right", ci = "left", estimate = "right"
  )
  pad_col <- function(col, w, side) {
    formatC(col, width = if (identical(side, "right")) w else -w, flag = "-")
  }
  render <- function(r) {
    cells <- vapply(names(block), function(nm) {
      pad_col(block[[nm]][r], widths[[nm]], just[[nm]])
    }, character(1))
    paste0("  ", paste(cells, collapse = "  "))
  }
  lines <- vapply(seq_len(nrow(block)), render, character(1))

  if (truncated) {
    # lines: [1]=header, [2..head_rows]=first (head_rows-1) rows, [last]=final row.
    n_head <- head_rows - 1L
    ell <- sprintf("  ... %d more rows ...", J - head_rows)
    lines <- c(lines[seq_len(1L + n_head)], ell, lines[[length(lines)]])
  }
  lines
}

#' Format a `gp_report_card` for the console
#'
#' Returns the plain-ASCII lines that [print.gp_report_card()] emits, reproducing
#' the Chapter 12 target box. Leads with the result: unit count, grade count and
#' composition (e.g. `3 (2/81/14)`), selected lambda, a HEAD of the ranked table
#' (and the final row when truncated), then the frozen INVARIANT 2 welfare
#' footnote and pointers. Recomputes nothing.
#'
#' The welfare note (invariant 2) describes grade 1 as the most-extreme theta and
#' is NOT hard-coded to one application: it flips between firms (most
#' discriminatory) and names (best-treated). Rows are described as "more
#' discriminatory" / "higher posterior contact rate", not as a contest.
#'
#' @param x A `gp_report_card`.
#' @param ... Unused.
#' @param width Optional integer rendering width.
#' @param max_rows Integer; how many ranked rows to show before an "... N more
#'   rows ..." line (default 8). Controls the HEAD length for large tables.
#' @return A `character` vector, one element per line.
#'
#' @examples
#' # The bundled fit's report card; the formatter only renders materialized rows.
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' card <- get_report_card(fit)
#' writeLines(format_gp_report_card_cli(card))
#'
#' @seealso [print.gp_report_card()], [gp_report_card()], [get_report_card()],
#'   [gp_plot_report_card()]
#' @family gradepath-cli
#' @export
format_gp_report_card_cli <- function(x, ..., width = NULL, max_rows = 8L) {
  width <- .gp_cli_width(width)
  tab <- x$table
  J <- nrow(tab)
  unit <- as.character(x$control$unit %gp_or% "unit")
  level <- x$control$interval_level %gp_or% 0.90
  composition <- .gp_grade_composition(x$grades)
  sel <- x$selected_lambda %gp_or% NA_real_

  header <- sprintf(
    "<gp_report_card>  unit: %s | %d rows | grades: %s | selected lambda: %s",
    unit, J, composition, .gp_fmt_num(sel, 2L)
  )
  subhead <- "sorted by Condorcet rank (endpoint lambda = 1; secondary key: id)"

  c(
    # Header + subhead are the RESULT line; never truncated (leading with the
    # result means the selected lambda etc. must always be visible).
    header,
    subhead,
    "",
    .gp_report_card_table_lines(tab, level, head_rows = max_rows),
    "",
    # ---- INVARIANT 2 (frozen): the firms<->names welfare-flip footnote --------
    .gp_cli_note(.gp_report_card_welfare_note(), width = width),
    .gp_cli_note(
      "full table: as.data.frame(card)   formatted CLI: format_gp_report_card_cli(card)",
      width = width
    )
  )
}

#' @rdname format_gp_report_card_cli
#' @param object A `gp_report_card`.
#' @return `summary()` returns a typed `gp_report_card_summary`.
#' @export
summary.gp_report_card <- function(object, ...) {
  units <- nrow(object$table)
  grade_count <- length(unique(object$grades))
  composition <- .gp_grade_composition(object$grades)
  dist <- sub("^[0-9]+ ", "", composition)   # the bare "(n1/n2/...)" distribution
  sel <- object$selected_lambda %gp_or% NA_real_
  level <- object$control$interval_level %gp_or% NA_real_

  # Interpretation: lead with the result + the FROZEN INVARIANT 2 welfare note
  # (reused verbatim via .gp_report_card_welfare_note(), never paraphrased).
  interp <- c(
    sprintf("%d units ranked into %d grades %s; grade 1 is the most extreme.",
            units, grade_count, dist),
    .gp_report_card_welfare_note()
  )
  if (is.finite(level)) {
    interp <- c(interp, sprintf(
      "Each row carries a %d%% credible interval for the unit's posterior theta.",
      round(100 * level)))
  }

  .gp_new_summary(
    "gp_report_card",
    units = units,
    grade_count = grade_count,
    composition = composition,
    selected_lambda = sel,
    interval_level = level,
    interpretation = interp,
    glossary = c(
      composition =
        "K (n1/n2/.../nK): number of grades, then the unit count in each grade from grade 1 up.",
      selected_lambda =
        "the grade-count selection knob (KRW baseline 0.25); larger lambda yields more, finer grades (smaller lambda pools units into fewer, coarser grades).",
      interval_level =
        "credible-interval level for the per-unit posterior intervals (e.g. 0.9 = 90% interval)."
    ),
    provenance = list(
      backend = object$control$backend %gp_or% NA_character_,
      sort_rule = "Condorcet rank at endpoint lambda; secondary key id"
    )
  )
}

#' @rdname format_gp_report_card_cli
#' @return `print()` returns `x` invisibly.
#' @export
print.gp_report_card <- function(x, ..., width = getOption("width"), max_rows = 8L) {
  cat(format_gp_report_card_cli(x, width = width, max_rows = max_rows), sep = "\n")
  invisible(x)
}

#' @rdname format_gp_report_card_cli
#' @param row.names,optional Passed through from [base::as.data.frame()];
#'   `row.names` resets the payload's row names.
#' @return `as.data.frame()` returns the 10-column payload (`x$table`), the
#'   object a researcher writes to CSV; this is the working target the print
#'   pointer advertises.
#' @export
as.data.frame.gp_report_card <- function(x, row.names = NULL, optional = FALSE, ...) {
  tbl <- x$table
  rownames(tbl) <- row.names
  tbl
}

# ===========================================================================
# gp_fit  (composite: one line per stage)
# ===========================================================================

#' Format a `gp_fit` for the console
#'
#' Returns the plain-ASCII lines that [print.gp_fit()] emits: a composite,
#' one-line-per-stage view leading with the result -- unit count, the selected
#' grade count and composition, the information-reliability metrics
#' `(1 - DR, tau-bar)` (read from the enriched grade-path summary when present),
#' and the backend. Then pointers to `summary()` / `get_report_card()`. Recomputes
#' nothing and never calls the solver: the metrics are read from materialized
#' fields, and shown as `NA` when the path summary is not enriched.
#'
#' @param x A `gp_fit`.
#' @param ... Unused.
#' @param width Optional integer rendering width.
#' @return A `character` vector, one element per line.
#'
#' @examples
#' # The pre-solved tiny fit bundled with the package (instant; no solve).
#' fit <- readRDS(system.file("extdata/examples/tiny_fit.rds", package = "gradepath"))
#' writeLines(format_gp_fit_cli(fit))
#'
#' @seealso [print.gp_fit()], [krw_report_card()], [get_report_card()]
#' @family gradepath-cli
#' @export
format_gp_fit_cli <- function(x, ..., width = NULL) {
  n_units <- length(x$ids)
  grades <- x$selected_grade$assignment$grade
  composition <- .gp_grade_composition(grades)
  sel <- x$grade_path$selection$selected_lambda %gp_or% NA_real_
  metrics <- .gp_path_selected_metrics(x$grade_path)
  backend <- x$control$backend %gp_or% x$grade_path$backend$name %gp_or% "NA"

  reliability_line <- if (is.null(metrics)) {
    "reliability: (1 - DR) = NA   tau-bar = NA"
  } else {
    sprintf(
      "reliability: (1 - DR) = %s   tau-bar = %s",
      .gp_fmt_num(metrics$reliability, 2L),
      .gp_fmt_num(metrics$tau_bar, 2L)
    )
  }

  title <- sprintf(
    "gp_fit  .  %d units  .  grades: %s  .  selected lambda = %s",
    n_units, composition, .gp_fmt_num(sel, 2L)
  )
  body <- c(
    sprintf("units      : %d", n_units),
    sprintf("grades     : %s", composition),
    reliability_line,
    sprintf("backend    : %s", as.character(backend))
  )
  c(
    .gp_cli_box(title, body, width = width),
    .gp_cli_note("summary(fit) for backend/selection details; get_report_card(fit) for the ranked table.", width = width)
  )
}

#' @rdname format_gp_fit_cli
#' @param object A `gp_fit`.
#' @return `summary()` returns a typed `gp_fit_summary`.
#' @export
summary.gp_fit <- function(object, ...) {
  grades <- object$selected_grade$assignment$grade
  metrics <- .gp_path_selected_metrics(object$grade_path)

  units <- length(object$ids)
  grade_count <- length(unique(grades))
  composition <- .gp_grade_composition(grades)
  dist <- sub("^[0-9]+ ", "", composition)   # the bare "(n1/n2/...)" distribution
  sel <- object$grade_path$selection$selected_lambda %gp_or% NA_real_
  reliability <- if (is.null(metrics)) NA_real_ else metrics$reliability # (1 - DR)
  tau_bar <- if (is.null(metrics)) NA_real_ else metrics$tau_bar

  # Interpretation: built from this object's values; NA-robust (each optional
  # sentence is gated on is.finite, so the prose never prints "NA" mid-sentence).
  interp <- sprintf("%d units sorted into %d grades %s at lambda = %s.",
                    units, grade_count, dist, .gp_fmt_num(sel, 2L))
  if (is.finite(reliability)) {
    interp <- c(interp, sprintf(paste0(
      "Grades agree with the pairwise (posterior) ranking on %d%% of comparable ",
      "pairs (reliability = 1 - discordance rate = %s); higher is cleaner."),
      round(100 * reliability), .gp_fmt_num(reliability, 3L)))
  }
  if (is.finite(tau_bar)) {
    interp <- c(interp, sprintf(paste0(
      "The grades carry tau-bar = %s of posterior Kendall rank agreement across ",
      "grade boundaries (in [-1, 1]); higher means more rank information."),
      .gp_fmt_num(tau_bar, 3L)))
  }

  .gp_new_summary(
    "gp_fit",
    units = units,
    grade_count = grade_count,
    composition = composition,
    selected_lambda = sel,
    reliability = reliability, # (1 - DR)
    tau_bar = tau_bar,
    interpretation = interp,
    glossary = c(
      reliability =
        "1 - discordance rate: share of comparable pairs the grades rank the same way as the posterior.",
      tau_bar =
        "posterior expected Kendall rank agreement over between-grade blocks, in [-1, 1]; higher = more rank information carried by the grading.",
      selected_lambda =
        "the grade-count selection knob (KRW baseline 0.25); larger lambda yields more, finer grades (smaller lambda pools units into fewer, coarser grades)."
    ),
    provenance = list(
      backend = object$control$backend %gp_or% object$grade_path$backend$name %gp_or% NA_character_,
      selection_rule = object$grade_path$selection$selection_rule %gp_or% NA_character_,
      frontier_source = if (is.null(metrics)) "not materialized" else "grade_path$summary (enriched)"
    )
  )
}

#' @rdname format_gp_fit_cli
#' @return `print()` returns `x` invisibly.
#' @export
print.gp_fit <- function(x, ..., width = getOption("width")) {
  cat(format_gp_fit_cli(x, width = width), sep = "\n")
  invisible(x)
}
