# vignettes/_setup.R
# Shared preamble for every gradepath vignette.
# Each vignette sources this as its first (include = FALSE) chunk.
#
# DESIGN OF DEFENSE (Definition of Done):
#   source()-ing this file MUST NOT error, even when knitr is not loaded and
#   when no data assets are installed yet. The example fit and the parity fits
#   are shipped in later build steps. Therefore every data loader below is a
#   FUNCTION that only errors WHEN CALLED, never at source time, and the knitr
#   block is guarded by requireNamespace().
#
#   Mirrored to inst/scripts/_setup.R; reachable from an installed package via
#     source(system.file("scripts/_setup.R", package = "gradepath"))
#   Keep the two copies in sync when editing.

# ----------------------------------------------------------------
# 1. knitr chunk options (only when knitr is available)
# ----------------------------------------------------------------

if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(
    echo       = TRUE,
    comment    = "#>",
    collapse   = TRUE,
    message    = FALSE,
    warning    = FALSE,
    fig.width  = 7,
    fig.height = 4.5,
    fig.retina = 2,
    dpi        = 144,
    fig.align  = "center",
    out.width  = "100%"
  )
}

# ----------------------------------------------------------------
# 2. Package-wide RNG seed (the KRW Monte Carlo seed)
# ----------------------------------------------------------------

set.seed(15238)

# ----------------------------------------------------------------
# 3. Attach packages quietly; set the default ggplot theme if present
# ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(gradepath)
  library(ggplot2)
})

# theme_gradepath() ships with the package; guard so sourcing survives a
# partially-built package where the verb may not yet be exported.
if (exists("theme_gradepath", mode = "function")) {
  theme_set(theme_gradepath())
}

# ----------------------------------------------------------------
# 4. Helper macros (each a function)
# ----------------------------------------------------------------

# Null/empty-coalescing operator: first non-empty argument wins.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Cross-link to another article by bare slug, e.g.
#   vlink("a2-the-grading-workflow", "the grading workflow")
vlink <- function(slug, text = NULL) {
  sprintf("[%s](%s.html)", text %||% slug, slug)
}

# Cache loader: return <dir>/<slug>.rds if it exists, else evaluate `expr`
# once, persist it, and return it. `expr` is captured lazily and forced here,
# so it is only computed on a cache miss:
#   fit <- load_or_compute("a2-fit", krw_report_card(input, "race", control))
load_or_compute <- function(slug, expr, dir = ".") {
  path <- file.path(dir, paste0(slug, ".rds"))
  if (file.exists(path)) {
    return(readRDS(path))
  }
  value <- force(expr)            # lazy promise: evaluated only on a miss
  saveRDS(value, path)
  value
}

# Load the bundled proven 97-firm parity fit for one demographic (race -> grades
# 2/81/14, gender -> 1/3/89/4). Shipped in a later step, so this only errors
# WHEN CALLED with no asset present -- never at source time.
gp_parity_fit <- function(demographic = c("race", "gender")) {
  demographic <- match.arg(demographic)
  p <- system.file(
    sprintf("extdata/cached/fit_%s_parity.rds", demographic),
    package = "gradepath"
  )
  if (!nzchar(p)) {
    stop("the bundled parity fit is not installed yet (shipped in a later ",
         "build step; rebuild assets and re-install).", call. = FALSE)
  }
  readRDS(p)
}

# Load the tiny pre-solved example fit for fast, license-free live examples.
# Shipped in a later step, so this only errors WHEN CALLED with no fixture
# present -- never at source time.
gp_example_fit <- function() {
  p <- system.file("extdata/examples/tiny_fit.rds", package = "gradepath")
  if (!nzchar(p)) {
    stop("the tiny example fit is not installed yet (shipped in a later ",
         "build step).", call. = FALSE)
  }
  readRDS(p)
}
