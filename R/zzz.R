# Package load/attach hooks for gradepath.
#
# Split the CRAN-safe way:
#   .onLoad   sets default options only -- no printing, no network, no solver
#             probe. Safe to run during R CMD check, vignette builds, and when
#             another package merely imports gradepath.
#   .onAttach prints user-facing guidance (the Gurobi prerequisite) via
#             packageStartupMessage(), which is suppressible with
#             suppressPackageStartupMessages() and never fires on a bare
#             requireNamespace("gradepath").

# Default run-time options, set only if the user has not already set them.
# `gradepath.backend` mirrors the gp_control() backend default (Chapters 9 and
# 13); gp_control() remains the single source of truth, and this option exists
# so internal helpers and the startup notice can read a backend default before
# a control object is constructed. The build/caching scripts read it via
# getOption("gradepath.backend", ...).
.gradepath_default_backend <- function() {
  env_backend <- Sys.getenv("GRADEPATH_BACKEND", unset = "")
  if (nzchar(env_backend)) {
    return(tolower(env_backend))
  }
  "gurobi"
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  defaults <- list(gradepath.backend = .gradepath_default_backend())
  toset <- !(names(defaults) %in% names(op))
  if (any(toset)) {
    options(defaults[toset])
  }
  invisible()
}

.onAttach <- function(libname, pkgname) {
  # Soft, non-fatal reference to the imported primitive supplier. ebrecipe is
  # a hard Imports dependency (the package will not load without it), so this
  # only matters for a broken/partial install; it also documents the boundary
  # at attach time. Kept as requireNamespace() rather than a version test so no
  # base-utils dependency is pulled in at scaffold stage.
  if (!requireNamespace("ebrecipe", quietly = TRUE)) {
    packageStartupMessage(
      "gradepath: the primitive supplier 'ebrecipe' (>= 0.5.0) is not ",
      "available; estimation will fail until it is installed."
    )
  }

  # Gurobi-prerequisite guidance (the drrank model). The default backend is
  # "gurobi", the same solver KRW used and a required prerequisite for the
  # default replication path, but only a Suggests dependency for packaging --
  # so the package loads, installs, and tests without a license. We do NOT
  # error here; the install-pointing error belongs to the solver verbs
  # (Chapter 10) at call time. We only inform, and only when the default
  # backend is selected and Gurobi is not detectable.
  backend <- getOption("gradepath.backend", "gurobi")
  if (identical(backend, "gurobi")) {
    have_gurobi <-
      requireNamespace("gurobi", quietly = TRUE) ||
      nzchar(Sys.which("gurobi_cl"))
    if (!have_gurobi) {
      packageStartupMessage(
        "gradepath: Gurobi is the strongly recommended default backend (KRW's ",
        "solver; fastest on the grade integer program), but no Gurobi ",
        "installation was detected\n",
        "  (neither the 'gurobi' R package nor 'gurobi_cl' on PATH).\n",
        "  Install Gurobi + a free academic license: ",
        "https://www.gurobi.com/academia/\n",
        "  License-free last resort (several times slower on a real 97-firm ",
        "solve): gp_control(backend = \"highs\")."
      )
    }
  }

  invisible()
}
