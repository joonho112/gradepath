registry_test_source_root <- function() {
  candidate_roots <- c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", ".."),
    file.path(getwd(), "..", "00_pkg_src", "gradepath"),
    file.path(getwd(), "..", "..", "00_pkg_src", "gradepath"),
    file.path(system.file(package = "gradepath"), ".."),
    file.path(system.file(package = "gradepath"), "..", "00_pkg_src", "gradepath")
  )
  candidate_roots <- unique(candidate_roots[file.exists(candidate_roots)])
  candidate_roots <- normalizePath(candidate_roots, winslash = "/", mustWork = TRUE)

  required <- c(
    "DESCRIPTION",
    file.path("data", "gp_registry.rda"),
    file.path("inst", "extdata", "registry", "paper_values.csv")
  )
  scores <- vapply(
    candidate_roots,
    function(root) {
      sum(file.exists(file.path(root, required))) +
        2L * file.exists(file.path(root, "data-raw", "make-registry.R"))
    },
    integer(1)
  )
  candidate_roots[[which.max(scores)]]
}

test_that("registry regeneration uses bundled source provenance", {
  root <- registry_test_source_root()
  source_path <- file.path("inst", "extdata", "registry", "paper_values.csv")
  source_abs <- file.path(root, source_path)
  expect_true(file.exists(source_abs))

  script_abs <- file.path(root, "data-raw", "make-registry.R")
  if (file.exists(script_abs)) {
    script <- readLines(script_abs, warn = FALSE)
    expect_false(any(grepl("/Users/joonholee", script, fixed = TRUE)))
    expect_false(any(grepl("krw-2024-companion/_companion", script, fixed = TRUE)))
  } else {
    expect_false(file.exists(script_abs))
  }

  load(file.path(root, "data", "gp_registry.rda"))
  src <- attr(gp_registry, "source", exact = TRUE)
  expect_named(src, c("name", "path", "sha256", "rows", "override_env"))
  expect_identical(src$path, source_path)
  expect_match(src$sha256, "^[0-9a-f]{64}$")

  paper_values <- utils::read.csv(
    source_abs,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  expect_equal(src$rows, nrow(paper_values))
  expect_identical(src$override_env, NA_character_)
  expect_identical(
    src$sha256,
    digest::digest(source_abs, file = TRUE, algo = "sha256")
  )
})

test_that("registry regeneration failures are actionable", {
  rscript <- file.path(R.home("bin"), "Rscript")
  root <- registry_test_source_root()
  script_abs <- file.path(root, "data-raw", "make-registry.R")
  source_abs <- file.path(root, "inst", "extdata", "registry", "paper_values.csv")

  if (file.exists(script_abs)) {
    missing_expr <- sprintf(
      "setwd(%s); source('data-raw/make-registry.R')",
      shQuote(root)
    )
    missing_source <- suppressWarnings(system2(
      rscript,
      c("-e", shQuote(missing_expr)),
      stdout = TRUE,
      stderr = TRUE,
      env = c("GRADEPATH_PAPER_VALUES_CSV=/definitely/missing/paper_values.csv"),
      wait = TRUE
    ))
    expect_true(!is.null(attr(missing_source, "status")))
    expect_true(any(grepl("GRADEPATH_PAPER_VALUES_CSV", missing_source, fixed = TRUE)))

    not_root_expr <- sprintf(
      "setwd(tempdir()); source(%s)",
      shQuote(script_abs)
    )
    not_root <- suppressWarnings(system2(
      rscript,
      c("-e", shQuote(not_root_expr)),
      stdout = TRUE,
      stderr = TRUE,
      wait = TRUE
    ))
    expect_true(!is.null(attr(not_root, "status")))
    expect_true(any(grepl("package root", not_root, fixed = TRUE)))
  } else {
    expect_true(file.exists(file.path(root, "DESCRIPTION")))
    expect_true(file.exists(source_abs))
    expect_false(file.exists(script_abs))
  }
})
