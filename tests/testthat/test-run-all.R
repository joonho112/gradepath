test_that("gp_run_all() scores named values and summarizes verdicts", {
  vals <- c(
    race_baseline_ngrades = 3,
    race_baseline_dr = 3.9,
    race_baseline_tau = 0.50
  )
  sc <- gp_run_all(vals)

  expect_s3_class(sc, "gp_scorecard")
  expect_identical(names(sc), c("checks", "table", "summary", "pass_rate", "provenance", "warnings"))
  expect_identical(sc$checks$status, c("PASS", "PASS", "FAIL"))
  expect_equal(sc$pass_rate, 2 / 3)
  total <- sc$summary[sc$summary$group == "TOTAL", ]
  expect_equal(total$pass, 2)
  expect_equal(total$fail, 1)
  expect_equal(total$unverified, 0)
})

test_that("gp_validate_targets() accepts data frames and fills missing requested targets", {
  replicated <- data.frame(
    id = "race_baseline_ngrades",
    replicated = 3,
    group = "M1 grades",
    stringsAsFactors = FALSE
  )
  sc <- gp_validate_targets(
    replicated,
    targets = c("race_baseline_ngrades", "race_baseline_dr")
  )

  expect_identical(sc$checks$id, c("race_baseline_ngrades", "race_baseline_dr"))
  expect_identical(sc$checks$status, c("PASS", "UNVERIFIED"))
  expect_identical(sc$checks$reason[[2]], "UNVERIFIED")
  expect_true(any(sc$summary$unverified > 0))
})

test_that("gp_run_all(NULL) returns requested targets as UNVERIFIED without heavy work", {
  testthat::local_mocked_bindings(
    gp_grade_path = function(...) stop("solver should not run in gp_run_all(NULL)"),
    krw_report_card = function(...) stop("monolith should not run in gp_run_all(NULL)"),
    .package = "gradepath"
  )

  sc <- gp_run_all(targets = c("race_baseline_ngrades", "race_baseline_dr"))
  expect_identical(sc$checks$status, c("UNVERIFIED", "UNVERIFIED"))
  expect_true(is.na(sc$pass_rate))
})

test_that("cache helper loads hits and reports stale metadata", {
  e <- new.env(parent = baseenv())
  sys.source(
    system.file("scripts", "build-cached-assets.R", package = "gradepath"),
    envir = e
  )
  dir <- tempfile("gp-cache-")
  counter <- new.env(parent = emptyenv())
  counter$count <- 0L
  first <- e$load_or_compute("toy", {
    counter$count <- counter$count + 1L
    list(value = counter$count)
  }, cache_dir = dir, seed = 1L, metadata = list(source_hash = "a"))
  second <- e$load_or_compute("toy", {
    counter$count <- counter$count + 1L
    list(value = counter$count)
  }, cache_dir = dir, seed = 1L, metadata = list(source_hash = "a"))
  first_value <- first
  second_value <- second
  attr(first_value, "build_metadata") <- NULL
  attr(second_value, "build_metadata") <- NULL
  expect_identical(first_value, list(value = 1L))
  expect_identical(second_value, list(value = 1L))
  expect_named(attr(second, "build_metadata"), c(
    "key",
    "seed",
    "built_at",
    "r_version",
    "gradepath_version",
    "backend",
    "solver_metadata",
    "source_hash",
    "extra"
  ))
  expect_match(attr(second, "build_metadata")$r_version, "^\\d+\\.\\d+\\.\\d+")
  expect_named(attr(second, "build_metadata")$solver_metadata, c(
    "backend_env",
    "gurobi_cl_available",
    # CCR-16: solver binary VERSION is part of the compatibility-checked
    # metadata so a Gurobi upgrade that changes MIP behaviour invalidates the
    # cache; it is NA_character_ when the solver is unavailable.
    "gurobi_version",
    "gurobi_r_package_available",
    "highs_package_available",
    "roi_package_available"
  ))
  expect_true(is.character(
    attr(second, "build_metadata")$solver_metadata$gurobi_version
  ))
  expect_length(
    attr(second, "build_metadata")$solver_metadata$gurobi_version,
    1L
  )
  expect_identical(counter$count, 1L)
  expect_error(
    e$load_or_compute("toy", 99, cache_dir = dir, seed = 1L,
                      metadata = list(source_hash = "b"), on_stale = "error"),
    regexp = "CACHE_STALE"
  )

  old <- list(value = 1L)
  attr(old, "build_metadata") <- list(key = "old", seed = 1L)
  expect_false(e$cache_metadata_compatible(
    old,
    e$cache_metadata("old", seed = 1L, extra = list(source_hash = "a"))
  ))
})

test_that("M1 cache dependency manifest is deterministic and portable", {
  e <- new.env(parent = baseenv())
  sys.source(
    system.file("scripts", "build-cached-assets.R", package = "gradepath"),
    envir = e
  )

  manifest <- e$cache_m1_dependency_manifest()
  expect_s3_class(manifest, "data.frame")
  expect_named(manifest, c("path", "category", "reason"))
  expect_identical(manifest$path, sort(manifest$path))
  expect_identical(manifest$path, unique(manifest$path))
  expect_false(any(grepl("^(/|~)", manifest$path)))
  expect_false(any(grepl("\\\\", manifest$path)))
  expect_true(all(nzchar(manifest$category)))
  expect_true(all(nzchar(manifest$reason)))

  candidate_roots <- c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", ".."),
    file.path(getwd(), "..", "00_pkg_src", "gradepath"),
    file.path(system.file(package = "gradepath"), ".."),
    file.path(system.file(package = "gradepath"), "..", "00_pkg_src", "gradepath")
  )
  candidate_roots <- unique(candidate_roots[file.exists(candidate_roots)])
  candidate_roots <- normalizePath(candidate_roots, winslash = "/", mustWork = TRUE)
  hit_counts <- vapply(
    candidate_roots,
    function(root) sum(file.exists(file.path(root, manifest$path))),
    integer(1)
  )
  package_root <- candidate_roots[[which.max(hit_counts)]]
  manifest_exists <- file.exists(file.path(package_root, manifest$path))
  missing_manifest_paths <- manifest$path[!manifest_exists]
  # CCR-18: the manifest must be fully present in source checkouts, but two
  # standard install shapes legitimately strip inputs: a built tarball drops
  # `data-raw/` (devtools::build()), and a binary/installed package drops
  # `tests/`. Tolerate ONLY those documented omissions -- any other missing
  # manifest path is a real portability regression.
  is_portably_strippable <- function(p) {
    p == "data-raw/make-registry.R" |
      grepl("^tests/testthat/.*\\.R$", p)
  }
  expect_true(
    all(manifest_exists) ||
      all(is_portably_strippable(missing_manifest_paths))
  )
  expect_identical(e$cache_m1_dependency_paths(), manifest$path)
  if (all(manifest_exists)) {
    expect_identical(
      e$cache_m1_source_hash(root = package_root),
      e$cache_source_hash(manifest$path, root = package_root)
    )
  } else {
    expect_identical(
      e$cache_source_hash(manifest$path, root = package_root, allow_missing = TRUE),
      e$cache_source_hash(
        manifest$path[manifest_exists],
        root = package_root
      )
    )
  }

  required <- c(
    "data/gp_registry.rda",
    "data/krw_firms.rda",
    "data-raw/make-registry.R",
    "R/harness-scorecard.R",
    "R/harness-check.R",
    "R/status.R",
    "R/grade.R",
    "R/grade-backends.R",
    "R/pairwise.R",
    "R/frontier.R",
    "R/report-card.R",
    "R/krw-report-card.R",
    "R/krw-report-card-targets.R",
    "R/krw-r2.R",
    # CCR-02: the beta-GMM / W-seam / deconvolution chain is LIVE-computed when
    # the M1 cache is built (t3_race_ni_beta via gp_w_seam()), so it must be in
    # the manifest the source hash watches.
    "R/seam-ebrecipe.R",
    "R/estimation-coupling.R",
    "R/estimation-input.R",
    "R/estimation-core.R",
    "R/deconvolution.R",
    "R/posterior-onelevel.R",
    "R/posterior-weights.R",
    "R/utils-math.R",
    # utils-validate.R defines the coercion/validation
    # primitives and the provenance stamp that gate the live acceptance solve's
    # producer_status (which a cached PASS depends on); zzz.R resolves the default
    # solver backend. Both drive the injected race+gender solve, so the source
    # hash must watch them.
    "R/utils-validate.R",
    "R/zzz.R",
    "inst/scripts/diagnose-race-n97-selected.R",
    "inst/scripts/diagnose-gender-deconvolution.R",
    "inst/scripts/diagnose-krw-r2.R",
    "inst/scripts/diagnose-report-card-targets.R",
    "inst/extdata/acceptance/m1-scorecard.csv",
    "inst/extdata/acceptance/m1-solver-evidence.csv",
    "inst/extdata/acceptance/m1-race-selected-diagnostic.csv",
    "inst/extdata/acceptance/m1-gender-deconvolution-diagnostic.csv",
    "inst/extdata/acceptance/m1-names-scope.csv",
    "inst/extdata/acceptance/m1-krw-r2-diagnostic.csv",
    "inst/extdata/acceptance/m1-report-card-targets.csv",
    "inst/extdata/registry/paper_values.csv",
    "tests/testthat/test-gurobi-robust.R",
    "tests/testthat/test-m1-acceptance.R",
    "tests/testthat/test-registry-regeneration.R"
  )
  expect_true(all(required %in% manifest$path))
})

test_that("mutating an estimation-chain file changes the M1 source hash (CCR-02)", {
  # CCR-02 regression: t3_race_ni_beta is LIVE-computed via gp_w_seam() when the
  # M1 cache is built, so the beta-GMM / W-seam / deconvolution chain is now in
  # the manifest. Prove the source hash is content-sensitive to that chain by
  # mutating R/estimation-coupling.R and confirming cache_m1_source_hash()
  # changes. The mutation happens on a TEMP COPY of the manifest sources, so no
  # tracked source file is ever touched and there is nothing to restore.
  e <- new.env(parent = baseenv())
  sys.source(
    system.file("scripts", "build-cached-assets.R", package = "gradepath"),
    envir = e
  )

  manifest_paths <- e$cache_m1_dependency_paths()
  expect_true("R/estimation-coupling.R" %in% manifest_paths)

  # Resolve the package source root that holds the most manifest files.
  candidate_roots <- c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", ".."),
    file.path(getwd(), "..", "00_pkg_src", "gradepath"),
    file.path(system.file(package = "gradepath"), ".."),
    file.path(system.file(package = "gradepath"), "..", "00_pkg_src", "gradepath")
  )
  candidate_roots <- unique(candidate_roots[file.exists(candidate_roots)])
  candidate_roots <- normalizePath(candidate_roots, winslash = "/", mustWork = TRUE)
  hit_counts <- vapply(
    candidate_roots,
    function(root) sum(file.exists(file.path(root, manifest_paths))),
    integer(1)
  )
  package_root <- candidate_roots[[which.max(hit_counts)]]

  mutated_rel <- "R/estimation-coupling.R"
  src_file <- file.path(package_root, mutated_rel)
  # In a binary/tarball install the live R source may be absent; the hash logic
  # is identical regardless of which manifest file is mutated, so skip cleanly.
  skip_if_not(file.exists(src_file), "estimation-coupling.R source not available")

  # Stage every available manifest file into an isolated temp root so the
  # mutation only ever touches the copy.
  staging <- tempfile("gp-mutate-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  present <- manifest_paths[file.exists(file.path(package_root, manifest_paths))]
  for (rel in present) {
    dest <- file.path(staging, rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(package_root, rel), dest, overwrite = TRUE)
  }

  hash_before <- e$cache_source_hash(present, root = staging)

  staged_file <- file.path(staging, mutated_rel)
  cat("\n# CCR-02 regression marker: semantic-drift probe\n",
      file = staged_file, append = TRUE)
  hash_after <- e$cache_source_hash(present, root = staging)

  expect_false(identical(hash_before, hash_after))

  # The original tracked source file must be byte-for-byte untouched.
  expect_true(file.exists(src_file))
})

test_that("cache source hashes are deterministic and content-sensitive", {
  e <- new.env(parent = baseenv())
  sys.source(
    system.file("scripts", "build-cached-assets.R", package = "gradepath"),
    envir = e
  )
  dir <- tempfile("gp-hash-")
  dir.create(dir)
  writeLines("alpha", file.path(dir, "a.txt"))
  writeLines("beta", file.path(dir, "b.txt"))

  paths <- c("a.txt", "b.txt")
  h1 <- e$cache_source_hash(paths, root = dir)
  h2 <- e$cache_source_hash(rev(paths), root = dir)
  expect_identical(h1, h2)

  writeLines("beta changed", file.path(dir, "b.txt"))
  h3 <- e$cache_source_hash(paths, root = dir)
  expect_false(identical(h1, h3))

  expect_error(
    e$cache_source_hash(c(paths, "missing.txt"), root = dir),
    regexp = "Missing cache dependency"
  )
})

test_that("cache source hash is locale-independent (radix order)", {
  # Regression for the locale-collation bug: cache_source_hash() must order the
  # manifest paths with method = "radix" so a cache built under one LC_COLLATE
  # (e.g. en_US) and verified under another (e.g. the C locale testthat and many
  # CRAN/CI checks impose) yields the SAME digest. Plain order() collates by
  # locale, which silently reorders the named hash vector across environments.
  e <- new.env(parent = baseenv())
  sys.source(
    system.file("scripts", "build-cached-assets.R", package = "gradepath"),
    envir = e
  )
  dir <- tempfile("gp-locale-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
  # Paths chosen to collate differently under C vs a UTF-8 locale: uppercase "R"
  # (ASCII 82) sorts before lowercase "data" (100) under C, but most UTF-8
  # collations fold case and order "data*" first.
  rel <- c("R/api.R", "data-raw/make.R", "data/x.txt", "tests/t.R")
  for (p in rel) {
    dir.create(file.path(dir, dirname(p)), recursive = TRUE, showWarnings = FALSE)
    writeLines(p, file.path(dir, p))
  }

  saved <- Sys.getlocale("LC_COLLATE")
  on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", saved)), add = TRUE)

  set_ok <- function(loc) {
    isTRUE(tryCatch(
      nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", loc))),
      error = function(e) FALSE
    ))
  }

  hashes <- character(0)
  for (loc in c("C", "en_US.UTF-8", "en_US")) {
    if (set_ok(loc)) {
      hashes[[loc]] <- e$cache_source_hash(rel, root = dir)
    }
  }
  # C must always be available; need it plus at least one other to be meaningful.
  expect_true("C" %in% names(hashes))
  if (length(hashes) >= 2L) {
    expect_identical(length(unique(unlist(hashes))), 1L)
  } else {
    skip("no non-C LC_COLLATE locale available to contrast")
  }
})

test_that("gp_validate_targets() extracts light M1 facts from a gp_fit", {
  ids <- paste0("u", seq_len(4))
  ctl <- gp_control(lambda_grid = c(0.25, 1), backend = "highs")
  pair <- new_gp_pairwise(
    ids = ids,
    matrix = matrix(c(
      0.5, 0.8, 0.8, 0.8,
      0.2, 0.5, 0.8, 0.8,
      0.2, 0.2, 0.5, 0.8,
      0.2, 0.2, 0.2, 0.5
    ), nrow = 4, byrow = TRUE, dimnames = list(ids, ids)),
    power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = 0.5, zero_floor = 1e-7),
    source = list(stage = "posterior", rule = "outer_product", assumption = "one_level_independence"),
    control = ctl,
    provenance = list()
  )
  path <- gp_grade_path(pair, control = ctl, selected_lambda = 0.25)
  selected <- gp_select_grade(path, 0.25)
  posterior <- new_gp_posterior(
    estimate = c(-1, -0.5, 0.5, 1),
    se = rep(0.1, 4),
    id = ids,
    label = ids,
    posterior_mean = c(-1, -0.5, 0.5, 1),
    posterior_sd = rep(0.1, 4),
    lower = c(-1.2, -0.7, 0.3, 0.8),
    upper = c(-0.8, -0.3, 0.7, 1.2),
    scale = "r",
    metadata = list()
  )
  estimates <- ebrecipe::eb_input(
    theta_hat = c(-1, -0.5, 0.5, 1),
    s = rep(0.1, 4),
    unit_id = ids
  )
  card <- gp_report_card(
    estimates,
    posterior = posterior,
    selected_grade = selected,
    grade_path = path
  )
  fit <- new_gp_fit(
    ids = ids,
    estimates = estimates,
    prior = new_gp_prior(
      support = c(-1, 0, 1),
      density = c(0.25, 0.5, 0.25),
      mean = 0,
      scale = "r",
      diagnostics = list(),
      metadata = list()
    ),
    posterior = posterior,
    precision_fit = NULL,
    pairwise = pair,
    grade_path = path,
    selected_grade = selected,
    report_card = card,
    control = ctl,
    # A `gp_fit` should always carry an explicit producer status; set it so the
    # base extraction asserts the producing stage was acceptance-ready ("OK").
    # (Fail-safe default of "UNVERIFIED" when absent is covered by a dedicated
    # regression test below.)
    provenance = list(demographic = "race", producer_status = "OK"),
    warnings = character(0)
  )
  fit <- validate_gp_fit(fit)

  sc <- gp_validate_targets(
    fit,
    targets = c("scale_n_firms_graded", "race_baseline_ngrades")
  )
  expect_identical(sc$checks$id, c("scale_n_firms_graded", "race_baseline_ngrades"))
  expect_identical(sc$checks$replicated[[1]], 4)
  expect_true(all(sc$checks$status %in% c("PASS", "FAIL")))

  for (status in c("SOLVER_GAP", "SOLVER_TIME_LIMIT")) {
    status_fit <- fit
    status_fit$provenance$producer_status <- status
    status_sc <- gp_validate_targets(
      status_fit,
      targets = c("scale_n_firms_graded", "race_baseline_ngrades")
    )
    expect_identical(status_sc$checks$status, c("UNVERIFIED", "UNVERIFIED"))
    expect_identical(status_sc$checks$reason, c(status, status))
    expect_identical(status_sc$checks$producer_status, c(status, status))
    expect_true(all(is.na(status_sc$checks$replicated)))
    expect_true(is.na(status_sc$pass_rate))
  }

  override_fit <- fit
  override_fit$provenance$producer_status <- "SOLVER_GAP"
  ok_override <- gp_validate_targets(
    override_fit,
    targets = c("scale_n_firms_graded", "race_baseline_ngrades"),
    producer_status = "OK"
  )
  expect_true(all(ok_override$checks$status %in% c("PASS", "FAIL")))

  r2_sc <- gp_validate_targets(fit, targets = "race_baseline_r2")
  expect_identical(r2_sc$checks$status, "UNVERIFIED")
})

# Shared toy `gp_fit` builder for the fail-safe / cross-column regressions below.
# `producer_status` controls only `provenance$producer_status`; when NULL the key
# is left absent so the harness must apply its own default.
.gp_test_toy_fit <- function(producer_status = "OK") {
  ids <- paste0("u", seq_len(4))
  ctl <- gp_control(lambda_grid = c(0.25, 1), backend = "highs")
  pair <- new_gp_pairwise(
    ids = ids,
    matrix = matrix(c(
      0.5, 0.8, 0.8, 0.8,
      0.2, 0.5, 0.8, 0.8,
      0.2, 0.2, 0.5, 0.8,
      0.2, 0.2, 0.2, 0.5
    ), nrow = 4, byrow = TRUE, dimnames = list(ids, ids)),
    power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = 0.5, zero_floor = 1e-7),
    source = list(stage = "posterior", rule = "outer_product", assumption = "one_level_independence"),
    control = ctl,
    provenance = list()
  )
  path <- gp_grade_path(pair, control = ctl, selected_lambda = 0.25)
  selected <- gp_select_grade(path, 0.25)
  posterior <- new_gp_posterior(
    estimate = c(-1, -0.5, 0.5, 1),
    se = rep(0.1, 4),
    id = ids,
    label = ids,
    posterior_mean = c(-1, -0.5, 0.5, 1),
    posterior_sd = rep(0.1, 4),
    lower = c(-1.2, -0.7, 0.3, 0.8),
    upper = c(-0.8, -0.3, 0.7, 1.2),
    scale = "r",
    metadata = list()
  )
  estimates <- ebrecipe::eb_input(
    theta_hat = c(-1, -0.5, 0.5, 1),
    s = rep(0.1, 4),
    unit_id = ids
  )
  card <- gp_report_card(
    estimates,
    posterior = posterior,
    selected_grade = selected,
    grade_path = path
  )
  provenance <- list(demographic = "race")
  if (!is.null(producer_status)) {
    provenance$producer_status <- producer_status
  }
  fit <- new_gp_fit(
    ids = ids,
    estimates = estimates,
    prior = new_gp_prior(
      support = c(-1, 0, 1),
      density = c(0.25, 0.5, 0.25),
      mean = 0,
      scale = "r",
      diagnostics = list(),
      metadata = list()
    ),
    posterior = posterior,
    precision_fit = NULL,
    pairwise = pair,
    grade_path = path,
    selected_grade = selected,
    report_card = card,
    control = ctl,
    provenance = provenance,
    warnings = character(0)
  )
  validate_gp_fit(fit)
}

test_that("a gp_fit with no provenance$producer_status routes to UNVERIFIED (fail-safe)", {
  # CCR-12 regression: a missing producer status must NOT default to the
  # fail-open "OK". A fit that does not carry an explicit, acceptance-ready
  # producer status is treated as UNVERIFIED rather than silently passing.
  fit <- .gp_test_toy_fit(producer_status = NULL)
  expect_null(fit$provenance$producer_status)

  sc <- gp_validate_targets(
    fit,
    targets = c("scale_n_firms_graded", "race_baseline_ngrades")
  )
  expect_identical(sc$checks$status, c("UNVERIFIED", "UNVERIFIED"))
  expect_identical(sc$checks$reason, c("UNVERIFIED", "UNVERIFIED"))
  expect_identical(sc$checks$producer_status, c("UNVERIFIED", "UNVERIFIED"))
  expect_true(all(is.na(sc$checks$replicated)))
  expect_true(is.na(sc$pass_rate))

  # An explicit OK override still wins, proving only the default changed.
  ok_sc <- gp_validate_targets(
    fit,
    targets = c("scale_n_firms_graded", "race_baseline_ngrades"),
    producer_status = "OK"
  )
  expect_true(all(ok_sc$checks$status %in% c("PASS", "FAIL")))

  # The internal helper itself defaults a missing fit status to UNVERIFIED.
  expect_identical(
    gradepath:::.gp_fit_producer_status(fit),
    "UNVERIFIED"
  )
  # ...but an explicitly OK-carrying fit still reports OK.
  ok_fit <- .gp_test_toy_fit(producer_status = "OK")
  expect_identical(
    gradepath:::.gp_fit_producer_status(ok_fit),
    "OK"
  )
})

test_that("gp_run_all() never marks a PASS row with a non-OK producer status (runtime)", {
  # CCR-13 regression: cross-column invariant asserted on a RUNTIME scorecard
  # object, not just the static bundled CSV. Whatever the producer status,
  # every PASS row must carry producer_status == "OK".
  invariant <- function(sc) {
    checks <- sc$checks
    sum(checks$status == "PASS" & checks$producer_status != "OK")
  }

  # (1) An OK fit is eligible for PASS/FAIL (its toy values FAIL the registry
  # band here), and every scored row carries producer_status == "OK".
  ok_fit <- .gp_test_toy_fit(producer_status = "OK")
  ok_sc <- gp_validate_targets(
    ok_fit,
    targets = c("scale_n_firms_graded", "race_baseline_ngrades")
  )
  expect_true(all(ok_sc$checks$status %in% c("PASS", "FAIL")))
  expect_true(all(ok_sc$checks$producer_status == "OK"))
  expect_identical(invariant(ok_sc), 0L)

  # (2) A non-OK fit must yield zero PASS rows at all.
  gap_fit <- .gp_test_toy_fit(producer_status = "SOLVER_GAP")
  gap_sc <- gp_validate_targets(
    gap_fit,
    targets = c("scale_n_firms_graded", "race_baseline_ngrades")
  )
  expect_false(any(gap_sc$checks$status == "PASS"))
  expect_identical(invariant(gap_sc), 0L)

  # (3) A fit missing its producer status (fail-safe) also yields no PASS rows.
  missing_fit <- .gp_test_toy_fit(producer_status = NULL)
  missing_sc <- gp_validate_targets(
    missing_fit,
    targets = c("scale_n_firms_graded", "race_baseline_ngrades")
  )
  expect_false(any(missing_sc$checks$status == "PASS"))
  expect_identical(invariant(missing_sc), 0L)

  # (4) Mixed named-value scorecard (PASS + FAIL) still upholds the invariant.
  mixed_sc <- gp_run_all(c(
    race_baseline_ngrades = 3,
    race_baseline_dr = 3.9,
    race_baseline_tau = 0.50
  ))
  expect_true(any(mixed_sc$checks$status == "PASS"))
  expect_identical(invariant(mixed_sc), 0L)
})
