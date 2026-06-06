# Tests for the accessor family and base generics over gp_fit.
#
# The valid gp_fit is built with the SAME construction as
# tests/testthat/test-class-fit.R (copied here so the test is self-contained):
# a real ebrecipe::eb_estimates when ebrecipe is installed, the native
# gp_prior / gp_posterior / gp_precision_fit shells, and the decision
# constructors. Run after sourcing the real R/ + accessors.R + testthat.

# ===========================================================================
# fixtures  (copied verbatim from tests/testthat/test-class-fit.R)
# ===========================================================================

GP_IDS  <- c("u1", "u2", "u3", "u4")
GP_J    <- length(GP_IDS)
GP_LAM  <- c(0.00, 0.25, 0.50, 1.00)
GP_SELL <- 0.25
GP_GRADES <- c(1L, 1L, 2L, 3L)   # the assignment at the selected lambda

make_estimates <- function(J = GP_J, ids = GP_IDS, krw_gmm = FALSE) {
  est <- seq_len(J) + 0.1
  se  <- rep(0.3, J)
  if (requireNamespace("ebrecipe", quietly = TRUE)) {
    e <- ebrecipe::eb_input(theta_hat = est, s = se, unit_id = ids)
    if (krw_gmm) {
      e$metadata$original_s <- se
    }
    return(e)
  }
  fields <- list(
    estimate = est, se = se, id = ids, label = ids,
    covariate = NULL, n = NULL, weight = NULL,
    metadata = list(), diagnostics = list()
  )
  if (krw_gmm) fields$original_s <- se
  structure(fields, class = c("eb_estimates", "list"))
}

make_prior <- function() {
  supp <- seq(-3, 3, length.out = 50)
  d    <- dnorm(supp); d <- d / sum(d)
  new_gp_prior(support = supp, density = d, mean = 0, scale = "r")
}

make_posterior <- function(ids = GP_IDS, J = GP_J) {
  pm <- seq_len(J) + 0.0
  new_gp_posterior(
    estimate = seq_len(J) + 0.1, se = rep(0.3, J),
    id = ids, label = ids,
    posterior_mean = pm, posterior_sd = rep(0.2, J),
    lower = pm - 0.5, upper = pm + 0.5, scale = "r"
  )
}

make_precision_fit <- function() {
  new_gp_precision_fit(
    parameters = list(model_form = "krw_4moment_2step", beta = 0.51, mu = 0.0),
    moments    = list(m_hat = c(0, 1, 0, 3)),
    diagnostics = list(J_stat = 1.2)
  )
}

make_pairwise <- function(ids = GP_IDS, J = GP_J) {
  m <- matrix(0.5, J, J); dimnames(m) <- list(ids, ids)
  for (i in seq_len(J)) for (j in seq_len(J)) if (i != j) {
    m[i, j] <- if (i < j) 0.6 else 0.4
  }
  new_gp_pairwise(
    ids = ids, matrix = m, power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = 0.5, zero_floor = 1e-7),
    source  = list(stage = "posterior", rule = "outer_product",
                   assumption = "one_level_independence"),
    control = gp_control()
  )
}

make_grade_fit <- function(ids = GP_IDS, lambda = GP_SELL, grades = GP_GRADES) {
  assignment <- data.frame(id = ids, grade = grades, stringsAsFactors = FALSE)
  new_gp_grade_fit(
    ids = ids, lambda = lambda, assignment = assignment,
    summary  = list(grade_count = length(unique(grades))),
    objective = list(value = 0.5),
    backend  = list(name = "gurobi"),
    control  = gp_control()
  )
}

make_grade_path <- function(ids = GP_IDS, lambdas = GP_LAM, sel = GP_SELL) {
  fits <- lapply(lambdas, function(l) {
    g <- if (isTRUE(all.equal(l, sel))) GP_GRADES else c(1L, 1L, 1L, 2L)
    make_grade_fit(ids = ids, lambda = l, grades = g)
  })
  summ <- data.frame(
    lambda = lambdas,
    grade_count = vapply(fits, function(f) f$summary$grade_count, integer(1))
  )
  new_gp_grade_path(
    ids = ids, lambda_grid = lambdas, fits = fits, summary = summ,
    backend = list(name = "gurobi"),
    selection = list(selected_lambda = sel, selection_rule = "default",
                     endpoint_lambda = lambdas[length(lambdas)]),
    control = gp_control()
  )
}

make_report_card <- function(ids = GP_IDS, J = GP_J, grades = GP_GRADES,
                             sel = GP_SELL, posterior = make_posterior()) {
  tab <- data.frame(
    id = ids, label = ids, grade = grades,
    sort_rank = seq_len(J),
    selected_lambda = rep(sel, J),
    posterior_mean = posterior$posterior_mean,
    lower = posterior$lower, upper = posterior$upper,
    estimate = seq_len(J) + 0.1, se = rep(0.3, J),
    stringsAsFactors = FALSE
  )
  new_gp_report_card(ids = ids, table = tab, selected_lambda = sel,
                     grades = grades, control = gp_control())
}

make_fit <- function(precision_rule = "none", krw_gmm = (precision_rule == "krw_gmm")) {
  post <- make_posterior()
  ctrl <- gp_control(precision_rule = precision_rule)
  pf   <- if (precision_rule == "none") NULL else make_precision_fit()
  path <- make_grade_path()
  selected <- path$fits[[which(abs(path$lambda_grid - path$selection$selected_lambda) < 1e-8)]]
  new_gp_fit(
    ids            = GP_IDS,
    estimates      = make_estimates(krw_gmm = krw_gmm),
    prior          = make_prior(),
    posterior      = post,
    precision_fit  = pf,
    pairwise       = make_pairwise(),
    grade_path     = path,
    selected_grade = selected,
    report_card    = make_report_card(posterior = post),
    control        = ctrl
  )
}

# ===========================================================================
# get_grades — integer vector NAMED by ids, from selected_grade
# ===========================================================================

test_that("get_grades returns an integer vector named by ids", {
  fit <- make_fit("none")
  g <- get_grades(fit)

  expect_type(g, "integer")
  expect_length(g, GP_J)
  expect_identical(names(g), GP_IDS)                 # names == canonical ids
  expect_identical(unname(g), GP_GRADES)             # values == selected grades
})

test_that("get_grades reads the SELECTED grade_fit assignment (aligned to ids)", {
  fit <- make_fit("none")
  g <- get_grades(fit)
  # the grade per id matches selected_grade$assignment, aligned by id
  asg <- fit$selected_grade$assignment
  expect_identical(
    g,
    setNames(as.integer(asg$grade[match(GP_IDS, asg$id)]), GP_IDS)
  )
})

test_that("get_grades aligns by id even when assignment row order is permuted", {
  fit <- make_fit("none")
  # Permute the assignment rows; the accessor must realign to canonical ids and
  # still return grades in ids order (selected_grade$ids stays canonical).
  perm <- c(4L, 1L, 3L, 2L)
  fit$selected_grade$assignment <- fit$selected_grade$assignment[perm, , drop = FALSE]
  g <- get_grades(fit)
  expect_identical(names(g), GP_IDS)
  expect_identical(unname(g), GP_GRADES)             # u1->1,u2->1,u3->2,u4->3
})

# ===========================================================================
# coef.gp_fit == get_grades
# ===========================================================================

test_that("coef.gp_fit equals get_grades", {
  fit <- make_fit("none")
  expect_identical(coef(fit), get_grades(fit))
  expect_type(coef(fit), "integer")
  expect_identical(names(coef(fit)), GP_IDS)
})

# ===========================================================================
# slot accessors return the right slot/type
# ===========================================================================

test_that("get_report_card returns the gp_report_card slot", {
  fit <- make_fit("none")
  rc <- get_report_card(fit)
  expect_identical(rc, fit$report_card)
  expect_s3_class(rc, "gp_report_card")
})

test_that("get_pairwise returns the gp_pairwise slot", {
  fit <- make_fit("none")
  pw <- get_pairwise(fit)
  expect_identical(pw, fit$pairwise)
  expect_s3_class(pw, "gp_pairwise")
})

test_that("get_prior returns the gp_prior slot (eb_prior-shaped)", {
  fit <- make_fit("none")
  pr <- get_prior(fit)
  expect_identical(pr, fit$prior)
  expect_s3_class(pr, "gp_prior")           # SHELL class; eb_* round-trip
})

test_that("get_posterior returns the gp_posterior slot (eb_posterior-shaped)", {
  fit <- make_fit("none")
  po <- get_posterior(fit)
  expect_identical(po, fit$posterior)
  expect_s3_class(po, "gp_posterior")       # SHELL class; eb_* round-trip
})

test_that("get_control returns the gp_control slot", {
  fit <- make_fit("none")
  ct <- get_control(fit)
  expect_identical(ct, fit$control)
  expect_s3_class(ct, "gp_control")
})

# ===========================================================================
# as.data.frame.gp_fit == report_card$table (row.names reset)
# ===========================================================================

test_that("as.data.frame.gp_fit returns the report-card table with reset row names", {
  fit <- make_fit("none")
  df <- as.data.frame(fit)
  expect_s3_class(df, "data.frame")

  expected <- fit$report_card$table
  rownames(expected) <- NULL
  expect_identical(df, expected)

  # exactly the report_card$table content (ignoring row names)
  expect_equal(df, get_report_card(fit)$table, ignore_attr = "row.names")
  expect_identical(nrow(df), GP_J)
})

# ===========================================================================
# non-gp_fit input errors cleanly (every accessor + base generics)
# ===========================================================================

test_that("get_* accessors reject non-gp_fit input cleanly", {
  bad <- list(a = 1)                      # plain list, not a gp_fit
  expect_error(get_grades(bad),      "gp_fit")
  expect_error(get_report_card(bad), "gp_fit")
  expect_error(get_pairwise(bad),    "gp_fit")
  expect_error(get_prior(bad),       "gp_fit")
  expect_error(get_posterior(bad),   "gp_fit")
  expect_error(get_control(bad),     "gp_fit")

  # also the bare-atomic case
  expect_error(get_grades(42),       "gp_fit")
  expect_error(get_grades("nope"),   "gp_fit")
})

test_that("the get_* error is a classed gradepath condition", {
  expect_error(get_grades(list()), class = "gradepath_error")
})

test_that("coef / as.data.frame on a non-gp_fit dispatch elsewhere (no gp_fit method hijack)", {
  # coef.gp_fit / as.data.frame.gp_fit must NOT fire for non-gp_fit inputs:
  # base dispatch should handle a plain data.frame / lm-less object normally.
  df <- data.frame(x = 1:3)
  expect_identical(as.data.frame(df), df)             # base method, untouched
})
