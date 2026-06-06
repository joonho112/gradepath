# Tests for the gp_fit composite constructor + EACH cross-slot check firing,
# plus the native gp_prior / gp_posterior / gp_precision_fit shells.
#
# `estimates` is a REAL ebrecipe::eb_estimates when ebrecipe is installed
# (built via eb_input); otherwise a clearly-marked faithful stand-in.

# ===========================================================================
# fixtures
# ===========================================================================

GP_IDS  <- c("u1", "u2", "u3", "u4")
GP_J    <- length(GP_IDS)
GP_LAM  <- c(0.00, 0.25, 0.50, 1.00)
GP_SELL <- 0.25
GP_GRADES <- c(1L, 1L, 2L, 3L)   # the assignment at the selected lambda

# --- a real eb_estimates (or a faithful stand-in) --------------------------

make_estimates <- function(J = GP_J, ids = GP_IDS, krw_gmm = FALSE) {
  est <- seq_len(J) + 0.1
  se  <- rep(0.3, J)
  if (requireNamespace("ebrecipe", quietly = TRUE)) {
    # REAL eb_input signature and output use theta_hat, s, unit_id.
    e <- ebrecipe::eb_input(theta_hat = est, s = se, unit_id = ids)
    if (krw_gmm) {
      e$original_s <- se
    }
    return(e)
  }
  # ---- faithful stand-in (ebrecipe not installed) ----
  # Mirrors the 9-field eb_estimates surface exactly; original_s top-level.
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

# --- decision-layer fixtures (use the constructors) -------------------------

make_pairwise <- function(ids = GP_IDS, J = GP_J) {
  m <- matrix(0.5, J, J); dimnames(m) <- list(ids, ids)
  # antisymmetric-ish off-diagonals in [0,1]; diagonal stays 0.5
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

# --- the assembled, valid composite -----------------------------------------

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
# native shells validate on their own
# ===========================================================================

test_that("gp_prior mirrors eb_prior and validates", {
  pr <- make_prior()
  expect_silent(validate_gp_prior(pr))
  expect_identical(names(pr), .gp_prior_fields)
  expect_s3_class(pr, "gp_prior")          # SHELL class; eb_* round-trip at estimation
  bad_scale <- pr; bad_scale$scale <- "theta"
  expect_error(validate_gp_prior(bad_scale), "scale")
  bad <- pr; bad$density <- bad$density[-1] # length mismatch
  expect_error(validate_gp_prior(bad), "density")
  bad_support <- pr; bad_support$support[2] <- bad_support$support[1]
  expect_error(validate_gp_prior(bad_support), "strictly increasing")
  bad_sum <- pr; bad_sum$density <- bad_sum$density * 2
  expect_error(validate_gp_prior(bad_sum), "sum to 1")
  bad_mean <- pr; bad_mean$mean <- bad_mean$mean + 1
  expect_error(validate_gp_prior(bad_mean), "sum\\(support \\* density\\)")
})

test_that("gp_posterior mirrors eb_posterior and validates", {
  po <- make_posterior()
  expect_silent(validate_gp_posterior(po))
  expect_identical(names(po), .gp_posterior_fields)
  expect_s3_class(po, "gp_posterior")      # SHELL class; eb_* round-trip at estimation
  bad <- po; bad$lower <- bad$posterior_mean + 1   # lower > mean
  expect_error(validate_gp_posterior(bad), "lower")
  bad_est <- po; bad_est$estimate[1] <- Inf
  expect_error(validate_gp_posterior(bad_est), "estimate.*finite|finite")
  bad_se <- po; bad_se$se[1] <- -0.1
  expect_error(validate_gp_posterior(bad_se), "se.*positive|positive")
  bad_sd <- po; bad_sd$posterior_sd[1] <- -0.1
  expect_error(validate_gp_posterior(bad_sd), "posterior_sd.*non-negative|non-negative")
  bad_ci <- po; bad_ci$upper[1] <- Inf
  expect_error(validate_gp_posterior(bad_ci), "finite or NA")
})

test_that("gp_posterior tolerates NA credible intervals (eb_posterior parity)", {
  po <- make_posterior()
  po$lower <- rep(NA_real_, GP_J); po$upper <- rep(NA_real_, GP_J)
  expect_silent(validate_gp_posterior(po))
})

test_that("gp_precision_fit requires model_form/beta/mu in parameters", {
  pf <- make_precision_fit()
  expect_silent(validate_gp_precision_fit(pf))
  expect_identical(names(pf), .gp_precision_fit_fields)
  bad <- pf; bad$parameters$beta <- NULL
  expect_error(validate_gp_precision_fit(bad), "beta|missing")
})

# ===========================================================================
# valid composite passes
# ===========================================================================

test_that("a well-formed gp_fit validates (precision_rule = none)", {
  fit <- make_fit("none")
  expect_silent(validate_gp_fit(fit))
  expect_identical(validate_gp_fit(fit), fit)
  expect_identical(names(fit), .gp_fit_fields)
  expect_length(.gp_fit_fields, 13L)              # EXACT slot count
  expect_s3_class(fit, "gp_fit")
  expect_s3_class(fit$estimates, "eb_estimates")  # composition: ebrecipe owns
  expect_null(fit$precision_fit)                  # NULL under "none"
})

test_that("a well-formed gp_fit validates (precision_rule = krw_gmm)", {
  fit <- make_fit("krw_gmm")
  expect_silent(validate_gp_fit(fit))
  expect_false(is.null(fit$precision_fit))
})

# ===========================================================================
# CROSS-SLOT CHECK (a): mismatched ids
# ===========================================================================

test_that("(a) mismatched ids across slots are rejected", {
  fit <- make_fit("none")
  fit$pairwise$ids <- c("x1", "x2", "x3", "x4")
  expect_error(validate_gp_fit(fit), "pairwise ids|ids")

  fit2 <- make_fit("none")
  fit2$ids <- c("z1", "z2", "z3", "z4")          # canonical ids drift
  expect_error(validate_gp_fit(fit2), "ids")

  fit3 <- make_fit("none")
  fit3$posterior$id <- c("p1", "p2", "p3", "p4") # posterior id drift
  expect_error(validate_gp_fit(fit3), "posterior id")
})

# ===========================================================================
# CROSS-SLOT CHECK (b): selected_grade != path member at selected_lambda
# ===========================================================================

test_that("(b) selected_grade not equal to path member is rejected", {
  fit <- make_fit("none")
  # change the selected_grade assignment so it no longer matches the
  # path member at selected_lambda
  fit$selected_grade$assignment$grade <- c(3L, 2L, 1L, 1L)
  expect_error(validate_gp_fit(fit), "path member")
})

test_that("(b) selected_lambda not on the grid is rejected", {
  fit <- make_fit("none")
  fit$grade_path$selection$selected_lambda <- 0.37   # not a grid point
  expect_error(validate_gp_fit(fit), "grid")
})

# ===========================================================================
# CROSS-SLOT CHECK (c): report <-> selected grade disagree
# ===========================================================================

test_that("(c) report_card grades disagreeing with selected_grade is rejected", {
  fit <- make_fit("none")
  fit$report_card$grades <- c(2L, 2L, 1L, 1L)        # != selected grades
  fit$report_card$table$grade <- c(2L, 2L, 1L, 1L)   # keep table self-consistent
  expect_error(validate_gp_fit(fit), "report_card grades")
})

test_that("(c) report_card selected_lambda disagreeing with path is rejected", {
  fit <- make_fit("none")
  fit$report_card$selected_lambda <- 0.50            # path selection is 0.25
  expect_error(validate_gp_fit(fit), "selected_lambda")
})

# ===========================================================================
# CROSS-SLOT CHECK (c'): report <-> posterior columns disagree
# ===========================================================================

test_that("(c') report_card posterior_mean disagreeing with gp_posterior fails", {
  fit <- make_fit("none")
  fit$report_card$table$posterior_mean <-
    fit$report_card$table$posterior_mean + 10        # diverge from posterior
  # keep the lower<=mean<=upper band valid so we isolate the cross-slot check
  fit$report_card$table$upper <- fit$report_card$table$posterior_mean + 0.5
  expect_error(validate_gp_fit(fit), "posterior_mean must agree")
})

test_that("(c') report_card lower/upper disagreeing with gp_posterior fails", {
  fit <- make_fit("none")
  fit$report_card$table$lower <- fit$report_card$table$lower - 5
  expect_error(validate_gp_fit(fit), "lower must agree")
})

# ===========================================================================
# CROSS-SLOT CHECK (d): precision_fit presence <-> precision_rule mismatch
# ===========================================================================

test_that("(d) precision_fit non-NULL under rule 'none' is rejected", {
  fit <- make_fit("none")
  fit$precision_fit <- make_precision_fit()          # should be NULL
  expect_error(validate_gp_fit(fit), "precision_fit must be NULL")
})

test_that("(d) precision_fit NULL under rule 'krw_gmm' is rejected", {
  fit <- make_fit("krw_gmm")
  # keep the slot PRESENT but NULL (`$<- NULL` would DELETE it and trip the
  # field-order check first); `["x"] <- list(NULL)` sets it to NULL in place.
  fit["precision_fit"] <- list(NULL)                  # should be non-NULL
  expect_error(validate_gp_fit(fit), "must be non-NULL")
})

# ===========================================================================
# CROSS-SLOT CHECK (e): krw_gmm missing original_s / parameter keys
# ===========================================================================

test_that("(e) krw_gmm with missing original_s is rejected", {
  fit <- make_fit("krw_gmm")
  # drop original_s from BOTH possible homes
  fit$estimates$original_s <- NULL
  if (!is.null(fit$estimates$metadata)) fit$estimates$metadata$original_s <- NULL
  expect_error(validate_gp_fit(fit), "original_s")
})

test_that("(e) krw_gmm with wrong-length original_s is rejected", {
  fit <- make_fit("krw_gmm")
  os <- rep(0.3, GP_J - 1L)                           # length J-1
  if (!is.null(fit$estimates$metadata)) {
    fit$estimates$metadata$original_s <- os
  } else {
    fit$estimates$original_s <- os
  }
  expect_error(validate_gp_fit(fit), "length J|same length")
})

test_that("(e) krw_gmm with precision_fit$parameters missing keys is rejected", {
  fit <- make_fit("krw_gmm")
  fit$precision_fit$parameters$mu <- NULL            # drop a required key
  expect_error(validate_gp_fit(fit), "mu|model_form|beta|missing")
})
