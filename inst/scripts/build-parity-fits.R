# Build (regenerate) the bundled parity fits shipped at
#   inst/extdata/cached/fit_race_parity.rds   (grades 2/81/14, beta 0.5095)
#   inst/extdata/cached/fit_gender_parity.rds (grades 1/3/89/4, beta 1.2554)
#
# These are the proven 97-firm Kline-Rose-Walters report-card fits that the
# applied vignettes load via gp_parity_fit() to show the published headline
# grades WITHOUT re-solving. Each fit is a full 97-firm grade integer-program
# solve, so this REQUIRES Gurobi (the default backend; a few minutes total) and
# is run only when the fits need to be regenerated -- the shipped RDS are
# otherwise stable build artifacts kept under version control.
#
# The solve policy matches the acceptance-mode injection in build-cached-assets.R
# (lambda = 0.25, mip_gap = 0): the published distributions are proven optimal.
# This script is deliberately SEPARATE from build-cached-assets.R so that
# regenerating the parity fits never perturbs that file's cache-source-hash
# manifest (which would risk a spurious CACHE_STALE on the M1 scorecard).
#
# Run from the package root:
#   Rscript inst/scripts/build-parity-fits.R

if (!file.exists("DESCRIPTION")) {
  stop("Run from the gradepath package root.", call. = FALSE)
}

# Load the package: devtools in a source checkout, else the installed package.
if (requireNamespace("devtools", quietly = TRUE)) {
  suppressMessages(devtools::load_all(".", quiet = TRUE))
} else {
  library(gradepath)
}

# gp_krw_gmm_input() is internal; reach it through the namespace (avoids a literal
# ':::' so R CMD check stays clean), mirroring build-cached-assets.R.
ns <- asNamespace("gradepath")
gmm_input <- get("gp_krw_gmm_input", envir = ns)

out_dir <- file.path("inst", "extdata", "cached")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ctrl <- gp_control(
  lambda_grid    = c(0.25, 1),
  backend        = "gurobi",
  precision_rule = "krw_gmm",
  time_limit     = 120,
  mip_gap        = 0
)

expected <- list(race = c(2L, 81L, 14L), gender = c(1L, 3L, 89L, 4L))

for (dem in c("race", "gender")) {
  message("Solving ", dem, " parity fit (97-firm grade IP, Gurobi) ...")
  inp <- gmm_input(dem)
  dat <- data.frame(
    theta_hat = inp$theta_hat, s = inp$s, unit_id = inp$unit_id,
    stringsAsFactors = FALSE
  )
  fit <- krw_report_card(
    dat, demographic = dem, control = ctrl,
    lambda = 0.25, acceptance_mode = TRUE
  )
  grades <- as.integer(table(get_grades(fit)))
  if (!identical(grades, expected[[dem]])) {
    stop(sprintf("%s: got grades %s, expected %s -- refusing to ship.",
                 dem, paste(grades, collapse = "/"),
                 paste(expected[[dem]], collapse = "/")), call. = FALSE)
  }
  path <- file.path(out_dir, sprintf("fit_%s_parity.rds", dem))
  saveRDS(fit, path)
  message("  wrote ", path, "  (grades ", paste(grades, collapse = "/"),
          ", beta ", round(get_prior(fit)$metadata$beta, 4), ")")
}

message("Done. Both parity fits regenerated and grade-verified.")
