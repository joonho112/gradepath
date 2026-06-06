# Diagnose KRW R2 targets from the companion `rsquared.do` definition.
#
# This script is intentionally inert when sourced by tests. To refresh the
# diagnostic artifact, execute from the package root with:
#
#   GRADEPATH_RUN_KRW_R2_DIAGNOSTIC=true \
#   Rscript inst/scripts/diagnose-krw-r2.R

diagnose_krw_r2_reference_root <- function() {
  candidates <- c(
    "../KRW-2024-companion-public",
    file.path(
      "log",
      "029_external-review-prep",
      "materials",
      "06_reference-code",
      "krw-replication-archive"
    ),
    file.path(
      "log",
      "029_external-review-prep",
      "materials",
      "06_reference-code",
      "krw-companion-public"
    )
  )
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "code", "rsquared.do")) &&
        file.exists(file.path(candidate, "data", "theta_estimates_race.csv")) &&
        file.exists(file.path(candidate, "dump", "ranking_results_log_dif_binary_race.csv"))) {
      return(candidate)
    }
  }
  stop("KRW R2 reference root was not found.", call. = FALSE)
}

diagnose_krw_r2_file_hash <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE))
  }
  NA_character_
}

diagnose_krw_r2_read <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

diagnose_krw_r2_row <- function(root,
                                demographic,
                                model,
                                target_id,
                                quantity,
                                value = c("r2_percent", "between_sd")) {
  value <- match.arg(value)
  theta_path <- file.path(root, "data", sprintf("theta_estimates_%s.csv", demographic))
  ranking_path <- file.path(root, "dump", sprintf("ranking_results_%s_%s.csv", model, demographic))
  theta <- diagnose_krw_r2_read(theta_path)
  ranking <- diagnose_krw_r2_read(ranking_path)
  calc <- gp_krw_r2(theta, ranking, scale = "percent")
  replicated <- if (identical(value, "r2_percent")) calc$r2 else calc$between_sd
  checked <- gp_check(target_id, replicated)

  data.frame(
    diagnostic_id = paste(demographic, model, value, sep = "_"),
    demographic = demographic,
    model = model,
    target_id = target_id,
    quantity = quantity,
    value = value,
    producer = "gp_krw_r2",
    source_script = file.path(root, "code", "rsquared.do"),
    source_script_sha256 = diagnose_krw_r2_file_hash(file.path(root, "code", "rsquared.do")),
    source_theta = theta_path,
    source_ranking = ranking_path,
    grade_col = calc$grade_col,
    n = calc$n,
    grade_count = calc$grade_count,
    overall_sd = calc$overall_sd,
    between_sd = calc$between_sd,
    r2_proportion = calc$r2_proportion,
    r2_percent = 100 * calc$r2_proportion,
    replicated = replicated,
    paper = checked$paper,
    delta = checked$delta,
    tolerance = checked$tol,
    unit = checked$unit,
    class = checked$class,
    milestone = checked$milestone,
    status = checked$status,
    reason = checked$reason %gp_or% "",
    notes = "Computed from KRW rsquared.do second-moment recipe; not produced by generic gp_r2().",
    stringsAsFactors = FALSE
  )
}

diagnose_krw_r2 <- function(output_dir = file.path("inst", "extdata", "acceptance")) {
  root <- diagnose_krw_r2_reference_root()
  out <- do.call(rbind, list(
    diagnose_krw_r2_row(
      root,
      demographic = "race",
      model = "log_dif_binary",
      target_id = "race_baseline_r2",
      quantity = "race baseline grade R2",
      value = "r2_percent"
    ),
    diagnose_krw_r2_row(
      root,
      demographic = "race",
      model = "log_dif_binary",
      target_id = "race_baseline_betweengrade_sd",
      quantity = "race between-grade SD",
      value = "between_sd"
    ),
    diagnose_krw_r2_row(
      root,
      demographic = "gender",
      model = "log_dif_binary",
      target_id = "gender_baseline_r2",
      quantity = "gender baseline grade R2",
      value = "r2_percent"
    ),
    diagnose_krw_r2_row(
      root,
      demographic = "race",
      model = "industry_rfe_binary",
      target_id = "race_industry_r2",
      quantity = "race-industry grade R2",
      value = "r2_percent"
    ),
    diagnose_krw_r2_row(
      root,
      demographic = "gender",
      model = "industry_rfe_binary",
      target_id = "gender_industry_r2",
      quantity = "gender-industry grade R2",
      value = "r2_percent"
    )
  ))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  path <- file.path(output_dir, "m1-krw-r2-diagnostic.csv")
  utils::write.csv(out, path, row.names = FALSE, na = "")
  path
}

if (identical(Sys.getenv("GRADEPATH_RUN_KRW_R2_DIAGNOSTIC"), "true")) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The `pkgload` package is required to run this diagnostic.", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE)
  message(diagnose_krw_r2())
}
