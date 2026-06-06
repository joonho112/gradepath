# =============================================================================
# krw-report-card-targets.R -- Table F5/F6 source-table target adapters
# -----------------------------------------------------------------------------
# These helpers read the KRW companion's published report-card tables and turn
# selected named cells into registry-scale replicated values. They do not claim
# that a live native gp_fit generated the table; they are source-truth adapters
# for the final published F5/F6 table cells.
# =============================================================================

.gp_krw_report_card_target_specs <- function() {
  data.frame(
    id = c(
      "f5_genuineparts_theta",
      "f5_genuineparts_postmean_baseline",
      "f5_genuineparts_postmean_industry",
      "f5_autonation_theta",
      "f5_autonation_condrank_industry",
      "f5_disney_theta",
      "f5_walmart_theta",
      "f5_charterspectrum_theta",
      "f5_charterspectrum_condrank_baseline",
      "f6_buildersfirstsource_theta",
      "f6_buildersfirstsource_postmean_baseline",
      "f6_buildersfirstsource_postmean_industry",
      "f6_lkq_theta",
      "f6_nationwide_theta",
      "f6_ascena_theta",
      "f6_ascena_condrank_baseline"
    ),
    table = c(rep("F5", 9L), rep("F6", 7L)),
    firm_pattern = c(
      "^Genuine Parts",
      "^Genuine Parts",
      "^Genuine Parts",
      "^AutoNation",
      "^AutoNation",
      "^Disney",
      "^Walmart",
      "^Charter",
      "^Charter",
      "^Builders FirstSource",
      "^Builders FirstSource",
      "^Builders FirstSource",
      "^LKQ",
      "^Nationwide",
      "^Ascena",
      "^Ascena"
    ),
    column = c(
      "log_dif",
      "post_mean_beta",
      "ind_post_mean_beta",
      "log_dif",
      "grade1_ind",
      "log_dif",
      "log_dif",
      "log_dif",
      "grade1",
      "log_dif",
      "post_mean_beta",
      "ind_post_mean_beta",
      "log_dif",
      "log_dif",
      "log_dif",
      "grade1"
    ),
    stringsAsFactors = FALSE
  )
}

.gp_krw_report_card_table <- function(x, name) {
  if (is.data.frame(x)) {
    return(x)
  }
  path <- .gradepath_validate_scalar_character(x, name)
  if (!file.exists(path)) {
    .gradepath_abort("File not found: %s.", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.gp_krw_report_card_cell <- function(table, firm_pattern, column, source_name) {
  required <- c("firm_name", column)
  miss <- setdiff(required, names(table))
  if (length(miss) > 0L) {
    .gradepath_abort("`%s` is missing column(s): %s.", source_name, paste(miss, collapse = ", "))
  }
  hit <- which(grepl(firm_pattern, table$firm_name, ignore.case = TRUE))
  if (length(hit) < 1L) {
    .gradepath_abort("No firm in `%s` matches `%s`.", source_name, firm_pattern)
  }
  if (length(hit) > 1L) {
    .gradepath_abort(
      "Expected exactly one firm in `%s` to match `%s`; found %d.",
      source_name,
      firm_pattern,
      length(hit)
    )
  }
  suppressWarnings(as.numeric(table[[column]][hit[[1L]]]))
}

.gp_krw_report_card_targets <- function(table_f5,
                                        table_f6,
                                        targets = NULL,
                                        producer_status = "OK",
                                        registry = NULL) {
  f5 <- .gp_krw_report_card_table(table_f5, "table_f5")
  f6 <- .gp_krw_report_card_table(table_f6, "table_f6")
  producer_status <- .gp_status_normalize(producer_status)
  specs <- .gp_krw_report_card_target_specs()
  if (!is.null(targets)) {
    targets <- as.character(targets)
    missing_specs <- setdiff(targets, specs$id)
    if (length(missing_specs) > 0L) {
      .gradepath_abort("Unknown F5/F6 target id(s): %s.", paste(missing_specs, collapse = ", "))
    }
    specs <- specs[match(targets, specs$id), , drop = FALSE]
  }

  rows <- vector("list", nrow(specs))
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, , drop = FALSE]
    source_table <- if (identical(spec$table, "F5")) f5 else f6
    source_name <- if (identical(spec$table, "F5")) "table_f5" else "table_f6"
    replicated <- .gp_krw_report_card_cell(
      source_table,
      firm_pattern = spec$firm_pattern,
      column = spec$column,
      source_name = source_name
    )
    checked <- gp_check(
      spec$id,
      replicated = replicated,
      producer_status = producer_status,
      registry = registry
    )
    rows[[i]] <- data.frame(
      id = spec$id,
      source_table = spec$table,
      firm_pattern = spec$firm_pattern,
      column = spec$column,
      replicated = replicated,
      producer_status = producer_status,
      group = "Report card",
      checked[, c(
        "quantity",
        "paper",
        "delta",
        "tol",
        "unit",
        "class",
        "milestone",
        "status",
        "reason"
      )],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  do.call(rbind, rows)
}
