# Diagnose KRW Table F5/F6 report-card target cells.
#
# This script is intentionally inert when sourced by tests. To refresh the
# diagnostic artifact, execute from the package root with:
#
#   GRADEPATH_RUN_REPORT_CARD_TARGETS_DIAGNOSTIC=true \
#   Rscript inst/scripts/diagnose-report-card-targets.R

diagnose_report_card_targets_reference_root <- function() {
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
    if (file.exists(file.path(candidate, "tables", "tableF5.csv")) &&
        file.exists(file.path(candidate, "tables", "tableF6.csv"))) {
      return(candidate)
    }
  }
  stop("KRW Table F5/F6 reference root was not found.", call. = FALSE)
}

diagnose_report_card_targets_file_hash <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE))
  }
  NA_character_
}

diagnose_report_card_targets <- function(output_dir = file.path("inst", "extdata", "acceptance")) {
  root <- diagnose_report_card_targets_reference_root()
  f5_path <- file.path(root, "tables", "tableF5.csv")
  f6_path <- file.path(root, "tables", "tableF6.csv")
  targets <- .gp_krw_report_card_targets(f5_path, f6_path)
  targets$producer <- "krw_report_card_table_cells"
  targets$source_table_path <- ifelse(
    targets$source_table == "F5",
    f5_path,
    f6_path
  )
  targets$source_table_sha256 <- vapply(
    targets$source_table_path,
    diagnose_report_card_targets_file_hash,
    character(1)
  )
  targets$notes <- paste(
    "Source-table adapter for published KRW Online Appendix Table F5/F6",
    "cells. Native live gp_fit report-card extraction remains gated by",
    "producer status and industry pipeline availability."
  )
  targets <- targets[, c(
    "id",
    "source_table",
    "source_table_path",
    "source_table_sha256",
    "firm_pattern",
    "column",
    "producer",
    "replicated",
    "producer_status",
    "group",
    "quantity",
    "paper",
    "delta",
    "tol",
    "unit",
    "class",
    "milestone",
    "status",
    "reason",
    "notes"
  )]
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  path <- file.path(output_dir, "m1-report-card-targets.csv")
  utils::write.csv(targets, path, row.names = FALSE, na = "")
  path
}

if (identical(Sys.getenv("GRADEPATH_RUN_REPORT_CARD_TARGETS_DIAGNOSTIC"), "true")) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The `pkgload` package is required to run this diagnostic.", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE)
  message(diagnose_report_card_targets())
}
