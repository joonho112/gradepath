#!/usr/bin/env Rscript
# Check or rebuild the M2 acceptance scorecard artifact.

.usage <- function() {
  paste(
    "Usage: Rscript inst/scripts/diagnose-m2-acceptance.R [--check|--write]",
    "",
    "Modes:",
    "  --check  Rebuild in memory and compare to inst/extdata/acceptance/m2-scorecard.csv (default).",
    "  --write  Overwrite inst/extdata/acceptance/m2-scorecard.csv with the rebuilt scorecard.",
    sep = "\n"
  )
}

.parse_mode <- function(args) {
  if (length(args) == 0L) return("check")
  if (any(args %in% c("-h", "--help"))) {
    cat(.usage(), "\n")
    quit(save = "no", status = 0L)
  }
  allowed <- c("--check", "--write")
  unknown <- setdiff(args, allowed)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown argument(s): %s\n%s",
                 paste(unknown, collapse = ", "), .usage()), call. = FALSE)
  }
  if (all(allowed %in% args)) {
    stop(sprintf("Use only one mode.\n%s", .usage()), call. = FALSE)
  }
  if ("--write" %in% args) "write" else "check"
}

.package_root <- function() {
  root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "DESCRIPTION"))) return(root)
  candidates <- c(".", "..", "../..")
  hits <- vapply(candidates, function(path) {
    file.exists(file.path(path, "DESCRIPTION")) &&
      dir.exists(file.path(path, "inst", "extdata", "acceptance"))
  }, logical(1))
  if (!any(hits)) {
    stop("Cannot find package root from current working directory.",
         call. = FALSE)
  }
  normalizePath(candidates[which(hits)[1L]], winslash = "/", mustWork = TRUE)
}

.load_gradepath <- function() {
  if (requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(devtools::load_all(".", quiet = TRUE))
  } else {
    library(gradepath)
  }
}

.write_scorecard_csv <- function(table, path) {
  utils::write.csv(table, path, row.names = FALSE, na = "")
}

.same_file_bytes <- function(a, b) {
  if (!file.exists(a) || !file.exists(b)) return(FALSE)
  identical(readBin(a, "raw", n = file.info(a)$size),
            readBin(b, "raw", n = file.info(b)$size))
}

.drift_summary <- function(expected_path, generated_path) {
  expected <- utils::read.csv(expected_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  generated <- utils::read.csv(generated_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
  reasons <- character()
  if (!identical(names(expected), names(generated))) {
    reasons <- c(reasons, "column names differ")
  }
  if (!identical(dim(expected), dim(generated))) {
    reasons <- c(reasons, sprintf("dimensions differ: expected %s, generated %s",
                                  paste(dim(expected), collapse = "x"),
                                  paste(dim(generated), collapse = "x")))
  }
  if (identical(dim(expected), dim(generated)) &&
      identical(names(expected), names(generated))) {
    mismatch <- expected != generated
    mismatch[is.na(mismatch)] <- FALSE
    if (any(mismatch)) {
      idx <- which(mismatch, arr.ind = TRUE)
      preview <- idx[seq_len(min(5L, nrow(idx))), , drop = FALSE]
      cells <- sprintf("row %d col %s", preview[, 1L],
                       names(expected)[preview[, 2L]])
      reasons <- c(reasons, sprintf("value mismatch at %s",
                                    paste(cells, collapse = "; ")))
    }
  }
  if (length(reasons) == 0L) {
    "byte-level difference only"
  } else {
    paste(reasons, collapse = "; ")
  }
}

mode <- .parse_mode(commandArgs(trailingOnly = TRUE))
root <- .package_root()
setwd(root)
.load_gradepath()

out <- file.path("inst", "extdata", "acceptance", "m2-scorecard.csv")
scorecard <- gp_m2_acceptance()
tmp <- tempfile("m2-scorecard-", fileext = ".csv")
on.exit(unlink(tmp), add = TRUE)
.write_scorecard_csv(scorecard$table, tmp)

message(sprintf("M2 acceptance diagnostic mode: %s", mode))
message(sprintf("Generated scorecard rows: %d", nrow(scorecard$table)))

if (identical(mode, "write")) {
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  before_same <- .same_file_bytes(out, tmp)
  .write_scorecard_csv(scorecard$table, out)
  if (before_same) {
    message(sprintf("Wrote %s (%d rows; byte-identical to prior artifact)",
                    out, nrow(scorecard$table)))
  } else {
    message(sprintf("Wrote %s (%d rows; artifact changed)",
                    out, nrow(scorecard$table)))
  }
  quit(save = "no", status = 0L)
}

if (!file.exists(out)) {
  stop(sprintf("Missing committed scorecard: %s. Run with --write to create it.",
               out), call. = FALSE)
}

if (.same_file_bytes(out, tmp)) {
  message(sprintf("Check OK: %s is byte-identical to regenerated scorecard.", out))
  quit(save = "no", status = 0L)
}

summary <- .drift_summary(out, tmp)
stop(sprintf(
  "M2 scorecard drift detected: %s. Run with --write to refresh %s.",
  summary, out
), call. = FALSE)
