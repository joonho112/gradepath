# Tests for gp_pairwise.
#
# Covers: a valid 3x3 gp_pairwise validates; every guard fires --
#   broken antisymmetry, diagonal != 0.5, exact off-diagonal 0, wrong slot
#   order, frozen cleanup-flag tampering, out-of-[0,1], dimname mismatch,
#   non-square, power != 0 -- plus the standalone W-orientation guard
#   (.gp_assert_w_row_stochastic): a transposed (M x J) W and a column-
#   normalized W are both rejected (frozen invariant 5).

# ---------------------------------------------------------------------------
# Helpers: build a valid cleaned 3x3 pairwise matrix + a valid gp_pairwise.
# ---------------------------------------------------------------------------

# A cleaned, antisymmetric, 0.5-diagonal, zero-floored 3x3 matrix on ids.
gpB_pairwise_matrix <- function(ids = c("a", "b", "c")) {
  # Strict upper triangle of raw probabilities, then mirror by 1 - x.
  m <- matrix(0.5, 3L, 3L, dimnames = list(ids, ids))
  m[1, 2] <- 0.70; m[2, 1] <- 0.30
  m[1, 3] <- 0.90; m[3, 1] <- 0.10
  m[2, 3] <- 0.55; m[3, 2] <- 0.45
  # zero-floor / cap to be safe (no exact zeros are present here, but mirror the
  # cleanup the real builder applies).
  off <- row(m) != col(m)
  m[off] <- pmax(pmin(m[off], 1), 1e-7)
  m
}

gpB_pairwise_cleanup <- function() {
  list(antisymmetry = TRUE, diagonal = 0.5, zero_floor = 1e-7)
}

gpB_pairwise_source <- function() {
  list(stage = "posterior", rule = "outer_product",
       assumption = "one_level_independence")
}

gpB_valid_pairwise <- function(ids = c("a", "b", "c")) {
  new_gp_pairwise(
    ids = ids,
    matrix = gpB_pairwise_matrix(ids),
    power = 0L,
    cleanup = gpB_pairwise_cleanup(),
    source = gpB_pairwise_source(),
    control = gp_control(),
    schema_version = .gradepath_schema_version,
    provenance = .gradepath_new_provenance(producer = "test"),
    warnings = character()
  )
}

# ---------------------------------------------------------------------------
# Valid object
# ---------------------------------------------------------------------------

test_that("a valid 3x3 gp_pairwise constructs and validates", {
  pw <- gpB_valid_pairwise()
  expect_s3_class(pw, "gp_pairwise")
  expect_identical(names(pw), .gp_pairwise_fields)
  expect_identical(dim(pw$matrix), c(3L, 3L))
  expect_true(all(abs(diag(pw$matrix) - 0.5) < 1e-12))
  expect_identical(pw$power, 0L)
  # idempotent: re-validating returns an equivalent object.
  expect_no_error(validate_gp_pairwise(pw))
})

# ---------------------------------------------------------------------------
# Guard: broken antisymmetry  (pi_ij + pi_ji != 1)
# ---------------------------------------------------------------------------

test_that("validate_gp_pairwise rejects broken antisymmetry", {
  pw <- gpB_valid_pairwise()
  # Perturb a single off-diagonal partner well beyond the ~1e-7 tolerance.
  pw$matrix[2, 1] <- 0.40   # now [1,2] + [2,1] = 0.70 + 0.40 = 1.10
  expect_error(
    validate_gp_pairwise(pw),
    regexp = "antisymmetry",
    class = "gradepath_error"
  )
})

# ---------------------------------------------------------------------------
# Guard: diagonal != 0.5  (invariant 6)
# ---------------------------------------------------------------------------

test_that("validate_gp_pairwise rejects a diagonal != 0.5 in the matrix", {
  pw <- gpB_valid_pairwise()
  pw$matrix[2, 2] <- 0.5001   # beyond 1e-8
  expect_error(
    validate_gp_pairwise(pw),
    regexp = "diagonal",
    class = "gradepath_error"
  )
})

test_that("validate_gp_pairwise rejects a cleanup$diagonal flag != 0.5 (frozen)", {
  pw <- gpB_valid_pairwise()
  pw$cleanup$diagonal <- 0.6
  expect_error(
    validate_gp_pairwise(pw),
    regexp = "diagonal.*0\\.5|0\\.5",
    class = "gradepath_error"
  )
})

# ---------------------------------------------------------------------------
# Guard: exact off-diagonal zero  (zero-floor)
# ---------------------------------------------------------------------------

test_that("validate_gp_pairwise rejects an exact off-diagonal zero", {
  pw <- gpB_valid_pairwise()
  # Inject an exact zero and keep antisymmetry consistent so THIS guard is what
  # fires (partner -> 1, which is allowed by the one-sided floor).
  pw$matrix[1, 3] <- 0
  pw$matrix[3, 1] <- 1
  expect_error(
    validate_gp_pairwise(pw),
    regexp = "zero",
    class = "gradepath_error"
  )
})

test_that("validate_gp_pairwise rejects a frozen zero_floor flag != 1e-7", {
  pw <- gpB_valid_pairwise()
  pw$cleanup$zero_floor <- 1e-6
  expect_error(
    validate_gp_pairwise(pw),
    regexp = "zero_floor|1e-7",
    class = "gradepath_error"
  )
})

# ---------------------------------------------------------------------------
# Guard: exact ones ARE allowed (one-sided cleanup -- Ch4 §cleanup)
# ---------------------------------------------------------------------------

test_that("validate_gp_pairwise allows an exact off-diagonal one", {
  pw <- gpB_valid_pairwise()
  pw$matrix[1, 3] <- 1          # exact one is a legitimate posterior verdict
  pw$matrix[3, 1] <- 1e-7       # partner floored, not exact zero
  expect_no_error(validate_gp_pairwise(pw))
})

# ---------------------------------------------------------------------------
# Guard: wrong slot order  (.gradepath_validate_named_fields)
# ---------------------------------------------------------------------------

test_that("validate_gp_pairwise rejects wrong slot order", {
  pw <- gpB_valid_pairwise()
  # Swap two slots: right field SET, wrong ORDER.
  reordered <- pw[c("matrix", "ids", "power", "cleanup", "source",
                    "control", "schema_version", "provenance", "warnings")]
  class(reordered) <- c("gp_pairwise", "list")
  expect_error(
    validate_gp_pairwise(reordered),
    regexp = "canonical order",
    class = "gradepath_error"
  )
})

test_that("validate_gp_pairwise rejects a missing slot", {
  pw <- gpB_valid_pairwise()
  dropped <- pw[setdiff(.gp_pairwise_fields, "power")]
  class(dropped) <- c("gp_pairwise", "list")
  expect_error(validate_gp_pairwise(dropped), class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# Guard: out-of-[0,1], non-square, dimname mismatch, power != 0
# ---------------------------------------------------------------------------

test_that("validate_gp_pairwise rejects entries outside [0, 1]", {
  pw <- gpB_valid_pairwise()
  pw$matrix[1, 2] <- 1.5
  pw$matrix[2, 1] <- -0.5
  expect_error(validate_gp_pairwise(pw), regexp = "\\[0, 1\\]",
               class = "gradepath_error")
})

test_that("validate_gp_pairwise rejects a non-square matrix", {
  pw <- gpB_valid_pairwise()
  pw$matrix <- matrix(0.5, nrow = 3L, ncol = 2L)
  expect_error(validate_gp_pairwise(pw), regexp = "square",
               class = "gradepath_error")
})

test_that("validate_gp_pairwise rejects dimnames inconsistent with ids", {
  pw <- gpB_valid_pairwise()
  rownames(pw$matrix) <- c("x", "y", "z")   # != ids
  expect_error(validate_gp_pairwise(pw), regexp = "row and column names",
               class = "gradepath_error")
})

test_that("validate_gp_pairwise rejects power != 0", {
  pw <- gpB_valid_pairwise()
  pw$power <- 2L
  expect_error(validate_gp_pairwise(pw), regexp = "power",
               class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# Guard: control must be a gp_control
# ---------------------------------------------------------------------------

test_that("validate_gp_pairwise rejects a non-gp_control control slot", {
  pw <- gpB_valid_pairwise()
  pw$control <- list(lambda_grid = c(0, 0.25, 1))  # bare list, not gp_control
  expect_error(validate_gp_pairwise(pw), class = "gradepath_error")
})

# ---------------------------------------------------------------------------
# W-orientation guard (frozen invariant 5 / Chapter 8 W1-W2).
# Not part of the gp_pairwise schema, but the named home of the orientation
# contract at the decision boundary -- test it directly.
# ---------------------------------------------------------------------------

test_that(".gp_assert_w_row_stochastic accepts a valid J x M row-stochastic W", {
  J <- 4L; M <- 5L
  W <- matrix(runif(J * M), J, M)
  W <- W / rowSums(W)                       # row-stochastic
  expect_no_error(.gp_assert_w_row_stochastic(W, J, M))
  expect_true(all(abs(rowSums(W) - 1) < 1e-12))
})

test_that(".gp_assert_w_row_stochastic rejects the v1 M x J transpose", {
  J <- 4L; M <- 5L
  W <- matrix(runif(J * M), J, M)
  W <- W / rowSums(W)
  Wt <- t(W)                                # M x J -- the v1 orientation
  expect_error(
    .gp_assert_w_row_stochastic(Wt, J, M),
    regexp = "J x M|transpose",
    class = "gradepath_error"
  )
})

test_that(".gp_assert_w_row_stochastic rejects a column-normalized (non-row-stochastic) W", {
  J <- 4L; M <- 5L
  W <- matrix(runif(J * M), J, M)
  W <- sweep(W, 2L, colSums(W), "/")        # columns sum to 1, rows do not
  # Right shape (J x M) but wrong normalization -> W2 fails.
  expect_error(
    .gp_assert_w_row_stochastic(W, J, M),
    regexp = "row-stochastic|sum to 1",
    class = "gradepath_error"
  )
})

test_that(".gp_assert_w_row_stochastic rejects negative / non-finite W", {
  J <- 3L; M <- 3L
  W <- matrix(1 / M, J, M)
  W[1, 1] <- -W[1, 1]
  expect_error(.gp_assert_w_row_stochastic(W, J, M), class = "gradepath_error")
})
