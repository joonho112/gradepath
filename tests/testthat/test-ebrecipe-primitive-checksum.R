# Primitive-compatibility checksum suite (hardens GP-DEC-07-A).
#
# A LOUD failure here => an upstream ebrecipe release drifted a reused primitive
# (same signature, changed numerics) — which the version pin and the boundary
# shape-assertion both miss. The remediation is the GP-DEC-07-A escape hatch:
# re-pin the last known-good ebrecipe, or vendor the affected primitive into R/,
# after auditing whether the drift is a fix to adopt (Ch7 sec-ch07-checksum).
#
# The golden digests live in inst/extdata/checksums/ebrecipe-primitives.rds,
# frozen by data-raw/make-checksums.R. The frozen-input + digest computation is
# shared via tests/testthat/helper-checksum.R (gp_compute_primitive_checksums()),
# so this test and the golden .rds can never drift apart.

test_that("reused ebrecipe primitives are numerically unchanged (checksum)", {
  skip_if_not_installed("ebrecipe")
  skip_if_not_installed("digest")

  golden_path <- system.file("extdata", "checksums", "ebrecipe-primitives.rds",
                             package = "gradepath")
  skip_if(golden_path == "" || !file.exists(golden_path),
          "golden checksum .rds not installed")

  golden <- readRDS(golden_path)
  current <- gp_compute_primitive_checksums()

  # Surface the pinned-vs-current ebrecipe version so a tripped check routes the
  # maintainer straight to the GP-DEC-07-A re-pin-or-vendor decision.
  pinned_ver  <- attr(golden, "build_metadata")$ebrecipe_version
  current_ver <- as.character(utils::packageVersion("ebrecipe"))
  ver_note <- paste0(" [pinned against ebrecipe ", pinned_ver,
                     "; currently resolving ", current_ver, "]")

  # Same primitive set, same order.
  expect_identical(names(current), names(golden))

  # Per-primitive comparison so a drift NAMES the offending primitive.
  for (nm in names(golden)) {
    expect_identical(
      current[[nm]], golden[[nm]],
      info = paste0("ebrecipe primitive '", nm,
                    "' drifted vs the pinned golden digest", ver_note,
                    " — see Ch7 sec-ch07-checksum (re-pin or vendor).")
    )
  }
})

test_that("the checksum suite covers every reused .gp_eb_* primitive", {
  skip_if_not_installed("ebrecipe")
  skip_if_not_installed("digest")
  # Guard against a primitive being added to the seam but forgotten in the
  # checksum (coverage drift). The seam wraps 9 numeric primitives + the
  # eb_input CALL (eb_input is exercised by test-seam-primitives, not hashed).
  current <- gp_compute_primitive_checksums()
  expect_setequal(
    names(current),
    c("spline_basis", "softmax_density_g", "softmax_density_log_g",
      "normal_mixture_matrix", "penalized_loglik", "row_log_sum_exp",
      "solve_alpha_T", "full_alpha", "pushforward_theta_mult",
      "pushforward_theta_add", "posterior_weights")
  )
})
