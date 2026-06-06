# Shared frozen-input + digest computation for the ebrecipe primitive checksum
# suite (hardens GP-DEC-07-A).
#
# This ONE definition is used by BOTH:
#   * tests/testthat/test-ebrecipe-primitive-checksum.R  (recompute + compare),
#   * data-raw/make-checksums.R                          (regenerate the golden .rds),
# so the shipped golden digests and the test can never drift apart.
#
# It runs each reused primitive THROUGH the `.gp_eb_*` seam wrappers
# (R/seam-ebrecipe.R) on a tiny, fully deterministic frozen input, and returns a
# named character vector of 10-decimal-rounded digests — one per primitive
# output. A behavioral drift in an upstream ebrecipe primitive (same signature,
# changed numerics) changes a digest and trips the test loudly, naming the
# primitive. Requires ebrecipe + digest (both Suggests at test time); callers
# guard with skip_if_not_installed().

# 10-decimal digest of a numeric object (matrix / vector / scalar).
# Digests `round(x, 10L)` DIRECTLY (not `as.numeric(x)`): `round` preserves the
# `dim` attribute for matrices and `digest` serializes it, so a pure
# transpose/reshape drift (same values, swapped dims) in an upstream primitive
# is caught too — the flatten form would miss it (QA-2.2 hardening).
.gp_checksum_hash10 <- function(x) {
  digest::digest(round(x, 10L))
}

# Build the frozen inputs and return the named digest vector. Internal seam
# wrappers (.gp_eb_*) are visible unqualified under both testthat (package
# namespace) and devtools::load_all().
gp_compute_primitive_checksums <- function() {
  if (!requireNamespace("ebrecipe", quietly = TRUE) ||
      !requireNamespace("digest", quietly = TRUE)) {
    stop("gp_compute_primitive_checksums() needs ebrecipe and digest.")
  }
  h <- .gp_checksum_hash10

  ## ---- frozen, deterministic inputs (no RNG) ------------------------------
  support   <- seq(-3, 3, length.out = 64L)
  basis     <- .gp_eb_spline_basis(support, n_knots = 5L)        # 64 x df matrix
  df        <- ncol(basis)
  alpha     <- rep(0, df)
  alpha_free <- rep(0, df - 1L)
  target_mean <- 0

  theta_hat <- c(-1.2, -0.5, 0.1, 0.7, 1.3)
  s         <- c(0.4, 0.5, 0.3, 0.6, 0.45)

  ## ---- per-primitive outputs through the seam wrappers --------------------
  # softmax returns list(g, log_g, Z, log_Z); we pin g + log_g (the density and
  # its log). Z/log_Z (normalizing constants) are intentionally omitted: g
  # already encodes the normalization, so a drift there shows up in g.
  sm   <- .gp_eb_softmax_density(basis, alpha)
  log_P <- .gp_eb_normal_mixture_matrix(theta_hat, s, support)  # M x J (log)
  pen  <- .gp_eb_penalized_loglik(alpha, basis, log_P, penalty = 0)
  rlse <- .gp_eb_row_log_sum_exp(log_P)
  sat  <- .gp_eb_solve_alpha_T(alpha_free, basis, support, target_mean)
  fa   <- .gp_eb_full_alpha(alpha_free, basis, support, target_mean)
  # one-level pushforward (invariant 8): single-support kernel. Pin BOTH branches
  # KRW uses: multiplicative (race, characteristic = "white") and additive
  # (gender, characteristic = "male") — the additive branch is the under-tested
  # one (cf. N10), so both are checksummed.
  # `.eb_pushforward_theta` returns list(support, g, density): `$g` is the
  # pass-through INPUT mixing density (echoes the argument — inert to pin), the
  # computed output is `$density` (the pushed one-level g_theta) — so we hash
  # `$density`. psi MUST be non-trivial (0.1, 0.5; the Ch7 checksum example):
  # with psi=(0,0) the multiplicative and additive transforms degenerate to the
  # SAME density (== input g), so the two branches would be indistinguishable.
  pf_mult <- .gp_eb_pushforward_theta(support, sm$g, s = 0.5,
                                      psi_1 = 0.1, psi_2 = 0.5,
                                      characteristic = "white")
  pf_add  <- .gp_eb_pushforward_theta(support, sm$g, s = 0.5,
                                      psi_1 = 0.1, psi_2 = 0.5,
                                      characteristic = "male")
  # posterior weights (the W primitive): object-based call (eb_estimates +
  # eb_prior), prior via eb_deconvolve. In ebrecipe 0.5.0 `.eb_posterior_weights()`
  # returns the J x M row-stochastic weight matrix DIRECTLY (a bare matrix, NOT a
  # list) — this exact matrix is what the GP-W-EXACT golden master compares
  # gradepath's native W against (max|W_gp - this| == 0).
  est  <- .gp_eb_input(theta_hat = theta_hat, s = s)
  prio <- ebrecipe::eb_deconvolve(est)
  Wmat <- .gp_eb_posterior_weights(est, prio)

  c(
    spline_basis              = h(basis),
    softmax_density_g          = h(sm$g),
    softmax_density_log_g      = h(sm$log_g),
    normal_mixture_matrix      = h(log_P),
    penalized_loglik           = h(pen),
    row_log_sum_exp            = h(rlse),
    solve_alpha_T              = h(sat),
    full_alpha                 = h(fa),
    pushforward_theta_mult     = h(pf_mult$density),
    pushforward_theta_add      = h(pf_add$density),
    posterior_weights          = h(Wmat)
  )
}
