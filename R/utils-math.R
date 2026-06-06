# Lean numeric helpers for gradepath v2.
#
# These are the small, broadly-reused numeric primitives the estimation core,
# the W-seam, the posterior, and the pairwise stage all reach for. They are kept
# here (not in utils-validate.R) so the validation file stays purely structural.
# The stable ROW log-sum-exp lives in utils-validate.R as `.gp_row_logsumexp()`
# (it is wanted by the validation/seam layer too); this file provides the SCALAR
# log-sum-exp and the normalization helpers that build on it.
#
# Nothing here is exported (Chapter 13 exports only the four verbs + monolith +
# gp_control + accessors). All carry the `.gp_*` internal prefix.

# Stable scalar log-sum-exp of a numeric vector x:
#   log(sum_i exp(x_i)) = m + log(sum_i exp(x_i - m)),  m = max_i x_i.
# Returns -Inf for an all -Inf vector (clean "zero mass"); rejects NA. This is
# the 1-D companion to `.gp_row_logsumexp()`; the deconvolution penalty grid and
# the GMM objective use the scalar form.
.gp_logsumexp <- function(x) {
  if (!is.numeric(x) || length(x) < 1L) {
    .gradepath_abort("`x` must be a non-empty numeric vector.")
  }
  if (anyNA(x)) {
    .gradepath_abort("`x` must not contain NA.")
  }

  m <- max(x)
  if (!is.finite(m)) {
    # all -Inf (or, defensively, +Inf): max already encodes the answer.
    return(m)
  }

  m + log(sum(exp(x - m)))
}

# Softmax of a numeric vector x: exp(x - logsumexp(x)), a probability vector
# summing to 1. Stable by construction. Used to turn log-weights into mixing
# weights without materializing an overflowing exp().
.gp_softmax <- function(x) {
  if (!is.numeric(x) || length(x) < 1L) {
    .gradepath_abort("`x` must be a non-empty numeric vector.")
  }
  if (anyNA(x)) {
    .gradepath_abort("`x` must not contain NA.")
  }

  lse <- .gp_logsumexp(x)
  if (!is.finite(lse)) {
    # No mass anywhere: fall back to a uniform vector rather than NaN.
    return(rep(1 / length(x), length(x)))
  }

  exp(x - lse)
}

# Safe normalization of a non-negative numeric weight vector to sum 1.
#
# Divides by the total mass; if the total is zero or non-finite (the degenerate
# row the W-seam guards against), it falls back to a
# uniform distribution and records nothing here -- the caller is responsible for
# the WEIGHT_DEGENERATE status (Chapter 13). Negative entries are rejected:
# weights are masses, and a negative mass is a bug upstream, not something to
# silently clamp.
#
# `w` must be a non-empty numeric vector with no NA and no negative entries.
# Returns a numeric vector the same length as `w` summing to 1 (to floating
# tolerance).
.gp_safe_normalize <- function(w) {
  if (!is.numeric(w) || length(w) < 1L) {
    .gradepath_abort("`w` must be a non-empty numeric vector.")
  }
  if (anyNA(w)) {
    .gradepath_abort("`w` must not contain NA.")
  }
  if (any(w < 0)) {
    .gradepath_abort("`w` must be non-negative (weights are masses).")
  }

  total <- sum(w)
  if (!is.finite(total) || total <= 0) {
    return(rep(1 / length(w), length(w)))
  }

  w / total
}

# Clamp a numeric vector into [lo, hi]. A tiny convenience used to enforce the
# pairwise zero-floor / probability range without scattering pmin/pmax pairs.
# `lo <= hi` required.
.gp_clamp <- function(x, lo, hi) {
  if (!is.numeric(x)) {
    .gradepath_abort("`x` must be numeric.")
  }
  lo <- .gradepath_validate_scalar_numeric(lo, "lo", finite = TRUE)
  hi <- .gradepath_validate_scalar_numeric(hi, "hi", finite = TRUE)
  if (lo > hi) {
    .gradepath_abort("`lo` must be <= `hi`.")
  }

  pmin(pmax(x, lo), hi)
}
