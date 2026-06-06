#' Two-level (group_fx == 1) six-moment cluster-robust GMM core
#'
#' @description
#' The two-level extension of the native beta-GMM core: a faithful port of the
#' KRW (2024) Matlab `group_fx == 1` branch (`moment_conditions.m`,
#' `get_moments_gmm.m`, `gmm_obj.m`, `estimate_lsqnonlin.m`). It fits the
#' four-parameter vector
#'
#'   delta = c(mu_raw, log_sigma_xi, beta, log_sigma_eta)
#'
#' (Matlab `sigma_psi = exp(delta(4))` is gradepath's `sigma_eta`; the archive
#' code-name `psi` IS the between-industry effect eta -- see Ch11 §2 "psi is eta")
#' against SIX moment conditions: the four one-level studentized-residual moments
#' `g_i = [r, r*s, r^2-1, (r^2-1)*s]` PLUS two between-industry (group) moments
#' built from the industry means `v_bar`, `s_bar` and the model-implied
#' between-industry mean/variance `mu_vbar`, `V_vbar`. The covariance of the
#' moment conditions is **cluster-robust** (clusters = industries):
#'   `V = (1/N) * sum_{i,j : d(i)==d(j)} g_i g_j'`.
#'
#' This file is a NEW source added; it does NOT edit
#' `R/estimation-core.R` (which is M1-cache-watched and implements only
#' `group_fx == 0`). It REUSES the existing internal helpers `gp_jacobian()`,
#' `.gradepath_abort()`, and `.gradepath_new_provenance()`, and MIRRORS the
#' structure of `gp_estimation_core()` (two-step efficient GMM, Omega1-vs-Omega2
#' split for J-stat-vs-sandwich, single factor of N). All functions are internal.
#'
#' @section PARITY (the published two-level numbers):
#' On KRW's actual GMM input (loaded via `gp_krw_gmm_input()`, the headerless
#' `[industry_d, theta_hat, s]` matrix, 97 firms / 19 SIC industries), this
#' two-step cluster-robust GMM reproduces the published two-level moments:
#' within-share ~ 0.366 (race) / 0.562 (gender), with sigma_eta ~ 0.53 and
#' sigma_xi ~ 0.11 (race), and a small chi-square(2) overidentification J-stat
#' (df = 6 moments - 4 params = 2). The one-level reduction is structurally
#' guaranteed: at `group_fx == 0` the two group moments and the cluster off-block
#' terms vanish and the construction collapses to the four-moment one-level core.
#'
#' @keywords internal
#' @noRd
NULL

# ---------------------------------------------------------------------------
# Two-level data matrix (with the REAL industry column)
# ---------------------------------------------------------------------------

#' Build the Matlab `[d, theta_hat, s]` data matrix from a two-level input
#'
#' @description
#' The two-level moment routines read the industry membership `d` (Matlab column
#' 1), unlike the one-level core which fills `d` with a constant placeholder.
#' This helper assembles `[d, theta_hat, s]` from a list that carries
#' `theta_hat`, `s`, and an integer `industry` membership vector -- exactly the
#' shape `gp_krw_gmm_input()` returns. The industry codes must be a contiguous
#' `1..J` labelling (the KRW input is already coded `1..19`); we re-key any
#' arbitrary integer/factor coding to dense `1..J` via `match()` so `dummyvar`
#' equivalents (the projection `P_D`) never index a gap.
#'
#' @param input A list with numeric `theta_hat`, `s`, and integer-codeable
#'   `industry`, all of equal length (e.g. `gp_krw_gmm_input()`).
#' @return A numeric `N x 3` matrix with columns `c("d", "theta_hat", "s")`,
#'   `d` a dense `1..J` industry index.
#' @keywords internal
#' @noRd
.gp_two_level_data_matrix <- function(input) {
  theta_hat <- as.numeric(input$theta_hat)
  s <- as.numeric(input$s)
  ind <- input$industry
  if (length(theta_hat) == 0L || length(s) == 0L) {
    .gradepath_abort(
      "`input` must carry non-empty `theta_hat` and `s` vectors.",
      class = "gradepath_validation_error"
    )
  }
  if (length(theta_hat) != length(s)) {
    .gradepath_abort(
      sprintf("`theta_hat` (%d) and `s` (%d) must be length-matched.",
              length(theta_hat), length(s)),
      class = "gradepath_validation_error"
    )
  }
  if (is.null(ind)) {
    .gradepath_abort(
      paste0("The two-level GMM (`group_fx == 1`) requires an `industry` ",
             "membership vector; `input$industry` is NULL. (Pass the output of ",
             "`gp_krw_gmm_input()`, which carries the 19-SIC industry coding.)"),
      class = "gradepath_validation_error"
    )
  }
  if (length(ind) != length(theta_hat)) {
    .gradepath_abort(
      sprintf("`industry` (%d) must be length-matched to `theta_hat` (%d).",
              length(ind), length(theta_hat)),
      class = "gradepath_validation_error"
    )
  }
  ## Re-key to dense 1..J (Matlab `dummyvar` assumes consecutive codes). The KRW
  ## input is already 1..19; match() is a no-op there but makes the helper safe.
  ## `method = "radix"` keeps the re-keying locale-independent (byte/code-point
  ## order), so the industry labelling -- and every downstream cluster-robust
  ## covariance built on it -- is reproducible across locales (CRAN/CI run under
  ## the C locale; this mirrors the M1 provenance-hash locale fix).
  d <- match(ind, sort(unique(ind), method = "radix"))
  if (max(d) < 2L) {
    .gradepath_abort(
      paste0("The two-level GMM (`group_fx == 1`) needs >= 2 industries to ",
             "identify the between-industry effect and form the cluster-robust ",
             "covariance; got ", max(d), " distinct industry code(s). (One-level ",
             "data should use the one-level core `gp_estimation_core()`.)"),
      class = "gradepath_validation_error"
    )
  }
  cbind(d = d, theta_hat = theta_hat, s = s)
}

# ---------------------------------------------------------------------------
# Industry-mean projection P_D (the `dummyvar` algebra of moment_conditions.m)
# ---------------------------------------------------------------------------

#' Per-firm industry means and sizes (the `P_D`, `n` carriers)
#'
#' @description
#' Reproduces the `dummyvar` block of `moment_conditions.m` L43-50 once per call.
#' In the Matlab, `D = dummyvar(d)`, `P_D = D*inv(D'*D)*D'` is the N x N
#' projection onto industry means (`(P_D x)_i` = mean of `x` over firm i's
#' industry), and `n = D*(sum(D,1)')` is the per-firm industry size `n_{k(i)}`.
#' The weighted carriers reduce algebraically to plain industry means:
#'
#'   w      = 1 / n_{k(i)}                              (L46-48)
#'   v_bar  = (P_D (w .* v_hat)) .* n = industry mean of v_hat   (L49)
#'   s_bar  = (P_D (w .* s))     .* n = industry mean of s       (L50)
#'
#' We return both a `proj()` closure (the literal `(P_D x)` industry-mean
#' projection, used so `V_vbar` is assembled verbatim from the Matlab
#' expressions) and the size vector `n`, the firm->industry index, and `J`.
#'
#' @param d Length-N integer (dense `1..J`) industry membership.
#' @return A list with `proj` (function: length-N x -> length-N industry means),
#'   `n` (length-N industry sizes `n_{k(i)}`), `idx` (firm->industry index),
#'   `J` (number of industries).
#' @keywords internal
#' @noRd
.gp_industry_projection <- function(d) {
  idx <- as.integer(d)
  J <- max(idx)
  sizes <- tabulate(idx, nbins = J)          # n_k per industry
  n <- sizes[idx]                            # n_{k(i)} per firm (Matlab `n`)
  ## (P_D x)_i = mean of x over firm i's industry. Vectorized industry mean.
  ## `rowsum` returns rows ordered by the sorted group key; divide by the matching
  ## industry sizes, then broadcast back to firms (`method = "radix"` keying keeps
  ## it locale-independent). This is the vectorized identity of the Matlab `P_D x`.
  grp_keys <- sort(unique(idx), method = "radix")
  proj <- function(x) {
    means_by_key <- as.numeric(rowsum(x, group = idx)) / sizes[grp_keys]
    means_by_key[match(idx, grp_keys)]
  }
  list(proj = proj, n = n, idx = idx, J = J)
}

# ---------------------------------------------------------------------------
# Two-level moment conditions (port of moment_conditions.m, group_fx == 1)
# ---------------------------------------------------------------------------

#' Six moment conditions for the two-level GMM (port of `moment_conditions.m`)
#'
#' @description
#' Forms the SIX cluster moment conditions for the two-level latent-signal model.
#' The first four are the one-level studentized-residual moments
#' `g_i = [r, r*s, r^2-1, (r^2-1)*s]` (`moment_conditions.m:38`); the last two are
#' the between-industry moments
#' `[ (v_bar - mu_vbar)^2 - V_vbar , ((v_bar - mu_vbar)^2 - V_vbar) * s_bar ]`
#' (`moment_conditions.m:60`), with `v_bar`, `s_bar` the industry means and
#' `mu_vbar`, `V_vbar` the model-implied between-industry mean/variance.
#'
#' RACE is multiplicative (`mu = exp(mu_raw)`, `v_hat = theta_hat / s^beta`,
#' `sigma^2 = sigma_xi^2 sigma_eta^2 + mu^2 sigma_eta^2 + sigma_xi^2`); GENDER is
#' additive (`v_hat = (theta_hat - mu)/s^beta`, `sigma^2 = sigma_eta^2 +
#' sigma_xi^2`). The RACE vs GENDER `V_vbar` formulas DIFFER and are ported
#' verbatim (`moment_conditions.m:52-58`):
#'
#'   RACE   (L54): V_vbar = (1+sigma_eta^2) sigma_xi^2 * (P_D(w^2) * n)
#'                          + sigma_eta^2 mu^2
#'                          + (P_D(w^2 * s^(2(1-beta))) * n)
#'                 mu_vbar = mu
#'   GENDER (L58): V_vbar = sigma_eta^2 + sigma_xi^2 / n + (P_D(s_v^2)) / n
#'                 mu_vbar = 0
#'
#' The covariance (`get_cov == TRUE`) is the CLUSTER-ROBUST form
#' `V = (1/N) sum_{i,j : d(i)==d(j)} g_i g_j'` (`moment_conditions.m:76-86`),
#' computed block-diagonally over industries (the only `i,j` pairs that survive
#' the `d(i)==d(j)` filter), which is the vectorized identity of the Matlab
#' double loop.
#'
#' @param delta Length-4 numeric `c(mu_raw, log_sigma_xi, beta, log_sigma_eta)`.
#' @param data Numeric N x 3 matrix `[d, theta_hat, s]` with a REAL industry
#'   column `d` (dense `1..J`).
#' @param characteristic `"race"` or `"gender"`.
#' @param get_cov Logical; if `TRUE`, also return the 6 x 6 cluster-robust
#'   covariance `Omega`.
#' @param proj Optional precomputed `.gp_industry_projection(data[,1])` (avoids
#'   rebuilding the projection on every Jacobian finite-difference call).
#' @return A list with `g` (length-6), `g_i` (N x 6), `Omega` (6 x 6 or `NULL`),
#'   `v_hat`, `s_v` (length-N), `r` (length-N), `v_bar`, `s_bar` (length-N).
#' @keywords internal
#' @noRd
gp_moment_conditions_2l <- function(delta, data, characteristic,
                                    get_cov = FALSE, proj = NULL) {
  ## Unpack data (Matlab column order [d, theta_hat, s]).
  d <- data[, 1L]
  theta_hat <- data[, 2L]
  s <- data[, 3L]
  N <- length(s)

  ## Unpack parameters: delta = [mu_raw, log_sigma_xi, beta, log_sigma_eta].
  mu <- delta[1L]
  sigma_xi <- exp(delta[2L])
  beta <- delta[3L]
  sigma_eta <- exp(delta[4L])              # Matlab sigma_psi = exp(delta(4))
  sigma <- sqrt(sigma_eta^2 + sigma_xi^2)

  if (identical(characteristic, "race")) {
    mu <- exp(mu)
    sigma <- sqrt((sigma_xi^2) * (sigma_eta^2) +
                    (mu^2) * (sigma_eta^2) +
                    (sigma_xi^2))
  } else if (!identical(characteristic, "gender")) {
    .gradepath_abort(
      sprintf("Unknown `characteristic`: '%s' (expected 'race' or 'gender').",
              characteristic),
      class = "gradepath_validation_error"
    )
  }

  ## Standardized residual r and the v_hat / s_v carriers.
  if (identical(characteristic, "race")) {
    v_hat <- theta_hat / (s^beta)
    s_v <- s^(1 - beta)
    r <- (v_hat - mu) / sqrt((sigma^2) + (s_v^2))
  } else {                                  # gender
    v_hat <- (theta_hat - mu) / (s^beta)
    s_v <- s^(1 - beta)
    r <- v_hat / sqrt((sigma^2) + (s_v^2))
  }

  ## Four base moments (moment_conditions.m:38).
  g_base <- cbind(r, r * s, r^2 - 1, (r^2 - 1) * s)

  ## ---- group-level moments (moment_conditions.m:41-60) --------------------
  if (is.null(proj)) proj <- .gp_industry_projection(d)
  P <- proj$proj                            # (P_D x)_i = industry mean of x
  n <- proj$n                               # n_{k(i)}
  w <- 1 / n                                # L46-48: w = 1/n_{k(i)}

  v_bar <- P(w * v_hat) * n                 # L49: industry mean of v_hat
  s_bar <- P(w * s) * n                     # L50: industry mean of s

  if (identical(characteristic, "race")) {
    mu_vbar <- mu                           # L53
    ## L54, verbatim:
    V_vbar <- (((1 + (sigma_eta^2)) * (sigma_xi^2)) * (P(w^2) * n)) +
      ((sigma_eta^2) * (mu^2)) +
      (P((w^2) * (s^(2 * (1 - beta)))) * n)
  } else {                                  # gender
    mu_vbar <- 0                            # L57
    ## L58, verbatim:
    V_vbar <- (sigma_eta^2) + ((sigma_xi^2) / n) + ((P(s_v^2)) / n)
  }

  ## Two appended group moments (moment_conditions.m:60).
  g_grp_resid <- ((v_bar - mu_vbar)^2) - V_vbar
  g_grp <- cbind(g_grp_resid, g_grp_resid * s_bar)

  g_i <- cbind(g_base, g_grp)
  colnames(g_i) <- c("r", "r_s", "r2m1", "r2m1_s", "grp", "grp_sbar")
  g <- colMeans(g_i)

  ## ---- cluster-robust covariance (moment_conditions.m:76-86) --------------
  ## V = (1/N) sum_{i,j : d(i)==d(j)} g_i g_j'. Only same-industry pairs survive,
  ## so the double loop is exactly a sum of per-industry rank-d blocks:
  ##   sum_{i,j in k} g_i g_j' = (sum_{i in k} g_i)(sum_{i in k} g_i)'.
  ## Hence V = (1/N) sum_k G_k G_k', G_k = colSums of g_i over industry k. This
  ## is the vectorized identity of the Matlab loop (verified against the literal
  ## double loop in the test file).
  Omega <- NULL
  if (isTRUE(get_cov)) {
    M <- ncol(g_i)
    Gk <- rowsum(g_i, group = proj$idx)     # J x M: per-industry moment sums
    Omega <- crossprod(Gk) / N              # (1/N) sum_k G_k G_k'  (M x M)
    dimnames(Omega) <- NULL
  }

  list(g = unname(g), g_i = unname(g_i), Omega = Omega,
       v_hat = v_hat, s_v = s_v, r = r,
       v_bar = unname(v_bar), s_bar = unname(s_bar))
}

# ---------------------------------------------------------------------------
# Two-level GMM-implied moments (port of get_moments_gmm.m, group_fx == 1)
# ---------------------------------------------------------------------------

#' GMM-implied two-level moments incl. within-share (port of `get_moments_gmm.m`)
#'
#' @description
#' Returns the natural-scale two-level moments (`get_moments_gmm.m`,
#' `group_fx == 1`). The coupling vector (`extra_moments == 0`) is:
#'   RACE   -> `[mu, sigma_xi, sigma_eta]`            (length 3)
#'   GENDER -> `[sigma_xi, sigma_eta]`                (length 2)
#' (`get_moments_gmm.m:57-68`: race prepends `mu = exp(mu_raw)`, then
#' `sigma_xi`, then `sigma_psi = sigma_eta`.)
#'
#' The `extra_moments == 1` (paper-reporting) path appends `[E_theta,
#' sd_theta, within_share]` (`get_moments_gmm.m:27-53`). The within-share is the
#' published two-level decomposition statistic:
#'   RACE   (L33): within_share = ((sigma_eta^2 + 1) sigma_xi^2) /
#'                  (sigma_eta^2 sigma_xi^2 + sigma_eta^2 mu^2 + sigma_xi^2)
#'   GENDER (L41): within_share = sigma_xi^2 / sigma^2
#'
#' @param delta Length-4 numeric `c(mu_raw, log_sigma_xi, beta, log_sigma_eta)`.
#' @param characteristic `"race"` or `"gender"`.
#' @param data Optional N x 3 matrix `[d, theta_hat, s]`; required when
#'   `extra_moments == 1` (E_theta / sd_theta read the `s` column).
#' @param extra_moments Integer `0` (default; coupling carriers) or `1`
#'   (paper-reporting moments, appends `[E_theta, sd_theta, within_share]`).
#' @return A numeric vector of GMM-implied moments.
#' @keywords internal
#' @noRd
gp_get_moments_2l <- function(delta, characteristic, data = NULL,
                              extra_moments = 0L) {
  mu <- delta[1L]
  sigma_xi <- exp(delta[2L])
  beta <- delta[3L]
  sigma_eta <- exp(delta[4L])              # Matlab sigma_psi
  sigma <- sqrt(sigma_eta^2 + sigma_xi^2)
  if (identical(characteristic, "race")) {
    mu <- exp(mu)
    sigma <- sqrt((sigma_xi^2) * (sigma_eta^2) +
                    (mu^2) * (sigma_eta^2) +
                    (sigma_xi^2))
  } else if (!identical(characteristic, "gender")) {
    .gradepath_abort(
      sprintf("Unknown `characteristic`: '%s'.", characteristic),
      class = "gradepath_validation_error"
    )
  }

  ## ---- extra (paper-reporting) moments incl. within_share -----------------
  m_extra <- numeric(0)
  if (identical(as.integer(extra_moments), 1L)) {
    if (is.null(data)) {
      .gradepath_abort(
        "`data` is required when `extra_moments == 1` (E_theta / sd_theta).",
        class = "gradepath_validation_error"
      )
    }
    s <- data[, 3L]
    if (identical(characteristic, "race")) {
      E_theta <- mu * mean(s^beta)                                  # L31
      E_theta_2 <- mean(s^(2 * beta)) * ((sigma^2) + (mu^2))        # L32
      within_share <-                                              # L33
        (((sigma_eta^2) + 1) * (sigma_xi^2)) /
          (((sigma_eta^2) * (sigma_xi^2)) +
             ((sigma_eta^2) * (mu^2)) +
             (sigma_xi^2))
    } else {                                # gender
      E_theta <- mu                                                # L39
      E_theta_2 <- (mu^2) + ((sigma^2) * mean(s^(2 * beta)))       # L40
      within_share <- (sigma_xi^2) / (sigma^2)                     # L41
    }
    ## get_moments_gmm.m:45,49: m_extra = [E_theta sd_theta]; for group_fx == 1
    ## within_share is appended -> [E_theta, sd_theta, within_share].
    m_extra <- c(E_theta, sqrt(E_theta_2 - (E_theta^2)), within_share)
  }

  ## ---- core coupling moments (get_moments_gmm.m:57-68) --------------------
  m <- numeric(0)
  if (identical(characteristic, "race")) {
    m <- mu                                 # exp(mu_raw)
  }
  m <- c(m, sigma_xi)                       # L62
  m <- c(m, sigma_eta)                      # L65 (group_fx == 1)
  m <- c(m, m_extra)                        # L68
  unname(m)
}

#' Within-share alone (convenience wrapper over `gp_get_moments_2l`)
#'
#' @description
#' Returns just the published within-share statistic at `delta`. RACE:
#' `((sigma_eta^2 + 1) sigma_xi^2) / (sigma_eta^2 sigma_xi^2 + sigma_eta^2 mu^2 +
#' sigma_xi^2)`; GENDER: `sigma_xi^2 / sigma^2` (`get_moments_gmm.m:33,41`). This
#' is closed-form in `delta` (no `data` dependence -- the `data` argument of
#' `gp_get_moments_2l` is read only for `E_theta`/`sd_theta`).
#'
#' @param delta Length-4 numeric `c(mu_raw, log_sigma_xi, beta, log_sigma_eta)`.
#' @param characteristic `"race"` or `"gender"`.
#' @return A length-1 numeric within-share in `(0, 1)`.
#' @keywords internal
#' @noRd
gp_within_share <- function(delta, characteristic) {
  mu <- delta[1L]
  sigma_xi <- exp(delta[2L])
  sigma_eta <- exp(delta[4L])
  if (identical(characteristic, "race")) {
    mu <- exp(mu)
    return(
      (((sigma_eta^2) + 1) * (sigma_xi^2)) /
        (((sigma_eta^2) * (sigma_xi^2)) +
           ((sigma_eta^2) * (mu^2)) +
           (sigma_xi^2))
    )
  }
  if (identical(characteristic, "gender")) {
    sigma <- sqrt(sigma_eta^2 + sigma_xi^2)
    return((sigma_xi^2) / (sigma^2))
  }
  .gradepath_abort(
    sprintf("Unknown `characteristic`: '%s'.", characteristic),
    class = "gradepath_validation_error"
  )
}

# ---------------------------------------------------------------------------
# Two-level start values
# ---------------------------------------------------------------------------

#' Two-level GMM start values (the Matlab reference start)
#'
#' @description
#' `estimate_lsqnonlin.m:67` uses the SAME 4-param start
#' `[-1.1598, -1.9752, 0.5193, -0.7595]` for BOTH race and gender in its
#' "one-step limited GMM" branch (when a SLURM/PBS array id is set). The
#' production `decons.m` run (no array id) instead draws 1000 `randn(4,1)`
#' multi-starts and keeps the best (`estimate_lsqnonlin.m:74-80`). gradepath
#' returns the documented fixed start as a deterministic anchor; the driver
#' `gp_two_level_gmm()` additionally runs a seeded multi-start (mirroring the
#' production path and the one-level gender hardening in `gp_estimation_core`).
#'
#' @param characteristic `"race"` or `"gender"` (the fixed start is shared, but
#'   the argument is validated for symmetry with `gp_gmm_start`).
#' @return A length-4 numeric start `c(mu_raw, log_sigma_xi, beta, log_sigma_eta)`.
#' @keywords internal
#' @noRd
gp_gmm_start_2l <- function(characteristic) {
  if (!characteristic %in% c("race", "gender")) {
    .gradepath_abort(
      sprintf("Unknown `characteristic`: '%s'.", characteristic),
      class = "gradepath_validation_error"
    )
  }
  ## estimate_lsqnonlin.m:67 (group_fx != 0 branch), shared race/gender.
  c(-1.1598, -1.9752, 0.5193, -0.7595)
}

# ---------------------------------------------------------------------------
# Minimize the two-level GMM objective J = N g' W g for one weighting
# ---------------------------------------------------------------------------

#' Minimize the two-level GMM objective (port of `gmm_obj.m` + the `lsqnonlin`
#' minimization, group_fx == 1)
#'
#' @description
#' Minimizes the scalar objective `J = N * g' W g` (`gmm_obj.m:10`, same form as
#' the one-level core) where `g` is the length-6 two-level moment vector. WEIGHT
#' convention is quadratic-form: step 1 `W = diag(6)`; step 2 `W = solve(Omega1)`
#' (the Matlab passes `chol(inv(Omega1))` to its residual form; `g' chol' chol g
#' = g' solve(Omega1) g`, the same quadratic form). The industry projection is
#' built once and reused across every objective / Jacobian evaluation.
#'
#' @param data Numeric N x 3 matrix `[d, theta_hat, s]`.
#' @param characteristic `"race"` or `"gender"`.
#' @param W The 6 x 6 GMM weight matrix (quadratic-form convention).
#' @param start Length-4 start `c(mu_raw, log_sigma_xi, beta, log_sigma_eta)`.
#' @param proj Optional precomputed `.gp_industry_projection(data[,1])`.
#' @param optimizer `"BFGS"` (default), `"Nelder-Mead"`, or `"CG"`.
#' @param reltol Relative convergence tolerance.
#' @param maxit Max iterations.
#' @return A list with `delta`, `objective`, `convergence`, `message`,
#'   `optimizer`.
#' @keywords internal
#' @noRd
gp_gmm_min_2l <- function(data, characteristic, W, start, proj = NULL,
                          optimizer = "BFGS", reltol = 1e-10, maxit = 1000L) {
  N <- nrow(data)
  if (is.null(proj)) proj <- .gp_industry_projection(data[, 1L])
  g_of <- function(delta) {
    gp_moment_conditions_2l(delta, data, characteristic,
                            get_cov = FALSE, proj = proj)$g
  }
  obj <- function(delta) {
    val <- tryCatch(
      as.numeric(N * crossprod(g_of(delta), W %*% g_of(delta))),
      error = function(e) NA_real_
    )
    ## Fence non-finite / erroring objective values so a pathological multistart
    ## draw cannot abort the optimizer mid-sweep; optim then steps away from it.
    if (!is.finite(val)) 1e10 else val
  }

  method <- if (optimizer %in% c("BFGS", "Nelder-Mead", "CG")) optimizer else "BFGS"
  fit <- stats::optim(par = start, fn = obj, method = method,
                      control = list(reltol = reltol, maxit = maxit))
  list(delta = fit$par, objective = fit$value,
       convergence = as.integer(fit$convergence),
       message = if (is.null(fit$message)) NA_character_ else fit$message,
       optimizer = method)
}

# ---------------------------------------------------------------------------
# Two-step efficient two-level GMM driver
# ---------------------------------------------------------------------------

#' Two-step efficient two-level GMM driver (port of `estimate_lsqnonlin.m`,
#' GMM section, group_fx == 1)
#'
#' @description
#' Runs the full two-step cluster-robust estimator on a two-level input and
#' returns the fit object plus the deconvolution coupling carriers and the
#' published within-share. Steps mirror `gp_estimation_core()` exactly, lifted to
#' 6 moments / 4 params / cluster-robust covariance:
#' \enumerate{
#'   \item Step 1 (identity weighting): minimize `N g' I g` from the documented
#'     Matlab start, then a seeded multi-start (the production `decons.m` path
#'     draws `randn(4,1)` x1000; we draw `n_starts` jittered starts, keep best).
#'   \item Form the CLUSTER-ROBUST `Omega1` at the step-1 optimum; set the GMM
#'     weight `W_gmm = solve(Omega1)` (Matlab `chol(inv(Omega1))`). A PD guard
#'     (`chol`) gives a clean abort on a singular/indefinite `Omega1`.
#'   \item Step 2 (efficient weighting): re-minimize `N g' W_gmm g` from delta1.
#'   \item Sandwich covariance `C = solve(G' Omega2^-1 G) / N` with ONE factor of
#'     N (invariant #7); `SE = sqrt(diag(C))`. `Omega2` is the cluster-robust
#'     covariance at delta2.
#'   \item Single-N J-statistic `J_stat = N g2' W_gmm g2`, df = 6 - 4 = 2; the
#'     J-stat weight is `W_gmm = solve(Omega1)` (the SAME weight that found
#'     delta2), NOT `solve(Omega2)` -- the Omega1-vs-Omega2 split is verbatim from
#'     the Matlab.
#'   \item Coupling carriers: `m_hat`, its Jacobian `dm`, and `V_m = dm C dm'`
#'     (`extra_moments == 0`), plus the published `within_share` and the full
#'     paper-reporting moment vector `m_hat_full` (`extra_moments == 1`).
#' }
#'
#' @param input A two-level input list (`theta_hat`, `s`, `industry`), e.g.
#'   `gp_krw_gmm_input()`.
#' @param control Optional `gp_control`; only `control$seed` is read (to seed the
#'   multi-start). The optimizer / tolerance / max-iter are explicit arguments.
#' @param characteristic `"race"` or `"gender"` (required).
#' @param n_starts Integer; seeded jittered starts in step 1 (`0L` to use only
#'   the documented start). Default `50L`.
#' @param optimizer Passed to `gp_gmm_min_2l`; default `"BFGS"`.
#' @param reltol Relative convergence tolerance for the inner minimizations.
#' @param maxit Max iterations for the inner minimizations.
#' @return A named list (the two-level GMM fit object) with: `beta`, `mu`,
#'   `sigma_xi`, `sigma_eta`, `within_share`, `delta`, `SE`, `J_stat`, `df`,
#'   `m_hat`, `m_hat_full`, `V_m`, `v_hat`, `s_v`, `v_bar`, `s_bar`,
#'   `characteristic`, plus diagnostics `Omega1`, `Omega2`, `W_gmm`, `W_chol`,
#'   `G`, `C`, `objective`, `convergence`, and a `provenance` stamp.
#' @keywords internal
#' @noRd
gp_two_level_gmm <- function(input, control = NULL, characteristic = NULL,
                             n_starts = 50L, optimizer = "BFGS",
                             reltol = 1e-10, maxit = 1000L) {
  if (is.null(characteristic)) {
    .gradepath_abort(
      paste0("`characteristic` is required ('race' or 'gender'): the input ",
             "carries no characteristic field."),
      class = "gradepath_validation_error"
    )
  }
  characteristic <- match.arg(characteristic, c("race", "gender"))
  if (is.null(control)) control <- gp_control()

  data <- .gp_two_level_data_matrix(input)
  N <- nrow(data)
  proj <- .gp_industry_projection(data[, 1L])
  n_moments <- 6L
  n_params <- 4L
  df <- n_moments - n_params            # 6 - 4 = 2

  step_min <- function(W, st) {
    gp_gmm_min_2l(data, characteristic, W, st, proj = proj,
                  optimizer = optimizer, reltol = reltol, maxit = maxit)
  }

  ## Cluster-robust Omega1 at a step-1 candidate, with its reciprocal condition
  ## number (used to reject degenerate step-1 corners whose cluster covariance is
  ## singular -- e.g. the divergent corner the documented gender start reaches
  ## under identity weighting; see the file-header PARITY note and the multistart
  ## rationale below). `rcond_min` is the well-conditioned threshold.
  rcond_min <- 1e-10
  omega1_of <- function(delta) {
    mc <- gp_moment_conditions_2l(delta, data, characteristic,
                                  get_cov = TRUE, proj = proj)
    Omega1 <- mc$Omega
    rc <- tryCatch(rcond(Omega1), error = function(e) 0)
    if (!is.finite(rc)) rc <- 0
    list(Omega1 = Omega1, rcond = rc)
  }

  ## ----- Step 1: identity weighting (documented start + seeded multi-start) --
  ## Production `decons.m` draws randn(4,1) x1000 under identity weighting and
  ## keeps the best (the SLURM-array branch's fixed start is only the "limited
  ## one-step" path; estimate_lsqnonlin.m:61-81). We mirror that: evaluate the
  ## documented start plus `n_starts` seeded jittered starts, and -- because the
  ## optimal weight requires an invertible CLUSTER Omega1 -- select the lowest-
  ## objective step-1 optimum *among those whose Omega1 is well-conditioned*.
  ## This is load-bearing for gender, whose documented start diverges under
  ## identity-weighted BFGS to a singular-Omega1 corner; the multistart finds the
  ## well-conditioned basin (faithful to the production randn-multistart).
  ## Deterministic: seed from control$seed, else 1234 (mirrors decons.m's
  ## `rng(1234)`, the production two-level seed). The recovered optimum is
  ## seed-invariant (verified across seeds); the seed only fixes the jitter draws.
  W1 <- diag(n_moments)
  start <- gp_gmm_start_2l(characteristic)

  seed_used <- if (!is.null(control$seed)) control$seed else 1234L
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed_used)
  sd_jit <- c(0.5, 0.5, 0.2, 0.5)         # mu_raw, log_sigma_xi, beta, log_sigma_eta
  starts <- list(start)
  if (n_starts > 0L) {
    for (k in seq_len(n_starts)) {
      starts[[length(starts) + 1L]] <- start + stats::rnorm(n_params, sd = sd_jit)
    }
  }

  fit1 <- NULL; Omega1 <- NULL          # best well-conditioned step-1
  fit1_any <- NULL                      # best overall (fallback diagnostics)
  for (st in starts) {
    fk <- tryCatch(step_min(W1, st), error = function(e) NULL)
    if (is.null(fk) || !is.finite(fk$objective)) next
    if (is.null(fit1_any) || fk$objective < fit1_any$objective) fit1_any <- fk
    ok <- omega1_of(fk$delta)
    if (ok$rcond >= rcond_min &&
        (is.null(fit1) || fk$objective < fit1$objective)) {
      fit1 <- fk
      Omega1 <- ok$Omega1
    }
  }

  if (is.null(fit1)) {
    rc_any <- if (is.null(fit1_any)) NA_real_ else omega1_of(fit1_any$delta)$rcond
    .gradepath_abort(
      paste0("No step-1 optimum yielded a well-conditioned cluster-robust ",
             "Omega1 (need rcond >= ", format(rcond_min), "; best seen ",
             format(rc_any, digits = 3), "). The optimal GMM weight (Matlab ",
             "`chol(inv(Omega1))`) is undefined. Try a larger `n_starts`, or ",
             "check for too few industries / collinear moments."),
      class = "gradepath_singular_error"
    )
  }
  delta1 <- fit1$delta

  ## ----- Optimal weight from the (pre-validated) CLUSTER-ROBUST Omega1 -------
  W_gmm <- tryCatch(
    solve(Omega1),
    error = function(e) {
      .gradepath_abort(
        paste0("Step-1 cluster-robust moment covariance Omega1 is singular; ",
               "the optimal GMM weight (Matlab `chol(inv(Omega1))`) is ",
               "undefined. (", conditionMessage(e), ")"),
        class = "gradepath_singular_error"
      )
    }
  )
  W_chol <- tryCatch(
    chol(W_gmm),
    error = function(e) {
      .gradepath_abort(
        paste0("inv(Omega1) is not positive definite; `chol(inv(Omega1))` ",
               "(the Matlab two-level GMM weight factor) fails. (",
               conditionMessage(e), ")"),
        class = "gradepath_singular_error"
      )
    }
  )

  ## ----- Step 2: efficient weighting ----------------------------------------
  fit2 <- step_min(W_gmm, delta1)
  delta2 <- fit2$delta

  ## ----- Moments + cluster covariance at the step-2 optimum -----------------
  mc2 <- gp_moment_conditions_2l(delta2, data, characteristic,
                                 get_cov = TRUE, proj = proj)
  g2 <- mc2$g
  Omega2 <- mc2$Omega

  ## ----- Sandwich covariance (single N; invariant #7) -----------------------
  g_fun <- function(d) {
    gp_moment_conditions_2l(d, data, characteristic,
                            get_cov = FALSE, proj = proj)$g
  }
  G <- gp_jacobian(g_fun, delta2)             # 6 x 4

  Omega2_inv <- tryCatch(
    solve(Omega2),
    error = function(e) {
      .gradepath_abort(
        paste0("Step-2 cluster-robust moment covariance Omega2 is singular; ",
               "the sandwich covariance is undefined. (",
               conditionMessage(e), ")"),
        class = "gradepath_singular_error"
      )
    }
  )
  bread <- t(G) %*% Omega2_inv %*% G          # G' Omega2^-1 G  (4 x 4)
  C <- solve(bread) / N                       # one N (invariant #7)
  C <- (C + t(C)) / 2                         # symmetrize FP asymmetry
  SE <- sqrt(diag(C))

  ## ----- Single-N J-statistic (invariant #7; df = 2) ------------------------
  ## Uses W_gmm = solve(Omega1) (the SAME weight that found delta2), NOT
  ## solve(Omega2). df = n_moments - n_params = 2.
  J_stat <- as.numeric(N * crossprod(g2, W_gmm %*% g2))

  ## ----- Coupling carriers (computed ONCE) ----------------------------------
  m_hat <- gp_get_moments_2l(delta2, characteristic, data = data,
                             extra_moments = 0L)
  m_fun <- function(d) {
    gp_get_moments_2l(d, characteristic, data = data, extra_moments = 0L)
  }
  dm <- gp_jacobian(m_fun, delta2)            # length(m_hat) x 4
  if (is.null(dim(dm))) dm <- matrix(dm, nrow = length(m_hat))
  V_m <- dm %*% C %*% t(dm)
  V_m <- (V_m + t(V_m)) / 2                   # symmetrize

  ## Full paper-reporting moments (incl. within_share) + the scalar within_share.
  m_hat_full <- gp_get_moments_2l(delta2, characteristic, data = data,
                                  extra_moments = 1L)
  within_share <- gp_within_share(delta2, characteristic)

  ## ----- Unpack natural-scale parameters ------------------------------------
  mu_raw <- delta2[1L]
  sigma_xi <- exp(delta2[2L])
  beta <- delta2[3L]
  sigma_eta <- exp(delta2[4L])
  mu <- if (identical(characteristic, "race")) exp(mu_raw) else mu_raw

  provenance <- .gradepath_new_provenance(
    step = "two-level:6-moment cluster GMM",
    optimizer = optimizer,
    reltol = reltol,
    maxit = maxit,
    n_moments = n_moments,
    n_params = n_params,
    df = df,
    n_starts = n_starts,
    seed = seed_used,
    characteristic = characteristic,
    jacobian_backend = .gp_jacobian_backend$get(),
    convergence_step1 = fit1$convergence,
    convergence_step2 = fit2$convergence
  )

  list(
    beta = beta,
    mu = mu,
    sigma_xi = sigma_xi,
    sigma_eta = sigma_eta,
    within_share = within_share,
    delta = unname(delta2),
    SE = unname(SE),
    J_stat = J_stat,
    df = df,
    m_hat = m_hat,
    m_hat_full = m_hat_full,
    V_m = V_m,
    industry = data[, 1L],
    v_hat = mc2$v_hat,
    s_v = mc2$s_v,
    v_bar = mc2$v_bar,
    s_bar = mc2$s_bar,
    characteristic = characteristic,
    ## diagnostics
    Omega1 = Omega1,
    Omega2 = Omega2,
    W_gmm = W_gmm,
    W_chol = W_chol,
    G = G,
    C = C,
    objective = fit2$objective,
    convergence = c(step1 = fit1$convergence, step2 = fit2$convergence),
    provenance = provenance
  )
}
