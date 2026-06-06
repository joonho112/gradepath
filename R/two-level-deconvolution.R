# =============================================================================
# two-level-deconvolution.R  --  native TWO-LEVEL deconvolution prior
# -----------------------------------------------------------------------------
# gradepath-NATIVE two-level (group_fx == 1) deconvolution: a faithful base-R
# port of the `group_fx == 1` path of Kline-Rose-Walters' deconvolution machinery
# (matlab_support/estimate_lsqnonlin.m + likelihood.m + get_model_moments.m +
# mean_xi.m). It fits the JOINT G_xi (x) G_eta penalized log-spline MLE with the
# eta-integration likelihood, the GMM-driven mu+-5sigma/7sigma support caps, the
# 2-D (c_xi, c_eta) penalty grid scored by the V_m^-1 GMM criterion, and the
# mean-constraint pins (E[eta] = 1 race / 0 gender; gender also pins E[xi] = 0).
#
# This is NOT a call to ebrecipe::eb_deconvolve. It REUSES ebrecipe's validated
# numerical primitives (via the .gp_eb_* seam wrappers) as internal
# helpers -- the ns(df=5) basis, the softmax exp-spline density, the
# mean-constraint root-solve (.eb_full_alpha), and the row-logsumexp -- and is
# otherwise a fully native orchestration of the two-level eta-integration loop,
# the 2-D penalty grid, the V_m^-1 scoring, and the caps (none of which has an
# ebrecipe host surface).
#
# It CONSUMES the two-level GMM fit (gp_two_level_gmm()): the coupling
# carriers `m_hat` and `V_m`, the standardized residual `v_hat`, and the
# `s_v = s^(1-beta)` scale; plus the industry membership off the two-level
# `gp_krw_gmm_input(...)` (`input$industry`).
#
# M1-SAFETY: the two-level deconvolution code lives in THIS new file. This file
# does not modify R/deconvolution.R or any M1-cache-watched file; it reuses
# existing internals by CALLING them (after devtools::load_all) and names every
# new function distinctly (gp_deconvolve_groups, gp_two_level_likelihood,
# .gp_2l_*). The theta-pushforward g_theta is a separate hook (a clean hook is
# left; NOT built here); the sandwich/delta/bootstrap inference is separate.
#
# MATLAB GROUND-TRUTH LINE CITATIONS:
#   likelihood.m  (THE core; group_fx == 1 path):
#     * unpack T_psi/M_psi for group_fx == 1 .................... L18-L21
#     * race xi free (all T_xi); gender pins E[xi] = 0 via fsolve  L25-L29
#     * group_fx == 1 pins E[eta] = mean_target (1 race / 0 gender) L31-L39
#     * softmax g_xi / g_psi = exp(Q a - logsumexp(Q a)) ........ L43, L49
#     * per-industry loop; per-eta-node inner loop ............. L52-L91
#     * conditional matrix P = (1/S) normpdf((Y - psi*xi)/S) race  L72-L74
#                          P = (1/S) normpdf((Y - psi - xi)/S) gen L75-L77
#     * L = P * g_xi; P_hat(l) = prod(L) ....................... L78, L81
#     * P_tilde(j) = P_hat' * g_psi ............................ L87
#     * penalized objective
#         logL = -(sum(log P_tilde) - c_xi*||a_xi|| - c_psi*||a_psi||) L95
#       (UNSQUARED L2 norms; sqrt(alpha' alpha)).
#   mean_xi.m: E = sum(supp .* g) / sum(g), g = softmax(Q alpha) (the pin target).
#   get_model_moments.m: mu = sum(supp g)/sum(g);
#       sigma = sqrt(sum((supp - mu)^2 g)/sum(g)) ............... L4-L9
#   estimate_lsqnonlin.m (the DRIVER, group_fx == 1):
#     * supp_pts = 200 (two-level); T_xi = T_psi = 5 (ncol Q) ... L15, L21-L23
#     * RACE  c_xi 0.08:0.005:0.12 (9) x c_psi 0.0025:0.0025:0.02 (8) = 72;
#         optimum c_xi = 0.105, c_psi = 0.0025 ................. L28-L29
#     * GENDER c_xi 0.01:0.0025:0.03 (9) x c_psi 0.0025:0.0025:0.03 (12) = 108;
#         optimum c_xi = 0.020, c_psi = 0.0025 ................. L38-L39
#     * empirical psi_hat (industry-weighted mean of v_hat) and
#         xi_hat (= v_hat/psi_hat race / v_hat - psi_hat gender)  L194-L219
#     * support caps (mu_xi = m_hat(1), sigma_xi = m_hat(2),
#         sigma_psi = m_hat(3) race / sigma_xi = m_hat(1),
#         sigma_psi = m_hat(2) gender):
#         RACE xi: min=0, max=min(max(max(xi_hat), mu+5s), mu+7s) L223-L242
#         RACE eta:min=0, max=min(max(max(psi_hat),1+5s),1+7s)
#         GENDER: symmetric +-5s/+-7s about 0 w/ empirical min/max L243-L262
#         nesting supp_min = max(min(min_1,min_2), min_cap);
#                 supp_max = min(max(max_1,max_2), max_cap) ..... L263-L266
#     * equispaced 200-pt supp_xi / supp_psi grids ............. L277-L288
#     * 2-D grid search: combs = combvec(c_xi_grid, c_psi_grid)  L350
#         per comb: spline MLE; m_hat_curr = get_model_moments of
#         g_xi (and g_psi(2) for eta); Js = (m_hat - m_hat_curr)'
#         inv(V_m) (m_hat - m_hat_curr); min-J wins ............ L358-L403
#     * final normalisation g_xi = g_xi / sum(g_xi) ............ L469, L474
#
# OUTPUT: a `gp_prior`-shaped object (the Step-1.5 shell) carrying support
# (= supp_xi), density (= g_xi probability masses summing to 1), mean, scale="r",
# with the two-level fields (support_eta, g_eta, c_xi, c_eta, J grid/score,
# group_fx = 1, characteristic, model_moments) riding in diagnostics/metadata
# (the eb_prior-shaped shell tolerates these advisory two-level fields; see
# R/class-objects-output.R "full-vs-shell").
# =============================================================================

#' Resolve the package abort helper whether sourced or namespaced
#'
#' Mirrors `.gp_decon_abort` (R/deconvolution.R) so this file works both under
#' `devtools::load_all()` (namespace) and when `source()`-d on top of it.
#'
#' @keywords internal
#' @noRd
.gp_2l_abort <- function(msg, class = "gradepath_error") {
  fn <- tryCatch(get(".gradepath_abort", envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (is.function(fn)) return(fn(msg, class = class))
  if (exists(".gradepath_abort", mode = "function")) {
    return(get(".gradepath_abort")(msg, class = class))
  }
  cnd <- structure(
    class = c(class, "gradepath_error", "error", "condition"),
    list(message = msg, call = NULL)
  )
  stop(cnd)
}

#' Resolve a package internal by name (namespace-first), erroring if absent
#' @keywords internal
#' @noRd
.gp_2l_get <- function(name) {
  fn <- tryCatch(get(name, envir = asNamespace("gradepath")),
                 error = function(e) NULL)
  if (!is.null(fn)) return(fn)
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  .gp_2l_abort(
    sprintf("Required gradepath internal '%s' not found.", name),
    class = "gradepath_internal_error"
  )
}

# ---------------------------------------------------------------------------
# Empirical component estimates (psi_hat / xi_hat) -- estimate_lsqnonlin.m:194-219
# ---------------------------------------------------------------------------

#' Industry-weighted empirical group/within component estimates
#'
#' @description
#' Port of `estimate_lsqnonlin.m:194-219` (the `group_fx == 1` block). The
#' industry-weighted mean of `v_hat` is the empirical group effect `psi_hat`
#' (per firm, constant within an industry, equal-weighted: `w_j = 1/n_j`):
#'   RACE   `psi_j = mean_{i in j}(v_hat_i) / mu`,  `xi_hat = v_hat / psi_hat`
#'   GENDER `psi_j = mean_{i in j}(v_hat_i)`,       `xi_hat = v_hat - psi_hat`
#' These feed ONLY the empirical `_1` cap terms (`max(psi_hat)`, `max(xi_hat)`,
#' and for gender also the mins); the GMM-implied `_2`/`_cap` terms come from
#' `m_hat`. `mu = m_hat(1)` for race (the natural-scale mean); gender has no mu.
#'
#' @param v_hat Length-N standardized residual carrier (`fit$v_hat`).
#' @param industry Length-N integer industry membership (re-keyed to dense 1..J).
#' @param characteristic "race" or "gender".
#' @param mu Natural-scale mean (race `m_hat(1)`; ignored for gender).
#' @return list(psi_hat, xi_hat), each length-N.
#' @keywords internal
#' @noRd
.gp_2l_components <- function(v_hat, industry, characteristic, mu = NULL) {
  idx <- as.integer(match(industry, sort(unique(industry), method = "radix")))
  J <- max(idx)
  sizes <- tabulate(idx, nbins = J)
  ## equal-weighted industry mean of v_hat (w_j = 1/n_j, so weighted sum == mean).
  grp_keys <- sort(unique(idx), method = "radix")
  vbar_by_key <- as.numeric(rowsum(v_hat, group = idx)) / sizes[grp_keys]
  vbar <- vbar_by_key[match(idx, grp_keys)]      # per-firm industry mean of v_hat

  if (identical(characteristic, "race")) {
    if (is.null(mu) || !is.finite(mu) || mu == 0) {
      .gp_2l_abort(
        "Race empirical components require a finite non-zero `mu` (= m_hat[1]).",
        class = "gradepath_validation_error"
      )
    }
    psi_hat <- vbar / mu                          # psi_j = industry-mean / mu
    xi_hat <- v_hat / psi_hat                      # xi_hat = v_hat / psi_hat
  } else {
    psi_hat <- vbar                                # psi_j = industry mean
    xi_hat <- v_hat - psi_hat                       # xi_hat = v_hat - psi_hat
  }
  list(psi_hat = psi_hat, xi_hat = xi_hat)
}

# ---------------------------------------------------------------------------
# GMM-driven support caps -- estimate_lsqnonlin.m:222-266 (two-level)
# ---------------------------------------------------------------------------

#' Two-level GMM-driven support caps for xi and eta (the mu+-5/7 sigma caps)
#'
#' @description
#' Faithful port of `estimate_lsqnonlin.m:222-266` for `group_fx == 1`.
#'
#' RACE (L223-242): `mu_xi = m_hat[1]`, `sigma_xi = m_hat[2]`,
#'   `sigma_psi = m_hat[3]`.
#'   xi : min_1=min_2=min_cap=0; max_1=max(xi_hat); max_2=mu_xi+5 sigma_xi;
#'        max_cap=mu_xi+7 sigma_xi.
#'   eta: min_1=min_2=min_cap=0; max_1=max(psi_hat); max_2=1+5 sigma_psi;
#'        max_cap=1+7 sigma_psi.
#' GENDER (L243-262): `sigma_xi = m_hat[1]`, `sigma_psi = m_hat[2]`, centered 0.
#'   xi : min_1=min(xi_hat); min_2=-5 sigma_xi; min_cap=-7 sigma_xi;
#'        max_1=max(xi_hat); max_2=5 sigma_xi; max_cap=7 sigma_xi.
#'   eta: min_1=min(psi_hat); min_2=-5 sigma_psi; min_cap=-7 sigma_psi;
#'        max_1=max(psi_hat); max_2=5 sigma_psi; max_cap=7 sigma_psi.
#' Nesting (L263-266): `supp_min = max(min(min_1, min_2), min_cap)`;
#'                     `supp_max = min(max(max_1, max_2), max_cap)`.
#'
#' @param m_hat Coupling moments (race `c(mu, sigma_xi, sigma_eta)`, length 3;
#'   gender `c(sigma_xi, sigma_eta)`, length 2) -- exactly the Step-7.1 layout.
#' @param psi_hat,xi_hat Empirical component estimates from
#'   `.gp_2l_components()`.
#' @param characteristic "race" or "gender".
#' @return list(xi = c(lo, hi), eta = c(lo, hi)).
#' @keywords internal
#' @noRd
.gp_2l_support_caps <- function(m_hat, psi_hat, xi_hat, characteristic) {
  if (identical(characteristic, "race")) {
    if (length(m_hat) < 3L) {
      .gp_2l_abort(
        sprintf("Two-level race caps need m_hat length 3 (mu, sigma_xi, sigma_eta); got %d.",
                length(m_hat)),
        class = "gradepath_validation_error"
      )
    }
    mu_xi_hat <- m_hat[1L]
    sigma_xi_hat <- m_hat[2L]
    sigma_psi_hat <- m_hat[3L]
    ## eta (psi) caps -- centered at 1.
    psi_min <- max(min(0, 0), 0)
    psi_max <- min(max(max(psi_hat), 1 + 5 * sigma_psi_hat), 1 + 7 * sigma_psi_hat)
    ## xi caps -- centered at mu_xi_hat, floored at 0.
    xi_min <- max(min(0, 0), 0)
    xi_max <- min(max(max(xi_hat), mu_xi_hat + 5 * sigma_xi_hat),
                  mu_xi_hat + 7 * sigma_xi_hat)
  } else if (identical(characteristic, "gender")) {
    if (length(m_hat) < 2L) {
      .gp_2l_abort(
        sprintf("Two-level gender caps need m_hat length 2 (sigma_xi, sigma_eta); got %d.",
                length(m_hat)),
        class = "gradepath_validation_error"
      )
    }
    sigma_xi_hat <- m_hat[1L]
    sigma_psi_hat <- m_hat[2L]
    ## eta (psi) caps -- symmetric about 0.
    psi_min <- max(min(min(psi_hat), -5 * sigma_psi_hat), -7 * sigma_psi_hat)
    psi_max <- min(max(max(psi_hat), 5 * sigma_psi_hat), 7 * sigma_psi_hat)
    ## xi caps -- symmetric about 0.
    xi_min <- max(min(min(xi_hat), -5 * sigma_xi_hat), -7 * sigma_xi_hat)
    xi_max <- min(max(max(xi_hat), 5 * sigma_xi_hat), 7 * sigma_xi_hat)
  } else {
    .gp_2l_abort(
      sprintf("Unknown `characteristic`: '%s' (expected 'race' or 'gender').",
              characteristic),
      class = "gradepath_validation_error"
    )
  }
  if (!(xi_min < xi_max) || !(psi_min < psi_max) ||
      any(!is.finite(c(xi_min, xi_max, psi_min, psi_max)))) {
    .gp_2l_abort(
      sprintf(paste0("Degenerate two-level support caps ",
                     "(xi: [%g, %g], eta: [%g, %g]); need lo < hi and finite."),
              xi_min, xi_max, psi_min, psi_max),
      class = "gradepath_validation_error"
    )
  }
  list(xi = c(lo = xi_min, hi = xi_max), eta = c(lo = psi_min, hi = psi_max))
}

# ---------------------------------------------------------------------------
# Equispaced support grids -- estimate_lsqnonlin.m:277-288
# ---------------------------------------------------------------------------

#' Equispaced support grid (linspace lo->hi inclusive)
#' @keywords internal
#' @noRd
.gp_2l_grid <- function(caps, supp_pts = 200L) {
  lo <- unname(caps[["lo"]]); hi <- unname(caps[["hi"]])
  seq(lo, hi, length.out = as.integer(supp_pts))
}

# ---------------------------------------------------------------------------
# Prior-implied moments -- get_model_moments.m (higher_moments == 0)
# ---------------------------------------------------------------------------

#' Mean + SD of a discrete density on a grid (get_model_moments.m:4-9)
#' @keywords internal
#' @noRd
.gp_2l_model_moments <- function(support, g) {
  sg <- sum(g)
  mu <- sum(support * g) / sg
  sigma <- sqrt(sum(((support - mu)^2) * g) / sg)
  c(mean = mu, sd = sigma)
}

#' Build the model moment vector m_hat_curr to compare against m_hat
#'
#' @description
#' Port of `estimate_lsqnonlin.m:370-377`. From the fitted densities:
#'   RACE   `m_hat_curr = [mu_xi, sigma_xi, sigma_eta]` (length 3) -- the MEAN
#'          and SD of g_xi, then the SD of g_eta (g_psi).
#'   GENDER `m_hat_curr = [sigma_xi, sigma_eta]` (length 2) -- the Matlab keeps
#'          only `m_hat_curr(2)` (the SD of g_xi; the mean is pinned to 0), then
#'          appends the SD of g_eta (g_psi).
#'
#' @param supp_xi,g_xi Within-firm support + density.
#' @param supp_eta,g_eta Between-industry support + density.
#' @param characteristic "race" or "gender".
#' @return Numeric vector matching the SHAPE of the two-level `m_hat`.
#' @keywords internal
#' @noRd
.gp_2l_model_moment_vector <- function(supp_xi, g_xi, supp_eta, g_eta,
                                       characteristic) {
  mom_xi <- .gp_2l_model_moments(supp_xi, g_xi)
  mom_eta <- .gp_2l_model_moments(supp_eta, g_eta)
  if (identical(characteristic, "race")) {
    ## [mu_xi, sigma_xi, sigma_eta]
    c(unname(mom_xi[["mean"]]), unname(mom_xi[["sd"]]), unname(mom_eta[["sd"]]))
  } else {
    ## gender: keep only m_hat_curr(2) (sigma_xi), then sigma_eta
    c(unname(mom_xi[["sd"]]), unname(mom_eta[["sd"]]))
  }
}

#' The coupling J for a two-level fit (estimate_lsqnonlin.m:378)
#'
#' `J = (m_hat - m_model)' solve(V_m) (m_hat - m_model)`.
#' @keywords internal
#' @noRd
.gp_2l_solve_V_m <- function(V_m, expected_dim = NULL) {
  if (!is.matrix(V_m) || nrow(V_m) != ncol(V_m)) {
    .gp_2l_abort(
      "`V_m` must be a square numeric matrix for the two-level penalty score.",
      class = "gradepath_validation_error"
    )
  }
  if (!is.null(expected_dim) && !identical(nrow(V_m), as.integer(expected_dim))) {
    .gp_2l_abort(
      sprintf("`V_m` dimension %d does not match the expected coupling dimension %d.",
              nrow(V_m), as.integer(expected_dim)),
      class = "gradepath_validation_error"
    )
  }
  rc <- tryCatch(rcond(V_m), error = function(e) 0)
  if (!is.finite(rc) || rc < 1e-12) {
    .gp_2l_abort(
      sprintf("`V_m` is numerically singular (rcond = %.3g); cannot score the two-level penalty grid with V_m^-1.", rc),
      class = "gradepath_singular_error"
    )
  }
  tryCatch(
    solve(V_m),
    error = function(e) {
      .gp_2l_abort(
        paste0("`V_m` could not be inverted for the two-level penalty score. (",
               conditionMessage(e), ")"),
        class = "gradepath_singular_error"
      )
    }
  )
}

#' @keywords internal
#' @noRd
.gp_2l_J <- function(m_hat, V_m, m_model, V_m_inv = NULL) {
  d <- as.numeric(m_hat) - as.numeric(m_model)
  if (is.null(V_m_inv)) V_m_inv <- .gp_2l_solve_V_m(V_m, expected_dim = length(d))
  as.numeric(crossprod(d, V_m_inv %*% d))
}

# ---------------------------------------------------------------------------
# The two-level eta-integration likelihood -- likelihood.m (group_fx == 1)
# ---------------------------------------------------------------------------

#' Two-level penalized negative log-likelihood (port of likelihood.m group_fx==1)
#'
#' @description
#' Faithful, HEAVILY VECTORIZED port of the `group_fx == 1` path of
#' `likelihood.m`. For free coefficient vectors `(alpha_xi_free, alpha_eta_free)`
#' it (a) rebuilds the FULL constrained coefficient vectors via the
#' mean-constraint root-solve (the `.eb_full_alpha` pin), (b) maps them to the
#' softmax densities `g_xi`, `g_eta`, (c) forms the per-industry eta-integrated
#' log-likelihood, and (d) returns the NEGATIVE penalized objective
#' `-(sum(log P_tilde) - c_xi ||alpha_xi|| - c_eta ||alpha_eta||)`
#' (likelihood.m:95; UNSQUARED L2 norms on the FULL constraint-completed vectors).
#'
#' Mean pins (likelihood.m:25-39):
#'   * RACE: `alpha_xi` is FREE (all T_xi coefficients; no xi mean constraint);
#'     `alpha_eta` is pinned so `E(eta) = 1`.
#'   * GENDER: `alpha_xi` is pinned so `E(xi) = 0` (the last xi coef is solved);
#'     `alpha_eta` is pinned so `E(eta) = 0`.
#'
#' Per-industry likelihood (likelihood.m:52-91), vectorized over firms x xi-nodes:
#'   For firm i and xi-node m, `P[i,m] = (1/s_v_i) phi((r_i - psi*supp_xi_m)/s_v_i)`
#'   (race; psi = supp_eta_l) or `phi((r_i - psi - supp_xi_m)/s_v_i)` (gender).
#'   `L[i] = sum_m P[i,m] g_xi[m]` (= P %*% g_xi). Per industry k and eta-node l:
#'   `log P_hat[k,l] = sum_{i in k} log L_l[i]` (the `prod(L)` over firms, in
#'   LOGS to avoid underflow). Then
#'   `log P_tilde[k] = logsumexp_l( log g_eta[l] + log P_hat[k,l] )`.
#'   `loglik = sum_k log P_tilde[k]`.
#'
#' We work entirely in LOGS: per eta-node we build the N x M_xi log-kernel
#'   `log P = -log(s_v) - 0.5 log(2 pi) - 0.5 ((r - psi*supp_xi)/s_v)^2`
#' and `log L_l[i] = logsumexp_m( log P[i,m] + log g_xi[m] )` (a row-logsumexp),
#' then sum within industry. This avoids the Matlab `prod(L)` underflow and the
#' nested R firm-loop; the only outer loop is over the M_eta eta-nodes.
#'
#' @param alpha_xi_free Free within-firm coefficients (length T_xi race;
#'   T_xi - 1 gender, the last pinned by the `E(xi) = 0` constraint).
#' @param alpha_eta_free Free between-industry coefficients (length T_eta - 1;
#'   the last pinned by the `E(eta) = mean_target` constraint).
#' @param r Length-N studentized residual carrier (`fit$v_hat`).
#' @param s_v Length-N scaled standard errors (`fit$s_v`).
#' @param industry_idx Length-N dense 1..J industry membership.
#' @param Q_xi,supp_xi Within-firm basis + grid.
#' @param Q_eta,supp_eta Between-industry basis + grid.
#' @param c_xi,c_eta Scalar penalties.
#' @param characteristic "race" or "gender".
#' @param return_densities If `TRUE`, also return `g_xi` / `g_eta` (for the
#'   moment scoring) alongside the objective.
#' @return If `return_densities` is FALSE, the scalar NEGATIVE penalized
#'   objective (to MINIMISE). Else a list(objective, g_xi, g_eta, loglik).
#' @keywords internal
#' @noRd
gp_two_level_likelihood <- function(alpha_xi_free, alpha_eta_free,
                                    r, s_v, industry_idx,
                                    Q_xi, supp_xi, Q_eta, supp_eta,
                                    c_xi, c_eta, characteristic,
                                    return_densities = FALSE) {
  eb_softmax <- .gp_2l_get(".gp_eb_softmax_density")
  eb_full    <- .gp_2l_get(".gp_eb_full_alpha")
  eb_rowlse  <- .gp_2l_get(".gp_eb_row_log_sum_exp")
  is_race <- identical(characteristic, "race")

  ## --- full (constraint-completed) coefficient vectors + L2 norms -----------
  ## RACE: alpha_xi free (all T_xi). GENDER: alpha_xi pinned E[xi] = 0.
  alpha_xi <- if (is_race) {
    as.numeric(alpha_xi_free)
  } else {
    eb_full(alpha_xi_free, Q_xi, supp_xi, target_mean = 0)
  }
  ## eta pinned: E[eta] = 1 (race) / 0 (gender).
  eta_target <- if (is_race) 1 else 0
  alpha_eta <- eb_full(alpha_eta_free, Q_eta, supp_eta, target_mean = eta_target)

  ## --- softmax densities (likelihood.m:43/49) -------------------------------
  sx <- eb_softmax(Q_xi, alpha_xi)
  se <- eb_softmax(Q_eta, alpha_eta)
  g_xi <- if (is.list(sx)) sx[["g"]] else sx
  g_eta <- if (is.list(se)) se[["g"]] else se
  log_g_xi <- log(g_xi)
  log_g_eta <- log(g_eta)

  N <- length(r)
  M_xi <- length(supp_xi)
  M_eta <- length(supp_eta)
  J <- max(industry_idx)

  ## --- per-eta-node log-likelihood (vectorized over firms x xi-nodes) -------
  ## logP_hat: J x M_eta, where logP_hat[k, l] = sum_{i in k} log L_l[i].
  log_const <- -log(s_v) - 0.5 * log(2 * pi)     # length N
  logP_hat <- matrix(0, nrow = J, ncol = M_eta)

  for (l in seq_len(M_eta)) {
    psi <- supp_eta[l]
    ## mean of the within-firm normal kernel at each xi-node, per firm:
    ##   race:   mu_im = psi * supp_xi_m   -> argument (r_i - psi*supp_xi_m)/s_v_i
    ##   gender: mu_im = psi + supp_xi_m   -> argument (r_i - psi - supp_xi_m)/s_v_i
    ## Build the N x M_xi standardized argument z, then log P = log_const - 0.5 z^2.
    if (is_race) {
      ## outer(r, psi*supp_xi, "-") = r_i - psi*supp_xi_m
      Mmean <- outer(r, psi * supp_xi, FUN = "-")        # N x M_xi
    } else {
      Mmean <- outer(r - psi, supp_xi, FUN = "-")         # N x M_xi
    }
    z <- Mmean / s_v                                       # recycle s_v down rows
    logP <- log_const - 0.5 * (z^2)                        # N x M_xi (log normal kernel)
    ## log L_l[i] = logsumexp_m( logP[i,m] + log_g_xi[m] )  (row-logsumexp)
    logL <- eb_rowlse(sweep(logP, 2L, log_g_xi, FUN = "+"))   # length N
    ## sum_{i in k} log L_l[i]  -> per-industry column l of logP_hat
    logP_hat[, l] <- as.numeric(rowsum(logL, group = industry_idx))
  }

  ## --- integrate out eta: log P_tilde[k] = logsumexp_l( log_g_eta[l] + logP_hat[k,l] )
  logP_tilde <- eb_rowlse(sweep(logP_hat, 2L, log_g_eta, FUN = "+"))   # length J
  loglik <- sum(logP_tilde)

  ## --- penalized objective (likelihood.m:95) --------------------------------
  ## UNSQUARED L2 norms on the FULL constraint-completed coefficient vectors.
  norm_xi <- sqrt(sum(alpha_xi^2))
  norm_eta <- sqrt(sum(alpha_eta^2))
  penalized <- loglik - c_xi * norm_xi - c_eta * norm_eta
  objective <- -penalized                          # NEGATIVE -> minimise

  if (!return_densities) return(objective)
  list(objective = objective, g_xi = g_xi, g_eta = g_eta, loglik = loglik)
}

# ---------------------------------------------------------------------------
# Single penalized joint MLE at one (c_xi, c_eta) node
# ---------------------------------------------------------------------------

#' Fit the joint G_xi (x) G_eta softmax MLE at one penalty node
#'
#' @description
#' Minimises `gp_two_level_likelihood()` over the FREE coefficient blocks at a
#' fixed `(c_xi, c_eta)` (the Matlab `fminunc` at one grid node;
#' `estimate_lsqnonlin.m:366`). The free-parameter layout matches the Matlab
#' `omega` packing (`likelihood.m:24-39`):
#'   * RACE:   `omega = c(alpha_xi(T_xi free), alpha_eta(T_eta - 1 free))`
#'             (alpha_xi all free; the LAST alpha_eta coef pinned `E(eta) = 1`).
#'   * GENDER: `omega = c(alpha_xi(T_xi - 1 free), alpha_eta(T_eta - 1 free))`
#'             (the LAST alpha_xi coef pinned `E(xi) = 0`; the LAST alpha_eta coef
#'             pinned `E(eta) = 0`).
#' We minimise with `stats::optim` (BFGS), the base-R analogue of `fminunc`.
#'
#' @return list(g_xi, g_eta, alpha_xi_free, alpha_eta_free, objective, loglik,
#'   converged, convergence_code, convergence_message).
#' @keywords internal
#' @noRd
.gp_2l_fit_node <- function(c_xi, c_eta, r, s_v, industry_idx,
                            Q_xi, supp_xi, Q_eta, supp_eta,
                            characteristic, start,
                            max_iter = 5000L, tol = 1e-8) {
  is_race <- identical(characteristic, "race")
  T_xi <- ncol(Q_xi)
  T_eta <- ncol(Q_eta)
  n_xi_free <- if (is_race) T_xi else (T_xi - 1L)
  n_eta_free <- T_eta - 1L
  n_free <- n_xi_free + n_eta_free

  if (length(start) != n_free) {
    .gp_2l_abort(
      sprintf("Two-level start length %d != expected free-param count %d.",
              length(start), n_free),
      class = "gradepath_validation_error"
    )
  }

  obj <- function(par) {
    a_xi <- par[seq_len(n_xi_free)]
    a_eta <- par[(n_xi_free + 1L):n_free]
    val <- tryCatch(
      gp_two_level_likelihood(a_xi, a_eta, r, s_v, industry_idx,
                              Q_xi, supp_xi, Q_eta, supp_eta,
                              c_xi, c_eta, characteristic,
                              return_densities = FALSE),
      error = function(e) NA_real_
    )
    ## Fence non-finite / erroring objective values so a pathological coefficient
    ## draw (e.g. a degenerate softmax or a failed mean-pin root-solve) cannot
    ## abort the optimizer mid-sweep; optim then steps away from it.
    if (!is.finite(val)) 1e10 else val
  }

  opt <- stats::optim(
    par = start, fn = obj, method = "BFGS",
    control = list(maxit = as.integer(max_iter), reltol = tol)
  )

  a_xi <- opt$par[seq_len(n_xi_free)]
  a_eta <- opt$par[(n_xi_free + 1L):n_free]
  final <- gp_two_level_likelihood(a_xi, a_eta, r, s_v, industry_idx,
                                   Q_xi, supp_xi, Q_eta, supp_eta,
                                   c_xi, c_eta, characteristic,
                                   return_densities = TRUE)

  list(g_xi = final$g_xi, g_eta = final$g_eta,
       alpha_xi_free = a_xi, alpha_eta_free = a_eta,
       objective = opt$value, loglik = final$loglik,
       converged = isTRUE(opt$convergence == 0L),
       convergence_code = as.integer(opt$convergence),
       convergence_message = if (is.null(opt$message)) NA_character_ else as.character(opt$message))
}

# ---------------------------------------------------------------------------
# 2-D penalty grid -- estimate_lsqnonlin.m:28-44, 350 (combvec)
# ---------------------------------------------------------------------------

#' Default two-level 2-D penalty grid (combvec(c_xi_grid, c_eta_grid))
#'
#' @description
#' Matches the Matlab two-level grids (`estimate_lsqnonlin.m:28-29, 38-39`):
#'   RACE   c_xi = seq(0.08, 0.12, by = 0.005)   (9) x
#'          c_eta = seq(0.0025, 0.02, by = 0.0025) (8) = 72 nodes.
#'   GENDER c_xi = seq(0.01, 0.03, by = 0.0025)  (9) x
#'          c_eta = seq(0.0025, 0.03, by = 0.0025) (12) = 108 nodes.
#' The Matlab `combvec` varies the FIRST argument (c_xi) fastest; we reproduce
#' that ordering (`expand.grid` varies its first column fastest) so the node
#' indexing matches the reference. Overridable through `control`.
#'
#' @param characteristic "race" or "gender".
#' @param control Optional list; reads `deconv2l_penalty_grid_race/_gender`
#'   (each a list(c_xi = <vec>, c_eta = <vec>)).
#' @return A data.frame with columns `c_xi`, `c_eta` (one row per node), c_xi
#'   varying fastest.
#' @keywords internal
#' @noRd
.gp_2l_penalty_grid <- function(characteristic, control = NULL) {
  override <- NULL
  if (!is.null(control)) {
    override <- if (identical(characteristic, "race")) {
      control$deconv2l_penalty_grid_race
    } else {
      control$deconv2l_penalty_grid_gender
    }
  }
  if (is.list(override) && length(override$c_xi) >= 1L &&
      length(override$c_eta) >= 1L) {
    cx <- as.numeric(override$c_xi)
    ce <- as.numeric(override$c_eta)
  } else if (identical(characteristic, "race")) {
    cx <- seq(0.08, 0.12, by = 0.005)
    ce <- seq(0.0025, 0.02, by = 0.0025)
  } else {
    cx <- seq(0.01, 0.03, by = 0.0025)
    ce <- seq(0.0025, 0.03, by = 0.0025)
  }
  ## combvec varies the FIRST arg fastest; expand.grid varies its 1st col fastest.
  g <- expand.grid(c_xi = cx, c_eta = ce, KEEP.OUT.ATTRS = FALSE,
                   stringsAsFactors = FALSE)
  g[, c("c_xi", "c_eta")]
}

# ---------------------------------------------------------------------------
# The native two-level deconvolution driver
# ---------------------------------------------------------------------------

#' Native two-level (group_fx == 1) deconvolution prior from a GMM fit
#'
#' @description
#' Takes the two-level GMM fit
#' (`gp_two_level_gmm()`) and the two-level input (`gp_krw_gmm_input()`) and
#' produces the JOINT mixing-density priors `G_xi` (within-firm) and `G_eta`
#' (between-industry) on the r/xi-scale, by penalized joint log-spline maximum
#' likelihood with the two penalties `(c_xi, c_eta)` selected to minimise the
#' coupling J (the `V_m^-1` GMM-moment fit of the prior-implied moments).
#' Faithful base-R port of the `group_fx == 1` branch of
#' `estimate_lsqnonlin.m` / `likelihood.m` (see file header for line citations).
#'
#' @section What this consumes (the coupling):
#' off `fit` -- `m_hat` (coupling moments: race `c(mu, sigma_xi, sigma_eta)`;
#' gender `c(sigma_xi, sigma_eta)`), `V_m` (their covariance), `v_hat` (the
#' standardized residual), `s_v` (= s^(1-beta)), `beta`, `characteristic`. Off
#' `input` -- the `industry` membership (re-keyed to dense 1..J).
#'
#' @param fit A two-level GMM fit (`gp_two_level_gmm()` output); must be
#'   two-level (`m_hat` length 3 race / 2 gender).
#' @param input Optional two-level input list (`theta_hat`, `s`, `industry`), e.g.
#'   `gp_krw_gmm_input()` -- read ONLY for `industry`. When omitted, the dense
#'   `fit$industry` carrier returned by `gp_two_level_gmm()` is used.
#' @param control Optional list; reads `deconv2l_supp_pts`, `deconv2l_n_knots`,
#'   `deconv2l_max_iter`, `tol`, `seed`, and the penalty-grid overrides. Falls
#'   back to the Matlab defaults when absent.
#' @param n_starts Positive integer; starts per penalty node. `n_starts = 1L`
#'   (DEFAULT) uses a single deterministic zeros-start (reproducible). `> 1L`
#'   ALSO fits from `n_starts - 1` seeded `rnorm` starts per node and keeps the
#'   lowest-objective fit (mirroring KRW's per-node `randn` restart,
#'   `estimate_lsqnonlin.m:356`); seeded once from `control$seed` (else 1234).
#' @param supp_pts Integer support-grid size (200 production; tests pass fewer).
#' @param penalties `NULL` (DEFAULT) runs the full 2-D V_m^-1 grid; or a
#'   `list(c_xi = , c_eta = )` to fit at a FIXED penalty and SKIP the grid
#'   (essential for bounded tests / the fixed-optimal-penalty parity solve).
#' @return A validated `gp_prior` with `support` (= supp_xi), `density` (= g_xi
#'   probability masses summing to 1), `mean`, `scale = "r"`, and the two-level
#'   fields in `diagnostics` (`support_eta`, `g_eta`, selected `c_xi`/`c_eta`,
#'   `J`, the J/penalty grid, `model_moments`, `caps`, `group_fx = 1`) and
#'   `metadata` (`characteristic`, `beta`, `alpha_xi_free`/`alpha_eta_free`,
#'   spline dims). The theta-pushforward `g_theta` is a separate hook (a clean
#'   hook; `diagnostics$g_theta` is left NULL here). Inference is separate.
#' @keywords internal
#' @noRd
gp_deconvolve_groups <- function(fit, input = NULL, control = NULL, n_starts = 1L,
                                 supp_pts = 200L, penalties = NULL) {
  ## --- validate n_starts ---------------------------------------------------
  if (length(n_starts) != 1L || !is.numeric(n_starts) ||
      !is.finite(n_starts) || n_starts < 1 ||
      n_starts != as.integer(n_starts)) {
    .gp_2l_abort(
      sprintf("`n_starts` must be a single positive integer (>= 1); got %s.",
              paste(format(n_starts), collapse = ", ")),
      class = "gradepath_validation_error"
    )
  }
  n_starts <- as.integer(n_starts)

  ## --- validate the GMM fit (two-level contract) ---------------------------
  if (is.null(fit$characteristic) || is.null(fit$m_hat) || is.null(fit$V_m) ||
      is.null(fit$v_hat) || is.null(fit$s_v)) {
    .gp_2l_abort(
      paste0("`fit` must be a Step-7.1 two-level GMM fit (gp_two_level_gmm()); ",
             "it must carry characteristic, m_hat, V_m, v_hat, s_v."),
      class = "gradepath_validation_error"
    )
  }
  characteristic <- fit$characteristic
  if (!characteristic %in% c("race", "gender")) {
    .gp_2l_abort(
      sprintf("Unknown `characteristic`: '%s'.", characteristic),
      class = "gradepath_validation_error"
    )
  }
  ## two-level shape: race m_hat length 3, gender length 2.
  expected_m <- if (identical(characteristic, "race")) 3L else 2L
  if (length(fit$m_hat) != expected_m) {
    .gp_2l_abort(
      sprintf(paste0("gp_deconvolve_groups() is two-level (group_fx == 1): ",
                     "expected m_hat length %d for %s, got %d. (A one-level fit ",
                     "should use gp_deconvolve().)"),
              expected_m, characteristic, length(fit$m_hat)),
      class = "gradepath_validation_error"
    )
  }
  raw_industry <- if (!is.null(input) && !is.null(input$industry)) {
    input$industry
  } else {
    fit$industry
  }
  if (is.null(raw_industry)) {
    .gp_2l_abort(
      paste0("The two-level deconvolution needs an industry membership vector. ",
             "Pass `input = gp_krw_gmm_input(...)`, or use a Step-7.1 fit that ",
             "carries `fit$industry`."),
      class = "gradepath_validation_error"
    )
  }

  ## --- carriers off the fit ------------------------------------------------
  v_hat <- as.numeric(fit$v_hat)               # the studentized residual r
  s_v   <- as.numeric(fit$s_v)
  m_hat <- as.numeric(fit$m_hat)
  V_m   <- fit$V_m
  industry <- raw_industry
  if (length(industry) != length(v_hat)) {
    .gp_2l_abort(
      sprintf("The industry membership vector (%d) must be length-matched to fit$v_hat (%d).",
              length(industry), length(v_hat)),
      class = "gradepath_validation_error"
    )
  }
  industry_idx <- as.integer(match(industry, sort(unique(industry), method = "radix")))
  if (!is.null(fit$industry)) {
    fit_industry <- as.integer(fit$industry)
    if (length(fit_industry) != length(industry_idx) ||
        !identical(fit_industry, industry_idx)) {
      .gp_2l_abort(
        paste0("`input$industry` does not match the dense industry ordering carried ",
               "by `fit$industry`; rebuild the Step-7.1 fit from the same ordered input."),
        class = "gradepath_validation_error"
      )
    }
  }
  V_m_inv <- .gp_2l_solve_V_m(V_m, expected_dim = expected_m)

  ## --- control knobs (Matlab defaults when absent) -------------------------
  supp_pts <- if (!is.null(control$deconv2l_supp_pts)) {
    as.integer(control$deconv2l_supp_pts)
  } else as.integer(supp_pts)
  n_knots <- if (!is.null(control$deconv2l_n_knots)) {
    as.integer(control$deconv2l_n_knots)
  } else 5L
  max_iter <- if (!is.null(control$deconv2l_max_iter)) {
    as.integer(control$deconv2l_max_iter)
  } else 5000L
  opt_tol <- if (!is.null(control$tol)) control$tol else 1e-8

  ## --- empirical components + GMM-driven caps + grids ----------------------
  mu_race <- if (identical(characteristic, "race")) m_hat[1L] else NULL
  comp <- .gp_2l_components(v_hat, industry, characteristic, mu = mu_race)
  caps <- .gp_2l_support_caps(m_hat, comp$psi_hat, comp$xi_hat, characteristic)
  supp_xi  <- .gp_2l_grid(caps$xi,  supp_pts = supp_pts)
  supp_eta <- .gp_2l_grid(caps$eta, supp_pts = supp_pts)

  ## --- bases via the ebrecipe primitive (ns df = 5, center + unit-norm) ----
  eb_basis <- .gp_2l_get(".gp_eb_spline_basis")
  Q_xi  <- eb_basis(supp_xi,  n_knots = n_knots)
  Q_eta <- eb_basis(supp_eta, n_knots = n_knots)
  T_xi <- ncol(Q_xi)
  T_eta <- ncol(Q_eta)

  is_race <- identical(characteristic, "race")
  n_xi_free <- if (is_race) T_xi else (T_xi - 1L)
  n_eta_free <- T_eta - 1L
  n_free <- n_xi_free + n_eta_free

  ## --- penalty grid (or a single fixed node) -------------------------------
  if (!is.null(penalties)) {
    if (!is.list(penalties) || is.null(penalties$c_xi) || is.null(penalties$c_eta)) {
      .gp_2l_abort(
        "`penalties` must be a list(c_xi = , c_eta = ) when supplied.",
        class = "gradepath_validation_error"
      )
    }
    if (length(penalties$c_xi) != 1L || length(penalties$c_eta) != 1L ||
        !is.finite(penalties$c_xi) || !is.finite(penalties$c_eta)) {
      .gp_2l_abort(
        "`penalties` must contain finite scalar `c_xi` and `c_eta` values.",
        class = "gradepath_validation_error"
      )
    }
    grid <- data.frame(c_xi = as.numeric(penalties$c_xi),
                       c_eta = as.numeric(penalties$c_eta))
    grid_is_fixed <- TRUE
  } else {
    grid <- .gp_2l_penalty_grid(characteristic, control)
    grid_is_fixed <- FALSE
  }
  nP <- nrow(grid)

  ## --- multistart draws (seeded; precomputed per node) ---------------------
  start0 <- rep(0, n_free)
  extra_starts <- if (n_starts > 1L) {
    seed <- if (!is.null(control$seed)) as.integer(control$seed) else 1234L
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
    set.seed(seed)
    lapply(seq_len(nP), function(.) {
      matrix(stats::rnorm((n_starts - 1L) * n_free), nrow = n_starts - 1L,
             ncol = n_free)
    })
  } else NULL

  fit_one <- function(c_xi, c_eta, start) {
    tryCatch(
      .gp_2l_fit_node(c_xi, c_eta, v_hat, s_v, industry_idx,
                      Q_xi, supp_xi, Q_eta, supp_eta,
                      characteristic, start = start,
                      max_iter = max_iter, tol = opt_tol),
      error = function(e) NULL
    )
  }
  fit_node <- function(c_xi, c_eta, k) {
    f <- fit_one(c_xi, c_eta, start0)
    if (!is.null(f) && (any(!is.finite(f$g_xi)) || any(!is.finite(f$g_eta)))) f <- NULL
    if (!is.null(extra_starts)) {
      offs <- extra_starts[[k]]
      for (j in seq_len(nrow(offs))) {
        fj <- fit_one(c_xi, c_eta, start0 + offs[j, ])
        if (is.null(fj) || any(!is.finite(fj$g_xi)) || any(!is.finite(fj$g_eta))) next
        if (is.null(f) || (is.finite(fj$objective) && fj$objective < f$objective)) {
          f <- fj
        }
      }
    }
    f
  }

  ## --- loop over penalties, score by coupling J ----------------------------
  J_vec    <- rep(NA_real_, nP)
  obj_vec  <- rep(NA_real_, nP)
  conv_vec <- rep(FALSE, nP)
  conv_code_vec <- rep(NA_integer_, nP)
  conv_msg_vec  <- rep(NA_character_, nP)
  best <- list(J = Inf, idx = NA_integer_)

  for (k in seq_len(nP)) {
    c_xi <- grid$c_xi[k]; c_eta <- grid$c_eta[k]
    f <- fit_node(c_xi, c_eta, k)
    if (is.null(f) || any(!is.finite(f$g_xi)) || any(!is.finite(f$g_eta))) next
    m_model <- .gp_2l_model_moment_vector(supp_xi, f$g_xi, supp_eta, f$g_eta,
                                          characteristic)
    Jk <- .gp_2l_J(m_hat, V_m, m_model, V_m_inv = V_m_inv)
    J_vec[k]    <- Jk
    obj_vec[k]  <- f$objective
    conv_vec[k] <- f$converged
    conv_code_vec[k] <- f$convergence_code
    conv_msg_vec[k]  <- f$convergence_message
    if (is.finite(Jk) && Jk < best$J) {
      best <- list(J = Jk, idx = k, g_xi = f$g_xi, g_eta = f$g_eta,
                   alpha_xi_free = f$alpha_xi_free, alpha_eta_free = f$alpha_eta_free,
                   objective = f$objective, loglik = f$loglik, m_model = m_model,
                   converged = f$converged, convergence_code = f$convergence_code,
                   convergence_message = f$convergence_message)
    }
  }

  if (is.na(best$idx)) {
    .gp_2l_abort(
      "Two-level deconvolution failed: no penalty node produced a finite fit.",
      class = "gradepath_estimation_error"
    )
  }
  if (!isTRUE(best$converged)) {
    warning(sprintf(
      "Selected two-level deconvolution fit did not converge (optim code %s): %s",
      best$convergence_code,
      if (is.na(best$convergence_message)) "no optimizer message" else best$convergence_message
    ), call. = FALSE)
  }

  ## --- final densities: normalise to probability masses (g/sum(g)) ---------
  ## estimate_lsqnonlin.m:469/474: g_xi = g_xi/sum(g_xi); g_psi = g_psi/sum(g_psi).
  density <- best$g_xi / sum(best$g_xi)
  g_eta <- best$g_eta / sum(best$g_eta)
  prior_mean <- sum(supp_xi * density)

  ## --- assemble the gp_prior ----------------------------------------------
  new_prior      <- .gp_2l_get("new_gp_prior")
  validate_prior <- .gp_2l_get("validate_gp_prior")

  diagnostics <- list(
    method            = "native-two-level-spline-softmax",
    group_fx          = 1L,
    c_xi              = grid$c_xi[best$idx],
    c_eta             = grid$c_eta[best$idx],
    penalty_index     = best$idx,
    J                 = best$J,
    log_likelihood    = best$loglik,
    objective         = best$objective,
    converged         = best$converged,
    convergence_code  = best$convergence_code,
    convergence_message = best$convergence_message,
    penalty_grid      = grid,
    J_grid            = J_vec,
    objective_grid    = obj_vec,
    converged_grid    = conv_vec,
    convergence_code_grid = conv_code_vec,
    convergence_message_grid = conv_msg_vec,
    model_moments     = best$m_model,         # [mu_xi, sigma_xi, sigma_eta] race
    coupling_m_hat    = m_hat,
    support_eta       = supp_eta,
    g_eta             = g_eta,
    caps              = caps,
    n_knots           = n_knots,
    spline_dim_xi     = T_xi,
    spline_dim_eta    = T_eta,
    grid_is_fixed     = grid_is_fixed,
    ## theta-pushforward hook: the native two-level theta-pushforward is built separately
    ## (gp_pushforward_theta, get_g_theta.m two-support port); left NULL here.
    g_theta           = NULL
  )
  metadata <- list(
    characteristic = characteristic,
    beta           = fit$beta,
    n_carriers     = length(v_hat),
    n_industries   = max(industry_idx),
    supp_pts       = supp_pts,
    alpha_xi_free  = best$alpha_xi_free,
    alpha_eta_free = best$alpha_eta_free
  )

  prior <- new_prior(
    support     = supp_xi,
    density     = density,
    mean        = prior_mean,
    scale       = "r",
    diagnostics = diagnostics,
    metadata    = metadata
  )
  validate_prior(prior)
}
