# =============================================================================
# Tool:        reach.hydro — PDM capacity distributions
# Description: Five soil-moisture capacity distributions from Moore (2007)
#              Appendices A-E, plus the Newton-Raphson c* solver (Appendix F).
#
#              Appendix A: Pareto       — closed-form Smax and c*
#              Appendix B: Rectangular  — closed-form Smax and c*
#              Appendix C: Exponential  — closed-form Smax and c*
#              Appendix D: Triangular   — closed-form Smax; c* via N-R
#              Appendix E: Log-Normal   — closed-form Smax; c* via N-R (App F)
#
#              The normal distribution has been removed — it is not in any
#              appendix of Moore (2007).
#
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-03-16 - JP: removed normal (not in paper); replaced
#                               uniroot with Newton-Raphson (Appendix F);
#                               renamed uniform->rectangular per Appendix B.
# Tier:        1
# References:
#   Moore, R.J. (2007). The PDM rainfall-runoff model.
#   Hydrol. Earth Syst. Sci., 11, 483-499.
# =============================================================================

# =============================================================================
# APPENDIX A — Pareto distribution
# F(c) = 1 - [(cmax-c)/(cmax-cmin)]^{1/(1+b)},  c in [cmin, cmax]
# Closed-form Smax and c*.
# =============================================================================

.pareto_cdf <- function(c, cmin, cmax, b) {
  c <- pmin(pmax(c, cmin), cmax)
  1 - ((cmax - c) / (cmax - cmin))^(1 / (1 + b))
}

.pareto_smax <- function(cmin, cmax, b) {
  cmin + (cmax - cmin) * b / (1 + b)
}

.pareto_cstar <- function(S, cmin, cmax, b) {
  Smax <- .pareto_smax(cmin, cmax, b)
  S    <- pmax(0, pmin(S, Smax))
  cs   <- cmax - (cmax - cmin) * (1 - S / Smax)^(1 / (1 + b))
  pmin(pmax(cs, cmin), cmax)
}

.pareto_runoff <- function(P, S, cmin, cmax, b) {
  if (P <= 0) return(0)
  cs     <- .pareto_cstar(S, cmin, cmax, b)
  fsat   <- .pareto_cdf(cs, cmin, cmax, b)
  avail  <- (cmax - cs) / (1 + b)
  runoff <- P * fsat + pmax(0, P - avail) * (1 - fsat)
  pmin(pmax(runoff, 0), P)
}

# =============================================================================
# APPENDIX B — Rectangular distribution
# F(c) = (c - cmin) / (cmax - cmin),  c in [cmin, cmax]
# Smax = (cmax + cmin) / 2
# Closed-form c*: c* = cmin + sqrt(2 * S * (cmax - cmin))
# =============================================================================

.rectangular_cdf <- function(c, cmin, cmax) {
  c <- pmin(pmax(c, cmin), cmax)
  (c - cmin) / (cmax - cmin)
}

.rectangular_smax <- function(cmin, cmax) {
  (cmax + cmin) / 2
}

.rectangular_cstar <- function(S, cmin, cmax) {
  Smax <- .rectangular_smax(cmin, cmax)
  S    <- pmax(0, pmin(S, Smax))
  cs   <- cmin + sqrt(2 * S * (cmax - cmin))
  pmin(pmax(cs, cmin), cmax)
}

.rectangular_runoff <- function(P, S, cmin, cmax) {
  if (P <= 0) return(0)
  cs     <- .rectangular_cstar(S, cmin, cmax)
  fsat   <- .rectangular_cdf(cs, cmin, cmax)
  avail  <- (cmax - cs) / 2
  runoff <- P * fsat + pmax(0, P - avail) * (1 - fsat)
  pmin(pmax(runoff, 0), P)
}

# =============================================================================
# APPENDIX C — Exponential distribution
# F(c) = 1 - exp(-c / cmax),  c >= 0   (cmax = mean capacity)
# Smax = cmax
# Closed-form c*: c* = -cmax * ln(1 - S/cmax)
# =============================================================================

.exponential_cdf <- function(c, cmax) {
  c <- pmax(c, 0)
  1 - exp(-c / cmax)
}

.exponential_smax <- function(cmax) {
  cmax
}

.exponential_cstar <- function(S, cmax) {
  S <- pmax(0, pmin(S, cmax * 0.9999))
  pmax(-cmax * log(1 - S / cmax), 0)
}

.exponential_runoff <- function(P, S, cmax) {
  if (P <= 0) return(0)
  cs     <- .exponential_cstar(S, cmax)
  fsat   <- .exponential_cdf(cs, cmax)
  avail  <- cmax * exp(-cs / cmax)
  runoff <- P * fsat + pmax(0, P - avail) * (1 - fsat)
  pmin(pmax(runoff, 0), P)
}

# =============================================================================
# APPENDIX D — Triangular distribution
# Symmetric triangle on [cmin, cmax] with peak at midpoint cm.
#
# PDF:
#   f(c) = 4(c-cmin)/(cmax-cmin)^2       cmin <= c <= cm
#   f(c) = 4(cmax-c)/(cmax-cmin)^2       cm < c <= cmax
#
# CDF:
#   F(c) = 2[(c-cmin)/(cmax-cmin)]^2     cmin <= c <= cm
#   F(c) = 1 - 2[(cmax-c)/(cmax-cmin)]^2 cm < c <= cmax
#
# Smax = (cmax + cmin) / 2  (triangular mean = midpoint by symmetry)
#
# c* solved via Newton-Raphson (Appendix F) since the cubic equation
# for the lower half has no simple closed form for arbitrary S.
# =============================================================================

.triangular_cdf <- function(c, cmin, cmax) {
  c   <- pmin(pmax(c, cmin), cmax)
  cm  <- (cmin + cmax) / 2
  rng <- cmax - cmin
  ifelse(
    c <= cm,
    2 * ((c - cmin) / rng)^2,
    1 - 2 * ((cmax - c) / rng)^2
  )
}

.triangular_smax <- function(cmin, cmax) {
  (cmax + cmin) / 2
}

.triangular_cstar <- function(S, cmin, cmax) {
  f <- function(c) .triangular_cdf(c, cmin, cmax)
  .newton_raphson_cstar(S, f, cmin, cmax)
}

.triangular_runoff <- function(P, S, cmin, cmax) {
  if (P <= 0) return(0)
  f    <- function(c) .triangular_cdf(c, cmin, cmax)
  Smax <- .triangular_smax(cmin, cmax)
  .generic_runoff(P, S, f, cmin, cmax, Smax)
}

# =============================================================================
# APPENDIX E — Log-Normal distribution
# F(c) = Phi((ln(c) - mu_lnc) / sigma_lnc),  c > 0
# Smax = exp(mu_lnc + sigma_lnc^2 / 2)  (exact lognormal mean)
# c* solved via Newton-Raphson (Appendix F).
# =============================================================================

.lognormal_cdf <- function(c, mu_lnc, sigma_lnc) {
  c <- pmax(c, 1e-12)
  pnorm((log(c) - mu_lnc) / sigma_lnc)
}

.lognormal_smax <- function(mu_lnc, sigma_lnc) {
  exp(mu_lnc + sigma_lnc^2 / 2)
}

.lognormal_cstar <- function(S, mu_lnc, sigma_lnc) {
  chi <- exp(mu_lnc + 5 * sigma_lnc)
  f   <- function(c) .lognormal_cdf(c, mu_lnc, sigma_lnc)
  .newton_raphson_cstar(S, f, 0, chi)
}

.lognormal_runoff <- function(P, S, mu_lnc, sigma_lnc) {
  if (P <= 0) return(0)
  chi  <- exp(mu_lnc + 5 * sigma_lnc)
  f    <- function(c) .lognormal_cdf(c, mu_lnc, sigma_lnc)
  Smax <- .numerical_smax(f, 0, chi)
  .generic_runoff(P, S, f, 0, chi, Smax)
}

# =============================================================================
# PUBLIC DISPATCH FUNCTIONS
# =============================================================================

#' CDF F(c) for the chosen capacity distribution
#'
#' @param c      Storage capacity values \[mm\].
#' @param dist   One of `"pareto"`, `"rectangular"`, `"exponential"`,
#'               `"triangular"`, `"lognormal"`. These correspond directly to
#'               Appendices A-E of Moore (2007).
#' @param params A `PdmParams` object or compatible named list.
#' @return Numeric vector of probabilities in \[0, 1\].
#' @references Moore (2007), Appendices A-E.
#' @export
capacity_cdf <- function(c, dist, params) {
  switch(dist,
    pareto      = .pareto_cdf(c, params$cmin, params$cmax, params$b),
    rectangular = .rectangular_cdf(c, params$cmin, params$cmax),
    exponential = .exponential_cdf(c, params$cmax),
    triangular  = .triangular_cdf(c, params$cmin, params$cmax),
    lognormal   = .lognormal_cdf(c, params$mu_lnc, params$sigma_lnc),
    stop(paste("Unknown distribution:", dist,
               "\nValid options: pareto, rectangular, exponential,",
               "triangular, lognormal"), call. = FALSE)
  )
}

#' Maximum basin storage Smax for the chosen distribution
#'
#' @inheritParams capacity_cdf
#' @return Smax \[mm\].
#' @references Moore (2007), Appendices A-E.
#' @export
capacity_smax <- function(dist, params) {
  switch(dist,
    pareto      = .pareto_smax(params$cmin, params$cmax, params$b),
    rectangular = .rectangular_smax(params$cmin, params$cmax),
    exponential = .exponential_smax(params$cmax),
    triangular  = .triangular_smax(params$cmin, params$cmax),
    lognormal   = .lognormal_smax(params$mu_lnc, params$sigma_lnc),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}

#' Critical capacity c* given current basin moisture S
#'
#' @param S Current soil moisture \[mm\].
#' @inheritParams capacity_cdf
#' @return c* \[mm\].
#' @references Moore (2007), Eq. 3 and Appendices A-F.
#' @export
capacity_cstar <- function(S, dist, params) {
  switch(dist,
    pareto      = .pareto_cstar(S, params$cmin, params$cmax, params$b),
    rectangular = .rectangular_cstar(S, params$cmin, params$cmax),
    exponential = .exponential_cstar(S, params$cmax),
    triangular  = .triangular_cstar(S, params$cmin, params$cmax),
    lognormal   = .lognormal_cstar(S, params$mu_lnc, params$sigma_lnc),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}

#' Direct runoff from rainfall P given current soil moisture S
#'
#' @param P Rainfall depth this timestep \[mm\].
#' @param S Current soil moisture \[mm\].
#' @inheritParams capacity_cdf
#' @return Direct runoff \[mm\], bounded to \[0, P\].
#' @references Moore (2007), Eq. 1 and Appendices A-F.
#' @export
capacity_runoff <- function(P, S, dist, params) {
  switch(dist,
    pareto      = .pareto_runoff(P, S, params$cmin, params$cmax, params$b),
    rectangular = .rectangular_runoff(P, S, params$cmin, params$cmax),
    exponential = .exponential_runoff(P, S, params$cmax),
    triangular  = .triangular_runoff(P, S, params$cmin, params$cmax),
    lognormal   = .lognormal_runoff(P, S, params$mu_lnc, params$sigma_lnc),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}
