# =============================================================================
# Tool:        reach.hydro — PDM capacity distributions
# Description: Six soil-moisture capacity distributions from Moore (2007).
#              Each provides: _cdf, _smax, _cstar, _runoff implementations.
#              Sections 2 and 3 of the original PDM script.
#              Internal (_) functions are unexported. Public dispatch
#              functions (capacity_*) are exported and documented below.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: ported from PDM standalone script
# Tier:        1
# Inputs:      Capacity c [mm], moisture S [mm], parameter lists
# Outputs:     F(c), Smax, c*, direct runoff [mm]
# Dependencies: pdm_numerical.R (internal helpers)
# =============================================================================

# =============================================================================
# APPENDIX A — Power / Pareto distribution (Moore 2007, Appendix A)
# F(c) = 1 - [(cmax - c) / (cmax - cmin)]^{1/(1+b)},  c in [cmin, cmax]
# Closed-form Smax and c* available.
# =============================================================================

.pareto_cdf <- function(c, cmin, cmax, b) {
  c <- pmin(pmax(c, cmin), cmax)
  1 - ((cmax - c) / (cmax - cmin))^(1 / (1 + b))
}

.pareto_smax <- function(cmin, cmax, b) {
  # Integral of [1-F(c)] from cmin to cmax (Moore 2007 Appendix A)
  cmin + (cmax - cmin) * b / (1 + b)
}

.pareto_cstar <- function(S, cmin, cmax, b) {
  Smax <- .pareto_smax(cmin, cmax, b)
  S    <- pmax(0, pmin(S, Smax))
  cs   <- cmax - (cmax - cmin) * (1 - S / (cmax - cmin))^(1 / (1 + b))
  pmin(pmax(cs, cmin), cmax)
}

.pareto_runoff <- function(P, S, cmin, cmax, b) {
  if (P <= 0) return(0)
  cs     <- .pareto_cstar(S, cmin, cmax, b)
  fsat   <- .pareto_cdf(cs, cmin, cmax, b)
  avail  <- (cmax - cs) / (1 + b)  # mean residual capacity on unsaturated part
  runoff <- P * fsat + pmax(0, P - avail) * (1 - fsat)
  pmin(pmax(runoff, 0), P)
}

# =============================================================================
# APPENDIX B — Uniform distribution (Moore 2007, Appendix B)
# F(c) = (c - cmin) / (cmax - cmin),  c in [cmin, cmax]
# Smax = (cmax + cmin) / 2
# Closed-form c*: c* = cmin + sqrt(2 * S * (cmax - cmin))
# =============================================================================

.uniform_cdf <- function(c, cmin, cmax) {
  c <- pmin(pmax(c, cmin), cmax)
  (c - cmin) / (cmax - cmin)
}

.uniform_smax <- function(cmin, cmax) {
  (cmax + cmin) / 2
}

.uniform_cstar <- function(S, cmin, cmax) {
  Smax <- .uniform_smax(cmin, cmax)
  S    <- pmax(0, pmin(S, Smax))
  cs   <- cmin + sqrt(2 * S * (cmax - cmin))
  pmin(pmax(cs, cmin), cmax)
}

.uniform_runoff <- function(P, S, cmin, cmax) {
  if (P <= 0) return(0)
  cs     <- .uniform_cstar(S, cmin, cmax)
  fsat   <- .uniform_cdf(cs, cmin, cmax)
  avail  <- (cmax - cs) / 2
  runoff <- P * fsat + pmax(0, P - avail) * (1 - fsat)
  pmin(pmax(runoff, 0), P)
}

# =============================================================================
# APPENDIX C — Exponential distribution (Moore 2007, Appendix C)
# F(c) = 1 - exp(-c / cmax),  c >= 0   (cmax = mean capacity = 1/lambda)
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
  Smax <- .exponential_smax(cmax)
  S    <- pmax(0, pmin(S, Smax * 0.9999))
  pmax(-cmax * log(1 - S / cmax), 0)
}

.exponential_runoff <- function(P, S, cmax) {
  if (P <= 0) return(0)
  cs     <- .exponential_cstar(S, cmax)
  fsat   <- .exponential_cdf(cs, cmax)
  avail  <- cmax * exp(-cs / cmax)  # mean residual on unsaturated fraction
  runoff <- P * fsat + pmax(0, P - avail) * (1 - fsat)
  pmin(pmax(runoff, 0), P)
}

# =============================================================================
# APPENDIX D — Generalised Logistic distribution (Moore 2007, Appendix D)
# Bounded logistic sigmoid rescaled to [cmin, cmax].
# b -> 0: approaches uniform; b large: sharp step at midpoint.
# No closed-form c* — uses numerical inversion.
# =============================================================================

.glogistic_cdf <- function(c, cmin, cmax, b) {
  cbar    <- (cmax + cmin) / 2
  b_sc    <- b / ((cmax - cmin) / 10)
  sigmoid <- function(x) 1 / (1 + exp(-x))
  raw     <- sigmoid(b_sc * (c - cbar))
  Flo     <- sigmoid(b_sc * (cmin - cbar))
  Fhi     <- sigmoid(b_sc * (cmax - cbar))
  pmin(pmax((raw - Flo) / (Fhi - Flo), 0), 1)
}

.glogistic_smax <- function(cmin, cmax, b) {
  f <- function(c) .glogistic_cdf(c, cmin, cmax, b)
  .numerical_smax(f, cmin, cmax)
}

.glogistic_cstar <- function(S, cmin, cmax, b) {
  f <- function(c) .glogistic_cdf(c, cmin, cmax, b)
  .numerical_cstar(S, f, cmin, cmax)
}

.glogistic_runoff <- function(P, S, cmin, cmax, b) {
  if (P <= 0) return(0)
  f    <- function(c) .glogistic_cdf(c, cmin, cmax, b)
  Smax <- .numerical_smax(f, cmin, cmax)
  .generic_runoff(P, S, f, cmin, cmax, Smax)
}

# =============================================================================
# APPENDIX E — Normal distribution (Moore 2007, Appendix E)
# Truncated to c >= 0; renormalised.
# No closed-form c* — uses numerical inversion.
# =============================================================================

.normal_cdf <- function(c, mu_c, sigma_c) {
  Flo <- pnorm(0, mu_c, sigma_c)
  raw <- pnorm(c, mu_c, sigma_c)
  pmin(pmax((raw - Flo) / (1 - Flo), 0), 1)
}

.normal_smax <- function(mu_c, sigma_c) {
  chi <- mu_c + 6 * sigma_c
  f   <- function(c) .normal_cdf(c, mu_c, sigma_c)
  .numerical_smax(f, 0, chi)
}

.normal_cstar <- function(S, mu_c, sigma_c) {
  chi <- mu_c + 6 * sigma_c
  f   <- function(c) .normal_cdf(c, mu_c, sigma_c)
  .numerical_cstar(S, f, 0, chi)
}

.normal_runoff <- function(P, S, mu_c, sigma_c) {
  if (P <= 0) return(0)
  chi  <- mu_c + 6 * sigma_c
  f    <- function(c) .normal_cdf(c, mu_c, sigma_c)
  Smax <- .numerical_smax(f, 0, chi)
  .generic_runoff(P, S, f, 0, chi, Smax)
}

# =============================================================================
# APPENDIX F — Log-Normal distribution (Moore 2007, Appendix F)
# F(c) = Phi((ln(c) - mu_lnc) / sigma_lnc),  c > 0
# Smax = exp(mu_lnc + sigma_lnc^2 / 2)  (exact lognormal mean)
# No closed-form c* — uses numerical inversion.
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
  .numerical_cstar(S, f, 0, chi)
}

.lognormal_runoff <- function(P, S, mu_lnc, sigma_lnc) {
  if (P <= 0) return(0)
  chi  <- exp(mu_lnc + 5 * sigma_lnc)
  f    <- function(c) .lognormal_cdf(c, mu_lnc, sigma_lnc)
  Smax <- .numerical_smax(f, 0, chi)
  .generic_runoff(P, S, f, 0, chi, Smax)
}

# =============================================================================
# SECTION 3 — Unified public dispatch functions
# Accept dist name (string) + PdmParams list; route to implementation above.
# =============================================================================

#' CDF F(c) for the chosen capacity distribution
#'
#' @param c      Storage capacity values \[mm\]
#' @param dist   Distribution name. One of `"pareto"`, `"uniform"`,
#'               `"exponential"`, `"glogistic"`, `"normal"`, `"lognormal"`.
#' @param params A `PdmParams` object or compatible named list.
#' @return Numeric vector of probabilities in \[0, 1\].
#' @references Moore (2007), Appendices A-F.
#' @export
capacity_cdf <- function(c, dist, params) {
  switch(dist,
    pareto      = .pareto_cdf(c, params$cmin, params$cmax, params$b),
    uniform     = .uniform_cdf(c, params$cmin, params$cmax),
    exponential = .exponential_cdf(c, params$cmax),
    glogistic   = .glogistic_cdf(c, params$cmin, params$cmax, params$b),
    normal      = .normal_cdf(c, params$mu_c, params$sigma_c),
    lognormal   = .lognormal_cdf(c, params$mu_lnc, params$sigma_lnc),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}

#' Maximum basin storage Smax for the chosen distribution
#'
#' @inheritParams capacity_cdf
#' @return Smax \[mm\].
#' @references Moore (2007), Appendices A-F.
#' @export
capacity_smax <- function(dist, params) {
  switch(dist,
    pareto      = .pareto_smax(params$cmin, params$cmax, params$b),
    uniform     = .uniform_smax(params$cmin, params$cmax),
    exponential = .exponential_smax(params$cmax),
    glogistic   = .glogistic_smax(params$cmin, params$cmax, params$b),
    normal      = .normal_smax(params$mu_c, params$sigma_c),
    lognormal   = .lognormal_smax(params$mu_lnc, params$sigma_lnc),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}

#' Critical capacity c* given current basin moisture S
#'
#' @param S      Current soil moisture \[mm\].
#' @inheritParams capacity_cdf
#' @return c* \[mm\].
#' @references Moore (2007), Eq. 3 and Appendices A-F.
#' @export
capacity_cstar <- function(S, dist, params) {
  switch(dist,
    pareto      = .pareto_cstar(S, params$cmin, params$cmax, params$b),
    uniform     = .uniform_cstar(S, params$cmin, params$cmax),
    exponential = .exponential_cstar(S, params$cmax),
    glogistic   = .glogistic_cstar(S, params$cmin, params$cmax, params$b),
    normal      = .normal_cstar(S, params$mu_c, params$sigma_c),
    lognormal   = .lognormal_cstar(S, params$mu_lnc, params$sigma_lnc),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}

#' Direct runoff from rainfall P given current soil moisture S
#'
#' @param P      Rainfall depth this timestep \[mm\].
#' @param S      Current soil moisture \[mm\].
#' @inheritParams capacity_cdf
#' @return Direct runoff \[mm\], bounded to \[0, P\].
#' @references Moore (2007), Eq. 1 and Appendices A-F.
#' @export
capacity_runoff <- function(P, S, dist, params) {
  switch(dist,
    pareto      = .pareto_runoff(P, S, params$cmin, params$cmax, params$b),
    uniform     = .uniform_runoff(P, S, params$cmin, params$cmax),
    exponential = .exponential_runoff(P, S, params$cmax),
    glogistic   = .glogistic_runoff(P, S, params$cmin, params$cmax, params$b),
    normal      = .normal_runoff(P, S, params$mu_c, params$sigma_c),
    lognormal   = .lognormal_runoff(P, S, params$mu_lnc, params$sigma_lnc),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}
