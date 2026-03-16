# =============================================================================
# Tool:        reach.hydro — PDM numerical helpers
# Description: Numerical integration and root-finding utilities used when
#              closed-form c* inversions do not exist (glogistic, normal,
#              lognormal distributions). Not exported; used internally by the
#              distribution implementations in pdm_distributions.R.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: ported from PDM standalone script (Moore 2007)
# Tier:        1
# Inputs:      CDF functions, integration bounds, moisture state S
# Outputs:     Smax (scalar), c* (scalar)
# Dependencies: base R only
# =============================================================================

# Numerically compute Smax = integral_clo^chi [1 - F(c)] dc
# Using the trapezoidal rule over n_grid equally-spaced points.
#
# @param cdf_fn  CDF function; signature f(c) -> [0,1]
# @param clo     Lower bound of integration (usually cmin or 0)
# @param chi     Upper bound of integration (usually cmax or a large quantile)
# @param n_grid  Number of grid points for trapezoidal integration
# @return Smax [mm]
.numerical_smax <- function(cdf_fn, clo, chi, n_grid = 500L) {
  cc   <- seq(clo, chi, length.out = n_grid)
  vals <- 1 - cdf_fn(cc)
  sum(diff(cc) * (vals[-n_grid] + vals[-1L]) / 2)
}

# Numerically invert S(c*) = integral_clo^c* [1-F(c)] dc = S
# Uses uniroot() on the cumulative trapezoidal integral.
# Falls back to chi (fully saturated) on any numerical failure.
#
# @param S       Target soil moisture [mm]; clamped to [0, Smax]
# @param cdf_fn  CDF function
# @param clo     Lower bound
# @param chi     Upper bound
# @param tol     Root-finding tolerance (default 1e-6)
# @return c* [mm]
.numerical_cstar <- function(S, cdf_fn, clo, chi, tol = 1e-6) {
  Smax <- .numerical_smax(cdf_fn, clo, chi)
  S    <- pmax(0, pmin(S, Smax))
  if (S <= 0)    return(clo)
  if (S >= Smax) return(chi)

  # Cumulative trapezoidal integral from clo to cs
  .S_of_cs <- function(cs) {
    n  <- 200L
    cc <- seq(clo, cs, length.out = n)
    sum(diff(cc) * ((1 - cdf_fn(cc[-n])) + (1 - cdf_fn(cc[-1L]))) / 2)
  }

  tryCatch(
    uniroot(function(cs) .S_of_cs(cs) - S,
            interval = c(clo, chi), tol = tol)$root,
    error = function(e) chi
  )
}

# Generic direct-runoff for any CDF (used when no closed-form exists).
# Saturation-excess runoff: fraction fsat already saturated contributes
# all rainfall as runoff; unsaturated fraction contributes excess above
# remaining mean capacity.
#
# @param P     Rainfall depth this timestep [mm]
# @param S     Current soil moisture [mm]
# @param cdf_fn  CDF function
# @param clo   Lower bound
# @param chi   Upper bound
# @param Smax  Pre-computed Smax [mm]
# @return Direct runoff [mm], clamped to [0, P]
.generic_runoff <- function(P, S, cdf_fn, clo, chi, Smax) {
  if (P <= 0) return(0)
  cs     <- .numerical_cstar(S, cdf_fn, clo, chi)
  fsat   <- cdf_fn(cs)
  avail  <- Smax - S
  runoff <- P * fsat + pmax(0, P - avail * (1 - fsat)) * (1 - fsat)
  pmin(pmax(runoff, 0), P)
}
