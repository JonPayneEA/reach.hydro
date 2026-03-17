# =============================================================================
# Tool:        reach.hydro — PDM numerical helpers
# Description: Newton-Raphson solver for c* as specified in Moore (2007)
#              Appendix F. Used when no closed-form inversion of
#              S(c*) = integral_clo^c* [1-F(c)] dc exists.
#
#              The Newton-Raphson iteration solves:
#                g(c*) = S(c*) - S = 0
#              where S(c*) = integral_clo^c* [1-F(c)] dc
#
#              By the fundamental theorem of calculus:
#                g'(c*) = dS/dc* = 1 - F(c*)
#
#              So the Newton-Raphson update is:
#                c*_{n+1} = c*_n - [S(c*_n) - S] / [1 - F(c*_n)]
#
#              The trapezoidal integral S(c*) is computed on an adaptive
#              grid. The iteration converges rapidly (typically 3-5 steps)
#              because g'(c*) = 1 - F(c*) >= 0 and the function is monotone.
#
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-03-16 - JP: replaced uniroot() with Newton-Raphson per
#                               Moore (2007) Appendix F.
# Tier:        1
# References:
#   Moore, R.J. (2007). The PDM rainfall-runoff model.
#   Hydrol. Earth Syst. Sci., 11, 483-499. Appendix F.
# =============================================================================

# Numerically compute Smax = integral_clo^chi [1 - F(c)] dc
# Using the trapezoidal rule on n_grid equally-spaced points.
.numerical_smax <- function(cdf_fn, clo, chi, n_grid = 500L) {
  cc   <- seq(clo, chi, length.out = n_grid)
  vals <- 1 - cdf_fn(cc)
  sum(diff(cc) * (vals[-n_grid] + vals[-1L]) / 2)
}

# Compute S(cs) = integral_clo^cs [1-F(c)] dc via trapezoidal rule
.numerical_s_of_cs <- function(cs, cdf_fn, clo, n = 200L) {
  if (cs <= clo) return(0)
  cc <- seq(clo, cs, length.out = n)
  sum(diff(cc) * ((1 - cdf_fn(cc[-n])) + (1 - cdf_fn(cc[-1L]))) / 2)
}

#' Newton-Raphson inversion of S(c*) = S (Moore 2007, Appendix F)
#'
#' Finds c* such that integral_clo^c*  \code{1-F(c)} dc = S using Newton-Raphson:
#'   c*_{n+1} = c*_n -  \code{S(c*_n) - S} /  \code{1 - F(c*_n)}
#'
#' @param S       Target soil moisture (mm); clamped to `[0, Smax]`.
#' @param cdf_fn  CDF function; signature f(c) -> `[0,1]`.
#' @param clo     Lower bound of capacity domain.
#' @param chi     Upper bound of capacity domain.
#' @param tol     Convergence tolerance. Default 1e-6.
#' @param max_iter Maximum Newton-Raphson iterations. Default 50.
#' @return c* (mm).
#' @keywords internal
.newton_raphson_cstar <- function(S, cdf_fn, clo, chi,
                                  tol = 1e-6, max_iter = 50L) {
  Smax <- .numerical_smax(cdf_fn, clo, chi)
  S    <- pmax(0, pmin(S, Smax))
  if (S <= 0)    return(clo)
  if (S >= Smax) return(chi)

  # Initial guess: linear interpolation
  cs <- clo + (chi - clo) * S / Smax

  for (i in seq_len(max_iter)) {
    S_cs  <- .numerical_s_of_cs(cs, cdf_fn, clo)
    resid <- S_cs - S
    if (abs(resid) < tol) break

    # Derivative: dS/dc* = 1 - F(c*)
    deriv <- 1 - cdf_fn(cs)

    # Avoid division by near-zero (at c* ~ chi, F(c*) ~ 1)
    if (abs(deriv) < 1e-12) break

    cs <- cs - resid / deriv
    cs <- pmax(clo, pmin(cs, chi))
  }

  cs
}

# Generic direct-runoff using Newton-Raphson c* inversion.
# Used for log-normal (the only distribution without a closed-form c*
# after removing the normal distribution).
.generic_runoff <- function(P, S, cdf_fn, clo, chi, Smax) {
  if (P <= 0) return(0)
  cs     <- .newton_raphson_cstar(S, cdf_fn, clo, chi)
  fsat   <- cdf_fn(cs)
  avail  <- Smax - S
  runoff <- P * fsat + pmax(0, P - avail * (1 - fsat)) * (1 - fsat)
  pmin(pmax(runoff, 0), P)
}
