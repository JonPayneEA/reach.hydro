# =============================================================================
# Tool:        reach.hydro — flood_frequency.R
# Description: Thin compatibility shim. The full flood frequency implementation
#              is in feh_vol3_statistical.R (FEH Vol 3 statistical methods),
#              feh_vol4_rainfall.R (DDF / design storm), feh_refh2.R (ReFH2),
#              and feh_fsr_uh.R (FSR unit hydrograph method).
#
#              fit_gev() and fit_glo() are kept here as convenience wrappers
#              so that existing call sites continue to work.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: replaced stub with shim to FEH Vol 3 methods
# Tier:        2
# =============================================================================

#' Fit a GEV distribution to annual maxima (convenience wrapper)
#'
#' Wrapper around [feh_single_site()] with `dist = "gev"`. For the full
#' FEH-compliant analysis use [feh_single_site()] directly.
#'
#' @param annual_max     Numeric vector of annual maximum flows.
#' @param return_periods Return periods [reach.hydro::years].
#' @return A `FehStatFit` object. See [feh_single_site()].
#' @seealso [feh_single_site()], [fit_glo()]
#' @export
fit_gev <- function(annual_max,
                    return_periods = c(2, 5, 10, 20, 50, 100, 200)) {
  feh_single_site(annual_max, dist = "gev",
                  return_periods = return_periods, n_boot = 0L)
}

#' Fit a GLO distribution to annual maxima (convenience wrapper)
#'
#' Wrapper around [feh_single_site()] with `dist = "glo"` (FEH-recommended
#' for GB). For the full FEH-compliant analysis use [feh_single_site()].
#'
#' @inheritParams fit_gev
#' @seealso [feh_single_site()], [fit_gev()]
#' @export
fit_glo <- function(annual_max,
                    return_periods = c(2, 5, 10, 20, 50, 100, 200)) {
  feh_single_site(annual_max, dist = "glo",
                  return_periods = return_periods, n_boot = 0L)
}

#' Return period flows from a fitted distribution (convenience wrapper)
#'
#' @param fit            A `FehStatFit` from [feh_single_site()], [fit_gev()],
#'                       or [fit_glo()].
#' @param return_periods Numeric vector of return periods [reach.hydro::years].
#' @return A `data.table` with columns `return_period_yr` and `flow`.
#' @export
return_period_flow <- function(fit,
                               return_periods = c(2, 5, 10, 20, 50, 100, 200)) {
  if (!inherits(fit, "FehStatFit"))
    stop("return_period_flow: expected a FehStatFit object.", call. = FALSE)
  probs  <- 1 - 1 / return_periods
  params <- fit$params
  flows  <- .qfun(probs, params, fit$dist) /
            .qfun(0.5, params, fit$dist) * fit$qmed
  data.table::data.table(
    return_period_yr = return_periods,
    flow             = round(flows, 3)
  )
}
