# =============================================================================
# Tool:        reach.hydro package-level documentation
# Description: Package entry point, imports, and onLoad hook
# Flode Module: reach.hydro (formerly flode.hydro)
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: initial skeleton
# Tier:        1
# Inputs:      N/A (package-level)
# Outputs:     N/A (package-level)
# Dependencies: data.table, collapse, checkmate
# =============================================================================

#' reach.hydro: Core hydrological calculations for the REACH framework
#'
#' Provides the PDM rainfall-runoff model (Moore 2007) with five selectable
#' capacity distributions, calibration utilities, goodness-of-fit metrics,
#' flow statistics, unit conversions, catchment aggregation, and flood
#' frequency tools.
#'
#' @section PDM rainfall-runoff model:
#' The main entry point is [pdm()]. Parameters are validated via [pdm_params()]
#' and [pdm_validate_params()]. Five soil-moisture capacity distributions are
#' supported: `"pareto"`, `"rectangular"`, `"exponential"`, `"triangular"`,
#' `"lognormal"` (Moore 2007, Appendices A-E).
#'
#' @section Calibration:
#' [calibrate_pdm()] optimises parameters against observed flow using
#' Nelder-Mead. [compare_distributions()] runs all five distributions on the
#' same forcing data and returns a summary `data.table`.
#'
#' @section Goodness-of-fit:
#' [nse()], [kge()], [pbias()], [far()] — standard hydrological metrics.
#'
#' @section Flow statistics:
#' [flow_stats()], [baseflow_index()], [flow_duration_curve()],
#' [annual_maxima()], [peaks_over_threshold()].
#'
#' @section References:
#' Moore, R.J. (2007). The PDM rainfall-runoff model.
#' *Hydrology and Earth System Sciences*, 11, 483-499.
#' <https://doi.org/10.5194/hess-11-483-2007>
#'
#' @docType package
#' @name reach.hydro-package
"_PACKAGE"

# Suppress R CMD CHECK notes for data.table's non-standard evaluation
utils::globalVariables(c(".SD", ".N", ".I", ".GRP", ".", ":="))

.onLoad <- function(libname, pkgname) {
  # S7's print.S7_object calls str() directly and never consults S7's
  # internal method table, so S7::method(print, ...) assignments are
  # silently ignored. Registering via registerS3method() puts the method
  # in R's own S3 dispatch table where it is found before print.S7_object.
  registerS3method("print", "reach.hydro::FlodeReFHParams",
                    .print_FlodeReFHParams, envir = asNamespace(pkgname))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "reach.hydro ", utils::packageVersion("reach.hydro"),
    " | PDM rainfall-runoff model (Moore 2007) | ",
    "Forecasting and Warning Team"
  )
}
