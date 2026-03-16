# =============================================================================
# Tool:        reach.hydro — catchment aggregation
# Description: Areal rainfall computation and weighted spatial aggregation
#              of catchment data. Uses data.table for all tabular operations.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: initial skeleton
# Tier:        2
# Inputs:      data.table of gauge/gridded values with area weights
# Outputs:     Numeric vector or data.table of spatially aggregated values
# Dependencies: data.table, collapse
# TODO:        Thiessen polygon weight computation; kriging interpolation
# =============================================================================

#' Compute areal rainfall from gauge network
#'
#' Calculates catchment-average rainfall from a set of point gauges using
#' pre-computed weights (e.g., Thiessen polygon areas, kriging weights).
#'
#' @param dt       A `data.table` with one column per gauge and one row per
#'                 timestep, containing rainfall depths \[mm\].
#' @param weights  Named numeric vector of weights summing to 1.0. Names must
#'                 match column names of `dt`. If `NULL`, equal weights are
#'                 applied to all columns.
#' @param na.rm    Logical. If `TRUE`, missing gauges are excluded from the
#'                 weighted mean (weights rescaled accordingly). Default `TRUE`.
#'
#' @return Numeric vector of catchment-average rainfall \[mm\], length = `nrow(dt)`.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' rain_dt <- data.table(G1 = c(2.1, 3.0, 0.0), G2 = c(1.8, NA, 0.5))
#' areal_rainfall(rain_dt, weights = c(G1 = 0.6, G2 = 0.4))
#' }
#'
#' @export
areal_rainfall <- function(dt, weights = NULL, na.rm = TRUE) {
  if (!data.table::is.data.table(dt))
    dt <- data.table::as.data.table(dt)

  gauges <- names(dt)
  n_ts   <- nrow(dt)

  if (is.null(weights)) {
    weights <- stats::setNames(rep(1 / length(gauges), length(gauges)), gauges)
  } else {
    checkmate::assert_numeric(weights, lower = 0, any.missing = FALSE)
    missing_g <- setdiff(names(weights), gauges)
    if (length(missing_g) > 0)
      stop(sprintf("areal_rainfall: weight names not in dt: %s",
                   paste(missing_g, collapse = ", ")), call. = FALSE)
    # Normalise weights in case they don't sum to exactly 1
    weights <- weights / sum(weights)
  }

  result <- numeric(n_ts)
  for (i in seq_len(n_ts)) {
    row_vals  <- unlist(dt[i, .SD, .SDcols = names(weights)])
    if (na.rm && any(is.na(row_vals))) {
      valid_w   <- weights[!is.na(row_vals)]
      valid_w   <- valid_w / sum(valid_w)
      result[i] <- sum(row_vals[!is.na(row_vals)] * valid_w)
    } else {
      result[i] <- sum(row_vals * weights, na.rm = FALSE)
    }
  }
  result
}

#' Weighted spatial mean of a catchment variable
#'
#' General-purpose weighted mean for any spatially distributed variable
#' (e.g., temperature, PET, land-cover fraction) across sub-catchments.
#'
#' @param values  Numeric vector of sub-catchment values.
#' @param weights Numeric vector of weights (e.g., sub-catchment areas).
#'               Does not need to sum to 1; normalised internally.
#' @param na.rm   Logical. Default `TRUE`.
#'
#' @return Weighted mean scalar.
#'
#' @export
catchment_weighted_mean <- function(values, weights, na.rm = TRUE) {
  checkmate::assert_numeric(values)
  checkmate::assert_numeric(weights, lower = 0)
  if (length(values) != length(weights))
    stop("catchment_weighted_mean: `values` and `weights` must be the same length.",
         call. = FALSE)
  if (na.rm) {
    v <- !is.na(values) & !is.na(weights)
    values  <- values[v]
    weights <- weights[v]
  }
  if (sum(weights) == 0) return(NA_real_)
  sum(values * weights) / sum(weights)
}
