# =============================================================================
# Tool:        reach.hydro — unit conversions
# Description: Conversion functions between common hydrological flow units.
#              All functions are vectorised and preserve NA values.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: initial skeleton
# Tier:        2
# Inputs:      Numeric flow vectors, catchment area [km2], timestep [hours]
# Outputs:     Numeric vectors in target units
# Dependencies: checkmate
# =============================================================================

#' Convert millimetres per timestep to cubic metres per second
#'
#' @param mm_per_ts Numeric vector of flow depth \[mm/timestep\].
#' @param area_km2  Catchment area \[km²\].
#' @param dt_hours  Timestep duration \[hours\]. Default 1 (hourly).
#'
#' @return Numeric vector of flow \[m³/s\].
#'
#' @examples
#' mm_to_m3s(5, area_km2 = 100, dt_hours = 1)
#'
#' @export
mm_to_m3s <- function(mm_per_ts, area_km2, dt_hours = 1) {
  checkmate::assert_numeric(mm_per_ts)
  checkmate::assert_number(area_km2, lower = 0)
  checkmate::assert_number(dt_hours, lower = 0)
  # mm/ts * km2 -> m3/s:  (mm/1000) * (km2 * 1e6) / (dt_hours * 3600)
  mm_per_ts * area_km2 * 1e6 / (1000 * dt_hours * 3600)
}

#' Convert cubic metres per second to millimetres per timestep
#'
#' @param m3s       Numeric vector of flow \[m³/s\].
#' @param area_km2  Catchment area \[km²\].
#' @param dt_hours  Timestep duration \[hours\]. Default 1 (hourly).
#'
#' @return Numeric vector of flow depth \[mm/timestep\].
#'
#' @export
m3s_to_mm <- function(m3s, area_km2, dt_hours = 1) {
  checkmate::assert_numeric(m3s)
  checkmate::assert_number(area_km2, lower = 0)
  checkmate::assert_number(dt_hours, lower = 0)
  m3s * 1000 * dt_hours * 3600 / (area_km2 * 1e6)
}

#' Alias: m3s_to_cumecs (identity — both are m³/s, provided for readability)
#'
#' @param m3s Numeric vector \[m³/s\].
#' @return Numeric vector \[cumecs = m³/s\], unchanged.
#' @export
m3s_to_cumecs <- function(m3s) {
  checkmate::assert_numeric(m3s)
  m3s
}
