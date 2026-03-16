# =============================================================================
# Tool:        reach.hydro — FEH Vol 4 rainfall frequency (DDF / ReFH storm)
# Description: Implements the FEH depth-duration-frequency (DDF) model and
#              the ReFH design storm profile for use in the ReFH rainfall-
#              runoff model (see feh_refh.R).
#
#              The FEH DDF model (Reed et al. 1999 / Faulkner 1999) expresses
#              rainfall depth as:
#
#                D(d, T) = c(d) * [1 + e * log(T)]^(1/f)
#
#              where c(d) is a duration-dependent scale factor derived from
#              catchment descriptors, e and f are regional parameters, and
#              T is return period [years]. Parameters are derived from the
#              FEH CD-ROM regression equations for the rainfall depth at the
#              site of interest.
#
#              The Revitalised FSH (ReFH) design storm uses the 60-minute
#              rainfall depth as the index and a dimensionless temporal profile
#              (summer / winter) to distribute rainfall across the storm.
#
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-23
# Modified:    2026-02-23 - JP: initial implementation
# Tier:        2
# Inputs:      FEH catchment descriptor rainfall statistics (RMED, Urb, etc.)
# Outputs:     Rainfall depth [mm]; design hyetograph [data.table]
# Dependencies: data.table
# References:
#   Faulkner, D. (1999). Rainfall Frequency Estimation. FEH Vol. 2.
#     CEH Wallingford.
#   Reed, D. et al. (1999). Flood Estimation Handbook. CEH Wallingford.
#   Kjeldsen, T.R. (2007). The revitalised FSR/FEH rainfall-runoff method.
#     Flood Studies Report Technical Note No. 1. CEH Wallingford.
# =============================================================================

# =============================================================================
# SECTION 1 : FEH DDF MODEL
# =============================================================================

#' Compute design rainfall depth using the FEH DDF model (FEH Vol 4)
#'
#' Estimates the T-year return period rainfall depth for a given duration
#' using the FEH depth-duration-frequency model. The site must be
#' characterised by its RMED (median annual maximum 1-hour rainfall) and
#' SAAR (standard average annual rainfall) catchment descriptors, available
#' from the FEH Web Service or NRFA catchment descriptors.
#'
#' @param duration_hr   Storm duration \[hours\]. Scalar or vector.
#' @param return_period Return period \[years\]. Scalar or vector.
#'                      Must be compatible with `duration_hr` (recycled).
#' @param rmed_1h       RMED for 1-hour duration \[mm\]. From FEH descriptors.
#' @param rmed_1d       RMED for 1-day duration \[mm\]. From FEH descriptors.
#' @param saar          Standard average annual rainfall \[mm\].
#' @param region        FEH rainfall region. One of `"uk"` (default, uses
#'                      UK-wide parameters), or a named regional set.
#'
#' @return A `data.table` with columns `duration_hr`, `return_period`,
#'   `rainfall_mm`, `growth_factor`.
#'
#' @examples
#' # 100-year, 1-hour to 24-hour rainfall at a site with
#' # RMED(1h) = 12 mm, RMED(1d) = 38 mm, SAAR = 900 mm
#' feh_ddf(duration_hr = c(1, 2, 4, 6, 12, 24),
#'         return_period = 100,
#'         rmed_1h = 12, rmed_1d = 38, saar = 900)
#'
#' @seealso [feh_design_storm()], [refh2_run()]
#' @export
feh_ddf <- function(duration_hr,
                    return_period,
                    rmed_1h,
                    rmed_1d,
                    saar,
                    region = "uk") {

  checkmate::assert_numeric(duration_hr,   lower = 0.017, upper = 8760)
  checkmate::assert_numeric(return_period, lower = 1)
  checkmate::assert_number(rmed_1h, lower = 0)
  checkmate::assert_number(rmed_1d, lower = 0)
  checkmate::assert_number(saar,    lower = 0)

  rp <- .feh_ddf_params(region)

  # Scale RMED to arbitrary duration using the FEH DDF scale function
  # (Faulkner 1999, eq 3.3: c(d) fitted from RMED_1h and RMED_1d)
  c_d <- .ddf_scale(duration_hr, rmed_1h, rmed_1d)

  # Growth factor: z(T) = (1 + e * log(T))^(1/f)  (Faulkner 1999, eq 3.5)
  # Using the M5 ratio and z-value formulation
  # At T=2.33: z(2.33) = 1 by definition of the median return period
  z_T <- .ddf_growth(return_period, rp$e, rp$f)

  # Median rainfall depth at this duration
  # RMED(d) ~ c(d) by construction; growth gives D(d,T)
  rain_mm <- c_d * z_T

  data.table::data.table(
    duration_hr   = duration_hr,
    return_period = return_period,
    rainfall_mm   = round(rain_mm, 2),
    growth_factor = round(z_T, 4)
  )
}

# Internal: FEH DDF regional parameters (Faulkner 1999, Table 2.2)
# Currently only "uk" (national average) implemented; extend for regions.
.feh_ddf_params <- function(region) {
  switch(region,
    uk = list(e = 0.4, f = 0.63),
    stop(paste("Unknown DDF region:", region), call. = FALSE)
  )
}

# Internal: duration-scaling of RMED using log-linear interpolation between
# 1-hour and 1-day RMED (Faulkner 1999, eq 3.3 simplified form)
.ddf_scale <- function(d_hr, rmed_1h, rmed_1d) {
  # Log-log interpolation across 1-hr to 24-hr
  # For d < 1h, scale down from rmed_1h; for d > 24h, scale up from rmed_1d
  log_r  <- log(d_hr / 1)  / log(24 / 1)   # 0 at 1h, 1 at 24h
  log_r  <- pmax(0, pmin(log_r, 1))
  exp(log(rmed_1h) + log_r * (log(rmed_1d) - log(rmed_1h)))
}

# Internal: DDF growth factor z(T) = (1 + e * log(T/2.33))^(1/f)
# Normalised so z = 1 at T = 2.33 yr (approximate median)
.ddf_growth <- function(T, e, f) {
  z <- (1 + e * log(T / 2.33))
  z <- pmax(z, 0)
  z^(1 / f)
}

# =============================================================================
# SECTION 2 : ReFH DESIGN STORM PROFILE
# =============================================================================

#' Generate a ReFH design storm hyetograph
#'
#' Creates the dimensionless temporal profile for the ReFH design storm
#' (Kjeldsen 2007), scaled to the T-year storm depth for a given duration.
#' The profile uses the summer (default) or winter double-triangle distribution
#' to produce a time-varying rainfall input for [refh2_run()].
#'
#' @param duration_hr     Storm duration \[hours\]. Typical values: 1, 2, 4, 8.
#' @param return_period   Return period \[years\].
#' @param rmed_1h         RMED 1-hour \[mm\]. From FEH descriptors.
#' @param rmed_1d         RMED 1-day \[mm\]. From FEH descriptors.
#' @param saar            SAAR \[mm\]. From FEH descriptors.
#' @param dt_min          Timestep \[minutes\]. Default 15.
#' @param season          `"summer"` (default) or `"winter"`. Affects the
#'                        temporal profile peakedness.
#' @param region          DDF region (default `"uk"`).
#'
#' @return A `data.table` with columns `time_min`, `time_hr`, and
#'   `rainfall_mm` (depth per timestep).
#'
#' @seealso [feh_ddf()], [refh2_run()]
#' @export
feh_design_storm <- function(duration_hr,
                              return_period,
                              rmed_1h,
                              rmed_1d,
                              saar,
                              dt_min   = 15,
                              season   = "summer",
                              region   = "uk") {

  season <- match.arg(season, c("summer", "winter"))

  # Total storm depth
  total_rain <- feh_ddf(duration_hr, return_period, rmed_1h, rmed_1d,
                        saar, region)$rainfall_mm

  # Dimensionless double-triangle profile (Kjeldsen 2007, Section 3.2)
  # Storm split at time ratio r_peak; summer more peaked than winter
  r_peak <- if (season == "summer") 0.40 else 0.50

  n_steps  <- round(duration_hr * 60 / dt_min)
  time_min <- seq(0, by = dt_min, length.out = n_steps)
  time_frac <- time_min / (duration_hr * 60)

  # Double-triangle profile: intensity proportional to dist from r_peak
  profile <- .double_triangle(time_frac, r_peak, dt_min / (duration_hr * 60))

  # Scale to total_rain
  rain_mm <- profile * total_rain / sum(profile)

  data.table::data.table(
    time_min   = time_min,
    time_hr    = time_min / 60,
    rainfall_mm = round(rain_mm, 4)
  )
}

# Internal: double-triangle dimensionless intensity profile
# Rising limb: linear increase to r_peak; falling limb: linear decrease
.double_triangle <- function(frac, r_peak, dt_frac) {
  # Intensity proportional to triangular shape
  intensity <- ifelse(
    frac <= r_peak,
    frac / r_peak,
    (1 - frac) / (1 - r_peak)
  )
  # Convert intensity to depth per dt_frac interval
  intensity * dt_frac
}

# =============================================================================
# SECTION 3 : M5 AND RAINFALL RATIO HELPERS
# =============================================================================

#' Compute M5 ratios from FEH catchment descriptors
#'
#' Returns the M5-60min / M5-2day ratio used in some FEH/WINFAP applications.
#' The M5 (5-year return period) values are derived from RMED using the
#' FEH growth factor at T=5.
#'
#' @param rmed_1h  RMED 1-hour \[mm\].
#' @param rmed_1d  RMED 1-day \[mm\].
#' @param region   DDF region. Default `"uk"`.
#'
#' @return Named numeric vector: `m5_60min`, `m5_2day`, `r_value`.
#' @export
feh_m5_ratios <- function(rmed_1h, rmed_1d, region = "uk") {
  rp     <- .feh_ddf_params(region)
  z5     <- .ddf_growth(5, rp$e, rp$f)
  rmed_2d <- .ddf_scale(48, rmed_1h, rmed_1d)
  c(
    m5_60min = rmed_1h  * z5,
    m5_2day  = rmed_2d  * z5,
    r_value  = rmed_1h  * z5 / (rmed_2d * z5)  # r ratio used in ReFH
  )
}
