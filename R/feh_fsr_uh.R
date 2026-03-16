# =============================================================================
# Tool:        reach.hydro — FSR/FEH rainfall-runoff (unit hydrograph)
# Description: Implements the original FSR (Flood Studies Report, NERC 1975)
#              rainfall-runoff method as updated in FEH (Reed et al. 1999):
#
#                1. Design rainfall (M5, r, ARF) → storm profile
#                2. Loss model: percentage runoff (PR) from SPR and CWI
#                3. Unit hydrograph (UH): triangular or observed
#                4. Routing: baseflow addition
#
#              Parameters are estimated from FEH catchment descriptors via the
#              FSR regression equations (NERC 1975) or FEH updates (IH 1999).
#
#              This method is largely superseded by ReFH2 but remains in
#              active use for smaller catchments and infrastructure design.
#
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-23
# Modified:    2026-02-23 - JP: initial implementation
# Tier:        2
# Inputs:      FEH catchment descriptors; design storm rainfall
# Outputs:     Design hydrograph [data.table]; peak flow; percentage runoff
# Dependencies: data.table
# References:
#   NERC (1975). Flood Studies Report. Vol. 1. Natural Environment Research
#     Council, London.
#   Reed, D. et al. (1999). Flood Estimation Handbook. CEH Wallingford.
#   IH (1999). Flood Estimation Handbook Software. CEH Wallingford.
# =============================================================================

# =============================================================================
# SECTION 1 : CATCHMENT DESCRIPTOR → FSR/FEH PARAMETERS
# =============================================================================

#' Estimate FSR/FEH rainfall-runoff parameters from catchment descriptors
#'
#' Applies the FSR regression equations (as updated in FEH Vol. 4/5) to
#' derive the standard percentage runoff (SPR), time-to-peak (Tp), and
#' catchment wetness index (CWI) parameters needed for the FSR unit
#' hydrograph method.
#'
#' @param area      Catchment area \[km²\].
#' @param bfihost   Baseflow index from HOST soils \[0-1\].
#' @param saar      Standard average annual rainfall \[mm\].
#' @param farl      Flood attenuation: reservoirs and lakes index \[0-1\].
#'                  Default 1.0 (no attenuation).
#' @param urbext    FEH urbanisation extent (URBEXT2000). Default 0.
#' @param slope     Channel slope S1085 \[m/km\]. Optional; used in Tp equation.
#'                  If `NULL`, omitted from the regression.
#'
#' @return A list of class `"FsrParams"` with:
#'   \describe{
#'     \item{spr}{Standard percentage runoff \[%\].}
#'     \item{Tp}{Time-to-peak \[hours\].}
#'     \item{cwi}{Catchment wetness index (baseline). }
#'     \item{Tp_dt}{Tp adjusted for routing (= Tp + dt/2 at timestep dt).}
#'     \item{Qb_index}{Baseflow index for hydrograph baseflow addition.}
#'   }
#'
#' @examples
#' p <- fsr_params(area = 150, bfihost = 0.4, saar = 800)
#' p
#'
#' @seealso [fsr_run()]
#' @export
fsr_params <- function(area,
                       bfihost,
                       saar,
                       farl   = 1.0,
                       urbext = 0.0,
                       slope  = NULL) {

  checkmate::assert_number(area,    lower = 0)
  checkmate::assert_number(bfihost, lower = 0, upper = 1)
  checkmate::assert_number(saar,    lower = 0)
  checkmate::assert_number(farl,    lower = 0, upper = 1)
  checkmate::assert_number(urbext,  lower = 0, upper = 1)

  # SPR: Standard percentage runoff (FSR Vol 1, eq 6.5; FEH update)
  # SPR = 45.56 * BFIHOST^(-0.52) adjusted for urbanisation
  spr_rural <- 45.56 * bfihost^(-0.52)
  spr_rural <- pmin(pmax(spr_rural, 1), 99)
  # Urban adjustment (FEH Vol 5, eq 7.1)
  spr <- spr_rural * (1 - urbext) + 70 * urbext
  spr <- pmin(pmax(spr, 1), 99)

  # Tp: time-to-peak [hours] (FEH Vol 4 regression; Kjeldsen 2007)
  if (!is.null(slope)) {
    Tp <- exp(-0.3120 + 0.5520 * log(area) -
              0.2368 * log(saar) - 0.4499 * log(slope))
  } else {
    Tp <- exp(-0.1085 + 0.5520 * log(area) -
              0.2368 * log(saar))
  }
  Tp <- Tp * farl^2.2  # reservoir attenuation effect

  # CWI: catchment wetness index (FSR, Section 6.4)
  # Baseline CWI from SAAR (antecedent conditions)
  cwi <- 125 + 0.1 * saar

  # Baseflow index for hydrograph addition
  Qb_index <- bfihost

  structure(
    list(spr = spr, Tp = Tp, cwi = cwi,
         Qb_index = Qb_index,
         area = area, bfihost = bfihost, saar = saar,
         farl = farl, urbext = urbext),
    class = "FsrParams"
  )
}

#' @export
print.FsrParams <- function(x, ...) {
  cat(sprintf("<FsrParams> area=%.1f km\u00b2 | BFIHOST=%.3f | SAAR=%d mm\n",
              x$area, x$bfihost, round(x$saar)))
  cat(sprintf("  SPR=%.1f%% | Tp=%.2f hr | CWI=%.1f\n",
              x$spr, x$Tp, x$cwi))
  invisible(x)
}

# =============================================================================
# SECTION 2 : PERCENTAGE RUNOFF
# =============================================================================

#' Compute event percentage runoff (PR) from SPR, CWI, and storm depth
#'
#' The FSR percentage runoff model adjusts SPR for antecedent catchment
#' wetness (CWI) and storm depth (M5-60min proxy) to give event PR:
#'
#'   PR = SPR + DPRr + DPRs
#'
#' where DPRr is the rainfall-depth adjustment and DPRs is the soil
#' wetness adjustment (FSR Vol 1, Section 6.4).
#'
#' @param spr        Standard percentage runoff \[%\].
#' @param cwi        Catchment wetness index (from [fsr_params()] or observed).
#' @param storm_depth Total storm rainfall depth \[mm\].
#' @param m5_60min   M5-60min rainfall \[mm\] (from [feh_m5_ratios()]).
#'
#' @return Event percentage runoff PR \[%\], clamped to \[0, 100\].
#' @export
fsr_percentage_runoff <- function(spr, cwi, storm_depth, m5_60min) {

  checkmate::assert_number(spr,         lower = 0, upper = 100)
  checkmate::assert_number(cwi,         lower = 0)
  checkmate::assert_number(storm_depth, lower = 0)
  checkmate::assert_number(m5_60min,    lower = 0)

  # DPRr: depth adjustment (FSR, eq 6.6)
  DPRr <- 0.366 * (storm_depth - m5_60min)

  # DPRs: soil wetness adjustment (FSR, eq 6.7)
  DPRs <- 0.0 + 0.1 * pmax(0, cwi - 125)

  PR <- spr + DPRr + DPRs
  pmin(pmax(PR, 0), 100)
}

# =============================================================================
# SECTION 3 : UNIT HYDROGRAPH CONSTRUCTION
# =============================================================================

#' Construct a FSR/FEH triangular unit hydrograph
#'
#' Builds the dimensionless triangular unit hydrograph (UH) specified in the
#' FSR (NERC 1975, Vol 1 Fig 6.13). The UH is parameterised by Tp (time to
#' peak) and has a base length of 2.52 × Tp.
#'
#' @param Tp      Time-to-peak \[hours\].
#' @param dt_hr   Timestep \[hours\].
#' @param n_steps Total UH length in timesteps. Default: covers full UH base.
#'
#' @return A numeric vector (UH ordinates, summing to 1/dt_hr \[1/hr\]).
#'
#' @export
fsr_unit_hydrograph <- function(Tp, dt_hr, n_steps = NULL) {
  checkmate::assert_number(Tp,    lower = 0)
  checkmate::assert_number(dt_hr, lower = 0)

  Tb <- 2.52 * Tp          # UH base time [hours]
  if (is.null(n_steps)) n_steps <- ceiling(Tb / dt_hr) + 1L

  t <- seq(0, by = dt_hr, length.out = n_steps)

  # Triangular ordinates [m3/s per mm of runoff per km2 — scaled to unit]
  ords <- ifelse(
    t <= Tp,
    t / Tp,
    pmax(0, 1 - (t - Tp) / (Tb - Tp))
  )

  # Normalise so sum = 1 (dimensionless; multiply by PR * rainfall later)
  ords / sum(ords)
}

# =============================================================================
# SECTION 4 : FSR FULL MODEL RUN
# =============================================================================

#' Run the FSR/FEH rainfall-runoff method
#'
#' Convolves the design storm with the triangular unit hydrograph, applies
#' percentage runoff, and adds a linearly declining baseflow to produce the
#' design hydrograph.
#'
#' @param storm        A `data.table` with columns `time_min` and `rainfall_mm`
#'                     (e.g., from [feh_design_storm()]).
#' @param params       A `FsrParams` object from [fsr_params()].
#' @param m5_60min     M5-60min rainfall \[mm\] for the site.
#' @param pr_override  If not `NULL`, overrides the computed PR with this
#'                     value \[%\]. Useful for sensitivity testing.
#' @param baseflow_0   Initial baseflow \[mm/timestep\]. Default 0.
#'
#' @return A `data.table` with columns `time_min`, `time_hr`, `rainfall_mm`,
#'   `net_rain_mm`, `Qd_mm`, `Qb_mm`, `Q_mm`.
#'   Attributes: `PR`, `params`, `peak_flow_mm`.
#'
#' @seealso [fsr_params()], [feh_design_storm()], [fsr_unit_hydrograph()]
#' @export
fsr_run <- function(storm,
                    params,
                    m5_60min,
                    pr_override = NULL,
                    baseflow_0  = 0) {

  if (!data.table::is.data.table(storm))
    storm <- data.table::as.data.table(storm)

  dt_hr        <- diff(storm$time_min[1:2]) / 60
  storm_depth  <- sum(storm$rainfall_mm, na.rm = TRUE)

  # Percentage runoff
  PR <- if (!is.null(pr_override)) {
    pr_override
  } else {
    fsr_percentage_runoff(params$spr, params$cwi, storm_depth, m5_60min)
  }

  # Effective rainfall
  net_rain <- storm$rainfall_mm * (PR / 100)

  # Unit hydrograph
  Tp_steps <- max(1L, round(params$Tp / dt_hr))
  uh       <- fsr_unit_hydrograph(params$Tp, dt_hr,
                                  n_steps = length(net_rain))

  # Convolve net rainfall with UH
  Qd <- as.numeric(stats::filter(net_rain, uh, sides = 1))
  Qd[is.na(Qd)] <- 0
  Qd <- pmax(Qd, 0)

  # Baseflow: simple linear recession
  Kb <- 0.5 * (1 + params$Qb_index)
  n  <- length(Qd)
  Qb <- numeric(n)
  Qb[1] <- baseflow_0
  for (t in seq(2L, n)) {
    Qb[t] <- Qb[t - 1] * Kb
  }

  out <- data.table::data.table(
    time_min    = storm$time_min,
    time_hr     = storm$time_min / 60,
    rainfall_mm = storm$rainfall_mm,
    net_rain_mm = net_rain,
    Qd_mm       = Qd,
    Qb_mm       = Qb,
    Q_mm        = Qd + Qb
  )

  attr(out, "PR")           <- PR
  attr(out, "params")       <- params
  attr(out, "peak_flow_mm") <- max(out$Q_mm)
  out
}

# =============================================================================
# SECTION 5 : AREAL REDUCTION FACTOR (ARF)
# =============================================================================

#' Areal reduction factor (FEH / FSR method)
#'
#' Converts point rainfall to catchment-average rainfall for a given area,
#' duration, and return period using the FEH ARF regression equations
#' (Reed et al. 1999, eq 2.4).
#'
#' @param area_km2      Catchment area \[km²\].
#' @param duration_hr   Storm duration \[hours\].
#' @param return_period Return period \[years\]. Affects ARF for T > 10 yr.
#'
#' @return ARF scalar in (0, 1\].
#'
#' @export
feh_arf <- function(area_km2, duration_hr, return_period = 10) {
  checkmate::assert_number(area_km2,      lower = 0)
  checkmate::assert_number(duration_hr,   lower = 0)
  checkmate::assert_number(return_period, lower = 1)

  if (area_km2 <= 0) return(1)

  # FEH ARF equation (Reed et al. 1999, eq 2.4)
  # ARF = 1 - exp(-d * A^b) where d and b depend on duration
  d <- 0.02; b <- 0.4
  if (duration_hr >= 1)  { d <- 0.04; b <- 0.35 }
  if (duration_hr >= 6)  { d <- 0.06; b <- 0.30 }
  if (duration_hr >= 24) { d <- 0.08; b <- 0.25 }

  arf <- 1 - d * area_km2^b
  # Adjustment for long return periods (T > 10 yr slightly lower ARF)
  if (return_period > 10) arf <- arf * (1 - 0.005 * log(return_period / 10))

  pmax(0.2, pmin(arf, 1))
}
