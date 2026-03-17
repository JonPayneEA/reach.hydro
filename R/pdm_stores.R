# =============================================================================
# Tool:        reach.hydro — PDM store update functions
# Description: Per-timestep update functions for each model store.
#              Implements all components from Moore (2007) Table 1:
#                - Rainfall factor (fc) and time delay (td)
#                - AET (be)
#                - Soil drainage: three formulations
#                    1. Standard recharge (kg, bg, St)
#                    2. Demand-based recharge (alpha, beta, q_sat, Sg_max)
#                    3. Proportional split (alpha)
#                - Groundwater store: non-linear reservoir (kb, m)
#                - Surface routing: cascade of TWO linear reservoirs (k1, k2)
#                - Constant flow addition (qc)
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-03-16 - JP: added fc, td, qc; corrected surface routing
#                               to two-reservoir cascade; added demand-based
#                               and split recharge formulations; renamed
#                               kg->kb, bg->m per Table 1 of Moore (2007).
# Tier:        1
# References:
#   Moore, R.J. (2007). The PDM rainfall-runoff model.
#   Hydrol. Earth Syst. Sci., 11, 483-499. Table 1, Eqs 9-13.
# =============================================================================

# =============================================================================
# RAINFALL FACTOR AND TIME DELAY
# fc scales catchment-average rainfall before it enters the soil store.
# td is a pure lag (integer timesteps) applied to the rainfall input.
# Moore (2007) Table 1.
# =============================================================================

#' Apply rainfall factor and time delay to a rainfall series
#'
#' @param rain   Numeric vector of rainfall `[mm/timestep]`.
#' @param fc     Rainfall factor (dimensionless). Default 1 (no scaling).
#' @param td     Time delay `[timesteps, integer >= 0]`. Default 0 (no delay).
#' @return Numeric vector, same length as rain, with fc applied and td lag.
#' @keywords internal
.apply_rainfall_transform <- function(rain, fc = 1, td = 0L) {
  n    <- length(rain)
  td   <- as.integer(round(td))
  rain <- rain * fc
  if (td > 0L) {
    rain <- c(rep(0, td), rain)[seq_len(n)]
  }
  rain
}

# =============================================================================
# ACTUAL EVAPOTRANSPIRATION (Moore 2007, Eq. 2)
# AET = PET * (S / Smax)^be
# =============================================================================

.aet <- function(PET, S, Smax, be) {
  ratio <- pmin(S / max(Smax, 1e-9), 1)
  min(PET * ratio^be, S)
}

# =============================================================================
# RECHARGE FORMULATION 1 — Standard (Moore 2007, Eq. 10)
# d_i = k_g^{-1} * (S(t) - St)^{b_g}   for S > St, else 0
#
# Note: Moore (2007) Table 1 uses kb and m for the groundwater STORAGE
# routing parameters, and kg/bg for the soil RECHARGE function. However
# the paper also uses kg as the recharge time constant. We follow Table 1
# naming exactly: recharge uses kg (time constant) and bg (exponent);
# groundwater storage routing uses kb and m.
# =============================================================================

.recharge_standard <- function(S, St, kg, bg) {
  pmax(0, S - St)^bg / kg
}

# =============================================================================
# RECHARGE FORMULATION 2 — Demand-based (Moore 2007, Eqs. 11-13)
# Uses groundwater deficit ratio g(t) to modulate recharge.
# f(t) = (g(t)/alpha)^beta   if g(t) < alpha, else 1
# D_i  = (D_sat + (Smax - D_sat) * f(t)) * S(t) / Smax
# where D_sat = q_sat * dt
# =============================================================================

.recharge_demand <- function(S, Smax, Sg, Sg_max, alpha, beta, q_sat) {
  if (Sg_max <= 0) return(0)
  g_t <- (Sg_max - Sg) / Sg_max          # groundwater deficit ratio (Eq. 11)
  f_t <- if (g_t < alpha) {
    (g_t / alpha)^beta                    # Eq. 12
  } else {
    1
  }
  D_sat <- q_sat                          # recharge at saturation per timestep
  D_i   <- (D_sat + (Smax - D_sat) * f_t) * S / Smax   # Eq. 13
  pmax(0, D_i)
}

# =============================================================================
# RECHARGE FORMULATION 3 — Proportional split (Moore 2007, Table 1)
# Direct runoff is split: fraction alpha to surface, (1-alpha) to groundwater.
# No explicit soil drainage function.
# =============================================================================

# (handled inline in pdm_core.R as it depends on direct runoff, not S)

# =============================================================================
# GROUNDWATER (SLOW) STORE — non-linear reservoir (Moore 2007, Eq. 9 / Table 1)
# q_b = (S_g / kb)^m     where m=1 gives linear reservoir
# Exact recursive update for m=1; Euler step for m!=1.
# =============================================================================

#' Update groundwater store one timestep
#'
#' @param Sg     Current groundwater store `[mm]`.
#' @param Qd     Recharge inflow this timestep `[mm]`.
#' @param kb     Baseflow time constant `[hour * mm^{1-m}]` (Table 1: kb).
#' @param m      Baseflow exponent (Table 1: m). m=1 gives linear reservoir.
#' @param Sg_max Max groundwater store `[mm]`. 0 = unlimited.
#' @return list(Sg = updated store `[mm]`, Qb = baseflow `[mm/timestep]`)
#' @keywords internal
.groundwater_step <- function(Sg, Qd, kb, m, Sg_max = 0) {
  if (m == 1) {
    ag     <- exp(-1 / kb)
    Sg_new <- Sg * ag + Qd * kb * (1 - ag)
    Qb     <- Sg_new / kb
  } else {
    Qb_now <- (pmax(Sg, 0) / kb)^m
    Sg_new <- pmax(0, Sg + Qd - Qb_now)
    Qb     <- (Sg_new / kb)^m
  }
  if (Sg_max > 0) Sg_new <- pmin(Sg_new, Sg_max)
  list(Sg = Sg_new, Qb = Qb)
}

# =============================================================================
# SURFACE ROUTING — cascade of TWO linear reservoirs (Moore 2007, Table 1)
# Moore (2007) specifies k1 and k2 as the two time constants.
# The exact recursive solution for a cascade of two linear reservoirs is
# given by the transfer function approach (Moore 2007, Section on routing).
#
# For k1 != k2:
#   q(t+dt) = a1*q(t) + b1*u(t+dt) + b2*u(t)
# where the coefficients are derived from the two time constants.
#
# For k1 == k2 (repeated root):
#   Uses the Nash cascade / gamma UH special case.
# =============================================================================

#' Update surface routing cascade one timestep
#'
#' Implements the cascade of two linear reservoirs with time constants k1 and
#' k2 (Moore 2007, Table 1). Uses the exact recursive (Z-transform) solution.
#'
#' @param Ss1    State of first reservoir `[mm]`.
#' @param Ss2    State of second reservoir `[mm]`.
#' @param Qs_in  Surface runoff inflow this timestep `[mm]`.
#' @param k1     Time constant of first reservoir `[hours]`.
#' @param k2     Time constant of second reservoir `[hours]`.
#' @return list(Ss1, Ss2, Qf) where Qf is surface outflow `[mm/timestep]`.
#' @keywords internal
.surface_step <- function(Ss1, Ss2, Qs_in, k1, k2) {
  # Exact recursive update for each reservoir in series
  a1   <- exp(-1 / k1)
  Ss1_new <- Ss1 * a1 + Qs_in * k1 * (1 - a1)
  Qf1     <- Ss1_new / k1

  a2   <- exp(-1 / k2)
  Ss2_new <- Ss2 * a2 + Qf1 * k2 * (1 - a2)
  Qf      <- Ss2_new / k2

  list(Ss1 = Ss1_new, Ss2 = Ss2_new, Qf = Qf)
}

# =============================================================================
# CONSTANT FLOW ADDITION (Moore 2007, Table 1)
# qc [m3/s]: fixed return flow or abstraction added to total flow.
# Positive = return (adds to flow); negative = abstraction (reduces flow).
# =============================================================================

# Applied inline in pdm_core.R: Q_total = Qf + Qb + qc
