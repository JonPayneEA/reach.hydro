# =============================================================================
# Tool:        reach.hydro — PDM store update functions
# Description: Per-timestep update functions for each model store.
#              Sections 4-7 of the original PDM script. All unexported;
#              called inside the main pdm() loop in pdm_core.R.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: ported from PDM standalone script
# Tier:        1
# Inputs:      Store states, fluxes, parameters
# Outputs:     Updated store states and fluxes (scalars)
# Dependencies: base R only
# =============================================================================

# =============================================================================
# SECTION 4 — Actual evapotranspiration
# AET = PET * (S / Smax)^be    (Moore 2007, Eq. 2)
# =============================================================================

# @param PET   Potential ET this timestep [mm]
# @param S     Current soil moisture [mm]
# @param Smax  Basin storage capacity [mm]
# @param be    AET exponent (default 5)
# @return AET [mm], clamped to [0, min(PET, S)]
.aet <- function(PET, S, Smax, be) {
  ratio <- pmin(S / max(Smax, 1e-9), 1)
  min(PET * ratio^be, S)
}

# =============================================================================
# SECTION 5 — Soil drainage
# D = max(0, S - St) / kg_soil  (linear above tension threshold)
# =============================================================================

# @param S        Current soil moisture [mm]
# @param St       Tension storage threshold [mm]
# @param kg_soil  Drainage time constant [timesteps] — note: distinct from
#                 the groundwater kg parameter; here the routing uses the same
#                 kg for simplicity, passed explicitly.
# @return Drainage flux [mm/timestep]
.soil_drainage <- function(S, St, kg_soil = 1) {
  pmax(0, S - St) / kg_soil
}

# =============================================================================
# SECTION 6 — Groundwater (slow) store
# dSg/dt = Qd - Qb
# Qb = (Sg / kg)^bg           (Moore 2007, Eq. 9)
# bg == 1: linear — exact recursive update; bg != 1: Euler step.
# Optional inhibited recharge near Sg_max (Moore 2007, Eq. 11).
# =============================================================================

# @param Sg     Current groundwater store [mm]
# @param Qd     Drainage inflow this timestep [mm]
# @param kg     Groundwater time constant [timesteps]
# @param bg     Groundwater exponent (1 = linear reservoir)
# @param Sg_max Max groundwater store [mm]; 0 = unlimited
# @return list(Sg = updated store [mm], Qb = baseflow [mm/timestep])
.groundwater_step <- function(Sg, Qd, kg, bg, Sg_max = 0) {
  # Inhibited recharge
  if (Sg_max > 0) {
    f_inh <- pmax(0, 1 - Sg / Sg_max)
    Qd    <- Qd * f_inh
  }

  if (bg == 1) {
    # Exact recursive solution for linear reservoir
    ag     <- exp(-1 / kg)
    Sg_new <- Sg * ag + Qd * kg * (1 - ag)
    Qb     <- Sg_new / kg
  } else {
    # Euler step for non-linear reservoir
    Qb_now <- (pmax(Sg, 0) / kg)^bg
    Sg_new <- pmax(0, Sg + Qd - Qb_now)
    Qb     <- (Sg_new / kg)^bg
  }

  if (Sg_max > 0) Sg_new <- pmin(Sg_new, Sg_max)

  list(Sg = Sg_new, Qb = Qb)
}

# =============================================================================
# SECTION 7 — Surface (fast) store
# Linear reservoir; exact recursive solution.
# =============================================================================

# @param Ss     Current surface store [mm]
# @param Qs_in  Surface runoff inflow this timestep [mm]
# @param ks     Surface store time constant [timesteps]
# @return list(Ss = updated store [mm], Qf = fast flow [mm/timestep])
.surface_step <- function(Ss, Qs_in, ks) {
  as     <- exp(-1 / ks)
  Ss_new <- Ss * as + Qs_in * ks * (1 - as)
  Qf     <- Ss_new / ks
  list(Ss = Ss_new, Qf = Qf)
}
