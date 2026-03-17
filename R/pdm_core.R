# =============================================================================
# Tool:        reach.hydro — PDM core model
# Description: Main pdm() entry point. Implements all components of Moore
#              (2007) Table 1:
#                - Rainfall factor (fc) and time delay (td)
#                - Probability-distributed soil moisture store (6 distributions)
#                - AET (be)
#                - Three recharge formulations: standard, demand, split
#                - Groundwater store: non-linear reservoir (kb, m)
#                - Surface routing: cascade of two linear reservoirs (k1, k2)
#                - Constant flow addition (qc)
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-03-16 - JP: full Table 1 implementation; corrected surface
#                               routing to two-reservoir cascade; added fc, td,
#                               qc; three recharge formulations; kb/m naming.
# Tier:        1
# References:
#   Moore, R.J. (2007). The PDM rainfall-runoff model.
#   Hydrol. Earth Syst. Sci., 11, 483-499.
# =============================================================================

#' Run the PDM rainfall-runoff model
#'
#' Implements the full Probability Distributed Model (PDM) of Moore (2007),
#' including all parameters from Table 1 of the paper.
#'
#' @param rain   Numeric vector of catchment-average rainfall \[mm/timestep\],
#'               OR a `Rainfall_Daily` / `Rainfall_15min` reach.io object.
#' @param pet    Numeric vector of potential evapotranspiration \[mm/timestep\],
#'               OR a reach.io HydroData object containing PET.
#' @param params A `PdmParams` object from [pdm_params()], or a named list.
#'               All Table 1 parameters are supported; see [pdm_params()] for
#'               the full reference.
#' @param dist   Capacity distribution (ignored if `params` is a `PdmParams`).
#'               One of `"pareto"` (default), `"rectangular"`, `"exponential"`,
#'               `"triangular"`, `"normal"`, `"lognormal"`.
#' @param S0     Initial soil moisture \[mm\]. Default: 50% of Smax.
#' @param Sg0    Initial groundwater store \[mm\]. Default: 0.
#' @param Ss10   Initial state of surface reservoir 1 \[mm\]. Default: 0.
#' @param Ss20   Initial state of surface reservoir 2 \[mm\]. Default: 0.
#'
#' @return A `ReachHydroResult` (data.table subclass) with columns:
#'   \describe{
#'     \item{rain}{Raw input rainfall \[mm/ts\]}
#'     \item{rain_eff}{Rainfall after fc scaling and td lag \[mm/ts\]}
#'     \item{pet}{Input PET \[mm/ts\]}
#'     \item{S}{Soil moisture \[mm\]}
#'     \item{cstar}{Critical capacity c* \[mm\]}
#'     \item{AET}{Actual evapotranspiration \[mm/ts\]}
#'     \item{SMD}{Soil moisture deficit = Smax - S \[mm\]}
#'     \item{direct_runoff}{Direct runoff from soil store \[mm/ts\]}
#'     \item{Qd}{Recharge to groundwater store \[mm/ts\]}
#'     \item{Qs}{Inflow to surface store \[mm/ts\]}
#'     \item{Qf}{Fast flow (surface cascade outflow) \[mm/ts\]}
#'     \item{Qb}{Baseflow (groundwater store outflow) \[mm/ts\]}
#'     \item{Q}{Total flow = Qf + Qb + qc \[mm/ts\]}
#'     \item{Sg}{Groundwater store state \[mm\]}
#'     \item{Ss1}{Surface reservoir 1 state \[mm\]}
#'     \item{Ss2}{Surface reservoir 2 state \[mm\]}
#'   }
#'
#' @references
#' Moore, R.J. (2007). The PDM rainfall-runoff model.
#' *Hydrology and Earth System Sciences*, 11, 483-499.
#' \doi{10.5194/hess-11-483-2007}
#'
#' @examples
#' set.seed(1)
#' rain <- pmax(0, rnorm(365, 2, 3))
#' pet  <- pmax(0, 2 + rnorm(365, 0, 0.3))
#' p    <- pdm_params(dist = "pareto", cmax = 350, b = 0.4,
#'                    fc = 1.0, td = 0L, kb = 150, m = 1, k1 = 5, k2 = 10)
#' res  <- pdm(rain, pet, params = p)
#' summary(res)
#'
#' @seealso [pdm_params()], [calibrate_pdm()], [compare_distributions()]
#' @export
pdm <- function(rain, pet,
                params = list(),
                dist   = "pareto",
                S0     = NULL,
                Sg0    = 0,
                Ss10   = 0,
                Ss20   = 0) {

  # ---- resolve reach.io HydroData inputs ------------------------------------
  if (.is_hydrodata(rain) || .is_hydrodata(pet)) {
    if (!.is_hydrodata(rain) || !.is_hydrodata(pet))
      stop("pdm(): if either `rain` or `pet` is a reach.io HydroData object, ",
           "both must be. Use as_pdm_input() to align and extract first.",
           call. = FALSE)
    inputs <- as_pdm_input(rain, pet)
    rain   <- inputs$rain
    pet    <- inputs$pet
  }

  # ---- resolve params -------------------------------------------------------
  if (inherits(params, "PdmParams")) {
    p    <- params
    dist <- p$dist
  } else {
    defaults <- list(
      dist = "pareto", fc = 1.0, td = 0L,
      cmin = 0, cmax = 400, b = 0.4,
      mu_c = 200, sigma_c = 80,
      mu_lnc = 5.0, sigma_lnc = 0.5,
      be = 5,
      recharge_type = "standard",
      kg = 200, bg = 1, St = 10,
      alpha = 0.4, beta = 1, q_sat = 2,
      kb = 200, m = 1, Sg_max = 0,
      k1 = 5, k2 = 5,
      qc = 0
    )
    p      <- modifyList(defaults, params)
    p$dist <- match.arg(dist, c("pareto", "rectangular", "exponential",
                                "triangular", "lognormal"))
    p$td   <- as.integer(round(p$td))
    class(p) <- "PdmParams"
    pdm_validate_params(p)
    dist <- p$dist
  }

  # ---- input checks ---------------------------------------------------------
  n <- length(rain)
  if (length(pet) != n) stop("`rain` and `pet` must be the same length.", call. = FALSE)
  if (n == 0L)          stop("`rain` is empty.", call. = FALSE)

  # ---- apply rainfall factor and time delay ---------------------------------
  rain_eff <- .apply_rainfall_transform(rain, fc = p$fc, td = p$td)

  # ---- pre-compute Smax and initial states ----------------------------------
  Smax <- capacity_smax(dist, p)
  if (is.null(S0)) S0 <- Smax * 0.5
  S0   <- pmax(0, pmin(S0,  Smax))
  Sg0  <- pmax(0, Sg0)
  Ss10 <- pmax(0, Ss10)
  Ss20 <- pmax(0, Ss20)

  # ---- pre-allocate output vectors (faster than set() in loop) --------------
  # TODO: This loop is the primary performance bottleneck for long high-
  # resolution records (e.g. 30 years at 15-minute timesteps ~ 1.05M
  # iterations). It is a strong candidate for Rcpp acceleration.
  # See GitHub issue #XX for tracking.
  #
  # Migration notes for when this is pursued:
  #   - The six distribution functions in pdm_distributions.R and the store
  #     update functions in pdm_stores.R are already isolated single-
  #     responsibility functions with no R-specific dependencies. They map
  #     directly to C++ with minimal refactoring.
  #   - The natural Rcpp interface is: pass all input vectors and scalar
  #     parameters in, return a named list of output vectors. The R wrapper
  #     (this function) handles validation, PdmParams construction, and
  #     coercion to ReachHydroResult — none of that changes.
  #   - The R-level API and all existing tests remain valid after the port.
  #   - As an interim step, the plain vector pre-allocation below (vs
  #     data.table::set() inside the loop) already gives a modest speedup.
  S_out  <- numeric(n); cstar_out <- numeric(n)
  AET_out <- numeric(n); SMD_out <- numeric(n)
  dr_out <- numeric(n); Qd_out <- numeric(n)
  Qs_out <- numeric(n); Qf_out <- numeric(n)
  Qb_out <- numeric(n); Q_out  <- numeric(n)
  Sg_out <- numeric(n); Ss1_out <- numeric(n); Ss2_out <- numeric(n)

  S   <- S0
  Sg  <- Sg0
  Ss1 <- Ss10
  Ss2 <- Ss20

  for (t in seq_len(n)) {
    P   <- max(0, rain_eff[t])
    PET <- max(0, pet[t])

    # AET
    Ea <- .aet(PET, S, Smax, p$be)

    # Direct runoff from probability-distributed soil store
    Qr <- capacity_runoff(P, S, dist, p)

    # Update soil moisture
    S <- pmax(0, pmin(S + P - Qr - Ea, Smax))

    # Recharge / routing split
    if (p$recharge_type == "split") {
      # Formulation 3: proportional split — no soil drainage
      Qd <- (1 - p$alpha) * Qr
      Qs <- p$alpha * Qr

    } else if (p$recharge_type == "demand") {
      # Formulation 2: demand-based recharge
      Qd <- .recharge_demand(S, Smax, Sg, p$Sg_max,
                             p$alpha, p$beta, p$q_sat)
      S  <- pmax(0, S - Qd)
      Qs <- Qr

    } else {
      # Formulation 1: standard recharge (default)
      Qd <- .recharge_standard(S, p$St, p$kg, p$bg)
      S  <- pmax(0, S - Qd)
      Qs <- Qr
    }

    # Groundwater store (kb, m)
    gw  <- .groundwater_step(Sg, Qd, p$kb, p$m, p$Sg_max)
    Sg  <- gw$Sg
    Qb  <- gw$Qb

    # Surface routing: cascade of two linear reservoirs (k1, k2)
    sf  <- .surface_step(Ss1, Ss2, Qs, p$k1, p$k2)
    Ss1 <- sf$Ss1
    Ss2 <- sf$Ss2
    Qf  <- sf$Qf

    # Total flow with constant flow addition (qc)
    Qt  <- Qf + Qb + p$qc

    S_out[t]   <- S
    cstar_out[t] <- capacity_cstar(S, dist, p)
    AET_out[t] <- Ea
    SMD_out[t] <- Smax - S
    dr_out[t]  <- Qr
    Qd_out[t]  <- Qd
    Qs_out[t]  <- Qs
    Qf_out[t]  <- Qf
    Qb_out[t]  <- Qb
    Q_out[t]   <- Qt
    Sg_out[t]  <- Sg
    Ss1_out[t] <- Ss1
    Ss2_out[t] <- Ss2
  }

  out <- data.table::data.table(
    rain          = rain,
    rain_eff      = rain_eff,
    pet           = pet,
    S             = S_out,
    cstar         = cstar_out,
    AET           = AET_out,
    SMD           = SMD_out,
    direct_runoff = dr_out,
    Qd            = Qd_out,
    Qs            = Qs_out,
    Qf            = Qf_out,
    Qb            = Qb_out,
    Q             = Q_out,
    Sg            = Sg_out,
    Ss1           = Ss1_out,
    Ss2           = Ss2_out
  )

  new_reach_hydro_result(out, params = p, Smax = Smax, dist = dist,
                         call = match.call())
}
