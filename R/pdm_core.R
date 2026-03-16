# =============================================================================
# Tool:        reach.hydro — PDM core model
# Description: Main pdm() entry point. Orchestrates the per-timestep model
#              loop, calling distribution, store, and ET functions.
#              Section 8 of the original PDM script, refactored to:
#                - Accept PdmParams objects (as well as plain lists)
#                - Return a ReachHydroResult (data.table subclass)
#                - Use data.table for allocation (fastverse convention)
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: ported + refactored for reach.hydro
# Tier:        1
# Inputs:      rain [mm/ts], pet [mm/ts], PdmParams or list, initial states
# Outputs:     ReachHydroResult (data.table subclass)
# Dependencies: data.table, pdm_classes.R, pdm_distributions.R, pdm_stores.R
# =============================================================================

#' Run the PDM rainfall-runoff model
#'
#' Implements the Probability Distributed Model (PDM) of Moore (2007).
#' The soil moisture store uses a selectable capacity distribution to compute
#' saturation-excess runoff. Drainage feeds a non-linear groundwater store
#' (baseflow) and a linear surface store (fast flow).
#'
#' @param rain   Numeric vector of catchment-average rainfall \[mm/timestep\],
#'               OR a `Rainfall_Daily` / `Rainfall_15min` object from reach.io.
#'               If a HydroData object is supplied, `pet` must also be a
#'               HydroData object; use [as_pdm_input()] for more control over
#'               quality filtering and alignment.
#' @param pet    Numeric vector of potential evapotranspiration \[mm/timestep\],
#'               OR any reach.io HydroData object containing PET values.
#'               Must be the same length as `rain` (after alignment).
#' @param params A `PdmParams` object (from [pdm_params()]) or a named list of
#'               parameter values. Missing values are filled with defaults.
#'               See [pdm_params()] for the full parameter reference.
#' @param dist   Capacity distribution. One of `"pareto"` (default),
#'               `"uniform"`, `"exponential"`, `"glogistic"`, `"normal"`,
#'               `"lognormal"`. Ignored if `params` is a `PdmParams` object
#'               (use `params$dist` in that case).
#' @param S0     Initial soil moisture \[mm\]. Default: 50% of Smax.
#' @param Sg0    Initial groundwater store \[mm\]. Default: 0.
#' @param Ss0    Initial surface store \[mm\]. Default: 0.
#'
#' @return A `ReachHydroResult` object (a `data.table` subclass) with columns:
#'   \describe{
#'     \item{rain}{Input rainfall \[mm/ts\]}
#'     \item{pet}{Input PET \[mm/ts\]}
#'     \item{S}{Soil moisture \[mm\]}
#'     \item{cstar}{Critical capacity c* \[mm\]}
#'     \item{AET}{Actual evapotranspiration \[mm/ts\]}
#'     \item{SMD}{Soil moisture deficit = Smax - S \[mm\]}
#'     \item{direct_runoff}{Saturation-excess runoff from soil store \[mm/ts\]}
#'     \item{Qd}{Drainage to groundwater store \[mm/ts\]}
#'     \item{Qs}{Inflow to surface store \[mm/ts\]}
#'     \item{Qf}{Fast (surface store) outflow \[mm/ts\]}
#'     \item{Qb}{Baseflow (groundwater store) outflow \[mm/ts\]}
#'     \item{Q}{Total flow = Qf + Qb \[mm/ts\]}
#'     \item{Sg}{Groundwater store state \[mm\]}
#'     \item{Ss}{Surface store state \[mm\]}
#'   }
#'   Attributes: `dist`, `Smax`, `params` (the resolved `PdmParams`).
#'
#' @references
#' Moore, R.J. (2007). The PDM rainfall-runoff model.
#' *Hydrology and Earth System Sciences*, 11, 483-499.
#' \doi{10.5194/hess-11-483-2007}
#'
#' @examples
#' # Minimal example with Pareto distribution (plain vectors)
#' set.seed(1)
#' rain <- pmax(0, rnorm(365, 2, 3))
#' pet  <- pmax(0, 2 + rnorm(365, 0, 0.3))
#' p    <- pdm_params(dist = "pareto", cmax = 350, b = 0.4)
#' res  <- pdm(rain, pet, params = p)
#' summary(res)
#'
#' # Using reach.io HydroData objects directly
#' \dontrun{
#' p      <- pdm_params(dist = "pareto", cmax = 350, b = 0.4)
#' res    <- pdm(rainfall_obj, pet_obj, params = p)
#'
#' # For more control over quality filtering and alignment:
#' inputs <- as_pdm_input(rainfall_obj, pet_obj)
#' res    <- pdm(inputs$rain, inputs$pet, params = p)
#' }
#'
#' @seealso [pdm_params()], [calibrate_pdm()], [compare_distributions()],
#'   [as_pdm_input()], [as_rain_input()]
#' @export
pdm <- function(rain, pet,
                params = list(),
                dist   = "pareto",
                S0     = NULL,
                Sg0    = 0,
                Ss0    = 0) {

  # ---- resolve params -------------------------------------------------------
  if (inherits(params, "PdmParams")) {
    p    <- params
    dist <- p$dist
  } else {
    # Merge user list over defaults then validate
    defaults <- list(
      cmin = 0, cmax = 400, b = 0.4,
      mu_c = 200, sigma_c = 80,
      mu_lnc = 5.0, sigma_lnc = 0.5,
      be = 5,
      St = 10, kg = 200, bg = 1,
      Sg_max = 0, ks = 10,
      use_split = FALSE, alpha = 0.4
    )
    p      <- modifyList(defaults, params)
    p$dist <- match.arg(
      dist,
      c("pareto", "uniform", "exponential", "glogistic", "normal", "lognormal")
    )
    class(p) <- "PdmParams"
    pdm_validate_params(p)
    dist <- p$dist
  }

  # ---- resolve reach.io HydroData inputs ------------------------------------
  if (.is_hydrodata(rain) || .is_hydrodata(pet)) {
    if (!.is_hydrodata(rain) || !.is_hydrodata(pet))
      stop(
        "pdm(): if either `rain` or `pet` is a reach.io HydroData object, ",
        "both must be. Supply plain numeric vectors or use as_pdm_input() ",
        "to align and extract the series first.",
        call. = FALSE)
    inputs <- as_pdm_input(rain, pet)
    rain   <- inputs$rain
    pet    <- inputs$pet
  }

  # ---- input checks ---------------------------------------------------------
  n <- length(rain)
  if (length(pet) != n) stop("`rain` and `pet` must be the same length.", call. = FALSE)
  if (n == 0L)          stop("`rain` is empty.", call. = FALSE)

  # ---- pre-compute Smax and initial state -----------------------------------
  Smax <- capacity_smax(dist, p)
  if (is.null(S0)) S0 <- Smax * 0.5

  S0  <- pmax(0, pmin(S0,  Smax))
  Sg0 <- pmax(0, Sg0)
  Ss0 <- pmax(0, Ss0)

  # ---- allocate output (data.table, pre-allocated for speed) ----------------
  out <- data.table::data.table(
    rain          = rain,
    pet           = pet,
    S             = NA_real_,
    cstar         = NA_real_,
    AET           = NA_real_,
    SMD           = NA_real_,
    direct_runoff = NA_real_,
    Qd            = NA_real_,
    Qs            = NA_real_,
    Qf            = NA_real_,
    Qb            = NA_real_,
    Q             = NA_real_,
    Sg            = NA_real_,
    Ss            = NA_real_
  )

  # ---- model loop -----------------------------------------------------------
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
  #   - As an interim step, replacing data.table::set() calls below with
  #     plain numeric vector assignment (filling pre-allocated vectors in the
  #     loop, then constructing the data.table once at the end) will give a
  #     modest speedup with no structural change.
  S  <- S0
  Sg <- Sg0
  Ss <- Ss0

  for (t in seq_len(n)) {
    P   <- max(0, rain[t])
    PET <- max(0, pet[t])

    # Actual ET (capped at available moisture)
    Ea <- .aet(PET, S, Smax, p$be)

    # Saturation-excess direct runoff
    Qr <- capacity_runoff(P, S, dist, p)

    # Update soil moisture after rainfall, runoff, and ET
    S <- pmax(0, pmin(S + P - Qr - Ea, Smax))

    # Route runoff: drainage function or proportional split
    if (isTRUE(p$use_split)) {
      Qd <- (1 - p$alpha) * Qr
      Qs <- p$alpha * Qr
    } else {
      Qd <- .soil_drainage(S, p$St)
      S  <- pmax(0, S - Qd)
      Qs <- Qr
    }

    # Groundwater store
    gw <- .groundwater_step(Sg, Qd, p$kg, p$bg, p$Sg_max)
    Sg <- gw$Sg
    Qb <- gw$Qb

    # Surface store
    sf <- .surface_step(Ss, Qs, p$ks)
    Ss <- sf$Ss
    Qf <- sf$Qf

    # Write to pre-allocated data.table by reference
    data.table::set(out, i = t, j = "S",             value = S)
    data.table::set(out, i = t, j = "cstar",         value = capacity_cstar(S, dist, p))
    data.table::set(out, i = t, j = "AET",           value = Ea)
    data.table::set(out, i = t, j = "SMD",           value = Smax - S)
    data.table::set(out, i = t, j = "direct_runoff", value = Qr)
    data.table::set(out, i = t, j = "Qd",            value = Qd)
    data.table::set(out, i = t, j = "Qs",            value = Qs)
    data.table::set(out, i = t, j = "Qf",            value = Qf)
    data.table::set(out, i = t, j = "Qb",            value = Qb)
    data.table::set(out, i = t, j = "Q",             value = Qf + Qb)
    data.table::set(out, i = t, j = "Sg",            value = Sg)
    data.table::set(out, i = t, j = "Ss",            value = Ss)
  }

  new_reach_hydro_result(out, params = p, Smax = Smax, dist = dist,
                         call = match.call())
}
