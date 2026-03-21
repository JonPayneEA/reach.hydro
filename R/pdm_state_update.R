# =============================================================================
# Tool:        reach.hydro — PDM state updating
# Description: State updating for the PDM model. Runs the model over an
#              assimilation window of recent observations to produce
#              calibrated initial conditions for forecasting, removing the
#              need for a burn-in / spin-up period.
#
#              Two methods are implemented:
#
#              "window"    — Run the PDM over the assimilation window from a
#                            cold start. The final states become the forecast
#                            initial conditions. Effective when the window is
#                            long enough for the stores to equilibrate
#                            (>= 30 days recommended; >= 365 days optimal).
#
#              "inversion" — As "window", then algebraically adjusts the
#                            groundwater store (Sg) and surface routing stores
#                            (Ss1, Ss2) at the forecast origin so that the
#                            implied instantaneous flow matches obs_q. The
#                            most operationally useful method.
#
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-03-21
# Modified:    2026-03-21
# Tier:        1
# References:
#   Moore, R.J. (2007). The PDM rainfall-runoff model.
#   Hydrol. Earth Syst. Sci., 11, 483-499.
#   See design/pdm-state-updating.md for the full design rationale.
# =============================================================================

#' Update PDM model states using recent observations
#'
#' Runs the PDM over a recent assimilation window to produce calibrated
#' initial conditions for a forecast run, removing the need for a burn-in
#' (spin-up) period. Optionally applies an algebraic state adjustment at the
#' forecast origin so that simulated flow matches observed flow before the
#' forecast begins.
#'
#' @param rain   Numeric vector of rainfall \[mm/timestep\] covering the
#'               assimilation window. The last element is at the forecast
#'               origin (time T).
#' @param pet    Numeric vector of PET \[mm/timestep\], same length as `rain`.
#' @param obs_q  Numeric vector of observed flow \[mm/timestep\], same length
#'               as `rain`. For `method = "inversion"`, `obs_q[length(obs_q)]`
#'               (the value at the forecast origin) is used to adjust the
#'               stores. For `method = "window"`, used only in diagnostics.
#' @param params A `PdmParams` object or named list. Must be the same
#'               parameter set used for the subsequent forecast [pdm()] call.
#' @param method State update method:
#'   \describe{
#'     \item{`"inversion"`}{(Default) Runs the PDM over the window, then
#'       algebraically adjusts the groundwater store `Sg` and surface stores
#'       `Ss1`/`Ss2` so that their implied instantaneous outflow matches
#'       `obs_q` at the forecast origin. Soil moisture `S` is not adjusted —
#'       it is determined by the rainfall history and is robust to the
#'       assimilation window length. Recommended for operational use.}
#'     \item{`"window"`}{Runs the PDM over the assimilation window from cold
#'       initial conditions and returns the final states without further
#'       adjustment. Simple and robust; effective for windows >= 30 days.}
#'   }
#' @param S0   Initial soil moisture for the window run \[mm\]. Default:
#'             50% of Smax (same as [pdm()]).
#' @param Sg0  Initial groundwater store for the window run \[mm\]. Default 0.
#' @param Ss10 Initial surface reservoir 1 for the window run \[mm\]. Default 0.
#' @param Ss20 Initial surface reservoir 2 for the window run \[mm\]. Default 0.
#'
#' @return A `PdmStateUpdate` object (named list) with elements:
#'   \describe{
#'     \item{`S0`}{Updated soil moisture \[mm\] — pass as `S0` to [pdm()].}
#'     \item{`Sg0`}{Updated groundwater store \[mm\] — pass as `Sg0`.}
#'     \item{`Ss10`}{Updated surface reservoir 1 \[mm\] — pass as `Ss10`.}
#'     \item{`Ss20`}{Updated surface reservoir 2 \[mm\] — pass as `Ss20`.}
#'     \item{`method`}{The method used (`"window"` or `"inversion"`).}
#'     \item{`hindcast`}{[ReachHydroResult] for the assimilation window run,
#'       before any inversion adjustment. Useful for diagnosing window length
#'       adequacy and model fit during the assimilation period.}
#'     \item{`diagnostics`}{Named list with values at the forecast origin:
#'       `Q_obs`, `Q_sim_window` (flow at end of window run),
#'       `Q_implied_update` (instantaneous flow implied by adjusted stores),
#'       `Sg_before`, `Sg_after`, `Ss1_before`, `Ss2_before`, `Ss1_after`,
#'       `Ss2_after`, `scale_factor`.}
#'   }
#'
#' @section Usage pattern:
#' ```r
#' # Assimilation window: last 30 days of observations
#' upd <- update_pdm_states(rain_window, pet_window, obs_q_window,
#'                          params = cal$params,
#'                          method = "inversion")
#'
#' # Updated states as initial conditions for the forecast
#' fcast <- pdm(rain_fcast, pet_fcast, params = cal$params,
#'              S0   = upd$S0,
#'              Sg0  = upd$Sg0,
#'              Ss10 = upd$Ss10,
#'              Ss20 = upd$Ss20)
#' ```
#'
#' @section Method details:
#' **`"window"`**: Runs `pdm()` over the supplied window from cold starts.
#' The window length determines how well the stores equilibrate. For daily
#' timesteps, 30 days is the practical minimum; 365 days gives well-spun-up
#' groundwater stores. For sub-daily timesteps, scale proportionally.
#'
#' **`"inversion"`**: After the window run, uses the observed flow at the
#' forecast origin to scale the routing stores. The decomposition preserves
#' the modelled baseflow fraction (BFI) — only the magnitude is corrected,
#' not the partition between fast and slow flow. This is appropriate when
#' the model structure is trusted but the absolute level of the stores may
#' drift from reality during the assimilation window.
#'
#' For `recharge_type = "demand"` or `"split"`, inversion still runs but
#' only adjusts the routing stores (Sg, Ss1, Ss2); the decomposition
#' preserves the simulated Qb/Qf partition in the same way.
#'
#' @section Soil moisture note:
#' Soil moisture `S` is not adjusted by the inversion. It is primarily
#' controlled by the recent rainfall/PET history (captured by the window
#' run) rather than by observed flow, and algebraic inversion of S from
#' Q is not tractable without additional information (the relationship
#' between S and Q passes through the runoff and recharge functions).
#'
#' @seealso [pdm()], [pdm_params()]
#' @export
update_pdm_states <- function(rain, pet, obs_q, params,
                               method = c("inversion", "window"),
                               S0   = NULL,
                               Sg0  = 0,
                               Ss10 = 0,
                               Ss20 = 0) {
  method <- match.arg(method)
  n      <- length(rain)

  if (length(pet)   != n) stop("`rain` and `pet` must be the same length.",   call. = FALSE)
  if (length(obs_q) != n) stop("`rain` and `obs_q` must be the same length.", call. = FALSE)
  if (n < 2L)             stop("`rain` must have at least 2 timesteps.",       call. = FALSE)
  if (anyNA(rain))        stop("`rain` contains NA values.",  call. = FALSE)
  if (anyNA(pet))         stop("`pet` contains NA values.",   call. = FALSE)

  # ---- Run PDM over the assimilation window ---------------------------------
  hindcast <- pdm(rain, pet, params = params,
                  S0 = S0, Sg0 = Sg0, Ss10 = Ss10, Ss20 = Ss20)

  # States and flow at the forecast origin (end of window)
  S_win   <- hindcast$S[n]
  Sg_win  <- hindcast$Sg[n]
  Ss1_win <- hindcast$Ss1[n]
  Ss2_win <- hindcast$Ss2[n]
  Qb_win  <- hindcast$Qb[n]
  Qf_win  <- hindcast$Qf[n]
  Q_win   <- hindcast$Q[n]
  Q_obs   <- obs_q[n]

  p <- attr(hindcast, "params")

  # Diagnostics common to both methods
  diag <- list(
    Q_obs              = Q_obs,
    Q_sim_window       = Q_win,
    Q_implied_update   = Q_win,   # overwritten by inversion if applied
    Sg_before          = Sg_win,
    Sg_after           = Sg_win,
    Ss1_before         = Ss1_win,
    Ss2_before         = Ss2_win,
    Ss1_after          = Ss1_win,
    Ss2_after          = Ss2_win,
    scale_factor       = 1
  )

  if (method == "window") {
    return(.new_state_update(S_win, Sg_win, Ss1_win, Ss2_win,
                             method, hindcast, diag))
  }

  # ---- Inversion: adjust routing stores to match Q_obs at forecast origin ---

  if (is.na(Q_obs) || Q_obs < 0) {
    warning(
      "obs_q at forecast origin is NA or negative; ",
      "returning window states without inversion adjustment.",
      call. = FALSE
    )
    return(.new_state_update(S_win, Sg_win, Ss1_win, Ss2_win,
                             method, hindcast, diag))
  }

  # Subtract constant flow before computing scale (qc does not come from stores)
  Q_win_stores <- max(Q_win  - p$qc, 0)
  Q_obs_stores <- max(Q_obs  - p$qc, 0)

  if (Q_win_stores > 0) {
    scale <- Q_obs_stores / Q_win_stores
  } else if (Q_obs_stores > 0) {
    # Model predicted near-zero flow but obs is non-zero.
    # Set stores to a fraction of their equilibrium values rather than scaling
    # from zero, which would leave them at zero.
    scale <- 1   # stores unchanged; inversion cannot be applied reliably
    warning(
      "Simulated flow at forecast origin is near zero but obs_q is positive. ",
      "Cannot invert stores reliably; returning window states unchanged.",
      call. = FALSE
    )
    return(.new_state_update(S_win, Sg_win, Ss1_win, Ss2_win,
                             method, hindcast, diag))
  } else {
    scale <- 1   # both zero — nothing to do
  }

  # -- Groundwater store adjustment -------------------------------------------
  # Qb = (Sg / kb)^m  =>  Sg = kb * Qb^(1/m)
  # Scale Qb by the same factor as total flow (preserves BFI).
  Qb_adj <- Qb_win * scale
  Sg_adj  <- if (Qb_adj > 0) p$kb * Qb_adj^(1 / p$m) else 0
  if (p$Sg_max > 0) Sg_adj <- min(Sg_adj, p$Sg_max)
  Sg_adj  <- max(Sg_adj, 0)

  # -- Surface routing store adjustment ---------------------------------------
  # Scale Ss1 and Ss2 proportionally to preserve their ratio.
  # Instantaneous outflow of the cascade: Qf ≈ Ss2 / k2.
  # A uniform scale on both stores produces a uniform scale on Qf.
  if (Qf_win > 0) {
    surf_scale <- (Qf_win * scale) / Qf_win   # = scale, kept explicit for clarity
    Ss1_adj    <- max(Ss1_win * surf_scale, 0)
    Ss2_adj    <- max(Ss2_win * surf_scale, 0)
  } else {
    Ss1_adj <- Ss1_win
    Ss2_adj <- Ss2_win
  }

  # -- Diagnostic: implied instantaneous flow from adjusted stores ------------
  # This is the "drain rate" of the adjusted stores (P = 0, PET = 0).
  # For m = 1 (linear): Qb_implied = Sg_adj / kb
  # For m != 1:         Qb_implied = (Sg_adj / kb)^m
  Qb_implied <- if (Sg_adj > 0) (Sg_adj / p$kb)^p$m else 0
  Qf_implied <- if (Ss2_adj > 0) Ss2_adj / p$k2 else 0
  Q_implied  <- Qb_implied + Qf_implied + p$qc

  diag$Q_implied_update <- Q_implied
  diag$Sg_after         <- Sg_adj
  diag$Ss1_after        <- Ss1_adj
  diag$Ss2_after        <- Ss2_adj
  diag$scale_factor     <- scale

  .new_state_update(S_win, Sg_adj, Ss1_adj, Ss2_adj, method, hindcast, diag)
}


# -----------------------------------------------------------------------------
# PdmStateUpdate S3 class
# -----------------------------------------------------------------------------

.new_state_update <- function(S0, Sg0, Ss10, Ss20, method, hindcast, diagnostics) {
  structure(
    list(
      S0          = S0,
      Sg0         = Sg0,
      Ss10        = Ss10,
      Ss20        = Ss20,
      method      = method,
      hindcast    = hindcast,
      diagnostics = diagnostics
    ),
    class = "PdmStateUpdate"
  )
}

#' @export
print.PdmStateUpdate <- function(x, ...) {
  d <- x$diagnostics
  cat("<PdmStateUpdate>\n")
  cat(sprintf("  Method          : %s\n", x$method))
  cat(sprintf("  Window length   : %d timesteps\n", nrow(x$hindcast)))
  cat(sprintf("  --- States at forecast origin ---\n"))
  cat(sprintf("  S0  (soil)      : %.2f mm\n",  x$S0))
  cat(sprintf("  Sg0 (gw store)  : %.2f mm",    x$Sg0))
  if (x$method == "inversion" && !is.null(d$Sg_before))
    cat(sprintf("  [was %.2f mm]", d$Sg_before))
  cat("\n")
  cat(sprintf("  Ss10 / Ss20     : %.2f / %.2f mm\n", x$Ss10, x$Ss20))
  cat(sprintf("  --- Flow at forecast origin ---\n"))
  cat(sprintf("  Q obs           : %.4f mm/ts\n", d$Q_obs))
  cat(sprintf("  Q sim (window)  : %.4f mm/ts\n", d$Q_sim_window))
  if (x$method == "inversion")
    cat(sprintf("  Q implied (upd) : %.4f mm/ts\n", d$Q_implied_update))
  invisible(x)
}
