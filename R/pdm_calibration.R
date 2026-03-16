# =============================================================================
# Tool:        reach.hydro — PDM calibration and distribution comparison
# Description: Nelder-Mead calibration of PDM parameters against observed
#              flow (Section 10), and multi-distribution comparison utility
#              (Section 11) from the original PDM script. Refactored to return
#              data.table outputs and accept PdmParams objects.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: ported + refactored for reach.hydro
# Tier:        2 (calibration is analytical, not operational)
# Inputs:      rain, pet, obs_q time series; parameter bounds; distribution
# Outputs:     Calibrated PdmParams + simulation; distribution summary table
# Dependencies: data.table, pdm_core.R, pdm_metrics.R
# =============================================================================

#' Calibrate PDM parameters using Nelder-Mead optimisation
#'
#' Optimises a chosen set of PDM parameters to maximise NSE or KGE against
#' observed flow, excluding a warm-up period from the objective function.
#'
#' @param rain     Numeric vector of rainfall \[mm/timestep\].
#' @param pet      Numeric vector of PET \[mm/timestep\].
#' @param obs_q    Numeric vector of observed flow \[mm/timestep\].
#' @param dist     Capacity distribution (default `"pareto"`).
#' @param warmup   Integer. Timesteps excluded from objective. Default 365.
#' @param par_init Named list of initial parameter values. Default calibrates
#'                 `cmax`, `b`, `St`, `kg`, `ks`.
#' @param par_lo   Named list of lower bounds (same names as `par_init`).
#' @param par_hi   Named list of upper bounds (same names as `par_init`).
#' @param fixed    Named list of parameters held constant during calibration.
#' @param metric   Objective function: `"nse"` (default) or `"kge"`.
#' @param maxit    Max iterations for Nelder-Mead. Default 3000.
#'
#' @return A named list:
#'   \describe{
#'     \item{params}{Calibrated `PdmParams` object.}
#'     \item{dist}{Distribution name.}
#'     \item{nse}{NSE on the evaluation period.}
#'     \item{kge}{KGE on the evaluation period.}
#'     \item{pbias}{Percent bias on the evaluation period.}
#'     \item{optim}{Raw `optim()` output.}
#'     \item{sim}{`ReachHydroResult` for the calibrated run.}
#'   }
#'
#' @examples
#' \dontrun{
#' cal <- calibrate_pdm(rain, pet, obs_q,
#'                      dist = "pareto",
#'                      par_init = list(cmax = 300, b = 0.4, kg = 150))
#' cal$nse
#' }
#'
#' @seealso [pdm()], [compare_distributions()]
#' @export
calibrate_pdm <- function(rain, pet, obs_q,
                          dist     = "pareto",
                          warmup   = 365L,
                          par_init = list(cmax = 300, b = 0.4, St = 10,
                                          kg = 150, ks = 8),
                          par_lo   = list(cmax = 50,  b = 0.05, St = 0,
                                          kg = 10,  ks = 1),
                          par_hi   = list(cmax = 800, b = 2.0,  St = 100,
                                          kg = 500, ks = 60),
                          fixed    = list(cmin = 0, be = 5, bg = 1,
                                          Sg_max = 0, use_split = FALSE,
                                          alpha = 0.4),
                          metric   = "nse",
                          maxit    = 3000L) {

  metric <- match.arg(metric, c("nse", "kge"))
  dist   <- match.arg(dist,
    c("pareto", "uniform", "exponential", "glogistic", "normal", "lognormal"))

  pnames <- names(par_init)
  lower  <- unlist(par_lo[pnames])
  upper  <- unlist(par_hi[pnames])
  theta0 <- unlist(par_init[pnames])
  idx    <- seq(warmup + 1L, length(rain))

  if (length(idx) < 10L)
    warning("Fewer than 10 timesteps in the evaluation period after warmup.",
            call. = FALSE)

  obj <- function(theta) {
    if (any(theta < lower) | any(theta > upper)) return(1e6)
    p_trial <- fixed
    for (i in seq_along(pnames)) p_trial[[pnames[i]]] <- theta[i]
    p_trial$dist <- dist
    class(p_trial) <- "PdmParams"
    res <- tryCatch(
      pdm(rain, pet, params = p_trial),
      error = function(e) NULL
    )
    if (is.null(res)) return(1e6)
    val <- switch(metric,
      nse = nse(obs_q[idx], res$Q[idx]),
      kge = kge(obs_q[idx], res$Q[idx])
    )
    if (is.na(val)) 1e6 else -val
  }

  opt <- stats::optim(
    theta0, obj, method  = "Nelder-Mead",
    control = list(maxit = maxit, reltol = 1e-7)
  )

  # Reconstruct calibrated PdmParams
  cal_list <- fixed
  for (i in seq_along(pnames)) cal_list[[pnames[i]]] <- opt$par[i]
  cal_list$dist <- dist
  class(cal_list) <- "PdmParams"
  pdm_validate_params(cal_list)

  sim <- pdm(rain, pet, params = cal_list)

  list(
    params = cal_list,
    dist   = dist,
    nse    = nse(obs_q[idx],   sim$Q[idx]),
    kge    = kge(obs_q[idx],   sim$Q[idx]),
    pbias  = pbias(obs_q[idx], sim$Q[idx]),
    optim  = opt,
    sim    = sim
  )
}

#' Compare all six capacity distributions on the same forcing data
#'
#' Runs [pdm()] with each of the six supported distributions and returns a
#' summary `data.table` with flow statistics and, optionally, performance
#' metrics against observed flow.
#'
#' @param rain        Numeric vector of rainfall \[mm/timestep\].
#' @param pet         Numeric vector of PET \[mm/timestep\].
#' @param obs_q       Optional numeric vector of observed flow. When supplied,
#'                    NSE, KGE, and PBIAS columns are added.
#' @param base_params Named list of shared parameters (applied to all runs).
#' @param warmup      Integer. Timesteps excluded from performance metrics.
#'
#' @return A `data.table` with one row per distribution and columns:
#'   `distribution`, `Smax_mm`, `mean_Q`, `mean_Qf`, `mean_Qb`, `BFI`,
#'   `mean_AET`, and (if `obs_q` supplied) `NSE`, `KGE`, `PBIAS`.
#'
#' @examples
#' \dontrun{
#' tbl <- compare_distributions(rain, pet, obs_q = obs_q)
#' tbl[order(-NSE)]
#' }
#'
#' @seealso [pdm()], [calibrate_pdm()]
#' @export
compare_distributions <- function(rain, pet,
                                  obs_q       = NULL,
                                  base_params = list(
                                    cmin = 0, cmax = 350, b = 0.4,
                                    mu_c = 175, sigma_c = 70,
                                    mu_lnc = 5.0, sigma_lnc = 0.5,
                                    St = 20, kg = 150, ks = 8
                                  ),
                                  warmup = 0L) {

  dists <- c("pareto", "uniform", "exponential",
             "glogistic", "normal", "lognormal")
  idx   <- seq(warmup + 1L, length(rain))

  rows <- lapply(dists, function(d) {
    p      <- base_params
    p$dist <- d
    class(p) <- "PdmParams"
    res    <- pdm(rain, pet, params = p)
    Smax   <- attr(res, "Smax")
    mean_Q <- mean(res$Q)

    row <- data.table::data.table(
      distribution = d,
      Smax_mm      = round(Smax,          1),
      mean_Q       = round(mean(res$Q),   3),
      mean_Qf      = round(mean(res$Qf),  3),
      mean_Qb      = round(mean(res$Qb),  3),
      BFI          = round(mean(res$Qb) / max(mean_Q, 1e-9), 3),
      mean_AET     = round(mean(res$AET), 3)
    )

    if (!is.null(obs_q)) {
      row[, NSE   := round(nse(obs_q[idx],   res$Q[idx]), 3)]
      row[, KGE   := round(kge(obs_q[idx],   res$Q[idx]), 3)]
      row[, PBIAS := round(pbias(obs_q[idx], res$Q[idx]), 1)]
    }
    row
  })

  data.table::rbindlist(rows)
}
