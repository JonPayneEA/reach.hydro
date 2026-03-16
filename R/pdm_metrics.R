# =============================================================================
# Tool:        reach.hydro — goodness-of-fit metrics
# Description: Standard hydrological performance metrics.
#              Section 9 of the original PDM script plus FAR.
#              All functions handle NA pairs and zero-denominator edge cases.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: ported + added FAR; input validation via checkmate
# Tier:        1
# Inputs:      Numeric vectors of observed and simulated values
# Outputs:     Scalar metric values
# Dependencies: checkmate
# =============================================================================

# Internal helper: strip paired NAs and check length
.metric_prep <- function(obs, sim, fn_name) {
  checkmate::assert_numeric(obs, any.missing = TRUE, min.len = 2L,
                            .var.name = "obs")
  checkmate::assert_numeric(sim, any.missing = TRUE, min.len = 2L,
                            .var.name = "sim")
  if (length(obs) != length(sim))
    stop(sprintf("%s: `obs` and `sim` must be the same length.", fn_name),
         call. = FALSE)
  v <- !is.na(obs) & !is.na(sim)
  if (sum(v) < 2L)
    stop(sprintf("%s: fewer than 2 non-NA pairs.", fn_name), call. = FALSE)
  list(obs = obs[v], sim = sim[v])
}

# =============================================================================
# Nash-Sutcliffe Efficiency
# NSE = 1 - sum((obs - sim)^2) / sum((obs - mean(obs))^2)
# Range: (-Inf, 1]; 1 = perfect; 0 = mean-flow benchmark.
# =============================================================================

#' Nash-Sutcliffe Efficiency
#'
#' @param obs Numeric vector of observed flow \[any consistent unit\].
#' @param sim Numeric vector of simulated flow, same length as `obs`.
#' @return NSE scalar. `NA` if denominator is zero (constant observed series).
#' @references Nash & Sutcliffe (1970). J. Hydrol. 10(3), 282-290.
#' @export
nse <- function(obs, sim) {
  d  <- .metric_prep(obs, sim, "nse")
  ss <- sum((d$obs - mean(d$obs))^2)
  if (ss == 0) return(NA_real_)
  1 - sum((d$obs - d$sim)^2) / ss
}

# =============================================================================
# Kling-Gupta Efficiency
# KGE = 1 - sqrt((r-1)^2 + (beta-1)^2 + (gamma-1)^2)
# r = Pearson correlation; beta = mean ratio; gamma = CV ratio.
# Range: (-Inf, 1]; 1 = perfect.
# =============================================================================

#' Kling-Gupta Efficiency
#'
#' @inheritParams nse
#' @return KGE scalar. `NA` if observed mean or CV is zero.
#' @references Gupta et al. (2009). J. Hydrol. 377(1-2), 80-91.
#' @export
kge <- function(obs, sim) {
  d     <- .metric_prep(obs, sim, "kge")
  mu_o  <- mean(d$obs)
  mu_s  <- mean(d$sim)
  if (mu_o == 0) return(NA_real_)
  cv_o  <- sd(d$obs) / mu_o
  if (cv_o == 0) return(NA_real_)
  r     <- stats::cor(d$obs, d$sim)
  beta  <- mu_s / mu_o
  gamma <- (sd(d$sim) / max(mu_s, 1e-12)) / cv_o
  1 - sqrt((r - 1)^2 + (beta - 1)^2 + (gamma - 1)^2)
}

# =============================================================================
# Percent bias
# PBIAS = 100 * sum(sim - obs) / sum(obs)
# =============================================================================

#' Percent bias
#'
#' @inheritParams nse
#' @return PBIAS \[%\]. Positive = over-prediction; negative = under-prediction.
#'   `NA` if sum of observed values is zero.
#' @export
pbias <- function(obs, sim) {
  d <- .metric_prep(obs, sim, "pbias")
  s <- sum(d$obs)
  if (s == 0) return(NA_real_)
  100 * sum(d$sim - d$obs) / s
}

# =============================================================================
# False Alarm Ratio
# FAR = FP / (TP + FP)
# Requires binary threshold exceedance vectors.
# =============================================================================

#' False Alarm Ratio for threshold exceedances
#'
#' @param obs_exceed Logical (or 0/1) vector: did observed flow exceed threshold?
#' @param sim_exceed Logical (or 0/1) vector: did simulated flow exceed threshold?
#' @return FAR scalar in \[0, 1\]. `NA` if no simulated exceedances exist.
#' @references WMO (2008). Forecast Verification: Issues, Methods and FAQ.
#' @export
far <- function(obs_exceed, sim_exceed) {
  checkmate::assert(
    checkmate::check_logical(obs_exceed),
    checkmate::check_integerish(obs_exceed, lower = 0, upper = 1)
  )
  checkmate::assert(
    checkmate::check_logical(sim_exceed),
    checkmate::check_integerish(sim_exceed, lower = 0, upper = 1)
  )
  if (length(obs_exceed) != length(sim_exceed))
    stop("far: `obs_exceed` and `sim_exceed` must be the same length.", call. = FALSE)

  TP <- sum( obs_exceed &  sim_exceed, na.rm = TRUE)
  FP <- sum(!obs_exceed &  sim_exceed, na.rm = TRUE)

  if ((TP + FP) == 0L) return(NA_real_)
  FP / (TP + FP)
}

# =============================================================================
# Convenience: compute all metrics at once
# =============================================================================

#' Compute all goodness-of-fit metrics in one call
#'
#' @inheritParams nse
#' @param warmup Integer. Number of leading timesteps to exclude. Default 0.
#' @return Named numeric vector: `nse`, `kge`, `pbias`.
#' @export
gof_metrics <- function(obs, sim, warmup = 0L) {
  idx <- seq(warmup + 1L, length(obs))
  c(
    nse   = nse(obs[idx],   sim[idx]),
    kge   = kge(obs[idx],   sim[idx]),
    pbias = pbias(obs[idx], sim[idx])
  )
}
