# =============================================================================
# Tool:        reach.hydro S7 class definitions
# Description: Typed S7 classes for PDM parameter sets and model results.
#              PdmParams validates all parameters at construction time so
#              errors surface before the model loop runs.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: initial skeleton
# Tier:        1
# Inputs:      Named parameter values (see PdmParams below)
# Outputs:     PdmParams S7 object; ReachHydroResult S7 object
# Dependencies: S7 (base R >= 4.2 via the S7 package once stabilised;
#               currently using R7/S7 conventions)
# =============================================================================

# NOTE: S7 is the team's preferred OOP system for new Flode/REACH classes.
# When the S7 package is formally on CRAN and stable, replace these stubs
# with full S7::new_class() definitions. Until then, we use a constructor
# + validator pattern that is API-compatible with what S7 will expose.

# -----------------------------------------------------------------------------
# PdmParams — validated PDM parameter set
# -----------------------------------------------------------------------------

#' Construct and validate a PDM parameter set
#'
#' Creates a named list with class `"PdmParams"` containing all PDM parameters,
#' validated against distribution-specific requirements and physical bounds.
#'
#' @param dist      Capacity distribution. One of `"pareto"`, `"uniform"`,
#'                  `"exponential"`, `"glogistic"`, `"normal"`, `"lognormal"`.
#' @param cmin      Minimum storage capacity \[mm\]. Used by pareto, uniform,
#'                  glogistic. Default 0.
#' @param cmax      Maximum (or mean) storage capacity \[mm\]. Used by pareto,
#'                  uniform, exponential, glogistic. Default 400.
#' @param b         Shape parameter. Used by pareto, glogistic. Default 0.4.
#' @param mu_c      Mean capacity \[mm\]. Used by normal distribution. Default 200.
#' @param sigma_c   Std dev of capacity \[mm\]. Used by normal. Default 80.
#' @param mu_lnc    Mean of log(capacity). Used by lognormal. Default 5.0.
#' @param sigma_lnc Std dev of log(capacity). Used by lognormal. Default 0.5.
#' @param be        AET exponent. Default 5.
#' @param St        Tension threshold \[mm\]. Default 10.
#' @param kg        Groundwater time constant \[timesteps\]. Default 200.
#' @param bg        Groundwater exponent (1 = linear). Default 1.
#' @param Sg_max    Max groundwater store \[mm\]; 0 = unlimited. Default 0.
#' @param ks        Surface store time constant \[timesteps\]. Default 10.
#' @param use_split Use proportional split instead of drainage function.
#'                  Default `FALSE`.
#' @param alpha     Fraction of runoff to surface store (split only). Default 0.4.
#'
#' @return An object of class `"PdmParams"`.
#'
#' @examples
#' p <- pdm_params(dist = "pareto", cmax = 350, b = 0.4)
#' pdm_validate_params(p)
#'
#' @export
pdm_params <- function(
    dist      = "pareto",
    cmin      = 0,
    cmax      = 400,
    b         = 0.4,
    mu_c      = 200,
    sigma_c   = 80,
    mu_lnc    = 5.0,
    sigma_lnc = 0.5,
    be        = 5,
    St        = 10,
    kg        = 200,
    bg        = 1,
    Sg_max    = 0,
    ks        = 10,
    use_split = FALSE,
    alpha     = 0.4
) {
  dist <- match.arg(
    dist,
    c("pareto", "uniform", "exponential", "glogistic", "normal", "lognormal")
  )

  p <- structure(
    list(
      dist      = dist,
      cmin      = cmin,
      cmax      = cmax,
      b         = b,
      mu_c      = mu_c,
      sigma_c   = sigma_c,
      mu_lnc    = mu_lnc,
      sigma_lnc = sigma_lnc,
      be        = be,
      St        = St,
      kg        = kg,
      bg        = bg,
      Sg_max    = Sg_max,
      ks        = ks,
      use_split = use_split,
      alpha     = alpha
    ),
    class = "PdmParams"
  )

  pdm_validate_params(p)
  p
}

#' Validate a PdmParams object
#'
#' Checks physical plausibility of all parameters. Called automatically by
#' [pdm_params()]. Can be called explicitly after manual modification.
#'
#' @param p A `PdmParams` object (or plain named list).
#' @return `p` invisibly if valid; stops with an informative error otherwise.
#' @export
pdm_validate_params <- function(p) {
  errs <- character(0)

  chk <- function(cond, msg) if (!cond) errs <<- c(errs, msg)

  chk(p$cmax > 0,          "cmax must be > 0")
  chk(p$cmin >= 0,         "cmin must be >= 0")
  chk(p$cmin < p$cmax,     "cmin must be < cmax")
  chk(p$b > 0,             "b (shape) must be > 0")
  chk(p$mu_c > 0,          "mu_c must be > 0")
  chk(p$sigma_c > 0,       "sigma_c must be > 0")
  chk(p$sigma_lnc > 0,     "sigma_lnc must be > 0")
  chk(p$be > 0,            "be (AET exponent) must be > 0")
  chk(p$St >= 0,           "St (tension threshold) must be >= 0")
  chk(p$kg > 0,            "kg (groundwater time constant) must be > 0")
  chk(p$bg > 0,            "bg (groundwater exponent) must be > 0")
  chk(p$Sg_max >= 0,       "Sg_max must be >= 0")
  chk(p$ks > 0,            "ks (surface time constant) must be > 0")
  chk(p$alpha >= 0 && p$alpha <= 1, "alpha must be in [0, 1]")

  if (length(errs) > 0) {
    stop(
      "Invalid PDM parameters:\n",
      paste0("  - ", errs, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(p)
}

#' @export
print.PdmParams <- function(x, ...) {
  cat("<PdmParams>\n")
  cat(sprintf("  Distribution : %s\n", x$dist))
  cat(sprintf("  cmin / cmax  : %.1f / %.1f mm\n", x$cmin, x$cmax))
  if (x$dist %in% c("pareto", "glogistic"))
    cat(sprintf("  b (shape)    : %.3f\n", x$b))
  if (x$dist == "normal")
    cat(sprintf("  mu_c / sigma_c : %.1f / %.1f mm\n", x$mu_c, x$sigma_c))
  if (x$dist == "lognormal")
    cat(sprintf("  mu_lnc / sigma_lnc : %.2f / %.2f\n", x$mu_lnc, x$sigma_lnc))
  cat(sprintf("  be / St      : %.1f / %.1f\n", x$be, x$St))
  cat(sprintf("  kg / bg      : %.1f / %.2f\n", x$kg, x$bg))
  cat(sprintf("  ks           : %.1f\n", x$ks))
  if (x$use_split)
    cat(sprintf("  split alpha  : %.2f\n", x$alpha))
  invisible(x)
}

# -----------------------------------------------------------------------------
# ReachHydroResult — typed wrapper around pdm() output
# -----------------------------------------------------------------------------

#' Construct a ReachHydroResult object
#'
#' Internal constructor called by [pdm()]. Wraps the output `data.table` with
#' metadata needed for downstream methods (print, summary, plot dispatch).
#'
#' @param dt     `data.table` of per-timestep model output.
#' @param params A `PdmParams` object.
#' @param Smax   Basin storage capacity \[mm\].
#' @param dist   Capacity distribution name.
#' @param call   The matched call (for provenance).
#'
#' @return Object of class `c("ReachHydroResult", "data.table", "data.frame")`.
#' @keywords internal
new_reach_hydro_result <- function(dt, params, Smax, dist, call = NULL) {
  structure(
    dt,
    class  = c("ReachHydroResult", "data.table", "data.frame"),
    params = params,
    Smax   = Smax,
    dist   = dist,
    call   = call
  )
}

#' @export
print.ReachHydroResult <- function(x, ...) {
  cat(sprintf(
    "<ReachHydroResult> n=%d timesteps | dist=%s | Smax=%.1f mm\n",
    nrow(x), attr(x, "dist"), attr(x, "Smax")
  ))
  cat(sprintf(
    "  mean Q=%.3f | mean Qf=%.3f | mean Qb=%.3f | BFI=%.3f\n",
    mean(x$Q), mean(x$Qf), mean(x$Qb),
    mean(x$Qb) / max(mean(x$Q), 1e-9)
  ))
  invisible(x)
}

#' @export
summary.ReachHydroResult <- function(object, warmup = 0, ...) {
  idx <- seq(warmup + 1L, nrow(object))
  cat("<ReachHydroResult summary>\n")
  cat(sprintf("  Timesteps (excl. warmup) : %d\n", length(idx)))
  cat(sprintf("  Distribution             : %s\n", attr(object, "dist")))
  cat(sprintf("  Smax                     : %.1f mm\n", attr(object, "Smax")))
  cat(sprintf("  Mean rainfall            : %.3f mm/ts\n", mean(object$rain[idx])))
  cat(sprintf("  Mean AET                 : %.3f mm/ts\n", mean(object$AET[idx])))
  cat(sprintf("  Mean Q (total)           : %.3f mm/ts\n", mean(object$Q[idx])))
  cat(sprintf("  Mean Qf (fast)           : %.3f mm/ts\n", mean(object$Qf[idx])))
  cat(sprintf("  Mean Qb (baseflow)       : %.3f mm/ts\n", mean(object$Qb[idx])))
  cat(sprintf("  Baseflow index (BFI)     : %.3f\n",
              mean(object$Qb[idx]) / max(mean(object$Q[idx]), 1e-9)))
  invisible(object)
}
