# =============================================================================
# Tool:        reach.hydro — PDM S3 class definitions
# Description: PdmParams (validated parameter set) and ReachHydroResult
#              (typed model output). All parameters from Moore (2007) Table 1.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-03-16 - JP: added fc, td, qc; k1/k2 surface routing;
#                               three recharge formulations; kb/m naming.
# Tier:        1
# =============================================================================

#' Construct and validate a PDM parameter set
#'
#' Creates a `PdmParams` object containing all parameters from Moore (2007)
#' Table 1. Parameters are validated at construction time.
#'
#' @section Rainfall transform:
#' \describe{
#'   \item{fc}{Rainfall factor \[-\]. Scales catchment-average rainfall before
#'     it enters the soil store. Default 1 (no scaling).}
#'   \item{td}{Time delay \[timesteps\]. Integer lag applied to rainfall input.
#'     Default 0.}
#' }
#'
#' @section Probability-distributed soil store:
#' \describe{
#'   \item{dist}{Capacity distribution: `"pareto"` (default), `"rectangular"`,
#'     `"exponential"`, `"triangular"`, `"normal"`, `"lognormal"`.}
#'   \item{cmin}{Minimum store capacity \[mm\]. Default 0.}
#'   \item{cmax}{Maximum store capacity \[mm\]. Default 400.}
#'   \item{b}{Exponent of Pareto distribution (also used by triangular as
#'     shape context). Default 0.4.}
#'   \item{mu_c}{Mean capacity \[mm\] — normal distribution only. Default 200.}
#'   \item{sigma_c}{Std dev of capacity \[mm\] — normal only. Default 80.}
#'   \item{mu_lnc}{Mean of log(capacity) — lognormal only. Default 5.0.}
#'   \item{sigma_lnc}{Std dev of log(capacity) — lognormal only. Default 0.5.}
#' }
#'
#' @section Evaporation:
#' \describe{
#'   \item{be}{Exponent in actual evaporation function. Default 5.}
#' }
#'
#' @section Recharge formulation (select one via `recharge_type`):
#' \describe{
#'   \item{recharge_type}{`"standard"` (default), `"demand"`, or `"split"`.}
#'   \item{kg}{Recharge time constant \[timesteps\] — standard only. Default 200.}
#'   \item{bg}{Exponent of recharge function \[-\] — standard only. Default 1.}
#'   \item{St}{Soil tension storage capacity \[mm\] — standard only. Default 10.}
#'   \item{alpha}{Groundwater deficit ratio threshold \[-\] — demand; or runoff
#'     split fraction to surface store — split. Default 0.4.}
#'   \item{beta}{Exponent in groundwater demand factor function \[-\] — demand
#'     only. Default 1.}
#'   \item{q_sat}{Maximum recharge rate \[mm/timestep\] — demand only.
#'     Default 2.}
#' }
#'
#' @section Groundwater storage routing:
#' \describe{
#'   \item{kb}{Baseflow time constant \[hour * mm^{1-m}\]. Default 200.}
#'   \item{m}{Exponent of baseflow non-linear storage \[-\]. m=1 gives linear
#'     reservoir. Default 1.}
#'   \item{Sg_max}{Maximum groundwater store \[mm\]. 0 = unlimited. Default 0.}
#' }
#'
#' @section Surface routing (cascade of two linear reservoirs):
#' \describe{
#'   \item{k1}{Time constant of first surface reservoir \[timesteps\].
#'     Default 5.}
#'   \item{k2}{Time constant of second surface reservoir \[timesteps\].
#'     Default 5.}
#' }
#'
#' @section Constant flow addition:
#' \describe{
#'   \item{qc}{Constant flow \[m³/s\] added to total flow at every timestep.
#'     Positive = return flow; negative = abstraction. Default 0.}
#' }
#'
#' @return An object of class `"PdmParams"`.
#'
#' @examples
#' # Standard recharge, Pareto distribution
#' p <- pdm_params(dist = "pareto", cmax = 350, b = 0.4,
#'                 fc = 1.0, td = 0L, kb = 150, m = 1, k1 = 5, k2 = 10)
#'
#' # Demand-based recharge
#' p2 <- pdm_params(recharge_type = "demand",
#'                  alpha = 0.3, beta = 2, q_sat = 1.5, Sg_max = 300)
#'
#' @export
pdm_params <- function(
    dist         = "pareto",
    fc           = 1.0,
    td           = 0L,
    cmin         = 0,
    cmax         = 400,
    b            = 0.4,
    mu_c         = 200,
    sigma_c      = 80,
    mu_lnc       = 5.0,
    sigma_lnc    = 0.5,
    be           = 5,
    recharge_type = "standard",
    kg           = 200,
    bg           = 1,
    St           = 10,
    alpha        = 0.4,
    beta         = 1,
    q_sat        = 2,
    kb           = 200,
    m            = 1,
    Sg_max       = 0,
    k1           = 5,
    k2           = 5,
    qc           = 0
) {
  dist          <- match.arg(dist, c("pareto", "rectangular", "exponential",
                                     "triangular", "lognormal"))
  recharge_type <- match.arg(recharge_type, c("standard", "demand", "split"))

  p <- structure(
    list(
      dist = dist, fc = fc, td = as.integer(round(td)),
      cmin = cmin, cmax = cmax, b = b,
      mu_c = mu_c, sigma_c = sigma_c,
      mu_lnc = mu_lnc, sigma_lnc = sigma_lnc,
      be = be,
      recharge_type = recharge_type,
      kg = kg, bg = bg, St = St,
      alpha = alpha, beta = beta, q_sat = q_sat,
      kb = kb, m = m, Sg_max = Sg_max,
      k1 = k1, k2 = k2,
      qc = qc
    ),
    class = "PdmParams"
  )

  pdm_validate_params(p)
  p
}

#' Validate a PdmParams object
#'
#' @param p A `PdmParams` object or plain named list.
#' @return `p` invisibly if valid; stops with an informative error otherwise.
#' @export
pdm_validate_params <- function(p) {
  errs <- character(0)
  chk  <- function(cond, msg) if (!cond) errs <<- c(errs, msg)

  chk(p$cmax > 0,           "cmax must be > 0")
  chk(p$cmin >= 0,          "cmin must be >= 0")
  chk(p$cmin < p$cmax,      "cmin must be < cmax")
  chk(p$b > 0,              "b must be > 0")
  chk(p$mu_c > 0,           "mu_c must be > 0")
  chk(p$sigma_c > 0,        "sigma_c must be > 0")
  chk(p$sigma_lnc > 0,      "sigma_lnc must be > 0")
  chk(p$be > 0,             "be must be > 0")
  chk(p$fc > 0,             "fc (rainfall factor) must be > 0")
  chk(p$td >= 0L,           "td (time delay) must be >= 0")
  chk(p$kb > 0,             "kb (baseflow time constant) must be > 0")
  chk(p$m > 0,              "m (baseflow exponent) must be > 0")
  chk(p$Sg_max >= 0,        "Sg_max must be >= 0")
  chk(p$k1 > 0,             "k1 (surface reservoir 1 time constant) must be > 0")
  chk(p$k2 > 0,             "k2 (surface reservoir 2 time constant) must be > 0")

  if (p$recharge_type == "standard") {
    chk(p$kg > 0,           "kg must be > 0 for standard recharge")
    chk(p$bg > 0,           "bg must be > 0 for standard recharge")
    chk(p$St >= 0,          "St must be >= 0")
  }
  if (p$recharge_type == "demand") {
    chk(p$alpha > 0 && p$alpha <= 1, "alpha must be in (0,1] for demand recharge")
    chk(p$beta > 0,         "beta must be > 0 for demand recharge")
    chk(p$q_sat > 0,        "q_sat must be > 0 for demand recharge")
    chk(p$Sg_max > 0,       "Sg_max must be > 0 for demand recharge")
  }
  if (p$recharge_type == "split") {
    chk(p$alpha >= 0 && p$alpha <= 1, "alpha must be in [0,1] for split recharge")
  }

  if (length(errs) > 0)
    stop("Invalid PDM parameters:\n",
         paste0("  - ", errs, collapse = "\n"), call. = FALSE)

  invisible(p)
}

#' @export
print.PdmParams <- function(x, ...) {
  cat("<PdmParams>\n")
  cat(sprintf("  Distribution    : %s\n", x$dist))
  cat(sprintf("  fc / td         : %.3f / %d timesteps\n", x$fc, x$td))
  cat(sprintf("  cmin / cmax     : %.1f / %.1f mm\n", x$cmin, x$cmax))
  if (x$dist %in% c("pareto"))
    cat(sprintf("  b (shape)       : %.3f\n", x$b))
  if (x$dist == "normal")
    cat(sprintf("  mu_c / sigma_c  : %.1f / %.1f mm\n", x$mu_c, x$sigma_c))
  if (x$dist == "lognormal")
    cat(sprintf("  mu_lnc/sigma_lnc: %.2f / %.2f\n", x$mu_lnc, x$sigma_lnc))
  cat(sprintf("  be              : %.1f\n", x$be))
  cat(sprintf("  Recharge type   : %s\n", x$recharge_type))
  if (x$recharge_type == "standard")
    cat(sprintf("  kg / bg / St    : %.1f / %.2f / %.1f mm\n",
                x$kg, x$bg, x$St))
  if (x$recharge_type == "demand")
    cat(sprintf("  alpha/beta/q_sat: %.2f / %.2f / %.2f\n",
                x$alpha, x$beta, x$q_sat))
  if (x$recharge_type == "split")
    cat(sprintf("  split alpha     : %.2f\n", x$alpha))
  cat(sprintf("  kb / m          : %.1f / %.2f\n", x$kb, x$m))
  cat(sprintf("  k1 / k2         : %.1f / %.1f timesteps\n", x$k1, x$k2))
  if (x$qc != 0)
    cat(sprintf("  qc              : %.4f m3/s\n", x$qc))
  invisible(x)
}

# -----------------------------------------------------------------------------
# ReachHydroResult — typed wrapper around pdm() output
# -----------------------------------------------------------------------------

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
  cat(sprintf("  Timesteps (excl. warmup) : %d\n",  length(idx)))
  cat(sprintf("  Distribution             : %s\n",  attr(object, "dist")))
  cat(sprintf("  Smax                     : %.1f mm\n", attr(object, "Smax")))
  cat(sprintf("  Mean rainfall (raw)      : %.3f mm/ts\n", mean(object$rain[idx])))
  cat(sprintf("  Mean rainfall (fc/td)    : %.3f mm/ts\n", mean(object$rain_eff[idx])))
  cat(sprintf("  Mean AET                 : %.3f mm/ts\n", mean(object$AET[idx])))
  cat(sprintf("  Mean Q (total)           : %.3f mm/ts\n", mean(object$Q[idx])))
  cat(sprintf("  Mean Qf (fast)           : %.3f mm/ts\n", mean(object$Qf[idx])))
  cat(sprintf("  Mean Qb (baseflow)       : %.3f mm/ts\n", mean(object$Qb[idx])))
  cat(sprintf("  Baseflow index (BFI)     : %.3f\n",
              mean(object$Qb[idx]) / max(mean(object$Q[idx]), 1e-9)))
  invisible(object)
}
