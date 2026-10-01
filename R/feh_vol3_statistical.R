# =============================================================================
# Tool:        reach.hydro — FEH Vol 3 statistical flood frequency
# Description: Implements the FEH statistical methods for flood frequency
#              estimation from AMAX series:
#                - Single-site GLO/GEV fitting by L-moments
#                - Pooled (regional) analysis using the index flood method
#                - POT (Peaks over Threshold) Poisson-GPD model
#                - Uncertainty via bootstrap
#                - Growth curve and design flood output
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-23
# Modified:    2026-02-23 - JP: initial implementation
# Tier:        2
# Inputs:      AMAX series [m3/s]; donor site pool; return periods
# Outputs:     FehStatFit objects; growth curves; return period flows
# Dependencies: data.table, feh_lmoments.R
# References:
#   Robson, A. & Reed, D. (1999). Statistical Procedures for Flood Frequency
#     Estimation. FEH Vol. 3. CEH Wallingford.
#   Kjeldsen, T.R. et al. (2008). Improving the FEH statistical procedures
#     for flood frequency estimation. CEH Report.
# =============================================================================

# =============================================================================
# SECTION 1 : SINGLE-SITE FITTING
# =============================================================================

#' Fit a flood frequency distribution to a single AMAX series (FEH Vol 3)
#'
#' Fits the GLO (default, FEH-recommended for GB) or GEV to a series of
#' annual maximum flows using L-moments. Returns growth curve quantiles and
#' optional bootstrap confidence intervals.
#'
#' The index flood (median QMED) is separated from the growth curve so that
#' pooled analyses can substitute a catchment descriptor-based QMED estimate.
#'
#' @param amax        Numeric vector of annual maximum flows \[m³/s\],
#'                    OR a `FlodeFlow_Daily` / `FlodeFlow_15min` object from reach.io
#'                    (annual maxima are extracted automatically via [as_amax()]).
#'                    At least 10 years recommended; 15+ for reliable estimates.
#' @param dist        Distribution: `"glo"` (default), `"gev"`, or `"gno"`.
#' @param return_periods Numeric vector of return periods \[years\].
#' @param n_boot      Number of bootstrap replicates for CIs. Set to 0 to
#'                    skip. Default 1000.
#' @param ci_level    Confidence level for bootstrap CIs. Default 0.90.
#' @param record_name Optional character label for the gauging station.
#'
#' @return A list of class `"FehStatFit"` containing:
#'   \describe{
#'     \item{dist}{Distribution name.}
#'     \item{params}{Fitted distribution parameters.}
#'     \item{lmom}{Sample L-moments.}
#'     \item{qmed}{Median annual maximum (index flood).}
#'     \item{growth_curve}{`data.table` of return period, growth factor,
#'           flow, and (if bootstrapped) CI bounds.}
#'     \item{n_years}{Record length.}
#'     \item{record_name}{Label.}
#'     \item{amax}{Input data.}
#'     \item{boot}{Bootstrap growth factor matrix (n_boot × n_rp), or NULL.}
#'   }
#'
#' @examples
#' set.seed(1)
#' amax <- rgev_sim(50, xi = 50, alpha = 20, k = -0.1)
#' fit  <- feh_single_site(amax, dist = "glo")
#' fit$growth_curve
#'
#' @seealso [feh_pooled()], [feh_pot()]
#' @export
feh_single_site <- function(amax,
                            dist         = "glo",
                            return_periods = c(2, 5, 10, 20, 50, 100, 200, 1000),
                            n_boot       = 1000L,
                            ci_level     = 0.90,
                            record_name  = NULL) {

  dist <- match.arg(dist, c("glo", "gev", "gno"))

  # ---- resolve reach.io FlodeHydroData input ------------------------------------
  if (.is_hydrodata(amax)) {
    record_name <- record_name %||% amax@period_name
    amax        <- as_amax(amax)
  }

  amax <- amax[!is.na(amax)]
  n    <- length(amax)
  if (n < 5L) stop("feh_single_site: need at least 5 AMAX values.", call. = FALSE)
  if (n < 10L) warning("feh_single_site: record length < 10 years; ",
                       "estimates are highly uncertain.", call. = FALSE)

  lmom   <- sample_lmom(amax)
  params <- .fit_by_dist(lmom, dist)
  qmed   <- stats::median(amax)

  gc <- .growth_curve(params, dist, return_periods, qmed)

  # Bootstrap CIs
  boot_mat <- NULL
  if (n_boot > 0L) {
    boot_mat <- .bootstrap_gc(amax, dist, return_periods, n_boot)
    alpha_ci <- (1 - ci_level) / 2
    gc[, ci_lo := apply(boot_mat, 2, stats::quantile, probs = alpha_ci)]
    gc[, ci_hi := apply(boot_mat, 2, stats::quantile, probs = 1 - alpha_ci)]
  }

  structure(
    list(dist = dist, params = params, lmom = lmom, qmed = qmed,
         growth_curve = gc, n_years = n,
         record_name = record_name %||% "unnamed",
         amax = amax, boot = boot_mat),
    class = "FehStatFit"
  )
}

# =============================================================================
# SECTION 2 : POOLED (REGIONAL) ANALYSIS — INDEX FLOOD METHOD
# =============================================================================

#' Pooled flood frequency analysis using the index flood method (FEH Vol 3)
#'
#' Implements the FEH pooling procedure:
#' 1. Standardise each donor AMAX series by its QMED (index flood).
#' 2. Fit a regional GLO/GEV growth curve to the pooled standardised series.
#' 3. Scale by the subject site's QMED (from gauged data or catchment
#'    descriptors) to give design floods.
#'
#' The pooling group should be assembled externally (e.g., from NRFA data
#' based on catchment similarity) and passed as a named list.
#'
#' @param subject_amax   Numeric vector of AMAX at the subject site \[m³/s\],
#'                       OR a `FlodeFlow_Daily` / `FlodeFlow_15min` reach.io object.
#'                       May be `NULL` if `subject_qmed` is supplied directly.
#' @param donor_list     Named list of numeric AMAX vectors for donor sites.
#'                       Names are used as labels in output.
#' @param subject_qmed   QMED for the subject site \[m³/s\]. If `NULL`,
#'                       estimated as the median of `subject_amax`.
#' @param dist           Distribution: `"glo"` (default) or `"gev"`.
#' @param return_periods Numeric vector of return periods \[years\].
#' @param n_boot         Bootstrap replicates. Default 1000.
#' @param ci_level       Confidence level. Default 0.90.
#' @param min_pooled_years Minimum total pooled record length \[station-years\]
#'                       for a warning. FEH guidance: 5× the target return
#'                       period. Default 500.
#'
#' @return A list of class `"FehPooledFit"` containing:
#'   \describe{
#'     \item{dist, regional_params, regional_lmom}{Regional fit.}
#'     \item{subject_qmed}{QMED used for the subject site.}
#'     \item{pooled_years}{Total station-years in pool.}
#'     \item{n_donors}{Number of donor sites.}
#'     \item{growth_curve}{`data.table` with return period, growth factor,
#'           design flood, and CI bounds.}
#'     \item{donor_summary}{`data.table` of per-donor QMED and record length.}
#'   }
#'
#' @seealso [feh_single_site()]
#' @export
feh_pooled <- function(subject_amax   = NULL,
                       donor_list,
                       subject_qmed   = NULL,
                       dist           = "glo",
                       return_periods = c(2, 5, 10, 20, 50, 100, 200, 1000),
                       n_boot         = 1000L,
                       ci_level       = 0.90,
                       min_pooled_years = 500L) {

  dist <- match.arg(dist, c("glo", "gev"))

  # ---- resolve reach.io FlodeHydroData input ------------------------------------
  if (!is.null(subject_amax) && .is_hydrodata(subject_amax)) {
    subject_amax <- as_amax(subject_amax)
  }

  # Subject QMED
  if (is.null(subject_qmed)) {
    if (is.null(subject_amax))
      stop("feh_pooled: supply either subject_amax or subject_qmed.", call. = FALSE)
    subject_qmed <- stats::median(subject_amax[!is.na(subject_amax)])
  }

  # Include subject in pool if provided
  all_sites <- donor_list
  if (!is.null(subject_amax))
    all_sites[["_subject"]] <- subject_amax

  # Standardise each site by its own QMED
  donor_summary <- data.table::rbindlist(lapply(names(all_sites), function(nm) {
    x    <- all_sites[[nm]][!is.na(all_sites[[nm]])]
    data.table::data.table(
      site    = nm,
      n_years = length(x),
      qmed    = stats::median(x)
    )
  }))

  pooled_years <- sum(donor_summary$n_years)
  if (pooled_years < min_pooled_years)
    warning(sprintf(
      "feh_pooled: only %d pooled station-years; FEH recommends >= %d for ",
      pooled_years, min_pooled_years,
      "reliable estimates at long return periods."), call. = FALSE)

  # Pool standardised series
  std_series <- unlist(lapply(names(all_sites), function(nm) {
    x    <- all_sites[[nm]][!is.na(all_sites[[nm]])]
    qmed <- stats::median(x)
    x / qmed
  }))

  regional_lmom   <- sample_lmom(std_series)
  regional_params <- .fit_by_dist(regional_lmom, dist)

  # Growth curve for subject site
  gc <- .growth_curve(regional_params, dist, return_periods,
                      qmed_index = 1)  # growth factors; scale below
  gc[, design_flood := growth_factor * subject_qmed]

  # Bootstrap CIs on growth factors
  if (n_boot > 0L) {
    boot_mat <- .bootstrap_pooled_gc(all_sites, dist, return_periods, n_boot)
    alpha_ci <- (1 - ci_level) / 2
    gc[, ci_lo := apply(boot_mat, 2, stats::quantile, probs = alpha_ci) * subject_qmed]
    gc[, ci_hi := apply(boot_mat, 2, stats::quantile, probs = 1 - alpha_ci) * subject_qmed]
  }

  structure(
    list(dist            = dist,
         regional_params = regional_params,
         regional_lmom   = regional_lmom,
         subject_qmed    = subject_qmed,
         pooled_years    = pooled_years,
         n_donors        = length(donor_list),
         growth_curve    = gc,
         donor_summary   = donor_summary),
    class = "FehPooledFit"
  )
}

# =============================================================================
# SECTION 3 : POT (PEAKS OVER THRESHOLD) — Poisson-GPD MODEL
# =============================================================================

#' POT flood frequency analysis using the Poisson-GPD model (FEH Vol 3)
#'
#' Fits a Poisson process model to peak counts and a Generalised Pareto
#' Distribution (GPD) to threshold exceedances. The annual exceedance
#' probability is combined via P(Q > q) = 1 - exp(-lambda * F_gpd(q)).
#'
#' @param peaks        Numeric vector of independent peak flows \[m³/s\],
#'                     OR a `FlodeFlow_Daily` / `FlodeFlow_15min` reach.io object
#'                     (peaks are extracted automatically via [as_pot()]).
#'                     Or use [peaks_over_threshold()] / [as_pot()] directly.
#' @param threshold    Threshold used to extract peaks \[m³/s\].
#' @param n_years      Total record length \[years\] for Poisson rate estimation.
#' @param return_periods Numeric vector of return periods \[years\].
#' @param n_boot       Bootstrap replicates. Default 500.
#' @param ci_level     Confidence level. Default 0.90.
#'
#' @return A list of class `"FehPotFit"` containing:
#'   \describe{
#'     \item{lambda}{Mean annual rate of threshold exceedances.}
#'     \item{gpd_params}{List with `sigma` (scale) and `k` (shape).}
#'     \item{growth_curve}{`data.table` with return period, flow, CI bounds.}
#'     \item{n_peaks, threshold, n_years}{Metadata.}
#'   }
#'
#' @seealso [peaks_over_threshold()], [feh_single_site()]
#' @export
feh_pot <- function(peaks,
                    threshold,
                    n_years,
                    return_periods = c(2, 5, 10, 20, 50, 100, 200, 1000),
                    n_boot         = 500L,
                    ci_level       = 0.90) {

  # ---- resolve reach.io FlodeHydroData input ------------------------------------
  if (.is_hydrodata(peaks)) {
    pot     <- as_pot(peaks, threshold = threshold, min_sep = 3L)
    n_years <- attr(pot, "n_years")
    peaks   <- pot$peak_flow
  }

  peaks <- peaks[!is.na(peaks)]
  peaks <- peaks[peaks > threshold]
  n_p   <- length(peaks)
  if (n_p < 5L) stop("feh_pot: fewer than 5 exceedances above threshold.", call. = FALSE)

  lambda     <- n_p / n_years
  exceedances <- peaks - threshold

  # GPD fit by L-moments (Hosking & Wallis 1987)
  gpd_params <- .gpd_fit_lmom(exceedances)

  gc <- .pot_growth_curve(gpd_params, lambda, threshold, return_periods)

  # Bootstrap CIs
  if (n_boot > 0L) {
    alpha_ci  <- (1 - ci_level) / 2
    boot_flows <- matrix(NA_real_, nrow = n_boot, ncol = length(return_periods))
    for (b in seq_len(n_boot)) {
      samp_exc  <- sample(exceedances, n_p, replace = TRUE)
      samp_lam  <- stats::rpois(1, n_p) / n_years
      samp_gpd  <- tryCatch(.gpd_fit_lmom(samp_exc), error = function(e) NULL)
      if (is.null(samp_gpd) || samp_lam <= 0) next
      boot_flows[b, ] <- .pot_growth_curve(samp_gpd, samp_lam, threshold,
                                           return_periods)$flow
    }
    gc[, ci_lo := apply(boot_flows, 2, stats::quantile, probs = alpha_ci, na.rm = TRUE)]
    gc[, ci_hi := apply(boot_flows, 2, stats::quantile, probs = 1 - alpha_ci, na.rm = TRUE)]
  }

  structure(
    list(lambda = lambda, gpd_params = gpd_params, growth_curve = gc,
         n_peaks = n_p, threshold = threshold, n_years = n_years),
    class = "FehPotFit"
  )
}

# =============================================================================
# SECTION 4 : PRINT / SUMMARY METHODS
# =============================================================================

#' @export
print.FehStatFit <- function(x, ...) {
  cat(sprintf("<FehStatFit> dist=%s | n=%d yrs | QMED=%.2f m\u00b3/s\n",
              toupper(x$dist), x$n_years, x$qmed))
  cat(sprintf("  L-CV=%.3f | L-skew=%.3f | L-kurt=%.3f\n",
              x$lmom$t2, x$lmom$t3, x$lmom$t4))
  cat("\nGrowth curve:\n")
  print(x$growth_curve[, .(return_period_yr, growth_factor,
                            flow = round(flow, 2))], ...)
  invisible(x)
}

#' @export
print.FehPooledFit <- function(x, ...) {
  cat(sprintf(
    "<FehPooledFit> dist=%s | %d donors | %d station-years | QMED=%.2f\n",
    toupper(x$dist), x$n_donors, x$pooled_years, x$subject_qmed))
  cat("\nDesign floods:\n")
  cols <- intersect(c("return_period_yr", "growth_factor", "design_flood",
                      "ci_lo", "ci_hi"), names(x$growth_curve))
  print(x$growth_curve[, .SD, .SDcols = cols], ...)
  invisible(x)
}

#' @export
print.FehPotFit <- function(x, ...) {
  cat(sprintf("<FehPotFit> lambda=%.2f/yr | threshold=%.2f | n_peaks=%d\n",
              x$lambda, x$threshold, x$n_peaks))
  cat(sprintf("  GPD sigma=%.3f | k=%.4f\n",
              x$gpd_params$sigma, x$gpd_params$k))
  cat("\nReturn period flows:\n")
  print(x$growth_curve, ...)
  invisible(x)
}

# =============================================================================
# SECTION 5 : INTERNAL HELPERS
# =============================================================================

# Dispatch to distribution-specific L-moment fitter
.fit_by_dist <- function(lmom, dist) {
  switch(dist,
    glo = .glo_fit_lmom(lmom),
    gev = .gev_fit_lmom(lmom),
    gno = .gno_fit_lmom(lmom),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}

# Quantile function dispatch
.qfun <- function(p, params, dist) {
  switch(dist,
    glo = qglo(p, params),
    gev = qgev(p, params),
    gno = qgno(p, params),
    stop(paste("Unknown distribution:", dist), call. = FALSE)
  )
}

# Build growth curve data.table from fitted params
# qmed_index: set to QMED for absolute flows, 1 for dimensionless growth factors
.growth_curve <- function(params, dist, return_periods, qmed_index) {
  probs   <- 1 - 1 / return_periods
  quants  <- .qfun(probs, params, dist)
  q_med_param <- .qfun(0.5, params, dist)  # median of fitted distribution
  gf <- quants / q_med_param               # growth factors

  data.table::data.table(
    return_period_yr = return_periods,
    growth_factor    = round(gf, 4),
    flow             = round(gf * qmed_index, 3)
  )
}

# Bootstrap growth curve (single-site)
.bootstrap_gc <- function(amax, dist, return_periods, n_boot) {
  n        <- length(amax)
  probs    <- 1 - 1 / return_periods
  boot_gf  <- matrix(NA_real_, nrow = n_boot, ncol = length(return_periods))
  for (b in seq_len(n_boot)) {
    s     <- sample(amax, n, replace = TRUE)
    lm    <- tryCatch(sample_lmom(s), error = function(e) NULL)
    if (is.null(lm)) next
    p     <- tryCatch(.fit_by_dist(lm, dist), error = function(e) NULL)
    if (is.null(p)) next
    qmed_b <- .qfun(0.5, p, dist)
    if (qmed_b <= 0) next
    boot_gf[b, ] <- .qfun(probs, p, dist) / qmed_b
  }
  boot_gf
}

# Bootstrap growth curve (pooled)
.bootstrap_pooled_gc <- function(all_sites, dist, return_periods, n_boot) {
  probs    <- 1 - 1 / return_periods
  site_nms <- names(all_sites)
  boot_gf  <- matrix(NA_real_, nrow = n_boot, ncol = length(return_periods))
  for (b in seq_len(n_boot)) {
    # Resample each site with replacement, then pool
    std_b <- unlist(lapply(site_nms, function(nm) {
      x    <- all_sites[[nm]][!is.na(all_sites[[nm]])]
      qmed <- stats::median(x)
      s    <- sample(x, length(x), replace = TRUE)
      s / qmed
    }))
    lm <- tryCatch(sample_lmom(std_b), error = function(e) NULL)
    if (is.null(lm)) next
    p  <- tryCatch(.fit_by_dist(lm, dist), error = function(e) NULL)
    if (is.null(p)) next
    qmed_b <- .qfun(0.5, p, dist)
    if (qmed_b <= 0) next
    boot_gf[b, ] <- .qfun(probs, p, dist) / qmed_b
  }
  boot_gf
}

# GPD fit by L-moments (Hosking & Wallis 1987)
# Parameterisation: x(F) = sigma/k * [1 - (1-F)^k], k != 0
.gpd_fit_lmom <- function(exceedances) {
  lm    <- sample_lmom(exceedances)
  # tau3 = (1 - 3k) / (3 - k)  =>  k = (1 - 3*tau3) / (3 - tau3)
  # (only first two L-moments needed for GPD: l1 = sigma/(1-k), l2 = sigma/((1-k)(2-k)))
  k     <- (1 - 3 * lm$t3) / (3 - lm$t3)
  sigma <- lm$l2 * (1 - k) * (2 - k)
  if (sigma <= 0) stop("GPD fit: negative scale parameter.", call. = FALSE)
  list(sigma = sigma, k = k)
}

# POT return period flows: Q(T) = threshold + GPD quantile at
# F = 1 - 1/(lambda * T)  via Poisson-GPD combination
.pot_growth_curve <- function(gpd_params, lambda, threshold, return_periods) {
  # Annual exceedance probability = 1 - exp(-lambda * (1-F_gpd))
  # => F_gpd = 1 - (-log(1 - 1/T)) / lambda
  aep      <- 1 / return_periods
  f_gpd    <- 1 - (-log(1 - aep)) / lambda
  f_gpd    <- pmax(0, pmin(1 - 1e-9, f_gpd))
  k        <- gpd_params$k
  sigma    <- gpd_params$sigma
  if (abs(k) < 1e-8) {
    excess <- -sigma * log(1 - f_gpd)
  } else {
    excess <- sigma / k * (1 - (1 - f_gpd)^k)
  }
  data.table::data.table(
    return_period_yr = return_periods,
    flow             = round(threshold + excess, 3)
  )
}

# Null coalescing operator (defined here to avoid importing rlang)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Simulate from GEV for examples / tests
#' @export
rgev_sim <- function(n, xi, alpha, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  u <- stats::runif(n)
  qgev(u, list(xi = xi, alpha = alpha, k = k))
}
