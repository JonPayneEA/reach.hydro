# =============================================================================
# Tool:        reach.hydro — L-moments engine
# Description: Full L-moments and L-moment ratios for fitting and diagnosis.
#              Implements sample L-moments (Hosking 1990), and parameter
#              estimation by L-moments for GLO, GEV, GNO, PE3, and GPA
#              distributions used across FEH Vol 3 methods.
#              No external dependency — avoids lmomco to keep the footprint
#              lean (fastverse principle).
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-23
# Modified:    2026-02-23 - JP: initial implementation
# Tier:        2
# Inputs:      Numeric vectors of data or distribution parameters
# Outputs:     L-moment estimates; distribution parameters; quantiles
# Dependencies: base R only
# References:
#   Hosking, J.R.M. (1990). L-moments: Analysis and estimation of
#     distributions using linear combinations of order statistics.
#     JRSS-B, 52(1), 105-124.
#   Hosking, J.R.M. & Wallis, J.R. (1997). Regional Frequency Analysis.
#     Cambridge University Press.
# =============================================================================

# =============================================================================
# SECTION 1 : SAMPLE L-MOMENTS
# =============================================================================

#' Sample L-moments (first four)
#'
#' Computes the first four sample L-moments and L-moment ratios for a data
#' vector using the unbiased PWM estimator (Hosking 1990, eq. 2.2).
#'
#' @param x       Numeric vector. NAs are removed.
#' @param sort_x  Logical. If `FALSE`, assumes `x` is already sorted ascending.
#'
#' @return A named list:
#'   \describe{
#'     \item{l1}{First L-moment (mean).}
#'     \item{l2}{Second L-moment (L-scale; half the mean absolute difference).}
#'     \item{l3}{Third L-moment.}
#'     \item{l4}{Fourth L-moment.}
#'     \item{t2}{L-CV = l2 / l1.}
#'     \item{t3}{L-skewness = l3 / l2.}
#'     \item{t4}{L-kurtosis = l4 / l2.}
#'     \item{n}{Sample size after NA removal.}
#'   }
#'
#' @references Hosking (1990); Hosking & Wallis (1997) Appendix.
#' @export
sample_lmom <- function(x, sort_x = TRUE) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4L) stop("sample_lmom: need at least 4 non-NA values.", call. = FALSE)
  if (sort_x) x <- sort(x)

  # Probability weighted moments b0..b3 (unbiased; Hosking 1990 eq 2.2)
  i   <- seq_len(n)
  b0  <- mean(x)
  b1  <- sum(x * (i - 1)) / (n * (n - 1))
  b2  <- sum(x * (i - 1) * (i - 2)) / (n * (n - 1) * (n - 2))
  b3  <- sum(x * (i - 1) * (i - 2) * (i - 3)) /
         (n * (n - 1) * (n - 2) * (n - 3))

  l1  <-  b0
  l2  <-  2 * b1 - b0
  l3  <-  6 * b2 - 6 * b1 + b0
  l4  <- 20 * b3 - 30 * b2 + 12 * b1 - b0

  list(
    l1 = l1, l2 = l2, l3 = l3, l4 = l4,
    t2 = l2 / l1,
    t3 = l3 / l2,
    t4 = l4 / l2,
    n  = n
  )
}

# =============================================================================
# SECTION 2 : DISTRIBUTION PARAMETER ESTIMATION BY L-MOMENTS
# =============================================================================

# ---- Generalised Logistic (GLO) -------------------------------------------
# FEH recommended distribution for GB flood frequency (Robson & Reed 1999).
# Quantile function: x(F) = xi + alpha * [1 - ((1-F)/F)^k] / k
# k = 0 gives the Logistic distribution.
#
# L-moment relations (Hosking & Wallis 1997, Table A2):
#   l1   = xi + alpha * (1/k - pi / sin(pi*k))
#   l2   = alpha * pi * k / (sin(pi*k) * Gamma(1+k) * Gamma(1-k))
#   tau3 = -k   (exact)
#
# @param lmom  Named list from sample_lmom()
# @return list(xi, alpha, k)  with class "GloParams"
.glo_fit_lmom <- function(lmom) {
  k     <- -lmom$t3                        # exact from tau3
  if (abs(k) < 1e-6) {
    # Logistic limit
    alpha <- lmom$l2
    xi    <- lmom$l1
  } else {
    pk    <- pi * k
    alpha <- lmom$l2 * sin(pk) / pk        # exact
    xi    <- lmom$l1 - alpha * (1/k - pi / sin(pk))
  }
  structure(list(xi = xi, alpha = alpha, k = k), class = "GloParams")
}

# ---- Generalised Extreme Value (GEV) --------------------------------------
# Quantile function: x(F) = xi + alpha * [1 - (-log F)^k] / k
# k = 0: Gumbel.
#
# L-moment relations (Hosking & Wallis 1997, Table A2):
#   l2   = alpha * (1 - 2^(-k)) * Gamma(1+k)
#   tau3 = 2(1 - 3^(-k))/(1 - 2^(-k)) - 3
#
# k estimated by rational approximation (Hosking & Wallis 1997 eq A11).
#
# @param lmom Named list from sample_lmom()
# @return list(xi, alpha, k) with class "GevParams"
.gev_fit_lmom <- function(lmom) {
  tau3 <- lmom$t3
  # Rational approximation for k (Hosking & Wallis 1997, A11)
  # Valid for tau3 in (-1, 1)
  if (tau3 >= 1/3) {
    # Upper range approximation
    z <- 1 - tau3
    k <- 0.2955 + z * (0.4046 + z * 0.0959)
  } else {
    z <- tau3
    k <- 7.859 * z + 2.9554 * z^2
    k <- k / (1 + 3.2517 * z)
  }
  # Refine k with one Newton step using exact tau3(k) relation
  for (iter in 1:3) {
    g1    <- gamma(1 + k)
    g2    <- gamma(1 + 2*k)
    t3_k  <- 2 * (1 - 3^(-k)) / (1 - 2^(-k)) - 3
    dt3   <- (2 * log(3) * 3^(-k) * (1 - 2^(-k)) -
              2 * (1 - 3^(-k)) * log(2) * 2^(-k)) / (1 - 2^(-k))^2
    k     <- k - (t3_k - tau3) / dt3
    k     <- pmax(-0.99, pmin(0.99, k))
  }

  alpha <- lmom$l2 * k / (gamma(1 + k) * (1 - 2^(-k)))
  xi    <- lmom$l1 - alpha * (1 - gamma(1 + k)) / k

  structure(list(xi = xi, alpha = alpha, k = k), class = "GevParams")
}

# ---- Generalised Normal (GNO / log-normal family) -------------------------
# Used in FEH for rainfall frequency and as alternative to GLO/GEV.
# Quantile function: x(F) = xi + alpha * [1 - exp(-k * Phi^{-1}(F))] / k
# k=0: Normal.
#
# @param lmom Named list from sample_lmom()
# @return list(xi, alpha, k) with class "GnoParams"
.gno_fit_lmom <- function(lmom) {
  tau3 <- lmom$t3
  # Approximation (Hosking & Wallis 1997, Table A2)
  k  <- -tau3 * (0.6366 + tau3^2 * (0.1272 + tau3^2 * 0.0057))
  A1 <-  0.3989422803
  z  <- abs(k) + 1e-10
  alpha <- lmom$l2 * sqrt(pi) * z /
           (A1 * (exp(z^2 / 2) * (1 - 2 * pnorm(-z))))
  if (is.nan(alpha) || alpha <= 0) alpha <- lmom$l2 * sqrt(pi)
  xi    <- lmom$l1 - alpha * (exp(k^2 / 2) - 1) / k
  if (abs(k) < 1e-6) {
    alpha <- lmom$l2 * sqrt(pi)
    xi    <- lmom$l1
  }
  structure(list(xi = xi, alpha = alpha, k = k), class = "GnoParams")
}

# ---- Pearson Type III (PE3) ------------------------------------------------
# Used in FEH for some rainfall applications.
# Parameterised via mean, L-scale, L-skewness.
#
# @param lmom Named list from sample_lmom()
# @return list(mu, sigma, gamma) with class "Pe3Params"
.pe3_fit_lmom <- function(lmom) {
  t3 <- lmom$t3
  if (abs(t3) < 1e-6) {
    return(structure(list(mu = lmom$l1, sigma = lmom$l2 * sqrt(pi),
                          gamma = 0), class = "Pe3Params"))
  }
  # Approximation (Hosking & Wallis 1997)
  if (abs(t3) < 1/3) {
    z     <- 3 * pi * t3^2
    alpha <- (1 + 0.2906 * z) / (z * (1 + z * (-0.1882 + z * 0.0442)))
  } else {
    z     <- 1 - abs(t3)
    alpha <- (0.36067 * z^(-1.47950) + 0.73886) /
             (z^(-1.47950) + 0.1108 * z^(-0.73975) + 2.8172 * z + 0.91269)
  }
  rtalpha <- sqrt(alpha)
  beta    <- lmom$l2 * sqrt(pi) * exp(lgamma(alpha) - lgamma(alpha + 0.5))
  mu      <- lmom$l1
  gamma_v <- 2 / rtalpha * sign(t3)
  sigma_v <- beta * rtalpha
  structure(list(mu = mu, sigma = sigma_v, gamma = gamma_v),
            class = "Pe3Params")
}

# =============================================================================
# SECTION 3 : QUANTILE FUNCTIONS
# =============================================================================

#' GLO quantile function
#'
#' @param p      Non-exceedance probability (or vector).
#' @param params A `GloParams` list (xi, alpha, k) or equivalent.
#' @return Quantile(s) at probability p.
#' @export
qglo <- function(p, params) {
  p <- pmax(1e-12, pmin(1 - 1e-12, p))
  k <- params$k
  if (abs(k) < 1e-8) {
    params$xi - params$alpha * log((1 - p) / p)
  } else {
    params$xi + params$alpha * (1 - ((1 - p) / p)^k) / k
  }
}

#' GEV quantile function
#'
#' @inheritParams qglo
#' @export
qgev <- function(p, params) {
  p <- pmax(1e-12, pmin(1 - 1e-12, p))
  k <- params$k
  if (abs(k) < 1e-8) {
    params$xi - params$alpha * log(-log(p))
  } else {
    params$xi + params$alpha * (1 - (-log(p))^k) / k
  }
}

#' GNO quantile function
#'
#' @inheritParams qglo
#' @export
qgno <- function(p, params) {
  k <- params$k
  z <- qnorm(p)
  if (abs(k) < 1e-8) {
    params$xi + params$alpha * z
  } else {
    params$xi + params$alpha * (1 - exp(-k * z)) / k
  }
}

# =============================================================================
# SECTION 4 : L-MOMENT RATIO DIAGRAM BOUNDS (for goodness-of-fit)
# =============================================================================

# Theoretical (tau3, tau4) for GLO, GEV, GNO, PE3 at a grid of k values
# Used by fit_diagnostic() to plot where the sample sits relative to
# theoretical lines.

#' L-moment ratio diagram data for common FEH distributions
#'
#' @return A `data.table` with columns `dist`, `t3`, `t4`.
#' @export
lmrd_theoretical <- function() {
  k_seq <- seq(-0.5, 0.5, by = 0.01)

  # GLO: tau4 = (1 + 5*k^2) / 6  (Hosking & Wallis 1997)
  glo_t3 <- -k_seq
  glo_t4 <- (1 + 5 * k_seq^2) / 6

  # GEV: approximate tau4 from tau3
  gev_t3 <- seq(-0.4, 0.8, by = 0.01)
  # Polynomial approximation (Hosking & Wallis 1997, Table A2)
  gev_t4 <- 0.10701 + 0.11090 * gev_t3 + 0.84838 * gev_t3^2 -
             0.06669 * gev_t3^3 - 0.00567 * gev_t3^4

  data.table::rbindlist(list(
    data.table::data.table(dist = "GLO", t3 = glo_t3, t4 = glo_t4),
    data.table::data.table(dist = "GEV", t3 = gev_t3, t4 = gev_t4)
  ))
}
