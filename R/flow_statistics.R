# =============================================================================
# Tool:        reach.hydro — flow statistics
# Description: Flow statistics for the REACH framework: summary statistics,
#              baseflow index, flow duration curves, annual maxima extraction,
#              and peaks-over-threshold (POT) series. All functions accept
#              numeric vectors or data.tables and return data.tables.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-02-23 - JP: initial skeleton
# Tier:        2
# Inputs:      Numeric flow vectors [any consistent unit]; optional date vector
# Outputs:     data.table of statistics or series
# Dependencies: data.table, collapse
# TODO:        Implement baseflow separation (Lyne-Hollick, Boughton-Eckhardt)
# =============================================================================

#' Summary flow statistics
#'
#' Computes a standard set of flow statistics for a time series, including
#' mean, median, percentile flows (Q5, Q10, Q50, Q90, Q95), coefficient of
#' variation, and skewness.
#'
#' @param flow   Numeric vector of flow values \[any unit\].
#' @param na.rm  Logical. Remove NAs before computing. Default `TRUE`.
#'
#' @return A one-row `data.table` with columns:
#'   `n`, `n_valid`, `mean`, `median`, `sd`, `cv`,
#'   `Q5`, `Q10`, `Q50`, `Q90`, `Q95`, `min`, `max`.
#'
#' @examples
#' flow_stats(c(1.2, 3.4, 0.5, 8.1, NA, 2.2))
#'
#' @export
flow_stats <- function(flow, na.rm = TRUE) {
  checkmate::assert_numeric(flow)
  v <- if (na.rm) flow[!is.na(flow)] else flow

  data.table::data.table(
    n       = length(flow),
    n_valid = length(v),
    mean    = collapse::fmean(v),
    median  = collapse::fmedian(v),
    sd      = collapse::fsd(v),
    cv      = collapse::fsd(v) / max(collapse::fmean(v), 1e-12),
    Q5      = collapse::fnth(v, 0.05),
    Q10     = collapse::fnth(v, 0.10),
    Q50     = collapse::fnth(v, 0.50),
    Q90     = collapse::fnth(v, 0.90),
    Q95     = collapse::fnth(v, 0.95),
    min     = min(v),
    max     = max(v)
  )
}

#' Baseflow index
#'
#' Computes the baseflow index (BFI = mean baseflow / mean total flow) from
#' a simulated baseflow and total flow series. For observed-only series,
#' a digital filter-based baseflow separation is planned (see TODO).
#'
#' @param Qb   Numeric vector of baseflow \[same unit as Q\].
#' @param Q    Numeric vector of total flow.
#' @param na.rm Logical. Default `TRUE`.
#'
#' @return BFI scalar in \[0, 1\].
#'
#' @export
baseflow_index <- function(Qb, Q, na.rm = TRUE) {
  checkmate::assert_numeric(Qb)
  checkmate::assert_numeric(Q)
  if (length(Qb) != length(Q))
    stop("baseflow_index: `Qb` and `Q` must be the same length.", call. = FALSE)
  v <- !is.na(Qb) & !is.na(Q)
  mean_Q <- collapse::fmean(Q[v])
  if (mean_Q <= 0) return(NA_real_)
  collapse::fmean(Qb[v]) / mean_Q
}

#' Flow duration curve
#'
#' Computes the flow duration curve (FDC) by ranking flow values and computing
#' exceedance probabilities.
#'
#' @param flow   Numeric vector of flow values.
#' @param n_bins Integer. Number of exceedance probability bins. Default 100.
#' @param na.rm  Logical. Default `TRUE`.
#'
#' @return A `data.table` with columns `exceedance_prob` \[0, 1\] and `flow`.
#'
#' @export
flow_duration_curve <- function(flow, n_bins = 100L, na.rm = TRUE) {
  checkmate::assert_numeric(flow)
  v <- if (na.rm) flow[!is.na(flow)] else flow
  v <- sort(v, decreasing = TRUE)
  n <- length(v)
  probs <- seq_len(n) / (n + 1)
  # Interpolate to n_bins
  target_probs <- seq(0, 1, length.out = n_bins)
  flow_interp  <- stats::approx(probs, v, xout = target_probs,
                                rule = 2)$y
  data.table::data.table(
    exceedance_prob = target_probs,
    flow            = flow_interp
  )
}

#' Extract annual maxima
#'
#' Extracts the maximum flow value for each water year from a dated flow series.
#'
#' @param flow       Numeric vector of flow values.
#' @param dates      `Date` or `POSIXct` vector, same length as `flow`.
#' @param water_year_start Integer month that starts the water year. Default 10
#'                         (October, consistent with UK hydrological practice).
#'
#' @return A `data.table` with columns `water_year` (integer) and `annual_max`.
#'
#' @export
annual_maxima <- function(flow, dates, water_year_start = 10L) {
  checkmate::assert_numeric(flow)
  checkmate::assert(
    checkmate::check_date(dates),
    checkmate::check_posixct(dates)
  )
  if (length(flow) != length(dates))
    stop("annual_maxima: `flow` and `dates` must be the same length.", call. = FALSE)

  dt <- data.table::data.table(flow = flow, date = as.Date(dates))
  dt[, month := data.table::month(date)]
  dt[, year  := data.table::year(date)]
  dt[, water_year := data.table::fifelse(
    month >= water_year_start, year, year - 1L
  )]
  dt[!is.na(flow), .(annual_max = max(flow)), by = water_year][
    order(water_year)
  ]
}

#' Peaks-over-threshold series
#'
#' Extracts independent flood peaks exceeding a threshold, with a minimum
#' separation window to ensure independence between events.
#'
#' @param flow        Numeric vector of flow values.
#' @param dates       `Date` or `POSIXct` vector, same length as `flow`.
#' @param threshold   Flow threshold. Peaks must exceed this value.
#' @param min_sep     Integer. Minimum number of timesteps between independent
#'                    peaks. Default 3.
#'
#' @return A `data.table` with columns `date`, `peak_flow`, `peak_index`.
#'
#' @export
peaks_over_threshold <- function(flow, dates, threshold, min_sep = 3L) {
  checkmate::assert_numeric(flow)
  checkmate::assert_number(threshold)
  if (length(flow) != length(dates))
    stop("peaks_over_threshold: `flow` and `dates` must be the same length.",
         call. = FALSE)

  n   <- length(flow)
  idx <- which(flow > threshold)
  if (length(idx) == 0L)
    return(data.table::data.table(date = as.Date(character(0)),
                                  peak_flow = numeric(0),
                                  peak_index = integer(0)))

  # Decluster: keep only the highest peak within each min_sep window
  peaks      <- integer(0)
  last_peak  <- -Inf

  for (i in idx) {
    if (i - last_peak >= min_sep) {
      # Find local maximum in the window around i
      win_start <- max(1L, i - min_sep + 1L)
      win_end   <- min(n,  i + min_sep - 1L)
      local_max <- win_start - 1L + which.max(flow[win_start:win_end])
      if (!(local_max %in% peaks)) {
        peaks     <- c(peaks, local_max)
        last_peak <- local_max
      }
    }
  }

  data.table::data.table(
    date       = as.Date(dates[peaks]),
    peak_flow  = flow[peaks],
    peak_index = peaks
  )[order(date)]
}
