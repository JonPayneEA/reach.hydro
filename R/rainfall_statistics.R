# =============================================================================
# Tool:        reach.hydro -- rainfall statistics
# Description: Rainfall analysis tools for the REACH framework: storm event
#              extraction, antecedent precipitation index (API), and empirical
#              intensity-duration-frequency (IDF) analysis.  All functions accept
#              numeric vectors and optional date/datetime vectors; they work at
#              any timestep (daily, hourly, 15-minute).
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-04-04
# Tier:        2
# Inputs:      Numeric rainfall vectors [mm/timestep]; optional date vector
# Outputs:     scalars, numeric vectors, or data.tables
# Dependencies: data.table, checkmate
# =============================================================================

#' Rainfall event extraction
#'
#' Identifies discrete storm events from a continuous rainfall series. Events
#' are separated by dry gaps of at least `min_dry` consecutive zero (or
#' near-zero) timesteps.  Each event is characterised by its duration, total
#' depth, peak intensity, and mean intensity.
#'
#' Duration is in timesteps when `dates` is `NULL`, in days when `dates` is a
#' `Date` vector, and in hours when `dates` is a `POSIXct` vector.
#'
#' @param rain      Numeric vector of rainfall depths per timestep. Must be
#'                  non-negative. Can also be a `reach.io` `FlodeRainfall_Daily` or
#'                  `FlodeRainfall_15min` FlodeHydroData object, in which case values and
#'                  datetimes are extracted automatically and `dates` is ignored.
#' @param dates     Optional `Date` or `POSIXct` vector the same length as
#'                  `rain`.  If supplied, `start_date` and `end_date` columns
#'                  are added to the output and duration is in real time units.
#'                  Ignored when `rain` is a FlodeHydroData object.
#' @param min_dry   Minimum number of consecutive dry timesteps (i.e. rain
#'                  below `dry_thresh`) required to separate two events.
#'                  Default `6L`.
#' @param min_depth Minimum total event depth to include in the output.  Events
#'                  with `total_depth < min_depth` are dropped.  Default `0`.
#' @param dry_thresh Threshold below which a timestep is considered dry.
#'                  Default `0.0` (strictly zero).
#'
#' @return A `data.table` with one row per event and columns:
#' \describe{
#'   \item{`start_idx`}{Integer index of the first wet timestep.}
#'   \item{`end_idx`}{Integer index of the last wet timestep.}
#'   \item{`duration`}{Event duration in timesteps (or days/hours if `dates`
#'     supplied).}
#'   \item{`total_depth`}{Total accumulated rainfall depth.}
#'   \item{`max_intensity`}{Maximum rainfall in any single timestep.}
#'   \item{`mean_intensity`}{Mean rainfall per wet timestep.}
#'   \item{`start_date`, `end_date`}{Present when `dates` is supplied.}
#' }
#' Returns a zero-row table if no qualifying events are found.
#'
#' @examples
#' rain <- c(0, 0, 2, 5, 3, 0, 0, 0, 0, 0, 0, 1, 4, 2, 0, 0)
#' rainfall_events(rain, min_dry = 4L)
#'
#' @export
rainfall_events <- function(rain, dates = NULL, min_dry = 6L,
                            min_depth = 0, dry_thresh = 0.0) {
  if (.is_hydrodata(rain)) {
    .assert_rainfall(rain, "rain")
    dates <- hydrodata_datetimes(rain)
    rain  <- hydrodata_values(rain)
  }
  checkmate::assert_numeric(rain, lower = 0, any.missing = FALSE, min.len = 1L)
  checkmate::assert_int(min_dry, lower = 1L)
  checkmate::assert_number(min_depth, lower = 0)
  checkmate::assert_number(dry_thresh, lower = 0)

  has_dates  <- !is.null(dates)
  is_posixct <- FALSE
  if (has_dates) {
    checkmate::assert(
      checkmate::check_date(dates),
      checkmate::check_posixct(dates)
    )
    if (length(dates) != length(rain))
      stop("rainfall_events: `rain` and `dates` must be the same length.",
           call. = FALSE)
    is_posixct <- inherits(dates, "POSIXct")
  }

  wet <- rain > dry_thresh
  n   <- length(rain)

  # Identify wet runs with gaps no larger than min_dry
  # Strategy: expand wet periods to bridge short dry gaps, then extract runs.
  expanded <- wet
  for (i in seq_len(n)) {
    if (!wet[i]) next
    # Look forward: if there is another wet timestep within min_dry, fill the gap
    j_end <- min(n, i + min_dry)
    if (any(wet[(i + 1L):j_end])) {
      expanded[i:(min(which(wet[(i + 1L):j_end]) + i, n))] <- TRUE
    }
  }

  # Extract contiguous wet runs from expanded
  starts  <- integer(0)
  ends    <- integer(0)
  in_evt  <- FALSE

  for (i in seq_len(n)) {
    if (expanded[i] && !in_evt) {
      starts  <- c(starts, i)
      in_evt  <- TRUE
    } else if (!expanded[i] && in_evt) {
      ends    <- c(ends, i - 1L)
      in_evt  <- FALSE
    }
  }
  if (in_evt) ends <- c(ends, n)

  empty_dt <- data.table::data.table(
    start_idx      = integer(0),
    end_idx        = integer(0),
    duration       = numeric(0),
    total_depth    = numeric(0),
    max_intensity  = numeric(0),
    mean_intensity = numeric(0)
  )
  if (length(starts) == 0L) return(empty_dt)

  dur  <- numeric(length(starts))
  tot  <- numeric(length(starts))
  mx   <- numeric(length(starts))
  mnI  <- numeric(length(starts))

  for (j in seq_along(starts)) {
    s  <- starts[j]
    e  <- ends[j]
    rr <- rain[s:e]
    tot[j] <- sum(rr)
    mx[j]  <- max(rr)
    n_wet  <- sum(rr > dry_thresh)
    mnI[j] <- if (n_wet > 0) tot[j] / n_wet else 0

    if (has_dates) {
      d1 <- as.POSIXct(dates[s])
      d2 <- as.POSIXct(dates[e])
      dt_units <- if (is_posixct) "hours" else "days"
      dur[j] <- as.numeric(difftime(d2, d1, units = dt_units)) + 1
    } else {
      dur[j] <- e - s + 1L
    }
  }

  out <- data.table::data.table(
    start_idx      = starts,
    end_idx        = ends,
    duration       = dur,
    total_depth    = tot,
    max_intensity  = mx,
    mean_intensity = mnI
  )

  if (has_dates) {
    out[, start_date := as.Date(dates[start_idx])]
    out[, end_date   := as.Date(dates[end_idx])]
  }

  out <- out[total_depth >= min_depth]
  out
}

#' Antecedent Precipitation Index (API)
#'
#' Computes the classic exponential-decay antecedent precipitation index:
#'
#' \deqn{API_t = k \cdot API_{t-1} + P_t}
#'
#' A high API indicates wet antecedent conditions; a low value indicates dry.
#' The decay factor `k` is interpreted per-timestep -- for daily data a typical
#' value is 0.85-0.95; for hourly data a smaller value (e.g. 0.99) may be more
#' appropriate to maintain the same effective memory window.
#'
#' @param rain  Numeric vector of rainfall depths per timestep. Must be
#'              non-negative. Can also be a `reach.io` `FlodeRainfall_Daily` or
#'              `FlodeRainfall_15min` FlodeHydroData object, in which case values are
#'              extracted automatically.
#' @param k     Decay factor per timestep, in (0, 1). Default `0.9`.
#'
#' @return Numeric vector of the same length as `rain`, giving the API at each
#'   timestep.
#'
#' @examples
#' rain <- c(0, 0, 10, 5, 0, 0, 0, 3, 0, 0)
#' api(rain, k = 0.9)
#'
#' @export
api <- function(rain, k = 0.9) {
  if (.is_hydrodata(rain)) {
    .assert_rainfall(rain, "rain")
    rain <- hydrodata_values(rain)
  }
  checkmate::assert_numeric(rain, lower = 0, any.missing = FALSE, min.len = 1L)
  checkmate::assert_number(k, lower = 0, upper = 1, finite = TRUE)

  n   <- length(rain)
  out <- numeric(n)
  out[1L] <- rain[1L]
  for (i in seq(2L, n)) {
    out[i] <- k * out[i - 1L] + rain[i]
  }
  out
}

#' Empirical intensity-duration-frequency table
#'
#' Extracts the maximum accumulated rainfall depth for each specified duration
#' from a continuous series, returning annual values when `dates` is supplied
#' or overall maxima otherwise.  This is a data-exploration tool; it does not
#' fit a statistical distribution (for fitted DDF see `feh_ddf()`).
#'
#' Durations must be specified in timestep units.  For example, with 15-minute
#' data pass `durations = c(4, 12, 24, 48, 96)` to get 1-hour, 3-hour, 6-hour,
#' 12-hour, and 24-hour maxima.  With daily data pass `durations = c(1, 3, 6)`
#' for 1-, 3-, and 6-day maxima.
#'
#' @param rain              Numeric vector of rainfall per timestep. Non-negative.
#'                          Can also be a `reach.io` `FlodeRainfall_Daily` or
#'                          `FlodeRainfall_15min` FlodeHydroData object, in which case
#'                          values and datetimes are extracted automatically and
#'                          `dates` is ignored.
#' @param durations         Integer vector of window widths (in timesteps) for
#'                          which to compute maxima. Default
#'                          `c(1L, 3L, 6L, 12L, 24L, 48L)`.
#' @param dates             Optional `Date` or `POSIXct` vector the same length
#'                          as `rain`.  If supplied, annual maxima are returned
#'                          grouped by water year. Ignored when `rain` is a
#'                          FlodeHydroData object.
#' @param water_year_start  Integer month that starts the water year. Default
#'                          `10L` (October, UK convention).
#'
#' @return A `data.table` with columns `duration` (integer, in timestep units)
#'   and `amax_depth`.  When `dates` is supplied, a `water_year` column is
#'   added and there is one row per `(duration, water_year)` combination.
#'   Rows are sorted by `duration` (then `water_year` if present).
#'
#' @examples
#' set.seed(42)
#' rain  <- pmax(0, rnorm(365, 2, 3))
#' idf_empirical(rain, durations = c(1L, 3L, 7L))
#'
#' @export
idf_empirical <- function(rain,
                          durations = c(1L, 3L, 6L, 12L, 24L, 48L),
                          dates = NULL,
                          water_year_start = 10L) {
  if (.is_hydrodata(rain)) {
    .assert_rainfall(rain, "rain")
    dates <- hydrodata_datetimes(rain)
    rain  <- hydrodata_values(rain)
  }
  checkmate::assert_numeric(rain, lower = 0, any.missing = FALSE, min.len = 1L)
  checkmate::assert_integerish(durations, lower = 1L, min.len = 1L)
  durations <- as.integer(durations)

  has_dates  <- !is.null(dates)
  if (has_dates) {
    checkmate::assert(
      checkmate::check_date(dates),
      checkmate::check_posixct(dates)
    )
    if (length(dates) != length(rain))
      stop("idf_empirical: `rain` and `dates` must be the same length.",
           call. = FALSE)
  }

  n   <- length(rain)

  # Rolling accumulation for each duration
  rows <- vector("list", length(durations))

  for (di in seq_along(durations)) {
    d    <- durations[di]
    acc  <- vapply(seq_len(n), function(i) {
      sum(rain[max(1L, i - d + 1L):i])
    }, numeric(1L))

    if (!has_dates) {
      rows[[di]] <- data.table::data.table(
        duration   = d,
        amax_depth = max(acc, na.rm = TRUE)
      )
    } else {
      dt <- data.table::data.table(acc = acc, date = as.Date(dates))
      dt[, month      := data.table::month(date)]
      dt[, year       := data.table::year(date)]
      dt[, water_year := data.table::fifelse(
        month >= as.integer(water_year_start), year, year - 1L
      )]
      wy_max <- dt[, .(amax_depth = max(acc, na.rm = TRUE)), by = water_year]
      wy_max[, duration := d]
      rows[[di]] <- wy_max[, .(duration, water_year, amax_depth)]
    }
  }

  out <- data.table::rbindlist(rows)
  if (has_dates) {
    data.table::setorder(out, duration, water_year)
  } else {
    data.table::setorder(out, duration)
  }
  out
}
