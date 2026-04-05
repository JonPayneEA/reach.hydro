# =============================================================================
# Tool:        reach.hydro — flow statistics
# Description: Flow statistics for the REACH framework: summary statistics,
#              baseflow index, flow duration curves, annual maxima extraction,
#              peaks-over-threshold (POT) series, baseflow separation, n-day
#              flow statistics, monthly statistics, flow deficit analysis, and
#              recession analysis. All functions accept numeric vectors and
#              return scalars, lists, or data.tables.
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# Modified:    2026-04-04 - add baseflow_separate, q_n_day, monthly_flow_stats,
#              flow_deficit, flow_recession; all functions work at any timestep
# Tier:        2
# Inputs:      Numeric flow vectors [any consistent unit]; optional date vector
# Outputs:     data.table of statistics or series
# Dependencies: data.table, collapse
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

# =============================================================================
# Baseflow separation
# =============================================================================

#' Baseflow separation using a digital filter
#'
#' Separates a flow series into baseflow and quickflow components using one of
#' two standard recursive digital filter methods. The filter can be applied in
#' multiple passes (forward/backward) to reduce phase shift.
#'
#' Filter coefficients are applied per-timestep and are consistent with the
#' literature for the specified timestep. For sub-daily data (e.g. 15-minute),
#' typical daily-calibrated alpha values (e.g. 0.925) will retain more
#' quickflow per step; recalibrate accordingly for sub-daily use.
#'
#' @param flow   Numeric vector of total flow values. Must be non-negative.
#'               Can also be a `reach.io` `Flow_Daily` or `Flow_15min` HydroData
#'               object, in which case values and datetimes are extracted
#'               automatically.
#' @param method One of `"lyne_hollick"` (default) or `"boughton_eckhardt"`.
#' @param alpha  Lyne-Hollick filter parameter in (0, 1). Default `0.925`.
#' @param k      Boughton-Eckhardt recession constant in (0, 1). Default `0.975`.
#' @param C      Boughton-Eckhardt direct runoff parameter in (0, 1). Default `0.1`.
#' @param passes Integer. Number of filter passes (alternating forward/backward).
#'               Default `3L` (forward-backward-forward).
#'
#' @return A list with:
#' \describe{
#'   \item{`baseflow`}{Numeric vector of separated baseflow, same length as `flow`.}
#'   \item{`quickflow`}{Numeric vector of quickflow (`flow - baseflow`).}
#'   \item{`bfi`}{Scalar baseflow index: `mean(baseflow) / mean(flow)`.}
#' }
#'
#' @examples
#' Q  <- c(2, 3, 8, 15, 10, 6, 4, 3, 2.5, 2)
#' bf <- baseflow_separate(Q, method = "lyne_hollick")
#' bf$bfi
#'
#' @export
baseflow_separate <- function(flow,
                              method = c("lyne_hollick", "boughton_eckhardt"),
                              alpha = 0.925, k = 0.975, C = 0.1,
                              passes = 3L) {
  method <- match.arg(method)
  if (.is_hydrodata(flow)) {
    .assert_flow(flow, "flow")
    flow <- hydrodata_values(flow)
  }
  checkmate::assert_numeric(flow, lower = 0, any.missing = FALSE, min.len = 2L)
  checkmate::assert_int(passes, lower = 1L)

  n <- length(flow)

  .lh_pass <- function(Q, forward) {
    if (!forward) Q <- rev(Q)
    qf <- numeric(n)
    qf[1L] <- 0
    for (i in seq(2L, n)) {
      qf[i] <- alpha * qf[i - 1L] +
        (1 + alpha) / 2 * (Q[i] - Q[i - 1L])
    }
    qb <- pmax(0, Q - pmax(0, qf))
    if (!forward) qb <- rev(qb)
    qb
  }

  .be_pass <- function(Q, forward) {
    if (!forward) Q <- rev(Q)
    qb <- numeric(n)
    qb[1L] <- Q[1L]
    a1 <- k / (1 - C)
    a2 <- C / (1 - C)
    for (i in seq(2L, n)) {
      qb[i] <- min(a1 * qb[i - 1L] + a2 * Q[i], Q[i])
    }
    if (!forward) qb <- rev(qb)
    qb
  }

  filter_fn <- if (method == "lyne_hollick") .lh_pass else .be_pass

  qb <- flow
  for (p in seq_len(passes)) {
    forward <- (p %% 2L == 1L)
    qb <- filter_fn(qb, forward)
  }
  qb <- pmin(qb, flow)

  mean_Q <- mean(flow)
  bfi    <- if (mean_Q > 0) mean(qb) / mean_Q else NA_real_

  list(
    baseflow  = qb,
    quickflow = flow - qb,
    bfi       = bfi
  )
}

# =============================================================================
# N-day flow statistics
# =============================================================================

#' N-day minimum or maximum flow
#'
#' Computes the rolling n-step minimum (or maximum) of a flow series, then
#' returns either the overall scalar or annual values grouped by water year.
#'
#' This is the standard Q7 dry-weather-flow index when `n = 7` and
#' `type = "min"`. For sub-daily data the rolling window is in timesteps, so
#' pass `n` in the appropriate units (e.g. `n = 7 * 96` for 7-day Q using
#' 15-minute data).
#'
#' @param flow              Numeric vector of flow values. Can also be a
#'                          `reach.io` `Flow_Daily` or `Flow_15min` HydroData
#'                          object, in which case values and datetimes are
#'                          extracted automatically and `dates` is ignored.
#' @param n                 Integer. Rolling window width in timesteps. Default `7L`.
#' @param type              `"min"` (default) or `"max"`.
#' @param dates             Optional `Date` or `POSIXct` vector the same length as
#'                          `flow`. If supplied, annual values are returned.
#'                          Ignored when `flow` is a HydroData object.
#' @param water_year_start  Integer month that starts the water year. Default `10L`.
#'
#' @return If `dates` is `NULL`: a scalar (overall rolling min/max).
#'   If `dates` is supplied: a `data.table` with columns `water_year` and
#'   `q_n_day`.
#'
#' @examples
#' Q <- c(5, 4, 3, 2, 1.5, 2, 3, 4, 5, 6)
#' q_n_day(Q, n = 3, type = "min")
#'
#' @export
q_n_day <- function(flow, n = 7L, type = c("min", "max"),
                    dates = NULL, water_year_start = 10L) {
  type <- match.arg(type)
  if (.is_hydrodata(flow)) {
    .assert_flow(flow, "flow")
    dates <- hydrodata_datetimes(flow)
    flow  <- hydrodata_values(flow)
  }
  checkmate::assert_numeric(flow, min.len = 1L)
  checkmate::assert_int(n, lower = 1L)

  nn  <- length(flow)
  fun <- if (type == "min") min else max

  # Rolling window: for each position i return fun over [i-n+1, i]
  roll <- vapply(seq_len(nn), function(i) {
    w <- flow[max(1L, i - n + 1L):i]
    if (all(is.na(w))) NA_real_ else fun(w, na.rm = TRUE)
  }, numeric(1L))

  if (is.null(dates)) {
    return(fun(roll, na.rm = TRUE))
  }

  checkmate::assert(
    checkmate::check_date(dates),
    checkmate::check_posixct(dates)
  )
  if (length(dates) != nn)
    stop("q_n_day: `flow` and `dates` must be the same length.", call. = FALSE)

  dt <- data.table::data.table(roll = roll, date = as.Date(dates))
  dt[, month      := data.table::month(date)]
  dt[, year       := data.table::year(date)]
  dt[, water_year := data.table::fifelse(
    month >= as.integer(water_year_start), year, year - 1L
  )]
  dt[!is.na(roll), .(q_n_day = fun(roll, na.rm = TRUE)), by = water_year][
    order(water_year)
  ]
}

# =============================================================================
# Monthly/seasonal flow statistics
# =============================================================================

#' Monthly flow statistics
#'
#' Computes summary statistics by calendar month from a dated flow series.
#' Works with any timestep (daily, hourly, 15-minute) provided a `dates` vector
#' is supplied.
#'
#' @param flow   Numeric vector of flow values. Can also be a `reach.io`
#'               `Flow_Daily` or `Flow_15min` HydroData object, in which case
#'               values and datetimes are extracted automatically and `dates`
#'               is ignored.
#' @param dates  `Date` or `POSIXct` vector, same length as `flow`. Ignored
#'               when `flow` is a HydroData object.
#'
#' @return A 12-row `data.table` with columns:
#'   `month` (1-12), `mean`, `median`, `Q10`, `Q90`, `max`.
#'   Rows are ordered by month number.
#'
#' @examples
#' set.seed(1)
#' dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 365)
#' flow  <- pmax(0, rnorm(365, 3, 1))
#' monthly_flow_stats(flow, dates)
#'
#' @export
monthly_flow_stats <- function(flow, dates = NULL) {
  if (.is_hydrodata(flow)) {
    .assert_flow(flow, "flow")
    dates <- hydrodata_datetimes(flow)
    flow  <- hydrodata_values(flow)
  }
  checkmate::assert_numeric(flow)
  if (is.null(dates))
    stop("monthly_flow_stats: `dates` is required when `flow` is a numeric vector.",
         call. = FALSE)
  checkmate::assert(
    checkmate::check_date(dates),
    checkmate::check_posixct(dates)
  )
  if (length(flow) != length(dates))
    stop("monthly_flow_stats: `flow` and `dates` must be the same length.",
         call. = FALSE)

  dt <- data.table::data.table(flow = flow, date = as.Date(dates))
  dt[, month := data.table::month(date)]

  dt[!is.na(flow),
    .(
      mean   = collapse::fmean(flow),
      median = collapse::fmedian(flow),
      Q10    = collapse::fnth(flow, 0.10),
      Q90    = collapse::fnth(flow, 0.90),
      max    = max(flow)
    ),
    by = month
  ][order(month)]
}

# =============================================================================
# Flow deficit analysis
# =============================================================================

#' Flow deficit analysis
#'
#' Identifies spells where flow falls below a threshold and characterises each
#' spell by its duration, total deficit volume, and maximum instantaneous
#' deficit.
#'
#' Duration is in timesteps if `dates` is `NULL`, or in days when `dates` is a
#' `Date` vector, or in hours when `dates` is `POSIXct`.
#'
#' @param flow       Numeric vector of flow values. Can also be a `reach.io`
#'                   `Flow_Daily` or `Flow_15min` HydroData object, in which
#'                   case values and datetimes are extracted automatically and
#'                   `dates` is ignored.
#' @param threshold  Flow threshold. Spells where `flow < threshold` are
#'                   identified.
#' @param dates      Optional `Date` or `POSIXct` vector the same length as
#'                   `flow`. If supplied, `start_date` and `end_date` columns
#'                   are added and duration is in real time units. Ignored when
#'                   `flow` is a HydroData object.
#'
#' @return A `data.table` with one row per deficit spell and columns:
#'   `start_idx`, `end_idx`, `duration`, `deficit_volume`, `max_deficit`.
#'   When `dates` is supplied: `start_date` and `end_date` columns are added.
#'   Returns a zero-row table if no deficit spells are found.
#'
#' @examples
#' Q <- c(5, 4, 2, 1, 0.5, 1.5, 3, 5, 4, 2, 0.8, 1.2, 6)
#' flow_deficit(Q, threshold = 3)
#'
#' @export
flow_deficit <- function(flow, threshold, dates = NULL) {
  if (.is_hydrodata(flow)) {
    .assert_flow(flow, "flow")
    dates <- hydrodata_datetimes(flow)
    flow  <- hydrodata_values(flow)
  }
  checkmate::assert_numeric(flow, min.len = 1L)
  checkmate::assert_number(threshold)

  has_dates  <- !is.null(dates)
  is_posixct <- FALSE
  if (has_dates) {
    checkmate::assert(
      checkmate::check_date(dates),
      checkmate::check_posixct(dates)
    )
    if (length(flow) != length(dates))
      stop("flow_deficit: `flow` and `dates` must be the same length.",
           call. = FALSE)
    is_posixct <- inherits(dates, "POSIXct")
  }

  n      <- length(flow)
  below  <- !is.na(flow) & flow < threshold
  starts <- integer(0)
  ends   <- integer(0)
  in_spell <- FALSE

  for (i in seq_len(n)) {
    if (below[i] && !in_spell) {
      starts   <- c(starts, i)
      in_spell <- TRUE
    } else if (!below[i] && in_spell) {
      ends     <- c(ends, i - 1L)
      in_spell <- FALSE
    }
  }
  if (in_spell) ends <- c(ends, n)

  if (length(starts) == 0L) {
    base <- data.table::data.table(
      start_idx      = integer(0),
      end_idx        = integer(0),
      duration       = numeric(0),
      deficit_volume = numeric(0),
      max_deficit    = numeric(0)
    )
    if (has_dates) {
      base[, start_date := as.Date(character(0))]
      base[, end_date   := as.Date(character(0))]
    }
    return(base)
  }

  duration <- numeric(length(starts))
  deficit  <- numeric(length(starts))
  maxdef   <- numeric(length(starts))

  for (j in seq_along(starts)) {
    s  <- starts[j]
    e  <- ends[j]
    dv <- threshold - flow[s:e]
    if (has_dates) {
      d1 <- as.POSIXct(dates[s])
      d2 <- as.POSIXct(dates[e])
      dt_units  <- if (is_posixct) "hours" else "days"
      duration[j] <- as.numeric(difftime(d2, d1, units = dt_units)) + 1
    } else {
      duration[j] <- e - s + 1L
    }
    deficit[j] <- sum(dv, na.rm = TRUE)
    maxdef[j]  <- max(dv,  na.rm = TRUE)
  }

  out <- data.table::data.table(
    start_idx      = starts,
    end_idx        = ends,
    duration       = duration,
    deficit_volume = deficit,
    max_deficit    = maxdef
  )

  if (has_dates) {
    out[, start_date := as.Date(dates[start_idx])]
    out[, end_date   := as.Date(dates[end_idx])]
  }

  out
}

# =============================================================================
# Recession analysis
# =============================================================================

#' Flow recession analysis
#'
#' Identifies recession limbs from a flow series and fits an exponential decay
#' constant `k` to each: `Q(t) = Q0 * exp(-t/k)`. Returns the median `k`
#' across all identified recessions and per-event statistics.
#'
#' A recession is identified as a run of at least `min_duration` consecutive
#' timesteps where `Q[t] / Q[t-1] <= min_ratio` (i.e. flow is non-increasing
#' within the tolerance set by `min_ratio`). The fitted `k` is in the same
#' time units as the timestep.
#'
#' @param flow          Numeric vector of flow values. Must be positive. Can
#'                      also be a `reach.io` `Flow_Daily` or `Flow_15min`
#'                      HydroData object, in which case values and datetimes
#'                      are extracted automatically and `dates` is ignored.
#' @param dates         Optional `Date` or `POSIXct` vector, same length as
#'                      `flow`. If supplied, `start_date` and `end_date` columns
#'                      are added to the per-event table. Ignored when `flow` is
#'                      a HydroData object.
#' @param min_duration  Minimum number of consecutive timesteps to qualify as a
#'                      recession. Default `5L`.
#' @param min_ratio     Maximum allowed `Q[t] / Q[t-1]` ratio for a timestep to
#'                      be classified as receding. Default `0.95`.
#'
#' @return A list with:
#' \describe{
#'   \item{`k`}{Median recession constant across all identified events (same
#'     time units as the input timestep). `NA` if no recessions found.}
#'   \item{`n_recessions`}{Integer count of identified recession events.}
#'   \item{`recessions`}{`data.table` with one row per event: `start_idx`,
#'     `end_idx`, `Q0`, `k_fit`. If `dates` supplied: `start_date`, `end_date`.}
#' }
#'
#' @examples
#' Q <- c(20, 15, 11, 8, 6, 5, 4, 3.5, 3, 2.8, 2.6, 2.5)
#' flow_recession(Q)$k
#'
#' @export
flow_recession <- function(flow, dates = NULL,
                           min_duration = 5L, min_ratio = 0.95) {
  if (.is_hydrodata(flow)) {
    .assert_flow(flow, "flow")
    dates <- hydrodata_datetimes(flow)
    flow  <- hydrodata_values(flow)
  }
  checkmate::assert_numeric(flow, lower = 0, min.len = 2L)
  checkmate::assert_int(min_duration, lower = 2L)
  checkmate::assert_number(min_ratio, lower = 0, upper = 1)

  has_dates <- !is.null(dates)
  if (has_dates) {
    checkmate::assert(
      checkmate::check_date(dates),
      checkmate::check_posixct(dates)
    )
    if (length(dates) != length(flow))
      stop("flow_recession: `flow` and `dates` must be the same length.",
           call. = FALSE)
  }

  n     <- length(flow)
  ratio <- c(NA_real_, flow[2:n] / pmax(flow[1:(n - 1L)], 1e-12))
  is_rec <- !is.na(ratio) & ratio <= min_ratio

  # Identify runs of recession timesteps
  starts <- integer(0)
  ends   <- integer(0)
  in_r   <- FALSE

  for (i in seq_len(n)) {
    if (is_rec[i] && !in_r) {
      starts <- c(starts, i - 1L)  # include the step before decline starts
      in_r   <- TRUE
    } else if (!is_rec[i] && in_r) {
      ends <- c(ends, i - 1L)
      in_r <- FALSE
    }
  }
  if (in_r) ends <- c(ends, n)

  # Keep only events meeting min_duration
  keep <- (ends - starts + 1L) >= min_duration
  starts <- starts[keep]
  ends   <- ends[keep]

  empty_dt <- data.table::data.table(
    start_idx = integer(0), end_idx = integer(0),
    Q0 = numeric(0), k_fit = numeric(0)
  )
  if (length(starts) == 0L)
    return(list(k = NA_real_, n_recessions = 0L, recessions = empty_dt))

  k_fits <- numeric(length(starts))
  Q0s    <- numeric(length(starts))

  for (j in seq_along(starts)) {
    s  <- starts[j]
    e  <- ends[j]
    Q0s[j] <- flow[s]
    seg    <- flow[s:e]
    t_vec  <- seq(0, e - s)
    # Fit log(Q) ~ t; k = -1/slope
    pos    <- seg > 0
    if (sum(pos) < 2L) {
      k_fits[j] <- NA_real_
      next
    }
    lm_fit   <- stats::lm.fit(cbind(1, t_vec[pos]), log(seg[pos]))
    slope    <- lm_fit$coefficients[2L]
    k_fits[j] <- if (!is.na(slope) && slope < 0) -1 / slope else NA_real_
  }

  out <- data.table::data.table(
    start_idx = starts,
    end_idx   = ends,
    Q0        = Q0s,
    k_fit     = k_fits
  )
  if (has_dates) {
    out[, start_date := as.Date(dates[start_idx])]
    out[, end_date   := as.Date(dates[end_idx])]
  }

  k_med <- stats::median(k_fits, na.rm = TRUE)

  list(
    k            = k_med,
    n_recessions = length(starts),
    recessions   = out
  )
}
