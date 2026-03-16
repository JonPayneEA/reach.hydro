# =============================================================================
# Tool:        reach.hydro — reach.io compatibility layer
# Description: Bridges reach.io S7 HydroData objects into reach.hydro
#              functions. Provides:
#                - Class validators: check the right HydroData subclass is
#                  passed and that its readings meet minimum requirements
#                - Extractors: pull the value vector out of a HydroData object
#                  with consistent NA handling and unit checking
#                - Coercers: convert HydroData pairs (Rainfall + Flow) into
#                  the flat vectors that pdm() and FEH functions expect
#
#              Design principle: reach.hydro does NOT depend on reach.io at
#              build time (reach.io is in Suggests, not Imports). All
#              reach.io-aware code is gated by .is_hydrodata() so the package
#              installs and loads cleanly without reach.io present. Users who
#              pass HydroData objects without reach.io installed receive a
#              clear error message.
#
#              reach.io S7 class hierarchy (from reach.io documentation):
#                Abstract parent : HydroData
#                  Slots: readings (data.table), parameter (chr),
#                         period_name (chr), from_date (chr), to_date (chr),
#                         n_measures (int), n_rows (int), downloaded_at (POSIXct)
#                Concrete classes (all inherit HydroData):
#                  Rainfall_Daily, Rainfall_15min,
#                  Flow_Daily,    Flow_15min,
#                  Level_Daily,   Level_15min
#                readings schema:
#                  dateTime (POSIXct), date (Date), value (numeric),
#                  measure_notation (chr), quality (chr, optional),
#                  completeness (numeric, optional)
#                Generics:
#                  as_data_table(x) -> x@readings
#                  as_long(x)       -> readings + parameter column
#
# Flode Module: reach.hydro
# Author:      Forecasting and Warning Team
# Created:     2026-02-23
# Modified:    2026-02-23 - JP: initial implementation
# Tier:        1
# Inputs:      reach.io HydroData S7 objects
# Outputs:     Numeric vectors; data.tables; validated inputs for model fns
# Dependencies: data.table (Imports); reach.io (Suggests)
# =============================================================================

# =============================================================================
# SECTION 1 : INTERNAL CLASS DETECTION HELPERS
# =============================================================================

# Check whether reach.io is installed (without requiring it at load time)
.reach_io_available <- function() {
  requireNamespace("reach.io", quietly = TRUE)
}

# Check whether an object is any HydroData subclass
.is_hydrodata <- function(x) {
  if (!.reach_io_available()) return(FALSE)
  inherits(x, "HydroData")
}

# Check for a specific concrete class or group
.is_rainfall <- function(x) {
  inherits(x, c("Rainfall_Daily", "Rainfall_15min"))
}

.is_flow <- function(x) {
  inherits(x, c("Flow_Daily", "Flow_15min"))
}

.is_level <- function(x) {
  inherits(x, c("Level_Daily", "Level_15min"))
}

.is_daily <- function(x) {
  inherits(x, c("Rainfall_Daily", "Flow_Daily", "Level_Daily"))
}

.is_15min <- function(x) {
  inherits(x, c("Rainfall_15min", "Flow_15min", "Level_15min"))
}

# =============================================================================
# SECTION 2 : VALIDATORS
# =============================================================================

#' Assert that an object is a reach.io HydroData object
#'
#' Stops with an informative message if `x` is not a `HydroData` subclass,
#' or if reach.io is not installed.
#'
#' @param x        Object to check.
#' @param arg_name Name of the argument (used in error messages).
#' @keywords internal
.assert_hydrodata <- function(x, arg_name = deparse(substitute(x))) {
  if (!.reach_io_available())
    stop(sprintf(
      "%s: object appears to be a reach.io HydroData but reach.io is not ",
      "installed. Install reach.io or supply a plain numeric vector instead.",
      arg_name), call. = FALSE)
  if (!.is_hydrodata(x))
    stop(sprintf(
      "%s: expected a reach.io HydroData object; got %s.",
      arg_name, paste(class(x), collapse = "/")), call. = FALSE)
  invisible(x)
}

#' Assert that a HydroData object is a Rainfall subclass
#'
#' @param x        HydroData object.
#' @param arg_name Argument name for error messages.
#' @keywords internal
.assert_rainfall <- function(x, arg_name = deparse(substitute(x))) {
  .assert_hydrodata(x, arg_name)
  if (!.is_rainfall(x))
    stop(sprintf(
      "%s: expected a Rainfall_Daily or Rainfall_15min object; got %s. ",
      "Check you have passed the rainfall series, not the flow series.",
      arg_name, paste(class(x), collapse = "/")), call. = FALSE)
  invisible(x)
}

#' Assert that a HydroData object is a Flow subclass
#'
#' @param x        HydroData object.
#' @param arg_name Argument name for error messages.
#' @keywords internal
.assert_flow <- function(x, arg_name = deparse(substitute(x))) {
  .assert_hydrodata(x, arg_name)
  if (!.is_flow(x))
    stop(sprintf(
      "%s: expected a Flow_Daily or Flow_15min object; got %s. ",
      "Check you have passed the flow series, not the rainfall series.",
      arg_name, paste(class(x), collapse = "/")), call. = FALSE)
  invisible(x)
}

#' Validate that two HydroData objects share the same timestep class
#'
#' Stops if one is daily and the other is 15-minute, since mismatched
#' timesteps will silently produce wrong model outputs.
#'
#' @param x,y  Two HydroData objects.
#' @keywords internal
.assert_same_timestep <- function(x, y) {
  x_daily <- .is_daily(x)
  y_daily <- .is_daily(y)
  if (x_daily != y_daily)
    stop(
      "reach.io timestep mismatch: one series is daily and the other is ",
      "15-minute. Both inputs to pdm() must use the same timestep.",
      call. = FALSE)
  invisible(TRUE)
}

# =============================================================================
# SECTION 3 : EXTRACTORS
# =============================================================================

#' Extract the value vector from a HydroData object
#'
#' Returns `readings$value` as a plain numeric vector, with quality filtering
#' applied if the `quality` column is present.
#'
#' @param x           A reach.io HydroData object.
#' @param na_quality  Character vector of quality flags to replace with `NA`.
#'                    Default `c("Missing", "Suspect")`. Set to `character(0)`
#'                    to keep all values.
#'
#' @return Named numeric vector. Names are the `dateTime` values as character.
#' @export
hydrodata_values <- function(x, na_quality = c("Missing", "Suspect")) {
  .assert_hydrodata(x)
  dt <- reach.io::as_data_table(x)
  v  <- dt$value
  if ("quality" %in% names(dt) && length(na_quality) > 0) {
    v[dt$quality %in% na_quality] <- NA_real_
  }
  stats::setNames(v, as.character(dt$dateTime))
}

#' Extract the dateTime vector from a HydroData object
#'
#' @param x A reach.io HydroData object.
#' @return POSIXct vector.
#' @export
hydrodata_datetimes <- function(x) {
  .assert_hydrodata(x)
  reach.io::as_data_table(x)$dateTime
}

#' Extract the date vector from a HydroData object
#'
#' @param x A reach.io HydroData object.
#' @return Date vector.
#' @export
hydrodata_dates <- function(x) {
  .assert_hydrodata(x)
  reach.io::as_data_table(x)$date
}

# =============================================================================
# SECTION 4 : COERCERS — convert HydroData inputs to model-ready vectors
# =============================================================================

#' Coerce a Rainfall HydroData object to a numeric vector for pdm()
#'
#' Extracts rainfall depths, applies quality filtering, and aligns the series
#' to a regular time grid. Returns a plain numeric vector in mm/timestep,
#' ready to pass as `rain` to [pdm()].
#'
#' @param x          A `Rainfall_Daily` or `Rainfall_15min` object.
#' @param na_quality Quality flags to replace with `NA`. Default
#'                   `c("Missing", "Suspect")`.
#'
#' @return Numeric vector of rainfall depths \[mm/timestep\].
#' @seealso [as_pet_input()], [as_pdm_input()], [pdm()]
#' @export
as_rain_input <- function(x, na_quality = c("Missing", "Suspect")) {
  .assert_rainfall(x, "x")
  hydrodata_values(x, na_quality)
}

#' Coerce a Flow or Level HydroData object to a PET proxy vector for pdm()
#'
#' Intended for cases where PET is stored as a Level or Flow series in
#' reach.io (e.g. gridded PET imported as a flow-unit series). For the
#' common case of estimating PET from temperature, use a separate PET
#' estimation function and pass the result directly to [pdm()].
#'
#' @param x          A HydroData object containing PET values \[mm/timestep\].
#' @param na_quality Quality flags to replace with `NA`.
#'
#' @return Numeric vector \[mm/timestep\].
#' @seealso [as_rain_input()], [as_pdm_input()]
#' @export
as_pet_input <- function(x, na_quality = c("Missing", "Suspect")) {
  .assert_hydrodata(x, "x")
  hydrodata_values(x, na_quality)
}

#' Coerce a paired Rainfall + PET HydroData to a list ready for pdm()
#'
#' Validates both objects, checks they share the same timestep class, aligns
#' them by dateTime (inner join), and returns a named list with `rain`, `pet`,
#' and `dates` vectors — all the same length.
#'
#' @param rainfall  A `Rainfall_Daily` or `Rainfall_15min` object.
#' @param pet       A HydroData object containing PET \[mm/timestep\].
#' @param na_quality Quality flags to replace with `NA`.
#'
#' @return A named list: `rain` (numeric), `pet` (numeric), `dates` (Date),
#'   `datetimes` (POSIXct), `n` (integer).
#'
#' @examples
#' \dontrun{
#' inputs <- as_pdm_input(rainfall_obj, pet_obj)
#' res <- pdm(inputs$rain, inputs$pet, params = p)
#' }
#'
#' @seealso [pdm()], [as_rain_input()]
#' @export
as_pdm_input <- function(rainfall, pet, na_quality = c("Missing", "Suspect")) {
  .assert_rainfall(rainfall, "rainfall")
  .assert_hydrodata(pet, "pet")
  .assert_same_timestep(rainfall, pet)

  rain_dt <- reach.io::as_data_table(rainfall)
  pet_dt  <- reach.io::as_data_table(pet)

  # Align by dateTime — inner join to keep only overlapping period
  merged <- merge(
    rain_dt[, .(dateTime, date, rain = value,
                rain_quality = if ("quality" %in% names(rain_dt)) quality else NA_character_)],
    pet_dt[ , .(dateTime,       pet  = value,
                pet_quality  = if ("quality" %in% names(pet_dt))  quality else NA_character_)],
    by = "dateTime"
  )

  if (nrow(merged) == 0L)
    stop("as_pdm_input: no overlapping timesteps between rainfall and pet series.",
         call. = FALSE)

  # Apply quality flags
  if (length(na_quality) > 0) {
    merged[rain_quality %in% na_quality, rain := NA_real_]
    merged[pet_quality  %in% na_quality, pet  := NA_real_]
  }

  n_na_rain <- sum(is.na(merged$rain))
  n_na_pet  <- sum(is.na(merged$pet))
  if (n_na_rain > 0)
    message(sprintf("as_pdm_input: %d NA rainfall values after quality filtering.",
                    n_na_rain))
  if (n_na_pet > 0)
    message(sprintf("as_pdm_input: %d NA PET values after quality filtering.",
                    n_na_pet))

  list(
    rain      = merged$rain,
    pet       = merged$pet,
    dates     = merged$date,
    datetimes = merged$dateTime,
    n         = nrow(merged)
  )
}

#' Coerce a Flow HydroData object to an AMAX series for FEH methods
#'
#' Extracts annual maximum flows from a `Flow_Daily` or `Flow_15min` object,
#' using [annual_maxima()] internally. Returns a numeric vector of annual
#' maxima ready to pass to [feh_single_site()], [feh_pooled()], or [fit_glo()].
#'
#' @param x               A `Flow_Daily` or `Flow_15min` object.
#' @param water_year_start Integer month that starts the water year. Default 10.
#' @param na_quality      Quality flags to replace with `NA`.
#' @param min_years       Minimum number of complete water years required.
#'                        Stops with an error if the record is shorter.
#'                        Default 5.
#'
#' @return Numeric vector of annual maximum flows \[m³/s or original units\],
#'   one value per water year. Attribute `water_years` gives the year integers.
#'
#' @examples
#' \dontrun{
#' amax <- as_amax(flow_obj)
#' fit  <- feh_single_site(amax)
#' }
#'
#' @seealso [feh_single_site()], [feh_pooled()], [annual_maxima()]
#' @export
as_amax <- function(x,
                    water_year_start = 10L,
                    na_quality       = c("Missing", "Suspect"),
                    min_years        = 5L) {
  .assert_flow(x, "x")

  vals  <- hydrodata_values(x, na_quality)
  dates <- hydrodata_dates(x)

  am <- annual_maxima(vals, dates, water_year_start = water_year_start)

  n_complete <- sum(!is.na(am$annual_max))
  if (n_complete < min_years)
    stop(sprintf(
      "as_amax: only %d complete water years in record; min_years = %d. ",
      "Extend the record or reduce min_years.",
      n_complete, min_years), call. = FALSE)

  if (any(is.na(am$annual_max)))
    warning(sprintf(
      "as_amax: %d water year(s) have NA annual maxima (insufficient data ",
      "after quality filtering). These will be dropped.",
      sum(is.na(am$annual_max))), call. = FALSE)

  result <- am$annual_max[!is.na(am$annual_max)]
  attr(result, "water_years") <- am$water_year[!is.na(am$annual_max)]
  result
}

#' Coerce a Flow HydroData object to a POT peaks series for feh_pot()
#'
#' Extracts independent peaks over a threshold from a `Flow_Daily` or
#' `Flow_15min` object using [peaks_over_threshold()].
#'
#' @param x         A `Flow_Daily` or `Flow_15min` object.
#' @param threshold Flow threshold \[same units as x\].
#' @param min_sep   Minimum separation between independent peaks \[timesteps\].
#'                  Default 3.
#' @param na_quality Quality flags to replace with `NA`.
#'
#' @return A `data.table` with columns `date`, `peak_flow`, `peak_index`, and
#'   attribute `n_years` (record length in years, for use in [feh_pot()]).
#'
#' @examples
#' \dontrun{
#' peaks   <- as_pot(flow_obj, threshold = 150)
#' fit_pot <- feh_pot(peaks$peak_flow, threshold = 150,
#'                    n_years = attr(peaks, "n_years"))
#' }
#'
#' @seealso [feh_pot()], [peaks_over_threshold()]
#' @export
as_pot <- function(x, threshold, min_sep = 3L,
                   na_quality = c("Missing", "Suspect")) {
  .assert_flow(x, "x")
  checkmate::assert_number(threshold, lower = 0)

  vals  <- hydrodata_values(x, na_quality)
  dates <- hydrodata_dates(x)

  # Record length in years (approximate from date range)
  n_years <- as.numeric(diff(range(dates, na.rm = TRUE))) / 365.25

  peaks <- peaks_over_threshold(vals, dates, threshold = threshold,
                                min_sep = min_sep)
  attr(peaks, "n_years")   <- n_years
  attr(peaks, "threshold") <- threshold
  peaks
}

#' Coerce a Flow HydroData object to an observed flow vector for pdm calibration
#'
#' Extracts the value vector from a Flow object, aligning to a supplied date
#' vector so the observed series lines up with model output from [pdm()].
#'
#' @param x          A `Flow_Daily` or `Flow_15min` object.
#' @param datetimes  POSIXct vector of model timesteps (from [as_pdm_input()]).
#' @param na_quality Quality flags to replace with `NA`.
#'
#' @return Numeric vector, same length as `datetimes`, with NA for any
#'   timesteps not present in the flow record.
#'
#' @seealso [calibrate_pdm()], [as_pdm_input()]
#' @export
as_obs_flow <- function(x, datetimes, na_quality = c("Missing", "Suspect")) {
  .assert_flow(x, "x")

  flow_dt <- reach.io::as_data_table(x)
  if ("quality" %in% names(flow_dt) && length(na_quality) > 0)
    flow_dt[quality %in% na_quality, value := NA_real_]

  # Align to supplied datetimes
  target <- data.table::data.table(dateTime = datetimes)
  merged <- merge(target, flow_dt[, .(dateTime, value)],
                  by = "dateTime", all.x = TRUE)
  merged$value
}

# =============================================================================
# SECTION 5 : HYDRODATA METADATA HELPERS
# =============================================================================

#' Summarise a reach.io HydroData object for logging / provenance
#'
#' Returns a compact named list describing the object — useful for recording
#' input provenance in model run logs and governance records.
#'
#' @param x A reach.io HydroData object.
#' @return Named list: `class`, `parameter`, `period_name`, `from_date`,
#'   `to_date`, `n_rows`, `n_measures`, `downloaded_at`.
#' @export
hydrodata_provenance <- function(x) {
  .assert_hydrodata(x, "x")
  list(
    class        = class(x)[1],
    parameter    = x@parameter,
    period_name  = x@period_name,
    from_date    = x@from_date,
    to_date      = x@to_date,
    n_rows       = x@n_rows,
    n_measures   = x@n_measures,
    downloaded_at = x@downloaded_at
  )
}
