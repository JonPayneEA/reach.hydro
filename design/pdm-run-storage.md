# PDM Run Storage Design

*Relates to [UPCOMING_FEATURES §5.5](../UPCOMING_FEATURES.md#55-pdm-run-storage--parquet--manifest-convention)*
*Governance tier: 2 (analytical) for calibration runs; 1 (operational) for historic and forecast runs*
*Owner: Deputy Director (Technology) | Steward: Lead Developer (G7)*

---

## Problem

PDM can be run continuously over long periods. A 30-year hindcast at 15-minute
timesteps produces ~1.05 million rows × 16 output columns per simulation. Calibration
jobs run thousands of such simulations internally. Without a standard storage
convention, runs are not reproducible, provenance is lost, and results cannot be
compared across catchments or time.

---

## Format: Apache Parquet

Apache Parquet (via the `arrow` package) is the recommended format for all
`ReachHydroResult` outputs.

**Rationale:**

| Property | Why it matters |
|---|---|
| Columnar storage | Reading only `Q` from a 30-year run reads ~1/16 of the data |
| Open standard | No R dependency to read it back; data outlives any tool (OAIS-aligned) |
| Compression | Temporal autocorrelation in time series gives ~10× reduction; a 30-year 15-min run lands at ~10–15 MB |
| Ecosystem | Readable by DuckDB, Python/pandas, and Power BI without conversion |

```r
# Write
arrow::write_parquet(as.data.frame(result), file.path(run_dir, "result.parquet"))

# Read back as data.table
arrow::read_parquet(path) |> data.table::as.data.table()
```

For runs longer than ~5 years, partition by water year to enable cheap
time-range queries:

```r
result_dt <- as.data.frame(result) |>
  dplyr::mutate(water_year = format(date, "%Y"))  # or derive from dateTime

arrow::write_dataset(result_dt, run_dir,
                     partitioning = "water_year",
                     format       = "parquet")
```

---

## Directory Structure

```
reach_runs/
  {catchment_id}/
    {run_type}/            # "historic" | "calibration" | "operational"
      {run_id}/
        result.parquet     # ReachHydroResult time series (or partitioned by water_year/)
        manifest.json      # provenance, params, run metadata
```

`run_id` is a deterministic SHA-256 hash of the parameter set and an input
fingerprint (e.g. `digest::digest(list(params, nrow(rain), range(dates)))`).
Re-running identical inputs produces the same `run_id`, which makes deduplication
and cache-checking straightforward.

---

## Manifest (Governance Record)

The `manifest.json` file is the primary governance artefact. It captures
everything needed to reproduce or audit a run without loading the result.

```json
{
  "run_id":          "a3f7c2d1...",
  "run_type":        "calibration",
  "catchment_id":    "25001",
  "tier":            1,
  "created_at":      "2026-03-21T09:14:00Z",
  "created_by":      "jpayne",
  "package_version": "0.2.0",
  "params": {
    "dist": "pareto", "cmax": 342.1, "b": 0.38,
    "St": 18.4, "kg": 163.0, "k1": 6.2, "k2": 5.0,
    "be": 5.0, "kb": 200.0, "m": 1.0
  },
  "Smax":          187.3,
  "dist":          "pareto",
  "period_start":  "1990-10-01",
  "period_end":    "2020-09-30",
  "n_timesteps":   1051920,
  "dt_minutes":    15,
  "input_provenance": {
    "class":        "Rainfall_15min",
    "parameter":    "rainfall",
    "period_name":  "Moorhouse 1990–2020",
    "downloaded_at":"2026-03-01T11:00:00Z"
  },
  "gof": {
    "nse": 0.84, "kge": 0.81, "pbias": -1.2
  }
}
```

Fields and their sources:

| Field | Source |
|---|---|
| `params` | `attr(result, "params")` — already stored on `ReachHydroResult` |
| `Smax`, `dist` | `attr(result, "Smax")`, `attr(result, "dist")` |
| `input_provenance` | `hydrodata_provenance()` — already implemented in `reach_io_compat.R` |
| `gof` | `gof_metrics(obs_q, result$Q, warmup)` — optional, added when `obs_q` supplied |
| `package_version` | `as.character(utils::packageVersion("reach.hydro"))` |

The `gof` block means you can query performance across all runs without loading
any result files.

---

## Calibration Runs

Nelder-Mead at 3,000 iterations runs 3,000 full-length simulations internally.
**Do not store intermediate trial runs.** Store only:

1. The **final calibrated `PdmParams`** in the manifest `params` block
2. The **final simulation** as `result.parquet`
3. Optionally, a `calibration_summary.json` with the raw `optim()` convergence
   message and final objective value

For the planned §1.2 multi-objective calibration, store the Pareto front as a
separate flat file (one row per non-dominated parameter set, no time series):

```
reach_runs/{catchment_id}/calibration/{run_id}/
  result.parquet           # simulation at the "best" Pareto solution
  pareto_front.parquet     # all non-dominated params + objective values
  manifest.json
```

---

## Proposed Package Helpers

Add `write_pdm_run()` and `read_pdm_run()` to a new file `R/pdm_io.R` (Tier 2).
These thin wrappers enforce the convention and auto-populate the manifest from
attributes already present on `ReachHydroResult`:

```r
write_pdm_run(result,
              path,
              catchment_id,
              run_type     = c("historic", "calibration", "operational"),
              obs_q        = NULL,
              created_by   = Sys.info()[["user"]],
              dt_minutes   = 15L,
              partition_by_year = FALSE)

read_pdm_run(path)
# Returns list(result = ReachHydroResult, manifest = list)
```

`write_pdm_run()` will:

- Derive `run_id` deterministically via `digest::digest()`
- Extract `params`, `Smax`, `dist` from `attr(result, ...)` (no user input needed)
- Call `hydrodata_provenance()` on the original reach.io inputs if supplied
- Compute `gof` if `obs_q` is provided
- Write `result.parquet` (optionally partitioned) + `manifest.json`
- Return the `run_id` invisibly

`read_pdm_run()` will:

- Read `manifest.json` and `result.parquet`
- Reconstruct the `ReachHydroResult` S3 class with `params`, `Smax`, and `dist`
  attributes re-attached

---

## Alignment with Planned Features

| Upcoming feature | Storage interaction |
|---|---|
| §1.1 Rcpp loop | No change to storage format — API is identical |
| §1.2 Multi-objective calibration | Add `pareto_front.parquet` alongside `result.parquet` |
| §1.4 Monte Carlo uncertainty | Ensemble stored as wide parquet (one column per member) or partitioned by `member_id` |
| §5.2 S7 class migration | `read_pdm_run()` reconstructs S7 class instead of S3 — manifest format unchanged |
| §5.3 reach.io HydroData round-trip | `write_pdm_run()` will optionally coerce to `Flow_Daily`/`Flow_15min` for storage via the reach.io data layer |

---

## What Not to Use

| Format | Reason |
|---|---|
| RDS | R-specific; no external tooling; no column-level reads; fails open-standard requirement |
| CSV | No compression; no schema; loses float64 precision; slow for 1M+ rows |
| One flat file per catchment | Defeats partitioning; makes retention management difficult |
| Embedded database (DuckDB in-process) | Adds operational dependency; Parquet + DuckDB is a better split if SQL queries are needed later |

---

*Last updated: 2026-03-21*
