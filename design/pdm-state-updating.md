# PDM State Updating Design

*Implemented in `R/pdm_state_update.R` | [UPCOMING_FEATURES §1.5](../UPCOMING_FEATURES.md)*
*Governance tier: 1 (operational)*
*Owner: Deputy Director (Technology) | Steward: Lead Developer (G7)*

---

## Problem

When `pdm()` is called from cold initial conditions (`S0 = Smax/2`, `Sg0 = 0`,
`Ss10 = 0`, `Ss20 = 0`), the model stores require a spin-up period — typically
one to two years — before they equilibrate and produce reliable flow estimates.
This is the burn-in (or warm-up) problem.

In operational flood forecasting there are two situations where this matters:

1. **Forecasting from arbitrary start dates** — a system that spins up on demand
   (e.g. triggered by a rain event) cannot afford a two-year warm-up.
2. **Historic/calibration runs** — the warmup period must be excluded from
   performance metrics (currently via the `warmup` argument to `gof_metrics()`),
   wasting data and obscuring performance near the start of the record.

State updating solves both problems by using recent observed data to drive the
model stores to a realistic state before the forecast begins.

---

## Model States

The PDM has four state variables at each timestep:

| State | Symbol | Description |
|---|---|---|
| Soil moisture | `S` [mm] | Controlled by rainfall, AET, and recharge — driven primarily by the recent rainfall/PET history |
| Groundwater store | `Sg` [mm] | Slow-response store; the primary control on baseflow |
| Surface reservoir 1 | `Ss1` [mm] | Fast-response routing — first linear reservoir |
| Surface reservoir 2 | `Ss2` [mm] | Fast-response routing — second linear reservoir |

`S` equilibrates relatively quickly (days to weeks) because it is directly forced
by rainfall and PET. `Sg` equilibrates slowly (months to years) because recharge
is small relative to storage.

---

## Methods

### Method 1: `"window"` (assimilation window run)

Run `pdm()` over a recent period of observed rain and PET. The states at the end
of the window become the initial conditions for the forecast.

```
recent observations
──────────────────►  pdm()  ──►  S(T), Sg(T), Ss1(T), Ss2(T)  ──►  forecast pdm()
[window: last N ts]              [initial conditions at time T]
```

**Window length guidelines** (daily timesteps):

| Window | Adequacy |
|---|---|
| < 14 days | Insufficient — groundwater store not equilibrated |
| 30 days | Practical minimum for operational use |
| 90 days | Good for most catchments |
| 365 days | Recommended for confidence in Sg; equivalent to a one-year burn-in |

For sub-daily timesteps, multiply thresholds proportionally (e.g. 15-min:
× 96 per day, so 30-day equivalent = 2,880 timesteps).

**When to use**: simple operational contexts where a long window (≥ 90 days) of
recent rain/PET observations is available and no observed flow is needed.

---

### Method 2: `"inversion"` (algebraic state adjustment)

As `"window"`, then adjusts `Sg`, `Ss1`, and `Ss2` at the forecast origin so
that the implied instantaneous outflow of the stores matches the observed flow
`Q_obs(T)`. This corrects any remaining drift in the routing stores after the
window run.

**Step-by-step:**

1. Run `pdm()` over the window → produces `Sg(T)`, `Ss1(T)`, `Ss2(T)`,
   `Qb(T)`, `Qf(T)`, `Q(T)`.

2. Compute the total store flow (excluding the constant term `qc`):

   ```
   Q_stores = Q(T) - qc
   Q_obs_stores = Q_obs(T) - qc
   scale = Q_obs_stores / Q_stores
   ```

3. Adjust the groundwater store, preserving the BFI partition:

   ```
   Qb_adj = Qb(T) × scale
   Sg_adj = kb × Qb_adj^(1/m)      [inversion of the reservoir equation]
   ```

4. Adjust the surface routing stores uniformly:

   ```
   Ss1_adj = Ss1(T) × scale
   Ss2_adj = Ss2(T) × scale
   ```

5. Soil moisture `S` is **not adjusted** — see note below.

**Diagnostic check** — the implied instantaneous outflow of the adjusted stores
(assuming zero new inflow):

```
Qb_implied = (Sg_adj / kb)^m
Qf_implied = Ss2_adj / k2
Q_implied  = Qb_implied + Qf_implied + qc
```

For the linear case (m = 1), `Q_implied ≈ Q_obs(T)` exactly (up to the BFI
decomposition). For non-linear groundwater (m ≠ 1), there is a small residual
because the Euler step used for non-linear routing introduces a bias; the match
is still close.

**When to use**: operational forecasting where a short window (30 days) is
available alongside observed flow at the forecast origin. The inversion corrects
for window-length inadequacy in the Sg estimate.

---

## Why Soil Moisture is Not Adjusted

Inverting `S` from `Q` is not tractable in closed form. The chain is:

```
S  →  c*(S)  →  Qr(P, S)  →  recharge(S)  →  Sg  →  Qb
```

This depends on the current rainfall `P` at time T, the distribution-specific
runoff function, and the recharge formulation — a multi-step non-linear map with
no clean inverse.

More importantly, `S` is primarily determined by the accumulated rainfall/PET
history. Given a reasonable window length (≥ 30 days), the window run will drive
`S` close to its true value automatically. Adjusting `Sg` and the surface stores
(which have the longer memory) gives most of the benefit.

---

## Usage

```r
library(reach.hydro)

# Calibrated params from a prior calibrate_pdm() run
cal <- calibrate_pdm(rain_hist, pet_hist, obs_q_hist, dist = "pareto")

# Assimilation window: last 30 days of recent data
n_window <- 30L
upd <- update_pdm_states(
  rain  = tail(rain_recent, n_window),
  pet   = tail(pet_recent,  n_window),
  obs_q = tail(obs_q_recent, n_window),
  params  = cal$params,
  method  = "inversion"
)

print(upd)
#> <PdmStateUpdate>
#>   Method          : inversion
#>   Window length   : 30 timesteps
#>   --- States at forecast origin ---
#>   S0  (soil)      : 187.43 mm
#>   Sg0 (gw store)  : 142.61 mm  [was 118.34 mm]
#>   Ss10 / Ss20     : 8.24 / 6.71 mm
#>   --- Flow at forecast origin ---
#>   Q obs           : 0.0412 mm/ts
#>   Q sim (window)  : 0.0341 mm/ts
#>   Q implied (upd) : 0.0408 mm/ts

# Use updated states as initial conditions for the forecast
fcast <- pdm(rain_fcast, pet_fcast, params = cal$params,
             S0   = upd$S0,
             Sg0  = upd$Sg0,
             Ss10 = upd$Ss10,
             Ss20 = upd$Ss20)
```

---

## Return Value

`update_pdm_states()` returns a `PdmStateUpdate` object:

| Field | Type | Description |
|---|---|---|
| `S0` | numeric | Updated soil moisture — pass as `S0` to `pdm()` |
| `Sg0` | numeric | Updated groundwater store |
| `Ss10` | numeric | Updated surface reservoir 1 |
| `Ss20` | numeric | Updated surface reservoir 2 |
| `method` | character | `"window"` or `"inversion"` |
| `hindcast` | `ReachHydroResult` | Full window run, before any inversion |
| `diagnostics` | list | Q_obs, Q_sim_window, Q_implied_update, Sg_before/after, Ss1/Ss2 before/after, scale_factor |

The `hindcast` object is useful for inspecting how well the model performed
during the assimilation window before the forecast starts.

---

## Integration with the Forecast Workflow

```
                    calibrate_pdm()
                          │
                          ▼
                      PdmParams  ─────────────────────────────────┐
                          │                                        │
                          ▼                                        │
              update_pdm_states()                                  │
            (assimilation window)                                  │
                          │                                        │
                          ▼                                        │
                   PdmStateUpdate                                  │
              (S0, Sg0, Ss10, Ss20)                                │
                          │                                        │
                          └──────────────────────────────►  pdm()
                                                           (forecast)
                                                                   │
                                                                   ▼
                                                          ReachHydroResult
```

---

## Planned Extensions

- **`"ensemble"` method** (future, §1.4): run an ensemble of window simulations
  with perturbed parameters, select the member whose end-states produce the best
  match to `Q_obs(T)`. Supports uncertainty-aware forecasting.

- **Multi-step nudging**: iterate the inversion over the last few timesteps of
  the window (rather than just the final point) for robustness to noisy
  observations.

- **reach.validate integration** (§4): feed `upd$hindcast` directly into
  `validate()` to report window fit quality as a pre-forecast health check.

- **write_pdm_run() integration** (§5.5): store the `PdmStateUpdate` alongside
  the forecast run in the Parquet/manifest convention, enabling retrospective
  audit of what states were used at each forecast issuance.

---

*Last updated: 2026-03-21*
