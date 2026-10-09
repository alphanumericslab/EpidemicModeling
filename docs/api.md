# API contracts

Author: **Reza Sameni · Emory University**.

## Shared numerical conventions

| Quantity | Both languages |
| --- | --- |
| Population compartments | Fractions; counts are divided by population before fitting |
| Rates | Inverse time; examples use days |
| States and observations | Variables × samples; Python single series may be 1-D |
| Covariances and gain | Variables × variables/observations × samples |
| SEIRP grid | K=round(duration/dt), initial included, last time (K−1)dt |
| Controlled SI grid | K samples, initial included |
| SI-alpha grid | K transitions, initial excluded |
| JHU first/threshold index | One-based in both languages; zero if absent |
| Regression start index | Python zero-based; MATLAB one-based |
| Noise | Standard-normal 3 × K array supplied to both for deterministic parity |
| Unconstrained boundaries | NaN entries in final state/covariance |

Floating-point linear algebra is compared with tolerances, not bitwise identity.
Do not interpret the legacy growth factor as a calibrated biological reproduction
number without an explicit generation-interval model.

## Compartment parameter glossary

| Parameter | Meaning |
| --- | --- |
| alpha_e, alpha_i | Transmission from exposed and infected individuals |
| kappa | Exposed-to-infected transition rate |
| rho | Exposed-to-recovered transition rate |
| beta | Infected removal/recovery rate |
| mu | Infected-to-passed mortality rate |
| gamma | Recovered-to-susceptible rate in SEIRP; contact-response speed in SI-alpha |
| beta_0, beta_s | Recovery rates before and after healthcare saturation |
| mu_0, mu_s | Mortality rates before and after saturation |
| sigma, i_0 | Saturation transition width and infected-fraction threshold |
| s0,e0,i0,r0,p0 | Initial population fractions |
| duration / MATLAB T | Requested simulation duration |
| count / MATLAB K | Number of samples or transitions, according to model |
| dt | Forward-Euler time step |
| u, u_max | Intervention schedule and maximum intervention levels |
| a, b | NPI influence vector and baseline contact-rate bias |
| alpha_min, alpha_max | Contact-rate bounds |
| s/i/alpha_noise_std | Process-noise standard deviations |

## SI-alpha EKF parameter struct / dictionary

`dt`, `beta`, `gamma`, `a`, `b`, `u_min`, `u_max`, `alpha_min`, `alpha_max`,
`obs_type` are required. `obs_type` is `NEWCASES` or `TOTALCASES`. Three-state
forward bounds use `s_min` and `i_min` (normally zero). Six-state control also
requires `epsilon`, `w`, `sigma`. `epsilon` mixes case and intervention objectives,
`w` is the intervention cost vector, and `sigma` is the slope used by the
historical approximate switching Jacobian. The SI Hessian callbacks are zero
even when `order=2`; the exponential filter has analytic nonzero Hessian terms.
Process-noise means are accepted but intentionally unused by the SI transition,
matching the source model; the exponential model does use its noise mean.

## Generic callback protocol

All callbacks are snake case. Their shapes follow the shared conventions above.
`k` is zero-based in Python and one-based in MATLAB. `params` is passed unchanged.

```text
state_hard_margins(state, params, k) -> bounded_state
obs_hard_margins(observation, params, k) -> bounded_observation
nlin_state_update(u, state, w_bar, params, k) -> (u_opt, next_state)
nlin_obs_update(u, state, v_bar, params, k) -> observation
state_jacobians(u, state, w_bar, params, k) -> (A, B)
obs_jacobian(u, state, v_bar, params, k) -> (C, D)
state_hessian_terms(u, state, P, w_bar, Q, params, k) -> (fs, Fsp, fw, Fwp)
obs_hessian_terms(u, state, P, v_bar, R, params, k) -> (gs, Gsp, gv, Gvp)
```

The generic result order is `u_opt, u_opt_smooth, s_minus, s_plus, s_smooth,
p_minus, p_plus, p_smooth, k_gain, innovations, rho`. Python returns a named tuple
that also supports positional unpacking. The older dedicated filter omits
`u_opt_smooth`. The Coder variant orders `p_minus, p_plus, k_gain` before
`s_smooth, p_smooth`; use `epidemic_modeling.codegen` for that alternate order.
The final generic smoothed control remains zero by historical convention;
`optimal_npi` explicitly fills its final policy from the final forward control.

## Growth parameter glossary

`new_cases` / MATLAB `NewCases`: regularly sampled incident counts. `wlen`:
integer window length; `generation_period`: lag in samples; `time_unit`:
positive sampling-scale parameter; `causal`: true for trailing windows, false
for symmetric windows of length `2*floor(wlen/2)+1`.

Log fitting uses slope per sample: `rt=exp(slope)`, `growth=slope/time_unit`.
Nonlinear fitting preserves the source's `n/time_unit` independent variable and
then returns `growth=fitted_rate/time_unit`, `rt=exp(fitted_rate)`. Therefore for
non-unit `time_unit` the two reported growth definitions differ. This is a
preserved historical convention rather than a newly standardized estimator.
`exp_fit=amplitude*rt` is a one-step forecast. Log-fitting unfitted endpoints
are amplitude=rt=fit=1 and growth=0. Nonlinear causal endpoints use lagged raw
counts; centered endpoints use raw counts. Log fit requires positive counts;
ratio zeros yield NaN/Inf; nonlinear zero windows use the current count and no
growth. Time-series zeros must not be silently replaced without documentation.

## Function reference

Python signatures below are executable API contracts. MATLAB uses positional
arguments and column vectors for intervention parameters. See MATLAB `help`
for its signatures; `prescribe_npi` requires a sixth positional `model_file`.
For the upgraded pipeline, model bundles use shared JSON schema version 1.

### models

#### `seirp`

```python
seirp(alpha_e, alpha_i, kappa, rho, beta, mu, gamma,
          s0, e0, i0, r0, p0, duration, dt)
```

Solve the SEIRP model using forward Euler.

Parameters
----------
alpha_e, alpha_i, kappa, rho, beta, mu, gamma : float or array_like
    Nonnegative rates per time unit, scalar or at least K-1 samples.
s0, e0, i0, r0, p0 : float
    Initial population fractions (susceptible, exposed, infected,
    recovered, and passed/deceased).
duration, dt : float
    Requested duration and step size in the same time unit.

Returns
-------
s, e, i, r, p : ndarray, shape (K,)
    K = round(duration/dt) samples, including the initial conditions.
    The final time is (K-1)*dt, preserving the original MATLAB convention.

Notes
-----
No clipping is applied: reduce dt if Euler produces negative fractions.
The sum of the compartments is conserved up to floating-point precision.

#### `seirp_saturated_resource`

```python
seirp_saturated_resource(alpha_e, alpha_i, kappa, rho, gamma,
                             s0, e0, i0, r0, p0, duration, dt,
                             beta_0, beta_s, mu_0, mu_s, sigma, i_0)
```

Solve SEIRP with a smooth transition to saturated healthcare rates.

``h = (tanh((i-i_0)/sigma)+1)/2`` interpolates recovery from beta_0 to
beta_s and mortality from mu_0 to mu_s. sigma must be positive; i_0 is
the prevalence threshold. Other inputs and outputs follow :func:`seirp`.

#### `si_controlled`

```python
si_controlled(alpha, beta, s0, i0, count, dt)
```

Integrate bounded SI fractions with a supplied contact-rate schedule.

alpha is scalar or length >= count-1; beta is the removal rate. Returns
two length-count arrays including the initial state. Each Euler update
is clipped to [0, 1], as in the original MATLAB solver.

#### `si_alpha_controlled`

```python
si_alpha_controlled(u, s0, i0, alpha0, u_max, alpha_min, alpha_max,
                        gamma, a, b, beta, s_noise_std, i_noise_std,
                        alpha_noise_std, count, dt, *, rng=None, noise=None)
```

Integrate NPI-driven SI-alpha dynamics, returning post-update samples.

u has shape (interventions, count); a and u_max have one element per
intervention. gamma is the contact-rate response speed. The equilibrium
contact rate is b + a @ (u_max-u). States and contact rate are clipped.
The returned three length-count arrays exclude the initial condition.

Pass ``noise`` of shape (3, count) containing standard-normal draws for
exact MATLAB/Python comparison. Otherwise ``rng`` is a NumPy Generator.
Zero noise standard deviations yield deterministic trajectories.

### growth

#### `exp_model`

```python
exp_model(params, t)
```

Evaluate ``amplitude * exp(growth * t)`` for a two-element parameter vector.

#### `rt_exp_fit_gen_ratios`

```python
rt_exp_fit_gen_ratios(new_cases, wlen, generation_period, time_unit)
```

Estimate lagged growth and its causal zero-padded moving average.

Returns (rt, growth, rt_smoothed, growth_smoothed), each length N.
growth[:generation_period] = 0; subsequent growth is log case ratio
divided by generation_period. rt = exp(growth*time_unit). Zero counts
propagate IEEE NaN/Inf, preserving the original definition.

#### `rt_exp_fit_log_lin_reg`

```python
rt_exp_fit_log_lin_reg(new_cases, wlen, time_unit, causal=True)
```

Fit rolling log cases using closed-form linear least squares.

Returns (rt, amplitude, growth, exp_fit). Causal windows end at each
sample; centered windows have 2*floor(wlen/2)+1 points. Unfitted ends
have rt=amplitude=exp_fit=1 and growth=0, matching MATLAB. rt=exp(slope),
growth=slope/time_unit, and exp_fit=amplitude*rt (one-step forecast).
Positive cases are required; zeros produce undefined log regressions.

#### `rt_exp_fit_nonlin_ls`

```python
rt_exp_fit_nonlin_ls(new_cases, wlen, time_unit, causal=True)
```

Fit rolling exponential cases by nonlinear least squares.

Inputs and return order follow :func:`rt_exp_fit_log_lin_reg`. The fit
uses n/time_unit as its independent variable, preserving the original
unusual scaling: growth = fitted_rate/time_unit and rt=exp(fitted_rate).
Causal unfitted amplitudes are delayed raw cases; centered ends retain
raw cases. Windows with zeros use the current count and zero growth.

### kalman

#### `generic_extended_kalman_filter`

```python
generic_extended_kalman_filter(u, x, handles, params, s_init, ps_init,
        s_final, ps_final, w_bar, v_bar, q_w, r_v, beta=1., gamma=1.,
        inv_monitor_len=21, order=1, *, covariance_update='joseph',
        switching_rho_epsilon=True)
```

Run an EKF and fixed-interval Rauch--Tung--Striebel smoother.

Parameters
----------
u, x : array_like, shape (controls, N), (observations, N)
    Inputs and observations. A column with any NaN observation skips
    correction. Model callbacks may replace NaN controls with optima.
handles : mapping or object
    Eight snake-case callbacks: state_hard_margins, obs_hard_margins,
    nlin_state_update, nlin_obs_update, state_jacobians, obs_jacobian,
    state_hessian_terms, obs_hessian_terms. Signatures match MATLAB,
    except time index k is zero-based. See docs/api.md.
params : object
    Passed unchanged to each callback.
s_init, ps_init : array_like
    Initial mean (m,) and covariance (m,m).
s_final, ps_final : array_like
    Final smoothing boundary values; NaN entries remain unconstrained.
w_bar, v_bar : array_like
    Process and observation noise means.
q_w, r_v : array_like
    Fixed covariance matrices, scalar variance series, or covariance
    stacks with shape (dimension, dimension, N).
beta, gamma : float
    Observation-noise adaptation factor [0,1] and covariance stability
    factor (0,1], respectively. beta=gamma=1 disables these adjustments.
inv_monitor_len : int
    Positive length of the rolling innovation monitor.
order : {1, 2}
    Linearized EKF or callback-supplied second-order corrections.

Returns
-------
filter_result
    Tuple of eleven arrays in the original generic MATLAB output order.
    Covariances and gains retain a final time axis; rho is squeezed.
    The final smoothed control is zero, preserving the legacy convention.

Notes
-----
Fixed R may adapt; time-dependent R does not. Joseph covariance updates
and symmetrization follow the generic original. The optional simple
update is used only by the two older dedicated filter ports.

#### `si_alpha_model_ekf`

```python
si_alpha_model_ekf(u,x,params,s_init,ps_init,s_final,ps_final,w_bar,v_bar,q_w,r_v,beta=1,gamma=1,inv_monitor_len=21,order=1)
```

Estimate three SI-alpha states; inputs/outputs follow the generic EKF contract.

#### `si_alpha_model_ekf_opt_controlled`

```python
si_alpha_model_ekf_opt_controlled(u,x,params,s_init,ps_init,s_final,ps_final,w_bar,v_bar,q_w,r_v,beta=1,gamma=1,inv_monitor_len=21,order=1)
```

Estimate six SI-alpha/costate states and fill NaN controls using phi > 0.

#### `si_alpha_model_backward_ekf`

```python
si_alpha_model_backward_ekf(u,x,params,s_init,ps_init,s_final,ps_final,w_bar,v_bar,q_w,r_v,beta=1,gamma=1,inv_monitor_len=21,order=1)
```

Run the three-state reverse-time SI filter and return chronological arrays.

#### `si_alpha_model_backward_ekf_opt_controlled`

```python
si_alpha_model_backward_ekf_opt_controlled(u,x,params,s_init,ps_init,s_final,ps_final,w_bar,v_bar,q_w,r_v,beta=1,gamma=1,inv_monitor_len=21,order=1)
```

Run six-state reverse-time filtering with bounded optimal NPI controls.

#### `new_case_ekf_estimator_with_optimal_npi`

```python
new_case_ekf_estimator_with_optimal_npi(u,x,params,s_init,ps_init,s_final,ps_final,w_bar,v_bar,q_w,r_v,beta=1,gamma=1,inv_monitor_len=21,order=1)
```

Port the older six-state estimator (simple covariance; switching phi >= 0).

Returns ten arrays: u_opt, s_minus, s_plus, s_smooth, p_minus, p_plus,
p_smooth, k_gain, innovations, rho, as in Tools' original function.

#### `rt_exp_fit_ekf`

```python
rt_exp_fit_ekf(x,s_init,params,w_bar,v_bar,ps_init,q_w,r_v,beta=1,gamma=1,inv_monitor_len=21,order=1)
```

Estimate cases and bounded growth with the original two-state EKF/EKS.

params = [time_scale, growth_memory, growth_saturation]; saturation > 0.
Returns nine arrays in the MATLAB order: s_minus, s_plus, p_minus,
p_plus, k_gain, s_smooth, p_smooth, innovations, rho. Missing observations
skip correction. Both first- and second-order filters are supported.

### data

#### `read_covid19_data`

```python
read_covid19_data(confirmed_datafile, death_datafile, recovered_datafile, region_list, min_cases)
```

Aggregate Johns Hopkins wide time-series tables by country substring.

Inputs are CSV paths with four metadata columns followed by aligned
date columns. Returns (total, infected, recovered, deceased, first_index,
threshold_index, num_days). Case arrays have shape (regions, days).
Indices are one-based for MATLAB parity, with 0 when never reached.

#### `read_oxford_data`

```python
read_oxford_data(path)
```

Read an Oxford/XPRIZE CSV, normalize region blanks, and sort daily rows.

Date accepts YYYYMMDD or ISO strings. CountryName, Date and ConfirmedCases
are required. Duplicate country/region/date observations are rejected
because they make daily fitting ambiguous. No online download is used.

#### `prepare_cases`

```python
prepare_cases(cumulative, window=7)
```

Return cleaned daily cases and their causal zero-padded moving average.

Daily cases begin with zero, negative revisions are clipped to zero,
and nonfinite differences become zero. This preprocessing uses only
present and past observations and is shared by both updated languages.

### npi

#### `npi_cost`

```python
npi_cost(new_cases, inputs, weights)
```

Return mean daily cases and mean elementwise weighted NPI intensity.

inputs is interventions-by-time; weights is a scalar, intervention
vector, or matching matrix. The denominator for NPI cost includes both
interventions and days, exactly as in the original MATLAB function.

#### `nonnegative_least_squares`

```python
nonnegative_least_squares(x,y,max_iter=10000,tolerance=1e-10)
```

Solve nonnegative least squares by matched cyclic coordinate descent.

Both languages use the same zero initialization, coordinate order,
relative stopping criterion, and iteration cap. This avoids differences
between MATLAB's lsqnonneg and SciPy solvers in the updated pipeline.

#### `default_si_params`

```python
default_si_params(u_max, *, epsilon=.5, weights=None)
```

Create reproducible SI-alpha parameters used by both implementations.

beta=-log(.01)/21, gamma=1/7, dt=1 day; other choices are documented
educational assumptions from the historical implementation.

#### `fit_npi_model`

```python
fit_npi_model(cumulative, inputs, population, u_max, *, regression_start=0)
```

Fit a two-pass EKF/NNLS SI-alpha model to one geographic time series.

cumulative is length N; inputs is (interventions,N); population is a
positive count. The first pass estimates alpha with zero NPI response;
NNLS regresses alpha on u_max-u, and the second pass refines alpha with
that response. A second NNLS uses the refined series (fixing the legacy
first-series reuse). Returns a serializable dict with parameters, final
state/covariance, coefficients and training length. See docs/migration.md.

#### `forecast_npi`

```python
forecast_npi(model, inputs)
```

Forecast daily cases and post-update states under a supplied NPI schedule.

Returns (cases, states), where states is (3,N). Cases use the pre-update
incidence N_population*s*i*alpha, avoiding a one-day shift between
forecasted cases and the transition producing those cases.

#### `optimal_npi`

```python
optimal_npi(model, horizon, weights, epsilon=.5, iterations=8)
```

Solve a bounded finite-horizon SI-alpha control problem by EKF/EKS shooting.

Six states include three costates. Repeated filtering/smoothing imposes
zero terminal costates; NaN controls select lower/upper bounds using
the paper's switching function. Returns (controls, cases, states).
This is a finite-iteration solver, not a claim of globally optimal convergence.

#### `train_npi_prescriptor`

```python
train_npi_prescriptor(start_date_str,end_date_str,data_file,geo_file,populations_file,included_ip,npi_maxes,trained_model_params_file, *, regression_start_date=None)
```

Train selected geographies from local Oxford data and save a JSON model bundle.

CountryName/RegionName identify rows in all three input CSVs.
Population2020 is required in the population table. Dates are inclusive
ISO strings. Missing interventions forward-fill then become zero.
Returns the same bundle written to trained_model_params_file.

#### `prescribe_npi`

```python
prescribe_npi(start_date_str,end_date_str,ip_file,costs_file,output_file, *, model_file,epsilon=.5)
```

Write XPRIZE-format daily prescriptions from a trained JSON bundle.

ip_file supplies selected CountryName/RegionName pairs and historic
interventions. costs_file has intervention weights in the bundle's
column names. model_file is explicit, replacing historical hardcoded
local paths. Returns the written DataFrame. No external service is used.

#### `train_predict_prescribe_npi`

```python
train_predict_prescribe_npi(npi_weights,human_npi_cost_factor,start_train_date_str,end_train_date_str,
        start_regression_date_str,end_predict_prescribe_date_str,data_file,geo_file,populations_file,
        included_ip,npi_mins,npi_maxes,trained_model_params_file)
```

Train and return per-geography fixed-policy and optimal-policy forecasts.

Dates and CSV arguments follow train_npi_prescriptor. npi_mins/maxes are
explicit bounds; output is a list of dicts with controls, cases, states,
and (human,NPI) costs. All numerical work is in src, outside notebooks.

#### `forecast_quality_assessment`

```python
forecast_quality_assessment(npi_weights,human_npi_cost_factor,start_train_date_str,end_train_date_str,
        start_regression_date_str,end_predict_prescribe_date_str,max_look_ahead_days,data_file,geo_file,
        populations_file,included_ip,npi_mins,npi_maxes,trained_model_params_file)
```

Evaluate held-out fixed-policy forecasts without fitting future observations.

Returns a table with country/region, number of held-out days, MAE, RMSE,
and bias. Actual held-out NPI schedules are known scenario inputs; case
data after the training endpoint are used only for scoring.

### layers

#### `exp_layer`

```python
exp_layer(x, alpha)
```

Evaluate the learnable exponential activation ``exp(alpha*x)``.

x and alpha must be broadcast-compatible; supply MATLAB Alpha weights
explicitly. The optional torch_exp_layer factory adds automatic gradients.

#### `my_tanh_layer`

```python
my_tanh_layer(x, alpha)
```

Evaluate ``alpha*tanh(x/alpha)`` with broadcast-compatible nonzero weights.

Pass the original MATLAB Alpha array to reproduce an existing layer.
A zero scale is rejected because the historical expression is undefined.

#### `torch_exp_layer`

```python
torch_exp_layer(alpha)
```

Create a trainable PyTorch module initialized with explicit exponential weights.

Requires the optional ``neural`` dependency. Input tensor and alpha must
be broadcast-compatible. No GPU, network, or random initialization is used.

#### `torch_my_tanh_layer`

```python
torch_my_tanh_layer(alpha)
```

Create a trainable PyTorch scaled-tanh module with explicit nonzero weights.

### spatial

#### `diffusion_2d`

```python
diffusion_2d(initial, diffusion, dt, spacing, steps)
```

Solve 2-D diffusion by explicit five-point finite differences.

initial is a finite (rows,columns) concentration grid. Periodic boundaries
conserve total mass. Returns (rows,columns,steps+1), including the initial
grid. The stability condition diffusion*dt/spacing**2 <= 1/4 is enforced.

#### `population_motion_2d`

```python
population_motion_2d(positions, velocities, dt, steps, box_size=1.)
```

Move agents with reflecting square boundaries using a triangular-wave map.

positions and velocities have shape (agents,2); initial positions lie
inside [0,box_size]. Returns positions (agents,2,steps+1). Reflection
handles arbitrarily many wall crossings in a single step.

### codegen

#### `state_hard_margins`

```python
state_hard_margins(s, params)
```

Clip s[0:2] to [0,1] and contact rate to its configured bounds.

#### `obs_hard_margins`

```python
obs_hard_margins(x, params)
```

Leave predicted observations unchanged (original Coder convention).

#### `nlin_state_update`

```python
nlin_state_update(u, s, w_bar, params)
```

Return optimal control and the next six-state SI-alpha/costate vector.

#### `nlin_obs_update`

```python
nlin_obs_update(u, s, v_bar, params)
```

Evaluate incidence plus v_bar; the standalone Coder model observes NEWCASES only.

#### `state_jacobians`

```python
state_jacobians(u, s, w_bar, params)
```

Return analytic state and process-noise Jacobians (6x6 each).

#### `obs_jacobian`

```python
obs_jacobian(u, s, v_bar, params)
```

Return observation and noise Jacobians (1x6 and 1x1).

#### `state_hessian_terms`

```python
state_hessian_terms(u, s, covariance, w_bar, q_w, params)
```

Return the four original zero-valued state/noise correction arrays.

#### `obs_hessian_terms`

```python
obs_hessian_terms(u, s, covariance, v_bar, r_v, params)
```

Return the four original zero-valued observation correction arrays.

#### `new_case_ekf_estimator_with_optimal_npi`

```python
new_case_ekf_estimator_with_optimal_npi(*args, **kwargs)
```

Return the dedicated estimator in MATLAB Coder's alternate output order.

Same arguments as the public dedicated estimator. Coder orders filtered
covariances and gain before smoothed mean/covariance; see docs/api.md.

### Geographic-table convenience reader

`read_geo_table(path)` is available in both languages. It preserves CSV column
names and normalizes missing RegionName entries to empty strings. MATLAB returns
a table; Python returns a pandas DataFrame. See `read_oxford_data` for date parsing
and duplicate-key validation.

The standalone Coder model differs from the Tools estimator: its observation
margin is a pass-through and it always uses NEWCASES observations. The Python
`codegen` module preserves these differences.
