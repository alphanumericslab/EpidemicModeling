# Migration and scope

The replacement root contains `src/`, `matlab/`, `notebooks/`, `tests/`, `data/`,
`figures/`, `docs/`, and `reports/`. Original files outside `update/` were not
edited. Copy the contents of `update/` or extract the replacement archive into
the intended repository root. Do not copy `.venv` or build/cache folders.

## Conversion coverage

All 22 original top-level MATLAB functions/classes have Python counterparts.
The original name-to-snake-case dictionary is in `name_map.json`.

| Original family | Updated implementation | Behavior |
| --- | --- | --- |
| SEIRP, SEIRPSaturatedResource | models.py; matched MATLAB functions | Original equations and Euler indexing; scalar rates added |
| SI_Controlled, SIalpha_Controlled | models.py; matched MATLAB functions | Original clipping and endpoints; explicit shared noise added |
| Rt_ExpFitGenRatios/LogLinReg/NonlinLS | growth.py; matched MATLAB functions | Original windows, endpoints, growth definitions preserved |
| Rt_ExpFitEKF | kalman.py; matched MATLAB function | First/second-order filter and smoother; zero linear-observation Hessian fixed |
| GenericExtendedKalmanFilter | kalman.py; matched MATLAB function | Callback EKF/EKS; covariance-stack handling fixed |
| Four SIAlphaModel variants | kalman.py; matched MATLAB functions | Forward/backward three/six-state dynamics and switching conventions preserved |
| NewCaseEKFEstimatorWithOptimalNPI | kalman.py and codegen.py; MATLAB equivalents | Simple covariance and phi>=0 convention preserved |
| NPICost, ReadCOVID19Data | npi.py, data.py; matched MATLAB functions | Original aggregation/cost semantics, robust missing-threshold indices |
| expLayer, MyTanhLayer | layers.py; snake-case MATLAB classes | Exact supplied-weight forward formulas; optional trainable PyTorch factories |
| TrainNPIPrescriptor | npi.py; train_npi_prescriptor.m | Redesigned matching deterministic pipeline |
| TrainPredictPrescribeNPI | npi.py; train_predict_prescribe_npi.m | Redesigned matching explicit result-returning pipeline |
| PrescribeNPI | npi.py; prescribe_npi.m | Explicit model path, shared JSON, reproducible prescription CSV |
| ForecastQualityAssessment | npi.py; forecast_quality_assessment.m | Held-out forecasts and explicit metrics without future-case fitting |
| MATLAB Coder standalone helpers | codegen.py; matlab/codegen | Analytic helper equations, NEWCASES-only observation and alternate ordering |

Nested SI callback equivalents are shared closures in `kalman.py`; forward and
reverse variants use the same analytic equations with opposite Euler signs.
The learnable NumPy functions take explicit weights rather than silently using
an unrelated random initializer. MATLAB retains its `Alpha` properties and
`predict` method required by Deep Learning Toolbox. Import optional trainable
Python modules via `torch_exp_layer` / `torch_my_tanh_layer`.

## Deliberate corrections and modernization

1. Generic covariance stacks are distinguished from fixed square matrices;
   second-order callbacks receive the current covariance slice.
2. Reverse-time results flip their actual final time dimension after `rho`
   is squeezed. Time-dependent noise covariance sequences are also reversed.
3. Terminal covariance assignments use paired indices rather than a Cartesian
   row/column submatrix.
4. The exponential model's linear observation has zero Hessian corrections;
   the historical dimensional mismatch is removed.
5. Smoothing uses a pseudoinverse in both languages for numerical robustness.
6. The original four large XPRIZE orchestration functions become a matching
   deterministic two-pass EKF/NNLS workflow in both languages. The refined
   regression uses refined data. Observation adaptation is disabled in this
   training routine for reproducibility; a positive variance floor
   and a finite contact-rate bound prevent degenerate synthetic experiments.
   The generic filters still expose the original adaptation options.
7. Pipeline I/O uses explicit paths, daily-gap validation, shared JSON models,
   returned forecast arrays and held-out metrics. The prescription solver
   uses finite-iteration EKF/EKS shooting with bounded policies and zero
   terminal costates. It does not certify a global optimum.
8. Daily incidence is calculated before each forecast transition, avoiding
   a one-day shift when using a solver that excludes the initial state.
9. Nonlinear zero-count fallback applies consistently to centered even-window
   fits as well as causal fits. Both languages reject invalid grids and
   nonpositive log-regression inputs explicitly.
10. Missing JHU threshold indices are zero in both languages. Date columns
   must align before aggregation.

The upgraded pipelines are **not literal reproductions of every historical
plotting branch, LASSO option, saved MAT-file format, or country-specific research
configuration**. Their new implementations match each other. Their original
source is preserved in `docs/historical_matlab/*_original.m.txt` for review.
Historical exploratory scripts are likewise retained as text, with snake-case
call sites, because many depended on external cloned datasets and toolboxes.
The seven executed notebooks and six MATLAB examples replace those
scripts as the supported class entry points. Optional LSTM exploratory training
scripts are preserved for provenance rather than presented as validated examples;
the custom layer formulas and trainable Python counterparts are provided.

Existing figures, EPS/FIG files, and videos are copied unchanged to `figures/`.
The historical CSV snapshots are copied unchanged to `data/`. No new research
results are attributed to the original papers. Figures displayed in notebooks
are labeled as historical illustrations when the new experiment differs.
