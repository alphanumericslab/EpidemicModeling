# Validation evidence

Validated on **October 9, 2026** with Python 3.12 on macOS arm64 and MATLAB R2024b.
The original repository files remain outside the replacement folder unchanged.

| Check | Result |
| --- | --- |
| Python tests | 35 passed |
| MATLAB function-based tests | 18 passed; 0 failed; 0 incomplete |
| Notebook execution | All 7 examples passed in separate fresh kernels |
| MATLAB examples | All 6 ran headlessly and rendered 10 figures |
| Actual MATLAB/Python numerical comparison | 124 arrays passed |
| Python package build | Wheel and source distribution built successfully |

The MATLAB tests exercised CSV training, JSON model loading, daily prescription
output, held-out scoring, JHU aggregation, conservation, known exponential growth,
scalar Kalman recursion, missing observations, NNLS, diffusion/reflection, and
learnable layer forward passes. Statistics and Machine Learning Toolbox and
Deep Learning Toolbox tests both ran successfully in this installation.

The Python suite additionally checks analytic SI/costate Jacobians against finite
differences, conversion coverage of all 22 original public functions/classes,
snake-case function names, documentation, and notebook imports/attribution.
All deterministic parity cases use shared supplied inputs and noise.

## Cross-language comparison tolerances

The baseline uses `rtol=2e-6`, `atol=2e-9`. Nonlinear regression uses `rtol=2e-4`,
`atol=2e-6` because MATLAB nlinfit and SciPy least_squares have different optimizer
internals. Two sensitive smoother outputs use documented absolute tolerances:
`legacy_s_smooth`: 1e-7; `si_backward_control_p_smooth`: 1e-6. The latter concerns
costate covariances in the reverse-time model; MATLAB and NumPy pseudoinverse
roundoff is amplified by the smoothing recursion. All other fields retain the
baseline threshold. The comparison report records each array's actual maximum
absolute difference and threshold. Matching is numerical, not bitwise.

Actual MATLAB results are stored in `reports/matlab_results.mat`; the detailed
comparison is in `reports/matlab_results.comparison.json`. The MATLAB execution
log is `reports/matlab_validation.log`. Notebook HTML previews are under
`reports/notebooks/`; MATLAB PNGs are under `reports/matlab_figures/`.

## Reproduce

```bash
python -m pytest
python scripts/execute_notebooks.py
```

```matlab
addpath('matlab'); root = setup_paths();
results = runtests(fullfile(root,'matlab','tests'));
assert(all([results.Passed]));
run_demo_checks;
export_parity_results;
```

```bash
python scripts/compare_matlab_python.py reports/matlab_results.mat
python -m build
python scripts/package_release.py
```

The optional PyTorch trainable factories require the `neural` extra and are not
part of the validated base environment; NumPy activation formulas and MATLAB
learnable layer forward passes are tested. MATLAB Coder compilation and the
historical LSTM training scripts are not claimed as validated. Upgraded
orchestration replaces historical plotting/configuration branches as explained
in `migration.md`.
