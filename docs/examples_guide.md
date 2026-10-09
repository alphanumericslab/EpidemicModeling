# Examples guide

**Reza Sameni · Emory University**

The notebooks demonstrate compartment models, growth estimation, Kalman filtering
and smoothing, intervention control, historical-data analysis, spatial models,
and neural activation layers. Each includes equations, parameter descriptions,
figures, numerical checks, further experiments, and source-paper references.

## Run the examples

Install `-e ".[notebooks,dev]"` in the environment used by Jupyter, then run the
notebook cells in order. The numerical routines are imported from
`src/epidemic_modeling`. HTML previews are available in `reports/notebooks/`.

MATLAB examples use the `demo_*` functions under `matlab/examples/`. Call
`setup_paths` before running them. `run_demo_checks` executes all six examples
without interactive windows and exports figures for review.

## Further experiments

- Change the Euler step and inspect conservation and convergence.
- Compare causal and centered growth estimates.
- Remove observation samples and inspect Kalman gain and smoothing behavior.
- Change intervention costs and inspect policy bounds and forecast tradeoffs.
- Evaluate historical forecasts at multiple held-out cutoffs.
- Inspect diffusion stability, reflecting boundaries, and supplied layer weights.

## Simulation recordings

The spatial notebook embeds both previous-run recordings: the M4V file and
the H.264 MP4 compressed from the larger AVI. Both are stored in `figures/`, with provenance in
`reports/previous_run_videos.json`.

## Troubleshooting

If imports fail, install the package in the notebook kernel's environment.
Preserve the repository's `figures/` and `data/` directories. Missing MATLAB
`nlinfit` or neural-layer classes require the respective optional toolboxes.
Reduce the time step if Euler produces negative compartment fractions; enforce
`D*dt/spacing**2 <= .25` for the explicit diffusion solver.
