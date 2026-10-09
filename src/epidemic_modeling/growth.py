"""Rolling growth estimates with the original MATLAB endpoint conventions.

Author: Reza Sameni, Emory University. See Sameni (2020, 2022).
"""

import numpy as np
from scipy.signal import lfilter
from scipy.optimize import least_squares


def _cases(new_cases, wlen, time_unit):
    """Validate a one-dimensional nonnegative case series and fitting grid."""
    x = np.asarray(new_cases, float).reshape(-1)

    if not np.all(np.isfinite(x)) or np.any(x < 0):
        raise ValueError("new_cases must be finite and nonnegative")

    if int(wlen) != wlen or wlen < 2 or time_unit <= 0:
        raise ValueError("wlen must be an integer >= 2 and time_unit positive")

    return x


def exp_model(params, t):
    """Evaluate ``amplitude * exp(growth * t)`` for a two-element parameter vector.

    params contains [amplitude, growth]. t contains sample times in the same
    time unit as growth. Returns an array with the shape of t.
    """

    return params[0] * np.exp(params[1] * np.asarray(t))


def rt_exp_fit_gen_ratios(new_cases, wlen, generation_period, time_unit):
    """Estimate lagged growth and its causal zero-padded moving average.

    Returns (rt, growth, rt_smoothed, growth_smoothed), each length N.
    growth[:generation_period] = 0; subsequent growth is log case ratio
    divided by generation_period. rt = exp(growth*time_unit). Zero counts
    propagate IEEE NaN/Inf, preserving the original definition.

    new_cases is a length-N nonnegative series. wlen is the averaging
    window in samples; generation_period is the lag in samples. time_unit
    sets the interval used to turn the growth estimate into a growth factor.
    """
    x = _cases(new_cases, wlen, time_unit)

    if int(generation_period) != generation_period or not 1 <= generation_period <= len(
        x
    ):
        raise ValueError("generation_period must be an integer in [1, N]")

    with np.errstate(divide="ignore", invalid="ignore"):
        growth = np.r_[
            np.zeros(generation_period),
            np.log(x[generation_period:] / x[:-generation_period]) / generation_period,
        ]

    smooth = lfilter(np.ones(wlen) / wlen, [1], growth)

    return np.exp(growth * time_unit), growth, np.exp(smooth * time_unit), smooth


def rt_exp_fit_log_lin_reg(new_cases, wlen, time_unit, causal=True):
    """Fit rolling log cases using closed-form linear least squares.

    Returns (rt, amplitude, growth, exp_fit). Causal windows end at each
    sample; centered windows have 2*floor(wlen/2)+1 points. Unfitted ends
    have rt=amplitude=exp_fit=1 and growth=0, matching MATLAB. rt=exp(slope),
    growth=slope/time_unit, and exp_fit=amplitude*rt (one-step forecast).
    Positive cases are required; zeros produce undefined log regressions.

    new_cases is a length-N daily case series. wlen is the fitting window
    in samples; time_unit is the time between samples. All four returned
    arrays have length N. causal=True uses only current and earlier cases.
    """
    x = _cases(new_cases, wlen, time_unit)

    if np.any(x <= 0):
        raise ValueError("log-linear fitting requires strictly positive cases")

    log_amplitude = np.zeros(len(x))
    slope = np.zeros(len(x))
    y = np.log(x)
    half = wlen // 2
    n = np.arange(-wlen + 1, 1) if causal else np.arange(-half, half + 1)

    en = n.mean()
    en2 = (n * n).mean()
    det = en2 - en * en

    # Fit only complete windows; preserve the original endpoint values.
    indices = range(wlen - 1, len(x)) if causal else range(half, len(x) - half)

    for k in indices:
        seg = y[k - wlen + 1 : k + 1] if causal else y[k - half : k + half + 1]
        log_amplitude[k] = (seg.mean() * en2 - (n * seg).mean() * en) / det
        slope[k] = ((n * seg).mean() - seg.mean() * en) / det

    amplitude = np.exp(log_amplitude)
    rt = np.exp(slope)

    return rt, amplitude, slope / time_unit, amplitude * rt


def rt_exp_fit_nonlin_ls(new_cases, wlen, time_unit, causal=True):
    """Fit rolling exponential cases by nonlinear least squares.

    Inputs and return order follow :func:`rt_exp_fit_log_lin_reg`. The fit
    uses n/time_unit as its independent variable, preserving the original
    unusual scaling: growth = fitted_rate/time_unit and rt=exp(fitted_rate).
    Causal unfitted amplitudes are delayed raw cases; centered ends retain
    raw cases. Windows with zeros use the current count and zero growth.

    new_cases is a length-N nonnegative case series. wlen is the fitting
    window in samples. Returns four length-N arrays. Use causal=True for
    forecasts that must not use later observations.
    """
    x = _cases(new_cases, wlen, time_unit)
    rate = np.zeros(len(x))
    amplitude = lfilter(np.r_[np.zeros(wlen - 1), 1.0], [1], x) if causal else x.copy()
    half = wlen // 2

    n = np.arange(-wlen + 1, 1) if causal else np.arange(-half, half + 1)

    # Fit only complete windows; preserve the original endpoint values.
    indices = range(wlen - 1, len(x)) if causal else range(half, len(x) - half)

    for k in indices:
        seg = x[k - wlen + 1 : k + 1] if causal else x[k - half : k + half + 1]

        if np.any(seg == 0):
            amplitude[k] = x[k]
            continue

        result = least_squares(
            lambda p: exp_model(p, n / time_unit) - seg,
            [x[k], 0],
            xtol=1e-6,
            ftol=1e-6,
            gtol=1e-6,
            max_nfev=250,
        )
        amplitude[k], rate[k] = result.x

    rt = np.exp(rate)

    return rt, amplitude, rate / time_unit, amplitude * rt
