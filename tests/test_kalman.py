"""Linear Kalman oracle, boundary conditions, missing data and analytic Jacobians."""

from types import SimpleNamespace
import numpy as np
import pytest
from epidemic_modeling.kalman import generic_extended_kalman_filter, _si_handles

import epidemic_modeling as em


def linear_handles():
    """Create scalar linear callbacks for an independent Kalman filter reference."""

    return dict(
        state_hard_margins=lambda s, p, k: s,
        obs_hard_margins=lambda x, p, k: x,
        nlin_state_update=lambda u, s, w, p, k: (u, s + w),
        nlin_obs_update=lambda u, s, v, p, k: s + v,
        state_jacobians=lambda u, s, w, p, k: (np.eye(1), np.eye(1)),
        obs_jacobian=lambda u, s, v, p, k: (np.eye(1), np.eye(1)),
        state_hessian_terms=lambda *args: (
            np.zeros(1),
            np.zeros((1, 1)),
            np.zeros(1),
            np.zeros((1, 1)),
        ),
        obs_hessian_terms=lambda *args: (
            np.zeros(1),
            np.zeros((1, 1)),
            np.zeros(1),
            np.zeros((1, 1)),
        ),
    )


def test_linear_filter_against_scalar_oracle():
    """Compare the generic filter with independent scalar Kalman equations."""

    x = np.array([1.0, 2.0, np.nan, 1.5])
    r = generic_extended_kalman_filter(
        np.zeros((1, 4)),
        x,
        linear_handles(),
        None,
        [0.0],
        [[1.0]],
        [np.nan],
        [[np.nan]],
        [0.0],
        0.0,
        0.1,
        0.5,
        1,
        1,
        3,
        1,
    )
    mean, var = 0.0, 1.0

    for k, value in enumerate(x):
        if np.isfinite(value):
            gain = var / (var + 0.5)
            mean += gain * (value - mean)
            var = (1 - gain) * var

        np.testing.assert_allclose(r.s_plus[0, k], mean)
        np.testing.assert_allclose(r.p_plus[0, 0, k], var)
        var += 0.1

    assert np.all(r.k_gain[:, :, 2] == 0)
    assert np.all(r.p_smooth <= r.p_plus + 1e-12)


def test_time_dependent_covariance_and_terminal_constraints():
    """Check time dependent covariance and terminal constraints."""

    r = generic_extended_kalman_filter(
        np.zeros((1, 4)),
        [1, 2, 3, 4],
        linear_handles(),
        None,
        [0.0],
        [[1.0]],
        [9.0],
        [[0.01]],
        [0.0],
        0.0,
        np.arange(1, 5) * 0.1,
        np.arange(1, 5) * 0.2,
        1,
        1,
        3,
        2,
    )
    assert r.s_smooth[0, -1] == 9 and r.p_smooth[0, 0, -1] == 0.01
    assert np.all(np.isfinite(r.s_plus))


@pytest.mark.parametrize(
    "controlled,backward", [(False, False), (False, True), (True, False), (True, True)]
)
def test_si_jacobians_match_finite_difference(controlled, backward):
    """Compare analytic SI/costate Jacobians with finite-difference derivatives."""

    p = SimpleNamespace(**em.default_si_params([3, 3]))
    p.a = np.array([0.06, 0.04])
    p.b = 0.02
    s = (
        np.array([0.8, 0.1, 0.3, 0.05, 0.02, 0.03])
        if controlled
        else np.array([0.8, 0.1, 0.3])
    )
    u = np.array([1.0, 2.0])

    handles = _si_handles(controlled, backward)
    a, b = handles.state_jacobians(u, s, np.zeros(len(s)), p, 0)
    numerical = np.zeros_like(a)
    delta = 1e-6

    for j in range(len(s)):
        plus = s.copy()
        minus = s.copy()
        plus[j] += delta
        minus[j] -= delta
        numerical[:, j] = (
            handles.nlin_state_update(u, plus, None, p, 0)[1]
            - handles.nlin_state_update(u, minus, None, p, 0)[1]
        ) / (2 * delta)

    np.testing.assert_allclose(a, numerical, atol=1e-8)


@pytest.mark.parametrize("order", [1, 2])
def test_exponential_filter_and_missing_observations(order):
    """Check exponential filter and missing observations."""

    truth = 25 * np.exp(0.02 * np.arange(40))
    x = truth + np.sin(np.arange(40))
    x[10:13] = np.nan
    r = em.rt_exp_fit_ekf(
        x,
        [25, 0.02],
        [1, 1, 0.2],
        [0, 0],
        0,
        np.diag([4, 0.001]),
        np.diag([1, 1e-5]),
        4,
        1,
        1,
        7,
        order,
    )
    assert all(np.all(np.isfinite(a)) for a in r)

    assert np.all(r[4][:, :, 10:13] == 0)
    assert np.mean(abs(r[5][0] - truth)) < 2


def test_all_si_variants():
    """Check all si variants."""

    p = em.default_si_params([3])
    p.update(a=np.array([0.08]), b=0.02)
    u = np.ones((1, 20))
    y = np.full(20, 0.0003)

    for controlled in (False, True):
        m = 6 if controlled else 3
        initial = np.r_[0.999, 0.001, 0.3, np.zeros(m - 3)]
        cov = np.diag(np.r_[1e-6, 1e-6, 0.01, np.ones(m - 3)])
        terminal = initial.copy()
        final_cov = cov.copy()

        args = [
            u,
            y,
            p,
            initial,
            cov,
            terminal,
            final_cov,
            np.zeros(m),
            0,
            np.eye(m) * 1e-8,
            1e-8,
            1,
            1,
            7,
            1,
        ]
        forward = (
            em.si_alpha_model_ekf_opt_controlled
            if controlled
            else em.si_alpha_model_ekf
        )
        backward = (
            em.si_alpha_model_backward_ekf_opt_controlled
            if controlled
            else em.si_alpha_model_backward_ekf
        )

        for fun in (forward, backward):
            r = fun(*args)
            assert r.s_plus.shape == (m, 20) and np.all(np.isfinite(r.s_plus))

        if controlled:
            old = em.new_case_ekf_estimator_with_optimal_npi(*args)
            assert len(old) == 10
