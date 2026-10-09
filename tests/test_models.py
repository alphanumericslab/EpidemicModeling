"""Independent invariants and analytic checks for deterministic compartment solvers."""

import numpy as np
import pytest
import epidemic_modeling as em


def test_seirp_initial_and_mass():
    """Check seirp initial and mass."""

    states = np.array(
        em.seirp(
            0.4, 0.3, 0.12, 0.04, 0.09, 0.002, 0, 0.998, 0.001, 0.001, 0, 0, 100, 0.1
        )
    )
    np.testing.assert_array_equal(states[:, 0], [0.998, 0.001, 0.001, 0, 0])
    np.testing.assert_allclose(states.sum(axis=0), 1, atol=1e-12)
    assert states.shape == (5, 1000) and states.min() >= 0
    assert np.all(np.diff(states[4]) >= 0)


def test_no_transmission_closed_form_euler():
    """Check no transmission closed form euler."""

    dt = 0.1
    count = 100
    states = np.array(
        em.seirp(0, 0, 0, 0, 0.1, 0.02, 0, 0.9, 0, 0.1, 0, 0, count * dt, dt)
    )
    expected = 0.1 * (1 - dt * 0.12) ** np.arange(count)
    np.testing.assert_allclose(states[2], expected)

    np.testing.assert_allclose(states[4], (0.1 - expected) * (0.02 / 0.12))


def test_saturated_equal_rates_reduces_to_seirp():
    """Check saturated equal rates reduces to seirp."""

    ordinary = em.seirp(
        0.4, 0.3, 0.12, 0.04, 0.09, 0.002, 0, 0.998, 0.001, 0.001, 0, 0, 20, 0.1
    )
    saturated = em.seirp_saturated_resource(
        0.4,
        0.3,
        0.12,
        0.04,
        0,
        0.998,
        0.001,
        0.001,
        0,
        0,
        20,
        0.1,
        0.09,
        0.09,
        0.002,
        0.002,
        0.01,
        0.04,
    )
    np.testing.assert_allclose(saturated, ordinary)


def test_si_closed_form_and_endpoint():
    """Check si closed form and endpoint."""

    s, i = em.si_controlled(0, 0.1, 0.9, 0.1, 100, 0.1)
    np.testing.assert_array_equal(s, np.full(100, 0.9))
    np.testing.assert_allclose(i, 0.1 * 0.99 ** np.arange(100))
    assert len(em.seirp(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0.25, 0.1)[0]) == 3


def test_shared_noise_si_alpha_and_equilibrium():
    """Check shared noise si alpha and equilibrium."""

    u = np.ones((2, 30))
    noise = np.zeros((3, 30))
    states = em.si_alpha_controlled(
        u,
        0.999,
        0.001,
        0.3,
        [3, 3],
        0,
        5,
        1 / 7,
        [0.06, 0.04],
        0.1,
        0.1,
        0,
        0,
        0,
        30,
        1,
        noise=noise,
    )
    np.testing.assert_allclose(states[2], 0.3)
    assert len(states[0]) == 30 and states[0][0] < 0.999

    with pytest.raises(ValueError):
        em.si_alpha_controlled(
            u,
            0.999,
            0.001,
            0.3,
            [3, 3],
            0,
            5,
            1 / 7,
            [0.06, 0.04],
            0.1,
            0.1,
            0,
            0,
            0,
            30,
            1,
            noise=np.zeros((3, 29)),
        )


@pytest.mark.parametrize("duration,dt", [(0, 0.1), (1, 0), (-1, 0.1), (1, -1)])
def test_bad_grids(duration, dt):
    """Check bad grids."""

    with pytest.raises(ValueError):
        em.seirp(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, duration, dt)
