"""Pipeline integration, data semantics, and spatial invariants."""

from pathlib import Path
import json
import numpy as np
import pandas as pd

import pytest
import epidemic_modeling as em


def synthetic_training():
    """Create a deterministic intervention schedule and case counts for fitting tests."""

    day = np.arange(50)
    u = np.vstack((np.where(day < 25, 0, 2), np.where(day < 35, 0, 2)))
    p = em.default_si_params([3, 3])
    p.update(a=np.array([0.06, 0.04]), b=0.02)
    truth = dict(
        population=1e6,
        params=p,
        state=[0.9999, 0.0001, 0.3],
        covariance=np.diag([1e-8, 1e-8, 0.01]),
    )

    cases, _ = em.forecast_npi(truth, u)

    return np.cumsum(cases), u


def test_nnls_kkt_conditions():
    """Check nonnegative coefficients and least-squares optimality conditions."""

    rng = np.random.default_rng(2)
    x = rng.uniform(size=(100, 4))
    y = x @ np.array([0.1, 0, 0.3, 0.2])
    a = em.nonnegative_least_squares(x, y)
    np.testing.assert_allclose(a, [0.1, 0, 0.3, 0.2], atol=1e-7)

    assert np.linalg.norm(x @ a - y) < 1e-6


def test_training_forecast_and_control_bounds():
    """Check training forecast and control bounds."""

    c, u = synthetic_training()
    model = em.fit_npi_model(c, u, 1e6, [3, 3])
    assert np.all(np.array(model["refined_coefficients"]) >= 0)
    json.dumps(model, allow_nan=False)
    controls, cases, states = em.optimal_npi(model, 15, [1, 1.5], 0.3)

    assert controls.shape == (2, 15) and states.shape == (3, 15)
    assert np.all((controls >= 0) & (controls <= 3)) and np.all(cases >= 0)
    assert cases[0] == pytest.approx(1e6 * np.prod(model["state"]))
    human, cost = em.npi_cost([10, 20], [[1, 2], [2, 3]], [1, 2])
    assert human == 15 and cost == 3.25


def test_data_reader_and_csv_pipeline(tmp_path):
    """Check data reader and csv pipeline."""

    c, u = synthetic_training()
    dates = pd.date_range("2020-01-01", periods=50)
    data = pd.DataFrame(
        {
            "CountryName": "Example",
            "RegionName": "",
            "Date": dates.strftime("%Y%m%d"),
            "ConfirmedCases": c,
            "NPI A": u[0],
            "NPI B": u[1],
        }
    )
    data_path = tmp_path / "data.csv"
    data.to_csv(data_path, index=False)

    geo = tmp_path / "geo.csv"
    pd.DataFrame({"CountryName": ["Example"], "RegionName": [""]}).to_csv(
        geo, index=False
    )
    pop = tmp_path / "pop.csv"
    pd.DataFrame(
        {"CountryName": ["Example"], "RegionName": [""], "Population2020": [1e6]}
    ).to_csv(pop, index=False)
    costs = tmp_path / "cost.csv"

    pd.DataFrame(
        {"CountryName": ["Example"], "RegionName": [""], "NPI A": [1], "NPI B": [1.5]}
    ).to_csv(costs, index=False)
    model_path = tmp_path / "model.json"
    bundle = em.train_npi_prescriptor(
        "2020-01-01",
        "2020-02-09",
        data_path,
        geo,
        pop,
        ["NPI A", "NPI B"],
        [3, 3],
        model_path,
    )
    assert len(bundle["models"]) == 1 and bundle["models"][0]["training_days"] == 40
    output = em.prescribe_npi(
        "2020-02-10",
        "2020-02-14",
        geo,
        costs,
        tmp_path / "out.csv",
        model_file=model_path,
    )

    assert len(output) == 5 and set(output.columns) == {
        "CountryName",
        "RegionName",
        "Date",
        "NPI A",
        "NPI B",
    }
    scores = em.forecast_quality_assessment(
        [1, 1.5],
        0.3,
        "2020-01-01",
        "2020-02-09",
        "2020-01-01",
        "2020-02-19",
        10,
        data_path,
        geo,
        pop,
        ["NPI A", "NPI B"],
        [0, 0],
        [3, 3],
        model_path,
    )
    assert scores.days.iloc[0] == 10 and np.isfinite(scores.rmse.iloc[0])
    results = em.train_predict_prescribe_npi(
        [1, 1.5],
        0.3,
        "2020-01-01",
        "2020-02-09",
        "2020-01-01",
        "2020-02-19",
        data_path,
        geo,
        pop,
        ["NPI A", "NPI B"],
        [0, 0],
        [3, 3],
        model_path,
    )
    assert len(results) == 1 and results[0]["fixed_cases"].shape == (10,)


def test_johns_hopkins_aggregation_and_missing_threshold(tmp_path):
    """Check johns hopkins aggregation and missing threshold."""

    files = []

    for name, values in [
        ("cases", [[0, 1, 3], [0, 2, 4]]),
        ("deaths", [[0, 0, 1], [0, 0, 0]]),
        ("recovered", [[0, 0, 0], [0, 1, 1]]),
    ]:
        frame = pd.DataFrame(
            {
                "Province/State": ["A", "B"],
                "Country/Region": ["Example", "Example"],
                "Lat": [0, 0],
                "Long": [0, 0],
                "1/1/20": [v[0] for v in values],
                "1/2/20": [v[1] for v in values],
                "1/3/20": [v[2] for v in values],
            }
        )
        path = tmp_path / f"{name}.csv"
        frame.to_csv(path, index=False)
        files.append(path)

    total, infected, recovered, deceased, first, threshold, days = em.read_covid19_data(
        *files, ["Example", "Absent"], 5
    )
    np.testing.assert_array_equal(total, [[0, 3, 7], [0, 0, 0]])
    np.testing.assert_array_equal(first, [2, 0])
    np.testing.assert_array_equal(threshold, [3, 0])
    assert days == 3

    np.testing.assert_array_equal(infected, total - recovered - deceased)


def test_preprocessing_causal():
    """Check preprocessing causal."""

    daily, smooth = em.prepare_cases([0, 4, 3, 8], 2)
    np.testing.assert_array_equal(daily, [0, 4, 0, 5])
    np.testing.assert_array_equal(smooth, [0, 2, 2, 2.5])


def test_diffusion_mass_nonnegative_and_stability():
    """Check diffusion mass nonnegative and stability."""

    x = np.zeros((21, 21))
    x[10, 10] = 1
    result = em.diffusion_2d(x, 1, 0.2, 1, 40)
    np.testing.assert_allclose(result.sum(axis=(0, 1)), 1)
    assert result.min() >= 0 and result[10, 10, -1] < 1

    with pytest.raises(ValueError):
        em.diffusion_2d(x, 1, 0.3, 1, 2)


def test_reflection_and_layers():
    """Check reflection and layers."""

    result = em.population_motion_2d([[0.2, 0.3]], [[10, -20]], 1, 3)
    assert np.all((result >= 0) & (result <= 1))
    np.testing.assert_allclose(result[:, :, 0], [[0.2, 0.3]])
    np.testing.assert_allclose(em.exp_layer([0, 1], 0.7), [1, np.exp(0.7)])
    np.testing.assert_allclose(
        em.my_tanh_layer([-1, 0, 1], 0.7), 0.7 * np.tanh(np.array([-1, 0, 1]) / 0.7)
    )
