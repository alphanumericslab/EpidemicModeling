"""Reproducible notebook execution, numerical parity fixtures, and release packaging."""

from pathlib import Path
import hashlib
import json
import sys

import zipfile
import numpy as np
from . import models, growth, kalman, npi, layers, spatial


def parity_results():
    """Compute deterministic results matching matlab/tests/export_parity_results.m."""
    out = {}
    out["seirp"] = np.array(
        models.seirp(
            0.4, 0.3, 0.12, 0.04, 0.09, 0.002, 0, 0.998, 0.001, 0.001, 0, 0, 10, 0.1
        )
    )
    out["saturated"] = np.array(
        models.seirp_saturated_resource(
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
            10,
            0.1,
            0.09,
            0.03,
            0.002,
            0.015,
            0.01,
            0.04,
        )
    )
    out["si"] = np.array(models.si_controlled(0.35, 0.1, 0.999, 0.001, 100, 0.1))

    u = np.vstack((np.linspace(0, 3, 40), np.linspace(3, 0, 40)))
    noise = np.sin(np.arange(120).reshape(3, 40, order="F"))
    out["si_alpha"] = np.array(
        models.si_alpha_controlled(
            u,
            0.999,
            0.001,
            0.3,
            [3, 3],
            0,
            5,
            1 / 7,
            [0.06, 0.04],
            0.02,
            0.1,
            0.001,
            0.001,
            0.01,
            40,
            0.1,
            noise=noise,
        )
    )
    cases = 25 * np.exp(0.035 * np.arange(40))
    out["ratios"] = np.array(growth.rt_exp_fit_gen_ratios(cases, 7, 5, 1))

    out["log_fit"] = np.array(growth.rt_exp_fit_log_lin_reg(cases, 7, 1))
    out["log_fit_centered"] = np.array(
        growth.rt_exp_fit_log_lin_reg(cases, 8, 1, False)
    )
    out["nonlinear_fit"] = np.array(growth.rt_exp_fit_nonlin_ls(cases, 7, 1))
    out["nonlinear_fit_centered"] = np.array(
        growth.rt_exp_fit_nonlin_ls(cases, 8, 1, False)
    )
    x = cases + 2 * np.sin(np.arange(40))

    x[15:18] = np.nan

    for order in (1, 2):
        r = kalman.rt_exp_fit_ekf(
            x,
            [25, 0.03],
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

        for name, value in zip(
            (
                "s_minus",
                "s_plus",
                "p_minus",
                "p_plus",
                "gain",
                "s_smooth",
                "p_smooth",
                "innovations",
                "rho",
            ),
            r,
        ):
            out[f"exp_{order}_{name}"] = value

    p = npi.default_si_params([3, 3])
    p.update(a=np.array([0.06, 0.04]), b=0.02)
    y = np.full(40, 0.0003)
    y[15:18] = np.nan
    common = [
        p,
        [0.999, 0.001, 0.3],
        np.diag([1e-6, 1e-6, 0.01]),
        np.full(3, np.nan),
        np.full((3, 3), np.nan),
        np.zeros(3),
        0,
        np.diag([1e-8, 1e-8, 1e-5]),
        1e-8,
        1,
        1,
        7,
        1,
    ]

    r = kalman.si_alpha_model_ekf(u, y, *common)

    for name, value in zip(r._fields, r):
        out["si_filter_" + name] = value

    q = np.repeat(np.diag([1e-8, 1e-8, 1e-5])[:, :, None], 40, axis=2)
    q *= np.linspace(1, 2, 40)
    variance = np.linspace(1e-8, 2e-8, 40)
    variable = common.copy()
    variable[7:9] = [q, variance]

    r = kalman.si_alpha_model_ekf(u, y, *variable)

    for name, value in zip(r._fields, r):
        out["si_variable_" + name] = value

    back = common.copy()
    back[3:5] = [[0.995, 0.002, 0.3], np.diag([1e-6, 1e-6, 0.01])]
    r = kalman.si_alpha_model_backward_ekf(u, y, *back)

    for name, value in zip(r._fields, r):
        out["si_backward_" + name] = value

    control_common = [
        p,
        [0.999, 0.001, 0.3, 0, 0, 0],
        np.diag([1e-6, 1e-6, 0.01, 1, 1, 1]),
        np.r_[np.full(3, np.nan), np.zeros(3)],
        np.full((6, 6), np.nan),
        np.zeros(6),
        0,
        np.diag([1e-8, 1e-8, 1e-5, 1e-4, 1e-4, 1e-4]),
        1e-8,
        1,
        1,
        7,
        1,
    ]
    r = kalman.si_alpha_model_ekf_opt_controlled(
        np.full_like(u, np.nan), np.full(40, np.nan), *control_common
    )

    for name, value in zip(r._fields, r):
        out["si_control_" + name] = value

    backwards = control_common.copy()
    backwards[3:5] = [
        [0.995, 0.002, 0.3, 0, 0, 0],
        np.diag([1e-6, 1e-6, 0.01, 1, 1, 1]),
    ]
    r = kalman.si_alpha_model_backward_ekf_opt_controlled(u, y, *backwards)

    for name, value in zip(r._fields, r):
        out["si_backward_control_" + name] = value

    legacy = kalman.new_case_ekf_estimator_with_optimal_npi(u, y, *control_common)

    for name, value in zip(
        (
            "u_opt",
            "s_minus",
            "s_plus",
            "s_smooth",
            "p_minus",
            "p_plus",
            "p_smooth",
            "k_gain",
            "innovations",
            "rho",
        ),
        legacy,
    ):
        out["legacy_" + name] = value

    out["costs"] = np.array(npi.npi_cost(cases, u, [1, 1.5]))
    out["exp_layer"] = layers.exp_layer(np.linspace(-2, 2, 21), 0.7)
    out["tanh_layer"] = layers.my_tanh_layer(np.linspace(-2, 2, 21), 0.7)
    initial = np.zeros((11, 11))
    initial[5, 5] = 1

    out["diffusion"] = spatial.diffusion_2d(initial, 1, 0.2, 1, 10)
    out["motion"] = spatial.population_motion_2d(
        [[0.2, 0.3], [0.7, 0.8]], [[0.14, 0.09], [-0.11, 0.07]], 0.1, 30
    )
    design = np.column_stack(
        (np.ones(40), np.arange(40) / 40, np.sin(np.arange(40)) ** 2)
    )
    out["nnls"] = npi.nonnegative_least_squares(
        design, design @ np.array([0.1, 0.2, 0.3])
    )
    train_u = np.vstack(
        (np.where(np.arange(50) < 25, 0, 2), np.where(np.arange(50) < 35, 0, 2))
    )

    truth = dict(
        population=1e6,
        params=p,
        state=[0.9999, 0.0001, 0.3],
        covariance=np.diag([1e-8, 1e-8, 0.01]),
    )
    train_cases, _ = npi.forecast_npi(truth, train_u)
    model = npi.fit_npi_model(np.cumsum(train_cases), train_u, 1e6, [3, 3])
    out["trained_coef"] = np.array(model["coefficients"])
    out["trained_refined"] = np.array(model["refined_coefficients"])

    out["trained_state"] = np.array(model["state"])
    out["trained_cov"] = np.array(model["covariance"])
    out["forecast"], out["forecast_states"] = npi.forecast_npi(
        model, np.repeat(train_u[:, -1, None], 10, axis=1)
    )
    out["optimal_controls"], out["optimal_cases"], out["optimal_states"] = (
        npi.optimal_npi(model, 10, [1, 1.5], 0.3)
    )
    from . import codegen

    state = np.array([0.8, 0.1, 0.3, 0.05, 0.02, 0.03])
    control = np.array([1.0, 2.0])
    out["coder_state_margin"] = codegen.state_hard_margins(state, p)
    out["coder_obs_margin"] = codegen.obs_hard_margins(np.array([-0.1]), p)
    out["coder_u"], out["coder_state"] = codegen.nlin_state_update(
        control, state, np.zeros(6), p
    )

    out["coder_obs"] = codegen.nlin_obs_update(control, state, 0, p)
    out["coder_a"], out["coder_b"] = codegen.state_jacobians(
        control, state, np.zeros(6), p
    )
    out["coder_c"], out["coder_d"] = codegen.obs_jacobian(control, state, 0, p)

    for key, value in zip(
        ("fs", "cs", "fw", "cw"),
        codegen.state_hessian_terms(
            control, state, np.eye(6), np.zeros(6), np.eye(6), p
        ),
    ):
        out["coder_" + key] = value

    for key, value in zip(
        ("gs", "gsp", "gv", "gvp"),
        codegen.obs_hessian_terms(control, state, np.eye(6), 0, 1, p),
    ):
        out["coder_" + key] = value

    return out


def compare_matlab_python(path, *, rtol=2e-6, atol=2e-9):
    """Compare actual MATLAB exported results against the matching Python cases.

    Missing fields or numerical mismatches raise AssertionError. Nonlinear
    regression uses looser tolerances because the optimizer internals differ.

    path is the MAT file produced by export_parity_results in MATLAB.
    Writes a JSON comparison report next to it and returns the number of
    compared arrays. Raises AssertionError if any comparison fails.
    """
    from scipy.io import loadmat

    actual = loadmat(path)
    expected = parity_results()
    failures = []

    skipped = []
    comparisons = {}

    for key, value in expected.items():
        if key not in actual:
            if key.startswith("nonlinear") and "optional_nonlinear_missing" in actual:
                skipped.append(key)
                continue

            failures.append(f"{key}: missing MATLAB output")
            continue

        left = np.squeeze(actual[key])
        right = np.squeeze(value)
        relative = 2e-4 if key.startswith("nonlinear") else rtol
        absolute = 2e-6 if key.startswith("nonlinear") else atol

        if key == "legacy_s_smooth":
            absolute = 1e-7

        if key == "si_backward_control_p_smooth":
            absolute = 1e-6

        passed = left.shape == right.shape and np.allclose(
            left, right, rtol=relative, atol=absolute, equal_nan=True
        )
        comparisons[key] = dict(
            passed=bool(passed),
            rtol=relative,
            atol=absolute,
            shape=list(right.shape),
            max_absolute_difference=(
                float(np.max(np.abs(left - right)))
                if left.shape == right.shape
                else None
            ),
        )

        if not passed:
            failures.append(f"{key}: numerical/shape mismatch")

    report = dict(
        compared=len(comparisons),
        skipped=skipped,
        failures=failures,
        arrays=comparisons,
    )
    Path(path).with_suffix(".comparison.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )

    if failures:
        raise AssertionError("\n".join(failures))

    print(f"MATLAB/Python parity passed for {len(expected)-len(skipped)} arrays.")

    if skipped:
        print("Optional nonlinear regression skipped: MATLAB nlinfit unavailable.")

    return len(expected) - len(skipped)


def execute_notebooks(root=None):
    """Execute each notebook in a fresh kernel and export self-contained HTML previews.

    root is the project directory, or is found from the working directory.
    Runs every notebook using the current Python interpreter, saves its
    outputs, and writes HTML previews to reports/notebooks. Execution
    errors stop the run. Returns None.
    """
    import nbformat
    from nbclient import NotebookClient
    from nbconvert import HTMLExporter
    from .plotting import repository_root

    root = Path(root) if root else repository_root()
    reports = root / "reports/notebooks"
    reports.mkdir(parents=True, exist_ok=True)
    # A temporary kernelspec guarantees the interpreter used to invoke this command.
    import os, tempfile

    with tempfile.TemporaryDirectory(prefix="epidemic-kernel-") as temporary:
        kernels = Path(temporary) / "kernels/epidemic_validation"
        kernels.mkdir(parents=True)
        (kernels / "kernel.json").write_text(
            json.dumps(
                dict(
                    argv=[
                        sys.executable,
                        "-m",
                        "ipykernel_launcher",
                        "-f",
                        "{connection_file}",
                    ],
                    display_name="Epidemic validation",
                    language="python",
                )
            )
        )
        previous = os.environ.get("JUPYTER_PATH")
        os.environ["JUPYTER_PATH"] = temporary + (
            os.pathsep + previous if previous else ""
        )

        try:
            for path in sorted((root / "notebooks").glob("*.ipynb")):
                print(f"Executing {path.name}", flush=True)
                nb = nbformat.read(path, as_version=4)
                NotebookClient(
                    nb,
                    timeout=180,
                    kernel_name="epidemic_validation",
                    resources={"metadata": {"path": str(root / "notebooks")}},
                ).execute()
                nbformat.write(nb, path)
                exporter = HTMLExporter()

                exporter.embed_images = True
                body, _ = exporter.from_notebook_node(nb)
                body = body.replace('href="../figures/', 'href="../../figures/')
                (reports / (path.stem + ".html")).write_text(body)
                print(
                    f'  passed: {sum(c.cell_type=="code" for c in nb.cells)} code cells',
                    flush=True,
                )
        finally:
            if previous is None:
                os.environ.pop("JUPYTER_PATH", None)
            else:
                os.environ["JUPYTER_PATH"] = previous


def package_release(root=None):
    """Write a drop-in archive and file-hash manifest, excluding generated build clutter.

    root is the project directory, or is found from the working directory.
    Writes a SHA-256 file manifest in reports and a ZIP next to root.
    Virtual environments, build products, and caches are excluded.
    Returns the archive path.
    """
    from .plotting import repository_root

    root = Path(root) if root else repository_root()
    archive = root.parent / "epidemic_modeling_update.zip"
    omitted = {
        ".venv",
        ".git",
        "__pycache__",
        ".pytest_cache",
        ".ipynb_checkpoints",
        "build",
        "dist",
    }

    files = [
        p
        for p in root.rglob("*")
        if p.is_file()
        and not any(
            part in omitted or part.endswith(".egg-info")
            for part in p.relative_to(root).parts
        )
        and p.name not in {"release_manifest.json", ".DS_Store"}
        and p.suffix not in {".pyc", ".zip"}
    ]
    manifest = {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(files)
    }
    (root / "reports/release_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    with zipfile.ZipFile(
        archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6
    ) as z:
        for p in files:
            z.write(p, p.relative_to(root))

        z.write(root / "reports/release_manifest.json", "reports/release_manifest.json")

    print(
        f"Created {archive} ({archive.stat().st_size/1e6:.1f} MB; {len(files)+1} files)"
    )

    return archive
