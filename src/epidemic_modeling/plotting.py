"""Plot and notebook presentation helpers; numerical models live in other modules."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def configure_plots():
    """Apply consistent, readable plotting defaults.

    Updates Matplotlib defaults for plots created afterwards. Call once
    before creating example figures. Returns None.
    """
    plt.rcParams.update(
        {
            "figure.figsize": (10, 4.8),
            "figure.dpi": 120,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.18,
            "font.size": 11,
            "axes.titlesize": 15,
            "axes.titleweight": "bold",
            "axes.labelsize": 11,
            "lines.linewidth": 2.3,
            "figure.constrained_layout.use": True,
            "axes.prop_cycle": plt.cycler(
                color=["#1565a8", "#df8a29", "#278574", "#b14362", "#7356a2"]
            ),
        }
    )


def repository_root(start=None):
    """Find the project root from a notebook, repo root, or installed checkout.

    start is a directory path, defaulting to the working directory.
    Searches parent directories for pyproject.toml and notebooks. Returns
    a pathlib.Path; raises FileNotFoundError if no project root is found.
    """
    start = Path(start or Path.cwd()).resolve()

    for path in (start, *start.parents):
        if (path / "pyproject.toml").is_file() and (path / "notebooks").is_dir():
            return path

    raise FileNotFoundError("Run from the repository root or its notebooks folder")


def plot_compartments(time, states, title="SEIRP dynamics"):
    """Plot five compartment trajectories with consistent labels and colors.

    time is a length-N vector in days. states contains five length-N
    series in S, E, I, R, P order. Returns the Matplotlib figure and axes
    so callers can adjust labels or save the plot.
    """
    fig, ax = plt.subplots()

    for values, label in zip(
        states, ["Susceptible", "Exposed", "Infected", "Recovered", "Passed"]
    ):
        ax.plot(time, values, label=label)

    ax.set(xlabel="Time (days)", ylabel="Population fraction", title=title)
    ax.legend(ncol=3)

    return fig, ax


def plot_growth(time, cases, fits):
    """Compare observed case counts with named growth-model forecasts.

    time and cases are length-N vectors. fits maps labels to the four
    arrays returned by the rolling regression functions. Returns a figure
    and two axes for case forecasts and estimated growth.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(time, cases, color="#263447", label="Observed cases")

    for name, result in fits.items():
        axes[0].plot(time, result[3], label=name)
        axes[1].plot(time, result[2], label=name)

    axes[0].set(
        xlabel="Day", ylabel="Cases / day", title="Observed cases and one-step fits"
    )
    axes[1].set(xlabel="Day", ylabel="Growth (per day)", title="Estimated growth")
    axes[1].axhline(0, color="black", lw=1)

    for ax in axes:
        ax.legend()

    return fig, axes
