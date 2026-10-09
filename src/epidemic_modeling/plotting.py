"""Plot and notebook presentation helpers; numerical models live in other modules."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def configure_plots():
    """Apply consistent, readable plotting defaults."""
    plt.rcParams.update({'figure.figsize':(10,4.8),'figure.dpi':120,'axes.spines.top':False,
        'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,'font.size':11,
        'axes.titlesize':15,'axes.titleweight':'bold','axes.labelsize':11,
        'lines.linewidth':2.3,'figure.constrained_layout.use':True,
        'axes.prop_cycle':plt.cycler(color=['#1565a8','#df8a29','#278574','#b14362','#7356a2'])})


def repository_root(start=None):
    """Find the project root from a notebook, repo root, or installed checkout."""
    start = Path(start or Path.cwd()).resolve()
    for path in (start,*start.parents):
        if (path/'pyproject.toml').is_file() and (path/'notebooks').is_dir(): return path
    raise FileNotFoundError("Run from the repository root or its notebooks folder")


def plot_compartments(time, states, title='SEIRP dynamics'):
    """Plot five compartment trajectories with consistent labels and colors."""
    fig,ax = plt.subplots()
    for values,label in zip(states,['Susceptible','Exposed','Infected','Recovered','Passed']): ax.plot(time,values,label=label)
    ax.set(xlabel='Time (days)',ylabel='Population fraction',title=title); ax.legend(ncol=3)
    return fig,ax


def plot_growth(time, cases, fits):
    """Compare observed case counts with named growth-model forecasts."""
    fig,axes = plt.subplots(1,2,figsize=(12,4.5))
    axes[0].plot(time,cases,color='#263447',label='Observed cases')
    for name,result in fits.items():
        axes[0].plot(time,result[3],label=name); axes[1].plot(time,result[2],label=name)
    axes[0].set(xlabel='Day',ylabel='Cases / day',title='Observed cases and one-step fits')
    axes[1].set(xlabel='Day',ylabel='Growth (per day)',title='Estimated growth'); axes[1].axhline(0,color='black',lw=1)
    for ax in axes: ax.legend()
    return fig,axes
