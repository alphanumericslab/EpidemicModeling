"""Deterministic spatial models with matching MATLAB implementations.

Author: Reza Sameni, Emory University. These refactor the historical spatial
scripts into reusable functions; they are educational diffusion/motion models.
"""
import numpy as np


def diffusion_2d(initial, diffusion, dt, spacing, steps):
    """Solve 2-D diffusion by explicit five-point finite differences.

    initial is a finite (rows,columns) concentration grid. Periodic boundaries
    conserve total mass. Returns (rows,columns,steps+1), including the initial
    grid. The stability condition diffusion*dt/spacing**2 <= 1/4 is enforced.
    """
    x = np.asarray(initial,float)
    ratio = diffusion*dt/spacing**2
    if x.ndim!=2 or not np.all(np.isfinite(x)) or not 0<=ratio<=.25 or dt<=0 or spacing<=0 or steps<0 or int(steps)!=steps:
        raise ValueError("finite 2-D grid, nonnegative integer steps, and stable positive time/space grid required")
    result = np.empty(x.shape+(steps+1,)); result[:,:,0] = x
    for k in range(steps):
        lap = np.roll(x,1,0)+np.roll(x,-1,0)+np.roll(x,1,1)+np.roll(x,-1,1)-4*x
        x = x+ratio*lap; result[:,:,k+1] = x
    return result


def population_motion_2d(positions, velocities, dt, steps, box_size=1.):
    """Move agents with reflecting square boundaries using a triangular-wave map.

    positions and velocities have shape (agents,2); initial positions lie
    inside [0,box_size]. Returns positions (agents,2,steps+1). Reflection
    handles arbitrarily many wall crossings in a single step.
    """
    p = np.asarray(positions,float); v = np.asarray(velocities,float)
    if p.ndim!=2 or p.shape[1]!=2 or p.shape!=v.shape or np.any(p<0) or np.any(p>box_size) or dt<=0 or box_size<=0 or steps<0 or int(steps)!=steps:
        raise ValueError("matching agent-by-2 arrays, bounded positions, positive dt/box, and integer steps required")
    phase = p[:,:,None]+v[:,:,None]*(dt*np.arange(steps+1))
    phase = np.mod(phase,2*box_size)
    return box_size-np.abs(phase-box_size)
