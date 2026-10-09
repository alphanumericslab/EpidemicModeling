function result = diffusion_2d(initial, diffusion, dt, spacing, steps)
    % DIFFUSION_2D Solve periodic 2-D diffusion by explicit finite differences.
    % result = diffusion_2d(initial,diffusion,dt,spacing,steps)
    % initial: rows-by-columns finite grid; diffusion: nonnegative coefficient;
    % dt/spacing: positive time/space steps; steps: nonnegative integer.
    % Output is rows-by-columns-by-(steps+1), including the initial grid.
    % Stability requires diffusion*dt/spacing^2 <= 1/4. Total mass is conserved.
    % Author: Reza Sameni | Emory University
    ratio = diffusion * dt / spacing^2;

    assert(ismatrix(initial) && all(isfinite(initial(:))) && ratio >= 0 && ratio <= .25 && ...
        dt > 0 && spacing > 0 && steps >= 0 && steps == fix(steps), ...
        'Invalid diffusion grid or unstable step.');
    x = initial;
    result = zeros([size(x) steps + 1]);
    result(:, :, 1) = x;

    for k = 1:steps
        lap = circshift(x, [1 0]) + circshift(x, [-1 0]) + circshift(x, [0 1]) + ...
            circshift(x, [0 -1]) - 4 * x;
        x = x + ratio * lap;
        result(:, :, k + 1) = x;
    end
end
