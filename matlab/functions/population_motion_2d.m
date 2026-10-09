function result = population_motion_2d(positions, velocities, dt, steps, box_size)
    % POPULATION_MOTION_2D Move agents with reflecting square boundaries.
    % result = population_motion_2d(positions,velocities,dt,steps,box_size)
    % positions/velocities: agents-by-2; initial positions lie in [0,box_size].
    % dt: positive time step; steps: nonnegative integer; box_size: positive
    % square side (default 1). Output is agents-by-2-by-(steps+1). Triangular-wave
    % reflection handles multiple wall crossings in one step.
    % Author: Reza Sameni | Emory University

    if nargin < 5
        box_size = 1;
    end

    assert(size(positions, 2) == 2 && isequal(size(positions), size(velocities)) && ...
        all(positions >= 0, 'all') && all(positions <= box_size, 'all') && dt > 0 && box_size > ...
        0 && steps >= 0 && steps == fix(steps), 'Invalid agent motion inputs.');
    phase = positions + velocities .* reshape(dt * (0:steps), 1, 1, []);
    result = box_size - abs(mod(phase, 2 * box_size) - box_size);
end
