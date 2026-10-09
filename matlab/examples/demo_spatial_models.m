function demo_spatial_models()
    % DEMO_SPATIAL_MODELS Example 06: periodic diffusion and reflected motion.
    % Author: Reza Sameni | Emory University
    initial = zeros(51);
    initial(26, 26) = 1;
    result = diffusion_2d(initial, 1, .2, 1, 100);
    figure('Color', 'w');
    tiledlayout(1, 3);

    for step = [0 20 100]
        nexttile;
        imagesc(result(:, :, step + 1));
        axis image;
        colorbar;

        title(sprintf('Step %d', step));
        xlabel('Column');
        ylabel('Row');
    end

    assert(max(abs(squeeze(sum(sum(result, 1), 2)) - 1)) < 1e-12);
    positions = [.2 .3; .7 .8; .4 .6];
    velocities = [.14 .09; -.11 .07; .08 -.13];
    paths = population_motion_2d(positions, velocities, .1, 180);
    figure('Color', 'w');

    hold on;

    for agent = 1:size(positions, 1)
        plot(squeeze(paths(agent, 1, :)), squeeze(paths(agent, 2, :)), 'LineWidth', 2);
    end

    axis equal;
    xlim([0 1]);
    ylim([0 1]);
    grid on;
    title('Reflecting agents');

    xlabel('x');
    ylabel('y');
end
