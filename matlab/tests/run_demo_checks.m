function run_demo_checks()
    % RUN_DEMO_CHECKS Run all MATLAB examples without interactive figure windows.
    % Save rendered example figures to reports/matlab_figures for visual review.
    % Author: Reza Sameni | Emory University
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root, 'matlab'));
    setup_paths();
    previous = get(groot, 'defaultFigureVisible');
    set(groot, 'defaultFigureVisible', 'off');

    cleanup = onCleanup(@() set(groot, 'defaultFigureVisible', previous));
    folder = fullfile(root, 'reports', 'matlab_figures');

    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    examples = {@demo_compartment_models, @demo_growth_estimation, @demo_kalman_estimation, ...
        @demo_intervention_control, @demo_historical_data, @demo_spatial_models};

    for k = 1:numel(examples)
        close all;
        examples{k}();
        figures = findall(groot, 'Type', 'figure');

        for j = 1:numel(figures)
            exportgraphics(figures(j), fullfile(folder, ...
                sprintf('example_%02d_figure_%02d.png', k, j)), 'Resolution', 120);
        end

        fprintf('Example %d passed (%d figures).\n', k, numel(figures));
    end

    close all;
end
