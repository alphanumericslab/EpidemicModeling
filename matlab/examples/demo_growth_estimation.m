function demo_growth_estimation()
    % DEMO_GROWTH_ESTIMATION Example 02: causal and centered exponential fits.
    % Uses deterministic sinusoidal reporting variability for repeatable examples.
    % Author: Reza Sameni | Emory University
    % Reference: Sameni (2020), arXiv:2003.11371.
    day = 0:89;
    exact = 25 * exp(.035 * day);
    cases = exact .* exp(.08 * sin(day));
    [~, ~, growth, fit] = rt_exp_fit_log_lin_reg(cases, 14, 1);
    [~, ~, centered] = rt_exp_fit_log_lin_reg(cases, 14, 1, false);

    [~, ~, ~, ratios] = rt_exp_fit_gen_ratios(cases, 7, 5, 1);
    figure('Color', 'w');
    tiledlayout(1, 2);
    nexttile;
    plot(day(14:end), [cases(14:end); fit(14:end)]', 'LineWidth', 2);

    grid on;
    xlabel('Day');
    ylabel('Cases/day');
    title('Cases and one-step forecast');
    legend('Observed', 'Log regression');

    nexttile;
    plot(day(15:end - 7), [growth(15:end - 7); centered(15:end - 7); ratios(15:end - 7)]', ...
        'LineWidth', 2);
    yline(.035, '--');
    grid on;
    xlabel('Day');

    ylabel('Growth/day');
    title('Information timing');
    legend('Causal', 'Centered', 'Ratio');
    [~, ~, known] = rt_exp_fit_log_lin_reg(exact, 14, 1);

    assert(max(abs(known(14:end) - .035)) < 1e-12);

    if exist('nlinfit', 'file') == 2
        [~, ~, nonlinear] = rt_exp_fit_nonlin_ls(cases, 14, 1);
        fprintf('Median nonlinear growth: %.4f\n', median(nonlinear(14:end)));
    end
end
