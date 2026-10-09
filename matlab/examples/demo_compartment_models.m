function demo_compartment_models()
    % DEMO_COMPARTMENT_MODELS Example 01: SEIRP, capacity saturation, controlled SI.
    % Run setup_paths first. All rates are per day and states are population fractions.
    % Author: Reza Sameni | Emory University
    % Reference: Sameni (2020), arXiv:2003.11371.
    dt = .1;
    duration = 180;
    initial = {.998, .001, .001, 0, 0};
    [s, e, i, r, p] = seirp(.4, .3, .12, .04, .09, .002, 0, initial{:}, duration, dt);
    time = (0:numel(s) - 1) * dt;

    states = [s; e; i; r; p];
    figure('Color', 'w');
    plot(time, states', 'LineWidth', 2);
    grid on;
    xlabel('Time (days)');

    ylabel('Population fraction');
    title('SEIRP dynamics');
    legend('Susceptible', 'Exposed', 'Infected', 'Recovered', 'Passed', 'Location', 'best');

    assert(max(abs(sum(states, 1) - 1)) < 1e-12);
    [~, ~, i_sat, ~, p_sat] = seirp_saturated_resource(.4, .3, .12, .04, 0, initial{:}, ...
        duration, dt, .09, .03, .002, .015, .01, .04);
    figure('Color', 'w');
    tiledlayout(1, 2);
    nexttile;

    plot(time, [i; i_sat]', 'LineWidth', 2);
    title('Active infections');
    xlabel('Day');
    ylabel('Fraction');
    grid on;

    legend('Adequate', 'Saturated');
    nexttile;
    plot(time, [p; p_sat]', 'LineWidth', 2);
    title('Passed compartment');
    xlabel('Day');

    ylabel('Fraction');
    grid on;
    legend('Adequate', 'Saturated');
    contact = [.35 * ones(1, 300) .08 * ones(1, 900)];
    [~, infected] = si_controlled(contact, .1, .999, .001, numel(contact), dt);

    figure('Color', 'w');
    plot((0:numel(contact) - 1) * dt, infected, 'LineWidth', 2);
    xline(30, '--');
    grid on;
    xlabel('Day');

    ylabel('Infected fraction');
    title('Contact-rate reduction');
end
