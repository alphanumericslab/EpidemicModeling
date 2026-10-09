function [cases, states] = forecast_npi(model, inputs)
    % FORECAST_NPI Forecast pre-update daily incidence under supplied interventions.
    % [cases,states] = forecast_npi(model,inputs)
    % model: fit_npi_model output; inputs: interventions-by-N within bounds.
    % cases: 1-by-N daily counts; states: 3-by-N post-update SI-alpha states.
    % Incidence uses the pre-update state to avoid a one-day forecast shift.
    % Author: Reza Sameni | Emory University
    p = model.params;
    initial = model.state(:);
    count = size(inputs, 2);

    assert(size(inputs, 1) == numel(p.a) && all(isfinite(inputs(:))) && all(inputs >= ...
        p.u_min(:), 'all') && all(inputs <= p.u_max(:), 'all'), 'Invalid forecast controls.');
    [s, i, alpha] = si_alpha_controlled(inputs, initial(1), initial(2), initial(3), ...
        p.u_max(:), p.alpha_min, p.alpha_max, p.gamma, p.a(:), p.b, p.beta, 0, 0, 0, count, ...
        p.dt, zeros(3, count));
    states = [s; i; alpha];
    before = [initial states(:, 1:end - 1)];
    cases = model.population * prod(before, 1);
end
