function s_k = state_hard_margins(s_k, params)
    % STATE_HARD_MARGINS Clip population fractions to [0,1] and contact rate to its bounds.
    %
    % Returns a six-element state column; costates remain unchanged.
    %
    % u is an intervention column; s_k contains [s; i; alpha] followed by
    % three costates. params supplies model fields from default_si_params,
    % plus epsilon, w and sigma. Noise inputs use the model state/observation
    % units. These callbacks keep the standalone Coder model conventions.
    %
    % Author: Reza Sameni | Emory University

    s_k(1) = min(1.0, max(0, s_k(1)));
    s_k(2) = min(1.0, max(0, s_k(2)));
    s_k(3) = min(params.alpha_max, max(params.alpha_min, s_k(3)));
end
