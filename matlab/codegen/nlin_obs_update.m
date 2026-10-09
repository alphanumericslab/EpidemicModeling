function x_k = nlin_obs_update(u, s_k, v_bar, params)
    % NLIN_OBS_UPDATE Return new-case incidence s*i*alpha plus the measurement-noise mean.
    %
    % Returns one scalar fraction; multiply by population to obtain counts.
    %
    % u is an intervention column; s_k contains [s; i; alpha] followed by
    % three costates. params supplies model fields from default_si_params,
    % plus epsilon, w and sigma. Noise inputs use the model state/observation
    % units. These callbacks keep the standalone Coder model conventions.
    %
    % Author: Reza Sameni | Emory University

    x_k = s_k(1) * s_k(2) * s_k(3) + v_bar;
end
