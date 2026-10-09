function [C, D] = obs_jacobian(u, s_k, v_bar, params)
    % OBS_JACOBIAN Return incidence and measurement-noise Jacobians.
    %
    % C is 1-by-6; D is the scalar noise map 1.
    %
    % u is an intervention column; s_k contains [s; i; alpha] followed by
    % three costates. params supplies model fields from default_si_params,
    % plus epsilon, w and sigma. Noise inputs use the model state/observation
    % units. These callbacks keep the standalone Coder model conventions.
    %
    % Author: Reza Sameni | Emory University

    C = [s_k(2) * s_k(3), s_k(1) * s_k(3), s_k(1) * s_k(2), 0, 0, 0];
    D = 1;
end
