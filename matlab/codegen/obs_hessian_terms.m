function [gs, Gsp, gv, Gvp] = obs_hessian_terms(u, s_k, Pk, v_bar, Rk, params)
    % OBS_HESSIAN_TERMS Return zero second-order observation/noise corrections.
    %
    % All four outputs are scalar zeros for this model. Covariance arguments
    % are accepted for the generic callback interface.
    %
    % u is an intervention column; s_k contains [s; i; alpha] followed by
    % three costates. params supplies model fields from default_si_params,
    % plus epsilon, w and sigma. Noise inputs use the model state/observation
    % units. These callbacks keep the standalone Coder model conventions.
    %
    % Author: Reza Sameni | Emory University

    gs = zeros(1);
    Gsp = zeros(1);

    gv = zeros(1);
    Gvp = zeros(1);

end
