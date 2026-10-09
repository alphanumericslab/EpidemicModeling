function [fs, Cs, fw, Cw] = state_hessian_terms(u, s_k, Pk, w_bar, Qk, params)
    % STATE_HESSIAN_TERMS Return the original zero second-order state/noise corrections.
    %
    % fs/fw are six-element columns; Cs/Cw are 6-by-6 matrices.
    % Covariance arguments are accepted for the generic callback interface.
    %
    % u is an intervention column; s_k contains [s; i; alpha] followed by
    % three costates. params supplies model fields from default_si_params,
    % plus epsilon, w and sigma. Noise inputs use the model state/observation
    % units. These callbacks keep the standalone Coder model conventions.
    %
    % Author: Reza Sameni | Emory University

    fs = zeros(6, 1);
    Cs = zeros(6);

    fw = zeros(6, 1);
    Cw = zeros(6);

end
