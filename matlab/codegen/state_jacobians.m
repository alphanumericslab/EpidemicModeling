function [A, B] = state_jacobians(u, s_k, w_bar, params)
    % STATE_JACOBIANS Return state-transition and process-noise Jacobians.
    %
    % A and B are 6-by-6 matrices. B is the identity noise map.
    %
    % u is an intervention column; s_k contains [s; i; alpha] followed by
    % three costates. params supplies model fields from default_si_params,
    % plus epsilon, w and sigma. Noise inputs use the model state/observation
    % units. These callbacks keep the standalone Coder model conventions.
    %
    % Author: Reza Sameni | Emory University

    A = zeros(6);
    A(1, 1) = 1 - params.dt * s_k(3) * s_k(2);
    A(1, 2) = -params.dt * s_k(3) * s_k(1);
    A(1, 3) = -params.dt * s_k(1) * s_k(2);

    A(2, 1) = params.dt * s_k(2) * s_k(3);
    A(2, 2) = 1 + params.dt * (s_k(1) * s_k(3) - params.beta);
    A(2, 3) = params.dt * s_k(1) * s_k(2);

    A(3, 3) = 1 - params.dt * params.gamma;
    phi = params.epsilon * params.w - params.gamma * s_k(6) * params.a;

    for kk = 1:length(u)

        if isnan(u(kk))

            if phi(kk) > -1.0 / params.sigma && phi(kk) < 1.0 / params.sigma
                A(3, 6) = A(3, 6) - params.gamma * params.dt * (params.sigma / ...
                    2) * params.a(kk) * (params.u_max(kk) - params.u_min(kk));
            end
        end
    end

    rho = s_k(4) - s_k(5) - (1 - params.epsilon);
    A(4, 2) = params.dt * s_k(3) * rho;
    A(4, 3) = params.dt * s_k(2) * rho;
    A(4, 4) = 1 + params.dt * s_k(2) * s_k(3);
    A(4, 5) = -params.dt * s_k(2) * s_k(3);

    A(5, 1) = params.dt * s_k(3) * rho;
    A(5, 3) = params.dt * s_k(1) * rho;
    A(5, 4) = params.dt * s_k(1) * s_k(3);
    A(5, 5) = 1 - params.dt * (s_k(1) * s_k(3) - params.beta);

    A(6, 1) = params.dt * s_k(2) * rho;
    A(6, 2) = params.dt * s_k(1) * rho;
    A(6, 4) = params.dt * s_k(1) * s_k(2);
    A(6, 5) = -params.dt * s_k(1) * s_k(2);
    A(6, 6) = 1 + params.dt * params.gamma;

    B = eye(6);
end
