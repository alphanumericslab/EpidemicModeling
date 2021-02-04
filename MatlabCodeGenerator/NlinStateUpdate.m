% Nonlinear state update
function [u, s_k_plus_one] = NlinStateUpdate(u, s_k, w_bar, params)

% if(isnan(u)) % Optimal control
%     u = zeros(length(params.u_max), 1);
%     phi = params.epsilon * params.w - params.gamma * s_k(6) * params.a;
%     u(phi >= 0) = params.u_min(phi >= 0);
%     u(phi < 0) = params.u_max(phi < 0);
% end

% Optimal control replaces nans (by design feature)
phi = params.epsilon * params.w - params.gamma * s_k(6) * params.a;
for kk = 1 : length(u)
    if(isnan(u(kk)))
        if(phi(kk) >= 0)
            u(kk) = params.u_min(kk);
        else
            u(kk) = params.u_max(kk);
        end
    end
end

rho = s_k(4) - s_k(5) - (1 - params.epsilon);

s_k_plus_one = zeros(6, 1);

% State equations
s_k_plus_one(1) = max(0.0, min(1.0, s_k(1) - params.dt * s_k(3) * s_k(1) * s_k(2)));
s_k_plus_one(2) = max(0.0, min(1.0, s_k(2) + params.dt * (s_k(3) * s_k(1) * s_k(2) - params.beta * s_k(2))));
s_k_plus_one(3) = max(params.alpha_min, min(params.alpha_max, s_k(3) + params.dt * (-params.gamma * s_k(3) + params.gamma * params.b + params.gamma * params.a'*(params.u_max - u))));

% Costate equations
s_k_plus_one(4) = s_k(4) + params.dt * rho * s_k(3) * s_k(2);
s_k_plus_one(5) = s_k(5) + params.dt * (rho * s_k(3) * s_k(1) + params.beta * s_k(5));
s_k_plus_one(6) = s_k(6) + params.dt * (rho * s_k(1) * s_k(2) + params.gamma * s_k(6));

end