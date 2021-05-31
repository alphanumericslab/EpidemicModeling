function [s, i, alpha] = SIalpha_Controlled(u, s0, i0, alpha0, u_max, alpha_min, alpha_max, gamma, a, b, beta, s_noise_std, i_noise_std, alpha_noise_std, K, dt)
%
% A numerical solver of the nonlinear susceptible-infected (SI) alpha model with
% controlled infection rate
%
% Reza Sameni, January 2021
% reza.sameni@gmail.com
%
% Modified:
% May 2021: initial condition removed from output sequences
% 
% The Open Source Electrophysiological Toolbox, version 3.14, January 2021
% Released under the GNU General Public License

s = zeros(1, K + 1);
i = zeros(1, K + 1);
alpha = zeros(1, K + 1);

s(1) = s0;
i(1) = i0;
alpha(1) = alpha0;

% State equations
for t = 1 : K
    s(t + 1) = max(0.0, min(1.0, s(t) - dt * (alpha(t) * s(t) * i(t) + randn * s_noise_std)));
    i(t + 1) = max(0.0, min(1.0, i(t) + dt * (alpha(t) * s(t) * i(t) - beta * i(t) + randn * i_noise_std)));
    alpha(t + 1) = max(alpha_min, min(alpha_max, alpha(t) + dt * (-gamma * alpha(t) + gamma * b + gamma * a'*(u_max - u(:, t)) + randn * alpha_noise_std)));
end

s = s(2 : end);
i = i(2 : end);
alpha = alpha(2 : end);
