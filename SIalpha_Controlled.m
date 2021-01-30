function [s, i, alpha] = SIalpha_Controlled(u, u_max, alpha_min, alpha_max, gamma, alpha0, beta, s0, i0, K, dt)
%
% A numerical solver of the nonlinear susceptible-infected (SI) alpha model with
% controlled infection rate
%
% Reza Sameni, January 2021
% reza.sameni@gmail.com
%
% The Open Source Electrophysiological Toolbox, version 3.14, January 2021
% Released under the GNU General Public License

s = zeros(1, K);
i = zeros(1, K);
alpha = zeros(1, K);

s(1) = s0;
i(1) = i0;
alpha(1) = alpha0;

% State equations
for t = 1 : K - 1
    s(t + 1) = max(0.0, min(1.0, s(t) - dt * alpha(t) * s(t) * i(t)));
    i(t + 1) = max(0.0, min(1.0, i(t) + dt * (alpha(t) * s(t) * i(t) - beta * i(t))));
    alpha(t + 1) = max(alpha_min, min(alpha_max, alpha(t) + dt * (-gamma * alpha(t) + gamma * b + gamma * a'*(u_max - u))));
end
