function [s, i] = SI_Controlled(alpha, beta, s0, i0, K, dt)
%
% A numerical solver of the nonlinear susceptible-infected (SI) model with
% controlled infection rate and time dependent alpha
%
% Reza Sameni, January 2021
% reza.sameni@gmail.com
%
% The Open Source Electrophysiological Toolbox, version 3.14, January 2021
% Released under the GNU General Public License

s = zeros(1, K);
i = zeros(1, K);

s(1) = s0;
i(1) = i0;

% State equations
for t = 1 : K - 1
    s(t + 1) = max(0.0, min(1.0, s(t) - dt * alpha(t) * s(t) * i(t)));
    i(t + 1) = max(0.0, min(1.0, i(t) + dt * (alpha(t) * s(t) * i(t) - beta * i(t))));
end
