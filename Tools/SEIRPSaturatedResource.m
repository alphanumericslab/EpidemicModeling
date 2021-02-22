function [s, e, i, r, p] = SEIRPSaturatedResource(alpha_e, alpha_i, kappa, rho, gamma, s0, e0, i0, r0, p0, T, dt, beta_0, beta_s, mu_0, mu_s, sigma, i_0)
%
% A numerical solver of the nonlinear susceptible-exposed-infected-recovered-passed (SEIRP) model for
% epidemic diseases under the saturated resources constraint
%
% Reza Sameni, March 2020
% reza.sameni@gmail.com
%
% The Open Source Electrophysiological Toolbox, version 3.14, March 2020
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/

K = round(T/dt);
s = zeros(1, K);
e = zeros(1, K);
i = zeros(1, K);
r = zeros(1, K);
p = zeros(1, K);

s(1) = s0;
e(1) = e0;
i(1) = i0;
r(1) = r0;
p(1) = p0;

for t = 1 : K - 1
    h = (tanh((i(t) - i_0)/sigma) + 1)/2;
    beta = (beta_s - beta_0)*h + beta_0;
    mu = (mu_s - mu_0)*h + mu_0;
    
    s(t + 1) = (-alpha_e(t) * s(t) * e(t) - alpha_i(t) * s(t) * i(t) + gamma(t) * r(t)) * dt + s(t);
    e(t + 1) = (alpha_e(t) * s(t) * e(t) + alpha_i(t) * s(t) * i(t) - kappa(t) * e(t) - rho(t) * e(t)) * dt + e(t);
    i(t + 1) = (kappa(t) * e(t) - beta * i(t) - mu * i(t))* dt + i(t);
    r(t + 1) = (beta * i(t) + rho(t) * e(t) - gamma(t) * r(t)) * dt + r(t);
    p(t + 1) = (mu * i(t)) * dt + p(t);
end


