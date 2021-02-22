% A simulation of multiple scenarios of the susceptible-exposed-infected-recovered-passed (SEIRP) model for
% epidemic diseases
%
% Reza Sameni, March 2020
% reza.sameni@gmail.com
%
% The Open Source Electrophysiological Toolbox, version 3.14, March 2020
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/

clear;
close all;
clc;

dt = 0.1; % simulation time step (in days)
scenario = 'A';

switch scenario
    case 'A' % Immunizing disease
        T = 50; % days
        K = round(T/dt);
        alpha_e = 0.65*ones(1, K);
        alpha_i = 0.005*ones(1, K);
        kappa = 0.05*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = 0.0*ones(1, K);
        N = 84.0e6;
    case 'B' % Non-immunizing disease
        T = 4000; % days
        K = round(T/dt);
        alpha_e = 0.65*ones(1, K);
        alpha_i = 0.005*ones(1, K);
        kappa = 0.05*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = 0.001*ones(1, K);
        N = 84.0e6;
    case 'C'
        T = 120; % days
        K = round(T/dt);
        alpha_e = 0.65*linspace(1, 0.01, K);
        alpha_i = 0.005*linspace(1, 0.01, K);
        kappa = 0.05*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = 0.001*ones(1, K);
        N = 84.0e6;
    case 'D'
        T = 4000; % days
        K = round(T/dt);
        alpha_e = 0.65*ones(1, K);
        alpha_i = 0.005*ones(1, K);
        kappa = 0.005*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = 0.001*ones(1, K);
        N = 84.0e6;
    case 'E'
        T = 4000; % days
        K = round(T/dt);
        alpha_e = 0.65*ones(1, K);
        alpha_i = 0.005*ones(1, K);
        kappa = 0.05*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = (1/365)*ones(1, K);
        N = 84.0e6;
end

e0 = 1/N;
s0 = 1 - e0;
[s, e, i, r, p] = SEIRP(alpha_e, alpha_i, kappa, rho, beta, mu, gamma, s0, e0, 0, 0, 0, T, dt);

% Find eigenvalues of the linearized system (for s(t) close to 1)
check_point = 30;
A = [alpha_e(check_point)-kappa(check_point)-rho(check_point), alpha_i(check_point), 0, 0 ; kappa(check_point), -beta(check_point)-mu(check_point), 0, 0 ; rho(check_point), beta(check_point), -gamma(check_point), 0 ; 0, mu(check_point), 0, 0];
C = [zeros(3, 1) eye(3)];
% check model observability
OB = obsv(A, C);
rnk = rank(OB)

% Find eigenvalues of the Jacobian of the system (for arbitrary s(t))
check_point = 30;
AA = zeros(4);
ss = 1 - e(check_point) - i(check_point) - r(check_point) - p(check_point);
AA(1,1) = alpha_e(check_point)*(ss - e(check_point)) - alpha_i(check_point)*i(check_point) - kappa(check_point) - rho(check_point);
AA(1,2) = alpha_i(check_point)*(ss - i(check_point)) - alpha_e(check_point)*e(check_point);
AA(1,3) = - alpha_e(check_point)*e(check_point) - alpha_i(check_point)*i(check_point);
AA(1,4) = - alpha_e(check_point)*e(check_point) - alpha_i(check_point)*i(check_point);
AA(2, 1) = kappa(check_point);
AA(2,2) = -beta(check_point) - mu(check_point);
AA(3,1) = rho(check_point);
AA(3,2) = beta(check_point);
AA(3,3) = -gamma(check_point);
AA(4,2) = mu(check_point);
OB2 = obsv(AA, C);
rnk2 = rank(OB2)


[V D] = eig(A)
lambda1 = 0;
lambda2 = -gamma(check_point);
delta = alpha_e(check_point) - kappa(check_point) - rho(check_point);
lambda3 = (delta - beta(check_point) - mu(check_point) + sqrt((beta(check_point) + mu(check_point) + delta)^2 + 4*kappa(check_point)*alpha_i(check_point))) / 2;
lambda4 = (delta - beta(check_point) - mu(check_point) - sqrt((beta(check_point) + mu(check_point) + delta)^2 + 4*kappa(check_point)*alpha_i(check_point))) / 2;
lambda = [lambda1, lambda2, lambda3, lambda4]
v1 = [0, 0, 0, 1]';
v2 = [0, 0, 1, 0]';
v3 = [1, (lambda3 - delta)/alpha_i(check_point), (rho(check_point)*alpha_i(check_point)+beta(check_point)*(lambda3-delta))/alpha_i(check_point)/(lambda3 + gamma(check_point)), mu(check_point)*(lambda3 - delta)/lambda3/alpha_i(check_point)]';
v3 = v3/sqrt(v3'*v3);
v4 = [1, (lambda4 - delta)/alpha_i(check_point), (rho(check_point)*alpha_i(check_point)+beta(check_point)*(lambda4-delta))/alpha_i(check_point)/(lambda4 + gamma(check_point)), mu(check_point)*(lambda4 - delta)/lambda4/alpha_i(check_point)]';
v4 = v4/sqrt(v4'*v4);

t = dt*(0 : K - 1);
ii = (e0/alpha_i(check_point))*(lambda3 - delta)*(lambda4 - delta)/(lambda3 - lambda4)*(exp(lambda4*t) - exp(lambda3*t));
ee = e0/(lambda3 - lambda4)*((lambda3 - delta)*exp(lambda4*t) + (delta - lambda4)*exp(lambda3*t));

figure;
hold on
plot(t, s, 'b', 'linewidth', 3);
plot(t, e, 'c', 'linewidth', 3);
plot(t, i, 'r', 'linewidth', 3);
plot(t, r, 'g', 'linewidth', 3);
plot(t, p, 'k', 'linewidth', 3);
grid
legend('S(t)', 'E(t)', 'I(t)', 'R(t)', 'P(t)');
ylabel('Population Ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');

figure;
hold on
plot(t, i*N, 'linewidth', 3);
plot(t, e*N, 'linewidth', 3);
plot(t, ii*N, 'linewidth', 1);
plot(t, ee*N, 'linewidth', 1);
grid
legend('I(t)', 'E(t)', 'II(t)', 'EE(t)');
% legend('I(t)', 'E(t)');
ylabel('Population ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
% title('The infection rate over time');

figure;
hold on
plot(e, i, 'linewidth', 1);
grid
ylabel('I');
xlabel('E');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
title('Phase plot');

tfinal = find(t >= 40, 1);
figure
hold on
% plot(t, e, 'linewidth', 3);
% plot(t, i, 'linewidth', 3);
plot(t(1 : tfinal), i(1 : tfinal)./e(1 : tfinal), 'linewidth', 3);
plot(t(1 : tfinal), ii(1 : tfinal)./ee(1 : tfinal), 'linewidth', 3);
set(gca, 'box', 'on');
set(gca, 'fontsize', 14)
title('I(t)/E(t)')
grid

ItoEratio_at_check_point = (lambda3 - delta)/alpha_i(check_point)
