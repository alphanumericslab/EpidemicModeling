% A comperative simulation of the susceptible-exposed-infected-recovered-passed (SEIRP) model for
% epidemic diseases under parameter changes
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
T = 150; % days
K = round(T/dt);
N = 84.0e6;
E0 = 1;

% parameter set 1
alpha_e = 0.6*ones(1, K);
alpha_i = 0.005*ones(1, K);
kappa = 0.05*ones(1, K);
rho = 0.08*ones(1, K);
beta = 0.1*ones(1, K);
mu = 0.02*ones(1, K);
gamma = 0.001*ones(1, K);
[s1, e1, i1, r1, p1] = SEIRP(alpha_e, alpha_i, kappa, rho, beta, mu, gamma, (N-E0)/N, E0/N, 0, 0, 0, T, dt);

% parameter set 2
QStartDay = round(30/dt);
QEndDay = round(90/dt);
alpha_e = [0.6*ones(1, QStartDay), 0.1*ones(1, (QEndDay - QStartDay)), 0.4*ones(1, K - QEndDay)];
alpha_i = [0.005*ones(1, QStartDay), 0.001*ones(1, (QEndDay - QStartDay)), 0.001*ones(1, K - QEndDay)];
kappa = 0.05*ones(1, K);
rho = 0.08*ones(1, K);
beta = 0.1*ones(1, K);
mu = 0.02*ones(1, K);
gamma = 0.001*ones(1, K);
[s2, e2, i2, r2, p2] = SEIRP(alpha_e, alpha_i, kappa, rho, beta, mu, gamma, (N-E0)/N, E0/N, 0, 0, 0, T, dt);

t = dt*(0 : K - 1);
% figure;
% hold on
% plot(t, s1, 'b', 'linewidth', 3);
% plot(t, e1, 'c', 'linewidth', 3);
% plot(t, i1, 'r', 'linewidth', 3);
% plot(t, r1, 'g', 'linewidth', 3);
% plot(t, p1, 'k', 'linewidth', 3);
% grid
% legend('s(t)', 'e(t)', 'i(t)', 'r(t)', 'p(t)');
% ylabel('Population Ratio');
% xlabel('days');
% set(gca, 'fontsize', 16)
% set(gca, 'box', 'on');

figure;
subplot(211);
plot(t, i1, 'linewidth', 3);
hold on;
plot(t, e1, 'linewidth', 3);
plot(t, i2, 'linewidth', 3);
plot(t, e2, 'linewidth', 3);
grid
legend('i1(t)', 'e1(t)', 'i2(t)', 'e2(t)');
ylabel('Population ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
subplot(212);
plot(t, p1, 'linewidth', 3);
hold on;
plot(t, p2, 'linewidth', 3);
grid
legend('p1(t)','p2(t)');
ylabel('Population ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
% title('The infection rate over time');