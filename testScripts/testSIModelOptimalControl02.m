% A script for testing optimal control of an SI model with
% non-pharmaceutical intervention plans as control inputs
%
% Corresponds to TWO-variable state equations in notes

close all
clear
clc

dt = 1; % Time unit
T = 300; % Time length
I0 = 1000; % initial infected seed
N = 100000; % total population
L = 12;
K = round(T/dt); % number of samples

u = zeros(L, K); % control inputs
% control input weights
ww = rand(L, 1); ww = L*ww/sum(ww);
w = diag(ww) * ones(L, K);
A = 0.01 * rand(L, 1);
u_min = zeros(L, 1);
u_max = round(4.0 * rand(L, 1));

alpha_min = 0.0;
alpha_max = 1.0;

s = zeros(1, K); % susceptibles
i = zeros(1, K); % infected
alpha = zeros(1, K); % infection rate
rho = zeros(1, K); % combined S costate (shadow state)
lambda2 = zeros(1, K); % I costate (shadow state)
H = zeros(1, K); % Hamiltonian
J0 = zeros(1, K); % human factor cost
J1 = zeros(1, K); % economic factor cost
J = zeros(1, K); % total cost

gamma = 0    *    1e-7;
beta = 1/7.0;
s(1) = (N-I0)/N;
i(1) = I0/N;
rho(1) = -1.0;
lambda2(1) = 0.0;

for t = 1 : K
    % Optimal control
    u(: , t) = round(max(u_min, min(u_max, gamma * w(:, t) ./ (-2 * A * s(t) * i(t) * rho(t)))));
    alpha(t) =  sum(A .* (u_max .^ 2 - u(:, t) .^ 2));
%     sum(A .* (u_max .^ 2 - u(:, t) .^ 2))
    sum(u(: , t))
%     alpha(t)
    
    % Hamiltonian
    H(t) = -rho(t) * alpha(t) * s(t) * i(t) - beta * lambda2(t) * i(t) + gamma * w(:, t)' * u(: , t);
    
    % Costs
    J0(t) = alpha(t) * s(t) * i(t);
    J1(t) = w(:, t)' * u(: , t);
    J(t) = J0(t) + gamma * J1(t);
    
    if(t < K)
        % Costate equations
        rho(t + 1) = rho(t) + dt * rho(t) * alpha(t) * (i(t) - s(t)) - dt * beta * lambda2(t);
        lambda2(t + 1) = lambda2(t) + dt * rho(t) * alpha(t) * s(t) + dt * beta * lambda2(t);
        
        % State equations
        s(t + 1) = max(0.0, min(1.0, s(t) - dt * alpha(t) * s(t) * i(t)));
        i(t + 1) = max(0.0, min(1.0, i(t) + dt * alpha(t) * s(t) * i(t) - dt * beta * i(t)));
    end
end

figure
subplot(511);
plot(s)
hold on
plot(i)
grid
legend('s', 'i');

subplot(512);
plot(alpha)
grid
legend('alpha');

subplot(513);
plot(H)
grid
legend('Hamiltonian');

subplot(514);
plot(J0)
hold on
plot(gamma * J1)
plot(J)
grid
legend('J0 (Human)', 'gamma * J1 (Economic)', 'J (Total)');

subplot(515);
plot(u');
grid
legend('inputs');
