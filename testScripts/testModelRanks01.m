% SEIRP Model rank tests:
% Reza Sameni
% April 2020

clear;
close all;
clc

alpha_e = 0.65;
alpha_i = 0.005;
kappa = 0.05;
rho = 0.08;
beta = 0.1;
mu = 0.02;
gamma_e = 0.001;
gamma_i = 0.001;
e = rand;
i = rand;
re = rand;
ri = rand;
p = rand;
s = 1 - e - i - re - ri - p;


% Find the Jacobian of the system
eta = alpha_e * e + alpha_i * i;
AA = zeros(5);
AA(1, 1) = alpha_e * s - (eta + kappa + rho);
AA(1, 2) = alpha_i * s - eta;
AA(1, 3) = -eta;
AA(1, 4) = -eta;
AA(1, 5) = -eta;
AA(2, 1) = kappa;
AA(2, 2) = -(mu + beta);
AA(3, 1) = rho;
AA(3, 3) = -gamma_e;
AA(4, 2) = beta;
AA(4, 4) = -gamma_i;
AA(5, 2) = mu;

% The observation matrix
CC = [0 1 0 0 0 ; 0 0 0 1 0 ; 0 0 0 0 1];

% check model observability
OB = obsv(AA, CC);
rnk = rank(OB)
