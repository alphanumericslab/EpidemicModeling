% state-parameter augmented SEIRP Model rank tests:
% Reza Sameni
% May 2020

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
epsilon_alpha_e = 1e-4;
epsilon_alpha_i = 1e-4;
epsilon_kappa = 1e-4;
epsilon_rho = 1e-4;
epsilon_beta = 1e-4;
epsilon_mu = 1e-4;
epsilon_gamma_e = 1e-4;
epsilon_gamma_i = 1e-4;

e = rand;
i = rand;
re = rand;
ri = rand;
p = rand;
s = 1 - e - i - re - ri - p;


% Find the Jacobian of the system
eta = alpha_e * e + alpha_i * i;
AA = zeros(13);
AA(1, 1) = alpha_e * s - (eta + kappa + rho);
AA(1, 2) = alpha_i * s - eta;
AA(1, 3) = -eta;
AA(1, 4) = -eta;
AA(1, 5) = -eta;
AA(1, 6) = e * s;
AA(1, 7) = i * s;
AA(1, 8) = -e;
AA(1, 9) = -e;

AA(2, 1) = kappa;
AA(2, 2) = -(mu + beta);
AA(2, 8) = e;
AA(2, 10) = -i;
AA(2, 11) = -i;

AA(3, 1) = rho;
AA(3, 3) = -gamma_e;
AA(3, 9) = e;
AA(3, 12) = -re;

AA(4, 2) = beta;
AA(4, 4) = -gamma_i;
AA(4, 10) = i;
AA(4, 13) = -ri;

AA(5, 2) = mu;
AA(5, 11) = i;
AA(6, 6) = epsilon_alpha_e;
AA(7, 7) = epsilon_alpha_i;
AA(8, 8) = epsilon_kappa;
AA(9, 9) = epsilon_rho;
AA(10, 10) = epsilon_beta;
AA(11, 11) = epsilon_mu;
AA(12, 12) = epsilon_gamma_e;
AA(13, 13) = epsilon_gamma_i;

% The observation matrix
CC = [0 1 0 0 0 ; 0 0 0 1 0 ; 0 0 0 0 1];
CC = [CC zeros(3, 8)];

% BB = [zeros(5, 8); eye(8)];
BB = [zeros(5, 2); eye(2) ; zeros(6, 2)];
% controllablevars = [1:5 12];
% BB = [zeros(3, 2); eye(2)];
% AA = AA([2:4 6:7], [2:4 6:7])

CC0 = [eye(5) zeros(5, 8)];
M = OutputCTRB(AA, BB, CC0);
ct_rank0 = rank(M)



CT = ctrb(AA, BB);
ct_rank = rank(CT)

% check model observability
observabalevars = [1:5 12];
OB = obsv(AA(observabalevars, observabalevars), CC(:, observabalevars));
ob_rank = rank(OB)
