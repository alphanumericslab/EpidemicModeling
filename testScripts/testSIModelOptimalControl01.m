% A script for testing optimal control of an SI model with
% non-pharmaceutical intervention plans as control inputs
%
% Corresponds to linear causal input dependency (THREE-variable state
% equations in Overleaf notes)

close all
clear
clc

dt = 0.1; % Time unit
T = 90; % Time length
N = 100000; % total population
I0 = 100; % initial infected seed popoulation
L = 12;
K = round(T/dt); % number of samples
u = zeros(L, K); % control inputs
% control input weights
ww = ones(L, 1); ww = L*ww/sum(ww);
w = diag(ww) * ones(L, K);
a = 0.01 * rand(L, 1); % input influence weight vector
u_min = zeros(L, 1); % minimum input values
u_max = round(4.0 * rand(L, 1)); % maximum input values
alpha_min = 0.0; % minimum alpha
alpha_max = 1.0; % maximum alpha
% % % epsilon = 0.0001; % [0, 1]: 0 neglects NPI cost and 1 neglects human cost!
gamma = 1/5.0; % input to contact influence rate (inverse time)
beta = 1/7.0; % recovery rate (inverse time)

max_instances_per_eps = 10;

figure
hold on
for multiobjective_ind = 1 : 100
    epsilon = rand;
    for instance = 1 : max_instances_per_eps
% % %         if(instance == max_instances_per_eps)
% % %             optimal_input = true;
% % %         else
% % %             optimal_input = false;
% % %         end
        
        s = zeros(1, K); % susceptibles
        i = zeros(1, K); % infected
        alpha = zeros(1, K); % infection rate
        rho = zeros(1, K); % auxiliary costate variable
        lambda1 = zeros(1, K); % S costate
        lambda2 = zeros(1, K); % I costate
        lambda3 = zeros(1, K); % alpha costate
        H = zeros(1, K); % Hamiltonian
        J0 = zeros(1, K); % human factor cost
        J1 = zeros(1, K); % economic factor cost
        J = zeros(1, K); % total cost
        s(1) = (N - I0) / N;
        i(1) = I0 / N;
        alpha(1) = 0.1;
        lambda1(1) = -1;
        lambda2(1) = 1;
        lambda3(1) = 1;
        
% % %         if(optimal_input) % optimal input
        if(instance == max_instances_per_eps) % optimal input
            % % %             for itr = 1 : 10 % number of iterations to fix the end points
            % time evolution
            for t = 1 : K
                % auxiliary variable
                rho(t) = lambda1(t) - lambda2(t) - (1 - epsilon);
                
                % Optimal control
                phi = epsilon * w(:, t) - gamma * lambda3(t) * a;
                u(phi >= 0 , t) = u_min(phi >= 0);
                u(phi < 0 , t) = u_max(phi < 0);
                
                % Hamiltonian
                H(t) =    +    (-rho(t) * alpha(t) * s(t) * i(t) - beta * lambda2(t) * i(t) + epsilon * w(:, t)' * u(: , t) + lambda3(t) * (-gamma * alpha(t) + gamma * a' * (u_max - u(:, t))));
                
                % Costs
                J0(t) = alpha(t) * s(t) * i(t);
                J1(t) = w(:, t)' * u(: , t);
                J(t) = (1 - epsilon) * J0(t) + epsilon * J1(t);
                
                if(t < K)
                    % Costate equations
                    lambda1(t + 1) = lambda1(t)       +          dt * rho(t) * alpha(t) * i(t);
                    lambda2(t + 1) = lambda2(t)       +          dt * (rho(t) * alpha(t) * s(t) + beta * lambda2(t));
                    lambda3(t + 1) = lambda3(t)       +          dt * (rho(t) * s(t) * i(t) + gamma * lambda3(t));
                    
                    % State equations
                    s(t + 1) = max(0.0, min(1.0, s(t) - dt * alpha(t) * s(t) * i(t)));
                    i(t + 1) = max(0.0, min(1.0, i(t) + dt * (alpha(t) * s(t) * i(t) - beta * i(t))));
                    alpha(t + 1) = max(alpha_min, min(alpha_max, alpha(t) + dt * (-gamma * alpha(t) + gamma * a'*(u_max - u(:, t)))));
                end
                % % %                     % fix end point conditions and repeat (shooting method)
                % % %                     lambda1(1) = lambda1(1) - lambda1(end);
                % % %                     lambda2(1) = lambda2(1) - lambda2(end);
                % % %                     lambda3(1) = lambda3(1) - lambda3(end);
                % % %                 end
            end
            plot(mean(J0), mean(J1), 'ro');
            %             plot(epsilon, mean(J), 'ro');
        else % random inputs
            % time evolution
            for t = 1 : K
                for jj = 1 : L
                    u(jj, t) = randi([u_min(jj), u_max(jj)]);
                end
                
                % Costs
                J0(t) = alpha(t) * s(t) * i(t);
                J1(t) = w(:, t)' * u(: , t);
                J(t) = (1 - epsilon) * J0(t) + epsilon * J1(t);
                
                if(t < K)
                    % State equations
                    s(t + 1) = max(0.0, min(1.0, s(t) - dt * alpha(t) * s(t) * i(t)));
                    i(t + 1) = max(0.0, min(1.0, i(t) + dt * alpha(t) * s(t) * i(t) - dt * beta * i(t)));
                    alpha(t + 1) = max(alpha_min, min(alpha_max, alpha(t) + dt * (-gamma * alpha(t) + gamma * a'*(u_max - u(:, t)))));
                end
            end
            plot(mean(J0), mean(J1), 'bo');
            %             plot(epsilon, mean(J), 'bo');
        end
    end
    disp(multiobjective_ind);
end
grid

t = (0 : K - 1) * dt;
figure

subplot(611);
plot(t, s)
hold on
plot(t, i)
plot(t, s.*i.*alpha)
grid
legend('s', 'i', 'n');

subplot(612);
plot(t, lambda1)
hold on
plot(t, lambda2)
plot(t, lambda3)
grid
legend('\lambda_1', '\lambda_2', '\lambda_3');

subplot(613);
plot(t, alpha)
grid
legend('alpha');

subplot(614);
plot(t, H)
grid
legend('Hamiltonian');

subplot(615);
plot(t, (1 - epsilon) * J0)
hold on
plot(t, epsilon * J1)
plot(t, J)
grid
legend('(1 - epsilon) * J0 (Human)', 'epsilon * J1 (Economic)', 'J (Total)');

subplot(616);
plot(t, u');
grid
legend('inputs');
