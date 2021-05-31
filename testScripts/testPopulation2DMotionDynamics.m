% A simulation of the diffusion PDE numeric solution for a population with
% accelerated motion with socializing/anti-socializing forces that keep
% them at an optimal distance.
%
% Reza Sameni
% Apr 2021

% function [p, v, a, j] = Population2DMotionDynamics(dt, T, environment_plan, p_inits, v_inits, a_inits, j_inits, eta_p_std, eta_v_std, eta_a_std, eta_j_std)
clear
close all
clc

n_subjects = 150; % number of subjects
dt = 0.1; % simulation time unit
T = round(250.0/dt); % total simulation time steps
dl = 0.3; % spatial resolution in meters (1 foot)

% min/max range of the simulated environment
x_min = 0; % m
y_min = 0; % m
x_max = 70.0; % m
y_max = 50.0; % m

v_max = 2.5; % max velocity (m/s)
a_max = 1.44; % max acceleration (m/s^2)
j_max = 0.6; % max jerk m/s^3 % https://www.mdpi.com/2079-9292/8/9/943/htm
eta_p_std = 0.1; % random position displacement std (m)
eta_v_std = 0.1/sqrt(2); % random velocity std (m/s)
eta_a_std = 0.01/sqrt(2); % random acceleration std (m/s^2)
eta_j_std = 0  *  0.1/sqrt(2); % random jerk std (m/s^3)

% Some refs for the above selections:
% https://www.researchgate.net/post/What_is_the_maximum_walking_acceleration_deceleration_over_a_very_short_time_period_eg_002_01_05_sec
% Most data I've seen is related to people's acceleration rates in specific situations. For example T. Korhonen collected data on pedestrian flow characteristics and observed that the acceleration distribution "bear a resemblance" to a normal distribution with an average of 0.68 m/s2 with a maximum of  1.44 m/s2
% In many models the acceleration parameter is calibrated to achieve a certain flow rate characteristic. I've seen ranges like 0.5 to 1.1 seconds for agents to achieve their maximum velocity.
% A 2016 paper by Peng Lin et al., see link below, does discuss acceleration in FDS+Evac and Usain Bolt's maximum acceleration which they quote as 3.09m/s2
% http://cpb.iphy.ac.cn/fileup/PDF/2016-3-034501.pdf
% http://dynardo.de/fileadmin/Material_Dynardo/dokumente/flyer/Verfahrenstechnik/RDO_Journal_02_2014_P_Evacuation.pdf

Lx = ceil((x_max - x_min)/dl); % environment X length
Ly = ceil((y_max - y_min)/dl); % environment Y length
environment_plan = zeros(Lx, Ly);
environment_plan(1, :) = 1;
environment_plan(end, :) = 1;
environment_plan(:, 1) = 1;
environment_plan(:, end) = 1;
environment_plan(round(10.0/dl):round(25.0/dl), round(10.0/dl)) = 1;
environment_plan(round(10.0/dl), round(10.0/dl):round(30.0/dl)) = 1;
% environment_plan(round(7.0/dl):round(25.0/dl), round(10.0/dl)) = 1;
% environment_plan(round(25.0/dl), round(10.0/dl) : round(30.0/dl)) = 1;

% p_inits = cat(2, randi([x_min x_max], n_subjects, 1), randi([y_min y_max], n_subjects, 1));
p_inits = cat(2, (x_max - x_min) * rand(n_subjects, 1) + x_min, (y_max - y_min) * rand(n_subjects, 1) + y_min);
v_inits = min(v_max/sqrt(2), 0.5 * randn(n_subjects, 2));
a_inits = min(a_max/sqrt(2), 0.1 * randn(n_subjects, 2));
j_inits = 0    *     min(j_max/sqrt(2), 0.075 * randn(n_subjects, 2));

N = size(p_inits, 1); % number of subjects in the environment
[obstacles_x_ind, obstacles_y_ind] = find(environment_plan);
obstacles_x = (obstacles_x_ind - 1) * dl + x_min;
obstacles_y = (obstacles_y_ind - 1) * dl + y_min;

p = zeros(N, 2, T); % positions of all subjects
v = zeros(N, 2, T); % velocities of all subjects
a = zeros(N, 2, T); % accelerations of all subjects
j = zeros(N, 2, T); % jerks of all subjects

p(:, :, 1) = p_inits;
v(:, :, 1) = v_inits;
a(:, :, 1) = a_inits;
j(:, :, 1) = j_inits;

marked_subjects = 1 : 2;
unmarked_subjects = 1 : N; unmarked_subjects(marked_subjects) = [];
E = 1.0; % magnitude of particles release by infected subjects through exhalation
D = 1e-2; % diffusion parameter
C = zeros(Lx, Ly); % particle concentration
init_x_ind = max(1, min(Lx, round((p_inits(marked_subjects, 1) - x_min)/dl)));
init_y_ind = max(1, min(Ly, round((p_inits(marked_subjects, 2) - y_min)/dl)));
C(init_x_ind, init_y_ind) = E;
alpha_x = D * dt / dl^2;
alpha_y = D * dt / dl^2;
alpha_xy = D * dt / (dl^2 + dl^2);
II = 2 : Lx - 1; % horizontal indexes
JJ = 2 : Ly - 1; % vertical indexes
% Check validity of the parameters
if((1 - 2 * alpha_x - 2 * alpha_y) < 0 || (1 - 2 * alpha_x - 2 * alpha_y - 2 * alpha_xy - 2 * alpha_xy) < 0)
    error('Stability condition for parameters not fulfilled. Make simulation time period smaller');
end


h = figure;
% figure_position = [100 100 540 400];
writerObj = VideoWriter('SocializingPopulation.avi');
writerObj.FrameRate = 10;%round(1/dt);
open(writerObj);

for t = 1 : T-1
    for n = 1 : N
        p(n, :, t + 1) = p(n, :, t) + dt * v(n, :, t) + eta_p_std * randn(1, 2);
        p(n, 1, t + 1) = max(x_min, min(x_max, p(n, 1, t + 1)));
        p(n, 2, t + 1) = max(y_min, min(y_max, p(n, 2, t + 1)));
        
        % Check if the pathway from t to t+1 has passed any obstacles
        xA = p(n, 1, t);
        xB = p(n, 1, t + 1);
        yA = p(n, 2, t);
        yB = p(n, 2, t + 1);
        
        %         obstacle_x_dist_from_path = obstacles_x - ((xA - xB) * obstacles_y + xB * yA - xA * yB) / (yA - yB);
        obstacle_y_dist_from_path = obstacles_y - ((yA - yB) * obstacles_x + yB * xA - yA * xB) / (xA - xB);
        
        %         I_path_obstacle_intersection = find(sqrt(obstacle_x_dist_from_path.^2 + obstacle_y_dist_from_path.^2) <= dl, 1, 'first');
        % I_path_obstacle_intersection = find(abs(obstacle_y_dist_from_path) <= dl & obstacle_y_dist_from_path <= y_max & obstacle_y_dist_from_path > y_min, 1, 'first');
        I_path_obstacle_intersection = find(abs(obstacle_y_dist_from_path) <= dl & obstacles_y <= max(yA, yB) & obstacles_y >= min(yA, yB) & obstacles_x <= max(xA, xB) & obstacles_x >= min(xA, xB), 1, 'first');
        
        if(environment_plan(max(1, min(Lx, round((p(n, 1, t + 1) - x_min)/dl))), max(1, min(Ly, round((p(n, 2, t + 1) - y_min)/dl)))) ...
                || ~isempty(I_path_obstacle_intersection)) % Undo the last move if crossed an obstacle
            p(n, :, t + 1) = p(n, :, t);
            %             v(n, :, t) = eta_v_std * randn(1, 2); % -v(n, :, t);
            %             a(n, :, t) = eta_a_std * randn(1, 2); % -a(n, :, t);
            %             j(n, :, t) = eta_j_std * randn(1, 2);
        end
        
        v(n, :, t + 1) = v(n, :, t) + dt * a(n, :, t) + eta_v_std * randn(1, 2);
        v(n, 1, t + 1) = max(-v_max, min(v_max, v(n, 1, t + 1)));
        v(n, 2, t + 1) = max(-v_max, min(v_max, v(n, 2, t + 1)));
        
        a(n, :, t + 1) = a(n, :, t) + dt * j(n, :, t) + eta_a_std * randn(1, 2);
        a(n, 1, t + 1) = max(-a_max, min(a_max, a(n, 1, t + 1)));
        a(n, 2, t + 1) = max(-a_max, min(a_max, a(n, 2, t + 1)));
        
        j(n, :, t + 1) = j(n, :, t) + eta_j_std * randn(1, 2);
        j(n, 1, t + 1) = max(-j_max, min(j_max, j(n, 1, t + 1)));
        j(n, 2, t + 1) = max(-j_max, min(j_max, j(n, 2, t + 1)));
    end
    
    % socializing/anti-socializing forces
    
    if(1)
        position_all_objects = [p(:, 1, t)' obstacles_x'; p(:, 2, t)' obstacles_y']';
        num_all_objects = size(position_all_objects, 1);
        
        object_interaction_matrix = zeros(num_all_objects);
        object_interaction_matrix(1 : n_subjects, 1 : n_subjects) = 10.0; % All subjects like to socialize
        object_interaction_matrix(n_subjects + 1: end, 1 : n_subjects) = -5.0; % People move away from physical obstacles
        object_interaction_matrix(1 : n_subjects, n_subjects + 1 : end) = -5.0; % Walls don't like people either!
        object_interaction_matrix(1 : num_all_objects + 1 : end) = 0; % set to zero the diagonal entries
        
        inter_subject_opt_dist = 6.0; % optimal social distance
        %         subject_obstacle_opt_dist = 3.0; % optimal obstacle distance
        no_further_effect_dist = 15.0;
        for n = 1 : N
            socializing_sign = zeros(num_all_objects, 1);
            
            subject_dist_vec_from_all_objects = position_all_objects(n * ones(1, num_all_objects), :) - position_all_objects;
            dist_vec_norm = vecnorm(subject_dist_vec_from_all_objects, 2, 2);
            
            %             socializing_sign(1 : n_subjects) = - sign(dist_vec_norm(1 : n_subjects) - inter_subject_opt_dist);
            socializing_sign(1 : n_subjects) = sign(inter_subject_opt_dist - dist_vec_norm(1 : n_subjects));
            socializing_sign(n_subjects + 1 : end) = -1;% No optimal distancing with the obstacles %dist_vec_norm(n_subjects + 1 : end) < subject_obstacle_opt_dist;
            
            socializing_sign(dist_vec_norm > no_further_effect_dist) = 0;
            
            Fx = (object_interaction_matrix(n, :)' .* socializing_sign .* subject_dist_vec_from_all_objects(:, 1)) ./ dist_vec_norm .^ 3;
            Fy = (object_interaction_matrix(n, :)' .* socializing_sign .* subject_dist_vec_from_all_objects(:, 2)) ./ dist_vec_norm .^3;
            %         socializing_forces(n, :, 1) = Fx;
            %         socializing_forces(n, :, 2) = Fy;
            a(n, :, t + 1) = a(n, :, t + 1) + [sum(Fx(isfinite(Fx))) sum(Fy(isfinite(Fy)))];
            a(n, 1, t + 1) = max(-a_max, min(a_max, a(n, 1, t + 1)));
            a(n, 2, t + 1) = max(-a_max, min(a_max, a(n, 2, t + 1)));
        end
    end
    
    
    x_ind = max(1, min(Lx, round((p(marked_subjects, 1, t) - x_min)/dl)));
    y_ind = max(1, min(Ly, round((p(marked_subjects, 2, t) - y_min)/dl)));
    indexes = sub2ind(size(C), x_ind, y_ind);
    C(indexes) = C(indexes) + E; % Exhale
    %     for kk = 1 : length(x_ind)
    %         C(x_ind(kk), y_ind(kk)) = C(x_ind(kk), y_ind(kk)) + E; % Exhale
    %     end
    
    %     % 4-point neighborhood update
    C(II, JJ) = (1 - 2 * alpha_x - 2 * alpha_y) * C(II, JJ) + alpha_x * C(II - 1, JJ) + alpha_x * C(II + 1, JJ) + alpha_y * C(II, JJ - 1) + alpha_y * C(II, JJ + 1);
    
    % 8-point neighborhood update
    C(II, JJ) = (1 - 2 * alpha_x - 2 * alpha_y - 2 * alpha_xy - 2 * alpha_xy) * C(II, JJ) + alpha_x * C(II - 1, JJ) + alpha_x * C(II + 1, JJ) + alpha_y * C(II, JJ - 1) + alpha_y * C(II, JJ + 1) +   ...
        alpha_xy * C(II - 1, JJ - 1) + alpha_xy * C(II + 1, JJ + 1) + alpha_xy * C(II - 1, JJ + 1) + alpha_xy * C(II + 1, JJ - 1);
    
    
    %     socializing_sign(1:20, 1: 20)
    %         plan = zeros(Lx, Ly);
    %         obstacles = sub2ind(obstacles_x, obstacles_y);
    %         plan(obstacles) = 1;
    plan = 0.9 * environment_plan;
    for n = 1 : N
        plan(max(1, min(Lx, round((p(n , 1, t) - x_min)/dl))), max(1, min(Lx, round((p(n , 2, t) - y_min) /dl)))) = 1;
    end
    figure(h);
    
    %     imshow(plan);
    
    subplot(121)
    plot(max(x_min, min(x_max, p(unmarked_subjects , 1, t))), max(y_min, min(y_max, p(unmarked_subjects , 2, t))), 'bo');
    hold on
    plot(max(x_min, min(x_max, p(marked_subjects , 1, t))), max(y_min, min(y_max, p(marked_subjects , 2, t))), 'ro');
    plot(obstacles_x, obstacles_y, 'k+', 'markersize', 0.5);
    hold off
    %     axis square
    axis equal
    axis tight
    xlim([x_min x_max])
    ylim([y_min y_max])
    xlabel('meters');
    ylabel('meters');
    %     title(['2D motion without social distancing policy, time = ' num2str((t+1)*dt) '(s)']);
    %     set(gca, 'fontsize', 14);
    title('Population motion');
    
    subplot(122)
    %     imagesc(C(:, end : -1 : 1)');
    %     imagesc(linspace(x_min, x_max, Lx), linspace(y_min, y_max, Ly), C(:, end : -1 : 1)');
    imagesc(linspace(x_min, x_max, Lx), linspace(y_min, y_max, Ly), C');
    xlim([x_min x_max])
    ylim([y_min y_max])
    set(gca,'YDir','normal')
    xlabel('meters');
    ylabel('meters');
    hold on
    plot(max(x_min, min(x_max, p(marked_subjects , 1, t))), max(y_min, min(y_max, p(marked_subjects , 2, t))), 'ro');
    %     axis square
    axis equal
    axis tight
    title('Particle concentrations');
    %     set(gca, 'fontsize', 14);
    sgtitle(['Population socializing at ' num2str(inter_subject_opt_dist) 'm optimal distance, t = ' num2str((t+1)*dt) '(s)'], 'fontsize', 18);
    
    %     pause(0.01)
    frame = getframe(h);
    writeVideo(writerObj, frame);
end
close(writerObj);
