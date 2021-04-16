% A simulation of the particle diffusion PDE numeric solution for a moving subject
% A comparison of two discretization methods
%
% Reza Sameni
% Apr 2021

clear
close all
clc

Lx = 51; % grid dimension 1
Ly = 51; % grid dimension 2
dx = 0.1; % in meters
dy = 0.1; % in meters

sd1 = 2; % SD of subject 1's random walk motion
sd2 = 1.5; % SD of subject 1's random walk motion
E1 = 0.1; % magnitude of particles release by subject 1 (through exhalation)
E2 = 0.2; % magnitude of particles release by subject 2 (through exhalation)

dt = 1/60.0; % simulation time unit
NT = round(10.0/dt + 1); % total simulation time iterations

D = 1e-1; % diffusion parameter

C1 = zeros(Lx, Ly); % method 1 concentration
C2 = zeros(Lx, Ly); % method 2 concentration
X1 = zeros(2, NT); % subject 1 coordinates
X2 = zeros(2, NT); % subject 2 coordinates

% Ininitialization
% first subject's initial point
x1_init = round(Lx/2);
y1_init = round(Lx/2);

% second subject's initial point
x2_init = round(Lx/7);
y2_init = round(Lx/5);
X1 = [x1_init ; y1_init];
X2 = [x2_init ; y2_init];

C1(x1_init, y1_init) = E1; % Point 1
C1(x2_init, y2_init) = E2; % Point 2

C2 = C1; % both methods start with the same initial conditions

alpha_x = D * dt / dx^2;
alpha_y = D * dt / dy^2;
alpha_xy = D * dt / (dx^2 + dy^2);

% Check validity of the parameters
if((1 - 2 * alpha_x - 2 * alpha_y) < 0 || (1 - 2 * alpha_x - 2 * alpha_y - 2 * alpha_xy - 2 * alpha_xy) < 0)
    error('Stability condition for parameters not fulfilled. Make simulation time period smaller');
end

h = figure;
i = 2 : Lx - 1; % horizontal indexes
j = 2 : Ly - 1; % vertical indexes

writerObj = VideoWriter('DiffusionVideoExample.avi');
writerObj.FrameRate = 60;
open(writerObj);

% Update all points using two methods
for t = 1 : NT - 1
    % Particles added by the moving subjects
    
    % 4-point neighborhood
    C1(X1(1), X1(2)) = C1(X1(1), X1(2)) + E1; % subject 1
    C1(X2(1), X2(2)) = C1(X2(1), X2(2)) + E2; % subject 2
    
    % 8-point neighborhood
    C2(X1(1), X1(2)) = C2(X1(1), X1(2)) + E1; % subject 1
    C2(X2(1), X2(2)) = C2(X2(1), X2(2)) + E2; % subject 2

    % Update subject positions for t + 1
    X1 = max(1, min(Lx, round(X1 + sd1 * (rand(2, 1) - 0.5))));              
    X2 = max(1, min(Ly, round(X2 + sd2 * (rand(2, 1) - 0.5))));
    
    % 4-point neighborhood update
    C1(i, j) = (1 - 2 * alpha_x - 2 * alpha_y) * C1(i, j) + alpha_x * C1(i - 1, j) + alpha_x * C1(i + 1, j) + alpha_y * C1(i, j - 1) + alpha_y * C1(i, j + 1);
    
    % 8-point neighborhood update
    C2(i, j) = (1 - 2 * alpha_x - 2 * alpha_y - 2 * alpha_xy - 2 * alpha_xy) * C2(i, j) + alpha_x * C2(i - 1, j) + alpha_x * C2(i + 1, j) + alpha_y * C2(i, j - 1) + alpha_y * C2(i, j + 1) +   ...
                      alpha_xy * C2(i - 1, j - 1) + alpha_xy * C2(i + 1, j + 1) + alpha_xy * C2(i - 1, j + 1) + alpha_xy * C2(i + 1, j - 1);
    
    figure(h);
    subplot(121)
    imagesc(C1);
    axis square
    title(['4 neighbor update, ITR = ', num2str(t)]);
    
    subplot(122)
    imagesc(C2);
    axis square
    title(['8 neighbor update, ITR = ', num2str(t)]);
    
    %pause(0.01);
    frame = getframe(h);
    writeVideo(writerObj,frame);
end

close(writerObj);
% movie(M)