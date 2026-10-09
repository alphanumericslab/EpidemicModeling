function [C, D] = obs_jacobian(u, s_k, v_bar, params)
% OBS_JACOBIAN Evaluate analytic observation and measurement-noise Jacobians.
C = [s_k(2)*s_k(3), s_k(1)*s_k(3), s_k(1)*s_k(2), 0 , 0, 0];
D = 1;
end
