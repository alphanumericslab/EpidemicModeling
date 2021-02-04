% Observation equation Jacobian
function [C, D] = ObsJacobian(u, s_k, v_bar, params)
C = [s_k(2)*s_k(3), s_k(1)*s_k(3), s_k(1)*s_k(2), 0 , 0, 0];
D = 1;
end