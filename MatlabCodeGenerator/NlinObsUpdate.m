% Nonlinear observation update
function x_k = NlinObsUpdate(u, s_k, v_bar, params)
x_k = s_k(1) * s_k(2) * s_k(3) + v_bar;
end