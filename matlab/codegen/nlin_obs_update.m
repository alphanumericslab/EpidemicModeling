function x_k = nlin_obs_update(u, s_k, v_bar, params)
% NLIN_OBS_UPDATE Evaluate the configured new-case or total-case observation.
x_k = s_k(1) * s_k(2) * s_k(3) + v_bar;
end
