function s_k = state_hard_margins(s_k, params)
% STATE_HARD_MARGINS Constrain population fractions and contact rate; preserve costates.
s_k(1) = min(1.0, max(0, s_k(1)));
s_k(2) = min(1.0, max(0, s_k(2)));
s_k(3) = min(params.alpha_max, max(params.alpha_min, s_k(3)));
end
