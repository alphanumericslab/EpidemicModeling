function [fs, Cs, fw, Cw] = StateHessianTerms(u, s_k, Pk, w_bar, Qk, params)
fs = zeros(6, 1);
Cs = zeros(6);

fw = zeros(6, 1);
Cw = zeros(6);

end