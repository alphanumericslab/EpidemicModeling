function [fs, Cs, fw, Cw] = state_hessian_terms(u, s_k, Pk, w_bar, Qk, params)
% STATE_HESSIAN_TERMS Standalone six-state SI-alpha helper for MATLAB Coder.
% Author: Reza Sameni | Emory University
% Shape and parameter contracts match docs/api.md; see Python codegen module.
fs = zeros(6, 1);
Cs = zeros(6);

fw = zeros(6, 1);
Cw = zeros(6);

end
