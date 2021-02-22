function [J0, J1] = NPICost(newcases, inputs, weights)
% Calculates the NPI cost for a given set of weights and bi-objective
% weights

% Human factor cost
J0 = mean(newcases);

% NPI costs
weighted_inputs = weights .* inputs;
J1 = mean(weighted_inputs(:));
