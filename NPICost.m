function [J0, J1, J] = NPICost(newcases, inputs, weights, epsilon)
% Calculates the NPI cost for a given set of weights and bi-objective
% weights

% Human factor cost
J0 = mean(newcases);

% NPI costs
weighted_inputs = weights .* inputs;
J1 = mean(weighted_inputs(:));

% Overall bi-objective cost
J = (1 - epsilon) * J0 + epsilon * J1;
