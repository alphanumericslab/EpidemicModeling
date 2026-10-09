function [human_cost, intervention_cost] = npi_cost(new_cases, inputs, weights)
    % NPI_COST Compute mean daily case burden and mean weighted intervention intensity.
    % [human_cost,intervention_cost] = npi_cost(new_cases,inputs,weights)
    % new_cases: case count vector; inputs: interventions-by-time; weights:
    % scalar, intervention column vector, or matching matrix. Intervention cost
    % averages over both days and interventions (the historical convention).
    % Author: Reza Sameni | Emory University

    if isvector(weights) && ~isscalar(weights)
        weights = weights(:);
    end

    human_cost = mean(new_cases(:));
    intervention_cost = mean(weights .* inputs, 'all');
end
