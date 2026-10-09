function demo_intervention_control()
% DEMO_INTERVENTION_CONTROL Example 04: training, forecasts and bounded control.
% Author: Reza Sameni | Emory University
% Reference: Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
population = 1e6; day = 0:99;
u = [2*(day>=35)-(day>=70);2*(day>=50)];
truth.population = population; truth.params = default_si_params([3;3]);
truth.params.a = [.06;.04]; truth.params.b = .02;
truth.state = [.9999;.0001;.32]; truth.covariance = diag([1e-8 1e-8 .01]);
[new_cases,~] = forecast_npi(truth,u); model = fit_npi_model(cumsum(new_cases),u,population,[3;3]);
disp(table(model.coefficients,model.refined_coefficients,'VariableNames',{'FirstPass','Refined'}));
horizon = 35; weights = [1;1.5]; fixed = repmat(u(:,end),1,horizon);
[fixed_cases,~] = forecast_npi(model,fixed); [controls,controlled_cases,~] = optimal_npi(model,horizon,weights,3e-5);
figure('Color','w'); tiledlayout(1,2);
nexttile; plot(0:horizon-1,[fixed_cases;controlled_cases]','LineWidth',2); grid on; title('Scenario forecasts'); xlabel('Forecast day'); ylabel('Cases/day'); legend('Continue policy','Control solver');
nexttile; stairs(0:horizon-1,controls','LineWidth',2); grid on; title('Bounded policy'); xlabel('Forecast day'); ylabel('NPI level'); legend('Intervention 1','Intervention 2');
assert(all(controls>=0,'all') && all(controls<=3,'all'));
end
