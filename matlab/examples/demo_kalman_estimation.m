function demo_kalman_estimation()
% DEMO_KALMAN_ESTIMATION Example 03: EKF/EKS and missing observations.
% Author: Reza Sameni | Emory University
% Reference: Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
day = 0:79; truth = 30*exp(.025*day); observed = truth+8*sin(day); observed(31:35) = nan;
[prior,filtered,~,~,gain,smoothed] = rt_exp_fit_ekf(observed,[30;.02],[1 1 .2],[0;0],0,diag([64 .001]),diag([4 1e-5]),64,1,1,14,1);
figure('Color','w'); tiledlayout(1,2);
nexttile; plot(day,[truth;filtered(1,:);smoothed(1,:)]','LineWidth',2); hold on; scatter(day,observed,15,'filled'); grid on; xlabel('Day'); ylabel('Cases'); title('A reporting gap'); legend('Truth','EKF','EKS','Observed');
nexttile; plot(day,[filtered(2,:);smoothed(2,:)]','LineWidth',2); grid on; xlabel('Day'); ylabel('Growth/day'); title('Hidden growth'); legend('EKF','EKS');
assert(all(gain(:,:,31:35)==0,'all')); assert(isequal(filtered(:,31:35),prior(:,31:35)));
end
