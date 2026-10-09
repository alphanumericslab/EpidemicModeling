function tests = test_kalman_filters()
% TEST_KALMAN_FILTERS Scalar Kalman oracle, missing data and second-order smoke tests.
setup_once([]);
tests = functiontests(localfunctions);
end
function setup_once(test_case)
% SETUP_ONCE Add function library.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
end
function handles = linear_handles()
% LINEAR_HANDLES Construct a one-state random-walk model.
handles.state_hard_margins = @(s,p,k) s; handles.obs_hard_margins = @(x,p,k) x;
handles.nlin_state_update = @linear_state_update; handles.nlin_obs_update = @(u,s,v,p,k) s+v;
handles.state_jacobians = @linear_jacobians; handles.obs_jacobian = @linear_jacobians;
handles.state_hessian_terms = @zero_hessians; handles.obs_hessian_terms = @zero_hessians;
end
function [u,next] = linear_state_update(u,s,w,p,k)
% LINEAR_STATE_UPDATE Random-walk mean transition.
next = s+w;
end
function [a,b] = linear_jacobians(varargin)
% LINEAR_JACOBIANS Unit state and noise Jacobians.
a = 1; b = 1;
end
function [a,b,c,d] = zero_hessians(varargin)
% ZERO_HESSIANS No second-order correction for a linear model.
a = 0; b = 0; c = 0; d = 0;
end
function test_linear_oracle(test_case)
% TEST_LINEAR_ORACLE Compare the generic filter with direct scalar recursion.
x = [1 2 nan 1.5]; [~,~,~,sp,~,~,pp,~,gain] = generic_extended_kalman_filter(zeros(1,4),x,linear_handles(),struct(),0,1,nan,nan,0,0,.1,.5,1,1,3,1);
mean_state = 0; variance = 1;
for k = 1:numel(x)
    if isfinite(x(k))
        g = variance/(variance+.5); mean_state = mean_state+g*(x(k)-mean_state); variance = (1-g)*variance;
    end
    verifyEqual(test_case,sp(k),mean_state,'AbsTol',1e-12); verifyEqual(test_case,pp(k),variance,'AbsTol',1e-12); variance = variance+.1;
end
verifyEqual(test_case,gain(:,:,3),0);
end
function test_variable_covariance_and_boundary(test_case)
% TEST_VARIABLE_COVARIANCE_AND_BOUNDARY Scalar variance series and terminal constraints.
[~,~,~,~,smooth,~,~,cov] = generic_extended_kalman_filter(zeros(1,4),1:4,linear_handles(),struct(),0,1,9,.01,0,0,(1:4)*.1,(1:4)*.2,1,1,3,2);
verifyEqual(test_case,smooth(end),9); verifyEqual(test_case,cov(end),.01);
end
function test_exponential_missing_both_orders(test_case)
% TEST_EXPONENTIAL_MISSING_BOTH_ORDERS Missing observations skip correction in both orders.
truth = 25*exp(.02*(0:39)); x = truth+sin(0:39); x(11:13) = nan;
for order = 1:2
    [~,~,~,~,gain,smooth] = rt_exp_fit_ekf(x,[25;.02],[1 1 .2],[0;0],0,diag([4 .001]),diag([1 1e-5]),4,1,1,7,order);
    verifyEqual(test_case,gain(:,:,11:13),zeros(2,1,3)); verifyLessThan(test_case,mean(abs(smooth(1,:)-truth)),2);
end
end
