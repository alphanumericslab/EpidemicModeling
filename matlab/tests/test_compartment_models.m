function tests = test_compartment_models()
% TEST_COMPARTMENT_MODELS Independent conservation, analytic and endpoint tests.
% Author: Reza Sameni | Emory University
 setup_once([]);
tests = functiontests(localfunctions);
end
function setup_once(test_case)
% SETUP_ONCE Add the function library.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'functions'));
end
function test_mass_conservation(test_case)
% TEST_MASS_CONSERVATION SEIRP preserves the total population fraction.
[s,e,i,r,p] = seirp(.4,.3,.12,.04,.09,.002,0,.998,.001,.001,0,0,100,.1);
verifyEqual(test_case,sum([s;e;i;r;p],1),ones(size(s)),'AbsTol',1e-12);
verifyGreaterThanOrEqual(test_case,min([s e i r p]),0);
end
function test_closed_form_removal(test_case)
% TEST_CLOSED_FORM_REMOVAL Euler decay has a known discrete solution.
[~,~,i,~,p] = seirp(0,0,0,0,.1,.02,0,.9,0,.1,0,0,10,.1);
expected = .1*(1-.1*.12).^(0:99);
verifyEqual(test_case,i,expected,'AbsTol',1e-12);
verifyEqual(test_case,p,(.1-expected)*(.02/.12),'AbsTol',1e-12);
end
function test_saturation_equal_rates(test_case)
% TEST_SATURATION_EQUAL_RATES Equal healthcare endpoints reduce to ordinary SEIRP.
[a,b,c,d,e] = seirp(.4,.3,.12,.04,.09,.002,0,.998,.001,.001,0,0,20,.1);
[f,g,h,i,j] = seirp_saturated_resource(.4,.3,.12,.04,0,.998,.001,.001,0,0,20,.1,.09,.09,.002,.002,.01,.04);
verifyEqual(test_case,[a;b;c;d;e],[f;g;h;i;j],'AbsTol',1e-12);
end
function test_si_decay(test_case)
% TEST_SI_DECAY Zero contact rate produces geometric infected-state decay.
[s,i] = si_controlled(0,.1,.9,.1,100,.1);
verifyEqual(test_case,s,.9*ones(1,100),'AbsTol',1e-12);
verifyEqual(test_case,i,.1*.99.^(0:99),'AbsTol',1e-12);
end
function test_shared_noise_equilibrium(test_case)
% TEST_SHARED_NOISE_EQUILIBRIUM Deterministic alpha remains at its equilibrium.
[s,~,alpha] = si_alpha_controlled(ones(2,30),.999,.001,.3,[3;3],0,5,1/7,[.06;.04],.1,.1,0,0,0,30,1,zeros(3,30));
verifyEqual(test_case,alpha,.3*ones(1,30),'AbsTol',1e-12); verifyLessThan(test_case,s(1),.999);
end
