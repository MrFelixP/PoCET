% This example shows the preparation (1), simulation (2), update (3), and
% simulation with inputs (4) of a stochastical system using PoCET.
% In (2), galerkin PCE is compared with Monte-Carlo-Simulation. 
% In (3) and (4), galerkin PCE is compared with collocation PCE.

clc, clear all, close all
addpath(genpath('../../../PoCET'));

%% (1.a) define system variables
%  For input instructions enter 'help PoCET' in your command window
states(1).name = 'x';
states(1).dist = 'uniform';
states(1).data = [4 5]; % lower and upper bound
states(1).rhs  = '-a^2*x^2 + b*u';

parameters(1).name = 'a';
parameters(1).dist = 'uniform'; 
parameters(1).data = [1 2]; % lower and upper bound

parameters(2).name = 'b';
parameters(2).dist = 'normal'; 
parameters(2).data = [1 0.01]; % mean and standard deviation

inputs(1).name = 'u';
inputs(1).rhs  = 'piecewise(ut,uv,t)';
inputs(1).ut = 1; % time of step
inputs(1).uv = 1; % height of step

outputs(1).name = 'y';
outputs(1).rhs  = '3*x^2';

%% (1.b) define simulation options
simoptions.tspan = [0 3];
simoptions.dt = 0.005;
simoptions.setup = odeset;
simoptions.solver = 'ode15s';

mc_samples = 1e3;
col_samples = 1e3;
pce_order = 3;

%% (1.c) compose the PCE system and write file for PCE-ODE and MC simulations
% mySystem = PoCETcompose(states, parameters, inputs, options);
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);

PoCETwriteFiles(sys,'ex1_ODE.m','ex1_OUT.m','ex1_MCODE.m','ex1_MCOUT.m')

%% (2.a) simulate system and compute moments
% run PCE simulation
results = PoCETsimGalerkin(sys,'ex1_ODE','ex1_OUT',simoptions);

% get moment matrices up to 4th moment
MomMats = PoCETmomentMatrices(sys,4);

% compute moments from simulation results
results.x.moments = PoCETcalcMoments(sys,MomMats,results.x.pcvals);
results.y.moments = PoCETcalcMoments(sys,MomMats,results.y.pcvals);

% run Monte-Carlo simulations
samples = PoCETsample(sys,'variables',mc_samples);
mcresults = PoCETsimMonteCarlo(sys,'ex1_MCODE','ex1_MCOUT',samples,simoptions,'method','moments');

%% (2.b) plot results of both simulations
figure(1)
subplot(2,1,1) % plot mean of state x
plot(results.time,results.x.moments(1,:),'r',mcresults.time,mcresults.x.moments(1,:),'b:')
hold on; title('Mean of state x')

subplot(2,1,2) % plot variance of state x
plot(results.time,results.x.moments(2,:),'r',mcresults.time,mcresults.x.moments(2,:),'b:')
hold on; title('Variance of state x')

figure(2)
subplot(2,1,1) % plot mean of output y
plot(results.time,results.y.moments(1,:),'r',mcresults.time,mcresults.y.moments(1,:),'b:')
hold on; title('Mean of output y')

subplot(2,1,2) % plot variance of output y
plot(results.time,results.y.moments(2,:),'r',mcresults.time,mcresults.y.moments(2,:),'b:')
hold on; title('Variance of output y')

%% (3) Change initial condition: x(0) = 2
% update initial condition for state x
sys = PoCETupdate(sys,'x',[6 6.5]);
results2 = PoCETsimGalerkin(sys,'ex1_ODE','ex1_OUT',simoptions);
results2.x.moments = PoCETcalcMoments(sys,MomMats,results2.x.pcvals);
results2.y.moments = PoCETcalcMoments(sys,MomMats,results2.y.pcvals);

% collocation method for deriving pce coefficients
colloc2 = PoCETsimCollocation(sys,'ex1_MCODE','ex1_MCOUT',col_samples,simoptions);
colloc2.x.moments = PoCETcalcMoments(sys,MomMats,colloc2.x.pcvals);
colloc2.y.moments = PoCETcalcMoments(sys,MomMats,colloc2.y.pcvals);

figure(1) % plot state x
subplot(2,1,1); plot(results.time,results2.x.moments(1,:),'g',colloc2.time,colloc2.x.moments(1,:),'m:')
subplot(2,1,2); plot(results.time,results2.x.moments(2,:),'g',colloc2.time,colloc2.x.moments(2,:),'m:')

figure(2) % plot output y
subplot(2,1,1); plot(results.time,results2.y.moments(1,:),'g',colloc2.time,colloc2.y.moments(1,:),'m:')
subplot(2,1,2); plot(results.time,results2.y.moments(2,:),'g',colloc2.time,colloc2.y.moments(2,:),'m:')

%% (4) Use higher step: uv = 2
% simulate with different input parameter uv
results3 = PoCETsimGalerkin(sys,'ex1_ODE','ex1_OUT',simoptions,'uv',2);
results3.x.moments = PoCETcalcMoments(sys,MomMats,results3.x.pcvals);
results3.y.moments = PoCETcalcMoments(sys,MomMats,results3.y.pcvals);

% collocation method with different input parameter uv
colloc3 = PoCETsimCollocation(sys,'ex1_MCODE','ex1_MCOUT',col_samples,simoptions,'uv',2);
colloc3.x.moments = PoCETcalcMoments(sys,MomMats,colloc3.x.pcvals);
colloc3.y.moments = PoCETcalcMoments(sys,MomMats,colloc3.y.pcvals);

figure(1)
subplot(2,1,1); plot(results.time,results3.x.moments(1,:),'k',colloc3.time,colloc3.x.moments(1,:),'y--')
legend('PCE (x(0)=0; uv=1)','MC (x(0)=0; uv=1)','PCE (x(0)=2; uv=1)','COL (x(0)=2; uv=1)','PCE (x(0)=4; uv=2)','COL (x(0)=4; uv=2)')
subplot(2,1,2); plot(results.time,results3.x.moments(2,:),'k',colloc3.time,colloc3.x.moments(2,:),'m--')

figure(2)
subplot(2,1,1); plot(results.time,results3.y.moments(1,:),'k',colloc3.time,colloc3.y.moments(1,:),'y--')
legend('PCE (x(0)=0; uv=1)','MC (x(0)=0; uv=1)','PCE (x(0)=2; uv=1)','COL (x(0)=2; uv=1)','PCE (x(0)=4; uv=2)','COL (x(0)=4; uv=2)')
subplot(2,1,2); plot(results.time,results3.y.moments(2,:),'k',colloc3.time,colloc3.y.moments(2,:),'m--')
