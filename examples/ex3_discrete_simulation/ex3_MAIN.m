% This example shows the declaration (1), simulation (2) and update (3) of
% a discrete-time stochastical system with five uncertainties (two inital 
% conditions, three parameters).

clear all, close all
addpath(genpath('../../../PoCET'));

%% (1) define states, parameters and simulation options
%  For input instructions enter 'help PoCET' in your command window
states(1).name = 'A';
states(1).dist = 'uniform';
states(1).data = [1.2 1.3];
states(1).rhs = 'p_1*A + 0.5*B';

states(2).name = 'B';
states(2).dist = 'uniform';
states(2).data = [0.9 1.3];
states(2).rhs = 'p_2*B + R';

parameters(1).name = 'p_1';
parameters(1).dist = 'uniform';
parameters(1).data = [0.3 0.6];

parameters(2).name = 'p_2';
parameters(2).dist = 'uniform';
parameters(2).data = [0.4 0.5];

parameters(3).name = 'R';
parameters(3).dist = 'uniform';
parameters(3).data = [0.2 0.4];

% define options
simoptions.tspan = [0 0.1];
simoptions.dt = 0.01;
simoptions.setup = odeset;
simoptions.solver = 'discrete-time'; % choosing 'discrete-time' as solver 
% is the only difference between continous- and discrete-time simulation

pce_order = 2;
mc_samples = 1e4;

% compose the PCE system and write files for galerkin-PCE and MC simulation
% mySystem = PoCETcompose(states,parameters,inputs,outputs,order);
sys = PoCETcompose(states,parameters,[],[],pce_order);
PoCETwriteFiles(sys,'ex3_RHS',[],'ex3_MCRHS');

%% (2) simulate system
% run galerkin-PCE simulation
results = PoCETsimGalerkin(sys,'ex3_RHS',[],simoptions);

% compute moments from simulation results
sys.MomMats = PoCETmomentMatrices(sys,4);
results.A.moments = PoCETcalcMoments(sys,sys.MomMats,results.A.pcvals);
results.B.moments = PoCETcalcMoments(sys,sys.MomMats,results.B.pcvals);

% run Monte-Carlo simulations
mcresults = PoCETsimMonteCarlo(sys,'ex3_MCRHS',[],mc_samples,simoptions,'method','moments');

% plot results of both simulations
figure(1) % central moments of state 1
subplot(2,2,1); plot(results.time,results.A.moments(1,:),'r',mcresults.time,mcresults.A.moments(1,:),'b:'); title('Mean')
subplot(2,2,2); plot(results.time,results.A.moments(2,:),'r',mcresults.time,mcresults.A.moments(2,:),'b:'); title('Variance')
subplot(2,2,3); plot(results.time,results.A.moments(3,:),'r',mcresults.time,mcresults.A.moments(3,:),'b:'); title('Skewness')
subplot(2,2,4); plot(results.time,results.A.moments(4,:),'r',mcresults.time,mcresults.A.moments(4,:),'b:'); title('Excess Kurtosis')

figure(2) % central moments of state 2
subplot(2,2,1); plot(results.time,results.B.moments(1,:),'r',mcresults.time,mcresults.B.moments(1,:),'b:'); title('Mean')
subplot(2,2,2); plot(results.time,results.B.moments(2,:),'r',mcresults.time,mcresults.B.moments(2,:),'b:'); title('Variance')
subplot(2,2,3); plot(results.time,results.B.moments(3,:),'r',mcresults.time,mcresults.B.moments(3,:),'b:'); title('Skewness')
subplot(2,2,4); plot(results.time,results.B.moments(4,:),'r',mcresults.time,mcresults.B.moments(4,:),'b:'); title('Excess Kurtosis')

%% (3) Updating initial conditions and parameter values
sys = PoCETupdate(sys,'A',[3 5],'p_2',[0.7 0.9]);

% run Galerkin-PCE simulation
results2 = PoCETsimGalerkin(sys,'ex3_RHS',simoptions);

% compute moments from simulation results (alternative function call)
results2 = PoCETcalcMoments(sys,sys.MomMats,results2);

% run Monte-Carlo simulations
mcresults2 = PoCETsimMonteCarlo(sys,'ex3_MCRHS',[],mc_samples,simoptions,'method','moments');

% plot results of both simulations
figure(1) % central moments of state 1
subplot(2,2,1); plot(results2.time,results2.A.moments(1,:),'r',mcresults2.time,mcresults2.A.moments(1,:),'b:'); title('Mean')
subplot(2,2,2); plot(results2.time,results2.A.moments(2,:),'r',mcresults2.time,mcresults2.A.moments(2,:),'b:'); title('Variance')
subplot(2,2,3); plot(results2.time,results2.A.moments(3,:),'r',mcresults2.time,mcresults2.A.moments(3,:),'b:'); title('Skewness')
subplot(2,2,4); plot(results2.time,results2.A.moments(4,:),'r',mcresults2.time,mcresults2.A.moments(4,:),'b:'); title('Excess Kurtosis')

figure(2) % central moments of state 2
subplot(2,2,1); plot(results2.time,results2.B.moments(1,:),'r',mcresults2.time,mcresults2.B.moments(1,:),'b:'); title('Mean')
subplot(2,2,2); plot(results2.time,results2.B.moments(2,:),'r',mcresults2.time,mcresults2.B.moments(2,:),'b:'); title('Variance')
subplot(2,2,3); plot(results2.time,results2.B.moments(3,:),'r',mcresults2.time,mcresults2.B.moments(3,:),'b:'); title('Skewness')
subplot(2,2,4); plot(results2.time,results2.B.moments(4,:),'r',mcresults2.time,mcresults2.B.moments(4,:),'b:'); title('Excess Kurtosis')



