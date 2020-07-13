% This example shows the (1) preparation and (2) simulation of an LTI 
% system using PoCET featuring a shorthand method to create the PoCET input
% structure for the states and the explicite Euler method for integration.
% Furthermore, part (3) shows different plotting techniques for Monte Carlo
% simulations (single trajectories, quantiles, and lower/upper bounds).

clc, clear all, close all
addpath(genpath('../../../PoCET'));

%% (1.a) define system variables
% define LTI system
x = str2sym('[x1;x2;x3]');
A = str2sym('[-R1,0,0; R1,-R2,0; 0,R2,-R3]'); % use commas!
B = sym([1,0;0,1;0,0]);
C = sym([0,0,1]);
D = sym([0,0]);
u = str2sym('[u1;u2]');

xdot = A*x + B*u;
yfun = C*x + D*u;

% define initial conditions
% rem: size(xdist) and size(xdata) have to equal size(x)
xdist = {'normal';'none';'none'};
xdata = {[3 0.1] ;    2 ;    0 };

% create input structs for PoCET from LTI system
% shorthand method using 'struct' and 'num2cell'
states = struct('name',num2cell(x),'rhs',num2cell(xdot),'dist',xdist,'data',xdata);
outputs = struct('name','y','rhs',yfun);

inputs(1).name = u(1);
inputs(1).rhs  = 'piecewise(ut1,uv1,t)';
inputs(1).ut1  = [0 0.5 1 1.5 2 2.5]; % time of step
inputs(1).uv1  = [1  0  1  0  1  0 ]; % height of step

inputs(2).name = u(2);
inputs(2).rhs  = '-uv2';
inputs(2).uv2  = 1; % height of step

parameters(1).name = 'R1';
parameters(1).dist = 'normal'; 
parameters(1).data = [1 0.01];

parameters(2).name = 'R2';
parameters(2).dist = 'normal'; 
parameters(2).data = [2 0.05];

parameters(3).name = 'R3';
parameters(3).dist = 'normal'; 
parameters(3).data = [1.5 0.02];

%% (1.b) define simulation options
simoptions.tspan = [0 3];
simoptions.dt = 0.005;
simoptions.setup = odeset;
simoptions.solver = 'explicit-euler';

mc_samples = 1e3;
pce_order = 2;

%% (1.c) compose the PCE system and write file for PCE-ODE and MC simulations
% mySystem = PoCETcompose(states, parameters, inputs, options);
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
MomMats = PoCETmomentMatrices(sys,2);
PoCETwriteFiles(sys,'ex2_ODE.m','ex2_OUT.m','ex2_MCODE.m','ex2_MCOUT.m')

%% (2.a) simulate system and compute moments
% run PCE simulation
results = PoCETsimGalerkin(sys,'ex2_ODE','ex2_OUT',simoptions);

% compute moments from simulation results and store also in results
results = PoCETcalcMoments(sys,MomMats,results);

% run Monte-Carlo simulations
samples = PoCETsample(sys,'variables',mc_samples);
mcresults = PoCETsimMonteCarlo(sys,'ex2_MCODE','ex2_MCOUT',samples,simoptions,'method','moments');

%% (2.b) plot results of both simulations
figure(1)
subplot(2,1,1) % plot mean of states
plot(results.time,results.x1.moments(1,:),'r') 
hold on;
plot(results.time,results.x2.moments(1,:),'b')
plot(results.time,results.x3.moments(1,:),'g')
title('Mean')
legend('x_1','x_2','x_3')

subplot(2,1,2) % plot variance of states
plot(results.time,results.x1.moments(2,:),'r') 
hold on;
plot(results.time,results.x2.moments(2,:),'b')
plot(results.time,results.x3.moments(2,:),'g')
title('Variance')
legend('x_1','x_2','x_3')

leg1 = sprintf('PCE (order %i)',pce_order);
leg2 = sprintf('MC (%i samples)',mc_samples);

figure(2)
subplot(1,2,1) % plot mean of output
plot(results.time,results.y.moments(1,:),'r',mcresults.time,mcresults.y.moments(1,:),'b:')
hold on; title('Mean of y')
legend(leg1,leg2)

subplot(1,2,2) % plot variance of output
plot(results.time,results.y.moments(2,:),'r',mcresults.time,mcresults.y.moments(2,:),'b:')
hold on; title('Variance of y')
legend(leg1,leg2)

%% (3.a) plotting Monte-Carlo trajectories
% runing a Monte Carlo simulation without specifying a method returns the
% trajectories of the states/outputs
mccomplete = PoCETsimMonteCarlo(sys,'ex2_MCODE','ex2_MCOUT',1000,simoptions);

figure(3)
subplot(3,1,1), plot(mccomplete.time,mccomplete.x1.mcvals(1:5,:)), title('x_1')
subplot(3,1,2), plot(mccomplete.time,mccomplete.x2.mcvals(1:5,:)), title('x_2')
subplot(3,1,3), plot(mccomplete.time,mccomplete.x3.mcvals(1:5,:)), title('x_3')

%% (3.b) plotting quantiles
% to plot quantiles use the auxiliary function plot_quantile

figure(4)
subplot(3,1,1), plot_quantile(mccomplete.time,mccomplete.x1.mcvals,{0.01 0.7},'r'), title('x_1')
subplot(3,1,2), plot_quantile(mccomplete.time,mccomplete.x2.mcvals,{0.01 0.7},'b'), title('x_2')
subplot(3,1,3), plot_quantile(mccomplete.time,mccomplete.x3.mcvals,{0.01 0.7},'g'), title('x_3')

%% (3.b) plotting lower and upper bounds
% using Monte Carlo simulation with method 'moments' will also return the
% worst cases in the resulting struct

figure(5)
subplot(3,1,1), plot(mcresults.time,mcresults.x1.lb,'r',mcresults.time,mcresults.x1.ub,'r'), hold on, title('x_1')
subplot(3,1,2), plot(mcresults.time,mcresults.x2.lb,'b',mcresults.time,mcresults.x2.ub,'b'), hold on, title('x_2')
subplot(3,1,3), plot(mcresults.time,mcresults.x3.lb,'g',mcresults.time,mcresults.x3.ub,'g'), hold on, title('x_3')

% using Monte Carlo simulations without any specified method you can simply
% use MATLAB's min() and max() functions
subplot(3,1,1), plot(mccomplete.time,min(mccomplete.x1.mcvals),'k:',mccomplete.time,max(mccomplete.x1.mcvals),'k:')
subplot(3,1,2), plot(mccomplete.time,min(mccomplete.x2.mcvals),'k:',mccomplete.time,max(mccomplete.x2.mcvals),'k:')
subplot(3,1,3), plot(mccomplete.time,min(mccomplete.x3.mcvals),'k:',mccomplete.time,max(mccomplete.x3.mcvals),'k:')
