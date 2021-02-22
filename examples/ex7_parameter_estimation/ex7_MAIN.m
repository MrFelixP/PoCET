% This is a simple example for parameter estimation with PoCET.
% Assume you have a measured data set and a dynamic model for the
% corresponding process. The goal is to fit the distribution of the
% system's state x to the measured data by finding optimal values for 
% the parameter a. 

clc, clear all, close all

%% (0) generate data to be fitted
data.samples = 40 + 2*randn(500,1); % generate data set 
data.mean = mean(data.samples); % mean of measured data 
data.std = std(data.samples); % standard deviation (std) of measured data 
data.time = 1; % time of the measurement

%% (1.a) define system to be fitted
%  For input instructions enter 'help PoCET' in your command window
states(1).name = 'x';
states(1).dist = 'none'; % assume no initial uncertainties
states(1).data = [100]; % initial value of x
states(1).rhs  = '-a*x'; % right-hand-side of ode (\dot{x} = rhs)

parameters(1).name = 'a';
parameters(1).dist = 'uniform'; % assumed distribution of uncertain parameter
parameters(1).data = [0.9 1.1]; % initial guesses of mean and std of a

%% (1.b) define simulation options
simoptions.tspan = [0 data.time]; 
simoptions.dt = 0.005;
simoptions.setup = odeset;
simoptions.solver = 'ode15s';

pce_order = 3;

%% (1.c) compose the PCE system and write file for PCE-ODE and MC simulations
% mySystem = PoCETcompose(states, parameters, inputs, options);
sys = PoCETcompose(states,parameters,[],[],pce_order);

PoCETwriteFiles(sys,'ex6_ODE.m')

%% (2.) estimate lb and ub using moments of x
% get moment matrices up to 2nd moment
MomMats = PoCETmomentMatrices(sys,2);

% optimize parameters of a (see below for definition of optfun_1)
% Note: the constraint [1 -1]*[lb; ub] <= -0.1 guarantees that the lower
% bound is always smaller than the uper bound
a_opt = fmincon(@(a)optfun(sys,data,a,MomMats,simoptions), parameters(1).data, [1 -1], -0.1);
disp(['Optimal paramters of parameter a found as    ' num2str(a_opt)])

%% (3.) plotting the results
figure(1)
histogram(data.samples,50,'Normalization','pdf') % plot initial data
hold on

% define normal distribution odf for plotting
ndfcn = @(x,mu,sd) exp(-(x-mu).^2 ./ (2*sd)) /(sqrt(2*sd*pi));

% simulate system with initial guesses for a
results_init = PoCETsimGalerkin(sys,'ex1_ODE',simoptions);
moments_init = PoCETcalcMoments(sys,MomMats,results_init.x.pcvals);
plot(30:0.5:50,ndfcn(30:0.5:50,moments_init(1,end),moments_init(2,end)),'LineWidth',2)

% update PCE system with optimal values for a
sys_opt = PoCETupdate(sys,'a',a_opt);
results_opt = PoCETsimGalerkin(sys_opt,'ex1_ODE',simoptions);
moments_opt = PoCETcalcMoments(sys_opt,MomMats,results_opt.x.pcvals);
plot(30:0.5:50,ndfcn(30:0.5:50,moments_opt(1,end),moments_opt(2,end)),'LineWidth',2)

legend('measured data','initial distribution','fitted distribution')

%% define optfun
% Note that fitting via only the first two moments is generally not a very
% robust procedure and might yield incorrect results since higher moments
% are disregarded. You could also adapt the optimization procoedure from
% example 4 (discrimination) which uses another distance measure for the
% involved distributions.
function value = optfun(sys,data,param_a,MomMats,simoptions)
 % update PCE system with current guess for a
 sys = PoCETupdate(sys,'a',param_a);

 % simulate the updated system
 sim = PoCETsimGalerkin(sys,'ex1_ODE',simoptions);

 % compute moments of x at the final time of the simulation 
 sim.x.moments = PoCETcalcMoments(sys,MomMats,sim.x.pcvals(:,end));
 
 % calculate quality of current guess 
 % (the 10* is for numeric normalization)
 value = (data.mean - sim.x.moments(1))^2 + 100*(data.std - sim.x.moments(2))^2;
end

