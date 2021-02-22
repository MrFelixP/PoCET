% This example briefly illustrates how PoCET can be used for solving
% stochastic model predictive control (SMPC) problems. In (1) the system
% is defined and transformed into a PCE system. (2) just entails a
% simulation of the autonomous system. The parameters of the MPC problem
% are defined in (3.1) and the actual "online" control is performed in
% (3.b). The value function is defined in the separate file 'ex6_optfun.m'
% while the nonlinear contraints are defined in 'ex6_confun.m'.
%
% !! Please note that (for currently unknown reasons) fmincon seems to
% have severe problems handling the constrained optimization problem so
% the current implementation in (3.b) only solves the unconstrained 
% problem using fminsearch. We encourage you to try different solvers.

clc, clear all, close all
addpath(genpath('../../../PoCET'));

%% (1.a) define system variables
%  For input instructions enter 'help PoCET' in your command window
states(1).name = 'x';
states(1).dist = 'uniform';
states(1).data = [40 45]; % lower and upper bound
states(1).rhs  = '-a*x + b*u';

parameters(1).name = 'a';
parameters(1).dist = 'uniform'; 
parameters(1).data = [4 5]; % lower and upper bound

parameters(2).name = 'b';
parameters(2).dist = 'uniform'; 
parameters(2).data = [9.9 10.1]; % mean and standard deviation

inputs(1).name = 'u';
inputs(1).rhs  = 'piecewise(ut,uv,t)';
inputs(1).ut = 1; % time of step
inputs(1).uv = 0; % height of step

%% (1.b) define simulation options
simopt.tspan = [0, 3];
simopt.dt = 0.1;
simopt.setup = odeset;
simopt.solver = 'ode15s';

mc_samples = 1e3;
col_samples = 1e3;
pce_order = 3;

%% (1.c) compose the PCE system and write file for expanded ODE
% mySystem = PoCETcompose(states, parameters, inputs, options);
sys = PoCETcompose(states,parameters,inputs,[],pce_order);
PoCETwriteFiles(sys,'ex6_ODE.m')

% get moment matrices up to 4th moment 
% (we use up tp 4th moment to fit a beta distribution later)
MomMats = PoCETmomentMatrices(sys,4);

%% (2) simulate autonomous system (uv=0) and plot results
% run PCE simulation
results = PoCETsimGalerkin(sys,'ex6_ODE',[],simopt,'uv',0);

% compute moments from simulation results
results.x.moments = PoCETcalcMoments(sys,MomMats,results.x.pcvals,'central',2);

figure(1)
subplot(2,1,1) % plot mean of state x
plot(results.time,results.x.moments(1,:),'r')
hold on; title('Mean of state x')

subplot(2,1,2) % plot variance of state x
plot(results.time,results.x.moments(2,:),'r')
hold on; title('Variance of state x')

%% (3.a) Set up MPC
% define parameters for the MPC
mpc_param.t_F = 20;    % total simulation time (in seconds)
mpc_param.t_H = 5;     % prediction horizon (in seconds)
mpc_param.dt = 1;    % stepsize for piecewise constant inputs

mpc_param.x_F = 75;    % target value for the MEAN of x
mpc_param.x_ub = 80;   % upper bound for state x
mpc_param.p_ub = 0.95; % maximum probability of violating upper bound for x

mpc_param.Q = 10;       % weighting for state cost

mpc_param.simopt = simopt; % use previous simulation options
mpc_param.simopt.tspan = [0 mpc_param.t_H]; % define prediction horizon
mpc_param.simopt.dt = mpc_param.dt; % evaluate state every dt seconds

% declare piecewise constant input signal
ut = 0; %0:mpc_param.dt:mpc_param.t_H; % set up time-vector for inputs
uv0 = 0;% zeros(size(ut)); % initial values for inputs

% change distribution type of state x to 'pce' (required for updating IC)
mpc_sys = sys;
mpc_sys.states(1).dist = 'pce';
mpc_sys.states(1).data = sys.states(1).pce;

%% (3.b) Perform MPC
% (re-)set initial condition and set-point
mpc_sys = PoCETupdate(mpc_sys,'x',sys.states(1).pce); % IC for x
mpc_param.x_F = 75; % initial set-point

% simulation variables
tt = 0; % global time in seconds
kk = 1; % current time index
mpc_param.simopt.tspan = [0, mpc_param.dt]; % time span of one step
uv_opt = zeros(1, mpc_param.t_F/mpc_param.dt + 1); % preallocation
x_opt = zeros(1, mpc_param.t_F/mpc_param.dt + 1);  % preallocation

% initialize plotting
figure(2); 
plot([0, mpc_param.t_F],[mpc_param.x_F, mpc_param.x_F],'k--'); hold on
plot([0, mpc_param.t_F],[mpc_param.x_ub, mpc_param.x_ub],'r')

while tt < mpc_param.t_F
 if tt == 10
  mpc_param.x_F = 50; % at tt = 10: change set point
 end
 
 % Solve unconstraint problem
 uv_opt_tt = fminsearch(@(uv)ex6_optfun(mpc_sys,mpc_param,MomMats,ut,uv),uv0);

 % fmincon seems to have problems solving the optimization problem;
 % unfortunately we don't know why at the moment but encourage you to try
 % different solvers
%  uv_opt_tt = fmincon(@(uv)ex6_optfun(mpc_sys,mpc_param,MomMats,ut,uv),uv0, ...
%                      [],[],[],[],-100*ones(size(ut)),100*ones(size(ut)),...
%                      @(uv)ex6_confun(mpc_sys,mpc_param,MomMats,ut,uv),options);
 
 % apply first input of the solution to the system
 uv_opt(kk) = uv_opt_tt(1);
 sim = PoCETsimGalerkin(mpc_sys,'ex6_ODE',[],simopt,'ut',0,'uv',uv_opt(kk));
 
 % calculate moments of x and plot mean
 moments = PoCETcalcMoments(mpc_sys,MomMats,sim.x.pcvals);
 figure(2); plot(tt+sim.time,moments(1,:),'b'); pause(0.01);
 
 % for visualization: fit and plot beta-4 distribution
 ablu = calcBeta4(moments(:,end));
 dx = linspace(ablu(3)-0.01,ablu(4)+0.01,100);
 dy = beta4(ablu,dx);
 figure(3); plot3(tt*ones(size(dy)),dx,dy,'b'); hold on; pause(0.01); 
 
 % update of initial condition for x (new initialization of MPC)
 mpc_sys = PoCETupdate(mpc_sys,'x',sim.x.pcvals(:,end));
          
 tt = tt + mpc_param.dt;
 kk = kk + 1;
end

figure(2); title('mean of state x'); grid on; hold off

figure(3); xlabel('time tt'); ylabel('state x'); zlabel('PDF \mu(x)');
title('fitted distribution of state x'); grid on; hold off
