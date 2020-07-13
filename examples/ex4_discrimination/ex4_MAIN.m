% This example demonstrates the discrimination of two model candidates
% using PoCET. The workflow is (1) defining the models, (2) simulating
% the undiscriminated models and (3) finding a discriminating input.
% In order to compare the PDFs a fourparametric beta-distribution is fitted
% to the moments.

clc, clear all, close all
addpath(genpath('../../../PoCET'));
parametersH = zeros(3,1);

%% (1) define model candidates, input, and options
% system 1
statesH(1).name = 'x_1'; % name
statesH(1).rhs = '(p_1+p_3)*(x_2-1)*x_1+(p_2+u)*x_2'; % right hand side of ODE
statesH(1).dist = 'beta4'; % initial distribution
statesH(1).data = [3 3 0.96 0.98]; % initial distribution parameters

statesH(2).name = 'x_2'; % name
statesH(2).rhs = 'p_1*(1-x_2)*x_1-(p_2+u)*x_2'; % right hand side of ODE
statesH(2).dist = 'beta4'; % initial distribution
statesH(2).data = [3 3 0.01 0.03]; % initial distribution parameters

for i = 1:3
 parametersH(i).name = ['p_' num2str(i)]; % name
 parametersH(i).dist = 'uniform'; % distribution
 parametersH(i).data = [0.9, 1.1]; % distribution parameters
end

% parametersH(2).name = 'p_2'; % name
% parametersH(2).dist = 'uniform'; % distribution
% parametersH(2).data = [0.9, 1.1]; % distribution parameters
% 
% parametersH(3).name = 'p_3'; % name
% parametersH(3).dist = 'uniform'; % distribution
% parametersH(3).data = [0.9, 1.1]; % distribution parameters

% system 2
statesM(1).name = 'x_1'; % name
statesM(1).dist = 'beta4'; % initial distribution 
statesM(1).data = [3 3 0.96 0.98]; % initial distribution parameters
statesM(1).rhs = 'p_1*(x_2-1)*x_1+(p_2+u1)*x_2';

statesM(2).name = 'x_2'; % name
statesM(2).dist = 'beta4'; % initial distribution 
statesM(2).data = [3 3 0.01 0.03]; % initial distribution parameters
statesM(2).rhs = 'p_1*(1-x_2)*x_1-(p_3+p_2+u1)*x_2';

parametersM(1).name = 'p_1'; % name
parametersM(1).dist = 'uniform'; % distribution
parametersM(1).data = [0.9, 1.15]; % distribution parameters

parametersM(2).name = 'p_2'; % name
parametersM(2).dist = 'uniform'; % distribution
parametersM(2).data = [0.9, 1.15]; % distribution parameters

parametersM(3).name = 'p_3'; % name
parametersM(3).dist = 'uniform'; % distribution
parametersM(3).data = [0.9, 1.15]; % distribution parameters

% define piecewise constant input
inputs(1).name = 'u'; % name
inputs(1).rhs  = 'piecewise(u_t,u_v,t)'; % right hand side (any MATLAB function possible)
inputs(1).u1_t = [0 0.5 1 1.5 2 2.5 3 3.5]; % vector of step times
inputs(1).u1_v = [0 0 0 0 0 0 0 0]; % initial vector of step sizes (preallocation)

% define simulation options
mc_samples = 1e3;


simoptions.setup = odeset;
simoptions.solver = 'ode15s';
simoptions.tspan = [0, 10];
simoptions.dt = 0.05;

% compose the PCE systems and write files for PCE-ODE simulations
pce_order = 2;
tic
sys_H = PoCETcompose(statesH,parametersH,inputs,[],pce_order); % calculate system expansion
sys_H.MomMats = PoCETmomentMatrices(sys_H,4); % calculate matrices for moment calculation
toc
keyboard
PoCETwriteFiles(sys_H,'ex4_ODE_H'); % write .m-function files for simulation

sys_M = PoCETcompose(statesM,parametersM,inputs,[],pce_order); % calculate system expansion
sys_M.MomMats = PoCETmomentMatrices(sys_M,4); % calculate matrices for moment calculation
PoCETwriteFiles(sys_M,'ex4_ODE_M'); % write .m-function files for simulation

%% (2) Galerkin PCE simulation of undiscriminated systems

nondisc_H = PoCETsimGalerkin(sys_H,'ex4_ODE_H',[],simoptions); 
nondisc_moments_H = PoCETcalcMoments(sys_H,sys_H.MomMats,nondisc_H.x_2.pcvals(:,end)); 
nondisc_beta_H = calcBeta4(nondisc_moments_H);

nondisc_M = PoCETsimGalerkin(sys_M,'ex4_ODE_M',[],simoptions);
nondisc_moments_M = PoCETcalcMoments(sys_M,sys_M.MomMats,nondisc_M.x_2.pcvals(:,end));
nondisc_beta_M = calcBeta4(nondisc_moments_M);


%% (3) Finding optimal discriminating input
% Because the problem is nonlinear and nonconvex, fmincon can only
% guarantee a feasible solution if the initial value is already feasible.
% In order to find such a feasible initial condition, first run fmincon
% with a objective function the always returns zero. Then run fmincon with
% an actual objective function using the solution from step one as the
% initial condition in order to find a cost-optimal discriminating input.
% (For step one the recommended algorithm is 'interior-point',
%  for step two the recommended algorithm is 'active-set'.)
%
% For model discrimination, the file 'ex4_mycon.m' can be used as a template 
% for the nonlinear constraint used in fmincon.

%  1. find any discriminating input
u0 = [1 .5 0 0 0 0 0 0];
ops = optimset(optimset,'Display','iter-detailed','MaxIter',25,'Algorithm','interior-point');
zerocon = @(x)0;
constr = @(x)ex4_mycon(x,sys_H,sys_M,simoptions);
u_feas = fmincon(zerocon,u0,[],[],[],[],zeros(size(u0)),5*ones(size(u0)),constr,ops);

%  find a cost-optimal discriminating input
u0 = [1 .5 0 0 0 0 0 0];
ops = optimset(optimset,'Display','iter-detailed','MaxIter',25,'Algorithm','active-set');
penalty = @(x) sqrt(x*eye(length(x))*x');
u_opt = fmincon(penalty,u_feas,[],[],[],[],zeros(size(u_feas)),5*ones(size(u_feas)),constr,ops);

%%
% Galerkin PCE simulation with discriminating input
disc_H = PoCETsimGalerkin(sys_H,'ex4_ODE_H',[],simoptions,'u1_v',u_opt);
disc_moments_H = PoCETcalcMoments(sys_H,sys_H.MomMats,disc_H.x_2.pcvals(:,end));
disc_beta_H = calcBeta4(disc_moments_H);

disc_M = PoCETsimGalerkin(sys_M,'ex4_ODE_M',[],simoptions,'u1_v',u_opt);
disc_moments_M = PoCETcalcMoments(sys_M,sys_M.MomMats,disc_M.x_2.pcvals(:,end));
disc_beta_M = calcBeta4(disc_moments_M);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeMCRHSfile(sys_H,'ex4_MCODE_H');
writeMCRHSfile(sys_M,'ex4_MCODE_M');

mcresultsH = PoCETsimMonteCarlo(sys_H,'ex4_MCODE_H',[],mc_samples,simoptions,'method','moments');
mcresultsM = PoCETsimMonteCarlo(sys_M,'ex4_MCODE_M',[],mc_samples,simoptions,'method','moments');

mcresultsH_d = PoCETsimMonteCarlo(sys_H,'ex4_MCODE_H',[],mc_samples,simoptions,'u1_v',u_opt,'method','moments');
mcresultsM_d = PoCETsimMonteCarlo(sys_M,'ex4_MCODE_M',[],mc_samples,simoptions,'u1_v',u_opt,'method','moments');

%%
minlb = min([nondisc_beta_H(3),nondisc_beta_M(3)])-0.001;
maxub = max([nondisc_beta_H(4),nondisc_beta_M(4)]);
x = linspace(minlb,maxub,1000);
y_H = beta4(nondisc_beta_H,x);
y_M = beta4(nondisc_beta_M,x);

figure(1)
subplot(2,1,1)
plot(x,y_H,'Color','b','Linewidth',1.5)
hold on
plot(x,y_M,'Color','r','Linewidth',1.5,'LineStyle','--')
legend({'$x_2^H$ (Henri kinetics)','$x_2^M$ (Michaelis-Menten kinetics)'},'Interpreter','latex','FontSize',15)
xlabel('Concentration of the substrate-enzyme-complex $x_2^*$','Interpreter','latex','FontSize',15)
ylabel('PDFs $\mu(x_2^*)$','Interpreter','latex','FontSize',15)
title('(a) PDFs before discrimination','Interpreter','latex','FontSize',15)
xlim([0.005 0.026])

minlb = min([disc_beta_H(3),disc_beta_M(3)]);
maxub = max([disc_beta_H(4),disc_beta_M(4)]);
x = linspace(minlb,maxub,1000);
y_H = beta4(disc_beta_H,x);
y_M = beta4(disc_beta_M,x);

figure(1)
subplot(2,1,2)
plot(x,y_H,'Color','b','Linewidth',1.5)
hold on
plot(x,y_M,'Color','r','Linewidth',1.5,'LineStyle','--')
% legend({'$x_2^H$ (Henri kinetics)','$x_2^M$ (Michaelis-Menten kinetics)'},'Interpreter','latex','FontSize',13)
xlabel('Concentration of the substrate-enzyme-complex $x_2^*$','Interpreter','latex','FontSize',15)
ylabel('PDFs $\mu(x_2^*)$','Interpreter','latex','FontSize',15)
xlim([0.0017 0.0275])
title('(b) PDFs after discrimination','Interpreter','latex','FontSize',15)

%%

figure(3)
subplot(2,1,1)
plot(mcresultsH.time,mcresultsH.x_2.lb,'Color','b','Linewidth',1.5)
hold on
plot(mcresultsM.time,mcresultsM.x_2.lb,'Color','r','Linewidth',1.5,'LineStyle','--')
plot(mcresultsH.time,mcresultsH.x_2.ub,'Color','b','Linewidth',1.5)
plot(mcresultsM.time,mcresultsM.x_2.ub,'Color','r','Linewidth',1.5,'LineStyle','--')
legend('System 1','System 2')
legend({'$x_2^H$ (Henri kinetics)','$x_2^M$ (Michaelis-Menten kinetics)'},'Interpreter','latex','FontSize',15)
ylabel('Concentration $x_2^*$','Interpreter','latex','FontSize',15)
xlabel('time $t$','Interpreter','latex','FontSize',15)
title('(a) Outer bounds before discrimination','Interpreter','latex','FontSize',15)

subplot(2,1,2)
plot(mcresultsH_d.time,mcresultsH_d.x_2.lb,'Color','b','Linewidth',1.5)
hold on
plot(mcresultsM_d.time,mcresultsM_d.x_2.lb,'Color','r','Linewidth',1.5,'LineStyle','--')
plot(mcresultsH_d.time,mcresultsH_d.x_2.ub,'Color','b','Linewidth',1.5)
plot(mcresultsM_d.time,mcresultsM_d.x_2.ub,'Color','r','Linewidth',1.5,'LineStyle','--')
ylabel('Concentration $x_2^*$','Interpreter','latex','FontSize',15)
xlabel('time $t$','Interpreter','latex','FontSize',15)
title('(b) Outer bounds after discrimination','Interpreter','latex','FontSize',15)

