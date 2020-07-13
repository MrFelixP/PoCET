% This example demonstrates the active fault detection for two fault 
% scenarios using the collocation PCE method. The workflow is (1) defining 
% the models, (2) simulating the nominal and fault scenarios and
% (3) finding a discriminating input.
% The PDFs are approximated using the PCE coefficients from the collocation
% method and a set of basis samples, which is evaluated beforehand.

clc, clear all, close all
addpath(genpath('../../../PoCET'));
fprintf(2,'NOTE: Upcoming warnings about unused variables can be ignored - still needs to be fixed!\n');

simoptions.tspan = [0 3000];
simoptions.dt = 50; % in continuous-time this only defines the stepwidth for plotting
simoptions.setup = odeset;
simoptions.solver = 'ode15s';

pce_order = 2;
col_samples = 30; % number of collocation points
n_samples = 5000; % number of basis evaluation samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) defining the scenarios

states(1).name = 'x1';
states(1).dist = 'none';
states(1).data = 0.2;
states(1).rhs = '(-c1*Sp*sign(x1-x3)*sqrt(2*g*abs(x1-x3))+u1)/A';

states(2).name = 'x2';
states(2).dist = 'none';
states(2).data = 0.15;
states(2).rhs = '(-c3*Sp*sign(x2-x3)*sqrt(2*g*abs(x2-x3))-c2*Sp*sqrt(2*g*x2)-Qf+u2)/A';

states(3).name = 'x3';
states(3).dist = 'none';
states(3).data = 0.1;
states(3).rhs = '(c1*Sp*sign(x1-x3)*sqrt(2*g*abs(x1-x3))-c3*Sp*sign(x3-x2)*sqrt(2*g*abs(x3-x2)))/A';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_org(1).name = 'c1';
param_org(1).dist = 'normal'; 
param_org(1).data = [1, sqrt(0.0025)];

param_org(2).name = 'c2';
param_org(2).dist = 'normal'; 
param_org(2).data = [1, sqrt(0.0025)];

param_org(3).name = 'c3';
param_org(3).dist = 'normal'; 
param_org(3).data = [1, sqrt(0.0025)];

param_org(4).name = 'Qf';
param_org(4).dist = 'none'; 
param_org(4).data =  0;

param_org(5).name = 'Sp';
param_org(5).dist = 'none';
param_org(5).data = 5e-5;

param_org(6).name = 'g';
param_org(6).dist = 'none'; 
param_org(6).data = 9.81;

param_org(7).name = 'A';
param_org(7).dist = 'none';
param_org(7).data = 0.0154;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_fauA(1).name = 'c1';
param_fauA(1).dist = 'normal'; 
param_fauA(1).data = [1, sqrt(0.0025)];

param_fauA(2).name = 'c2';
param_fauA(2).dist = 'normal'; 
param_fauA(2).data = [1, sqrt(0.0025)];

param_fauA(3).name = 'c3';
param_fauA(3).dist = 'normal'; 
param_fauA(3).data = [1, sqrt(0.0025)];

param_fauA(4).name = 'alph';
param_fauA(4).dist = 'normal'; 
param_fauA(4).data =  [0.6, 4e-4];

param_fauA(5).name = 'Sp';
param_fauA(5).dist = 'none';
param_fauA(5).data = 5e-5;

param_fauA(6).name = 'g';
param_fauA(6).dist = 'none'; 
param_fauA(6).data = 9.81;

param_fauA(7).name = 'A';
param_fauA(7).dist = 'none';
param_fauA(7).data = 0.0154;

param_fauA(8).name = 'Qf';
param_fauA(8).dist = 'none'; 
param_fauA(8).data =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_fauB(1).name = 'c1';
param_fauB(1).dist = 'normal'; 
param_fauB(1).data = [1, sqrt(0.0025)];

param_fauB(2).name = 'c2';
param_fauB(2).dist = 'normal'; 
param_fauB(2).data = [1, sqrt(0.0025)];

param_fauB(3).name = 'c3';
param_fauB(3).dist = 'normal'; 
param_fauB(3).data = [1, sqrt(0.0025)];

param_fauB(4).name = 'r';
param_fauB(4).dist = 'normal'; 
param_fauB(4).data =  [0.002, 1e-6];

param_fauB(5).name = 'Sp';
param_fauB(5).dist = 'none';
param_fauB(5).data = 5e-5;

param_fauB(6).name = 'g';
param_fauB(6).dist = 'none'; 
param_fauB(6).data = 9.81;

param_fauB(7).name = 'A';
param_fauB(7).dist = 'none';
param_fauB(7).data = 0.0154;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputs(1).name = 'y';
outputs(1).rhs  = 'x3';

inputs(1).name = 'u1';
inputs(1).rhs  = 'piecewise(ut1,uv1,t)';
inputs(1).ut1 = [0 1000 2000]; % time of step
inputs(1).uv1 = [1    1    1]*1e-5; % height of step

inputs(2).name = 'u2';
inputs(2).rhs  = 'piecewise(ut2,uv2,t)';
inputs(2).ut2 = [0 1000 2000]; % time of step
inputs(2).uv2 = [1    1    1]*4e-5; % height of step

lb_in = 1e-5; % lower input bound
ub_in = 1e-4; % upper input bound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Composing nominal system model...\n');
sys_org = PoCETcompose(states,param_org,inputs,outputs,pce_order);
writeMCRHSfile(sys_org,'ex5_ODE_org.m');
writeMCOUTfile(sys_org,'ex5_OUT_org.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputs(1).rhs  = 'alph*piecewise(ut1,uv1,t)';

fprintf('\nComposing fault scenario A...\n');
sys_fauA = PoCETcompose(states,param_fauA,inputs,outputs,pce_order);
writeMCRHSfile(sys_fauA,'ex5_ODE_fauA.m');
writeMCOUTfile(sys_fauA,'ex5_OUT_fauA.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputs(1).rhs  = 'piecewise(ut1,uv1,t)';
inputs(3).name = 'Qf';
inputs(3).rhs  = 'c2*pi*r^2*sqrt(2*g*x2)';

fprintf('\nComposing fault scenario B...\n');
sys_fauB = PoCETcompose(states,param_fauB,inputs,outputs,pce_order);
writeMCRHSfile(sys_fauB,'ex5_ODE_fauB.m');
writeMCOUTfile(sys_fauB,'ex5_OUT_fauB.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nEvaluating basis samples...\n')
basis_evals_org  = PoCETsample(sys_org, 'basis',n_samples);
basis_evals_fauA = PoCETsample(sys_fauA,'basis',n_samples);
basis_evals_fauB = PoCETsample(sys_fauB,'basis',n_samples);

fprintf('\nEvaluating collocation samples...\n')
col.org  = PoCETsample(sys_org, 'basis',col_samples);
col.fauA = PoCETsample(sys_fauA,'basis',col_samples);
col.fauB = PoCETsample(sys_fauB,'basis',col_samples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) simulating the scenarios
fprintf('\nSimulating undiscriminated scenarios...\n')
Org_PC  = PoCETsimCollocation(sys_org ,'ex5_ODE_org' ,'ex5_OUT_org' ,col.org,simoptions,'method','final');
FauA_PC = PoCETsimCollocation(sys_fauA,'ex5_ODE_fauA','ex5_OUT_fauA',col.fauA,simoptions,'method','final');
FauB_PC = PoCETsimCollocation(sys_fauB,'ex5_ODE_fauB','ex5_OUT_fauB',col.fauB,simoptions,'method','final');
fprintf('\n')

Org_MC_y  = basis_evals_org' *Org_PC.y.pcvals;
FauA_MC_y = basis_evals_fauA'*FauA_PC.y.pcvals;
FauB_MC_y = basis_evals_fauB'*FauB_PC.y.pcvals;

lb = min([Org_MC_y; FauA_MC_y; FauB_MC_y]);
ub = max([Org_MC_y; FauA_MC_y; FauB_MC_y]);

EDGES = linspace(lb,ub,100); % Intervals used to generate the histograms
N_Org_MC  = histc(Org_MC_y ,EDGES); % histogram of Model 1 (constructed with MC simulations) at the end point of simulation
N_FauA_MC = histc(FauA_MC_y,EDGES);
N_FauB_MC = histc(FauB_MC_y,EDGES);

% plot output y
figure(1)
plot(EDGES,N_Org_MC,'b',EDGES,N_FauA_MC,'r',EDGES,N_FauB_MC,'k')
title('undiscriminated case')
legend('Nominal Scenario','Fault Scenario A','Fault Scenario B')

pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) finding discriminating input
% optimize input
fprintf('Performing the optimization...\n')
vzero  = @(x)0;
value  = @(x)sqrt(x*eye(length(x))*x');
constr = @(x)ex5_objfun(x,sys_org,sys_fauA,sys_fauB,col,basis_evals_org,basis_evals_fauA,basis_evals_fauB,simoptions);
ops = optimset(optimset,'Display','iter-detailed','MaxIter',2,'Algorithm','interior-point');

uv_0 = [inputs(1).uv1 inputs(2).uv2];
uv_opt = fmincon(constr,uv_0,[],[],[],[],lb_in*ones(size(uv_0)),ub_in*ones(size(uv_0)),[],ops);

uv1 = uv_opt(1:numel(uv_opt)/2);
uv2 = uv_opt(numel(uv_opt)/2+1:end);

% collocation method for deriving pce coefficients
Org_PC  = PoCETsimCollocation(sys_org, 'ex5_ODE_org' ,'ex5_OUT_org' ,col.org,simoptions,'uv1',uv1,'uv2',uv2,'method','final');
FauA_PC = PoCETsimCollocation(sys_fauA,'ex5_ODE_fauA','ex5_OUT_fauA',col.fauA,simoptions,'uv1',uv1,'uv2',uv2,'method','final');
FauB_PC = PoCETsimCollocation(sys_fauB,'ex5_ODE_fauB','ex5_OUT_fauB',col.fauB,simoptions,'uv1',uv1,'uv2',uv2,'method','final');

Org_MC_y  = basis_evals_org'*Org_PC.y.pcvals;
FauA_MC_y = basis_evals_fauA'*FauA_PC.y.pcvals;
FauB_MC_y = basis_evals_fauB'*FauB_PC.y.pcvals;

lb = min([Org_MC_y; FauA_MC_y; FauB_MC_y]);
ub = max([Org_MC_y; FauA_MC_y; FauB_MC_y]);

EDGES = linspace(lb,ub,100); % Intervals used to generate the histograms
N_Org_MC  = histc(Org_MC_y ,EDGES);
N_FauA_MC = histc(FauA_MC_y,EDGES);
N_FauB_MC = histc(FauB_MC_y,EDGES);

% plot output y
figure(2)
plot(EDGES,N_Org_MC,'b',EDGES,N_FauA_MC,'r',EDGES,N_FauB_MC,'k')
title('discriminated case')
legend('Nominal Scenario','Fault Scenario A','Fault Scenario B')


