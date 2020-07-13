function [obj] = ex5_objfun(x,sys_org,sys_fauA,sys_fauB,col,basis_evals_org,basis_evals_fauA,basis_evals_fauB,simoptions)
uv1 = x(1:numel(x)/2);
uv2 = x(numel(x)/2+1:end);

% get PoCET coefficients
fprintf('Nominal Scenario - ')
Org_PC  = PoCETsimCollocation(sys_org ,'ex5_ODE_org' ,'ex5_OUT_org' ,col.org,simoptions,'uv1',uv1,'uv2',uv2,'method','final');
fprintf('Fault Scenario A - ')
FauA_PC = PoCETsimCollocation(sys_fauA,'ex5_ODE_fauA','ex5_OUT_fauA',col.fauA,simoptions,'uv1',uv1,'uv2',uv2,'method','final');
fprintf('Fault Scenario B - ')
FauB_PC = PoCETsimCollocation(sys_fauB,'ex5_ODE_fauB','ex5_OUT_fauB',col.fauB,simoptions,'uv1',uv1,'uv2',uv2,'method','final');

% get samples
Org_MC_y  = basis_evals_org'*Org_PC.y.pcvals;
FauA_MC_y = basis_evals_fauA'*FauA_PC.y.pcvals;
FauB_MC_y = basis_evals_fauB'*FauB_PC.y.pcvals;

lb = min([Org_MC_y; FauA_MC_y; FauB_MC_y]);
ub = max([Org_MC_y; FauA_MC_y; FauB_MC_y]);

EDGES = linspace(lb,ub,100); % Intervals used to generate the histograms
N_Org_MC  = histc(Org_MC_y,EDGES); % histogram of Model 1 (constructed with MC simulations) at the end point of simulation
N_FauA_MC = histc(FauA_MC_y,EDGES); % histogram of Model 2 (constructed with MC simulations) at the end point of simulation
N_FauB_MC = histc(FauB_MC_y,EDGES); % histogram of Model 2 (constructed with MC simulations) at the end point of simulation

BC_OA = sum(sqrt((N_Org_MC/sum(N_Org_MC)).*(N_FauA_MC/sum(N_FauA_MC))));
BC_OB = sum(sqrt((N_Org_MC/sum(N_Org_MC)).*(N_FauB_MC/sum(N_FauB_MC))));
BC_AB = sum(sqrt((N_FauA_MC/sum(N_FauA_MC)).*(N_FauB_MC/sum(N_FauB_MC))));

obj = -(sqrt(1-BC_OA)+sqrt(1-BC_OB)+sqrt(1-BC_AB)); % Objective is to maximzie the distance between the two distributions
% ceq = 0;

fprintf(repmat('\b',1,3*85))
end