%% define MPC value-function optfun
% The MPC objective is to minimize the distance between a desired target
% state x_F (defined in mpc_param) and the mean of the state x at the end 
% of the prediction horizon. To this end, the following steps are taken:
%  (1) simulate the system with the current input sequence uv
%  (2) calculate the trajectory mean of x
%  (3) calculate the value

function value = ex6_optfun(sys,mpc_param,MomMats,ut,uv)
 x_F = mpc_param.x_F;
 Q = mpc_param.Q;
 simopt = mpc_param.simopt;

 % (1) simulate the system with current input sequence
 sim = PoCETsimGalerkin(sys,'ex6_ODE',[],simopt,'ut',ut,'uv',uv);

 % compute mean of x 
 mean = PoCETcalcMoments(sys,MomMats,sim.x.pcvals,[],1);
 
 % get value of current input signal
 value = (Q*(mean(end)-x_F)).^2;
end