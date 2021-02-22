%% define nonlinear (chance-)constraints
% The constraint to be enforced is P(x>x_ub)<1-p_ub, i.e. the probability
% of NOT violating an upper bound x_ub for x should not exeed p_ub at any 
% time instant. Note that the constraints are enforced for each time 
% instant individually, i.e. the constraint is, more accurately, given by
% P(max_t{x(t)>x_ub})<1-p_ub.  The procedure here is as follows:
%  (1) simulate the system with the current input signal.
%  (2) calculate the trajectories of the first 4 central moments from the pce
%  (3) fit a 4-parameter beta distribution at each time instance
%  (4) calculate the postitions xpc(t) of the p_ub percentiles
%  (5) enforce contraints c = xpc(t)-p_ub <= 0 for all t

function [c,ceq] = ex6_confun(sys,mpc_param,MomMats,ut,uv)
 x_ub = mpc_param.x_ub;
 p_ub = mpc_param.p_ub;
 simopt = mpc_param.simopt;

 % (1) simulate the system with current input sequence
 sim = PoCETsimGalerkin(sys,'ex6_ODE',[],simopt,'ut',ut,'uv',uv);

 % (2) compute moments of x 
 moments = PoCETcalcMoments(sys,MomMats,sim.x.pcvals(:,2:end));
 
 % (3) fit a Beta-4-distribution to moments
 ablu = zeros(1,4,size(moments,2));
 for i = 1:size(moments,2)
  ablu(1,:,i) = calcBeta4(moments(:,i));
 end
 
 % (4) calculate positions of p_ub upper percentiles
 xpc = calcPercentileBeta4(ablu,p_ub);
 
 % (5) constraints
 ceq = 0;
 c = xpc - x_ub;
end