function out = piecewise(tvec,vvec,t)
% out = piecewise(tvec,vvec,t) 
% returns 
%  0          for             t <  tvec(1)
%  vvec(i)    for  tvec(i) >= t >  tvec(i+1)
%  vvec(end)  for             t >= tvec(end)
% 
%  'tvec' is a vector containing the time-points, at which the 
%  input value changes.
%
%  'vvec' is a vector containing the values, the function takes 
%  at the specified timepoints.
%
%  't' is the timepoint, at which the function is evaluated.
%
% In simulations, specify tvec and vvec as input parameters and
% piecewise(tvec,vvec,t) as the input function, where t will be 
% interpreted as the simulation time.

 dv = diff([0 vvec]);
 out = dv*(ones(size(t'))*tvec<=t'*ones(size(tvec)))';

end