function dH = calcMDCbeta(ABLU1,ABLU2)
% dH = calcMDCbeta(ABLU1,ABLU2)
% returns the value of the integral over the product of two four-parameter 
% beta distributions defined by ABLU1 and ABLU2 with
% ABLU = [ALPHA, BETA, LOWER_BOUND, UPPER_BOUND].
% NOTE: Since    dH = 0  <=>  the distributions don't overlap at all
%       and      dH = 1  <=>  the distributions are identical
%       this value can be used as a measure of the overlap or similarity of
%       both distrubtions.

 nvals = size(ABLU1,3);
 dH_vals = ones(1,nvals);
 for j = 1:nvals
  low1 = ABLU1(1,3,j);
  upp1 = ABLU1(1,4,j);
  low2 = ABLU2(1,3,j);
  upp2 = ABLU2(1,4,j);
  if upp1 > low2 && upp2 > low1
   dx = linspace(max(low1,low2),min(upp1,upp2),1e3);
   dy1 = beta4(ABLU1(1,:,j),dx);
   dy2 = beta4(ABLU2(1,:,j),dx);
   dH_vals(j) = trapz(dx,sqrt(dy1.*dy2));
  else
   dH_vals(j) = 0;
  end
 end
%  dH = sum(dH_vals)*(1/nvals);
%  dH = prod(dH_vals)^(1/(nvals-1));
 dH = sum(dH_vals);
end