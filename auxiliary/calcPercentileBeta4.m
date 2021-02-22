function XU = calcPercentileBeta4(ABLU,PU)
% XU = calcPercentileBeta4(ABLU,PU)
% returns the position XU of the PU'th percentile of a beta-4 distributed 
% variable X, i.e. P(X>XU) < 1-PU, where the distribution of X is defined
% via the parameters ABLU with
% ABLU = [ALPHA, BETA, LOWER_BOUND, UPPER_BOUND].

 nvals = size(ABLU,3);
 dx_steps = 1e3;
 XU = zeros(1,nvals);
 for j = 1:nvals
  low = ABLU(1,3,j);
  upp = ABLU(1,4,j);
  dx = linspace(low,upp,dx_steps);
  dy = beta4(ABLU(1,:,j),dx);
  for k = 2:dx_steps
   q_k = trapz(dx(1:k),dy(1:k));
   if q_k > PU
    XU(j) = dx(k-1);
    break;
   end
  end
 end
end