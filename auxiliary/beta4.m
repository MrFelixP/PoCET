function Y = beta4(ABLU,X)
% Y = beta4(ABLU, X)
% returns the values of the four-parameter beta distribution defined by
% ABLU = [ALPHA, BETA, LOWER_BOUND, UPPER_BOUND] evaluated at points X. 

 alf = ABLU(1);
 bet = ABLU(2);
 low = ABLU(3);
 upp = ABLU(4);
 Y = (X-low).^(alf).*(upp-X).^(bet)./((upp-low).^(alf+bet+1)*beta(alf+1,bet+1));
 Y(X<low) = 0;
 Y(X>upp) = 0;
 Y(isinf(Y)) = max(Y(isfinite(Y)));
end