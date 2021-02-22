function ablu = calcBeta4(mu_xi) 
% ABLU = calcBeta4(mu) 
% returns an 1-by-4 array containing the beta-distibution-parameters alpha
% and beta based on the central moments mu of some random variable xi as
% well as the lower and upper bounds of its support using the method of
% moments for the four-parameter beta distribution.
% Scource:
% http://en.wikipedia.org/wiki/Beta_distribution#Four_unknown_parameters
%
% The input mu is an array containing the first four central moments mu 
% of a random variable x i, i. e. its mean, variance, skewness, and excess 
% curtosis.

 ablu = [1 1 0 1];    % preallocate output
 
% mu_xi(4) = -abs(mu_xi(4)); 
% The method is conditional on 
%   mu_xi(3)^2-2 < mu_xi(4) < 1.5*mu_xi(3)^2 
% and since |mu_xi(3)| is usually very small positive values of 
% mu_xi(4) are very likely to violate this condition. Forcing it to be
% negative is more a brute force way to ensure the procedure doesn't
% collapse, although mathematically it is... highly quenstionable.

% The rest of the script are simply the formulas given in the scource above
 v = 3*((mu_xi(4)-mu_xi(3)^2+2)/(1.5*mu_xi(3)^2-mu_xi(4)));
 if abs(mu_xi(3)) < 1e-5 % skewness is numerically 0
  ablu(1) = v/2;
  ablu(2) = ablu(1);
 else
     assert(1+16*(v+1)/((v+2)^2*mu_xi(3)^2) >= 0, 'Error fitting Beta.');
     alfbet = 1/sqrt(1+16*(v+1)/((v+2)^2*mu_xi(3)^2));
     alfbet = [v/2*(1+alfbet) v/2*(1-alfbet)];
  if mu_xi(3) < 0
   ablu(1) = max(alfbet);
   ablu(2) = min(alfbet);
  else
   ablu(1) = min(alfbet);
   ablu(2) = max(alfbet);
  end
 end
 assert( (v+2)^2*mu_xi(3)^2+16*(v+1) >= 0, 'Error fitting Beta.');
 lu = 0.5*sqrt(mu_xi(2))*sqrt((v+2)^2*mu_xi(3)^2+16*(v+1));
 ablu(3) = mu_xi(1)-ablu(1)*lu/v;
 ablu(4) = lu+ablu(3);
 
% correction af alpha and beta to match definition of jacobi polynomials
 ablu(1) = ablu(1) - 1;
 ablu(2) = ablu(2) - 1;
 end