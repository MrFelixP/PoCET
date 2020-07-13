function [X W] = gauss_legendre_quadrature(N)
% [X W] = gauss_legendre_quadrature(N)
%
% Find the Gauss-Legendre abscissae and weights.
%
% Arguments:
%  N - The number of abscissae and weights to return.
%
% Return Values:
%  X - A column vector containing the abscissae.
%  W - A column vector containing the corresponding weights.
%
% Gauss-Legendre quadrature approximates definite integrals of the form
%
%     \\int^{-\\infty}_{\\infty} dx W(x) f(x)
%
% where
%
%     W(x) = 1
%
% with the sum
%
%     \\sum_{n=1}^{N} w_{n} f(x_{n}).
%
% This function returns the set of abscissae and weights
%
%     {x_{n}, w_{n}}^{N}_{n=1}
%
% for performing this calculation given N, the number of abscissae.
% These abscissae correspond to the zeros of the Nth Hermite
% polynomial.  It can be shown that such integration is exact when f(x)
% is a polynomial of maximum order 2N-1.
%
% The procedure in this calculation is taken more or less directly from
%
% @BOOK{ press-etal-1992a,
%  AUTHOR   = { Press, William  H.   and 
%               Flannery, Brian  P.  and
%               Teukolsky, Saul  A.  and
%               Vetterling, William  T. },
%   ISBN      = {0521431085},
%	MONTH     = {October},
%	PUBLISHER = {{Cambridge University Press}},
%	TITLE     = {Numerical Recipes in C : The Art of Scientific Computing},
%	YEAR      = {1992}
% }
%

EPS = 3.0e-14; % precision
PIAP = 3.141592654; % sufficient approximation of pi
MAXIT = 10; % maximum number of loops

% allocate the return values
X = zeros([N 1]);
W = zeros([N 1]);

for i=1:(N+1)/2
  
  % good guesses at initial values for specific roots
  z = cos(PIAP*(i-0.25)/(N+0.5));
  
  for iter=1:MAXIT+1
    p1 = 1;
    p2 = 0.0;
    
    for j=1:N
      p3 = p2;
      p2 = p1;
      p1 = ((2*j-1)*z*p2-(j-1)*p3)/j;
    end
    
    % the derivative
    pp = N*(z*p1-p2)/(z*z-1);
    
    % newton step
    z1 = z;
    z  = z1 - p1/pp;
    
    if abs(z-z1) <= EPS
      break;
    end    
  end
  
  if iter == MAXIT+1
    fprintf('Too many iterations in hermquad.\n');
  end
  
  X(i)     = z;
  X(N+1-i) = -z;
  W(i)     = 2/((1-z*z)*pp*pp);
  W(N+1-i) = W(i);

end
