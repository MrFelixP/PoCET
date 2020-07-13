function [X W] = gauss_jacobi_quadrature(N,alf,bet)
% [X W] = gauss_jacobi_quadrature(N)
%
% Find the Gauss-Jacobi abscissae and weights.
%
% Arguments:
%  N - The number of abscissae and weights to return.
%
% Return Values:
%  X - A column vector containing the abscissae.
%  W - A column vector containing the corresponding weights.
%
% Gauss-Hermite quadrature approximates definite integrals of the form
%
%     \\int^{-\\infty}_{\\infty} dx W(x) f(x)
%
% where
%
%     W(x) = \\exp( - x^2 )
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
%	AUTHOR    = { Press, William  H.   and 
%                 Flannery, Brian  P.  and
%                 Teukolsky, Saul  A.  and
%                 Vetterling, William  T. },
%   ISBN      = {0521431085},
%	MONTH     = {October},
%	PUBLISHER = {{Cambridge University Press}},
%	TITLE     = {Numerical Recipes in C : The Art of Scientific Computing},
%	YEAR      = {1992}
% }
%

EPS = 3.0e-14; % precision
MAXIT = 10; % maximum number of loops
alfbet = alf+bet;

% allocate the return values
X = zeros([N 1]);
W = zeros([N 1]);

for i=1:N
 % good guesses at initial values for specific roots
 if i == 1
  an = alf/N;
  bn = bet/N;
  r1 = (1+alf)*(2.78/(4+N*N)+0.768*an/N);
  r2 = 1+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
  z  = 1-r1/r2;
 elseif i == 2
  r1 = (4.1+alf)/((1+alf)*(1+0.156*alf));
  r2 = 1+0.06*(N-8)*(1+0.12*alf)/N;
  r3 = 1+0.012*bet*(1+0.25*abs(alf))/N;
  z  = z-(1-z)*r1*r2*r3;
 elseif i == 3
  r1 = (1.67+0.28*alf)/(1+0.37*alf);
  r2 = 1+0.22*(N-8)/N;
  r3 = 1+8*bet/((6.28+bet)*N*N);
  z  = z-(X(1)-z)*r1*r2*r3;
 elseif i == N-1
  r1 = (1+0.235*bet)/(0.766+0.119*bet);
  r2 = 1/(1+0.639*(N-4)/(1+0.71*(N-4)));
  r3 = 1/(1+20*alf/((7.5+alf)*N*N));
  z  = z+(z-X(N-3))*r1*r2*r3;
 elseif i == N
  r1 = (1+0.37*bet)/(1.67+0.28*bet);
  r2 = 1/(1+0.22*(N-8)/N);
  r3 = 1/(1+8*alf/((6.28+alf)*N*N));
  z  = z+(z-X(N-2))*r1*r2*r3;
 else
  z  = 3*X(i-1)-3*X(i-2)+X(i-3);
 end

 for iter = 1:MAXIT+1
  temp = 2+alfbet;
  p1 = (alf-bet+temp*z)/2;
  p2 = 1;
  for j=2:N
   p3 = p2;
   p2 = p1;
   temp = 2*j+alfbet;
   a = 2*j*(j+alfbet)*(temp-2);
   b = (temp-1)*(alf*alf-bet*bet+temp*(temp-2)*z);
   c = 2*(j-1+alf)*(j-1+bet)*temp;
   p1 = (b*p2-c*p3)/a;
  end

  % the derivative
  pp = (N*(alf-bet-temp*z)*p1+2*(N+alf)*(N+bet)*p2)/(temp*(1-z*z));

  % newton step
  z1 = z;
  z  = z1-p1/pp;

  if abs(z-z1) <= EPS
   break;
  end    
 end

 if iter == MAXIT+1
  fprintf('Too many iterations in hermquad.\n');
 end

 X(i) = z;
 W(i) = exp(gammaln(alf+N)+gammaln(bet+N)-gammaln(N+1)-gammaln(N+alfbet+1))*temp*2^alfbet/(pp*p2);
end
