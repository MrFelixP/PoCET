function is0 = iszero(PSI)
% PoCET internal function. Do not call directly!

N_xi = size(PSI,2);
is0 = 0;
for i=1:N_xi
 pows_i = PSI(:,i);
 pows_i = pows_i(pows_i~=0);
 
 % if only one polynomial contains xi_i the integral will be zero
 if numel(pows_i) == 1 
  is0 = 1; return;
  
 % if there are exactly two polynomials containing xi_i but with different
 % powers the integral will be zero
 elseif numel(pows_i) == 2 && pows_i(1) ~= pows_i(2) %diff(pows_i)~=0
  is0 = 1; return;
 
 % if the total power of xi_i is odd the integral will be zero
 elseif mod(sum(pows_i),2)
  is0 = 1; return;
 end
end
end

%% Explanation for developers 
% iszero(PHI) checks if a projection of the input vectors is zero. 
% PSI is expected to be a matrix whose rows contain the polynomial orders 
% for each polynomial. 
% If     <PSI(1,:),PSI(2,:),...,PSI(n,:)> == 0     iszero(PSI) returns 1 
% If     <PSI(1,:),PSI(2,:),...,PSI(n,:)> ~= 0     iszero(PSI) returns 0
%
% These vectors PSI{n} have to be vectors of integers of the same length 
% and contain the orders of their corresponding univariate polynomials.
%
% Example: For two polynomials
%   PSI_1(xi_1,xi_2,xi_3) = PHI_1(xi_1)*PHI_3(xi_2)*PHI_0(xi_3)
%   PSI_2(xi_1,xi_2,xi_3) = PHI_0(xi_1)*PHI_2(xi_2)*PHI_1(xi_3)
% the expected input matrix is
%   PSI = [1 3 0 ; 0 2 1]