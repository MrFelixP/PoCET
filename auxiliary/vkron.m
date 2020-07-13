function X = vkron(A,B)
% X = vkron(A,B) 
% returns the Kronecker product of two dense vectors A and B of dimensions 
% M x N with either M = 1 or N = 1, and K x L with either K = 1 or L = 1, 
% respectively. The result is an M*K x N*L matrix in which the (m,n)-th 
% element is defined as A(m)*B(n). 
 
 [M,N] = size(A);
 [K,L] = size(B);
 
 if N == 1 && L == 1 % A and B are column vectors
  X = B*A';
  X = X(:);
 elseif M == 1 && K == 1 % A and B are row vectors
  X = B'*A;
  X = (X(:))';
 elseif N == 1 && K == 1 % A is column vector and B is row vector
  X = A*B;
 elseif M == 1 && L == 1 % A is row vector and B is column vector
  X = B*A;
 end
 
end