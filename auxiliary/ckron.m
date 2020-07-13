function X = ckron(A,B)
% X = ckron(A,B) 
% returns the Kronecker product of two dense column vectors A and B of 
% dimensions M x 1 and N x 1, respectively. The result is an M*N x 1 column
% vector in which the ((n-1)*m+n)-th element is defined as A(m)*B(n). 

 X = B*A';
 X = X(:);
end