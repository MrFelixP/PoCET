function X = rkron(A,B)
% X = fastkron(A,B) 
% returns the Kronecker product of two dense row vectors A and B of 
% dimensions 1 x M and 1 x N, respectively. The result is an 1 x M*N row
% vector in which the ((n-1)*m+n)-th element is defined as A(m)*B(n). 

 X = B'*A; 
 X = (X(:))';
end