function fout = ex1_MCOUT(t,X_in,PAR,ut,uv)

 x = X_in(1,:);
 
 u = piecewise(ut,uv,t);
 
 y = 3*x.^2;
 
 fout = [y];
end