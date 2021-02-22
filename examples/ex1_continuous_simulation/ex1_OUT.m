function fout = ex1_OUT(t,X,PoCETsys,ut,uv)
 M = PoCETsys.coeff_matrices;

 x = X(0*PoCETsys.pce.options.n_phi+1:1*PoCETsys.pce.options.n_phi);
 
 u = piecewise(ut,uv,t);
 
 xx = ckron(x,x);
 
 y = 3*M.one_O2*xx;
 
 fout = [y];
end