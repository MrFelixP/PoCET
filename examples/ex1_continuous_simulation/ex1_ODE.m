function dXdt = ex1_ODE(t,X,PoCETsys,ut,uv)
 M = PoCETsys.coeff_matrices;

 x = X(0*PoCETsys.pce.options.n_phi+1:1*PoCETsys.pce.options.n_phi);
 
 u = piecewise(ut,uv,t);
 
 xx = ckron(x,x);
 
 ddt_x = - M.aa_O2*xx + u*M.b_O0;
 
 dXdt = [ddt_x];
end