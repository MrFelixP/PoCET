function dXdt = ex1_MCODE(t,X,PAR,ut,uv)

 x = X(1);
 
 u = piecewise(ut,uv,t);
 
 ddt_x = -PAR.a^2*x^2 + PAR.b*u;

 dXdt = [ddt_x];
end