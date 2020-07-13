function PHI_rj = calc_legendre_polynomials(rj,P,qp)
% PoCET internal function. Do not call directly!

 PHI_rj = zeros(P+1,qp);
 
 PHI_rj(1,:) = 1;
 PHI_rj(2,:) = rj;
 for j = 3:P+1
  PHI_rj(j,:) = ((2*j-3)*rj.*PHI_rj(j-1,:)-(j-2)*PHI_rj(j-2,:))/(j-1);
 end
end