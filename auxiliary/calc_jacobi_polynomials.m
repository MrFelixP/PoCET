function PHI_rj = calc_jacobi_polynomials(rj,alf,bet,P,qp)
% PoCET internal function. Do not call directly!

 PHI_rj = zeros(P+1,qp);
 
 PHI_rj(1,:) = 1;
 PHI_rj(2,:) = (alf-bet+(2+alf+bet)*rj)/2;
 for j = 3:P+1
  temp = 2*(j-1)+alf+bet;
  a = 2*(j-1)*(j-1+alf+bet)*(temp-2);
  b = (temp-1)*(alf*alf-bet*bet+temp*(temp-2)*rj);
  c = 2*(j-2+alf)*(j-2+bet)*temp;
  PHI_rj(j,1:qp) = (b.*PHI_rj(j-1,:)-c.*PHI_rj(j-2,:))./a;
 end
end