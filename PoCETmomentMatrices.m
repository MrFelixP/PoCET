function MomMats = PoCETmomentMatrices(PoCETsys, M)
% MomMats = PoCETmomentMatrices(PoCETsys, M) 
% returns a struct array of size (M-1)-by-1 containing the matrices
% required to calculate the second to M-th moments of a variable from its
% PCE-coefficients. 
% NOTE: computation times increase exponentially with M!

% check if enough quadrature points
if PoCETsys.pce.options.hpo < M
 PoCETsys.pce.options.hpo = M;
 PoCETsys.quadrature = PoCETquadRules(PoCETsys);
end

% calculate matrices
MomMats = struct('val',cell(M,1),'phi',cell(M,1));
for i = 1:M
 [MomMats(i).val, MomMats(i).phi] = ...
  calc_MomentMatrx(PoCETsys.pce.options, ...
                   PoCETsys.quadrature.PHI_rj,...
                   PoCETsys.quadrature.wj,i);
end
end