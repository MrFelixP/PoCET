function quadrules = PoCETquadRules(PoCETsys)
% q = PoCETquadRules(PoCETsys) 
% returns a struct containing the roots, weights, normalized weights,
% polynomials, and normalized polynomials required for integration via
% gaussian quadrature. All polynomials are already evaluated at the roots.

opts = PoCETsys.pce.options;
state_names = {PoCETsys.states(:).name};
if ~isempty(PoCETsys.parameters(1).name)
 param_names = {PoCETsys.parameters(:).name};
else
 param_names = {};
end

qp = ceil(opts.hpo*opts.order/2 + 0.5);

% preallocation
rj = zeros(opts.n_xi,qp);
wj = zeros(opts.n_xi,qp);
PHI_rj = zeros(opts.order+1,qp,opts.n_xi);

% calculate gaussian quadrature roots & weights and evaluate polynomials
for i = 1:opts.n_xi
 statenum = find(ismember(state_names,PoCETsys.pce.vars{i}));
 paramnum = find(ismember(param_names,PoCETsys.pce.vars{i}));
 if ~isempty(statenum)
  curvar = PoCETsys.states(statenum);
 else
  curvar = PoCETsys.parameters(paramnum);
 end

 if strcmp(curvar.dist,'normal')
  [rj(i,1:qp), wj(i,1:qp)] = gauss_hermite_quadrature(qp);
  wj(i,:) = wj(i,:)/sqrt(2*pi);
  PHI_rj(:,:,i) = calc_hermite_polynomials(rj(i,:),opts.order,qp);
 elseif strcmp(curvar.dist,'uniform')
  [rj(i,1:qp), wj(i,1:qp)] = gauss_legendre_quadrature(qp);
  wj(i,:) = wj(i,:)/2;
  PHI_rj(:,:,i) = calc_legendre_polynomials(rj(i,:),opts.order,qp);
 elseif strcmp(curvar.dist,'beta') || strcmp(curvar.dist,'beta4')
  alf = curvar.data(1);
  bet = curvar.data(2);
  [rj(i,1:qp), wj(i,1:qp)] = gauss_jacobi_quadrature(qp,bet-1,alf-1);
  wj(i,:) = wj(i,:)/(2^(alf+bet-1)*beta(alf,bet));
  PHI_rj(:,:,i) = calc_jacobi_polynomials(rj(i,:),bet-1,alf-1,opts.order,qp);
 end
end

quadrules.rj = rj;
quadrules.wj = wj;
quadrules.PHI_rj = PHI_rj;
end