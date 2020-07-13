function C = calc_CoefMatrx(opts,project_on,mult_order,PHI_rj,wj)
% PoCET internal function. Do not call directly!

psimap = get_PSImap(opts.n_xi,opts.order); % list of all multi-polynomials
combis = combs_rep(size(psimap,1),mult_order+1); % list of unique polynomial-multiplications
nonzero = zeros(size(combis,1),1);

for i = 1:size(combis,1)
 psi = psimap([project_on combis(i,:)],:); % get polynomials to be integrated
 nonzero(i) = ~iszero(psi); % 
end
combis = combis(nonzero~=0,:);

c_num = ones(size(combis,1),1);
for j = 1:size(combis,1)
 psi = psimap([project_on combis(j,:)],:);
 for k = 1:opts.n_xi     % cycle through all random variables (RVs)
  xi_i = psi(:,k);       % get polynomial orders of current RV
  c_num_vals = wj(k,:);  % initialize gauss quadrature by setting weights
  for l = 1:length(xi_i) % cycle through polynomial orders of current RV
   c_num_vals = c_num_vals.*PHI_rj(xi_i(l)+1,:,k); % 
  end
  c_num_vals = sum(c_num_vals);
  c_num(j) = c_num(j)*c_num_vals;
 end
end

c_den = ones(opts.n_phi,1); % denominator <\PSI_{i_j}\PSI_{i_j}> constant in each row
for i = 1:opts.n_phi
 for j = 1:opts.n_xi
  xi_j = psimap(i,j); % get order of xi_j in i'th base polynomial
  den_vals = PHI_rj(xi_j+1,:,j).^2*wj(j,:)';
  c_den(i) = c_den(i) * den_vals;
 end
end

C = sparse(opts.n_phi,opts.n_phi^mult_order);
idxvec = opts.n_phi.^[mult_order-1:-1:0];
for i = 1:size(combis,1)
 cperms = uperm(combis(i,:));
 for j = 1:size(cperms,1)
  row = cperms(j,end);
  col = (cperms(j,1:end-1)-1)*idxvec' + 1;
  c_ij = c_num(i)/c_den(row);
  if c_ij > opts.quadtol
   C(row,col) = c_ij;
  end
 end
end
% C = sparse(C);
end

function CR = combs_rep(N,K)
% Subfunction multichoose:  combinations with replacement.
% cr = @(N,K) prod((N):(N+K-1))/(prod(1:K)); Number of rows.

M = double(N);  % Single will give us trouble on indexing.
WV = ones(1,K,class(N));  % This is the working vector.
mch = prod((M:(M+K-1)) ./ (1:K));  % Pre-allocation.
CR = ones(round(mch),K,class(N));

for ii = 2:mch
 if WV(K) == N
  cnt = K-1;  % Work backwards in WV.

  while WV(cnt) == N
   cnt = cnt-1;  % Work backwards in WV.
  end

  WV(cnt:K) = WV(cnt) + 1;  % Fill forward.
 else
  WV(K) = WV(K)+1;   % Keep working in this group.
 end

 CR(ii,:) = WV;
end
end

function p = uperm(a)
[u, ~, J] = unique(a);
p = u(up(J, length(a)));
end % uperm

function p = up(J, n)
ktab = histc(J,1:max(J));
l = n;
p = zeros(1, n);
s = 1;
for i=1:length(ktab)
    k = ktab(i);
    c = nchoosek(1:l, k);
    m = size(c,1);
    [t, ~] = find(~p.');
    t = reshape(t, [], s);
    c = t(c,:)';
    s = s*m;
    r = repmat((1:s)',[1 k]);
    q = accumarray([r(:) c(:)], i, [s n]);
    p = repmat(p, [m 1]) + q;
    l = l - k;
end
end