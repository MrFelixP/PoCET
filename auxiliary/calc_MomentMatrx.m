function [val, phi] = calc_MomentMatrx(opts,PHI_rj,wj,moment_order)
% PoCET internal function. Do not call directly!

 mo = moment_order;
 psimap = get_PSImap(opts.n_xi,opts.order);
 
 if mo == 1
  c_num = 1;
  for k = 1:opts.n_xi
   c_num = c_num*(wj(k,:)*PHI_rj(1,:,k)');
  end
  val = c_num;
  phi = [0 0];
  
 elseif mo == 2
  combis = [1:opts.n_phi; 1:opts.n_phi]';
  c_num = get_num(combis,psimap,wj,PHI_rj,opts.n_xi);
  val = c_num;
  phi = combis;
 
 else % mo > 2
  combis = combs_rep(size(psimap,1),mo);
  nonzero = zeros(size(combis,1),1);
  for i = 1:size(combis,1)
   psi = psimap(combis(i,:),:);
   nonzero(i) = ~iszero(psi);
  end
  combis = combis(nonzero~=0,:);

  c_num = get_num(combis,psimap,wj,PHI_rj,opts.n_xi);
  
  val = zeros(factorial(mo)*size(combis,1),1);
  phi = zeros(factorial(mo)*size(combis,1),mo);
  k = 0;
  
  for i = 1:size(combis,1)
   cperms = uperm(combis(i,:));
   for j = 1:size(cperms,1)
    k = k+1;
    if c_num(i) > opts.quadtol
     val(k) = c_num(i);
     phi(k,:) = combis(i,:);
    end
   end
  end
  
  val(val==0) = [];
  phi(~any(phi,2),:) = [];
 end
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

function c_num = get_num(combis,psimap,wj,PHI_rj,n_xi)
 c_num = ones(size(combis,1),1);
 for j = 1:size(combis,1)
  psi = psimap(combis(j,:),:);
  for k = 1:n_xi
   xi_i = psi(:,k);
   c_num_vals = wj(k,:);
   for l = 1:length(xi_i)
    c_num_vals = c_num_vals.*PHI_rj(xi_i(l)+1,:,k);
   end
   c_num_vals = sum(c_num_vals);
   c_num(j) = c_num(j)*c_num_vals;
  end
 end
end