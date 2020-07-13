function basis = PoCETbasis(PoCETsys)
% basis = PoCETbasis(PoCETsys)
% returns an array of function handles containing the basis functions of
% the polynomial chaos expanded system PoCETsys.

opts = PoCETsys.pce.options;
state_names = {PoCETsys.states.name};
if ~isempty(PoCETsys.parameters(1).name)
 param_names = {PoCETsys.parameters.name};
else
 param_names = {};
end
basis = sym(1);
basis2 = [];

for i = 1:opts.n_xi
 x = sym(['xi_' num2str(i)]);
 statenum = find(ismember(state_names,PoCETsys.pce.vars{i}));
 paramnum = find(ismember(param_names,PoCETsys.pce.vars{i}));
 if ~isempty(statenum)
  curvar = PoCETsys.states(statenum);
 else
  curvar = PoCETsys.parameters(paramnum);
 end
 
 % construct orthogonal basis
 if strcmpi(curvar.dist,'normal')
  % normal distribution
  tmpBASE = [sym(1) x];
  for j = 3:opts.order+1
   tmpBASE(j) = x*tmpBASE(j-1) - (j-2)*tmpBASE(j-2);
  end
 elseif strcmpi(curvar.dist,'beta') || strcmpi(curvar.dist,'beta4')
  % beta2 and beta4 distributions
  alf = curvar.data(2)-1;
  bet = curvar.data(1)-1;
  tmpBASE = [sym(1) (alf-bet+(2+alf+bet)*x)/2];
  for j = 3:opts.order+1
   temp = 2*(j-1)+alf+bet;
   a = 2*(j-1)*(j-1+alf+bet)*(temp-2);
   b = (temp-1)*(alf*alf-bet*bet+temp*(temp-2)*x);
   c = 2*(j-2+alf)*(j-2+bet)*temp;
   tmpBASE(j) = (b*tmpBASE(j-1)-c*tmpBASE(j-2))/a;
  end
 elseif strcmpi(curvar.dist,'uniform')
  % uniform distribution
  tmpBASE = [sym(1) x];
  for j = 3:opts.order+1
   tmpBASE(j) = ((2*j-3)*x*tmpBASE(j-1) - (j-2)*tmpBASE(j-2))/(j-1);
  end
 elseif strcmpi(curvar.dist,'matrix')
  % everything is fine
 elseif strcmpi(curvar.dist,'pce')
  % everything is fine
 else
  % error
  error('Unknown distribution ''%s''.',curvar.dist);
 end
 
 %     basis = kron(basis,tmpBASE);
 basis2 = [basis2,reshape(tmpBASE,numel(tmpBASE),1)];
end

PSImap = get_PSImap(PoCETsys.pce.options.n_xi,PoCETsys.pce.options.order);
for i = 1:size(PSImap,1),
 basis(i) = sym(1);
 for j = 1:opts.n_xi,
  basis(i) = basis(i)*basis2(PSImap(i,j)+1,j);
 end
 basis(i) = expand(basis(i));
end

% % basis = unique(expand(basis));
% basis = expand(basis);
% degrees = zeros(size(basis));
%
% for i=1:numel(basis)
%     degrees(i) = feval(symengine,'degree',basis(i));
% end
%
% basis = basis(degrees<=opts.order);
%
%% convert into anonymous function
basisStr = '[';
for i = 1:numel(basis),
 basisStr = [basisStr,char(basis(i)),';'];
end
varsStr = '';
for i = 1:numel(PoCETsys.randomVariables),
 varsStr = [varsStr,'xi_',num2str(i),','];
end
if isempty(varsStr)
 varsStr = 'xi ';
end
eval(['basis = ',['@(',varsStr(1:end-1),')',basisStr(1:end-1),'];']]);

end