function [mats raw] = PoCETcoeffMatrices(PoCETsys)
% mats = PoCETcoeffMatrices(PoCETsys) 
% returns the coefficient matrices for the system PoCETsys which has been
% composed using PoCETcompose. 
%
% [mats raw] = PoCETcoeffMatrices(PoCETsys) 
% also returns the 'raw' coefficient matrices, i. e. the coefficient
% matrices without including the parameters' values.

% parameters = {PoCETsys.parameters(:).name};
if ~isempty(PoCETsys.parameters(1).name)
 parameters = {PoCETsys.parameters.name};
else
 parameters = {};
end
params = {};
sortparams = {};

%% analyze ode for each state
for i=1:numel(PoCETsys.states)
 analyze_terms(PoCETsys.ode(i).terms);
end

for i=1:numel(PoCETsys.outputs)
 analyze_terms(PoCETsys.out(i).terms);
end

%% calculate coefficient matrix
[~, uidx] = unique(sortparams);
structvals = params(1:2,uidx);
rawnames = params(1,uidx);
rawvals = params(2,uidx);
orders = params(4,uidx);
pars = params(3,uidx);

for i = 1:numel(uidx)
 [structvals{2,i},rawvals{i}]  = ...
  calc_MultiCoefMatrx(PoCETsys,structvals{2,i},pars{i},orders{i},...
  PoCETsys.quadrature.PHI_rj,PoCETsys.quadrature.wj);
end
structvals = structvals(:)';
mats = struct(structvals{:});

raw = cell2struct(rawvals,rawnames,2);
 
function analyze_terms(terms)
 for j=1:numel(terms)
  % find only valid parameters
  [~,zIDX] = find(terms(j).param_index);
  tmpIDX = terms(j).param_index(zIDX);
  tmpDEG = terms(j).param_degrees(zIDX);

  % find uique parameter combinations and get degrees of parameters
  [~,sIDX] = sort(tmpIDX);
  tmpIDX = tmpIDX(sIDX);
  tmpDEG = tmpDEG(sIDX);

  % get order of state polynomial
  tmpORD = sum(terms(j).state_degrees);

  tmpPAR = [];
  tmpXI  = [];

  % build names for parameter combinations
  for k=1:numel(tmpIDX)
   tmpPAR = [tmpPAR repmat(tmpIDX(k),1,tmpDEG(k))]; 
  end

  % check if parameters are associated with random variables
  for k=1:numel(tmpPAR)
   xi_num = find(ismember(PoCETsys.pce.vars,parameters(tmpPAR(k))));
   if ~isempty(xi_num)
    tmpXI = [tmpXI xi_num];
   else
    tmpXI = [tmpXI 0];
   end
  end

  if ~isempty(tmpPAR)
   params = [params {[parameters{tmpPAR} '_O' num2str(tmpORD)]; tmpXI; tmpPAR; tmpORD}];
   sortparams = [sortparams {[parameters{tmpPAR} ' - ' num2str(tmpXI) ' - ' num2str(tmpORD)]}];
  else
   params = [params {['one_O' num2str(tmpORD)]; []; tmpPAR; tmpORD}];
   sortparams = [sortparams {['0 - 0 - ' num2str(tmpORD)]}];
  end
 end
end

end

%  for j=1:numel(PoCETsys.ode(i).terms)
%   % find only valid parameters
%   [~,zIDX] = find(PoCETsys.ode(i).terms(j).param_index);
%   tmpIDX = PoCETsys.ode(i).terms(j).param_index(zIDX);
%   tmpDEG = PoCETsys.ode(i).terms(j).param_degrees(zIDX);
% 
%   % find uique parameter combinations and get degrees of parameters
%   [~,sIDX] = sort(tmpIDX);
%   tmpIDX = tmpIDX(sIDX);
%   tmpDEG = tmpDEG(sIDX);
% 
%   % get order of state polynomial
%   tmpORD = sum(PoCETsys.ode(i).terms(j).state_degrees);
% 
%   tmpPAR = [];
%   tmpXI  = [];
% 
%   % build names for parameter combinations
%   for k=1:numel(tmpIDX)
%    tmpPAR = [tmpPAR repmat(tmpIDX(k),1,tmpDEG(k))]; 
%   end
% 
%   % check if parameters are associated with random variables
%   for k=1:numel(tmpPAR)
%    xi_num = find(ismember(PoCETsys.pce.vars,parameters(tmpPAR(k))));
%    if ~isempty(xi_num)
%     tmpXI = [tmpXI xi_num];
%    else
%     tmpXI = [tmpXI 0];
%    end
%   end
% 
%   if ~isempty(tmpPAR)
%    params = [params {[parameters{tmpPAR} '_O' num2str(tmpORD)]; tmpXI; tmpIDX; tmpORD}];
%    sortparams = [sortparams {[parameters{tmpPAR} ' - ' num2str(tmpXI) ' - ' num2str(tmpORD)]}];
%   else
%    params = [params {['one_O' num2str(tmpORD)]; []; tmpIDX; tmpORD}];
%    sortparams = [sortparams {['0 - 0 - ' num2str(tmpORD)]}];
%   end
%  end