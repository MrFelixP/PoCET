function samples = PoCETsample(PoCETsys, type, varargin)
% samples = PoCETsample(PoCETsys, 'variables', N)
% returns a struct containing N samples of all initial conditions and 
% parameters in the system specified by PoCETsys, which is a struct created
% with PoCETcompose.
%
% samples = PoCETsample(PoCETsys, 'basis', N)
% returns a N_PHI-by-N matrix containing the values of the polynomial chaos 
% basis functions evaluated at N samples of the random variables of the
% system PoCETsys. N_PHI is the number of terms used for the PCE.
%
% samples = PoCETsample(PoCETsys, 'PCE', N, pcvals)
% returns a N-by-SIZE(pcvals,2) matrix containing samples that are
% calculated from the polynomial chaos coefficients 'pcvals'.
%
% samples = PoCETsample(PoCETsys, 'PCE', [], pcvals, basis)
% returns a SIZE(basis,2)-by-SIZE(pcvals,2) matrix containing samples that
% are calculated from the polynomial chaos coefficients 'pcvals' and the
% evaluated basis functions 'basis'.

%% init and deal with input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin >= 2 && strcmpi(type,'variables'),
 if nargin == 2,
  samples = sampleRandomVariables(PoCETsys,1);
 else
  samples = sampleRandomVariables(PoCETsys,varargin{1});
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif nargin >= 2 && strcmpi(type,'basis'),
 if nargin == 2,
  tmpSamples = sampleStandardRandomVariables(PoCETsys,1);
%   tmpSamples = [tmpSamples.states;tmpSamples.parameters];
  tmpSamples = tmpSamples.xi;
 else
  tmpSamples = sampleStandardRandomVariables(PoCETsys,varargin{1});
%   tmpSamples = [tmpSamples.states;tmpSamples.parameters];
  tmpSamples = tmpSamples.xi;
 end
 tmpSamples2 = zeros(PoCETsys.pce.options.n_phi,size(tmpSamples,2));
 for i = 1:size(tmpSamples,2),
  tmp = arrayfun(@(x){x},tmpSamples(:,i));
%   keyboard
  if isempty(tmp)
   tmp = {1};
  end
  tmpSamples2(:,i) = PoCETsys.basis(tmp{:});
 end
 samples = tmpSamples2;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif nargin >= 2 && strcmpi(type,'PCE'),
 if nargin == 4,
  basisSamples = PoCETsample(PoCETsys,'basis',varargin{1});
 else
  basisSamples = varargin{3};
 end
 coeffVec = varargin{2};
 samples = basisSamples'*coeffVec;
else
 error('Wrong inputs.');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function samples = sampleRandomVariables(PoCETsys, numSamples)

% sample variables
fprintf('Sampling random variables:   0.0%% finished\n');
 
%  show progress
samples = struct('states',zeros(numel(PoCETsys.states),numSamples),...
                 'parameters',zeros(numel(PoCETsys.parameters),numSamples),...
                 'xi',zeros(numel(PoCETsys.pce.vars),numSamples));
loops = (numel(PoCETsys.states)+numel(PoCETsys.parameters))*numSamples;
k=0;
for vv = {'states','parameters'},
 vars = vv{1};
 for i=1:numel(PoCETsys.(vars)),
  for j=1:numSamples,
   k=k+1;
   if floor(mod(k,loops/20)) == 0,
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6s%% finished',sprintf('%3.1f',k/loops*100));
   end
   switch lower(PoCETsys.(vars)(i).dist),
    case 'normal',
     samples.(vars)(i,j) = PoCETsys.(vars)(i).data(1) + PoCETsys.(vars)(i).data(2)*randn;
    case 'uniform',
     samples.(vars)(i,j) = PoCETsys.(vars)(i).data(1) + (PoCETsys.(vars)(i).data(2) - PoCETsys.(vars)(i).data(1))*rand;
    case 'beta',
     samples.(vars)(i,j) = betarnd(PoCETsys.(vars)(i).data(1),PoCETsys.(vars)(i).data(2));
    case 'beta4',
     samples.(vars)(i,j) = PoCETsys.(vars)(i).data(3) + (PoCETsys.(vars)(i).data(4) - PoCETsys.(vars)(i).data(3))*betarnd(PoCETsys.(vars)(i).data(1),PoCETsys.(vars)(i).data(2));
    case {'none','dirac'}
     samples.(vars)(i,j) = PoCETsys.(vars)(i).data(1);
    otherwise
     error('Wrong input.');
   end
  end
  if any(strcmp(PoCETsys.pce.vars,PoCETsys.(vars)(i).name))
   samples.xi(strcmp(PoCETsys.pce.vars,PoCETsys.(vars)(i).name),:) = samples.(vars)(i,:);
  end
 end
end

fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function samples = sampleStandardRandomVariables(PoCETsys, numSamples)

% convert into standard random variables and then call sampleRandomVariables
for vv = {'states','parameters'},
 vars = vv{1};
 for i=1:numel(PoCETsys.(vars)),
  switch lower(PoCETsys.(vars)(i).dist),
   case 'normal',
    PoCETsys.(vars)(i).data = [+0,+1];
   case 'uniform',
    PoCETsys.(vars)(i).data = [-1,+1];
   case 'beta',
    PoCETsys.(vars)(i).dist = 'beta4';
    PoCETsys.(vars)(i).data(3) = -1;
    PoCETsys.(vars)(i).data(4) = 1;
   case 'beta4',
    PoCETsys.(vars)(i).data(3) = -1;
    PoCETsys.(vars)(i).data(4) = 1;
   case {'none','dirac'}
    PoCETsys.(vars)(i).data = [1];
   otherwise
    error('Wrong input.');
  end
 end
end
samples = sampleRandomVariables(PoCETsys, numSamples);
end
