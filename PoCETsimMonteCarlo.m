function MCsim = PoCETsimMonteCarlo(PoCETsys, MCRHSfile, MCOUTfile, samples, simoptions, varargin)
% MCsim = PoCETsimMonteCarlo(PoCETsys, 'MCRHSfile', 'MCOUTfile', samples, simoptions)
% simulates the system PoCETsys and returns a struct containing the
% Monte-Carlo trajectories for all states and outputs in PoCETsys.
%
% - PoCETsys is the system struct returned by PoCETcompose
% - MCRHSfile is a string containing the filename of the states'
%   right-hand-sides, which can be created using PoCETwriteFiles
% - MCOUTfile is a string containing the filename of the output'
%   right-hand-side, which can be created using PoCETwriteFiles
% - samples is either an integer defining the number of Monte-Carlo-samples
%   or a struct containig samples for all variables, which can be created
%   using PoCETsample(PoCETsys, 'variables', N)
% - simoptions is a struct with fields 
%   .setup   -  odeset-struct
%   .solver 	-  string containing the solver's name; default is 'ode15s'
%   .tspan   -  vector [t_0, t_end] or [t_0, t_1, ..., t_end]
%   .dt      -  double defining the stepwidth of the solver's output;
%               if tspan is a vector of more than two timepoints, setting
%               dt will overwrite these values
% 
% MCsim = PoCETsimMonteCarlo(..., 'method', 'final')
% returns a struct only containing the values of all states and outputs at 
% time simoptions.tspan(end)
%
% MCsim = PoCETsimMonteCarlo(..., 'method', 'moments')
% returns a struct containing the central moments of all states and outputs
% which are calculated recursively during the simulation in order to save
% memory
%
% MCsim = PoCETsimMonteCarlo(..., 'in1', data1, ..., 'inN', dataN)
% simulates the system with inputs 'in1' to 'inN'. The inputs' names and
% data have to match those specified in PoCETsys.inputs (see the examples
% for more details on working with inputs).

%% process input and check for errors
% get number of states
nSTATES = numel(PoCETsys.states);

% check sample specification
if nargin < 2
 error('Not enough input arguments.')
elseif nargin == 2 || nargin == 3
 loops = 1;
 samples = PoCETsample(PoCETsys,'variables',loops);
elseif nargin > 3
 if isnumeric(samples)
  loops = floor(samples);
  samples = PoCETsample(PoCETsys,'variables',loops);
 elseif isstruct(samples)
  samplefields = fieldnames(samples);
  if ismember('states',samplefields) && ismember('parameters',samplefields)
   chkStates = nSTATES == size(samples.states,1);
   chkParams = numel(PoCETsys.parameters) == size(samples.parameters,1);
   if chkStates && chkParams
    loops = size(samples.states,2);
   elseif ~chkStates
    error('Number of given state-samples does not match number of states!')
   else
    error('Number of given parameter-samples does not match number of parameters!')
   end
  else
   error('Provided samples-struct does not match expected specifications!')
  end
 else
  error('Input argument ''samples'' is expected to be an integer or a struct of samples created with PoCETsample!')
 end
end

% check method
if nargin < 6
 method = 'complete';
end

% check simulation options and inputs
if nargin < 5
 simoptions = struct('setup',odeset,'solver','ode45','tspan',[0 5]);
 invars = []; invals = [];
else
 if isstruct(simoptions)
  specs = fieldnames(simoptions);
  if ~ismember('setup',specs);  simoptions.setup =  odeset;  end
  if ~ismember('solver',specs); simoptions.solver = 'ode45'; end
  if ~ismember('tspan',specs);  simoptions.tspan =  [0, 5];  end
 elseif isempty(simoptions)
  simoptions = struct('setup',odeset,'solver','ode45','tspan',[0 5]);
 else
  error('Input argument simoptions is expected to be a struct.')
 end
 invars = varargin(1:2:end);
 invals = varargin(2:2:end);
 if any(strcmp(invars,'method'))
  method = invals(strcmp(invars,'method'));
 else
  method = 'complete';
 end
 if numel(invars)~=numel(invals)
  error('Additional function inputs must come in pairs. Enter ''help simulatePoCETsys'' for more details.')
 end
 if ~all(cellfun(@(x)ischar(x),invars))
  error('Wrong specification of system inputs. Enter ''help simulatePoCETsys'' for more details.')
 end
end

% evaluate input parameters
IPAR = fieldnames(PoCETsys.inputs);
writeIPAR = '';
for i = 1:numel(IPAR) 
 if ~(strcmpi(IPAR{i},'rhs') || strcmpi(IPAR{i},'name'))
  if ismember(IPAR{i},invars)
   parstring = [IPAR{i} '=[' num2str(invals{strcmp(invars,IPAR{i})}) '];'];
  else
   parstring = [IPAR{i} '=[' num2str([PoCETsys.inputs(:).(IPAR{i})]) '];'];
  end
  eval(parstring);
  writeIPAR = [writeIPAR ',' IPAR{i}];
 end
end

% check output file and get number of outputs
if ~isempty(MCOUTfile)
 if exist(MCOUTfile,'file');
  OUT = eval(sprintf('@(t,x,PAR)%s(t,x,PAR%s)',MCOUTfile,writeIPAR));
  nOUT = numel(PoCETsys.outputs);
 else
  error('Output file %s.m does not exist.',MCOUTfile);
 end
else
 OUT = @(t,x,p)[];
 nOUT = 0;
 chkout = cellfun(@(x)~isempty(x),{PoCETsys.outputs(:).name});
 if any(chkout)
  fprintf(2,'! WARNING: No output file specified! Output %s cannot be evaluated!\n',PoCETsys.outputs(chkout).name);
 end
end

% check time specifications and create time-vector
if isfield(simoptions,'dt')
 if simoptions.dt ~= 0
  tspan = simoptions.tspan(1):simoptions.dt:simoptions.tspan(end);
 else
  tspan = simoptions.tspan;
 end
else
 tspan = simoptions.tspan;
end
steps = numel(tspan);
% 
% if steps == 1
%  error('Invalid time-span definition! Please define at least one start and final time point!')
% end

%% preallocation
%-------------------------------------------------------------------------%
if strcmpi(method,'moments')
 for i = 1:nSTATES
  MCsim.(PoCETsys.states(i).name).moments = zeros(4,steps);
 end
 x(1:nSTATES) = ...
    struct('lb',+inf(1,steps),'ub',-inf(1,steps),'m1',zeros(1,steps), ...
   'sqsum',zeros(1,steps),'cusum',zeros(1,steps),'qusum',zeros(1,steps));
 for i = 1:nOUT
  MCsim.(PoCETsys.outputs(i).name).moments = zeros(4,steps);
 end
 y(1:nOUT) = ...
    struct('lb',+inf(1,steps),'ub',-inf(1,steps),'m1',zeros(1,steps), ...
   'sqsum',zeros(1,steps),'cusum',zeros(1,steps),'qusum',zeros(1,steps));
%-------------------------------------------------------------------------%
elseif strcmpi(method,'final')
 for i = 1:nSTATES
  MCsim.(PoCETsys.states(i).name).mcvals = zeros(loops,1);
 end
 for i = 1:nOUT
  MCsim.(PoCETsys.outputs(i).name).mcvals = zeros(loops,1);
 end
%-------------------------------------------------------------------------%
elseif strcmpi(method,'complete')
 for i = 1:nSTATES
  MCsim.(PoCETsys.states(i).name).mcvals = zeros(loops,steps);
 end
 for i = 1:nOUT
  MCsim.(PoCETsys.outputs(i).name).mcvals = zeros(loops,steps);
 end
end

%% Simulation
fprintf('Performing Monte-Carlo simulation of %i samples:   0.0%% finished\n',loops);
for j = 1:loops
 
%  show progresss
 if floor(mod(j,loops/20)) == 0,
  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6s%% finished',sprintf('%3.1f',j/loops*100));
 end

 % get samples
 x0 = samples.states(:,j);
 for i = 1:numel(PoCETsys.parameters),
  PAR.(PoCETsys.parameters(i).name) = samples.parameters(i,j);
 end
 
 % re-evaluate RHS-function (necessary in each loop to match ode syntax)
 RHS = eval(sprintf('@(t,x)%s(t,x,PAR%s)',MCRHSfile,writeIPAR)); 
 
 % do simulation
 if numel(tspan) > 1
  if strcmp(simoptions.solver,'discrete-time')
   mc = zeros(numel(x0),steps); % preallocation (already right dimensions)
   mc(:,1) = x0; % set initial conditions
   for k = 1:steps-1  % simulation
    mc(:,k+1) = RHS(tspan(k),mc(:,k)); % disc-time-sim
   end
  elseif strcmp(simoptions.solver,'explicit-euler')
   mc = zeros(numel(x0),steps); % preallocation (already right dimensions)
   mc(:,1) = x0; % set initial conditions
   for k = 1:steps-1  % simulation
    mc(:,k+1) = mc(:,k) + simoptions.dt*RHS(tspan(k),mc(:,k));
   end
  else
   [~,mc] = feval(simoptions.solver,RHS,tspan,x0,simoptions.setup); % cont-time-sim
   mc = mc'; % transpose output to match dimensions in following steps
  end
 else
  mc = x0;
 end
 
 % evaluate output
 mcout = OUT(tspan,mc,PAR);
 
%-------------------------------------------------------------------------%
 if strcmpi(method,'moments')
  for i = 1:nSTATES
   x(i).lb(x(i).lb>mc(i,:)) = mc(i,x(i).lb>mc(i,:));
   x(i).ub(x(i).ub<mc(i,:)) = mc(i,x(i).ub<mc(i,:));
   x(i).m1 = (x(i).m1*(j-1)+mc(i,:))/j;
   x(i).sqsum = x(i).sqsum + mc(i,:).^2;
   x(i).cusum = x(i).cusum + mc(i,:).^3;
   x(i).qusum = x(i).qusum + mc(i,:).^4;
  end
  for i = 1:nOUT
   y(i).lb(y(i).lb>mcout(i,:)) = mcout(i,y(i).lb>mcout(i,:));
   y(i).ub(y(i).ub<mcout(i,:)) = mcout(i,y(i).ub<mcout(i,:));
   y(i).m1 = (y(i).m1*(j-1) + mcout(i,:))/j;
   y(i).sqsum = y(i).sqsum + mcout(i,:).^2;
   y(i).cusum = y(i).cusum + mcout(i,:).^3;
   y(i).qusum = y(i).qusum + mcout(i,:).^4;
  end
%-------------------------------------------------------------------------%
 elseif strcmpi(method,'final')
  for i = 1:nSTATES
   MCsim.(PoCETsys.states(i).name).mcvals(j) = mc(i,end);
  end
  for i = 1:nOUT
   MCsim.(PoCETsys.outputs(i).name).mcvals(j) = mcout(i,end);
  end 
%-------------------------------------------------------------------------%
 elseif strcmpi(method,'complete')
  for i = 1:nSTATES
   MCsim.(PoCETsys.states(i).name).mcvals(j,:) = mc(i,:);
  end
  for i = 1:nOUT
   MCsim.(PoCETsys.outputs(i).name).mcvals(j,:) = mcout(i,:);
  end
 end
end

% calculate moments if method is 'moments'
if strcmpi(method,'moments')
 for i = 1:nSTATES
  m1 = x(i).m1;
  m2 = (x(i).sqsum-loops*m1.^2)/(loops-1);
  m3 = (x(i).cusum-3*m1.*x(i).sqsum+2*loops*m1.^3)./(loops*sqrt(m2).^3);
  m4 = (x(i).qusum-4*x(i).cusum.*m1+6*x(i).sqsum.*m1.^2-3*loops*m1.^4)./(loops*sqrt(m2).^4)-3;
  
  MCsim.(PoCETsys.states(i).name).moments(1,:) = m1;
  MCsim.(PoCETsys.states(i).name).moments(2,:) = m2;
  MCsim.(PoCETsys.states(i).name).moments(3,:) = m3;
  MCsim.(PoCETsys.states(i).name).moments(4,:) = m4;
  
  MCsim.(PoCETsys.states(i).name).lb = x(i).lb;
  MCsim.(PoCETsys.states(i).name).ub = x(i).ub;
 end
 for i = 1:nOUT
  m1 = y(i).m1;
  m2 = (y(i).sqsum-loops*m1.^2)/(loops-1);
  m3 = (y(i).cusum-3*m1.*y(i).sqsum+2*loops*m1.^3)./(loops*sqrt(m2).^3);
  m4 = (y(i).qusum-4*y(i).cusum.*m1+6*y(i).sqsum.*m1.^2-3*loops*m1.^4)./(loops*sqrt(m2).^4)-3;
  
  MCsim.(PoCETsys.outputs(i).name).moments(1,:) = m1;
  MCsim.(PoCETsys.outputs(i).name).moments(2,:) = m2;
  MCsim.(PoCETsys.outputs(i).name).moments(3,:) = m3;
  MCsim.(PoCETsys.outputs(i).name).moments(4,:) = m4;
  
  MCsim.(PoCETsys.outputs(i).name).lb = y(i).lb;
  MCsim.(PoCETsys.outputs(i).name).ub = y(i).ub;
 end
end

MCsim.time = tspan;
fprintf('\n');
end