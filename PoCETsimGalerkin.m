function results = PoCETsimGalerkin(PoCETsys, RHSfile, OUTfile, simoptions, varargin)
% PoCETsim = PoCETsimGalerkin(PoCETsys, RHSfile, OUTfile) 
% simulates the system PoCETsys and returns a struct containing the PCE
% coefficients for all states and outputs specified in PoCETsys 
% The coefficients are calculated via the Galerkin method. 
% 
% - PoCETsys is the system struct returned by PoCETcompose
% - RHSfile is a string containing the filename of the states'
%   right-hand-sides, which can be created using PoCETwriteFiles
% - OUTfile is a string containing the filename of the system's output,
%   which can be created using PoCETwriteFiles
% - simoptions is a struct with fields 
%   .setup   -  odeset-struct
%   .solver 	-  string containing the solver's name; default is 'ode15s'
%   .tspan   -  vector [t_0, t_end] or [t_0, t_1, ..., t_end]
%   .dt      -  double defining the stepwidth of the solver's output;
%               if tspan is a vector of more than two timepoints, setting
%               dt will overwrite these values
%
% PoCETsim = PoCETsimGalerkin(..., 'in1', data1, ..., 'inN', dataN) 
% simulates the system with inputs 'in1' to 'inN'. The inputs' names and
% data have to match those specified in PoCETsys.inputs (see examples for
% more details on working with inputs).

%% process input
if nargin == 2
 simoptions = struct('setup',odeset,'solver','ode45','tspan',[0 5]);
 OUTfile = []; invars = []; invals = [];
elseif nargin == 3
 if ischar(OUTfile) || isempty(OUTfile)
  simoptions = struct('setup',odeset,'solver','ode45','tspan',[0 5]);
 elseif isstruct(OUTfile) && any(ismember({'setup' 'solver' 'tspan' 'dt'},fieldnames(OUTfile)))
  simoptions = OUTfile; OUTfile = [];
 else
  error('Invalid input %s. Enter ''help PoCETsimGalerkin'' for expected input specifications.',inputname(3))
 end
 invars = []; invals = [];
else
 if isempty(simoptions)
  simoptions = struct('setup',odeset,'solver','ode45','tspan',[0 5]);
 elseif ~isstruct(simoptions)
  error('Input argument simoptions is expected to be a struct!')
 end
 invars = varargin(1:2:end);
 invals = varargin(2:2:end);
 
 if numel(invars)~=numel(invals)
  error('Wrong specification of system inputs. Enter ''help PoCETsimGalerkin'' for more details.')
 end
 if ~all(cellfun(@(x)ischar(x),invars))
  error('Wrong specification of system inputs. Enter ''help PoCETsimGalerkin'' for more details.')
 end
end

% if some solver-specifications are not set use default values
specs = fieldnames(simoptions);
if ~ismember('setup',specs);  simoptions.setup  = odeset;  end
if ~ismember('solver',specs); simoptions.solver = 'ode45'; end
if ~ismember('tspan',specs);  simoptions.tspan  = [0, 5];  end

% check, if all input parameters are set
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

if steps == 1
 error('Invalid time-span definition! Please define at least one start and final time point!')
end

n_phi = PoCETsys.pce.options.n_phi;

%% set initial conditions
x0 = zeros(numel(PoCETsys.states)*n_phi,1);
for i = 1:numel(PoCETsys.states)
 x0((i-1)*n_phi+1:i*n_phi) = PoCETsys.states(i).pce;
end

%% simulation
if ~exist(RHSfile,'file');
 error('Right-hand-side-file %s.m does not exist.',RHSfile);
else
 RHS = eval(sprintf('@(t,x)%s(t,x,PoCETsys%s)',RHSfile,writeIPAR));
end
 
if strcmp(simoptions.solver,'discrete-time') % discrete time simulation
 x = zeros(numel(x0),steps); % preallocation
 x(:,1) = x0;
 for k = 1:steps-1  % simulation
  x(:,k+1) = RHS(tspan(k),x(:,k));
 end
elseif strcmp(simoptions.solver,'explicit-euler') 
 x = zeros(numel(x0),steps); % preallocation
 x(:,1) = x0;
 for k = 1:steps-1  % simulation
  x(:,k+1) = x(:,k) + simoptions.dt*RHS(tspan(k),x(:,k)); % continuous time simulation with euler method
 end
else
 [~,x] = feval(simoptions.solver,RHS,tspan,x0,simoptions.setup); % continuous time simulation 
 x = x';
end

for i = 1:numel(PoCETsys.states)
 results.(PoCETsys.states(i).name).pcvals = x((i-1)*n_phi+1:i*n_phi,:);
end
results.time = tspan;

if ~isempty(OUTfile)
 if exist(OUTfile,'file');
  OUT = eval(sprintf('@(t,x)%s(t,x,PoCETsys%s)',OUTfile,writeIPAR));
  y = zeros(numel(PoCETsys.outputs)*n_phi,steps);
  for k = 1:steps
   y(:,k) = OUT(tspan(k),x(:,k));
  end
  for i = 1:numel(PoCETsys.outputs)
   results.(PoCETsys.outputs(i).name).pcvals = y((i-1)*n_phi+1:i*n_phi,:);
  end
 else
  error('Output-file %s.m does not exist.',OUTfile);
 end
end
end

%% Codesnippet for alternative setting of inital conditions
%% (moment-independent; requires different base polynomials!)
% for i = 1:numel(PoCETsys.states)
%     sim.x(i).pcvals = zeros(n_phi,steps); % contains all expansion values of state i
%     switch lower(PoCETsys.states(i).dist),
%         case 'uniform'
%             sim.x(i).pcvals(1,1)   = PoCETsys.states(i).moments(1);
%             if numel(PoCETsys.states(i).data) > 1
%                 sim.x(i).pcvals(i+1,1) = (PoCETsys.states(i).data(2)-PoCETsys.states(i).data(1))/2;
%             end            
%         case 'normal'
%             sim.x(i).pcvals(1,1)   = PoCETsys.states(i).moments(1);
%             if numel(PoCETsys.states(i).data) > 1
%                 sim.x(i).pcvals(i+1,1) = PoCETsys.states(i).moments(2);
%             end
%         otherwise
%             error('@Felix: please fix!');
%     end
%     if numel(PoCETsys.states(i).data) > 1
%         sim.x(i).pcvals(i+1,1) = PoCETsys.states(i).moments(2);
%     end
% end