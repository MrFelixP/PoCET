function results = PoCETsimCollocation(PoCETsys, MCRHSfile, MCOUTfile, BasisSamples, simoptions, varargin)
% PoCETsim = PoCETsimCollocation(PoCETsys, 'MCRHSfile', 'MCOUTfile', BasisSamples, simoptions) 
% simulates the system PoCETsys and returns a struct containing the PCE
% coefficients for all states and outputs specified in PoCETsys 
% The coefficients are calculated via the stochastic collocation method. 
% 
% - PoCETsys is the system struct returned by PoCETcompose
% - MCRHSfile is a string containing the filename of the states'
%   right-hand-sides, which can be created using PoCETwriteFiles
% - MCOUTfile is a string containing the filename of the system's output,
%   which can be created using PoCETwriteFiles
% - samples is either an integer defining the number of collocation-points
%   or a matrix containig basis samples, which can be created using
%   PoCETsample(PoCETsys, 'basis', N)
% - simoptions is a struct with fields 
%   .setup   -  odeset-struct
%   .solver 	-  string containing the solver's name; default is 'ode15s'
%   .tspan   -  vector [t_0, t_end] or [t_0, t_1, ..., t_end]
%   .dt      -  double defining the stepwidth of the solver's output;
%               if tspan is a vector of more than two timepoints, setting
%               dt will overwrite these values
%
% PoCETsim = PoCETsimCollocation(..., 'method', 'final')
% returns a struct only containing the PCE coefficient values of all states
% and outputs at time simoptions.tspan(end)
%
% PoCETsim = PoCETsimCollocation(..., 'in1', data1, ..., 'inN', dataN) 
% simulates the system with inputs. The inputs' names and data have to
% match those specified in PoCETsys.inputs (see examples for more details on
% working with inputs).

if ~isempty(MCOUTfile)
 if exist(MCOUTfile,'file');
  nOUT = numel(PoCETsys.outputs);
 else
  error('Output file %s.m does not exist.',MCOUTfile);
 end
else
 nOUT = 0;
end

% sampling
if max(size(BasisSamples)) == 1
 nsamples = BasisSamples;
 BasisSamples = PoCETsample(PoCETsys,'basis',nsamples);
 samples = struct('states',zeros(numel(PoCETsys.states),nsamples), ...
                  'parameters',zeros(numel(PoCETsys.parameters),nsamples));
 for i = 1:numel(PoCETsys.states)
  samples.states(i,:) = PoCETsys.states(i).pce'*BasisSamples;
 end
 for i = 1:numel(PoCETsys.parameters)
  samples.parameters(i,:) = PoCETsys.parameters(i).pce'*BasisSamples;
 end
elseif any(size(BasisSamples) == PoCETsys.pce.options.n_phi)
 [ro co] = size(BasisSamples);
 if ro == PoCETsys.pce.options.n_phi
  nsamples = co;
 else
  nsamples = ro;
 end
 samples = struct('states',zeros(numel(PoCETsys.states),nsamples), ...
                  'parameters',zeros(numel(PoCETsys.parameters),nsamples));
 for i = 1:numel(PoCETsys.states)
  samples.states(i,:) = PoCETsys.states(i).pce'*BasisSamples;
 end
 for i = 1:numel(PoCETsys.parameters)
  samples.parameters(i,:) = PoCETsys.parameters(i).pce'*BasisSamples;
 end
else
 error('Input argument ''BasisSamples'' is expected to be an integer or a matrix of basis samples created with PoCETsample(PoCETsys, ''basis'', N)!')
end

% simulation
mc_vals = PoCETsimMonteCarlo(PoCETsys,MCRHSfile,MCOUTfile,samples,simoptions,varargin{:});

% calculation of pce coefficients
basis_pinv = pinv(BasisSamples);
for i = 1:numel(PoCETsys.states)
 results.(PoCETsys.states(i).name).pcvals = ...
 basis_pinv'*mc_vals.(PoCETsys.states(i).name).mcvals;
end
for i = 1:nOUT
 results.(PoCETsys.outputs(i).name).pcvals = ...
 basis_pinv'*mc_vals.(PoCETsys.outputs(i).name).mcvals;
end
results.time = mc_vals.time;
end