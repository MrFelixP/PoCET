function [PoCETsys] = PoCETcompose(states, parameters, inputs, outputs, order)
% PoCETsys = PoCETcompose(states, parameters, inputs, outputs) 
% returns a struct containing the Polynomial Chaos Expansion of order 2 of 
% the system defined by the states-, parameters-, and inputs-structs used 
% as input arguments. See 'help PoCET' for more information on how to define
% a system.
% 
% PoCETcompose(states, parameters, inputs, outputs, order) additionally sets 
% the order of the Polynomial Chaos Expansion.
%
% NOTE: PoCETcompose(s)      does the same as  PoCETcompose(s,[],[],[],2)
%       PoCETcompose(s,p)    does the same as  PoCETcompose(s,p,[],[],2)
%       PoCETcompose(s,p,i)  does the same as  PoCETcompose(s,p,i,[],2)
%
% In order to simulate the composed system, first create its expanded
% right-hand-sides using PoCETwriteFiles(PoCETsys).

% check input
if nargin == 1, parameters = []; inputs = []; outputs = []; order = 2;
 elseif nargin == 2, inputs = []; outputs = []; order = 2;
 elseif nargin == 3, outputs = []; order = 2;
 elseif nargin == 4, order = 2;
 elseif nargin > 5, error('Too many input arguments!')
end

% set dummy structs, if input is empty
if isempty(parameters), parameters = struct('name','','dist','','data',[],'moments',[]); end
if isempty(inputs), inputs = struct('name','','rhs','','vars',''); end
if isempty(outputs), outputs = struct('name','','rhs','','vars',''); end
if isempty(order), order = 2; end

for i = 1:numel(states)
 states(i).name = char(states(i).name);
 states(i).rhs  = char(states(i).rhs);
end
for i = 1:numel(parameters)
 parameters(i).name = char(parameters(i).name);
end
for i = 1:numel(inputs)
 inputs(i).name = char(inputs(i).name);
 inputs(i).rhs = char(inputs(i).rhs);
end
for i = 1:numel(outputs)
 outputs(i).name = char(outputs(i).name);
 outputs(i).rhs = char(outputs(i).rhs);
end

fprintf('Analyzing input...\n')
% create cells with variable names for later comparisons
state_names = {states.name}; 
param_names = {parameters.name};
input_names = {inputs.name};
output_names = {outputs.name};
all_names = [state_names, param_names, input_names, output_names];
all_names = all_names(~cellfun('isempty',all_names));

%% now do the real work
% preallocate system struct
PoCETsys.states = states;
PoCETsys.parameters = parameters;
PoCETsys.inputs = inputs;
PoCETsys.outputs = outputs;
PoCETsys.pce.vars = {};
PoCETsys.pce.pars = {};
PoCETsys.pce.dep_vars = {};

% check if states' ICs are distributed
 for i = 1:numel(states)
  if ~any([strcmp(states(i).dist,'none'),strcmp(states(i).dist,'pce')]);
   PoCETsys.pce.vars = [PoCETsys.pce.vars {states(i).name}];
  end
 end

% check if parameters are distributed
for i = 1:numel(parameters) 
 if ~any([strcmp(parameters(i).dist,'none'),strcmp(parameters(i).dist,'pce')]) && ~isempty(parameters(i).name);
  PoCETsys.pce.vars = [PoCETsys.pce.vars {parameters(i).name}];
  PoCETsys.pce.pars = [PoCETsys.pce.pars {parameters(i).name}];
 end
end

% Collect PCE options
PoCETsys.pce.options.order = order;
PoCETsys.pce.options.quadtol = 1e-8;
PoCETsys.pce.options.n_xi = numel(PoCETsys.pce.vars);
PoCETsys.pce.options.n_phi = ...
     factorial(PoCETsys.pce.options.n_xi+PoCETsys.pce.options.order) / ...
    (factorial(PoCETsys.pce.options.n_xi)*factorial(PoCETsys.pce.options.order));

InCheckMsg = checkInput(PoCETsys);

if ~isempty(InCheckMsg)% if there is a wrong input
 fprintf(2,[InCheckMsg '\n']);
 error('Invalid input: Please check your input file or enter ''help PoCET'' for expected input specifications!');
end

%% Analyse ODEs and output functions
VarCheckMsg = '';
for i=1:numel(PoCETsys.states)
 if ischar(states(i).rhs)
  PoCETsys.ode(i) = eq_decompose(states(i).rhs,'s');
 else
  error('RHS of state %s is of class %s but expected to be of class char or sym!',PoCETsys.states(i).name,class(states(i).rhs))
 end
end

for i=1:numel(outputs)
 if ischar(outputs(i).rhs)
  PoCETsys.out(i) = eq_decompose(outputs(i).rhs,'o');
 else
  error('RHS of output %s is of class %s but expected to be of class char or sym!',PoCETsys.outputs(i).name,class(outputs(i).rhs))
 end
end

for i=numel(inputs):-1:1
 if ischar(inputs(i).rhs)
  tmp.in(i) = eq_decompose(inputs(i).rhs,'i');
 else
  error('RHS of input %s is of class %s but expected to be of class char or sym!',PoCETsys.inputs(i).name,class(inputs(i).rhs))
 end
end

hpo = 0;
for i = 1:numel(PoCETsys.ode)
 for j = 1:numel(PoCETsys.ode(i).terms)
  hpo = max([hpo, sum(PoCETsys.ode(i).terms(j).state_degrees)]);
 end
end

for i = 1:numel(PoCETsys.out)
 for j = 1:numel(PoCETsys.out(i).terms)
  hpo = max([hpo, sum(PoCETsys.out(i).terms(j).state_degrees)]);
 end
end

PoCETsys.pce.options.hpo = hpo + 1;

%% check if all variables are used & warn user if not
vars_used = unique([PoCETsys.ode.vars PoCETsys.out.vars tmp.in.vars]);
% names_defined = [state_names param_names input_names];
names_used = cell(1,numel(vars_used));
for i = 1:numel(vars_used)
 names_used{i} = char(vars_used(i));
end
for i = 1:numel(param_names)
 cur = param_names{i};
 if ~ismember(cur,names_used) && ~isempty(cur)
  fprintf(2,'! WARNING: Parameter %s is unused!\n',cur);
 end
end
for i = 1:numel(input_names)
 cur = input_names{i};
 if ~ismember(cur,names_used) && ~isempty(cur)
  fprintf(2,'! WARNING: Input variable %s is unused!\n',cur);
 end
end

if ~isempty(VarCheckMsg)
 fprintf(2,[VarCheckMsg '\n']);
 error('Invalid input: Please check your input file or enter ''help PoCET'' for expected input specifications!');
end

%% orthogonal basis
fprintf('Determining orthogonal basis...\n')
if ~isempty(PoCETsys.pce.vars)
 varsStr = cell(1,numel(PoCETsys.pce.vars));
 for i = 1:numel(PoCETsys.pce.vars),
  varsStr{i} = ['xi_',num2str(i)];
 end
 PoCETsys.randomVariables = sym(varsStr);
else
 PoCETsys.randomVariables = [];
end
PoCETsys.basis = PoCETbasis(PoCETsys);
%% quadrature rules
fprintf('Calculating quadrature rules...\n')
PoCETsys.quadrature = PoCETquadRules(PoCETsys);

%% initial conditions
for vv = {'states','parameters'},
 vars = vv{1};
 for i=1:numel(PoCETsys.(vars)),
  xi_i = find(ismember(PoCETsys.pce.vars,PoCETsys.(vars)(i).name));
  PoCETsys.(vars)(i).pce = zeros(PoCETsys.pce.options.n_phi,1); % contains pce-coefficients of variable
  switch lower(PoCETsys.(vars)(i).dist),
   case 'pce'
    if numel(PoCETsys.(vars)(i).data) == PoCETsys.pce.options.n_phi;
     PoCETsys.(vars)(i).pce = PoCETsys.(vars)(i).data;
    else
     error('Size of PCE coefficients specified for %s %s does not match size (%i) of PCE basis.',vars(1:end-1),PoCETsys.(vars)(i).name,PoCETsys.pce.options.n_phi)
    end
   case 'normal',
    PoCETsys.(vars)(i).pce(1) = PoCETsys.(vars)(i).data(1);
    PoCETsys.(vars)(i).pce(xi_i+1) = PoCETsys.(vars)(i).data(2);
   case 'uniform',
    PoCETsys.(vars)(i).pce(1) = (PoCETsys.(vars)(i).data(2)+PoCETsys.(vars)(i).data(1))/2;
    PoCETsys.(vars)(i).pce(xi_i+1) = (PoCETsys.(vars)(i).data(2)-PoCETsys.(vars)(i).data(1))/2;
   case 'beta',
    alf = PoCETsys.(vars)(i).data(1);
    bet = PoCETsys.(vars)(i).data(2);
    PoCETsys.(vars)(i).pce(1) = alf/(alf+bet);
    PoCETsys.(vars)(i).pce(xi_i+1) = 1/(alf+bet);
   case 'beta4',
    alf = PoCETsys.(vars)(i).data(1);
    bet = PoCETsys.(vars)(i).data(2);
    low = PoCETsys.(vars)(i).data(3);
    upp = PoCETsys.(vars)(i).data(4);
    PoCETsys.(vars)(i).pce(1) = low + (upp-low)*alf/(alf+bet);
    PoCETsys.(vars)(i).pce(xi_i+1) = (upp-low)/(alf+bet);
   case {'none','dirac'}
    PoCETsys.(vars)(i).pce(1) = PoCETsys.(vars)(i).data(1);
   otherwise
    error('Invalid distribution ''%s'' for %s %s.',PoCETsys.(vars)(i).dist,vars(1:end-1),PoCETsys.(vars)(i).name)
  end
 end
end

%% coefficient matrices
if all([PoCETsys.ode.ispoly PoCETsys.out.ispoly])
 fprintf('Calculating coefficient matrices...\n')
 [PoCETsys.coeff_matrices PoCETsys.coeff_matrices.raw] = PoCETcoeffMatrices(PoCETsys);
end

%% output
fprintf('\nPoCET finished composing the expanded system.\n')
nvars = numel(PoCETsys.pce.vars);
if nvars == 1
 fprintf('%i random variables has been identified:\n',nvars)
else
 fprintf('%i random variables have been identified:\n',nvars)
end
for i = 1:nvars
 fprintf(' xi_%i = %s\n',i,PoCETsys.pce.vars{i})
end
fprintf('\nTo simulate the system, first write it to file using PoCETwriteFiles and then run PoCETsimGalerkin.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = degreeSym(in,var)
    out = double(feval(symengine,'degree',sym(in),sym(var)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eq_dec = eq_decompose(eq_string,varflag)
 % prepare ode decomposition
 
 eq_dec.terms = struct('monomial',[],'coefficient',[],...
                       'xi_index',[],'xi_degrees',[],...
                       'state_index',[],'state_degrees',[],...
                       'param_index',[],'param_degrees',[],...
                       'input_index',[],'input_degrees',[],...
                       'output_index',[],'output_degrees',[]); 
 eq_dec.ispoly = true;
 eq_dec.vars = [];
 eq_dec.equation = [];
 eq_dec.coeffs = [];
 eq_dec.mons = [];    
 
 if isempty(eq_string)
  % if eq_string is empty: return default values
  return;
 end
 
 [eq_dec.ispoly,eq_dec.vars] = ispolynomial(eq_string,all_names,varflag);
 
 switch varflag
  case 's'
   varname = ['state ' state_names{i}];
   in_fields = {};
  case 'o'
   varname = ['output ' output_names{i}];
   in_fields = {};
  case 'i'
   varname = ['input ' input_names{i}];
   in_fields = fieldnames(inputs)';
   in_fields(strcmp(in_fields,'name')) = [];
   in_fields(strcmp(in_fields,'rhs')) = [];
   in_fields = [in_fields 't'];
 end
 
 if ~isempty(eq_dec.vars) % not polynomial but convertible to sym
  if any(~ismember(eq_dec.vars,[all_names, in_fields])),
   vars = find(~ismember(eq_dec.vars,all_names));
   for jj = 1:numel(vars)
    VarCheckMsg = [VarCheckMsg 'Variable ' char(eq_dec.vars(vars(jj))) ...
        ' in RHS of ' varname ' is undefined!\n'];
   end
  end
 end
 
 if eq_dec.ispoly
  
  dum_var = sym(all_names); % create symbolic dummy variables
  dum_str = [' ' eq_string ' ']; % create dummy equation string with extra blankspace
  
  for i_var = 1:numel(all_names)
   dum_str = regexprep(dum_str,['(\W)' all_names{i_var} '(\W)'],['$1dum_var(' num2str(i_var) ')$2']);
  end
  
  eq_dec.equation = str2syms(dum_str); % convert to symbolic
  if ~isempty(eq_dec.vars)
   [eq_dec.coeffs, eq_dec.mons] = coeffs(eq_dec.equation,eq_dec.vars); % decompose terms
  else
   eq_dec.coeffs = coeffs(eq_dec.equation,eq_dec.vars); % decompose terms
   eq_dec.mons = [];
  end
  for kk = 1:numel(eq_dec.mons), % cycle through terms
   eq_dec.terms(kk).monomial = char(eq_dec.mons(kk)); % convert monomials back to char
   eq_dec.terms(kk).coefficient = double(eq_dec.coeffs(kk)); % convert coefficients to double
   vars = symvar(eq_dec.mons(kk)); % get single variables in current term
   for jj = 1:numel(vars), % cycle through variables and chek if state, parameter, or input
    if ~isempty(find(vars(jj) == state_names,1)), % is state?
     eq_dec.terms(kk).state_index(jj) = find(vars(jj) == state_names);
     eq_dec.terms(kk).state_degrees(jj) = degreeSym(eq_dec.mons(kk),vars(jj));
    elseif ~isempty(find(vars(jj) == param_names,1)) % is parameter?
     eq_dec.terms(kk).param_index(jj) = find(vars(jj) == param_names);
     eq_dec.terms(kk).param_degrees(jj) = degreeSym(eq_dec.mons(kk),vars(jj));
    elseif ~isempty(find(vars(jj) == input_names,1)) % is input?
     eq_dec.terms(kk).input_index(jj) = find(vars(jj) == input_names);
     eq_dec.terms(kk).input_degrees(jj) = degreeSym(eq_dec.mons(kk),vars(jj));
    elseif ~isempty(find(vars(jj) == output_names,1)) % is output?
     eq_dec.terms(kk).output_index(jj) = find(vars(jj) == output_names);
     eq_dec.terms(kk).output_degrees(jj) = degreeSym(eq_dec.mons(kk),vars(jj));
    end
    if ~isempty(find(vars(jj) == PoCETsys.pce.vars,1)) % check if variable is random
     eq_dec.terms(kk).xi_index(jj)   = find(vars(jj) == PoCETsys.pce.vars);
     eq_dec.terms(kk).xi_degrees(jj) = degreeSym(eq_dec.mons(kk),vars(jj));
    end
   end
  end
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [is_poly,sym_var] = ispolynomial(eq_str,all_names,varflag)
 
 is_poly = false;
 sym_var = [];

 dum_var = sym(all_names); % create symbolic dummy variables
 dum_str = [' ' eq_str ' ']; % create dummy equation string with extra blankspaces
 
 for i_var = 1:numel(all_names)
  dum_str = regexprep(dum_str,['(\W)' all_names{i_var} '(\W)'],['$1dum_var(' num2str(i_var) ')$2']);
 end
 
 try 
  sym_eq = str2syms(dum_str); % convert equation to symbolic expression
  sym_eq = expand(sym_eq); % get rid of brackets
  sym_var = symvar(sym_eq); % get single variables in expression
  
  try
   sym_coef = coeffs(sym_eq,sym_var); % get coefficients of all polynomial terms
  catch % try to extract coeffs
   if ~strcmp(varflag,'i')
    fprintf('Nonpolynomial nonlinearity detected. Use collocation method only!\n');
   end
  end
  
  % if this is possible, nonpolynomial nonlinearities are considered as 
  % coefficients of polynomial terms; if not, error will be caught
  try
   % try, if all cofficients are numeric <=> expression contains no nonpolynomial nonlinearities 
   double(sym_coef); % if all coefficients are convertible to double, expression is polynomial
   is_poly = true;
  catch % try to convert coeffs to double
   if ~strcmp(varflag,'i')
    fprintf('Nonpolynomial nonlinearity detected. Use collocation method only!\n');
   end
  end
 catch % try to convert expression to symbolic
  if ~strcmp(varflag,'i')
   fprintf('Nonpolynomial nonlinearity detected. Use collocation method only!\n');
  end
 end
end
end