% Using the Polynomial Chaos Expansion Toolbox (PoCET)
% ---------------------------------------------------
%% Step 1: defining a system
% There are only four different input types: states, parameters, inputs,
% and outputs. Each type is defined via a struct-array variable. 
%
% States are defined in a struct containing the following fields:
%  .name - String containing the state's name.
%         
%  .dist - String containing the distribution of the initial condition.
%          If the initial condition is not distributed use 'none'.
%          Valid inputs are
%          'none' | 'normal' | 'uniform' | 'beta' | 'beta4' | 'pce'
%
%  .data - Vector containing the distribution's specification. The expected 
%          input depends on the specified distribution:
%          'none'   : [value]
%          'normal' : [mean, standard_deviation]
%          'uniform': [lower_bound, upper_bound]
%          'beta'   : [alpha, beta]
%          'beta4'  : [alpha, beta, lower_bound, upper_bound]
%          'pce'    : array with n_phi entries containing the state's 
%                     pce-values at time 0. Only possible if the state
%                     itself is not random (see example 5).
%
%  .rhs  - String containing the right-hand-side of the state's ODE.
%
%
% Parameters are defined in a struct containing the following fields:
%  .name - String containing the parameter's name. 
%
%  .dist - String containing the distribution of the parameter. If the
%          parameter is not distributed use 'none'. Valid inputs are
%          'none' | 'normal' | 'uniform' | 'beta' | 'beta4' | 'pce'
%
%  .data - Vector containing the distribution's specification. The expected 
%          input depends on the specified distribution:
%          'none'   : [value]
%          'normal' : [mean, standard_deviation]
%          'uniform': [lower_bound, upper_bound]
%          'beta'   : [alpha, beta]
%          'beta4'  : [alpha, beta, lower_bound, upper_bound] 
%          'pce'    : array with n_phi entries containing the parameter's 
%                     pce-values. Only possible if the parameter itself
%                     not random (see example 5).
% 
%
% Inputs are defined in a struct containing the following fields:
%  .name - String containing the input's name. 
%
%  .rhs  - string contatining the right-hand-side of the input equation.
%          It can contain arbitrary scalar MATLAB-functions. 
%          All inputs parameters used in rhs have to be defined in an 
%          additional field (see example 1)
%
%
% Ouputs are defined in a struct containing the following fields:
%  .name - String containing the ouputs's name. It has to start with a
%          letter and can only contain underscores as special characters.
%
%  .rhs  - string contatining the right-hand-side of the ouput equation.
%
%
% Options for solving the ODEs can be set in a struct with fields: 
%  .setup   -  odeset-struct
%
%  .solver 	-  string containing the solver's name; default is 'ode15s'
%
%  .tspan   -  vector [t_0, t_end] or [t_0, t_1, ..., t_end]
%
%  .dt      -  double defining the stepwidth of the solver's output;
%              if tspan is a vector of more than two timepoints, setting
%              dt will overwrite these values
%
% NOTE: The names of all variables have to be unique, start with a
% letter and can only contain letters, numbers, and underscores.
%
%% Step 2: composing the system struct
% Almost every function of this toolbox requires a PoCET-system-struct as an
% input. This structure can be created using the function PoCETcompose with
% the above mentioned structs as input arguments (see 'help PoCETcompose'
% or any example for more information).
%
%% Step 3: simulating the system
% Before a simulation can be done, you have to write all system equations
% to function files using PoCETwriteFiles.
% Afer that, PoCET offers three different ways to simulate a system defined
% by its PoCET-system-struct: 
% 
% PoCETsimGalerkin simulates the system using the Galerkin-method and
% returns the PCE-coefficient of the system's time-dependent variables.
% This method works for polynomial systems ONLY! If you work with arbitrary
% nonlinearities, use collocation method (see example 4) or transform the
% system into a polynomial representation.
% 
% PoCETsimCollocation simulates the system using the stoachastic collocation
% method and returns the PCE-coefficient of the system's time-dependent
% variables.
% 
% PoCETsimMonteCarlo does a Monte-Carlo-simulation of the system and returns
% the trajectories or moments (however specified) of the system's
% time-dependent variables.
% 
% See the individual help files for more information.
%
%% Step 4: visualizing the results
% After retrieving the PCE coefficients of the time-dependent variables
% in Step 3, you can calculate their the moment using the Galerkin method. 
% These two functions will help you with that: 
% 
% PoCETmomentMatrices evaluates the integrals that occur during moment
% calculation. These integrals are independent of the PCE coefficients
% and, hence, only have to be evaluated once for the system.
%
% PoCETcalcMoments calculates the moments of the time depending variables
% from their PCE coefficients and the system's moment matrices. 
%
% See the individual help files for more information.

help PoCET

