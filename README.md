# PoCET
A Polynomial Chaos Expansion Toolbox for MATLAB created by the 
[Automatic Control & System Dynamics Lab](https://www.tu-chemnitz.de/etit/control/index.php.en "ACSD Lab") 
at Chemnitz University of Technology.

PoCET was designed for the automatic generation of polynomial chaos expansions (PCE) for linear and nonlinear 
dynamic systems with time-invariant stochastic parameters or initial conditions. It features
<ul>
  <li>a simple syntax and usability,</li>
  <li>built-in handling of Gaussian, uniform, and beta probability density functions,</li>
  <li>projection and collocation-based calculation of PCE coefficients,</li>
  <li>routines for the calculation of stochastic moments from PCE coefficients,</li>
  <li>several routines for system simulation and visualization of results,</li>
  <li>as well as a variety of introductory and instructive examples.</li>
</ul>

Some use-cases include
<ul>
  <li>uncertainty analysis</li>
  <li>parameter estimation</li>
  <li>state estimation and prediction</li>
  <li>optimal experimental design</li>
  <li>(active) fault diagnosis</li>
  <li>stochastic nonlinear MPC</li>
</ul>

PoCET was officially released and introduced on the 21st IFAC World Congress 2020. 
A preprint of the corresponding paper, which includes more detailed information about its usage and features, 
is [available on arXiv](https://arxiv.org/abs/2007.05245 "arXiv").

For more information, see <https://www.tu-chemnitz.de/etit/control/research/PoCET/>.

When publishing results gained by using PoCET use the following citation:
F. Petzke, A. Mesbah, and S. Streif. PoCET: a Polynomial Chaos Expansion Toolbox for Matlab. In Proc. 21st IFAC World Congress, 2020.