function out = PoCETcalcMoments(PoCETsys, MomMats, pcvals, type, M)
% moments = PoCETcalcMoments(PoCETsys, MomMats, pcvals) 
% calculates the central moments of a variable from its PCE
% coefficients pcvals which are returned by PoCETsimGalerkin or
% PoCETsimCollocation, respectively. Only the first M=SIZE(MomMats,1)
% moments are calculated. 
% If pcvals is a matrix of PCE coefficients, PoCETcalcMoments will return a 
% M-by-size(pcvals,2) matrix containing where each row contains the 
% moment of corresponding order (1st row -> 1st central moment, ...).
% If pcvals is a struct, PoCETcalcMoments returns a struct containing the
% moments of all variables in pcvals.
% 
% moments = PoCETcalcMoments(PoCETsys, MomMats, pcvals, 'raw')
% returnes the raw moments instead.
%
% moments = PoCETcalcMoments(PoCETsys, MomMats, pcvals, ..., M)
% returnes moments up to order M with M <= SIZE(MomMats,1).
%
% Example: 
%  PoCETcalcMoments(PoCETsys, MomMats, results.x1.pcvals(:,end), 'raw', 2) 
%  returns the first two raw moments of x1 at the end of the simulation.

if nargin < 3
 error('Not enough input arguments.');
elseif nargin == 3
 type = 'central';
 M = size(MomMats,1);
elseif nargin == 4
 M = size(MomMats,1);
elseif nargin == 5
 if isempty(type)
  type = 'central';
 end
elseif nargin > 5
 error('Too many input arguments.');
end

if isstruct(pcvals)
 out = pcvals;
 vars = fieldnames(out);
 vars = vars(~strcmp(vars,'time'));
 for i = 1:numel(vars)
  if isfield(out.(vars{i}),'pcvals')
   out.(vars{i}).moments = calcMoments(MomMats,out.(vars{i}).pcvals,type,M);
  else
   fprintf(2,'! WARNING: No valid results for variable %s found!\n',vars{i});
  end
 end
else
 % check dimensions of coefficient matrix
 [r,c] = size(pcvals);
 if r == PoCETsys.pce.options.n_phi
  xi = pcvals;
 elseif c == PoCETsys.pce.options.n_phi
  xi = pcvals';
 else
  error('\nSize of input argument ''pcvals'' is wrong:\nThe size of at least one dimension is expected to equal %i.',PoCETsys.pce.options.n_phi);
 end
 out = calcMoments(MomMats,xi,type,M);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function moments = calcMoments(MomMats,xi,type,M)
steps = size(xi,2);
raw = zeros(M,steps);
raw(1,:) = xi(1,:).*MomMats(1).val;
if M >= 2
 raw(2,:) = sum(xi.^2.*(MomMats(2).val*ones(1,steps)),1);
end
if M >= 3
 for iM = 3:M
  val = MomMats(iM).val;
  phi = MomMats(iM).phi;
  for k = 1:steps
   raw(iM,k) = val'*prod(reshape(xi(phi,k),[size(phi),1]),2);
  end
 end
end

moments = zeros(M,steps);
switch lower(type)
 case 'raw'
  moments = raw;
 case 'central'
  moments(1,:) = raw(1,:); 
  if M >= 2; moments(2,:) = raw(2,:)-raw(1,:).^2; end
  if M >= 3; moments(3,:) = ...
     (2*raw(1,:).^3 - 3*raw(1,:).*raw(2,:) + raw(3,:))./sqrt(raw(2,:)-raw(1,:).^2).^3; end
  if M >= 4; moments(4,:) = ...
     (-3*raw(1,:).^4 + 6*raw(1,:).^2.*raw(2,:) - 4*raw(1,:).*raw(3,:) + raw(4,:))./(raw(2,:)-raw(1,:).^2).^2 - 3; end
end
end