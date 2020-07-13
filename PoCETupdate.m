function [PoCETsys] = PoCETupdate(PoCETsys, varargin)
% PoCETsys = PoCETupdate(PoCETsys, 'var1', data1,..., 'varN', dataN)
% returns the system PoCETsys with updated/modified uncertainties of the 
% initial conditions or parameter values as specified by the input 
% arguments. The input values data1 to dataN have to match the distribution
% of the corresponding variable in PoCETsys.
% 
% NOTE: PoCETupdate cannot change the system's dynamics, the order of the 
%       polynomial chaos expansion, or the disbribution of a variable.
%       Use PoCETcompose instead!

update_vars = varargin(1:2:end);
update_vals = varargin(2:2:end);

if numel(update_vars)~=numel(update_vals)
 error('Wrong specification of update variables. Enter ''help PoCETupdate'' for more details.')
end
if ~all(cellfun(@(x)ischar(x),update_vars))
 error('Wrong specification of update variables. Enter ''help PoCETupdate'' for more details.')
end

system_states = {PoCETsys.states(:).name};
system_params = {PoCETsys.parameters(:).name};

% fprintf('Analyzing input...\n')
CheckMsg = checkUpdate(PoCETsys,update_vars,update_vals);
error(CheckMsg)

for i = 1:numel(update_vars)
 if ismember(update_vars{i},system_states)
  PoCETsys.states(strcmp(update_vars{i},system_states)).data = update_vals{i};
 elseif ismember(update_vars{i},system_params)
  PoCETsys.parameters(strcmp(update_vars{i},system_params)).data = update_vals{i};
 else
  fprintf('Unknown input variable %s.',update_vars{i})
 end
end

% PoCETinputProcessing % evaluate moments (actually nowhere used!)

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
     if strcmp(vars,'states')
      error('Size of PCE coefficients specified for state %s does not match size (%i) of PCE basis.',PoCETsys.(vars)(i).name,PoCETsys.pce.options.n_phi)
     else
      error('Size of PCE coefficients specified for parameter %s does not match size (%i) of PCE basis.',PoCETsys.(vars)(i).name,PoCETsys.pce.options.n_phi)
     end
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
    if strcmp(vars,'states')
     error('Invalid distribution ''%s'' for state %s.',PoCETsys.(vars)(i).dist,PoCETsys.(vars)(i).name)
    else
     error('Invalid distribution ''%s'' for parameter %s.',PoCETsys.(vars)(i).dist,PoCETsys.(vars)(i).name)
    end
  end
 end
end

if any(ismember(update_vars,system_params))
 matnames = fieldnames(PoCETsys.coeff_matrices.raw);
 for i = 1:numel(matnames)
  C = zeros(size(PoCETsys.coeff_matrices.raw.(matnames{i}){3,1}));
  for j = 1:size(PoCETsys.coeff_matrices.raw.(matnames{i}),2)
   pars   = PoCETsys.coeff_matrices.raw.(matnames{i}){1,j};
   coefs  = PoCETsys.coeff_matrices.raw.(matnames{i}){2,j};
   tmpMAT = PoCETsys.coeff_matrices.raw.(matnames{i}){3,j};
   tmpCOEF = 1;
   for k = 1:numel(pars)
    tmpCOEF = tmpCOEF * PoCETsys.parameters(pars(k)).pce(coefs(k));
   end
   C = C + tmpCOEF * tmpMAT;
  end
  PoCETsys.coeff_matrices.(matnames{i}) = C;
 end
end
%% output
fprintf('PoCET finished updating the expanded system.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function PoCETinputProcessing
% 
% for ii = 1:numel(PoCETsys.states)
%     if ~isempty(PoCETsys.states(ii).data);
%         d = PoCETsys.states(ii).data;
%         if strcmp(PoCETsys.states(ii).dist,'normal')
%             PoCETsys.states(ii).moments = [d(1) d(2)];
%         elseif strcmp(PoCETsys.states(ii).dist,'uniform')
%             PoCETsys.states(ii).moments = [(d(1)+d(2))/2, ...
%                 (d(2)-d(1))/sqrt(12)];
%         elseif strcmp(PoCETsys.states(ii).dist,'beta')
%             PoCETsys.states(ii).moments = [d(1)/(d(1)+d(2)), ...
%                 d(1)*d(2)/((d(2)+d(1))^2*(d(2)+d(1)+1))];
%         elseif strcmp(PoCETsys.states(ii).dist,'beta4')
%             PoCETsys.states(ii).moments = [(d(1)*d(4)+d(2)*d(3))/(d(1)+d(2)), ...
%                 d(1)*d(2)*(d(4)-d(3))^2/((d(2)+d(1))^2*(d(2)+d(1)+1))];
%         elseif strcmp(PoCETsys.states(ii).dist,'none')
%             PoCETsys.states(ii).moments = d(1);
%         end
%     end
% end
% 
% for ii = 1:numel(PoCETsys.parameters)
%     if ~isempty(PoCETsys.parameters(ii).data);
%         d = PoCETsys.parameters(ii).data;
%         if strcmp(PoCETsys.parameters(ii).dist,'normal')
%             PoCETsys.parameters(ii).moments = [d(1) d(2)];
%         elseif strcmp(PoCETsys.parameters(ii).dist,'uniform')
%             PoCETsys.parameters(ii).moments = [(d(1)+d(2))/2, ...
%                 (d(2)-d(1))/sqrt(12)];
%         elseif strcmp(PoCETsys.parameters(ii).dist,'beta')
%             PoCETsys.parameters(ii).moments = [d(1)/(d(1)+d(2)), ...
%                 d(1)*d(2)/((d(2)+d(1))^2*(d(2)+d(1)+1))];
%         elseif strcmp(PoCETsys.parameters(ii).dist,'beta4')
%             PoCETsys.parameters(ii).moments = [(d(1)*d(4)+d(2)*d(3))/(d(1)+d(2)), ...
%                 d(1)*d(2)*(d(4)-d(3))^2/((d(2)+d(1))^2*(d(2)+d(1)+1))];
%         elseif strcmp(PoCETsys.parameters(ii).dist,'none')
%             PoCETsys.parameters(ii).moments = d(1);
%         end
%     end
% end
% 
% % for i = 1:numel(inputs)
% %     if strcmpi(inputs(i).type,'piecewiseC')
% %         tmpTIME = inputs(i).time;
% %         tmpVALS = [inputs(i).vals(1) diff(inputs(i).vals)];
% %         if any(diff(tmpTIME)<=0);
% %          errormsg = sprintf('Time-specification of input %s is incorrect: Time-points have to be distinct and non-decreasing!',inputs(i).name);
% %          error(errormsg);
% %         end
% %         ifun = '';
% %         
% %         for k = 1:numel(tmpTIME)
% %          ifun = [ifun num2str(tmpVALS(k)) '*(t>=' num2str(tmpTIME(k)) ') + '];
% %         end
% %         ifun = ifun(1:end-2);
% %         PoCETsys.inputs(i).func = ifun;
% % %         tmpTIME = inputs(i).time;
% % %         mindiff = min(diff(inputs(i).time));
% % %         if tmpTIME(1) == odeoptions.tspan(1); tmpTIME = tmpTIME(2:end); end
% % %         if tmpTIME(end) == odeoptions.tspan(end); tmpTIME = tmpTIME(1:end-1); end
% % %         addvec = repmat([-mindiff*1e-2 mindiff*1e-2],1,numel(tmpTIME));
% % %         tmpTIME = [tmpTIME; tmpTIME];
% % %         tmpTIME = [odeoptions.tspan(1) tmpTIME(:)'+addvec odeoptions.tspan(end)];
% % %         PoCETsys.inputs(i).time = tmpTIME;
% %     end
% % end
% end

end