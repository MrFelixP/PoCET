function CheckMsg = checkUpdate(PoCETsystem, update_vars, update_vals)
% PoCET internal function. Do not call directly!

CheckMsg = [];
system_states = {PoCETsystem.states(:).name};
system_params = {PoCETsystem.parameters(:).name};
for i = 1:numel(update_vars)
 if ismember(update_vars{i},system_params)
  idx = find(strcmp(update_vars{i},system_params));
  if strcmp(PoCETsystem.parameters(idx).dist,'uniform')
   if numel(update_vals{i})~=2
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed: Data does not match distribution.\n'];
   elseif update_vals{i}(1) > update_vals{i}(2)
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed: The lower bound has to be smaller than the upper bound!\n'];
   end
  elseif strcmp(PoCETsystem.parameters(idx).dist,'normal')
   if numel(update_vals{i})~=2
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed: Data does not match distribution. \n'];
   elseif update_vals{i}(2) <= 0
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed: The standard deviation cannot be smaller or equal to zero!\n'];
   end
  elseif strcmp(PoCETsystem.parameters(idx).dist,'beta')
    if numel(update_vals{i})~=2
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed: Data does not match distribution. \n'];
   elseif update_vals{i}(1)<0 || update_vals{i}(2)<0
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed - alpha and beta have to be bigger than zero!\n'];
   end
  elseif strcmp(PoCETsystem.parameters(idx).dist,'beta4')
   if numel(update_vals{i})~=4
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed: Data does not match distribution.\n'];
   elseif update_vals{i}(1)<0 || update_vals{i}(2)<0
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed - alpha and beta have to be bigger than zero!\n'];
   elseif update_vals{i}(3) > update_vals{i}(4)
    CheckMsg = [CheckMsg 'Data update of parameter ' update_vars{i} ...
                   ' failed: The lower bound has to be smaller than the upper bound!\n'];
   end
  end
 elseif ismember(update_vars{i},system_states)
  idx = find(strcmp(update_vars{i},system_states));
  if strcmp(PoCETsystem.states(idx).dist,'uniform')
   if numel(update_vals{i})~=2
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed: Data does not match distribution.\n'];
   elseif update_vals{i}(1) > update_vals{i}(2)
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed: The lower bound needs to be smaller than the upper bound!\n'];
   end
  elseif strcmp(PoCETsystem.states(idx).dist,'normal')
   if numel(update_vals{i})~=2
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed: Data does not match distribution. \n'];
   elseif update_vals{i}(2) <= 0
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed: The standard deviation cannot be smaller or equal to zero!\n'];
   end
  elseif strcmp(PoCETsystem.states(idx).dist,'beta')
    if numel(update_vals{i})~=2
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed: Data does not match distribution. \n'];
   elseif update_vals{i}(1)<0 || update_vals{i}(2)<0
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed - alpha and beta need to be bigger than zero!\n'];
   end
  elseif strcmp(PoCETsystem.states(idx).dist,'beta4')
   if numel(update_vals{i})~=4
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed: Data does not match distribution.\n'];
   elseif update_vals{i}(1)<0 || update_vals{i}(2)<0
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed - alpha and beta need to be bigger than zero!\n'];
   elseif update_vals{i}(3) > update_vals{i}(4)
    CheckMsg = [CheckMsg 'Data update of state ' update_vars{i} ...
                   ' failed: The lower bound needs to be smaller than the upper bound!\n'];
   end
  end
 else
  CheckMsg = [CheckMsg 'Unknown variable ' update_vars{i} '.\n']; 
 end
end

if ~isempty(CheckMsg)
 CheckMsg = sprintf(CheckMsg);
end

end