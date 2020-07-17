function CheckMsg = checkInput(PoCETsys)
% PoCET internal function. Do not call directly!

vv = {'states','parameters','outputs','inputs'};
CheckMsg = '';
nvars1 = numel(PoCETsys.states) + numel(PoCETsys.parameters) + numel(PoCETsys.outputs) + numel(PoCETsys.inputs);
nvars2 = numel(unique({PoCETsys.states.name PoCETsys.parameters.name PoCETsys.outputs.name PoCETsys.inputs.name}));
if nvars1 ~= nvars2 % double variable names?
 %% check for double names within each catergory
 for vi = 1:numel(vv)
  varsi = vv{vi};
  if numel(unique({PoCETsys.(varsi).name})) ~= numel(PoCETsys.(varsi))
   [n,~,i2] = unique({PoCETsys.(varsi).name});
   u_id = accumarray(i2(:),(1:length(i2))',[],@(x) {sort(x)});
   for j = 1:numel(u_id)
    if numel(u_id{j})>1
     CheckMsg = sprintf([CheckMsg '%s%s %i and %i are both named ''%s''!\n'],upper(varsi(1)),varsi(2:end),u_id{j}(1),u_id{j}(2),n{j});
    end
   end
  end
  for vj = vi+1:numel(vv)
   %% check for double names in different categories
   varsj = vv{vj};
   [n,i1,i2] = intersect({PoCETsys.(varsi).name}, {PoCETsys.(varsj).name});
   if ~isempty(n) && ~isempty(n{:})
    for i = 1:numel(n)
     CheckMsg = sprintf([CheckMsg '%s%s %i and %s %i are both named ''%s''!\n'],upper(varsi(1)),varsi(2:end-1),i1(i),varsj(1:end-1),i2(i),n{i});
    end
   end
  end
 end
end

%% check if data matches distribution
for vi = 1:2
 vars = vv{vi};
 for i=1:numel(PoCETsys.(vars)),
  switch lower(PoCETsys.(vars)(i).dist),
   case 'uniform'
    if numel(PoCETsys.(vars)(i).data)~=2
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s does not match distribution!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    elseif PoCETsys.(vars)(i).data(1) > PoCETsys.(vars)(i).data(2)
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s is incorrect: The lower bound has to be smaller than the upper bound!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    end
   case 'normal'
    if numel(PoCETsys.(vars)(i).data)~=2
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s does not match distribution!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    elseif PoCETsys.(vars)(i).data(2) <= 0
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s is incorrect: The standard deviation has to be greater than zero!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    end
   case 'beta'
     if numel(PoCETsys.(vars)(i).data)~=2
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s does not match distribution!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    elseif PoCETsys.(vars)(i).data(1)<0 || PoCETsys.(vars)(i).data(2)<0
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s is incorrect - alpha and beta have to be greater than zero!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    end
   case 'beta4'
    if numel(PoCETsys.(vars)(i).data)~=4
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s does not match distribution!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    elseif PoCETsys.(vars)(i).data(1)<0 || PoCETsys.(vars)(i).data(2)<0
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s is incorrect - alpha and beta have to be greater than zero!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    elseif PoCETsys.(vars)(i).data(3) > PoCETsys.(vars)(i).data(4)
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s is incorrect: The lower bound has to be smaller than the upper bound!\n'],vars(1:end-1),PoCETsys.(vars)(i).name);
    end
   case 'pce'
    if numel(PoCETsys.(vars)(i).data) ~= PoCETsys.pce.options.n_phi;
     CheckMsg = sprintf([CheckMsg 'Data input of %s %s is incorrect: Number of PCE coefficients does not match size (%i) of PCE basis.\n'],vars(1:end-1),PoCETsys.(vars)(i).name,PoCETsys.pce.options.n_phi);
    end
  end
 end
end
end