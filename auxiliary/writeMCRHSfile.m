function writeMCRHSfile(PoCETsys,filename)
% writeMCRHSfile(PoCETsys, 'filename') 
% writes a function-file filename.m containing the states' scalar
% right-hand-sides of the system defined in PoCETsys. The file will be
% created in the current directory. This function-file is required for a
% simulation using PoCETsimCollocation or PoCETsimMonteCarlo.

states = {PoCETsys.states.name};
params = {PoCETsys.parameters.name};
inputs = {PoCETsys.inputs.name};
writeIC   = '\n ';
writeIFUN = '\n ';
writeIPAR = '';
writeODE  = '\n';
writeDXDT = '\n dXdt = [';

if length(filename)>1
 if ~strcmp(filename(end-1:end),'.m')
  filename = [filename '.m'];
 end
else
 filename = [filename '.m'];
end

fprintf('Writing MCODE file ''%s'' ...\n', filename);

for i=1:numel(states)
 writeIC  = [writeIC states{i} ' = X(' num2str(i) ');\n '];
 writeODE = [writeODE ' ddt_' states{i} ' = ' PoCETsys.states(i).rhs ';\n'];
 writeDXDT = [writeDXDT 'ddt_' states{i} '; '];
end

for i = 1:numel(params)
 writeODE  = regexprep(writeODE ,['\<' params{i} '\>'],['PAR.' params{i}]);
end

writeDXDT(end-1:end+2) = '];\n';

writeODE = strrep(writeODE,'*;',';');
writeODE = strrep(writeODE,' 1**',' ');
writeODE = strrep(writeODE,' 1*',' ');
writeODE = strrep(writeODE,'-1*M','- M');
writeODE = strrep(writeODE,'**','*');
writeODE = strrep(writeODE,'^1','');
writeODE = strrep(writeODE,'+ -','-');
writeODE = strrep(writeODE,'* ',' ');

for i = 1:numel(inputs)
 if ~isempty(inputs{i})
  writeIFUN = [writeIFUN inputs{i} ' = ' PoCETsys.inputs(i).rhs ';\n '];
 end
end

IPAR = fieldnames(PoCETsys.inputs);
for i = 1:numel(IPAR)
 if ~(strcmpi(IPAR{i},'rhs') || strcmpi(IPAR{i},'name'))
  writeIPAR = [writeIPAR ',' IPAR{i}];
 end
end

for i = 1:numel(params)
 writeIFUN = regexprep(writeIFUN,['\<' params{i} '\>'],['PAR.' params{i}]);
end

writeIFUN = strrep(writeIFUN,' -1*',' - ');
writeIFUN = strrep(writeIFUN,' 1*',' ');
writeIFUN = strrep(writeIFUN,'+ -','- ');
writeIFUN = strrep(writeIFUN,' ;',';');
writeIFUN = strrep(writeIFUN,'  ',' ');

fid = fopen(filename,'w');

PoCETmcode = ...
    sprintf('function dXdt = %s(t,X,PAR%s)\n%send',...
             filename(1:end-2), writeIPAR, ...
             [writeIC writeIFUN writeODE writeDXDT]);

fprintf(fid,PoCETmcode);

fclose(fid);

end