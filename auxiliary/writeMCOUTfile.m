function writeMCOUTfile(PoCETsys, filename)
% writeMCOUTfile(PoCETsys, 'filename') 
% writes a function-file filename.m containing the outputs' scalar
% right-hand-sides of the system defined in PoCETsys. The file will be
% created in the current directory. This function-file is required for a
% simulation using PoCETsimCollocation or PoCETsimMonteCarlo.

states = {PoCETsys.states.name};
params = {PoCETsys.parameters.name};
inputs = {PoCETsys.inputs.name};
outputs = {PoCETsys.outputs.name};
writeIC   = '\n ';
writeIFUN = '\n ';
writeIPAR = '';
writeODE  = '\n ';
writeFOUT = '\n fout = [';

if length(filename)>1
 if ~strcmp(filename(end-1:end),'.m')
  filename = [filename '.m'];
 end
else
 filename = [filename '.m'];
end

fprintf('Writing MCODE file ''%s'' ...\n', filename);

for i=1:numel(states)
 writeIC  = [writeIC states{i} ' = X_in(' num2str(i) ',:);\n '];
end

for i=1:numel(outputs)
 writeODE = [writeODE outputs{i} ' = ' PoCETsys.outputs(i).rhs ';\n '];
 writeFOUT = [writeFOUT outputs{i} '; '];
end

for i = 1:numel(params)
 writeODE  = regexprep(writeODE ,['\<' params{i} '\>'],['PAR.' params{i}]);
end

writeFOUT(end-1:end+2) = '];\n';

writeODE = strrep(writeODE,'^','.^');
writeODE = strrep(writeODE,'..','.');
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

PoCETmcout = ...
    sprintf('function fout = %s(t,X_in,PAR%s)\n%send',...
             filename(1:end-2), writeIPAR, ...
             [writeIC writeIFUN writeODE writeFOUT]);

fprintf(fid,PoCETmcout);

fclose(fid);

end