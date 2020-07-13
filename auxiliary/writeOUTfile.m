function writeOUTfile(PoCETsys,filename)
% writeOUTfile(PoCETsys, 'filename') 
% writes a function-file filename.m containing the outputs' expanded
% right-hand-sides of the system defined in PoCETsys. The file will be
% created in the current directory. This function-file is required for a
% simulation using PoCETsimGalerkin.

states = {PoCETsys.states.name};
parameters = {PoCETsys.parameters.name};
inputs = {PoCETsys.inputs.name};
outputs = {PoCETsys.outputs.name};
writeIC   = '\n ';
writeKRON = '\n ';
writeIFUN = '\n ';
writeIPAR = '';
writeODE  = '\n ';
writeFOUT = '\n fout = [';
terms  = {};

if length(filename)>1
 if ~strcmp(filename(end-1:end),'.m'), 
  filename = [filename '.m']; 
 end
else
 filename = [filename '.m'];
end

fprintf('Writing output-file ''%s'' ...\n', filename);

%% generate input functions
for i = 1:numel(states)
 writeIC  = [writeIC states{i} ' = X(' num2str(i-1) '*PoCETsys.pce.options.n_phi+1:' ...
             num2str(i) '*PoCETsys.pce.options.n_phi);\n '];
end


for i=1:numel(outputs)
 writeODE = [writeODE outputs{i} ' = '];

 for j=1:numel(PoCETsys.out(i).terms)
  tmpTERM = {};
  tmpPARAM = {};
  tmpINPUT = [];

  %% states
  [~,nIDX] = find(PoCETsys.out(i).terms(j).state_index);
  tmpIDX = PoCETsys.out(i).terms(j).state_index(nIDX);
  tmpDEG = PoCETsys.out(i).terms(j).state_degrees(nIDX);

  [~,zIDX] = find(tmpIDX);
  tmpIDX = tmpIDX(zIDX);
  tmpDEG = tmpDEG(zIDX);

  [~,sIDX] = sort(tmpIDX);
  tmpIDX = tmpIDX(sIDX);
  tmpDEG = tmpDEG(sIDX);

  tmpXI   = [];

  for k=1:numel(tmpIDX)
   tmpXI = [tmpXI repmat(tmpIDX(k),1,tmpDEG(k))]; 
  end

  if ~isempty(tmpXI)
   tmpTERM = [tmpTERM {[states{tmpXI}]; tmpXI}];
  else
   tmpTERM = [tmpTERM {''; tmpXI}];
  end

  terms = [terms tmpTERM];

  %% parameters
  [~,nIDX] = find(PoCETsys.out(i).terms(j).param_index);
  tmpIDX = PoCETsys.out(i).terms(j).param_index(nIDX);
  tmpDEG = PoCETsys.out(i).terms(j).param_degrees(nIDX);

  [~,sIDX] = sort(tmpIDX);
  tmpIDX = tmpIDX(sIDX);
  tmpDEG = tmpDEG(sIDX);

  tmpORD = sum(PoCETsys.out(i).terms(j).state_degrees);
  tmpXI   = [];

  for k=1:numel(tmpIDX)
   tmpXI = [tmpXI repmat(tmpIDX(k),1,tmpDEG(k))]; 
  end

  if ~isempty(tmpXI)
   tmpPARAM = [tmpPARAM {['M.' parameters{tmpXI} '_O' num2str(tmpORD)]; tmpXI}];
  else
   tmpPARAM = [tmpPARAM {['M.one_O' num2str(tmpORD)]; tmpXI}];
  end

  %% inputs
  [~,nIDX] = find(PoCETsys.out(i).terms(j).input_index);
  tmpIDX = PoCETsys.out(i).terms(j).input_index(nIDX);
  tmpDEG = PoCETsys.out(i).terms(j).input_degrees(nIDX);

  [~,sIDX] = sort(tmpIDX);
  tmpIDX = tmpIDX(sIDX);
  tmpDEG = tmpDEG(sIDX);

  for k=1:numel(tmpIDX)
   tmpINPUT = [tmpINPUT inputs{tmpIDX(k)} '^' num2str(tmpDEG(k)) '*'];
  end

  %% add current term to ODE-file
  writeODE = [writeODE num2str(PoCETsys.out(i).terms(j).coefficient) ...
              '*' tmpINPUT tmpPARAM{1,1} '*' tmpTERM{1,1} ' + '];
 end
 writeODE(end-2:end+1) = ';\n ';

 writeFOUT = [writeFOUT outputs{i} '; '];
end

writeFOUT(end-1:end+2) = '];\n';

writeODE = strrep(writeODE,'^1','');
writeODE = strrep(writeODE,' 1**',' ');
writeODE = strrep(writeODE,' 1*',' ');
writeODE = strrep(writeODE,' -1*',' - ');
writeODE = strrep(writeODE,'**','*');
writeODE = strrep(writeODE,'*;',';');
writeODE = strrep(writeODE,' ;',';');
writeODE = strrep(writeODE,'* +',' +');
writeODE = strrep(writeODE,'+ -','-');

%% generate input parameters & functions
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

writeIFUN = strrep(writeIFUN,' -1*',' - ');
writeIFUN = strrep(writeIFUN,' 1*',' ');
writeIFUN = strrep(writeIFUN,'+ -','- ');
writeIFUN = strrep(writeIFUN,' ;',';');
writeIFUN = strrep(writeIFUN,'  ',' ');

%% generate kronecker products
[~,uIDX] = unique(terms(1,:)); % find unique state-multiplication combinations
uterms = terms(:,uIDX);
if strcmp(uterms{1,1},'')
 uterms = uterms(:,2:end);
end

for k=1:size(uterms,2)
 if numel(uterms{2,k}) == 2
  writeKRON = [writeKRON uterms{1,k} ...
               ' = ckron(' states{uterms{2,k}(1)} ...
                        ',' states{uterms{2,k}(2)} ');\n '];
 elseif numel(uterms{2,k}) > 2
  tmpKRON1 = [];
  tmpKRON2 = [];
  for l=1:numel(uterms{2,k})-1
   tmpKRON1 = [tmpKRON1 'ckron(' states{uterms{2,k}(l)} ','];
   tmpKRON2 = [tmpKRON2 ')'];
  end
  writeKRON = [writeKRON uterms{1,k} ' = ' tmpKRON1 ...
               states{uterms{2,k}(end)} tmpKRON2 ';\n '];
 end
end

fid = fopen(filename,'w');

PoCETode = ...
 sprintf('function fout = %s(t,X,PoCETsys%s)\n M = PoCETsys.coeff_matrices;\n%send',...
  filename(1:end-2),writeIPAR,...
  [writeIC writeIFUN writeKRON writeODE writeFOUT]);

fprintf(fid,PoCETode);

fclose(fid);

end