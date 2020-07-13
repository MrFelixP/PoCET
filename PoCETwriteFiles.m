function PoCETwriteFiles(PoCETsys,RHSfile,OUTfile,MCRHSfile,MCOUTfile)
% PoCETwriteFiles(PoCETsys, RHSfile, OUTfile, MCRHSfile, MCOUTfile) 
% writes .m-files containing all dynamic- and output-equations of the 
% Polynomial Chaos Expansion of the system defined in PoCETsys. The files
% will be created in the current directory.
%
% PoCETwriteFiles(PoCETsys) creates all files using default filenames. 
%
% EXAMPLES:
%  PoCETwriteFiles(PoCETsys) 
%  creates files RHS.m, OUT.m, MCRHS.m, and MCOUT.m.
%
%  PoCETwriteFiles(PoCETsys,'myRHS.m',[],'','myMCOUT')
%  creates dynamics-file myRHS.m and output-file myMCOUT.m
%
%  PoCETwriteFiles(PoCETsys,[],'sysOUT')
%  creates output-file sysOUT.m

if nargin == 1
 RHSfile = 'RHS.m'; 
 OUTfile = 'OUT.m'; 
 MCRHSfile = 'MCRHS.m'; 
 MCOUTfile = 'MCOUT.m';
end

if exist('RHSfile','var') && ~isempty(RHSfile) && ischar(RHSfile), 
 writeRHSfile(PoCETsys,RHSfile); 
end

if exist('MCRHSfile','var') && ~isempty(MCRHSfile) && ischar(MCRHSfile), 
 writeMCRHSfile(PoCETsys,MCRHSfile); 
end

if ~isempty(PoCETsys.outputs(1).name)
 if exist('OUTfile','var') && ~isempty(OUTfile) && ischar(OUTfile), 
  writeOUTfile(PoCETsys,OUTfile); 
 end
 if exist('MCOUTfile','var') && ~isempty(MCOUTfile) && ischar(MCOUTfile), 
  writeMCOUTfile(PoCETsys,MCOUTfile); 
 end
end
end