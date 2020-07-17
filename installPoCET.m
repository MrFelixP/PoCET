function installPoCET(mode,varargin)
% installPoCET This is an "install" function that checks 
% consistency of the toolbox and adds needed path entries. 
% 
% To only add paths, run the script using the option 'path'.  
% To remove path definitions, use 'removePaths'.
% 
% installPoCET()
% installPoCET('paths')
% installPoCET('removePaths')
%
% Note this script smust be ran after an update of the toolbox.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear global workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
evalin('base','clear all');
evalin('base','clear classes');
evalin('base','clear all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fullInstall = 1;
pathInstall = 0;
webDocuInstall = 0;
docuInstall = 0;
if nargin == 0,
elseif nargin == 1,
    mode = upper(mode);
    if strcmpi(mode,'paths') || strcmpi(mode,'path'),
        fullInstall = 0;
        pathInstall = 1;
    elseif strcmpi(mode,'removePaths') || strcmpi(mode,'removePath'),
        removePathDefs();    
        return;
    else
        error('Incorrect input arguments.');
    end        
else
    error('Incorrect input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set license texts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
licenseText1 = {...
    'PoCET (Polynomial Chaos Expansion Toolbox) for MATLAB',...
    '  by [Petzke et al., 2020]'};
licenseText2 = {...
    '*  See also "https://www.tu-chemnitz.de/etit/control/research/PoCET/" where',...
    '   you can always download the newest version and can get support.',...
    '*  Please give credits to the authors of the toolbox and cite the reference',...
    '   given above in your talks and published work.'};
licenseText3 = {'Copyright (c) 2020 Stefan Streif, stefan.streif@etit.tu-chemnitz.de', ...
    'Licensed under the EUPL-1.2-or-later'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up desktop and show message
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
curDate=clock;
diary off
diary(['installPoCET_' datestr(curDate,'dd-mm-yy_HH-MM') '.log']);
disp('+--------------------------------------------------------------------------------+');
disp('| INSTALLATION                                                                   |');
disp('+--------------------------------------------------------------------------------+');
for i = 1:length(licenseText1),
    fprintf('| %-78s |\n',licenseText1{i});
end
fprintf('+--------------------------------------------------------------------------------+\n');
for i = 1:length(licenseText3),
    fprintf('| %-78s |\n',licenseText3{i});
end
fprintf('+--------------------------------------------------------------------------------+\n');
q = input('Do you accept the licensing conditions? (y/n [n]): ','s');
if ~strcmpi('Y',q),
    disp('Aborting installation.');
    return
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set path for the toolbox components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% search for existing installations of the toolbox
if fullInstall || pathInstall,
    disp('+--------------------------------------------------------------------------------+');
    % check for previous installations
    disp('Searching for previous installation of PoCET...');
    removePathDefs();    
    % add new paths
    disp('Adding search paths for PoCET functions...');
    neededPaths = {'auxiliary'};    
    if checkIsPoCETRootFolder(neededPaths),
        % Apply path settings and save path
        PATH_TOOLBOX = pwd;
        warning off;
        addpath(PATH_TOOLBOX);
        for needPath = neededPaths,
            if ispc,
                addpath(genpath([PATH_TOOLBOX,'\',needPath{1}]));
            else %ISUNIX,
                addpath(genpath([PATH_TOOLBOX,'/',needPath{1}]));
            end
        end
        result = savepath;
        warning on;
        if result == 1,
            disp(' ');
            disp('Your MATLAB installation does not allow the saving of the updated MATLAB path variable.');
            disp('This means that you have to add the PoCET folder and all its subfolders to the MATLAB');
            disp('path each time you start MATLAB or run this script using the option installPoCET(''path''),');
            disp('which does that without testing external tools that might be required or other toolboxes.');        
        end
    else
        fprintf(2,'Error! Necessary folders are missing, please make sure you unpacked the toolbox\npreserving the structure and run the script from the toolbox root folder.\n');
        return;
    end
    disp(' ');
%     input('Press any key to continue. ','s');
    disp('Press any key to continue.')
    pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for prerequisite toolboxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fullInstall,
    disp('+--------------------------------------------------------------------------------+');
    pathProblems = 0;
    try
        requiredNotFound = 0;
        notFound = 0;
        requiredToolbox(1).name    = 'MATLAB';
        requiredToolbox(1).version = [7,10];
        requiredToolbox(1).msg = 'PoCET WILL NOT WORK without it, aborting installation.';
        requiredToolbox(1).requiredNotFound = 1;
        %
        requiredToolbox(2).name    = 'Symbolic Math Toolbox';
        requiredToolbox(2).version = [5,0];
        requiredToolbox(2).msg = 'PoCET WILL NOT WORK without it, aborting installation.';
        requiredToolbox(2).requiredNotFound = 1;
        %
        requiredToolbox(3).name    = 'Optimization Toolbox';
        requiredToolbox(3).version = [5,0];
        requiredToolbox(3).msg = 'LIMITED FUNCTIONALITY of PoCET: cannot run examples "ex3_discrimination" and "ex4_fault_detection".';
        requiredToolbox(3).requiredNotFound = 0;
        %
        fprintf('Determining versions of installed toolboxes (this might take few seconds) ...\n\n');
        verS = ver;
        for i = 1:length(requiredToolbox),
            fprintf('Checking "%s" installation ... ',requiredToolbox(i).name);
            %%% all other toolboxes
            idx = regexp(requiredToolbox(i).name,{verS(:).Name});
            idx = find(~cellfun(@isempty,idx));
            if isempty(idx),
                fprintf(2,'\b NOT FOUND or version TOO OLD!\n%s\n',requiredToolbox(i).msg);
                notFound = 1;
                requiredNotFound = requiredNotFound | requiredToolbox(i).requiredNotFound;
                if requiredNotFound,
                    error(' ');
                end
            else
                installedVer = zeros(0,2);
                for j = 1:length(idx),
                    if regexp(verS(idx(j)).Version,'\.'),                            
                        tmp = regexp(verS(idx(j)).Version,'\.','split');
                        installedVer(j,1:length(tmp)) = cellfun(@(x) str2double(x), tmp);
                    else
                        installedVer(j,1) = str2double(verS(idx(j)).Version);
                        installedVer(j,2) = 0;
                    end
                end
                if size(installedVer,1) > 1,
                    [installedVer,tmp] = sortrows(installedVer);
                    installedVer = installedVer(end,:);
                    idx = idx(tmp(end));
                end        
                % note: currently we only check major and first minor
                % version number, i.e. if installed-version = 6.2.1,
                % then we will only consider 6.2
                if any(isnan(installedVer)) || ...
                   ((installedVer(1) < requiredToolbox(i).version(1)) || ...
                    ((installedVer(1) == requiredToolbox(i).version(1)) && ...
                     (installedVer(2) < requiredToolbox(i).version(2)))),
                    fprintf(2,'\b Version TOO OLD or UNKNOWN!\n%s\n',requiredToolbox(i).msg);
                    requiredNotFound = requiredNotFound | requiredToolbox(i).requiredNotFound;
                    notFound = 1;
                    if requiredNotFound,
                        error(' ');
                    end
                else
                    fprintf('\b OK! Installed version: %s\n',verS(idx).Version);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if notFound,
            fprintf(2,'\nNote: PoCET can be used, but some functions for advanced tasks will not work. \nCheck help texts of the affected functions to find out more.\n');
        end
        disp(' ');
%         input('Press any key to continue. ','s');
        disp('Press any key to continue.')
        pause
    catch
        if pathProblems,
            fprintf(2,'\nUNKNOWN ERROR or PATH NOT FOUND! Make sure to run "installPoCET" from the toolbox root folder. Installation aborted.\n',requiredToolbox(i).msg);
            disp(' ');
            return;
        end
        if ~requiredNotFound,
            fprintf(2,'\nUNKNOWN INSTALLATION OR MATLAB PROBLEMS! Installation aborted.\n',requiredToolbox(i).msg);
            disp(' ');
            return;
        end
        disp(' ');
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating Documentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fullInstall || docuInstall || webDocuInstall,
    if fullInstall || docuInstall,
        disp('+--------------------------------------------------------------------------------+');
        disp('Creating documentation...')
        if ispc
            if ~exist('documentation', 'dir')
                mkdir('documentation')
            end
            addpath([PATH_TOOLBOX,'\documentation']);
            if ~exist('documentation\help', 'dir')
                mkdir('documentation\help')
            end
            addpath([PATH_TOOLBOX,'\documentation\help']);
        else
            if ~exist('documentation', 'dir')
                mkdir('documentation')
            end
            addpath([PATH_TOOLBOX,'/documentation']);
            if ~exist('documentation/help', 'dir')
                mkdir('documentation/help')
            end
            addpath([PATH_TOOLBOX,'/documentation/help']);
        end
         
        warning('off','MATLAB:dispatcher:nameConflict');
        try
            createDocumentation(0);
        catch ME
%             rethrow(ME)
            fprintf(2,'Problems during creation of documentation and help texts. PoCET should work nevertheless.\n');
        end
        warning('on','MATLAB:dispatcher:nameConflict');

%         if ispc
%             builddocsearchdb([PATH_TOOLBOX,'\documentation\help']);
%         else %ISUNIX,
%             builddocsearchdb([PATH_TOOLBOX,'/documentation/help']);
%         end
        fprintf('\nNote: You might have to restart Matlab once in order to \nrefresh help paths and to use the PoCET documentation.\n')
    end
%     if webDocuInstall,
%         disp('+--------------------------------------------------------------------------------+');
%         disp('Creating documentation for homepage...')
%         warning('off','MATLAB:dispatcher:nameConflict');
%         createDocumentation(1);
%         warning('on','MATLAB:dispatcher:nameConflict');
%     end
    disp(' ');
%     input('Press any key to continue. ','s');
    disp('Press any key to continue.')
    pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output license information, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diary off;
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('+--------------------------------------------------------------------------------+\n');
for i = 1:length(licenseText1),
    fprintf('| %-78s |\n',licenseText1{i});
end
fprintf('+--------------------------------------------------------------------------------+\n');
for i = 1:length(licenseText2),
    fprintf('| %-78s |\n',licenseText2{i});
end
fprintf('+--------------------------------------------------------------------------------+\n');
for i = 1:length(licenseText3),
    fprintf('| %-78s |\n',licenseText3{i});
end
fprintf('+--------------------------------------------------------------------------------+\n');
fprintf('Installation finished.\n\nFor more information on the installation process check log-file:\n');
if ispc,
    fprintf('\t"%s\\installPoCET_%s.log"\n',pwd,datestr(curDate,'dd-mm-yy_HH-MM'));
else %ISUNIX,
    fprintf('\t"%s/installPoCET_%s.log"\n',pwd,datestr(curDate,'dd-mm-yy_HH-MM'));
end
fprintf('\n\nThe documentation for Matlab 2012b or newer can be accessed via\nHelp->Supplemental Software.\n');
fprintf('\nTo start, please have a look at the examples in \n\t"%s\\examples".\n',pwd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createDocumentation(forWebsite)

% directories where we look for m-files
directories = {'.', 'auxiliary/'};

% the html-files go here
help_dir = 'documentation/help';
toolbox_dir = pwd;

for dir_idx = 1:length(directories)
    directory = directories{dir_idx};
    dir_listing = dir(directory);
    for list_idx = 1:length(dir_listing)
        entry = dir_listing(list_idx);
        if ~entry.isdir && strcmp(entry.name(end-1:end), '.m')
            if forWebsite
                generateSHTMLDoc(entry.name, fullfile(toolbox_dir, help_dir))
            else
                generateHTMLDoc(entry.name, fullfile(toolbox_dir, help_dir))
            end
        end
    end
end

% Create the html file which automatically redirects to the PoCET web site.
redirect_html_file = fopen(fullfile(toolbox_dir, help_dir, '_redirect.html'), 'w');
redirect_html_code = '<html><head><meta http-equiv="refresh" content="0; url=http://ifatwww.et.uni-magdeburg.de/syst/PoCET/"></head><body></body></html>';
fprintf(redirect_html_file, '%s', redirect_html_code);
fclose(redirect_html_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function removePathDefs()
    if ispc,
        pts1 = regexp(path,';','split');
    else %ISUNIX,
        pts1 = regexp(path,':','split');
    end
    %%% ANTON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 'PoCET' IS NOT always used as the name of the folder!
    % Using more sophisticated way to determine previous installations
    PoCETfolder=which('installPoCET','-all');
    for ind = 1:length(PoCETfolder)
        if ispc,
            slashind=strfind(PoCETfolder{ind},'\');
        else %ISUNIX,
            slashind=strfind(PoCETfolder{ind},'/');
        end
        PoCETfolder{ind}=PoCETfolder{ind}(1:slashind(end)-1);
        pts2=strfind(pts1,PoCETfolder{ind});
        pti = cellfun(@(x) ~isempty(x), pts2);
        if any(pti),
            fprintf(2,'Found existing installation of PoCET. Removing paths to existing installation ...\n');
        end
        for i = 1:length(pti),
            if pti(i),
                rmpath(pts1{i});
            end
        end
    end
    % Changing directory to the one including installPoCET
    isInPoCET=strcmp(cd,PoCETfolder);
    if length(PoCETfolder)==1 && isInPoCET==0
        fprintf('Changing current directory to the toolbox root:\n%s\n',PoCETfolder{1});
        cd(PoCETfolder{1});
    elseif length(PoCETfolder)>1 && any(isInPoCET)==0
        fprintf(2,'Multiple choices of PoCET toolbox folder.\nPlease enter the number of the correct choice:\n');
        for ind = 1:length(PoCETfolder)
            fprintf('%d: %s\n',ind,PoCETfolder{ind});
        end
        while 1
            toolboxinp=input(sprintf('Enter the number (1-%d)>> ',length(PoCETfolder)),'s');
            [toolboxnum, status]=str2num(toolboxinp);
            if status && length(toolboxnum)==1 && toolboxnum>=1 && toolboxnum<=length(PoCETfolder)
                fprintf('Changing current directory to the toolbox root:\n%s\n',PoCETfolder{toolboxnum});
                cd(PoCETfolder{toolboxnum});
                break
            else
                fprintf(2,'Wrong input! Please try again...\n');
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = checkIsPoCETRootFolder(neededPaths)
    try
        if nargin == 0,    
            neededPaths = {'auxiliary'};
%             neededPaths = {'auxiliary','documentation'};
        end
        checkedNeededPaths = cellfun(@(x) length(dir(x))>=1, neededPaths);
        ok = all(checkedNeededPaths);
    catch
        fprintf(2,'Error checking PoCET root folder.\n');
    end
end
