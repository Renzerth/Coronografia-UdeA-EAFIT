%[analysDir,toolsDir,dataDir,outDir] = f(analysFldr,toolsFldr,dataFlrd,outFlrd,filemanag,subscripts)
%% Add all the directories of the algoritm
% If they don't exist, they get created
addpath(strcat('..',pathSep,toolsFldr,pathSep,filemanag));
analysDir = f_makeParentFolder(0,analysFldr); % Store script directory
% analysDir = pwd; % Same functionality
toolsDir = f_makeParentFolder(1,toolsFldr); % Tools is 1 folder above
dataDir = f_makeParentFolder(1,dataFlrd); % Store data directory
outDir = f_makeParentFolder(1,outFlrd); % Store output directory
addpath(genpath(toolsDir)); % Add all folders in the functions folder
addpath(subscripts); % Adds all the sub programs of the algorithm
% Restore back default paths, type: restoredefaultpath