function [analysDir,toolsDir,dataDir,outDir] = f_addDirectories(analysFldr,toolsFldr,dataFlrd,outFlrd)
%% Add all the directories of the algoritm
% If they don't exist, they get created
dirabove = 1; % All the folders are at maximum one above the current one
analysDir = f_makeParentFolder(dirabove,analysFldr); % Adds the directory
% analysDir = pwd; % Same functionality
toolsDir = f_makeParentFolder(dirabove,toolsFldr); % Adds the directory
dataDir = f_makeParentFolder(dirabove,dataFlrd); % Adds the directory
outDir = f_makeParentFolder(dirabove,outFlrd); % Adds the directory
addpath(genpath(toolsDir)); % Add all folders inside the tools folder
% Restore back default paths, type: restoredefaultpath
end