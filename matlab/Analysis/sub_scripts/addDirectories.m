%% Add all the directories of the algoritm
analysDir = pwd; cd ..; % Store script directory
cd(toolsFldr); toolsDir = pwd; cd ..; % Store function directory
cd(dataFlrd); dataDir = pwd;  cd ..; % Store data directory
cd(outFlrd); outDir = pwd; cd ..; % Store output directory
addpath(genpath(toolsDir)); cd(analysDir); % Add all folders in functions
% restore back default paths, type: restoredefaultpath