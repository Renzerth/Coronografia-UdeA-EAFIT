function [folderPath] = f_makeParentFolder(TGTParentLevel,targetFolder,pathSep)
% Retrieves the path of a specified folder name TGTParetLevel above
% Inputs:
%  -TGTParentLevel: folder to search with a level above the current one, i.e.
%  from [0,1,2,3,...]
%  -targetFolder: folder name that one wants to search or to create. String

%% Retrieve current folder jerarchy
currentDirectory = pwd;
directoryElements = strsplit(currentDirectory,pathSep);
ResPath = strjoin(directoryElements(1:end-TGTParentLevel),pathSep);
%% Create folder if not present
folderPath = strcat(ResPath,pathSep,targetFolder);
if exist(folderPath,'dir') == 0
  mkdir(ResPath,targetFolder);
end
end