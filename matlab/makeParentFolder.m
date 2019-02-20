function [folderPath] = makeParentFolder(TGTParentLevel,targetFolder)
%% Retrieve current folder jerarchy
currentDirectory = pwd;
directoryElements = strsplit(currentDirectory,'\');
ResPath = strjoin(directoryElements(1:end-TGTParentLevel),'\');
%% Create folder if not present
folderPath = strcat(ResPath,'\',targetFolder);
if exist(folderPath,'dir') == 0
  mkdir(ResPath,targetFolder);
end
end