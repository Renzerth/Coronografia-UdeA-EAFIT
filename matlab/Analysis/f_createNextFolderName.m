function [numberedFolder] = f_createNextFolderName(dirPath,folderName,dirDelim,pathSep,createFold)
% Creates the next possible folder name for a given "folderName" in a
% specific directory
%
% Inputs:
%  dirPath: path of the folders that one wants to count
%  folderName: name of the root folder
%  dirDelimiter: splits the rootfolder and the distinct character between
%                folders (a string)
%  pathSep: path separation
%  createFolder: creates or not the folder (boolean)
%
%  Output:
%   numberedFolder: with the currentCounter (extended name)
%
currentCounter = f_addFolderCount(dirPath,folderName,dirDelim);
numberedFolder = strcat(folderName,dirDelim,num2str(currentCounter));
if createFold
  if currentCounter > 0 % Folder already exist
      mkdir(strcat(dirPath,pathSep,numberedFolder));
  else % folder doesn't exist
      mkdir(strcat(dirPath,pathSep,folderName));
  end
end
end