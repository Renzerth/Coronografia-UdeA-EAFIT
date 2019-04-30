function [folderCounter] = f_addFolderCount(dirPath,folderName,dirDelimiter)
% Returns the number of folders that have the same rootname "folderName"
% and that may differ by a delimiter "dirDelimiter" and something after
% this delimiter (normally a number)
%
% Inputs:
%  dirPath: path of the folders that one wants to count
%  folderName: name of the root folder that shouldn't have the dirDelimiter
%  dirDelimiter: splits the rootfolder and the distinct character between
%                folders (a string)
%
% Outputs:
%  folderCounter: number of folders that have the same rootname "folderName"
%
% Note: "FolderName" cannot contain an equal character as specified by
%       dirDelimiter

dirInfo = dir(dirPath);
uniqueTokens = cellfun(@(string) strtok(string,dirDelimiter), ...
              {dirInfo.name},'UniformOutput',false); % Get all the posible
                                                     % folder coincidences
folderCounter = sum([dirInfo(ismember(uniqueTokens,folderName)).isdir]);
% Count the valid foldernames coincidences to be added to the counter
end