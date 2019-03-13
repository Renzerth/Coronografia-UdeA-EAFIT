function [folderCounter] = addFolderCount(dirPath,folderName,dirDelimiter)
% FolderName cannot contain an equal character as specified by dirDelimiter
dirInfo = dir(dirPath);
uniqueTokens = cellfun(@(string) strtok(string,dirDelimiter),{dirInfo.name},'UniformOutput',false); % Get all posible folder coincidences
folderCounter = sum([dirInfo(ismember(uniqueTokens,folderName)).isdir]); % Count valid foldernames coincidences to be added to the counter
end