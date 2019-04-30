function numberedFolder = f_createNextFolderName(dirPath,folderName, ...
                                               dirDelim,pathSep,createFold)
% Creates the next possible folder name for a given "folderName" in a
% specific directory
%
% Inputs:
%  dirPath: path of the folders that one wants to count
%  folderName: name of the root folder
%  dirDelimiter: splits the rootfolder and the distinct character between
%                folders (a string)
%  pathSep: path separation
%  createFold: creates or not the folder (boolean)
%
%  Output:
%   numberedFolder: with the currentCounter (extended name)
%   numberedFolder = folderName if the folder doesn't exist,
%   so currentCounter = 0
%
currentCounter = f_addFolderCount(dirPath,folderName,dirDelim);

if createFold % Case where one wants to create a folder
  if currentCounter > 0 % Folder already exist
      numberedFolder = strcat(folderName,dirDelim,num2str(currentCounter));
      mkdir(strcat(dirPath,pathSep,numberedFolder));
  else % folder doesn't exist: currentCounter = 0
      mkdir(strcat(dirPath,pathSep,folderName));
      numberedFolder = folderName; % Root name is maintained
  end
elseif currentCounter > 1 % Cases where one wants to retrieve the current folder's name with the maximum counter
  numberedFolder = strcat(folderName,dirDelim,num2str(currentCounter-1)); % This is the existing number of folders
else % currentCounter = 1 and therefore this function returns folderName
   numberedFolder = folderName; % Root name is maintained
end

end