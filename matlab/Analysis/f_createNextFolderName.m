function f_createNextFolderName(dirPath,folderName,dirDelim)
% Creates the next possible folder name for a given "folderName" in a
% specific directory
%
% Inputs:
%  dirPath: path of the folders that one wants to count
%  folderName: name of the root folder
%  dirDelimiter: splits the rootfolder and the distinct character between
%                folders (a string)
%
currentCounter = f_addFolderCount(dirPath,folderName,dirDelim);
if currentCounter > 0 % Folder already exist
    mkdir(strcat(folderName,'_',num2str(currentCounter)))
else % folder doesn't exist
    mkdir(folderName)
end

end