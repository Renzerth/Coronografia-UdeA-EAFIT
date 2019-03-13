folderName = 'DATA';
currentCounter = addFolderCount(pwd,folderName,'_');
if currentCounter > 0
    mkdir(strcat(folderName,'_',num2str(currentCounter)))
else
    mkdir(folderName)
end