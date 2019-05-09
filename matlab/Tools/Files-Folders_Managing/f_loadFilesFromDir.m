function [foundFiles] = f_loadFilesFromDir(fileFormat)
% loadFilesFromDir loads files of specified fileFormat
% from selected directory.
%   Inputs:  fileFormat (string)
%
%   Outputs: images: Cell - Structure with all pictures
%
%  Version: 1.2 for Matlab R2014b
%
%  Author: Juan Jose Cadavid Munoz - EAFIT UNIVERSITY
%
%  Date: 08/02/2015
%
%% Get Directory
folder_name = uigetdir;
if ~ischar(folder_name); error('Folder was not specified.'); end
FolderDataCll = struct2cell(dir(fullfile(folder_name, ['*',fileFormat])));
FileList = FolderDataCll(1,1:end)';
Numfiles = length(FileList);
%% File existance verification
if Numfiles==0
    error('No files found');
end
%% Initializing
foundFiles = cell(1,Numfiles);
for index=1:Numfiles
    file_location = strcat(folder_name,'/',FileList{index,1});
    foundFiles = load(file_location);
end
end