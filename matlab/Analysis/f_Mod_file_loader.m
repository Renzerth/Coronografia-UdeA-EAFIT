function [images,Numfiles,folder_name,FileList] = f_Mod_file_loader(format,Flag)
% file_loader loads all image of format '.jpg' or '.png'
% in current directory. All files must have a sequential designation 
% of 'image_1', 'image_2'...'image_n' in order to be read.
%
%   Inputs:  None
%
%   Outputs: images: Cell - Structure with all pictures
%                Numfiles: Number of files with same format
%
%  Version: 1.2 for Matlab R2014b
%
%  Author: Juan Jos? Cadavid Mu?oz - EAFIT UNIVERSITY
%
%  Date: 08/02/2015
%
%% Retrieve Folder information
folder_name = uigetdir;
if ~ischar(folder_name); error('Folder was not specified.'); end
FolderDataStr = dir(fullfile(folder_name, ['*',format]));
FolderDataCll = struct2cell(FolderDataStr);
FileList = FolderDataCll(1,1:end)';
Numfiles = length(FileList);
%% File existance verification
if Numfiles==0
    error('No image files found');
end
%% Initializing
images = cell(1,Numfiles);
%% Load images
for index=1:Numfiles
    file_location = strcat(folder_name,'\',FileList{index,1});
    images{index} = imread(file_location);
end
%% Channel Conversion
if Flag
    for index=1:Numfiles
        if size(size(images{index}),2) ==3
            images{index} = double(rgb2gray(images{index}));
        else
            images{index} = double(images{index});
        end
    end
end
end