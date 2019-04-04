function SingleFrame = f_CaptureImage(vid,dataDir,filename,imgformat, ...
                                    pathSep,dirDelim,snapsfldr,previewBool)
                                  
% Inputs:
%   
%
% Outputs:



%% Target directory
createFoldSnap = 1; % Always creates a new folder
numberedFolderSnap = f_createNextFolderName(dataDir,snapsfldr, ...
                                          dirDelim,pathSep,createFoldSnap);    
snapsdir = strcat(dataDir,pathSep,numberedFolderSnap); % tests directory           

%% Get frame data
SingleFrame = getsnapshot(vid); % Take a picture
% isfile(filename) MATLAB 2017; exist(filename,'file'): MATLAB 2015
% imwrite(variables,directory+filename+extension)
fullFramePath = strcat(snapsdir,pathSep,filename,'.',imgformat); % Snapshot
                                                                 % path
if exist(fullFramePath, 'file') == 2 
  newfileName = inputdlg(['The file name already exists, try another ' ...
  'one (otherwise, the file may be overwritten): ']);
  imwrite(SingleFrame,strcat(snapsdir,pathSep,newfileName{1},'.',imgformat));
else
  imwrite(SingleFrame,fullFramePath); % File does not exist yet and will be
                                      % written.
end

%% Show frame
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium

% The numbers after 'position' were empirically obtained
figure('units','normalized','position',[1/10 1/10 1/3 1/2]);
imagesc(SingleFrame); % normalized
colorbar; title(['Camera image: ' filename]);

% The numbers after 'position' were empirically obtained
figure('units','normalized','position',[5/10 1/10 1/3 1/2]);
imagesc(log(im2double(SingleFrame))); 
colorbar; title(['LOG camera image: ' filename])

%% Open preview
if previewBool
  preview(vid);
end

end