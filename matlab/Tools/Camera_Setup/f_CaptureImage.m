function SingleFrame = f_CaptureImage(vid,dataDir,filename,imgformat, ...
                                      pathSep,snapsfldr)

%% Target directory
snapsdir = [dataDir pathSep snapsfldr]; % tests directory
ax = exist(snapsdir, 'dir'); % 7 if folder exists, 0 if not
if ax ~= 7 % Create a folder if it doesn't exist
    mkdir(snapsfldr);
end


%% Get frame data
SingleFrame = getsnapshot(vid); % Take a picture
% isfile(filename) MATLAB 2017
% imwrite(variables,directory+filename+extension)
fullFramePath = [snapsdir pathSep filename imgformat]; % Snapshot path
if exist(fullFramePath, 'file') == 2 
  newfileName = inputdlg(['The file name already exists, try another ' ...
  'one (otherwise, the file may be overwritten): ']);
  imwrite(SingleFrame,[snapsdir pathSep newfileName{1} imgformat]);
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
end