function f_ImageCapture(vid,dataDir,filename,imgformat)

%% Target directory
cd(dataDir); % Move to the data directory
snapsfldr = 'TestSnapshots'; % Tests folder
ax = exist(snapsfldr, 'dir'); % 7 if folder exists, 0 if not
if ax ~= 7 % Create a folder if it doesn't exist
    mkdir(snapsfldr);
end
cd(snapsfldr); % Go to tests folder

%% Get frame data
SingleFrame = getsnapshot(vid); % Take a picture

if exist([filename imgformat], 'file') == 2 % isfile(filename) MATLAB 2017
  % warning('The file name already exists'); % File exists.
  newfileName = inputdlg(['The file name already exists, try another ' ...
  'one (otherwise, the file may be overwritten): ']);
  imwrite(SingleFrame,[newfileName{1} imgformat]);
else
  imwrite(SingleFrame,[filename imgformat]); % File does not exist yet and 
                                             % will be written.

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