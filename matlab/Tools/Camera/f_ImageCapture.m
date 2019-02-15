function f_ImageCapture(vid,dataDir,filename)

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
figure; imagesc(SingleFrame); title(['Snapshot: ' filename]);
if isfile(filename) 
 % warning('The file name already exists'); % File exists.
 newfileName = inputdlg(['The file name already exists, try another one' ...
 '(otherwise, the file may be overwritten): ']);
 imwrite(SingleFrame,newfileName{1});
else
  imwrite(SingleFrame,filename,'fmt','png'); % File does not exist yet and 
                                             % will be written.
end