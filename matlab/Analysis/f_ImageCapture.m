function f_ImageCapture(vid,dataDir,filename)

%% Target directory
cd(dataDir); 
snapsfldr = 'TestSnapshots';
ax = exist(snapsfldr, 'dir'); % 7 if folder exists, 0 if not
if ax ~= 7 % Create a folder if it doesn't exist
    mkdir(snapsfldr);
end
cd(snapsfldr);
%% Get frame data
SingleFrame = getsnapshot(vid);

if isfile(filename) 
 % warning('The file name already exists'); % File exists.
 newfileName = inputdlg(['The file name already exists, try another one' ...
 '(otherwise, the file may be overwritten): ']);
 imwrite(SingleFrame,newfileName{1});
else
  imwrite(SingleFrame,filename,'fmt','png');  % File does not exist.
end