function SingleFrame = f_CaptureImage(vid,dataDir,filename,imgformat, ...
                                    pathSep,infoDelim,dirDelim,snapsfldr,previewBool)
                                  
% Inputs:
%   
%
% Outputs:

%% Open preview
if previewBool
  preview(vid);
  pause(1); % In order for the preview to start safely
end

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
                                                                 
%% Loop condition (for taking several snapshots)
loopCondition = 0; % Variable initialization
i = 1; % counter to store the snapshot's name

while loopCondition == 0
  
  %% Save the snapshot
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
  h1 = figure('units','normalized','position',[1/10 1/10 1/3 1/2]);
  imagesc(SingleFrame); % normalized
  colorbar; title(['Camera image: ' filename]);

  % The numbers after 'position' were empirically obtained
  h2 = figure('units','normalized','position',[5/10 1/10 1/3 1/2]);
  imagesc(log(im2double(SingleFrame))); 
  colorbar; title(['LOG camera image: ' filename])
  
  %% Histogram
  h3 = figure; imhist(SingleFrame); % Shows a histogram of the snapshot
  title('Histogram of the current snapshot');
  xlabel('Gray levels'); ylabel('Counts');

  %% Ask to leave figures open or not
  answer = input('Do you want to take another snapshot? y/n: ','s'); % returned always as a string
  if strcmp(answer,'y') % Compare string. Loop repeats
    loopCondition = 0;  
    i = i + 1;
  elseif strcmp(answer,'n') % Loop finishes
    loopCondition = 1;    
  else % Loop finishes
    loopCondition = 1; 
    warning('Try a valid command, the program will assume a "no" as an answer')
  end
  
  %% Close current figures
  close(h1); close(h2); close(h3);
  
end

end