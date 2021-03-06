function [SingleFrame,circShiftX,circShiftY] = f_CaptureImage(vid,dataDir,filename,imgformat, ...
pathSep,infoDelim,dirDelim,snapsfldr,previewBool,loghist,camera, ...
cameraPlane,exposure,format,fps,maskFig,plotMask)
                                  
% Inputs:
%   
%
% Outputs:

%% Target directory
createFoldSnap = 1; % Always creates a new folder
numberedFolderSnap = f_createNextFolderName(dataDir,snapsfldr, ...
                                          dirDelim,pathSep,createFoldSnap);    
snapsdir = strcat(dataDir,pathSep,numberedFolderSnap); % tests directory           
                                                    
%% Loop condition (for taking several snapshots)
loopCondition = 0; % Variable initialization
i = 1; % counter to store the snapshot's name

circShiftX = 0; % Initialization
circShiftY = 0; % Initialization

while loopCondition == 0  
  %% Open preview
  % In general, it is faster to open the preview first and then to 
  % implement the snapshot
  if previewBool
    % Original preview: previewHandle = preview(vid); 
    [N,D] = rat(exposure);
    tit = strcat('Preview of the',{' '},cameraPlane,{' '}, ...
    'camera (',camera,')',{' '},'[Exposure:',{' '},num2str(N),'/', ...
    num2str(D),';',{' '},'format:',{' '},num2str(format),';',{' '}, ...
    'fps:',{' '},num2str(fps),']');
    previewHandle = f_CustomPreview(vid,tit);
    wait(vid); % Waits until vid is not running or logging
  end

  %% Get frame data
  SingleFrame = getsnapshot(vid); % Take a picture

  %% Save the snapshot
  % imwrite(variables,directory+filename+extension)
  imwrite(SingleFrame,strcat(snapsdir,pathSep,filename,infoDelim, ...
                             num2str(i),'.',imgformat));

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
  [histImg,~] = imhist(SingleFrame);
  if loghist == 0 % truncated 
    histImg(histImg>1000) = 1000; % Truncates all the values above 1000
    h3 = figure; bar(histImg); % Shows a histogram of the snapshot
    ylabel('Counts');
  else % loghist == 1 % Logarithmic
    % SingleFrame = im2double(SingleFrame);
    % h3 = figure; imhist(log10(SingleFrame));
    h3 = figure; bar(histImg);
    set(gca,'yscale','log')
    ylabel('Counts (logscale)');
  end
  
  % Another option: log of the hist
  title('Histogram of the current snapshot');
  xlabel('Bins (gray levels)'); xlim([0 255]);
  
  %% Don't continue until the preview is closed
  
  
  %% 
  disp('Close the preview in order to continue the program"s execution.');
  if plotMask == 2 && ishandle(maskFig{1})% SLM plot
    % circShiftXp, circShiftYp : p: previous shift
    [circShiftXp, circShiftYp] = f_addSliderPositioning(maskFig{1},maskFig{2},previewHandle);
    circShiftX = circShiftXp + circShiftX;
    circShiftY = circShiftYp + circShiftY;
  end
  
  while ishandle(previewHandle)
    pause(0.1); % In order for the preview to start safely
    % Delay for user response: dead time of the processor
    % Allow idle instead of processing
  end
  
  % OLD method:
  % disp('Press any valid key to close the preview (find it!)');
  % pause;

  %% Ask to leave figures open or not
  answer1 = input('Do you want to take another snapshot? y/n: ','s');
  % returned always as a string  
  switch answer1
    case 'y' % Compare string. Loop repeats
    loopCondition = 0;  
    i = i + 1;
    case 'n' % Loop finishes
    loopCondition = 1;    
    otherwise % Loop finishes
    loopCondition = 1; 
    warning('Try a valid command, the program will assume a "no" as an answer');
  end
  
  %% Close current figures (if they exist)
  if ishandle(h1)
    close(h1);
  end
  if ishandle(h2)
    close(h2); 
  end
  if ishandle(h3)
    close(h3);
  end
  if ishandle(maskFig{1}) && ~strcmp(answer1,'y')
    close(maskFig{1});
  end
  
end

%% Close preview if open
if strcmp(vid.previewing,'on') 
  closepreview(vid); % Gets closed in case it was already opened
end

end