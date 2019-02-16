function ImageCapture(filename)
%% Get video source
vid = videoinput('tisimaq_r2013', 1, 'RGB24 (1024x768)');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
%% Set Camera properties
src.Brightness = 0;
src.Contrast = 0;
src.Exposure = 0.0556;
src.GainAuto = 'Off';
src.Gain = 34;
%% Target directory
MakeFolder = 0;
if MakeFolder
    TGTParentLevel = 1;
    targetFolder = '\Image_Data';
    currentDirectory = pwd;
    directoryElements = strsplit(currentDirectory,'\');
    targetDirectory = strjoin(directoryElements(1:end-TGTParentLevel),'\');
    ResPath = strcat(targetDirectory,targetFolder);
    ResFolder = date;
    
    if exist(strcat(ResPath,ResFolder),'dir') == 0
        mkdir(ResPath,ResFolder);
    end
    
end
%% Get frame data
SingleFrame = getsnapshot(vid);
imwrite(SingleFrame,filename);
%% Clean up
delete(vid);
% exit;
end