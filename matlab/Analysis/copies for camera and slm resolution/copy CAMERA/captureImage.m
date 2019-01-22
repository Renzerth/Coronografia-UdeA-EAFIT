function [SingleFrame] = captureImage(filename,MakeFolderFlag)
%% Set video source
camera = 'DMK23U445';
exposure = 1/40;
format = 'Y800 (1024x768)';
[vid,~] = CAMERA.selectCamera(camera,exposure,format);
vid.FramesPerTrigger = 1;
%% Make Target directory at CWD
if MakeFolderFlag
  targetFolder = '\Image_Data';
  FolderParentLevel = 1;
  TOOLS.makeParentFolder(FolderParentLevel,targetFolder);
end
%% Get frame data
SingleFrame = getsnapshot(vid);
imwrite(SingleFrame,filename);
%% Clean up
delete(vid);
end