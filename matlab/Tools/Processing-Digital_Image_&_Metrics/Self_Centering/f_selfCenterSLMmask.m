function [coorShiftY,coorShiftX,systemPupilPixelSize,mainDataCenter, ...
mainDataRadius] = f_selfCenterSLMmask(detectorPixelPitch, ...
systemPupilSIze,screenIndex,vid,tit)

%% Set Reference
close all;
previewHandler = f_CustomPreview(vid,tit);
% Old version: preview(vid);
wait(vid); % Waits until vid is not running or logging
referenceData = getsnapshot(vid); % Input reference, must be gray 'Y800'
                                  % format.

%% Lyot Intensity Feedback coordinates
[mainDataCenter, mainDataRadius,~] = ...
                          f_findCircleShapedIntensity(referenceData,false);

%% System Pixel Size
[systemPupilPixelSize] = f_computePupilPixelSize(mainDataRadius, ...
                         detectorPixelPitch,systemPupilSIze);

%% Mask Generation
shiftX = 0; shiftY =0; TC = 2; enablechange = true;
[X,Y,aspectRatio,monitorSize] = f_MakeScreenCoords(screenIndex, ...
                                                   enablechange);

%% Polar coordinates with the Aspect Ratio scaling
scaledY = Y/aspectRatio;
[angularTranstion,rhoRadius] = cart2pol(X - shiftX,scaledY - shiftY);

%% Mask generation
MaskProportion = 0.1; % Relative proportion of the size of the mask. Ref: 0.1
projectionMask = mat2gray((rhoRadius <= MaskProportion/aspectRatio).* ...
                 angle(exp(1i*TC*(angularTranstion + pi/TC))));

%% Find Vortex center coordinates
[circShiftY,circShiftX] = f_centerMaskToSpot(referenceData, ...
          projectionMask,mainDataCenter,mainDataRadius,monitorSize,TC,vid);
coorShiftX = monitorSize(1)/2 + circShiftX;
coorShiftY = monitorSize(2)/2 + circShiftY;

%% Close preview
if strcmp(vid.previewing,'on') 
  closepreview(vid); % Gets closed in case it was already opened
end
if ishandle(previewHandler)
  close(previewHandler);  
end
end