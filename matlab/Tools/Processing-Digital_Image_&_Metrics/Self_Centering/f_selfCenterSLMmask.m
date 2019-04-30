function [coorShiftX,coorShiftY,systemPupilPixelSize,mainDataCenter, ...
mainDataRadius] = f_selfCenterSLMmask(detectorPixelPitch, ...
systemPupilSIze,screenIndex,vid,tit)
% systemPupilSIze: used to calculate statistics

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
shiftX = 0; shiftY =0; TC = 10; enablechange = true;
[X,Y,aspectRatio,monitorSize] = f_MakeScreenCoords(screenIndex, ...
                                                   enablechange);
maskGen = @(rho,theta,r0,TC) mat2gray(rho <= r0).*angle(exp(1i*TC*(theta...
                                                                + pi/TC)));

%% Polar coordinates with the Aspect Ratio scaling
scaledY = Y/aspectRatio;
[angularTranstion,rhoRadius] = cart2pol(X - shiftX,scaledY - shiftY);

%% Mask generation
pupilRatio = 0.20; % Relative proportion of the size of the mask. Ref: 0.1
vortexTC10 = maskGen(rhoRadius,angularTranstion,pupilRatio/aspectRatio,10);
vortexTC02 = maskGen(rhoRadius,angularTranstion,pupilRatio/aspectRatio,2);
projectionMasks = {vortexTC10,vortexTC02};

%% Find Vortex center coordinates
[coorRelShiftX,coorRelShiftY] = f_centerMaskToSpot(referenceData, ...
          projectionMasks,mainDataCenter,mainDataRadius,monitorSize,TC,vid);
        
coorShiftX = monitorSize(1)/2 + coorRelShiftX;
coorShiftY = monitorSize(2)/2 + coorRelShiftY;

%% Close preview
if strcmp(vid.previewing,'on') 
  closepreview(vid); % Gets closed in case it was already opened
end
if ishandle(previewHandler)
  close(previewHandler);  
end
end