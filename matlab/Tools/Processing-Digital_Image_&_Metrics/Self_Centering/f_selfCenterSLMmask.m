function [coorShiftY,coorShiftX,systemPupilPixelSize,mainDataCenter, ...
mainDataRadius] = f_selfCenterSLMmask(detectorPixelPitch, ...
systemPupilSIze,screenIndex,vid)

%% Set Reference
close all;
preview(vid);
referenceData = getsnapshot(vid); % Input reference, must be gray 'Y800' Format.

%% Lyot Intensity Feedback coordinates
[mainDataCenter, mainDataRadius,~] = f_findCircleShapedIntensity(referenceData,false);

%% System Pixel Size
[systemPupilPixelSize] = f_computePupilPixelSize(mainDataRadius, ...
                         detectorPixelPitch,systemPupilSIze);

%% Mask Generation
shiftX = 0; shiftY =0; TC = 10; enablechange = true;
[X,Y,aspectRatio,monitorSize] = f_MakeScreenCoords(screenIndex,enablechange);

%% Polar coordinates with the Aspect Ratio scaling
scaledY = Y/aspectRatio;
[angularTranstion,rhoRadius] = cart2pol(X - shiftX,scaledY - shiftY);

%% Mask generation
projectionMask = mat2gray((rhoRadius <= 0.25/aspectRatio).* ...
                 angle(exp(1i*TC*(angularTranstion + pi/TC))));

%% Find Vortex center coordinates
[circShiftY,circShiftX] = f_centerMaskToSpot(referenceData,projectionMask, ...
                  mainDataCenter,mainDataRadius,monitorSize,TC,vid);
coorShiftX = monitorSize(1)/2 + circShiftX;
coorShiftY = monitorSize(2)/2 + circShiftY;
%% Clear resources
clear vid;
end