function [shiftY,shiftX,systemPupilPixelSize,mainDataCenter, ...
mainDataRadius] = f_selfCenterSLMmask(detectorPixelPitch, ...
systemPupilSIze,screenIndex,vid)

%% Set Reference
close all;
preview(vid);
referenceData = getsnapshot(vid); % Input reference, must be gray 'Y800' Format.

%% Lyot Intensity Feedback coordinates
[mainDataCenter, mainDataRadius,~] = f_findCircleShapedIntensity(referenceData,false);

%% System Pixel Size
[systemPupilPixelSize] = f_computePupilPixelSize(mainDataRadius, detectorPixelPitch, systemPupilSIze);

%% Mask Generation
shiftX = 0; shiftY =0; TC = 10;
[X,Y,rhoRadius, monitorSize] = f_MakeScreenCoords(screenIndex);
aspectRatio = monitorSize(1)/monitorSize(2);
[angularTranstion,~] = cart2pol(X - shiftX,Y - shiftY);
projectionMask = mat2gray((rhoRadius <= 0.25/aspectRatio).*angle(exp(1i*TC*(angularTranstion + pi/TC))));

%% Compute initial Radial Profile
[shiftY, shiftX] = centerMaskToSpot(referenceData, projectionMask, mainDataCenter, mainDataRadius, monitorSize, TC, vid);

%% Clear resources
clear vid;
end